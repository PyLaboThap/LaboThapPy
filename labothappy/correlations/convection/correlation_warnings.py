# -*- coding: utf-8 -*-
"""
Categories of the warnings raised by the correlations (heat transfer and pressure drop).

All three derive from UserWarning: for the user nothing changes (same messages, same filters).
They let the automatic choice of the dispatchers (choose_htc_correlation / choose_dp_correlation)
know WHY a correlation warned when it is tried at a given state:

    AssumedParameterWarning : a parameter the correlation needs was not given and a value was
                              assumed (default, substitute key, simplified geometry)
    FluidWarning            : the fluid is not one the correlation was established for
    ValidityRangeWarning    : a condition or a geometry value is outside the validity range of the
                              correlation (Re, Pr, x, G, q, p, angle, diameter, regime, ...)
"""


class AssumedParameterWarning(UserWarning):
    """A needed parameter was missing: an assumed value was used."""


class FluidWarning(UserWarning):
    """The fluid is not one the correlation was established for."""


class ValidityRangeWarning(UserWarning):
    """A condition or a geometry value is outside the validity range of the correlation."""


def choose_by_trial(names, evaluate, label):
    """Automatic choice among candidate correlations, by trying each one at a representative state.

    Each candidate is evaluated once (evaluate(name) -> value) and the warnings it raises are
    recorded. Ranking (first = best), as decided by Basile (2026-11-25):
      0. the correlation must run: an error or a non-finite / negative value discards it
      1. parameters: no AssumedParameterWarning (every parameter it needs was given)
      2. validity range: no ValidityRangeWarning
      3. fluid: no FluidWarning
    The criteria are relaxed in the reverse order (fluid first, then range, then parameters).
    Between candidates of the same rank, the order of `names` is kept (the caller puts the
    previous rule-based choice first).

    Arguments
    ---------
    names    : candidate names, in priority order
    evaluate : function name -> value (heat transfer coefficient or pressure drop)
    label    : text naming the case, for the warning

    Returns
    -------
    (name, rejected) : chosen name (None if no candidate runs) and a dictionary
                       name -> reason for the candidates not chosen.
    """
    import math
    import warnings

    ranked = []                     # (rank tuple, position, name)
    rejected = {}
    for position, name in enumerate(names):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                value = evaluate(name)
            except Exception as e:
                rejected[name] = f"error: {type(e).__name__}"
                continue
        if value is None or not math.isfinite(value) or value < 0:
            rejected[name] = f"value {value!r}"
            continue
        # Criteria not met (True = problem), in the order of importance
        assumed = any(issubclass(w.category, AssumedParameterWarning) for w in caught)
        out_of_range = any(issubclass(w.category, ValidityRangeWarning) for w in caught)
        other_fluid = any(issubclass(w.category, FluidWarning) for w in caught)
        ranked.append(((assumed, out_of_range, other_fluid), position, name))

    if not ranked:
        return None, rejected
    ranked.sort()
    best_rank, _, best = ranked[0]
    for rank, _, name in ranked[1:]:
        texts = ("assumed parameter", "out of range", "other fluid")
        problems = [text for flag, text in zip(rank, texts) if flag]
        rejected[name] = ", ".join(problems) if problems else "lower priority"

    # Warn when the best candidate does not meet every criterion
    relaxed = [text for flag, text in zip(best_rank, ("assumed parameter", "out of validity range",
                                                       "fluid not tested")) if flag]
    if relaxed:
        warnings.warn(f"{label}: automatic choice '{best}' ({', '.join(relaxed)}): no candidate "
                      "meets every criterion.", stacklevel=3)
    return best, rejected


def parameters_by_removal(geom, cond, evaluate, skip=()):
    """Parameters a correlation needs, found by REMOVING them and evaluating it again.

    The correlation is evaluated once with every parameter (reference value). A removal "fails"
    when the evaluation raises an error (KeyError, ...) or a NEW AssumedParameterWarning (a
    default value is then used).
      1. Each key of geom and cond is removed alone:
           failure                    -> REQUIRED (no other given key can replace it)
           value changed, no failure  -> OPTIONAL (used, but not needed to run)
      2. The other keys are removed one after the other, CUMULATIVELY, starting from the end of
         the dictionaries (hex_MB puts the computed keys after the user ones): a key whose
         removal then fails or changes the value is kept -> ALTERNATIVE (needed, but another
         given key could give the same quantity, e.g. Dh / D_t_i / D_ch; one of them is needed).
    REQUIRED + ALTERNATIVE is a minimal set of keys giving the same value.
    The test sees everything the dispatcher and the correlation read, even through copies of
    geom / cond.

    Arguments
    ---------
    geom     : geometry dictionary of the side
    cond     : condition dictionary of the state
    evaluate : function (geom, cond) -> value (heat transfer coefficient or pressure drop)
    skip     : keys not tested (routing keys such as HTX_Type, flow_state, side)

    Returns
    -------
    dict {'geom': {'required': [...], 'alternative': [...], 'optional': [...]},
          'cond': {'required': [...], 'alternative': [...], 'optional': [...]}}
    (keys in the order of the dictionaries)
    """
    import math
    import warnings

    def trial(trial_geom, trial_cond):
        """(failed, value) of one evaluation."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                value = evaluate(trial_geom, trial_cond)
            except Exception:
                return True, None
        assumed = {str(w.message) for w in caught
                   if issubclass(w.category, AssumedParameterWarning)}
        return bool(assumed - reference_assumed), value

    def same_value(value):
        if value == reference:
            return True
        if not (math.isfinite(value) and math.isfinite(reference)):
            return False
        return math.isclose(value, reference, rel_tol=1e-9, abs_tol=0.0)

    # Reference evaluation with every parameter (an error here is raised to the caller)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        reference = evaluate(geom, cond)
    reference_assumed = {str(w.message) for w in caught
                         if issubclass(w.category, AssumedParameterWarning)}

    result = {}
    for part in ('geom', 'cond'):
        result[part] = {'required': [], 'alternative': [], 'optional': []}
    tested = [('geom', key) for key in geom if key not in skip]
    tested += [('cond', key) for key in cond if key not in skip]

    # 1. Each key removed alone
    unused = []
    for part, key in tested:
        trial_geom = dict(geom)
        trial_cond = dict(cond)
        if part == 'geom':
            del trial_geom[key]
        else:
            del trial_cond[key]
        failed, value = trial(trial_geom, trial_cond)
        if failed:
            result[part]['required'].append(key)
        elif not same_value(value):
            result[part]['optional'].append(key)
        else:
            unused.append((part, key))

    # 2. Keys not used alone: removed cumulatively (from the end), kept if still needed
    kept_geom = dict(geom)
    kept_cond = dict(cond)
    for part, key in reversed(unused):
        trial_geom = dict(kept_geom)
        trial_cond = dict(kept_cond)
        if part == 'geom':
            del trial_geom[key]
        else:
            del trial_cond[key]
        failed, value = trial(trial_geom, trial_cond)
        if failed or not same_value(value):
            result[part]['alternative'].append(key)
        else:
            kept_geom = trial_geom
            kept_cond = trial_cond

    # Order of the dictionaries
    for part, data in (('geom', geom), ('cond', cond)):
        order = list(data)
        result[part]['alternative'].sort(key=order.index)
    return result
