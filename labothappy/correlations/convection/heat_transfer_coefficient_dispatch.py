# -*- coding: utf-8 -*-
# Adapt the import paths to where the dispatchers live.
from labothappy.correlations.convection.BPHE_HTC import BPHE_HTC
from labothappy.correlations.convection.STHE_HTC import STHE_HTC
from labothappy.correlations.convection.pipe_HTC import pipe_HTC
from labothappy.correlations.convection.PCHE_HTC import PCHE_HTC
from labothappy.correlations.convection.TFHE_HTC import TFHE_HTC
from labothappy.correlations.convection.curved_pipe_HTC import curved_pipe_HTC
from labothappy.correlations.convection.SPHE_HTC import SPHE_HTC
from labothappy.correlations.convection.bare_tube_bank_HTC import bare_tube_bank_HTC

_VALID_STATES = ('SP', 'EVAP', 'COND')
_VALID_SIDES = ('shell', 'tube')

# Values of geom['finned'] that mean "no fins" (bare tubes)
_NOT_FINNED = (0, False, '0', 'no', 'No', 'False', 'false')


def _is_bare_tube_bank(htx, side, geom):
    """True for the outer side of a TFHE whose tubes have no fins (geom['finned'] = 0 / False):
    bank of bare tubes in crossflow, routed to bare_tube_bank_HTC (2026-11-27)."""
    return htx in ('TFHE', 'Tube&Fins') and side == 'fin' and geom.get('finned', 1) in _NOT_FINNED


def _choices(correlation, state):
    """A dict {state: name} goes straight through; a plain name applies to the current state."""
    return correlation if isinstance(correlation, dict) else {state: correlation}


def heat_transfer_coefficient(AS, geom, cond, correlation):
    """cond is a dictionary built by the model (_cond_dict): it holds G and flow_state,
    and every field of CellConditions by name.

    BPHE : G, p, x, q, m_dot; mu, k, Pr (+ mu_w) for single phase.
    STHE : cond['side'] (or geom['side']) selects the side of the shell & tube exchanger.
      'shell' : SP: m_dot (WHOLE side), mu, k, Pr, cp (+ mu_w optional);
                EVAP / COND: m_dot, p, x, q, T_w;
                geom holds the shell & tube parameters (D_t_o, PR, D_s, ...).
      'tube'  : SP / EVAP / COND, inside the tubes (pipe correlations). G (mass flux in ONE tube,
                or m_dot of the whole side), p, x, q, T_w, mu, k, Pr (+ mu_w, rho, cp, h_in, h_out
                for sCO2); geom holds D_t_i (or D_t_o, t_t) and L_t.
    PCHE : SP (PCHE correlations 'Straight', 'Lee', 'Bennett_Chen', or any pipe correlation on the
           equivalent tube D = D_h) / EVAP / COND (pipe correlations). G (mass flux in ONE channel),
           mu, k, Pr, T (+ T_in, p_in, p_other_in for 'Bennett_Chen'); geom holds D_ch (or g_f, d_f
           for ch_section 'rectangular'), L_ch, alpha, type_channel ('Straight' / 'Zigzag' / 'S_fin')
           and the equivalent tube (Dh, D_t_i, L_t, n_t).
    HCHE : coil side SP 'Sigalotti' / EVAP 'Zhao','Nariai' / COND 'Zakeralhoseini'; geom: Dh, R_curv (or D_coil, P_coil), L.
    SCHE : same as HCHE (curved pipe correlations on the elliptical channel, D = D_h).
    SPHE : SP 'Picon_Nunez','Rajavel','Curved_rect_duct','Rect_duct'; geom: Dh, R_curv, b (channel gap), H, L.
    TFHE ('Tube&Fins') : cond['side'] (or geom['side']) = 'fin' (gas side, single phase: Pismennyi, REDUCED
           coefficient referred to the total outer area) or 'tube' (pipe correlations). Fin side: m_dot (WHOLE
           side), rho, mu, k, Pr; geom holds the tube and fin geometry (D_f_o, t_f, f_density, k_f, w, ...).
           geom['finned'] = 0 / False: outer side = bank of BARE tubes in crossflow (bare_tube_bank_HTC:
           SP 'Zukauskas' / 'ESDU_73031', EVAP 'Shah', COND 'Rose_Kern' / 'Nusselt_single_tube')."""
    htx = geom.get('HTX_Type')

    if htx in ('BPHE', 'STHE', 'PCHE', 'TFHE', 'Tube&Fins', 'HCHE', 'SCHE', 'SPHE'):
        state = cond.get('flow_state')
        if state not in _VALID_STATES:
            raise ValueError(f"cond['flow_state'] = {state!r}: expected one of {_VALID_STATES}.")
        choices = _choices(correlation, state)

        if htx == 'BPHE':
            return float(BPHE_HTC(AS, geom, cond, state, choices))

        if htx == 'PCHE':
            return float(PCHE_HTC(AS, geom, cond, state, choices))

        if htx in ('HCHE', 'SCHE'):                     # HCHE: coil side only (shell side not implemented); SCHE: both sides
            return float(curved_pipe_HTC(AS, geom, cond, state, choices))

        if htx == 'SPHE':
            return float(SPHE_HTC(AS, geom, cond, state, choices))

        if htx in ('TFHE', 'Tube&Fins'):      # cond['side'] / geom['side'] = 'fin' or 'tube'
            if _is_bare_tube_bank(htx, cond.get('side', geom.get('side')), geom):
                return float(bare_tube_bank_HTC(AS, geom, cond, state, choices))   # finned = 0
            return float(TFHE_HTC(AS, geom, cond, state, choices))

        side = cond.get('side', geom.get('side'))
        if side == 'shell':
            return float(STHE_HTC(AS, geom, cond, state, choices))   # SP, EVAP, COND
        if side == 'tube':
            return float(pipe_HTC(AS, geom, cond, state, choices))
        raise ValueError(f"STHE: cond['side'] (or geom['side']) = {side!r}: expected one of {_VALID_SIDES}.")

    raise NotImplementedError(f"heat_transfer_coefficient: HTX_Type {htx!r} is not handled "
                              "by the three-level scheme (the legacy registry is disabled).")
    
    


# =============================================================================
# NAME OF THE AUTOMATIC CHOICE (used by hex_MB to fix one correlation per phase)
# =============================================================================

def auto_htc_correlation(AS, geom, cond):
    """Name of the correlation that heat_transfer_coefficient(AS, geom, cond, None) would choose
    for the state in cond, WITHOUT computing the coefficient.

    Purpose: the model (hex_MB) calls it once per side and per phase, at a representative state of
    the phase, and then gives this name explicitly to every cell of the phase. The automatic rules
    (select_* functions of each dispatcher) are not duplicated: they are called here.

    Arguments
    ---------
    AS   : CoolProp AbstractState set at the state of cond (h, p)
    geom : geometry dictionary of the side (as for heat_transfer_coefficient)
    cond : condition dictionary (as for heat_transfer_coefficient): flow_state 'SP', 'EVAP' or
           'COND', and the keys read by the select_* functions (G, mu, p, x, q, T_w, ...)

    Returns
    -------
    str : name of the correlation, valid as an explicit choice for heat_transfer_coefficient.
    """
    # Selection functions of the dispatchers (imported here to keep this helper self-contained)
    from labothappy.correlations.convection.BPHE_HTC import (
        select_plate_HTC, select_plate_evaporation_HTC, select_plate_condensation_HTC)
    from labothappy.correlations.convection.STHE_HTC import (
        select_shell_HTC, select_shell_evaporation_HTC, select_shell_condensation_HTC)
    from labothappy.correlations.convection.pipe_HTC import (
        select_pipe_HTC, select_pipe_evaporation_HTC, select_pipe_condensation_HTC)
    from labothappy.correlations.convection.PCHE_HTC import select_PCHE_HTC
    from labothappy.correlations.convection.PCHE_HTC import _equivalent_tube as pche_equivalent_tube
    from labothappy.correlations.convection.TFHE_HTC import select_fin_HTC
    from labothappy.correlations.convection.curved_pipe_HTC import (
        select_evap_HTC as select_curved_evap_HTC)
    from labothappy.correlations.convection.curved_pipe_HTC import _AUTO as CURVED_AUTO
    from labothappy.correlations.convection.SPHE_HTC import _AUTO as SPHE_AUTO

    htx = geom.get('HTX_Type')
    state = cond.get('flow_state')
    if state not in _VALID_STATES:
        raise ValueError(f"cond['flow_state'] = {state!r}: expected one of {_VALID_STATES}.")
    side = cond.get('side', geom.get('side'))

    # Tube side of STHE / TFHE, and equivalent tubes: in-tube (pipe) rules
    pipe_rules = {'SP': select_pipe_HTC,
                  'EVAP': select_pipe_evaporation_HTC,
                  'COND': select_pipe_condensation_HTC}

    if htx == 'BPHE':
        if state == 'SP':
            return select_plate_HTC(geom, cond)
        if state == 'EVAP':
            return select_plate_evaporation_HTC(geom, cond)
        return select_plate_condensation_HTC(geom, cond)

    if htx == 'STHE':
        if side == 'tube':
            return pipe_rules[state](AS, geom, cond)
        if state == 'SP':
            return select_shell_HTC(geom, cond)
        if state == 'EVAP':
            return select_shell_evaporation_HTC(geom, cond)
        return select_shell_condensation_HTC(geom, cond)

    if htx in ('TFHE', 'Tube&Fins'):
        if side == 'tube':
            return pipe_rules[state](AS, geom, cond)
        if _is_bare_tube_bank(htx, side, geom):
            from labothappy.correlations.convection.bare_tube_bank_HTC import (
                select_bare_tube_bank_HTC)
            return select_bare_tube_bank_HTC(geom, cond, state)
        return select_fin_HTC(geom, cond)       # finned side: same correlation for every state

    if htx == 'PCHE':
        if state == 'SP':
            return select_PCHE_HTC(geom, cond)
        # Two phase: equivalent tube; straight pipe rules without a radius of curvature,
        # curved pipe rules with one (same routing as PCHE_HTC._pipe)
        tube = pche_equivalent_tube(geom)
        if 'R_curv' not in geom:
            return pipe_rules[state](AS, tube, cond)
        if state == 'EVAP':
            return select_curved_evap_HTC(AS, tube, cond)
        return CURVED_AUTO['COND']

    if htx in ('HCHE', 'SCHE'):
        if state == 'EVAP':
            return select_curved_evap_HTC(AS, geom, cond)
        return CURVED_AUTO[state]

    if htx == 'SPHE':
        return SPHE_AUTO[state]

    raise NotImplementedError(f"auto_htc_correlation: HTX_Type {htx!r} is not handled.")


def _htc_case_tables(htx, side, state, curved, bare=False):
    """Registry tables of the heat transfer correlations of a case, in priority order, as
    (dispatcher module, table {name: function}).
    htx: HTX_Type; side: 'shell' / 'tube' (STHE), 'fin' / 'tube' (TFHE), None otherwise;
    state: 'SP', 'EVAP', 'COND'; curved: True if the channel has a radius of curvature;
    bare: True for the outer side of a TFHE without fins (bank of bare tubes)."""
    from labothappy.correlations.convection import BPHE_HTC, STHE_HTC, pipe_HTC, PCHE_HTC, TFHE_HTC
    from labothappy.correlations.convection import curved_pipe_HTC, SPHE_HTC
    from labothappy.correlations.convection import bare_tube_bank_HTC as bare_bank

    pipe_tables = {'SP': pipe_HTC._PIPE_SP_HTC, 'EVAP': pipe_HTC._PIPE_EVAP_HTC,
                   'COND': pipe_HTC._PIPE_COND_HTC}

    if htx == 'BPHE':
        return [(BPHE_HTC, {'SP': BPHE_HTC._BPHE_SP_HTC, 'EVAP': BPHE_HTC._BPHE_EVAP_HTC,
                            'COND': BPHE_HTC._BPHE_COND_HTC}[state])]
    if htx == 'STHE' and side == 'tube':
        return [(pipe_HTC, pipe_tables[state])]
    if htx == 'STHE':
        return [(STHE_HTC, {'SP': STHE_HTC._STHE_SP_HTC, 'EVAP': STHE_HTC._STHE_EVAP_HTC,
                            'COND': STHE_HTC._STHE_COND_HTC}[state])]
    if htx in ('TFHE', 'Tube&Fins') and side == 'tube':
        return [(pipe_HTC, pipe_tables[state])]
    if htx in ('TFHE', 'Tube&Fins') and bare:
        return [(bare_bank, bare_bank._TABLES[state])]
    if htx in ('TFHE', 'Tube&Fins'):
        return [(TFHE_HTC, TFHE_HTC._TFHE_FIN_HTC)]
    if htx == 'PCHE' and state == 'SP':
        return [(PCHE_HTC, PCHE_HTC._PCHE_SP_HTC)]
    if htx == 'PCHE':
        # two phase: PCHE correlations, then the tube correlations of the equivalent tube
        # (curved pipe ones if the channel has a radius of curvature, straight pipe ones otherwise)
        own = {'EVAP': PCHE_HTC._PCHE_EVAP_HTC, 'COND': PCHE_HTC._PCHE_COND_HTC}[state]
        if curved:
            return [(PCHE_HTC, own), (curved_pipe_HTC, curved_pipe_HTC._HCHE[state])]
        return [(PCHE_HTC, own), (pipe_HTC, pipe_tables[state])]
    if htx in ('HCHE', 'SCHE'):
        return [(curved_pipe_HTC, curved_pipe_HTC._HCHE[state])]
    if htx == 'SPHE':
        return [(SPHE_HTC, SPHE_HTC._SPHE[state])]
    raise NotImplementedError(f"htc_candidates: HTX_Type {htx!r} is not handled.")


def htc_candidates(AS, geom, cond):
    """Names of every registered heat transfer correlation usable for the case of cond
    (heat exchanger type, side, flow state), in priority order: the rule-based automatic choice
    (auto_htc_correlation) first, then the others in the order of the registry tables.
    A function registered under several names (aliases) appears once."""
    htx = geom.get('HTX_Type')
    state = cond.get('flow_state')
    side = cond.get('side', geom.get('side'))
    bare = _is_bare_tube_bank(htx, side, geom)
    tables = [table for _, table in _htc_case_tables(htx, side, state, 'R_curv' in geom, bare)]

    # Rule-based choice first, then the registry order, each function once
    names = [auto_htc_correlation(AS, geom, cond)]
    functions = []
    for table in tables:
        for name, function in table.items():
            if name in names or function in functions:
                continue
            names.append(name)
            functions.append(function)
    return names


def choose_htc_correlation(AS, geom, cond):
    """Automatic choice of the heat transfer correlation for the state of cond, by TRYING every
    registered candidate (htc_candidates) at this state: first the correlations whose parameters
    are all given, then the ones inside their validity range, then the ones established for the
    fluid (relaxed in the reverse order; see correlation_warnings.choose_by_trial).

    Arguments: as heat_transfer_coefficient (AS set at the state of cond).
    Returns: name of the correlation (a valid explicit choice for heat_transfer_coefficient).
    Raises ValueError if no candidate can be computed at this state."""
    import CoolProp.CoolProp as CP
    from labothappy.correlations.convection.correlation_warnings import choose_by_trial

    # The correlations may move AS (saturation properties): it is set back before every trial
    h, p = AS.hmass(), AS.p()

    def evaluate(name):
        AS.update(CP.HmassP_INPUTS, h, p)
        return heat_transfer_coefficient(AS, geom, cond, name)

    side = cond.get('side', geom.get('side', ''))
    label = f"{geom.get('HTX_Type')} {side} {cond.get('flow_state')} HTC"
    names = htc_candidates(AS, geom, cond)
    chosen, rejected = choose_by_trial(names, evaluate, label)
    AS.update(CP.HmassP_INPUTS, h, p)
    if chosen is None:
        raise ValueError(f"{label}: no candidate correlation can be computed ({rejected}).")
    return chosen


# Keys that only route the call (type, side, state): not tested by the required-parameter search
_ROUTING_KEYS = ('HTX_Type', 'side', 'flow_state', 'phase')


def htc_required_parameters(AS, geom, cond, correlation=None):
    """Parameters needed by a heat transfer correlation at the state of cond.

    Each key of geom and cond is removed in turn and the correlation is evaluated again
    (see correlation_warnings.parameters_by_removal): a key whose removal makes the evaluation
    fail, or makes the correlation assume a value, is required; a key whose removal only
    changes the value is optional. Routing keys (_ROUTING_KEYS) are not tested.

    Arguments
    ---------
    AS          : CoolProp AbstractState set at the state of cond (h, p)
    geom, cond  : as for heat_transfer_coefficient
    correlation : name of the correlation (or dict {state: name}); None or 'auto' = the
                  automatic choice (choose_htc_correlation)

    Returns
    -------
    dict {'correlation': name,
          'geom': {'required': [...], 'optional': [...]},
          'cond': {'required': [...], 'optional': [...]}}
    """
    import CoolProp.CoolProp as CP
    from labothappy.correlations.convection.correlation_warnings import parameters_by_removal

    if correlation is None or correlation == 'auto':
        correlation = choose_htc_correlation(AS, geom, cond)

    # The correlations may move AS (saturation properties): it is set back before every trial
    h, p = AS.hmass(), AS.p()

    def evaluate(trial_geom, trial_cond):
        AS.update(CP.HmassP_INPUTS, h, p)
        return heat_transfer_coefficient(AS, trial_geom, trial_cond, correlation)

    try:
        result = parameters_by_removal(geom, cond, evaluate, skip=_ROUTING_KEYS)
    finally:
        AS.update(CP.HmassP_INPUTS, h, p)
    result['correlation'] = correlation
    return result

def htc_correlation_parameters(htx, side, state, name=None, curved=False, bare=False):
    """Model parameters that a heat transfer correlation needs IN ADDITION to the parameters of the
    geometry (tables _REQUIRED_PARAMETERS of the dispatchers of each type), for one case.

    Arguments
    ---------
    htx    : HTX_Type ('BPHE', 'STHE', 'TFHE', 'PCHE', 'SPHE', 'SCHE', 'HCHE')
    side   : 'shell' / 'tube' (STHE), 'fin' / 'tube' (TFHE), None for the other types
    state  : 'SP', 'EVAP' or 'COND'
    name   : correlation name; None or 'auto' = automatic choice: the parameters of the
             registered correlation of the case that needs the FEWEST are returned (with the
             trial choice, the solve can always fall back on it)
    curved : True if the channel has a radius of curvature (PCHE two phase: curved pipe tables)
    bare   : True for the outer side of a TFHE without fins (geom['finned'] = 0): bare tube bank

    Returns
    -------
    tuple of parameter names. A name that is not registered for the case returns () (the
    dispatcher raises the error at the first evaluation).
    """
    tables = _htc_case_tables(htx, side, state, curved, bare)
    if name is None or name == 'auto':
        fewest = None
        for module, table in tables:
            for registered in table:
                needed = module.required_parameters(registered)
                if fewest is None or len(needed) < len(fewest):
                    fewest = needed
        return tuple(fewest or ())
    for module, table in tables:
        if name in table:
            return tuple(module.required_parameters(name))
    return ()
