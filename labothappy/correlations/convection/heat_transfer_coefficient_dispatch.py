# -*- coding: utf-8 -*-
# Adapt the import paths to where the dispatchers live.
from labothappy.correlations.convection.BPHE_HTC import BPHE_HTC
from labothappy.correlations.convection.STHE_HTC import STHE_HTC
from labothappy.correlations.convection.pipe_HTC import pipe_HTC
from labothappy.correlations.convection.PCHE_HTC import PCHE_HTC
from labothappy.correlations.convection.TFHE_HTC import TFHE_HTC
from labothappy.correlations.convection.curved_pipe_HTC import curved_pipe_HTC
from labothappy.correlations.convection.SPHE_HTC import SPHE_HTC

_VALID_STATES = ('SP', 'EVAP', 'COND')
_VALID_SIDES = ('shell', 'tube')


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
           mu, k, Pr, T (+ T_in, p_in, p_other_in for 'Bennett_Chen'); geom holds D_c, L_c, alpha,
           type_channel ('Straight' / 'Zigzag') and the equivalent tube (Dh, D_t_i, L_t, n_t).
    HCHE : coil side SP 'Sigalotti' / EVAP 'Zhao','Nariai' / COND 'Zakeralhoseini'; geom: Dh, R_curv (or D_coil, P_coil), L.
    SCHE : same as HCHE (curved pipe correlations on the elliptical channel, D = D_h).
    SPHE : SP 'Picon_Nunez','Rajavel','Curved_rect_duct','Rect_duct'; geom: Dh, R_curv, ch_w, H (or ch_b), L.
    TFHE ('Tube&Fins') : cond['side'] (or geom['side']) = 'fin' (gas side, single phase: Pismennyi, REDUCED
           coefficient referred to the total outer area) or 'tube' (pipe correlations). Fin side: m_dot (WHOLE
           side), rho, mu, k, Pr; geom holds the tube and fin geometry (D_f_o, t_f, f_density, k_f, w, ...)."""
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
