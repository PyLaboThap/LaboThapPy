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
    
    