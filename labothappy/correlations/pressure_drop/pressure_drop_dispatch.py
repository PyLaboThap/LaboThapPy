# Adapt the import paths to where the dispatchers live.
from labothappy.correlations.pressure_drop.BPHX_DP import BPHE_DP, BPHX_port_DP
from labothappy.correlations.pressure_drop.STHE_DP import STHE_DP
from labothappy.correlations.pressure_drop.pipe_DP import pipe_DP
from labothappy.correlations.pressure_drop.TFHE_DP import TFHE_DP
from labothappy.correlations.pressure_drop.curved_pipe_DP import curved_pipe_DP
from labothappy.correlations.pressure_drop.SPHE_DP import SPHE_DP
from labothappy.correlations.pressure_drop.PCHE_DP import PCHE_DP

_VALID_STATES = ('SP', 'EVAP', 'COND')
_VALID_SIDES = ('shell', 'tube')
_CURVED_PIPE_TYPES = ('SCHE', 'SPCE', 'HCHE')       # spiral channel / helical coil: curved pipes


def pressure_drop(AS, geom, cond, correlation):
    """Frictional pressure drop [Pa] of one side of one cell, over the FULL flow length of the side
    (the model scales it to the cell). Port losses are not included: see port_pressure_drop.

    cond is a dictionary built by the model (_cond_dict): it holds G, phase and flow_state, and every
    field of CellConditions by name; for a single-phase pressure drop it also holds mu, D / rho (density).

    The flow state of the pressure drop is cond['flow_state'] ('SP', 'EVAP', 'COND'), except for
    a 'vapor-wet' cell: its bulk is superheated vapour, so its friction is single phase ('SP').
    `correlation` is a name for the current flow state, None / 'auto', or a dict
    {'SP': name, 'EVAP': name, 'COND': name}.

    BPHE : plate correlations (BPHE_DP).
    STHE : cond['side'] (or geom['side']) = 'shell' (SP, EVAP, COND: STHE_DP) or 'tube' (pipe_DP).
    TFHE : cond['side'] (or geom['side']) = 'fin' (gas side: Pismennyi bank pressure drop,
           TFHE_DP) or 'tube' (pipe_DP).
    PCHE : PCHE_DP (SP zigzag: 'Lee_Jiang', auto for geom['type_channel'] = 'Zigzag'; otherwise pipe_DP
           names on the equivalent tube: geom Dh / D_t_i, L, K; cond G of ONE channel); no port losses.
    SPHE : SPHE_DP (spiral plate: SP 'Picon_Nunez', 'Curved_rect_duct', 'Rect_duct', 'Sigalotti', 'Ito';
           two phase: curved pipe correlations on D_h and the cell R_curv).
    SCHE (SPCE), HCHE : curved_pipe_DP (curved pipes: SP 'Sigalotti', 'Ito'; EVAP 'Cioncolini', 'Ruffell';
           COND 'Fsadni', 'Cioncolini', 'Ruffell'; or pipe names on the straight equivalent tube)."""
    htx = geom.get('HTX_Type')

    if htx in ('BPHE', 'STHE', 'PCHE', 'TFHE', 'SPHE') + _CURVED_PIPE_TYPES:
        state = 'SP' if cond.get('phase') == 'vapor-wet' else cond.get('flow_state')
        if state not in _VALID_STATES:
            raise ValueError(f"Pressure drop: flow state {state!r}, expected one of {_VALID_STATES}.")
        choices = correlation if isinstance(correlation, dict) else {state: correlation}

        if htx == 'BPHE':
            return float(BPHE_DP(AS, geom, cond, state, choices))

        if htx == 'PCHE':      # PCHE dispatcher: zigzag 'Lee_Jiang', else pipe_DP on the equivalent tube
            return float(PCHE_DP(AS, geom, cond, state, choices))

        if htx == 'SPHE':      # spiral plate: rectangular curved channel (G of one channel, R_curv of the cell)
            return float(SPHE_DP(AS, geom, cond, state, choices))

        if htx in _CURVED_PIPE_TYPES:      # curved pipe on the equivalent tube D = D_h (R_curv of the cell)
            return float(curved_pipe_DP(AS, geom, cond, state, choices))

        if htx == 'TFHE':      # cond['side'] / geom['side'] = 'fin' or 'tube'
            return float(TFHE_DP(AS, geom, cond, state, choices))

        side = cond.get('side', geom.get('side'))
        if side == 'shell':
            return float(STHE_DP(AS, geom, cond, state, choices))
        if side == 'tube':
            return float(pipe_DP(AS, geom, cond, state, choices))
        raise ValueError(f"STHE: cond['side'] (or geom['side']) = {side!r}: expected one of {_VALID_SIDES}.")

    raise NotImplementedError(f"pressure_drop: heat exchanger type {htx!r} is not handled yet by the "
                              "three-level scheme ('BPHE', 'STHE', 'PCHE', 'TFHE', 'SPHE', 'SCHE', 'HCHE').")


def port_pressure_drop(geom, cond):
    """Port / manifold pressure loss [Pa] of one side, for the WHOLE exchanger (not per cell).
    cond must hold m_dot and D (density, evaluated at the inlet state), or G_port; geom may hold
    D_port and n_passes. Returns 0 (with a warning) when no port size is available.
    All other types (STHE, PCHE, TFHE, SPHE, SCHE, HCHE): 0 (nozzle, header and head losses are not modelled)."""
    if geom.get('HTX_Type') == 'BPHE':
        return float(BPHX_port_DP(geom, cond))
    return 0.0


# =============================================================================
# NAME OF THE AUTOMATIC CHOICE (used by hex_MB to fix one correlation per phase)
# =============================================================================

def auto_dp_correlation(AS, geom, cond):
    """Name of the correlation that pressure_drop(AS, geom, cond, None) would choose for the state
    in cond, WITHOUT computing the pressure drop.

    Purpose: the model (hex_MB) calls it once per side and per phase, at a representative state of
    the phase, and then gives this name explicitly to every cell of the phase. The automatic rules
    (select_* functions of each dispatcher) are not duplicated: they are called here.

    Arguments
    ---------
    AS   : CoolProp AbstractState set at the state of cond (h, p)
    geom : geometry dictionary of the side (as for pressure_drop)
    cond : condition dictionary (as for pressure_drop): flow_state 'SP', 'EVAP' or 'COND' and phase
           ('vapor-wet' cells use the single-phase pressure drop, as in pressure_drop)

    Returns
    -------
    str : name of the correlation, valid as an explicit choice for pressure_drop.
    """
    # Selection functions of the dispatchers (imported here to keep this helper self-contained)
    from labothappy.correlations.pressure_drop.BPHX_DP import (
        select_plate_DP, select_plate_evaporation_DP, select_plate_condensation_DP)
    from labothappy.correlations.pressure_drop.STHE_DP import select_shell_DP
    from labothappy.correlations.pressure_drop.pipe_DP import (
        select_pipe_DP, select_pipe_two_phase_DP)
    from labothappy.correlations.pressure_drop.PCHE_DP import select_PCHE_DP

    htx = geom.get('HTX_Type')
    # Same flow state as pressure_drop: a 'vapor-wet' cell has a single-phase pressure drop
    if cond.get('phase') == 'vapor-wet':
        state = 'SP'
    else:
        state = cond.get('flow_state')
    if state not in _VALID_STATES:
        raise ValueError(f"Pressure drop: flow state {state!r}, expected one of {_VALID_STATES}.")
    side = cond.get('side', geom.get('side'))

    if htx == 'BPHE':
        if state == 'SP':
            return select_plate_DP(geom, cond)
        if state == 'EVAP':
            return select_plate_evaporation_DP(geom, cond)
        return select_plate_condensation_DP(geom, cond)

    if htx == 'STHE':
        if side == 'tube':
            if state == 'SP':
                return select_pipe_DP(AS, geom, cond)
            return select_pipe_two_phase_DP(AS, geom, cond)
        if state == 'SP':
            return select_shell_DP(geom, cond)
        return 'Ishihara'                               # STHE_two_phase_DP automatic choice

    if htx == 'TFHE':
        if side == 'tube':
            if state == 'SP':
                return select_pipe_DP(AS, geom, cond)
            return select_pipe_two_phase_DP(AS, geom, cond)
        return 'Pismennyi'                              # TFHE_fin_DP automatic choice

    if htx == 'PCHE':
        if state == 'SP':
            name = select_PCHE_DP(geom, cond)
            if name is None:                            # straight channel: pipe rule
                name = select_pipe_DP(AS, geom, cond)
            return name
        return select_pipe_two_phase_DP(AS, geom, cond)

    if htx == 'SPHE':
        if state == 'SP':
            return 'Picon_Nunez'                        # SPHE_DP automatic choice
        return select_pipe_two_phase_DP(AS, geom, cond)

    if htx in _CURVED_PIPE_TYPES:
        if state == 'SP':
            return 'Sigalotti'                          # curved_pipe_DP automatic choice
        return select_pipe_two_phase_DP(AS, geom, cond)

    raise NotImplementedError(f"auto_dp_correlation: heat exchanger type {htx!r} is not handled.")


def dp_candidates(AS, geom, cond):
    """Names of every registered pressure drop correlation usable for the case of cond
    (heat exchanger type, side, flow state; a 'vapor-wet' cell is single phase), in priority
    order: the rule-based automatic choice (auto_dp_correlation) first, then the others in the
    order of the registry tables. A function registered under several names appears once."""
    from labothappy.correlations.pressure_drop import BPHX_DP, STHE_DP, pipe_DP, PCHE_DP, TFHE_DP
    from labothappy.correlations.pressure_drop import curved_pipe_DP, SPHE_DP

    htx = geom.get('HTX_Type')
    if cond.get('phase') == 'vapor-wet':
        state = 'SP'
    else:
        state = cond.get('flow_state')
    side = cond.get('side', geom.get('side'))
    pipe_tables = {'SP': pipe_DP._PIPE_SP_DP, 'EVAP': pipe_DP._PIPE_EVAP_DP,
                   'COND': pipe_DP._PIPE_COND_DP}

    # Registry tables of the case
    if htx == 'BPHE':
        tables = [{'SP': BPHX_DP._BPHX_SP_DP, 'EVAP': BPHX_DP._BPHX_EVAP_DP,
                   'COND': BPHX_DP._BPHX_COND_DP}[state]]
    elif htx in ('STHE', 'TFHE') and side == 'tube':
        tables = [pipe_tables[state]]
    elif htx == 'STHE' and state == 'SP':
        tables = [STHE_DP._STHE_SP_DP]
    elif htx == 'STHE':
        tables = [STHE_DP._STHE_2P_DP]
    elif htx == 'TFHE':
        tables = [TFHE_DP._TFHE_FIN_DP]
    elif htx == 'PCHE' and state == 'SP':
        tables = [PCHE_DP._PCHE_SP_DP, pipe_tables['SP']]
    elif htx == 'PCHE':
        tables = [PCHE_DP._PCHE_2P_DP, pipe_tables[state]]
    elif htx == 'SPHE' and state == 'SP':
        tables = [SPHE_DP._SPHE_SP, curved_pipe_DP._CURVED['SP']]
    elif htx in ('SPHE',) + _CURVED_PIPE_TYPES and state != 'SP':
        tables = [curved_pipe_DP._CURVED[state], pipe_tables[state]]
    elif htx in _CURVED_PIPE_TYPES:
        tables = [curved_pipe_DP._CURVED['SP']]
    else:
        raise NotImplementedError(f"dp_candidates: heat exchanger type {htx!r} is not handled.")

    # Rule-based choice first, then the registry order, each function once
    names = [auto_dp_correlation(AS, geom, cond)]
    functions = []
    for table in tables:
        for name, function in table.items():
            if name in names or function in functions:
                continue
            names.append(name)
            functions.append(function)
    return names


def choose_dp_correlation(AS, geom, cond):
    """Automatic choice of the pressure drop correlation for the state of cond, by TRYING every
    registered candidate (dp_candidates) at this state: first the correlations whose parameters
    are all given, then the ones inside their validity range, then the ones established for the
    fluid (relaxed in the reverse order; see correlation_warnings.choose_by_trial).

    Arguments: as pressure_drop (AS set at the state of cond).
    Returns: name of the correlation (a valid explicit choice for pressure_drop).
    Raises ValueError if no candidate can be computed at this state."""
    import CoolProp.CoolProp as CP
    from labothappy.correlations.convection.correlation_warnings import choose_by_trial

    # The correlations may move AS (saturation properties): it is set back before every trial
    h, p = AS.hmass(), AS.p()

    def evaluate(name):
        AS.update(CP.HmassP_INPUTS, h, p)
        return pressure_drop(AS, geom, cond, name)

    side = cond.get('side', geom.get('side', ''))
    label = f"{geom.get('HTX_Type')} {side} {cond.get('flow_state')} DP"
    names = dp_candidates(AS, geom, cond)
    chosen, rejected = choose_by_trial(names, evaluate, label)
    AS.update(CP.HmassP_INPUTS, h, p)
    if chosen is None:
        raise ValueError(f"{label}: no candidate correlation can be computed ({rejected}).")
    return chosen
