# Adapt the import paths to where the dispatchers live.
from labothappy.correlations.pressure_drop.BPHX_DP import BPHE_DP, BPHX_port_DP
from labothappy.correlations.pressure_drop.STHE_DP import STHE_DP
from labothappy.correlations.pressure_drop.pipe_DP import pipe_DP
from labothappy.correlations.pressure_drop.TFHE_DP import TFHE_DP

_VALID_STATES = ('SP', 'EVAP', 'COND')
_VALID_SIDES = ('shell', 'tube')


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
    STHE : cond['side'] (or geom['side']) = 'shell' (single phase only: STHE_DP) or 'tube' (pipe_DP).
    TFHE ('Tube&Fins') : cond['side'] (or geom['side']) = 'fin' (gas side: Pismennyi bank pressure drop,
           TFHE_DP) or 'tube' (pipe_DP).
    PCHE : pipe_DP on the equivalent tube (geom D_h / Tube_ID, L, K; cond G of ONE channel); no port losses."""
    htx = geom.get('HTX_Type')

    if htx in ('BPHE', 'STHE', 'PCHE', 'TFHE', 'Tube&Fins'):
        state = 'SP' if cond.get('phase') == 'vapor-wet' else cond.get('flow_state')
        if state not in _VALID_STATES:
            raise ValueError(f"Pressure drop: flow state {state!r}, expected one of {_VALID_STATES}.")
        choices = correlation if isinstance(correlation, dict) else {state: correlation}

        if htx == 'BPHE':
            return float(BPHE_DP(AS, geom, cond, state, choices))

        if htx == 'PCHE':      # same as the pipe pressure drop, on the equivalent tube (D_h, L_c * n_series, K, G of one channel)
            return float(pipe_DP(AS, geom, cond, state, choices))

        if htx in ('TFHE', 'Tube&Fins'):      # cond['side'] / geom['side'] = 'fin' or 'tube'
            return float(TFHE_DP(AS, geom, cond, state, choices))

        side = cond.get('side', geom.get('side'))
        if side == 'shell':
            return float(STHE_DP(AS, geom, cond, state, choices))
        if side == 'tube':
            return float(pipe_DP(AS, geom, cond, state, choices))
        raise ValueError(f"STHE: cond['side'] (or geom['side']) = {side!r}: expected one of {_VALID_SIDES}.")

    raise NotImplementedError(f"pressure_drop: heat exchanger type {htx!r} is not handled yet by the "
                              "three-level scheme ('BPHE', 'STHE', 'PCHE' and 'Tube&Fins').")


def port_pressure_drop(geom, cond):
    """Port / manifold pressure loss [Pa] of one side, for the WHOLE exchanger (not per cell).
    cond must hold m_dot and D (density, evaluated at the inlet state), or G_port; geom may hold
    D_port and n_passes. Returns 0 (with a warning) when no port size is available.
    STHE, PCHE, Tube&Fins: 0 (nozzle, header and head losses are not modelled)."""
    if geom.get('HTX_Type') == 'BPHE':
        return float(BPHX_port_DP(geom, cond))
    return 0.0

