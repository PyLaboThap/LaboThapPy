# -*- coding: utf-8 -*-
"""
Tube & fins (TFHE) pressure drop dispatcher: side ('fin' / 'tube') -> flow state -> correlation,
(AS, geom, cond) formalism.

side 'fin'  : finned side (gas): 'Pismennyi' (also 'Pismennyi_Annular' / 'Pismennyi_Square': same formula). Two-phase / 'vapor-wet' cells (EVAP / COND) reuse it with
              the bulk properties of the cell (homogeneous mixture for a two-phase cell), phase change ignored, warning.
side 'tube' : SP / EVAP / COND inside the tubes: pipe correlations (pipe_DP).
Pressure drops are those of the FULL flow length of the side (the model scales them to the cells).
"""
import warnings
from labothappy.correlations.pressure_drop.TFHE_DP_correlation_bank import TFHE_pismennyi_DP
from labothappy.correlations.pressure_drop.pipe_DP import pipe_DP
from labothappy.correlations.convection.TFHE_HTC import _ensure_gas_props

# same formula for both fin shapes: the names of the HTC dispatcher are accepted so that one set of names serves both
_TFHE_FIN_DP = {'Pismennyi': TFHE_pismennyi_DP,
                'Pismennyi_Annular': TFHE_pismennyi_DP,
                'Pismennyi_Square': TFHE_pismennyi_DP}


def TFHE_fin_DP(AS, geom, cond, flow_state, choices=None):
    if flow_state not in ('SP', 'EVAP', 'COND'):
        raise ValueError(f"TFHE finned side: flow state {flow_state!r}, expected 'SP', 'EVAP' or 'COND'.")
    if flow_state != 'SP':
        warnings.warn(f"TFHE finned side: no {flow_state} pressure drop correlation, the gas correlation is used "
                      f"for a {cond.get('phase')!r} cell (phase change ignored).", stacklevel=2)
    cond = _ensure_gas_props(AS, cond, flow_state)
    choice = (choices or {}).get(flow_state) or (choices or {}).get('SP')
    if choice in (None, 'auto'):
        choice = 'Pismennyi'
    if choice not in _TFHE_FIN_DP:
        raise ValueError(f"Unknown finned-side pressure drop correlation {choice!r}: {sorted(_TFHE_FIN_DP)} or None/'auto'.")
    return _TFHE_FIN_DP[choice](AS, geom, cond)


def TFHE_DP(AS, geom, cond, flow_state, choices=None):
    side = cond.get('side', geom.get('side'))
    if side == 'fin':
        return TFHE_fin_DP(AS, geom, cond, flow_state, choices)
    if side == 'tube':
        return pipe_DP(AS, geom, cond, flow_state, choices)
    raise ValueError(f"TFHE: cond['side'] (or geom['side']) = {side!r}: expected 'fin' or 'tube'.")
    
    