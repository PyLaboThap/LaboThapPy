# -*- coding: utf-8 -*-
"""
Curved pipe (helical coil, spiral channel) pressure drop dispatcher: flow state (SP / EVAP / COND) -> correlation.

SP   : 'Sigalotti' (default, your curved_pipe_1P_DP), 'Ito', or a pipe name (straight equivalent tube).
EVAP / COND : the pipe two-phase correlations of pipe_DP ('Friedel', 'MSH', 'Choi', 'Lockhart_Martinelli',
       'Homogeneous'; None / 'auto': Friedel or MSH as chosen by the pipe dispatcher) on the straight equivalent
       tube D = D_h, without curvature effect (no warning: it is the default).
       The curved two-phase correlations stay available by name: 'Cioncolini', 'Ruffell', 'Fsadni' (COND; your
       curved_pipe_condensing_DP).
"""
import warnings
from labothappy.correlations.pressure_drop.curved_pipe_DP_correlation_bank import (
    curved_pipe_ito_SP_DP, curved_pipe_sigalotti_SP_DP, curved_pipe_cioncolini_boiling_DP, curved_pipe_ruffell_boiling_DP,
    curved_pipe_fsadni_condensing_DP, D_over_2R)
from labothappy.correlations.pressure_drop.pipe_DP import pipe_DP
from labothappy.correlations.convection.correlation_warnings import ValidityRangeWarning

_CURVED = {
    'SP': {'Sigalotti': curved_pipe_sigalotti_SP_DP, 'Ito': curved_pipe_ito_SP_DP},
    'EVAP': {'Cioncolini': curved_pipe_cioncolini_boiling_DP, 'Ruffell': curved_pipe_ruffell_boiling_DP},
    'COND': {'Fsadni': curved_pipe_fsadni_condensing_DP, 'Cioncolini': curved_pipe_cioncolini_boiling_DP,
             'Ruffell': curved_pipe_ruffell_boiling_DP},
}


def _equivalent_tube(geom):
    g = dict(geom)
    D = g.get('Dh') or g.get('D_t_i')
    g.update(D_t_i=D, Dh=D, n_passes=1)
    g.setdefault('L_t', g.get('L'))
    g.setdefault('n_t', 1)
    return g


def curved_pipe_DP(AS, geom, cond, flow_state, choices=None):
    if flow_state not in _CURVED:
        raise ValueError(f"Curved pipe DP: flow_state {flow_state!r}: expected 'SP', 'EVAP' or 'COND'.")
    choice = (choices or {}).get(flow_state)
    table = _CURVED[flow_state]
    if choice in (None, 'auto') and flow_state == 'SP':
        choice = 'Sigalotti'
    elif choice in (None, 'auto'):
        return pipe_DP(AS, _equivalent_tube(geom), cond, flow_state, {})      # pipe two-phase, auto
    if choice in table:
        return table[choice](AS, geom, cond)
    if flow_state == 'SP':
        warnings.warn(f"Curved pipe DP: pipe correlation {choice!r} on the straight equivalent tube (no curvature).",
                      ValidityRangeWarning)
    return pipe_DP(AS, _equivalent_tube(geom), cond, flow_state, {flow_state: choice})
