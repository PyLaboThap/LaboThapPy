# -*- coding: utf-8 -*-
"""
SPHE (spiral plate) pressure drop dispatcher, (AS, geom, cond) formalism.

SP   : 'Picon_Nunez' (default, your single_phase_spiral_plate_DP: Minton), 'Curved_rect_duct', 'Rect_duct' (straight rectangular duct), 'Sigalotti' / 'Ito' (circular curved pipe on Dh),
       or any pipe correlation name on the equivalent tube.
EVAP / COND : pipe two-phase correlations on D_h (default: auto Friedel / MSH; no curvature effect);
       'Cioncolini', 'Ruffell', 'Fsadni' by name (curved_pipe_DP).
"""
from labothappy.correlations.pressure_drop.SPHE_DP_correlation_bank import (SPHE_curved_rect_SP_DP, SPHE_rect_SP_DP,
                                                                          SPHE_picon_nunez_SP_DP)
from labothappy.correlations.pressure_drop.curved_pipe_DP import curved_pipe_DP

_SPHE_SP = {'Picon_Nunez': SPHE_picon_nunez_SP_DP, 'Curved_rect_duct': SPHE_curved_rect_SP_DP, 'Rect_duct': SPHE_rect_SP_DP}

def SPHE_DP(AS, geom, cond, flow_state, choices=None):
    choice = (choices or {}).get(flow_state)
    if flow_state == 'SP':
        if choice in (None, 'auto'):
            choice = 'Picon_Nunez'
        if choice in _SPHE_SP:
            return _SPHE_SP[choice](AS, geom, cond)
    return curved_pipe_DP(AS, geom, cond, flow_state, choices)       # 'Ito', pipe names, two phase

