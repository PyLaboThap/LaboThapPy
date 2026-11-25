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


# =============================================================================
# PARAMETERS NEEDED BY EACH CORRELATION IN ADDITION TO THE GEOMETRY
# =============================================================================

# Model parameters that each registered correlation needs IN ADDITION to the
# parameters of the heat exchanger geometry (hex_geometry required_parameters): keys
# that the geometry does not compute. () = the geometry parameters are enough.
# Found on 2026-11-25 by removing the keys one at a time at hex_MB states
# (correlation_warnings.parameters_by_removal) and by reading the code.
_REQUIRED_PARAMETERS = {
    'Picon_Nunez': (),
    'Curved_rect_duct': (),
    'Rect_duct': (),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"SPHE_DP: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
