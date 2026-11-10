# -*- coding: utf-8 -*-
"""
SPHE (spiral plate) HTC dispatcher, (AS, geom, cond) formalism.

SP   : spiral plate correlations only: 'Picon_Nunez' (default), 'Rajavel' (water), 'Curved_rect_duct', 'Rect_duct'.
EVAP : helical pipe correlations on Dh / R_curv: 'Zhao' (default), 'Nariai'.      (needs cond q)
COND : helical pipe correlation: 'Zakeralhoseini' (R1234yf).
"""
import warnings
from labothappy.correlations.convection.SPHE_HTC_correlation_bank import (
    SPHE_picon_nunez_HTC, SPHE_rajavel_HTC, SPHE_curved_rect_duct_HTC, SPHE_rect_duct_HTC)
from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import (
    curved_pipe_zhao_boiling_HTC, curved_pipe_nariai_boiling_HTC, curved_pipe_zakeralhoseini_condensing_HTC)

_SPHE = {
    'SP':   {'Picon_Nunez': SPHE_picon_nunez_HTC, 'Rajavel': SPHE_rajavel_HTC,
             'Curved_rect_duct': SPHE_curved_rect_duct_HTC, 'Rect_duct': SPHE_rect_duct_HTC},
    'EVAP': {'Zhao': curved_pipe_zhao_boiling_HTC, 'Nariai': curved_pipe_nariai_boiling_HTC},
    'COND': {'Zakeralhoseini': curved_pipe_zakeralhoseini_condensing_HTC},
}
_AUTO = {'SP': 'Picon_Nunez', 'EVAP': 'Zhao', 'COND': 'Zakeralhoseini'}


def SPHE_HTC(AS, geom, cond, flow_state, choices=None):
    if flow_state not in _SPHE:
        raise ValueError(f"SPHE HTC: flow_state {flow_state!r}: expected 'SP', 'EVAP' or 'COND'.")
    choice = (choices or {}).get(flow_state)
    if choice in (None, 'auto'):
        choice = _AUTO[flow_state]
    if choice not in _SPHE[flow_state]:
        raise ValueError(f"SPHE HTC: unknown {flow_state} correlation {choice!r}; choose from {list(_SPHE[flow_state])}.")
    if choice == 'Zakeralhoseini':
        warnings.warn("SPHE: Zakeralhoseini condensation correlation established for R1234yf in helical tubes.")
    return _SPHE[flow_state][choice](AS, geom, cond)

