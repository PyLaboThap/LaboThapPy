# -*- coding: utf-8 -*-
"""
Tube bank HTC dispatcher (outer side of banks of bare tubes in crossflow):
flow state -> correlation, (AS, geom, cond) formalism.
Created 2026-11-26 (replaces the old convection/tube_bank_htc.py).

  'SP'   : 'Zukauskas' (= None / 'auto')            tube_bank_zukauskas_SP_HTC
           'ESDU_73031'                             tube_bank_esdu_73031_SP_HTC
  'EVAP' : 'Shah' (= None / 'auto')                 tube_bank_shah_boiling_HTC (Shah 2017, bundles)
  'COND' : 'Rose_Kern' (= None / 'auto')            tube_bank_rose_kern_condensing_HTC
                                                    (single tube x Kern inundation, N = rows
                                                    crossed by the flow)
           'Nusselt_single_tube'                    tube_nusselt_film_condensing_HTC (one tube)
The critical heat flux of crossflow boiling is available as
bare_tube_bank_HTC_correlation_bank.tube_lienhard_eichhorn_CHF (not an HTC, not registered).

Geometry and conditions: see bare_tube_bank_HTC_correlation_bank (TFHE geometry names).
"""
from labothappy.correlations.convection.bare_tube_bank_HTC_correlation_bank import (
    tube_bank_zukauskas_SP_HTC, tube_nusselt_film_condensing_HTC, tube_bank_esdu_73031_SP_HTC)
from labothappy.correlations.convection.bare_tube_bank_HTC_correlation_bank import (
    tube_bank_shah_boiling_HTC, tube_bank_rose_kern_condensing_HTC)

_TUBE_BANK_SP_HTC = {
    'Zukauskas': tube_bank_zukauskas_SP_HTC,
    'ESDU_73031': tube_bank_esdu_73031_SP_HTC,
}
_TUBE_BANK_EVAP_HTC = {
    'Shah': tube_bank_shah_boiling_HTC,
}
_TUBE_BANK_COND_HTC = {
    'Rose_Kern': tube_bank_rose_kern_condensing_HTC,
    'Nusselt_single_tube': tube_nusselt_film_condensing_HTC,
}
_TABLES = {'SP': _TUBE_BANK_SP_HTC, 'EVAP': _TUBE_BANK_EVAP_HTC, 'COND': _TUBE_BANK_COND_HTC}
_AUTO = {'SP': 'Zukauskas', 'EVAP': 'Shah', 'COND': 'Rose_Kern'}


def select_bare_tube_bank_HTC(geom, cond, flow_state):
    """Automatic choice: SP -> 'Zukauskas' ('ESDU_73031' is an explicit choice), EVAP -> 'Shah',
    COND -> 'Rose_Kern' ('Nusselt_single_tube' is an explicit choice)."""
    return _AUTO[flow_state]


def bare_tube_bank_HTC(AS, geom, cond, flow_state, choices=None):
    """Heat transfer coefficient of the outer side of a bare tube bank [W/m2/K].
    flow_state : 'SP', 'EVAP' or 'COND'.
    choices    : dict {'SP': name, 'EVAP': name, 'COND': name}, a name, or None / 'auto'."""
    if flow_state not in _TABLES:
        raise ValueError(f"bare_tube_bank_HTC: flow state {flow_state!r}, expected 'SP', "
                         "'EVAP' or 'COND'.")
    table = _TABLES[flow_state]
    if isinstance(choices, dict):
        choice = choices.get(flow_state)
    else:
        choice = choices
    if choice in (None, 'auto'):
        choice = select_bare_tube_bank_HTC(geom, cond, flow_state)
    if choice not in table:
        raise ValueError(f"Unknown tube bank {flow_state} correlation {choice!r}: {sorted(table)} "
                         "or None/'auto'.")
    return table[choice](AS, geom, cond)


# Model parameters that each registered correlation needs IN ADDITION to the geometry
# (D_t_o, L_t, n_t, n_rows, n_passes, P_V, P_H, theta). () = none.
_REQUIRED_PARAMETERS = {
    'Zukauskas': (),
    'ESDU_73031': (),
    'Shah': (),
    'Rose_Kern': (),
    'Nusselt_single_tube': (),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"bare_tube_bank_HTC: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
