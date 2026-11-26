# -*- coding: utf-8 -*-
"""
Tube bank pressure drop dispatcher (outer side of banks of bare tubes in crossflow): flow state ->
correlation, (AS, geom, cond) formalism. Created 2026-11-26 (replaces the pressure drop functions of
the old convection/tube_bank_htc.py).

  'SP'            : 'Zukauskas' (= None / 'auto')  tube_bank_zukauskas_SP_DP (whole bank, all rows)
  'EVAP' / 'COND' : 'Ishihara' (= None / 'auto')   tube_bank_ishihara_two_phase_DP (Ishihara 1980
                    multiplier on the Zukauskas single-phase DP)

Geometry and conditions: see bare_tube_bank_DP_correlation_bank (TFHE geometry names).
"""
from labothappy.correlations.pressure_drop.bare_tube_bank_DP_correlation_bank import (
    tube_bank_zukauskas_SP_DP, tube_bank_ishihara_two_phase_DP)

_TUBE_BANK_SP_DP = {
    'Zukauskas': tube_bank_zukauskas_SP_DP,
}
_TUBE_BANK_2P_DP = {
    'Ishihara': tube_bank_ishihara_two_phase_DP,
}
_TABLES = {'SP': _TUBE_BANK_SP_DP, 'EVAP': _TUBE_BANK_2P_DP, 'COND': _TUBE_BANK_2P_DP}
_AUTO = {'SP': 'Zukauskas', 'EVAP': 'Ishihara', 'COND': 'Ishihara'}


def select_bare_tube_bank_DP(geom, cond, flow_state='SP'):
    """Automatic choice: SP -> 'Zukauskas', EVAP / COND -> 'Ishihara'."""
    return _AUTO[flow_state]


def bare_tube_bank_DP(AS, geom, cond, flow_state, choices=None):
    """Pressure drop of the outer side of a bare tube bank [Pa] (whole bank).
    flow_state : 'SP', 'EVAP' or 'COND'.
    choices    : dict {'SP': name, 'EVAP': name, 'COND': name}, a name, or None / 'auto'."""
    if flow_state not in _TABLES:
        raise ValueError(f"bare_tube_bank_DP: flow state {flow_state!r}, expected 'SP', "
                         "'EVAP' or 'COND'.")
    table = _TABLES[flow_state]
    if isinstance(choices, dict):
        choice = choices.get(flow_state)
    else:
        choice = choices
    if choice in (None, 'auto'):
        choice = select_bare_tube_bank_DP(geom, cond, flow_state)
    if choice not in table:
        raise ValueError(f"Unknown tube bank {flow_state} pressure drop correlation {choice!r}: "
                         f"{sorted(table)} or None/'auto'.")
    return table[choice](AS, geom, cond)


# Model parameters that each registered correlation needs IN ADDITION to the geometry. () = none.
_REQUIRED_PARAMETERS = {
    'Zukauskas': (),
    'Ishihara': (),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"bare_tube_bank_DP: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
