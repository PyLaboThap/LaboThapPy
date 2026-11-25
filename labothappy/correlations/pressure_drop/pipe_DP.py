# -*- coding: utf-8 -*-
"""
In-pipe (tube-side) pressure drop dispatcher: flow state (SP / EVAP / COND) -> correlation.
"""

import warnings

from labothappy.correlations.pressure_drop.pipe_DP_correlation_bank import (
    pipe_churchill_SP_DP, pipe_swamee_jain_SP_DP, pipe_haaland_SP_DP, pipe_konakov_SP_DP, pipe_petukhov_SP_DP,
    pipe_cheng_CO2_SP_DP,
    pipe_friedel_boiling_DP, pipe_msh_boiling_DP, pipe_choi_boiling_DP, pipe_homogeneous_boiling_DP,
    pipe_lockhart_martinelli_boiling_DP,
    pipe_friedel_condensing_DP, pipe_msh_condensing_DP, pipe_choi_condensing_DP, pipe_homogeneous_condensing_DP,
    pipe_lockhart_martinelli_condensing_DP,
    pipe_xu_fang_condensing_DP, pipe_chisholm_DP, pipe_gronnerud_boiling_DP,
    pipe_kim_mudawar_condensing_DP)
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props

_PIPE_SP_DP = {
    'Churchill': pipe_churchill_SP_DP,
    'Swamee_Jain': pipe_swamee_jain_SP_DP,
    'Haaland': pipe_haaland_SP_DP,
    'Konakov': pipe_konakov_SP_DP,
    'Petukhov': pipe_petukhov_SP_DP,
    'Cheng_CO2': pipe_cheng_CO2_SP_DP,
}
_PIPE_EVAP_DP = {
    'Friedel': pipe_friedel_boiling_DP,
    'MSH': pipe_msh_boiling_DP,
    'Choi': pipe_choi_boiling_DP,
    'Lockhart_Martinelli': pipe_lockhart_martinelli_boiling_DP,
    'Homogeneous': pipe_homogeneous_boiling_DP,
    'Chisholm': pipe_chisholm_DP,
    'Gronnerud': pipe_gronnerud_boiling_DP,
}
_PIPE_COND_DP = {
    'Friedel': pipe_friedel_condensing_DP,
    'MSH': pipe_msh_condensing_DP,
    'Choi': pipe_choi_condensing_DP,
    'Lockhart_Martinelli': pipe_lockhart_martinelli_condensing_DP,
    'Homogeneous': pipe_homogeneous_condensing_DP,
    'Xu_Fang': pipe_xu_fang_condensing_DP,
    'Chisholm': pipe_chisholm_DP,
    'Kim_Mudawar': pipe_kim_mudawar_condensing_DP,    # mini / micro-channels (D 0.07-6.2 mm)
}

def _cheng_available():
    try:
        from labothappy.correlations.properties.supercritical import (  # noqa: F401
            pseudo_boiling_bounds, pseudo_critical_temperature_CO2)
        return True
    except ImportError:
        return False

def select_pipe_DP(AS, geom, cond):
    """Single-phase friction (heuristic):
      - supercritical CO2 inside the Cheng window (D 8-12 mm, G 497-1346 kg/m2/s, p 7.5-23.5 MPa; needs cond h)
        -> Cheng_CO2 (same window as the Cheng HTC)
      - otherwise Churchill (laminar to turbulent, roughness).
    Swamee_Jain, Haaland (turbulent), Konakov, Petukhov (turbulent, smooth) are explicit choices only."""
    p = cond.get('p')
    if p is not None and 'h' in cond and AS.fluid_names()[0] in ('CarbonDioxide', 'CO2') and p > AS.p_critical():
        D = geom['Dh'] if 'Dh' in geom else geom['D_t_o'] - 2 * geom['t_t']
        G = cond.get('G', 0.0)
        if 8e-3 <= D <= 12e-3 and 496.7 <= G <= 1346.2 and 7.5e6 <= p <= 23.51e6 and _cheng_available():
            return 'Cheng_CO2'
    return 'Churchill'


def select_pipe_two_phase_DP(AS, geom, cond):
    """Two-phase friction (heuristic):
      - mu_l / mu_v < 1000 and a surface tension available -> Friedel (the range it is recommended for)
      - otherwise -> MSH (Muller-Steinhagen & Heck: no surface tension, wide database).
    Choi (needs the quality change of the pipe), Lockhart_Martinelli and Homogeneous are explicit choices."""
    s = _sat_props(AS, cond['p'], cond)
    if s['mu_l'] / s['mu_v'] < 1000 and s['sigma'] is not None:
        return 'Friedel'
    return 'MSH'


def _dispatch(table, select, AS, geom, cond, choice, label):
    if choice in (None, 'auto'):
        choice = select(AS, geom, cond)
    elif choice not in table:
        raise ValueError(f"Unknown {label} correlation {choice!r}: {sorted(table)} or None/'auto'.")
    return table[choice](AS, geom, cond)


def pipe_SP_DP(AS, geom, cond, choice=None):
    return _dispatch(_PIPE_SP_DP, select_pipe_DP, AS, geom, cond, choice, 'single-phase pressure drop')


def pipe_boiling_DP(AS, geom, cond, choice=None):
    return _dispatch(_PIPE_EVAP_DP, select_pipe_two_phase_DP, AS, geom, cond, choice, 'evaporation pressure drop')


def pipe_condensing_DP(AS, geom, cond, choice=None):
    return _dispatch(_PIPE_COND_DP, select_pipe_two_phase_DP, AS, geom, cond, choice, 'condensation pressure drop')


_STATE_DISPATCH = {'SP': pipe_SP_DP, 'EVAP': pipe_boiling_DP, 'COND': pipe_condensing_DP}


def pipe_DP(AS, geom, cond, flow_state, choices=None):
    """Frictional pressure drop of a side inside tubes, over the full flow length geom['L'] [Pa].
    flow_state : 'SP', 'EVAP' or 'COND'.   choices : dict {'SP': name, 'EVAP': name, 'COND': name} or None / 'auto'."""
    if flow_state not in _STATE_DISPATCH:
        raise ValueError(f"Unknown flow state {flow_state!r}: expected one of {sorted(_STATE_DISPATCH)}.")
    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Churchill'}}, got {type(choices).__name__}.")
    unknown = set(choices) - set(_STATE_DISPATCH)
    if unknown:
        raise ValueError(f"Unknown keys in choices: {sorted(unknown)} (expected {sorted(_STATE_DISPATCH)}).")
    choice = choices.get(flow_state)
    if not choice or (isinstance(choice, str) and choice.lower() == 'auto'):
        choice = None
    return _STATE_DISPATCH[flow_state](AS, geom, cond, choice)


# =============================================================================
# PARAMETERS NEEDED BY EACH CORRELATION IN ADDITION TO THE GEOMETRY
# =============================================================================

# Model parameters that each registered correlation needs IN ADDITION to the
# parameters of the heat exchanger geometry (hex_geometry required_parameters): keys
# that the geometry does not compute. () = the geometry parameters are enough.
# Found on 2026-11-25 by removing the keys one at a time at hex_MB states
# (correlation_warnings.parameters_by_removal) and by reading the code.
_REQUIRED_PARAMETERS = {
    'Churchill': (),
    'Swamee_Jain': (),
    'Haaland': (),
    'Konakov': (),
    'Petukhov': (),
    'Cheng_CO2': (),
    'Friedel': (),
    'MSH': (),
    'Choi': (),
    'Lockhart_Martinelli': (),
    'Homogeneous': (),
    'Chisholm': (),
    'Gronnerud': (),
    'Xu_Fang': (),
    'Kim_Mudawar': (),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"pipe_DP: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
