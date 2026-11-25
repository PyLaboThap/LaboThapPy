# -*- coding: utf-8 -*-
"""
Shell-side HTC dispatcher of shell & tube heat exchangers, same style as BPHX_HTC_dispatcher:
HX type -> flow state (SP / EVAP / COND) -> correlation, with the (AS, geom, cond) formalism.
Implemented: SP (Kern, Donohue, Bell-Delaware) and EVAP (Cooper, Mostinski pool boiling, no
tube-bundle factor).
Added 2026-11-23: EVAP 'Shah' (Shah 2017, tube bundle, automatic when the Kern geometry is known)
and COND 'Rose_Kern' (Bonneau et al. 2019 recommendation).
"""

#%% SINGLE PHASE : SHELL-SIDE HTC

from labothappy.correlations.convection.STHE_HTC_correlation_bank import (
    STHE_kern_shell_HTC, STHE_donohue_shell_HTC, STHE_bell_delaware_shell_HTC)

_STHE_SP_HTC = {
    'Kern': STHE_kern_shell_HTC,
    'Donohue': STHE_donohue_shell_HTC,
    'Bell_Delaware': STHE_bell_delaware_shell_HTC,
}

# Geometry keys needed by each correlation (on top of the common ones: m_dot, mu, Pr, k/cp in cond)
_KERN_KEYS = ('D_t_o', 'PR', 'D_s', 'B', 'theta')
_DONOHUE_KEYS = _KERN_KEYS + ('BC', 'n_t')
_BELL_KEYS = _DONOHUE_KEYS + ('D_OTL', 'delta_TB', 'delta_BS', 'N_strips',
                              'L_t', 'TS_t', 'B_i', 'B_o')
_REQUIRED = {'Kern': _KERN_KEYS, 'Donohue': _DONOHUE_KEYS, 'Bell_Delaware': _BELL_KEYS}


def _missing(geom, name):
    return [k for k in _REQUIRED[name] if k not in geom]

def select_shell_HTC(geom, cond):
    """Automatic choice of the shell-side HTC correlation (heuristic): the most complete
    correlation the available geometry allows.
      - full baffle / clearance / end-space data -> Bell_Delaware (leakage, bypass, window,
                                                    unequal end-spacing corrections)
      - baffle cut and tube count known          -> Donohue (window flow included, no leakage terms)
      - otherwise                                -> Kern (ideal crossflow only, needs tube bank data)
    """
    if not _missing(geom, 'Bell_Delaware'):
        return 'Bell_Delaware'
    if not _missing(geom, 'Donohue'):
        return 'Donohue'
    return 'Kern'

def STHE_SP_HTC(AS, geom, cond, choice=None):
    """Single-phase shell-side heat transfer coefficient [W/m2/K].
    choice: None/'auto' or the name of a correlation.

    cond : m_dot (flow of the whole shell side) [kg/s], mu, Pr, k, cp; mu_w optional (mu_w = mu if absent)
    geom : see _REQUIRED; n_parallel (optional) divides the flow rate."""
    if choice in (None, 'auto'):
        choice = select_shell_HTC(geom, cond)
    elif choice not in _STHE_SP_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_STHE_SP_HTC)} or None/'auto'.")
    missing = _missing(geom, choice)
    if missing:
        raise KeyError(f"{choice} shell HTC: missing geometry keys {missing}.")
    return _STHE_SP_HTC[choice](AS, geom, cond)

#%% TWO-PHASE : SHELL-SIDE BOILING (pool boiling on the tubes, no tube-bundle factor)

from labothappy.correlations.convection.STHE_HTC_correlation_bank import STHE_mostinski_boiling_HTC
from labothappy.correlations.convection.BPHE_HTC_correlation_bank import BPHE_cooper_boiling_HTC
from labothappy.correlations.convection.STHE_HTC_correlation_bank import (
    STHE_shah_bundle_boiling_HTC, STHE_rose_kern_condensing_HTC)

_STHE_EVAP_HTC = {
    'Cooper': BPHE_cooper_boiling_HTC,          # Cooper (1984) pool boiling, needs p, q (+ roughness)
    'Mostinski': STHE_mostinski_boiling_HTC,    # Mostinski (1963) reduced-pressure method, p, q
    'Shah': STHE_shah_bundle_boiling_HTC,       # Shah (2017) tube bundle, m_dot, p, x, q + Kern geometry
}

def select_shell_evaporation_HTC(geom, cond):
    """Automatic choice of the shell-side boiling HTC: Cooper (1984), built on a wide multi-fluid
    pool boiling database. Both correlations are single-tube pool boiling: no tube-bundle factor.
    Shah (2017) tube-bundle correlation when the Kern geometry and m_dot, x are known."""
    if 'm_dot' in cond and 'x' in cond and not _missing(geom, 'Kern'):
        return 'Shah'
    return 'Cooper'

def STHE_boiling_HTC(AS, geom, cond, choice=None):
    """Shell-side boiling heat transfer coefficient [W/m2/K] (outer tube surface).
    choice: None/'auto' or the name of a correlation of _STHE_EVAP_HTC.
    cond : p (saturation pressure) [Pa], q (local heat flux on the outer surface) [W/m2]
    geom : roughness [m] (optional, Cooper)."""
    if choice in (None, 'auto'):
        choice = select_shell_evaporation_HTC(geom, cond)
    elif choice not in _STHE_EVAP_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_STHE_EVAP_HTC)} or None/'auto'.")
    return _STHE_EVAP_HTC[choice](AS, geom, cond)

#%% TWO-PHASE : SHELL-SIDE CONDENSATION (horizontal tube bundle)

_STHE_COND_HTC = {
    'Rose_Kern': STHE_rose_kern_condensing_HTC,  # Rose/Fujii/Nusselt single tube x Kern N^(-1/6)
}

def select_shell_condensation_HTC(geom, cond):
    """Automatic choice of the shell-side condensation HTC: Rose (single tube) with Kern
    inundation, as recommended by Bonneau et al. (2019)."""
    return 'Rose_Kern'

def STHE_condensation_HTC(AS, geom, cond, choice=None):
    """Shell-side condensation heat transfer coefficient [W/m2/K] (outer tube surface).
    choice: None/'auto' or the name of a correlation of _STHE_COND_HTC.
    cond : m_dot [kg/s], p [Pa], x [-], T_w [K]; geom : Kern geometry (+ N_rows optional)."""
    if choice in (None, 'auto'):
        choice = select_shell_condensation_HTC(geom, cond)
    elif choice not in _STHE_COND_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_STHE_COND_HTC)} or None/'auto'.")
    return _STHE_COND_HTC[choice](AS, geom, cond)

#%% GLOBAL DISPATCHER

_STATE_DISPATCH = {
    'SP': STHE_SP_HTC,
    'EVAP': STHE_boiling_HTC,
    'COND': STHE_condensation_HTC,
}

def STHE_HTC(AS, geom, cond, flow_state, choices=None):
    """Shell-side HTC of a shell & tube heat exchanger.
    flow_state : 'SP', 'EVAP' or 'COND' ('vapor-wet' etc. must be mapped to 'SP' by the caller if wanted)
    choices    : dict {'SP': name}. None, empty, 'auto' (or empty / 'auto' for a state) = automatic."""
    if flow_state not in _STATE_DISPATCH:
        raise NotImplementedError(f"STHE shell-side HTC: flow state {flow_state!r} not implemented "
                                  f"(available: {sorted(_STATE_DISPATCH)}).")

    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Kern'}}, "
                        f"got {type(choices).__name__}.")

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
# Bell_Delaware: leakage and bypass geometry (KeyError without it).
_REQUIRED_PARAMETERS = {
    'Kern': (),
    'Donohue': (),
    'Bell_Delaware': ('D_OTL', 'delta_TB', 'delta_BS', 'N_strips', 'TS_t', 'B_i', 'B_o'),
    'Cooper': (),
    'Mostinski': (),
    'Shah': (),
    'Rose_Kern': (),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"STHE_HTC: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
