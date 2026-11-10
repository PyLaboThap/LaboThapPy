# -*- coding: utf-8 -*-
"""
Shell-side HTC dispatcher of shell & tube heat exchangers, same style as BPHX_HTC_dispatcher:
HX type -> flow state (SP / EVAP / COND) -> correlation, with the (AS, geom, cond) formalism.
Only single-phase shell-side flow is implemented for now.
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
_KERN_KEYS = ('Tube_OD', 'pitch_ratio', 'Shell_ID', 'central_spacing', 'tube_layout')
_DONOHUE_KEYS = _KERN_KEYS + ('Baffle_cut', 'n_tubes')
_BELL_KEYS = _DONOHUE_KEYS + ('D_OTL', 'clear_TB', 'clear_BS', 'N_strips',
                              'Tube_L', 'Tubesheet_t', 'inlet_spacing', 'outlet_spacing')
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

#%% GLOBAL DISPATCHER

_STATE_DISPATCH = {
    'SP': STHE_SP_HTC,
}

def STHE_HTC(AS, geom, cond, flow_state, choices=None):
    """Shell-side HTC of a shell & tube heat exchanger.
    flow_state : 'SP' ('vapor-wet' etc. must be mapped to 'SP' by the caller if wanted)
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