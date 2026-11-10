# -*- coding: utf-8 -*-
"""
Shell-side pressure drop dispatcher of shell & tube heat exchangers (single phase only for now),
same style as STHE_HTC_dispatcher. The automatic choice mirrors select_shell_HTC, so that the HTC and
the pressure drop of the shell side come from the same family of correlations.
"""

from labothappy.correlations.pressure_drop.STHE_DP_correlation_bank import (
    STHE_kern_shell_DP, STHE_donohue_shell_DP, STHE_bell_delaware_shell_DP)
from labothappy.correlations.convection.STHE_HTC import select_shell_HTC, _missing

_STHE_SP_DP = {
    'Kern': STHE_kern_shell_DP,
    'Donohue': STHE_donohue_shell_DP,
    'Bell_Delaware': STHE_bell_delaware_shell_DP,
}


def select_shell_DP(geom, cond):
    return select_shell_HTC(geom, cond)


def STHE_SP_DP(AS, geom, cond, choice=None):
    """Single-phase shell-side pressure drop [Pa] (whole shell side)."""
    if choice in (None, 'auto'):
        choice = select_shell_DP(geom, cond)
    elif choice not in _STHE_SP_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_STHE_SP_DP)} or None/'auto'.")
    missing = _missing(geom, choice)
    if missing:
        raise KeyError(f"{choice} shell pressure drop: missing geometry keys {missing}.")
    return _STHE_SP_DP[choice](AS, geom, cond)


_STATE_DISPATCH = {'SP': STHE_SP_DP}


def STHE_DP(AS, geom, cond, flow_state, choices=None):
    """Shell-side pressure drop of a shell & tube heat exchanger. flow_state: 'SP' (only one implemented)."""
    if flow_state not in _STATE_DISPATCH:
        raise NotImplementedError(f"STHE shell-side pressure drop: flow state {flow_state!r} not implemented "
                                  f"(available: {sorted(_STATE_DISPATCH)}).")
    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Kern'}}, got {type(choices).__name__}.")
    unknown = set(choices) - set(_STATE_DISPATCH)
    if unknown:
        raise ValueError(f"Unknown keys in choices: {sorted(unknown)} (expected {sorted(_STATE_DISPATCH)}).")
    choice = choices.get(flow_state)
    if not choice or (isinstance(choice, str) and choice.lower() == 'auto'):
        choice = None
    return _STATE_DISPATCH[flow_state](AS, geom, cond, choice)

