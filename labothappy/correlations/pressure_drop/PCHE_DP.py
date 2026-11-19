# -*- coding: utf-8 -*-
"""
PCHE pressure drop dispatcher: flow state (SP / EVAP / COND) -> correlation, (AS, geom, cond) formalism.

geom['type_channel'] : 'Straight' (default) or 'Zigzag'.

SP    : PCHE correlation 'Lee_Jiang' (zigzag channel: Lee et al. 2024 for Re <= 2000, Jiang et al.
        2018 for Re >= 3000, linear blend in between; PCHE_DP_correlation_bank.py)
        or any pipe correlation name ('Churchill', 'Haaland', ...), evaluated by pipe_DP on the
        equivalent tube (Tube_ID = D_h, L = channel length, G of one channel).
        auto : Zigzag -> 'Lee_Jiang' ; Straight -> pipe_DP auto choice.
EVAP / COND : pipe_DP on the equivalent tube (pipe two-phase correlations, auto as pipe_DP).

Returns the frictional pressure drop [Pa] over the full channel length geom['L'].
(Replaces the former content of this file, an unused copy of pipe_DP.py; PCHE was routed to pipe_DP.)
"""

from labothappy.correlations.pressure_drop.PCHE_DP_correlation_bank import (
    PCHE_lee_jiang_zigzag_SP_DP)
from labothappy.correlations.pressure_drop.pipe_DP import pipe_DP

_PCHE_SP_DP = {
    'Lee_Jiang': PCHE_lee_jiang_zigzag_SP_DP,
}

_VALID_STATES = ('SP', 'EVAP', 'COND')


def select_PCHE_DP(geom, cond):
    """Zigzag channel -> 'Lee_Jiang'; straight channel -> None (pipe_DP auto choice)."""
    if str(geom.get('type_channel', 'Straight')).lower() == 'zigzag':
        return 'Lee_Jiang'
    return None


def PCHE_SP_DP(AS, geom, cond, choice=None):
    """Single-phase pressure drop [Pa]: PCHE correlation, or pipe correlation on the equivalent tube."""
    if choice in (None, 'auto'):
        choice = select_PCHE_DP(geom, cond)
    if choice in _PCHE_SP_DP:
        return _PCHE_SP_DP[choice](AS, geom, cond)
    # Pipe correlation name (or None = pipe auto choice); pipe_DP raises on an unknown name
    return pipe_DP(AS, geom, cond, 'SP', {'SP': choice})


def PCHE_DP(AS, geom, cond, flow_state, choices=None):
    """Frictional pressure drop of one PCHE side, over the full channel length geom['L'] [Pa].
    flow_state : 'SP', 'EVAP' or 'COND'.
    choices    : dict {'SP': name, 'EVAP': name, 'COND': name}, or None / 'auto'."""
    if flow_state not in _VALID_STATES:
        raise ValueError(f"PCHE DP: flow state {flow_state!r}: expected one of {_VALID_STATES}.")
    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Lee_Jiang'}}, "
                        f"got {type(choices).__name__}.")
    choice = choices.get(flow_state)
    if isinstance(choice, str) and choice.lower() == 'auto':
        choice = None
    if flow_state == 'SP':
        return PCHE_SP_DP(AS, geom, cond, choice)
    # Two phase: pipe correlations on the equivalent tube
    return pipe_DP(AS, geom, cond, flow_state, {flow_state: choice})
