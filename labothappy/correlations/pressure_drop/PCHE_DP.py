# -*- coding: utf-8 -*-
"""
PCHE pressure drop dispatcher: flow state (SP / EVAP / COND) -> correlation, (AS, geom, cond) formalism.

geom['type_channel'] : 'Straight' (default), 'Zigzag' or 'S_fin' (S-shaped fins).

SP    : PCHE correlation 'Lee_Jiang' (zigzag channel: Lee et al. 2024 for Re <= 2000, Jiang et al.
        2018 for Re >= 3000, linear blend in between; PCHE_DP_correlation_bank.py)
        or any pipe correlation name ('Churchill', 'Haaland', ...), evaluated by pipe_DP on the
        equivalent tube (D_t_i = D_h, L = channel length, G of one channel).
        'Ngo' (zigzag fins 52 deg, sCO2, 3500 < Re < 22000, Ngo et al. 2007): explicit choice only.
        'Ngo_S_fin' (S-shaped fins 52 deg, sCO2, 3500 < Re < 23000, Ngo et al. 2007).
        'Kim_No_2012' (zigzag, laminar Re < 2500, helium, angle 5-45 deg, Kim & No 2012; optional
        geom['P_zz']), 'Kim_No_2011' (zigzag 15 deg, laminar, helium / water),
        'Liu_straight' (straight, sCO2 cooling, 3600 < Re < 36500, Liu et al. 2020),
        'Meshram' (zigzag, sCO2, 5000 < Re < 32000, Meshram et al. 2016): explicit choice only.
        Zigzag, explicit choice only (2026-11-23): 'Katz' (37 deg, sCO2 / helium), 'Kim_SG_2016'
        (32.5 / 40 deg, sCO2 CFD), 'Baik' (32.5 deg, CO2 / water, modified diameter),
        'Saeed' (40 deg, sCO2 / water), 'Chen_2016' (15 deg, helium).
EVAP / COND (2026-11-23): 'Bae' (CO2 near the critical point, zigzag 32.5 deg, Bae et al. 2019:
        Baik gas-only f x HEM Cicchitti multiplier x form-loss correction): explicit choice only.
        auto : Zigzag -> 'Lee_Jiang' ; S_fin -> 'Ngo_S_fin' ; Straight -> pipe_DP auto choice.
EVAP / COND : pipe_DP on the equivalent tube (pipe two-phase correlations, auto as pipe_DP).

Returns the frictional pressure drop [Pa] over the full channel length geom['L'].
(Replaces the former content of this file, an unused copy of pipe_DP.py; PCHE was routed to pipe_DP.)
"""

from labothappy.correlations.pressure_drop.PCHE_DP_correlation_bank import (
    PCHE_lee_jiang_zigzag_SP_DP)
from labothappy.correlations.pressure_drop.PCHE_DP_correlation_bank import PCHE_ngo_zigzag_SP_DP
from labothappy.correlations.pressure_drop.PCHE_DP_correlation_bank import PCHE_ngo_S_fin_SP_DP
from labothappy.correlations.pressure_drop.PCHE_DP_correlation_bank import (
    PCHE_kim_no_2012_zigzag_SP_DP, PCHE_kim_no_2011_zigzag_SP_DP, PCHE_liu_straight_SP_DP,
    PCHE_meshram_zigzag_SP_DP)
from labothappy.correlations.pressure_drop.PCHE_DP_correlation_bank import (
    PCHE_katz_zigzag_SP_DP, PCHE_kim_sg_zigzag_SP_DP, PCHE_baik_zigzag_SP_DP,
    PCHE_saeed_zigzag_SP_DP, PCHE_chen_zigzag_SP_DP, PCHE_bae_two_phase_DP)
from labothappy.correlations.pressure_drop.pipe_DP import pipe_DP

_PCHE_SP_DP = {
    'Lee_Jiang': PCHE_lee_jiang_zigzag_SP_DP,
    'Ngo': PCHE_ngo_zigzag_SP_DP,      # zigzag fins 52 deg, sCO2 (explicit choice only)
    'Ngo_S_fin': PCHE_ngo_S_fin_SP_DP,  # S-shaped fins 52 deg, sCO2 (auto for 'S_fin')
    'Kim_No_2012': PCHE_kim_no_2012_zigzag_SP_DP,  # zigzag laminar, helium (explicit only)
    'Kim_No_2011': PCHE_kim_no_2011_zigzag_SP_DP,  # zigzag 15 deg laminar (explicit only)
    'Liu_straight': PCHE_liu_straight_SP_DP,       # straight, sCO2 cooling (explicit only)
    'Meshram': PCHE_meshram_zigzag_SP_DP,          # zigzag, sCO2 (explicit only)
    'Katz': PCHE_katz_zigzag_SP_DP,                # 37 deg, sCO2 / helium (explicit only)
    'Kim_SG_2016': PCHE_kim_sg_zigzag_SP_DP,       # 32.5 / 40 deg, sCO2 CFD (explicit only)
    'Baik': PCHE_baik_zigzag_SP_DP,                # 32.5 deg, CO2 / water (explicit only)
    'Saeed': PCHE_saeed_zigzag_SP_DP,              # 40 deg, sCO2 / water, Darcy (explicit only)
    'Chen_2016': PCHE_chen_zigzag_SP_DP,           # 15 deg, helium (explicit only)
}

_VALID_STATES = ('SP', 'EVAP', 'COND')

# Two-phase PCHE pressure drop (EVAP and COND, explicit choice only)
_PCHE_2P_DP = {'Bae': PCHE_bae_two_phase_DP}   # CO2 near the critical point (Bae et al. 2019)


def select_PCHE_DP(geom, cond):
    """Zigzag channel -> 'Lee_Jiang'; S-shaped fins -> 'Ngo_S_fin';
    straight channel -> None (pipe_DP auto choice)."""
    type_channel = str(geom.get('type_channel', 'Straight')).lower()
    if type_channel == 'zigzag':
        return 'Lee_Jiang'
    if type_channel == 's_fin':
        return 'Ngo_S_fin'
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
    # Two phase: PCHE two-phase correlation if named (explicit only), else pipe correlations
    if choice in _PCHE_2P_DP:
        return _PCHE_2P_DP[choice](AS, geom, cond)
    return pipe_DP(AS, geom, cond, flow_state, {flow_state: choice})


# =============================================================================
# PARAMETERS NEEDED BY EACH CORRELATION IN ADDITION TO THE GEOMETRY
# =============================================================================

# Model parameters that each registered correlation needs IN ADDITION to the
# parameters of the heat exchanger geometry (hex_geometry required_parameters): keys
# that the geometry does not compute. () = the geometry parameters are enough.
# Found on 2026-11-25 by removing the keys one at a time at hex_MB states
# (correlation_warnings.parameters_by_removal) and by reading the code.
# alpha: bending angle of the zigzag / S-fin channels (error without it).
# D_ch: Baik and Bae assume Dh without it (AssumedParameterWarning).
_REQUIRED_PARAMETERS = {
    'Lee_Jiang': ('alpha',),
    'Ngo': ('alpha',),
    'Ngo_S_fin': ('alpha',),
    'Kim_No_2012': ('alpha',),
    'Kim_No_2011': ('alpha',),
    'Liu_straight': (),
    'Meshram': (),
    'Katz': ('alpha',),
    'Kim_SG_2016': ('alpha',),
    'Baik': ('alpha', 'D_ch'),
    'Saeed': ('alpha',),
    'Chen_2016': ('alpha',),
    'Bae': ('D_ch',),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"PCHE_DP: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
