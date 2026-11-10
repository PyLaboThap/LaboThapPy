import numpy as np
import warnings

#%% SINGLE PHASE : FRICTION DP

from labothappy.correlations.pressure_drop.BPHX_DP_correlation_bank import (
    chevron_deg, BPHX_sinnott_SP_DP, BPHX_martin_SP_DP, BPHX_kumar_SP_DP,
    BPHX_thonon_SP_DP, BPHX_muley_manglik_SP_DP)

_BPHX_SP_DP = {
    'Sinnott': BPHX_sinnott_SP_DP,
    'Kumar': BPHX_kumar_SP_DP,
    'Thonon': BPHX_thonon_SP_DP,
    'Muley_Manglik': BPHX_muley_manglik_SP_DP,
    'Martin': BPHX_martin_SP_DP,
}
_PERP_CONVENTION = {'Kumar', 'Thonon'}   # angle from the perpendicular (to be verified)

def select_plate_DP(geom, cond):
    """Automatic choice of the single-phase friction correlation (heuristic).
    The rules are those of select_plate_HTC, so that HTC and friction DP of a
    side are by default taken from the same family of correlations."""
    if 'beta' not in geom:
        return 'Sinnott'                          # only one that needs no angle
    beta = chevron_deg(geom['beta'])              # from the flow direction (Martin's convention)
    beta_perp = 90.0 - beta                       # Kumar / Thonon tables
    Re = cond['G'] * geom['Dh'] / cond['mu']
    if Re >= 1000 and 30 <= beta <= 60 and 'phi' in geom:
        return 'Muley_Manglik'
    if geom.get('duty') in ('evaporator', 'condenser') and 50 <= Re <= 15000 \
            and 30 <= beta_perp <= 75:
        return 'Thonon'
    if Re < 100 and 30 <= beta_perp <= 65:
        return 'Kumar'
    return 'Martin'

def BPHX_SP_DP(AS, geom, cond, choice=None):
    """Single-phase friction pressure drop of one side [Pa], over the full flow length.
    choice: None/'auto' or the name of a correlation.

    cond : G, mu, D (density)
    geom : Dh, L, beta (angle from the main flow direction), plus n_passes and D_port
           (optional) and phi (Martin, Muley-Manglik).

    Unlike the HTC functions, the friction functions of Kumar and Thonon do NOT convert
    the angle themselves: it is converted here, to the perpendicular convention, on a copy of geom."""
    if choice in (None, 'auto'):
        choice = select_plate_DP(geom, cond)
    elif choice not in _BPHX_SP_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_SP_DP)} or None/'auto'.")
    g = geom
    if choice in _PERP_CONVENTION:
        g = dict(geom)
        g['beta'] = np.radians(90.0 - chevron_deg(geom['beta']))
    return _BPHX_SP_DP[choice](AS, g, cond)

#%% TWO-PHASE : BOILING

from labothappy.correlations.pressure_drop.BPHX_DP_correlation_bank import BPHX_lockhart_martinelli_DP

_BPHX_EVAP_DP = {
    'Lockhart_Martinelli': BPHX_lockhart_martinelli_DP,
}

def select_plate_evaporation_DP(geom, cond):
    """Automatic choice of the evaporation friction correlation (heuristic).
      - only Lockhart-Martinelli is available for the moment (separated flow model, Chisholm form;
        extrapolated from smooth tubes, with Martin's friction factor for each phase)."""
    return 'Lockhart_Martinelli'

def BPHX_boiling_DP(AS, geom, cond, choice=None):
    """Evaporation friction pressure drop of one side [Pa], over the full flow length.
    choice: None/'auto' or the name of a correlation.

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-], p (saturation pressure) [Pa]
    geom : Dh, L, beta; optional lm_C and lm_plate_friction (Lockhart-Martinelli).
    Acceleration, gravity and port losses are not included."""
    if choice in (None, 'auto'):
        choice = select_plate_evaporation_DP(geom, cond)
    elif choice not in _BPHX_EVAP_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_EVAP_DP)} or None/'auto'.")
    return _BPHX_EVAP_DP[choice](AS, geom, cond)

#%% TWO-PHASE : CONDENSING

from labothappy.correlations.pressure_drop.BPHX_DP_correlation_bank import BPHX_han_condensing_DP

_BPHX_COND_DP = {
    'Han': BPHX_han_condensing_DP,
    'Lockhart_Martinelli': BPHX_lockhart_martinelli_DP,
}

def select_plate_condensation_DP(geom, cond):
    """Automatic choice of the condensation friction correlation (heuristic).
      - p_co and beta known, inside Han's data      -> Han (plate-specific: corrugation pitch and angle;
        range (G 13-34 kg/m2/s, 45-70 deg from the      data of R410A / R22 only), same rule as the
        flow direction)                                 Han condensation HTC
      - otherwise                                   -> Lockhart-Martinelli (separated flow model)
    """
    if 'p_co' in geom and 'beta' in geom:
        theta = chevron_deg(geom['beta'])         # from the main flow direction
        if 13 <= cond['G'] <= 34 and 45 <= theta <= 70:
            return 'Han'
    return 'Lockhart_Martinelli'

def BPHX_condensing_DP(AS, geom, cond, choice=None):
    """Condensation friction pressure drop of one side [Pa], over the full flow length.
    choice: None/'auto' or the name of a correlation.

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-], p (saturation pressure) [Pa]
    geom : Dh, L, beta for both; p_co for Han; optional lm_C and lm_plate_friction (Lockhart-Martinelli).
    Acceleration, gravity and port losses are not included."""
    if choice in (None, 'auto'):
        choice = select_plate_condensation_DP(geom, cond)
    elif choice not in _BPHX_COND_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_COND_DP)} or None/'auto'.")
    return _BPHX_COND_DP[choice](AS, geom, cond)

#%% PORT DP

def BPHX_port_DP(geom, cond):
    """
    Port / manifold pressure loss of ONE side of a plate heat exchanger [Pa],
    for the whole exchanger (inlet + outlet ports together, not per cell).

        dP_port = 1.3 * G_port^2 / (2 rho) * n_passes      (Sinnott)

    G_port : cond['G_port'] if given, else m_dot / (pi D_port^2 / 4) with
             geom['D_port'] [m] (port diameter) and cond['m_dot'] (side flow).
    rho    : cond['D'], evaluated by the caller (inlet state is a good choice).
    Without port information the loss is NOT estimated: a warning is raised
    and 0 is returned.
    """
    if 'G_port' in cond:
        G_port = cond['G_port']
    elif 'D_port' in geom:
        G_port = cond['m_dot'] / (np.pi * geom['D_port']**2 / 4)
    else:
        warnings.warn("BPHX_port_DP: no port size ('D_port' or 'G_port'): port loss set to 0.")
        return 0.0
    return 1.3 * G_port**2 / (2 * cond['D']) * geom.get('n_passes', 1)

#%% GLOBAL DISPATCHER

_STATE_DISPATCH_DP = {
    'SP': BPHX_SP_DP,
    'EVAP': BPHX_boiling_DP,
    'COND': BPHX_condensing_DP,
}

def BPHE_DP(AS, geom, cond, flow_state, choices=None):
    """Perte de charge par frottement d'un côté d'un échangeur à plaques brasées [Pa],
    sur toute la longueur d'écoulement (les pertes aux ports sont à part : BPHX_port_DP).
    flow_state : 'SP', 'EVAP' ou 'COND'
    choices    : dict {'SP': nom, 'EVAP': nom, 'COND': nom}.
                 None, vide, 'auto' (ou valeur vide / 'auto' pour un état) = sélection automatique."""
    if flow_state not in _STATE_DISPATCH_DP:
        raise ValueError(f"Unknown flow state {flow_state!r}: expected one of {sorted(_STATE_DISPATCH_DP)}.")

    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}                                  # None, {}, '', 0, False, 'auto' -> tout en auto
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Martin'}}, "
                        f"got {type(choices).__name__}.")

    unknown = set(choices) - set(_STATE_DISPATCH_DP)
    if unknown:
        raise ValueError(f"Unknown keys in choices: {sorted(unknown)} (expected {sorted(_STATE_DISPATCH_DP)}).")

    choice = choices.get(flow_state)
    if not choice or (isinstance(choice, str) and choice.lower() == 'auto'):
        choice = None                                 # -> sélection automatique du dispatcher de l'état
    return _STATE_DISPATCH_DP[flow_state](AS, geom, cond, choice)