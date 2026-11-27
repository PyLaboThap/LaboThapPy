import numpy as np
import warnings

#%% SINGLE PHASE : FRICTION DP

from labothappy.correlations.pressure_drop.BPHE_DP_correlation_bank import (
    chevron_deg, BPHE_sinnott_SP_DP, BPHE_martin_SP_DP, BPHE_kumar_SP_DP,
    BPHE_thonon_SP_DP, BPHE_muley_manglik_SP_DP)
from labothappy.correlations.convection.correlation_warnings import AssumedParameterWarning

_BPHE_SP_DP = {
    'Sinnott': BPHE_sinnott_SP_DP,
    'Kumar': BPHE_kumar_SP_DP,
    'Thonon': BPHE_thonon_SP_DP,
    'Muley_Manglik': BPHE_muley_manglik_SP_DP,
    'Martin': BPHE_martin_SP_DP,
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

def BPHE_SP_DP(AS, geom, cond, choice=None):
    """Single-phase friction pressure drop of one side [Pa], over the full flow length.
    choice: None/'auto' or the name of a correlation.

    cond : G, mu, D (density)
    geom : Dh, L, beta (angle from the main flow direction), plus n_passes and D_port
           (optional) and phi (Martin, Muley-Manglik).

    Unlike the HTC functions, the friction functions of Kumar and Thonon do NOT convert
    the angle themselves: it is converted here, to the perpendicular convention, on a copy of geom."""
    if choice in (None, 'auto'):
        choice = select_plate_DP(geom, cond)
    elif choice not in _BPHE_SP_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHE_SP_DP)} or None/'auto'.")
    g = geom
    if choice in _PERP_CONVENTION:
        g = dict(geom)
        g['beta'] = np.radians(90.0 - chevron_deg(geom['beta']))
    return _BPHE_SP_DP[choice](AS, g, cond)

#%% TWO-PHASE : BOILING

from labothappy.correlations.pressure_drop.BPHE_DP_correlation_bank import (
    BPHE_lockhart_martinelli_DP, BPHE_amalfi_boiling_DP)
from labothappy.correlations.pressure_drop.BPHE_DP_correlation_bank import BPHE_hsieh_lin_boiling_DP

_BPHE_EVAP_DP = {
    'Lockhart_Martinelli': BPHE_lockhart_martinelli_DP,
    'Amalfi': BPHE_amalfi_boiling_DP,
    'Hsieh_Lin': BPHE_hsieh_lin_boiling_DP,          # R410A data (Hsieh & Lin 2002)
}

def select_plate_evaporation_DP(geom, cond):
    """Automatic choice of the evaporation friction correlation (heuristic), same rule as
    select_plate_evaporation_HTC so that HTC and DP come from the same paper:
      - chevron angle in geom  -> Amalfi (plate-specific, Amalfi et al. 2016 Part 2)
      - otherwise              -> Lockhart-Martinelli (separated flow model, Chisholm form;
        extrapolated from smooth tubes, with Martin's friction factor for each phase)."""
    if 'beta' in geom:
        return 'Amalfi'
    return 'Lockhart_Martinelli'

def BPHE_boiling_DP(AS, geom, cond, choice=None):
    """Evaporation friction pressure drop of one side [Pa], over the full flow length.
    choice: None/'auto' or the name of a correlation.

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-], p (saturation pressure) [Pa]
    geom : Dh, L, beta; optional lm_C and lm_p_f (Lockhart-Martinelli).
    Amalfi reads geom['beta'] directly (from the main flow direction, deg or rad).
    Acceleration, gravity and port losses are not included."""
    if choice in (None, 'auto'):
        choice = select_plate_evaporation_DP(geom, cond)
    elif choice not in _BPHE_EVAP_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHE_EVAP_DP)} or None/'auto'.")
    return _BPHE_EVAP_DP[choice](AS, geom, cond)

#%% TWO-PHASE : CONDENSING

from labothappy.correlations.pressure_drop.BPHE_DP_correlation_bank import BPHE_han_condensing_DP
from labothappy.correlations.pressure_drop.BPHE_DP_correlation_bank import (
    BPHE_kuo_condensing_DP, BPHE_yan_lio_lin_condensing_DP)

_BPHE_COND_DP = {
    'Han': BPHE_han_condensing_DP,
    'Lockhart_Martinelli': BPHE_lockhart_martinelli_DP,
    'Kuo': BPHE_kuo_condensing_DP,                   # R410A data (Kuo et al. 2005)
    'Yan_Lio_Lin': BPHE_yan_lio_lin_condensing_DP,   # R134a data (Yan, Lio & Lin 1999)
}

def select_plate_condensation_DP(geom, cond):
    """Automatic choice of the condensation friction correlation (heuristic).
      - Lambda and beta known, inside Han's data      -> Han (plate-specific: corrugation pitch and angle;
        range (G 13-34 kg/m2/s, 45-70 deg from the      data of R410A / R22 only), same rule as the
        flow direction)                                 Han condensation HTC
      - otherwise                                   -> Lockhart-Martinelli (separated flow model)
    """
    if 'Lambda' in geom and 'beta' in geom:
        theta = chevron_deg(geom['beta'])         # from the main flow direction
        if 13 <= cond['G'] <= 34 and 45 <= theta <= 70:
            return 'Han'
    return 'Lockhart_Martinelli'

def BPHE_condensing_DP(AS, geom, cond, choice=None):
    """Condensation friction pressure drop of one side [Pa], over the full flow length.
    choice: None/'auto' or the name of a correlation.

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-], p (saturation pressure) [Pa]
    geom : Dh, L, beta for both; Lambda for Han; optional lm_C and lm_p_f (Lockhart-Martinelli).
    Acceleration, gravity and port losses are not included."""
    if choice in (None, 'auto'):
        choice = select_plate_condensation_DP(geom, cond)
    elif choice not in _BPHE_COND_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHE_COND_DP)} or None/'auto'.")
    return _BPHE_COND_DP[choice](AS, geom, cond)

#%% PORT DP

def BPHE_port_DP(geom, cond):
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
        warnings.warn("BPHE_port_DP: no port size ('D_port' or 'G_port'): port loss set to 0.",
                      AssumedParameterWarning)
        return 0.0
    return 1.3 * G_port**2 / (2 * cond['D']) * geom.get('n_passes', 1)

#%% GLOBAL DISPATCHER

_STATE_DISPATCH_DP = {
    'SP': BPHE_SP_DP,
    'EVAP': BPHE_boiling_DP,
    'COND': BPHE_condensing_DP,
}

def BPHE_DP(AS, geom, cond, flow_state, choices=None):
    """Perte de charge par frottement d'un côté d'un échangeur à plaques brasées [Pa],
    sur toute la longueur d'écoulement (les pertes aux ports sont à part : BPHE_port_DP).
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


# =============================================================================
# PARAMETERS NEEDED BY EACH CORRELATION IN ADDITION TO THE GEOMETRY
# =============================================================================

# Model parameters that each registered correlation needs IN ADDITION to the
# parameters of the heat exchanger geometry (hex_geometry required_parameters): keys
# that the geometry does not compute. () = the geometry parameters are enough.
# Found on 2026-11-25 by removing the keys one at a time at hex_MB states
# (correlation_warnings.parameters_by_removal) and by reading the code.
_REQUIRED_PARAMETERS = {
    'Sinnott': (),
    'Kumar': (),
    'Thonon': (),
    'Muley_Manglik': (),
    'Martin': (),
    'Lockhart_Martinelli': (),
    'Amalfi': (),
    'Han': (),
    'Hsieh_Lin': (),
    'Kuo': (),
    'Yan_Lio_Lin': (),
}


def required_parameters(name):
    """Model parameters that the correlation `name` needs in addition to the geometry
    parameters (tuple, empty if none). Raises KeyError for a name that is not registered here."""
    if name not in _REQUIRED_PARAMETERS:
        raise KeyError(f"BPHE_DP: unknown correlation {name!r}; registered: "
                       f"{list(_REQUIRED_PARAMETERS)}.")
    return _REQUIRED_PARAMETERS[name]
