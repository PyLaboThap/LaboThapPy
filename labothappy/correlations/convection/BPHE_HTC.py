# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:51:41 2026

@author: basil
"""

#%% SINGLE PHASE : HTC

from labothappy.correlations.convection.BPHE_HTC_correlation_bank import (
    chevron_deg, BPHE_sinnott_SP_HTC, BPHE_martin_SP_HTC, BPHE_kumar_SP_HTC,
    BPHE_thonon_SP_HTC, BPHE_muley_manglik_SP_HTC)

_BPHX_SP_HTC = {
    'Sinnott': BPHE_sinnott_SP_HTC,
    'Kumar': BPHE_kumar_SP_HTC,
    'Thonon': BPHE_thonon_SP_HTC,
    'Muley_Manglik': BPHE_muley_manglik_SP_HTC,
    'Martin': BPHE_martin_SP_HTC,
}
# These need cond['mu_w']; Martin treats it as optional.
_NEEDS_MU_W = {'Sinnott', 'Kumar', 'Muley_Manglik'}

def select_plate_HTC(geom, cond):
    """Automatic choice of the single-phase HTC correlation (heuristic).
    The rules mirror select_plate_friction_DP, so that HTC and friction DP of a
    side are by default taken from the same family of correlations."""
    if 'beta' not in geom:
        return 'Sinnott'                          # only one that needs no angle
    beta = chevron_deg(geom['beta'])              # from the flow direction (Martin's convention)
    beta_perp = 90.0 - beta                       # Kumar / Thonon tables
    Re = cond['G'] * geom['Dh'] / cond['mu']
    if Re >= 1000 and 30 <= beta <= 60 and 'phi' in geom:
        return 'Muley_Manglik'                    # best in its validated window, uses phi
    if geom.get('duty') in ('evaporator', 'condenser') and 50 <= Re <= 15000 \
            and 30 <= beta_perp <= 75:
        return 'Thonon'                           # designed for evaporators / condensers
    if Re < 100 and 30 <= beta_perp <= 65:
        return 'Kumar'                            # explicit laminar / transition rows
    return 'Martin'                               # continuous, any angle, any Re

def BPHE_SP_HTC(AS, geom, cond, choice=None):
    """Single-phase heat transfer coefficient of one side [W/m2/K].
    choice: None/'auto' or the name of a correlation.

    geom['beta'] stays in the library convention (angle from the main flow
    direction): Kumar and Thonon convert it to their own convention INSIDE
    their functions, so it must not be converted here (unlike the friction DP
    functions, whose dispatcher converts before the call).
    If cond has no 'mu_w' and the correlation needs it, the wall viscosity
    correction is set to 1 (mu_w = mu), on a copy of cond."""
    if choice in (None, 'auto'):
        choice = select_plate_HTC(geom, cond)
    elif choice not in _BPHX_SP_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_SP_HTC)} or None/'auto'.")
    c = cond
    if choice in _NEEDS_MU_W and 'mu_w' not in cond:
        c = dict(cond)
        c['mu_w'] = cond['mu']
            
    return _BPHX_SP_HTC[choice](AS, geom, c)

#%% TWO-PHASE : BOILING 

from labothappy.correlations.convection.BPHE_HTC_correlation_bank import (
    BPHE_cooper_boiling_HTC, BPHE_amalfi_boiling_HTC)

_BPHX_EVAP_HTC = {
    'Cooper': BPHE_cooper_boiling_HTC,
    'Amalfi': BPHE_amalfi_boiling_HTC,
}

def select_plate_evaporation_HTC(geom, cond):
    """Automatic choice of the evaporation HTC correlation (heuristic).
      - no chevron angle in geom        -> Cooper (pool boiling, needs no plate geometry;
                                           no quality dependence: rough estimate)
      - otherwise                       -> Amalfi (plate-specific, quality and mass flux
                                           dependent, built on a large multi-fluid database)
    """
    if 'beta' not in geom:
        return 'Cooper'
    return 'Amalfi'
 
def BPHE_boiling_HTC(AS, geom, cond, choice=None):
    """Evaporation heat transfer coefficient of one side [W/m2/K].
    choice: None/'auto' or the name of a correlation.
 
    cond : p (saturation pressure) [Pa], q (local heat flux) [W/m2] for both;
           G (mass flux in ONE channel) and x (mean quality) also for Amalfi.
    geom : Dh and beta (angle from the main flow direction, rad or deg) for Amalfi;
           roughness (optional) for Cooper.
    Amalfi needs the vapour viscosity when Bd >= 4: if it is not available for the
    fluid it raises an error, and choice='Cooper' must then be given explicitly."""
    if choice in (None, 'auto'):
        choice = select_plate_evaporation_HTC(geom, cond)
    elif choice not in _BPHX_EVAP_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_EVAP_HTC)} or None/'auto'.")
    return _BPHX_EVAP_HTC[choice](AS, geom, cond)

#%% TWO-PHASE : CONDENSING 

from labothappy.correlations.convection.BPHE_HTC_correlation_bank import (
    BPHE_shah_condensing_HTC, BPHE_han_condensing_HTC, BPHE_longo_condensing_HTC)
  
_BPHX_COND_HTC = {
    'Shah': BPHE_shah_condensing_HTC,
    'Han': BPHE_han_condensing_HTC,
    'Longo': BPHE_longo_condensing_HTC,
}
 
#%% CONDENSATION : HTC
 
def select_plate_condensation_HTC(geom, cond):
    """Automatic choice of the condensation HTC correlation (heuristic).
      - enlargement factor phi known                -> Shah (wide database, gravity and
                                                       forced-convection regimes)
      - p_co and beta known, inside Han's data      -> Han (plate-specific: corrugation pitch
        range (G 13-34 kg/m2/s, 45-70 deg from the      and angle; data of R410A / R22 only)
        flow direction)
      - otherwise                                   -> Longo (ACHP form: needs only Dh, G, x, p)
    """
    if 'phi' in geom:
        return 'Shah'
    if 'p_co' in geom and 'beta' in geom:
        theta = chevron_deg(geom['beta'])         # from the main flow direction
        if 13 <= cond['G'] <= 34 and 45 <= theta <= 70:
            return 'Han'
    return 'Longo'
 
def BPHE_condensing_HTC(AS, geom, cond, choice=None):
    """Condensation heat transfer coefficient of one side [W/m2/K].
    choice: None/'auto' or the name of a correlation.
 
    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-],
           p (saturation pressure) [Pa]
    geom : Dh for all; phi (enlargement factor) for Shah (default 1.17, with a warning);
           p_co (corrugation pitch) and beta (angle from the main flow direction,
           rad or deg) for Han.
    geom['beta'] stays in the library convention: Han converts it to its own
    convention inside its function, so it is not converted here."""
    if choice in (None, 'auto'):
        choice = select_plate_condensation_HTC(geom, cond)
    elif choice not in _BPHX_COND_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_COND_HTC)} or None/'auto'.")
    return _BPHX_COND_HTC[choice](AS, geom, cond)

#%% GLOBAL DISPATCHER

_STATE_DISPATCH = {
    'SP': BPHE_SP_HTC,
    'EVAP': BPHE_boiling_HTC,
    'COND': BPHE_condensing_HTC,
}

def BPHE_HTC(AS, geom, cond, flow_state, choices=None):
    """HTC d'un côté d'un échangeur à plaques brasées.
    flow_state : 'SP', 'EVAP' ou 'COND'
    choices    : dict {'SP': nom, 'EVAP': nom, 'COND': nom}.
                 None, vide, 'auto' (ou valeur vide / 'auto' pour un état) = sélection automatique."""
                 
    if flow_state not in _STATE_DISPATCH:
        raise ValueError(f"Unknown flow state {flow_state!r}: expected one of {sorted(_STATE_DISPATCH)}.")

    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}                                  # None, {}, '', 0, False, 'auto' -> tout en auto
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Martin'}}, "
                        f"got {type(choices).__name__}.")

    unknown = set(choices) - set(_STATE_DISPATCH)
    if unknown:
        raise ValueError(f"Unknown keys in choices: {sorted(unknown)} (expected {sorted(_STATE_DISPATCH)}).")

    choice = choices.get(flow_state)
    if not choice or (isinstance(choice, str) and choice.lower() == 'auto'):
        choice = None                                 # -> sélection automatique du dispatcher de l'état
    return _STATE_DISPATCH[flow_state](AS, geom, cond, choice)
