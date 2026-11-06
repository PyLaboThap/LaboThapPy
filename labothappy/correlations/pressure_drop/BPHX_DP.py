
import numpy as np
import warnings
from labothappy.correlations.pressure_drop.plate_DP_correlation_bank import (chevron_deg, sinnott_plate_DP, martin_plate_DP,
                                                                             kumar_plate_DP, thonon_plate_DP, muley_manglik_plate_DP)

#%% SINGLE PHASE

_BPHX_FRICTION_SP_DP = {
    'Sinnott': sinnott_plate_DP,
    'Kumar': kumar_plate_DP,
    'Thonon': thonon_plate_DP,
    'Muley_Manglik': muley_manglik_plate_DP,
    'Martin': martin_plate_DP,
}
_PERP_CONVENTION = {'Kumar', 'Thonon'}   # angle from the perpendicular (to be verified)

#%% PORT DP

def BHPX_port_DP(geom, cond):
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
        warnings.warn("plate_port_DP: no port size ('D_port' or 'G_port'): port loss set to 0.")
        return 0.0
    return 1.3 * G_port**2 / (2 * cond['D']) * geom.get('n_passes', 1)

#%% SINGLE PHASE : FRICTION DP

def select_plate_friction_DP(geom, cond):
    """Automatic choice of the friction correlation (heuristic)."""
    if 'beta' not in geom:
        return 'Sinnott'
    beta = chevron_deg(geom['beta'])             # from the flow direction
    beta_perp = 90.0 - beta
    Re = cond['G'] * geom['Dh'] / cond['mu']
    if Re >= 1000 and 30 <= beta <= 60 and 'phi' in geom:
        return 'Muley_Manglik'
    if geom.get('duty') in ('evaporator', 'condenser') and 50 <= Re <= 15000 \
            and 30 <= beta_perp <= 75:
        return 'Thonon'
    if Re < 100 and 30 <= beta_perp <= 65:
        return 'Kumar'
    return 'Martin'

def BHPX_friction_SP_DP(AS, geom, cond, choice=None):
    """Friction pressure drop over the full flow length [Pa]. choice: None/'auto' or a name."""
    if choice in (None, 'auto'):
        choice = select_plate_friction_DP(geom, cond)
    elif choice not in _BPHX_FRICTION_SP_DP:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_BPHX_FRICTION_SP_DP)} or None/'auto'.")
    g = geom
    if choice in _PERP_CONVENTION:
        g = dict(geom)
        g['beta'] = np.radians(90.0 - chevron_deg(geom['beta']))
    return _BPHX_FRICTION_SP_DP[choice](AS, g, cond)

#%% TWO PHASE : EVAPORATION



