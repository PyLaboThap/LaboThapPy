"""
dp_corr_distributor.py

Single entry point to compute a pressure drops from any correlation,
for any heat exchanger/pipe geometry:

    dp = compute_dp(AS, m_dot, geom, dp_correlation, void_fraction_model)

The model only provides the fluid state (AS, at the mean state of the zone),
the mass flow rate and the geometry of one side. Each adapter below translates
these into the specific inputs of its correlation, including the translation
from one geometry type to another when a correlation is used outside the
geometry it was made for.

author: Elise
"""

import CoolProp.CoolProp as CP

from labothappy.correlations.pressure_drop.pipe_DP import (
    friction_factor_churchill,
    pressure_drop_muller_steinhagen_heck,
    pressure_drop_friedel
)
from labothappy.correlations.pressure_drop.plate_DP import (
    dp_martin_plate,
    dp_lockhart_martinelli
)
from labothappy.correlations.properties.dimensionless import compute_reynolds
from labothappy.correlations.properties.two_phase import get_saturated_phase_properties

from correlations.properties.thermal_conductivity import conducticity_R1233zd

import math
G_GRAVITY = 9.81  # m/s²
# ===============================================================
# Entry point
# ===============================================================

def compute_dp(AS, m_dot, geom, L, dp_correlation, void_fraction_model):
    """
    Pressure drop computation [Pa]

    AS          : CoolProp AbstractState, set at the MEAN state of the zone
    m_dot       : total mass flow rate on this side [kg/s]
    geom        : geometry of this side (PlateGeometry, TubeGeometry)
    L : Tube length
    dp_correlation : name of the correlation (key of _ADAPTERS)
    void_fraction_model
    """
    try:
        adapter = DP_CORRELATIONS[dp_correlation]
    except KeyError:
        raise ValueError(f"Unknown htc correlation '{dp_correlation}'. "
                         f"Available: {list(DP_CORRELATIONS)}") from None

    # Adapters may update AS (saturated or wall properties): save the zone
    # state and always restore it, so the caller's AS is never modified.
    p, h = AS.p(), AS.hmass()
    try:
        return adapter(AS, m_dot, geom, L, void_fraction_model)
    finally:
        AS.update(CP.HmassP_INPUTS, h, p)



# ===============================================================
# Adapters: all have the same signature (AS, m_dot, geom)
# ===============================================================

# ---------------- Single-phase ----------------

def _adapter_churchill(AS, m_dot, geom, L, void_fraction_model):
    """Darcy friction factor using Churchill correlation. Straight pipes originally."""

    p, h = AS.p(), AS.hmass()
    fluid = AS.fluid_names()[0]
    try:
        mu = AS.viscosity()
    except:
        try:
            mu = CP.PropsSI('V', 'P', p, 'H', h, fluid)
        except:
            raise RuntimeError(f"Problem computing viscosity at P={p}, h={h}")
    rho = AS.rhomass()

    # Mass flux through ONE channel
    G = m_dot / (geom['A_cs_channel'] * geom['n_channels'])
    v = G / rho # Mean velocity [m/s]

    Re = compute_reynolds(geom['D_h'], mu, rho, v)

    f = friction_factor_churchill(geom['roughness'], geom['D_h'], Re)

    # Frictional pressure drop (Darcy-Weisbach)
    dp= f * (L / geom['D_h']) * (rho * v ** 2 / 2.0)

    # Gravitational pressure drop (hydrostatic)
    # NEGLECTED FOR THE PHE SO FAR!!
    # theta = 
    # dP_gravity = G_GRAVITY * rho * L * math.sin(math.radians(theta))
    return dp


def _adapter_martin(AS, m_dot, geom, L, void_fraction_model):
    """Martin (VDI): chevron plates only"""
    if geom['chevron_angle'] is None:
        raise ValueError("Martin needs chevron_angle")

    # Single-phase properties at the current state
    p, h = AS.p(), AS.hmass()
    fluid = AS.fluid_names()[0]
    # Get other properties with safe guard in case the property doesn't exist
    try:
        mu = AS.viscosity()
        # print('mu', mu)
    except:
        try:
            mu = CP.PropsSI('V', 'P', p, 'H', h, fluid)
        except:
            raise RuntimeError(f"Problem computing viscosity at P={p}, h={h}")
    rho = AS.rhomass()

    # Mass flux through ONE channel
    G_ch = m_dot / (geom['A_cs_channel'] * geom['n_channels'])

    dp = dp_martin_plate(G_ch, rho, mu, geom['D_h'], geom['chevron_angle'], L)

    return dp

# ---------------- Two-phase ----------------

def _adapter_msh(AS, m_dot, geom, L, void_fraction_model):
    """Darcy friction factor using Churchill correlation. Straight pipes originally."""

    p, h = AS.p(), AS.hmass()
    fluid = AS.fluid_names()[0]

    props = get_saturated_phase_properties(AS)
    x = props["x"]
    rho_l = props["rho_l"]
    rho_v = props["rho_v"]
    mu_l = props["mu_l"]
    mu_v = props["mu_v"]
    sigma = props["sigma"]
    # Mass flux through ONE channel
    G = m_dot / (geom['A_cs_channel'] * geom['n_channels'])

    dp = pressure_drop_muller_steinhagen_heck(G, x, rho_l, rho_v, mu_l, mu_v, geom['D_h'], L, geom['roughness'])
    
    return dp

def _adapter_lockhart_martinelli(AS, m_dot, geom, L, void_fraction_model):
    """Darcy friction factor using Churchill correlation. Straight pipes originally."""

    p, h = AS.p(), AS.hmass()
    fluid = AS.fluid_names()[0]

    props = get_saturated_phase_properties(AS)
    x = props["x"]
    rho_l = props["rho_l"]
    rho_v = props["rho_v"]
    mu_l = props["mu_l"]
    mu_v = props["mu_v"]
    sigma = props["sigma"]
    # Mass flux through ONE channel
    G = m_dot / (geom['A_cs_channel'] * geom['n_channels'])

    dp = dp_lockhart_martinelli(G, x, rho_l, rho_v, mu_l, mu_v, geom['D_h'], L,
                                   C=4.67, chevron_angle=geom['chevron_angle'])

    # NOTE: C=4.67 according by Claesson (verifiy this based on ACHP)
    
    return dp

DP_CORRELATIONS = {
    'Churchill':     _adapter_churchill,
    'Martin':        _adapter_martin,
    'MSH':           _adapter_msh,
    'Lockhart-Martinelli': _adapter_lockhart_martinelli,
}


# dp_lockhart_martinelli_2phase(G, x, rho_l, rho_v, mu_l, mu_v, D_h, L,
#                                   C=None, chevron_angle=None):
# # ---------------- Condensation ----------------

# def _adapter_longo_cond(AS, m_dot, geom, q_flux, T_wall):
#     """Longo: condensation in brazed plates. Uses D_h = 2b, without Phi."""
#     if geom['corrugation_amplitude'] is None:
#       raise ValueError("Longo condensation needs corrugation_amplitude.")

#     x = AS.Q()                                   # read BEFORE _sat_props updates AS = x_avg over the cell
#     p = AS.p()                                   # read before updating AS

#     # Saturated liquid (bubble point for a mixture)
#     AS.update(CP.PQ_INPUTS, p, 0.0)
#     fluid = AS.fluid_names()[0]
#     rho_L = AS.rhomass()
#     try:
#         mu_L = AS.viscosity()
#     except:
#         try:
#             mu_L = CP.PropsSI('viscosity', 'P', p, 'Q', 0, fluid)
#         except:
#             raise RuntimeError(f"Problem computing viscosity at P={p}, Q=0")
#     cp_L = AS.cpmass()

#     if fluid == 'R1233ZDE' or fluid == 'R1233zd(E)':
#         T = AS.T()
#         k_L = conducticity_R1233zd(T, p)
#     else:
#         try:
#             k_L = AS.conductivity()
#         except:
#             try:
#                 k_L = CP.PropsSI('L', 'P', p, 'Q', 0, fluid)
#             except:
#                 raise RuntimeError(f"Problem computing conductivity at P={p}, Q=0")

#     # Saturated vapour (dew point for a mixture)
#     AS.update(CP.PQ_INPUTS, p, 1.0)
#     rho_V = AS.rhomass()
#     # Mass flux through ONE channel
#     G_ch = m_dot / (geom['A_cs_channel'] * geom['n_channels'])
#     htc = htc_longo_condensation(x=x, G=G_ch, D_h=geom['D_h'], rho_L=rho_L, rho_V=rho_V, mu_L=mu_L, cp_L=cp_L, k_L=k_L)

#     return htc


# HTC_CORRELATIONS = {
#     'Martin':     _adapter_martin,
#     'Cooper':     _adapter_cooper,
#     'Longo_cond': _adapter_longo_cond,
# }
