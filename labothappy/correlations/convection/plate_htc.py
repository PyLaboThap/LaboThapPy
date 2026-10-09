"""
author: Elise
"""

import numpy as np
from scipy.optimize import fsolve
import warnings

import CoolProp.CoolProp as CP


# ---------------------------------------------------------------
# Single-phase heat transfer coefficient correlations
# ---------------------------------------------------------------

def htc_martin_plate_1phase(G_ch, mu, cp, k, D_h, chevron_angle, mu_wall=None):
    """
    Single-phase heat transfer coefficient and pressure drop in a chevron plate
    heat exchanger channel, after H. Martin, VDI Heat Atlas, chapter N6.

    Inputs
    ------
    G_ch          : mass flux in ONE channel [kg/s] (= m_dot_side / (n_channels_side*A_cs))
    rho           : density [kg/m^3]
    mu            : dynamic viscosity at bulk temperature [Pa.s]
    cp            : specific heat capacity [J/kg/K]
    k             : thermal conductivity [W/m/K]
    D_h           : hydraulic diameter of one channel, 2 * b / enlargement_factor [m]
    chevron_angle : corrugation angle measured FROM THE MAIN FLOW DIRECTION [rad]
    mu_wall       : dynamic viscosity at wall temperature [Pa.s], optional

    Output
    -------
    htc [-]
    """
    # Dimensionless numbers
    Pr = cp * mu / k
    Re = G_ch * D_h / mu
    print('Re martin', Re)

    # Friction factors for phi = 0 (zeta0) and phi = 90 deg (zeta1_0)
    if Re < 2000:
        zeta0 = 64 / Re
        zeta1_0 = 597 / Re + 3.85
    else:
        zeta0 = (1.8 * np.log(Re) - 1.5) ** (-2)
        zeta1_0 = 39 / Re ** 0.289

    # Martin's fitting constants (a, b, c in the original paper)
    a_M, b_M, c_M = 3.8, 0.18, 0.36
    zeta1 = a_M * zeta1_0

    # Martin's equation 18: combined friction factor
    rhs = (np.cos(chevron_angle) / np.sqrt(b_M * np.tan(chevron_angle) + c_M * np.sin(chevron_angle) + zeta0 / np.cos(chevron_angle))
           + (1 - np.cos(chevron_angle)) / np.sqrt(zeta1))
    zeta = 1 / rhs**2

    # Hagen number
    Hg = zeta * Re**2 / 2
    
    # Viscosity correction (set to 1 if the wall viscosity is not given)
    visc_corr = (mu / mu_wall) ** (1 / 6) if mu_wall is not None else 1.0
    
    # Nusselt number and heat transfer coefficient
    c_q, q = 0.122, 0.374
    Nu = c_q * Pr ** (1 / 3) * visc_corr * (2 * Hg * np.sin(2 * chevron_angle)) ** q
    htc = Nu * k / D_h
    return htc


# ---------------------------------------------------------------
# Evaporating heat transfer coefficient correlations
# ---------------------------------------------------------------

def htc_cooper_pool_boiling(p, p_crit, q_flux, M_molar, roughness=None):
    """
    Nucleate pool boiling heat transfer coefficient, after M.G. Cooper (1984),
    "Heat flow rates in saturated nucleate boiling - A wide-ranging examination
    using reduced properties", Advances in Heat Transfer, Vol. 16, pp. 157-239.

    Inputs
    ------
    p         : saturation pressure [Pa]
    p_crit    : critical pressure of the fluid [Pa]
    q_flux    : heat flux through the wall [W/m^2]
    M_molar   : molar mass of the fluid [kg/mol] (CoolProp convention)
    roughness : surface roughness [m], optional. If None or 0, Cooper's
                recommended value of 1 micrometre is used.

    Output
    ------
    h : nucleate boiling heat transfer coefficient [W/m^2/K]
    """
    # Unit conversions to what Cooper's correlation expects
    p_red = p / p_crit                       # reduced pressure [-]
    M_gmol = M_molar * 1e3                   # [kg/mol] -> [g/mol]
    if roughness is None or roughness <= 0:
        Rp_um = 1.0                          # Cooper's default when unknown [um]
    else:
        Rp_um = roughness * 1e6              # [m] -> [um]

    htc = 55 * p_red ** (0.12 - 0.2 * np.log10(Rp_um)) * (-np.log10(p_red)) ** (-0.55) * q_flux ** 0.67 * M_gmol ** (-0.5)

    return htc


# ---------------------------------------------------------------
# Condensing heat transfer coefficient correlations
# ---------------------------------------------------------------

def htc_longo_condensation(x, G, D_h, rho_L, rho_V, mu_L, cp_L, k_L):
    """
    Condensation heat transfer coefficient in a brazed plate heat exchanger,
    after Longo (form as implemented in ACHP), based on the equivalent
    Reynolds number of Akers.

    Inputs
    ------
    x     : vapour quality, mean over the cell [-]
    G     : mass flux in ONE channel [kg/m^2/s] (= m_dot_side / (n_channels_side * A_cs))
    D_h   : hydraulic diameter [m] (Longo uses 2 * corrugation_depth, without the enlargement factor)
    rho_L : saturated liquid density [kg/m^3]
    rho_V : saturated vapour density [kg/m^3]
    mu_L  : saturated liquid dynamic viscosity [Pa.s]
    cp_L  : saturated liquid specific heat capacity [J/kg/K]
    k_L   : saturated liquid thermal conductivity [W/m/K]

    Output
    ------
    h : condensation heat transfer coefficient [W/m^2/K]
    """
    Pr_L = cp_L * mu_L / k_L

    # Equivalent Reynolds number (Akers): two-phase flow converted to an all-liquid flow
    Re_eq = G * ((1 - x) + x * np.sqrt(rho_L / rho_V)) * D_h / mu_L

    if Re_eq < 1750:
        Nu = 60 * Pr_L ** (1 / 3)
    else:
        if Re_eq > 3000:
            warnings.warn(f"Longo condensation: Re_eq = {Re_eq:.0f} > 3000, extrapolating.")
        Nu = ((75 - 60) / (3000 - 1750) * (Re_eq - 1750) + 60) * Pr_L ** (1 / 3)

    return Nu * k_L / D_h
