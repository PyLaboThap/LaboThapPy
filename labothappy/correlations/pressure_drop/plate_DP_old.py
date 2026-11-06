# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 09:55:59 2026

@author: Basile
"""

import numpy as np

#%% CONDENSATION

def han_BPHEX_DP(mu, G, Dh, beta, pitch_co, rho_v, rho_l, L_v, N_cp, m_dot, D_p): 
    """
    Inputs
    ------
    mu       : viscosity [Pa*s]
    rho_l    : Liquid density [kg/m^3]
    rho_v    : Vapor density [kg/m^3]
    G        : Mass flux [kg/(m^2 * s)]
    Dh       : Plate spacing [m]
    pitch_co : Corrugated pitch [m]
    beta     : Chevron angle [°]
    L_v      : Vertical length between fluid ports [m]
    N_cp     : Number of canals [-]
    m_dot    : Flowrate [kg/s]
    D_p      : Port diameter [m]
    
    Outputs
    -------
    h_cond : Condensation HTC [W/(m^2*K)]
    Nu     : Nusselt Number [-]
    DP_tot : Total Pressure Drop [Pa]
    
    Reference
    ---------
    "The caracteristics of condensation in brazed plate heat exchangers with different chevron angles", Journal of the Korean Physical Society, 2003
    Han et. al
    
    """
    
    # Preliminary calculations
    theta = np.pi/2 - beta
    g = 9.81 # gravity acceleration
    
    Re_eq = G*Dh/mu
    
    # Pressure drop
    Ge3 = 3521.1*(pitch_co/Dh)**(4.17)*(theta)**(-7.75)
    Ge4 = -1.024*(pitch_co/Dh)**(0.0925)*(theta)**(-1.3)
    
    f = Ge3*Re_eq**Ge4
    
    rho = rho_v + rho_l
    
    # Two phase related pressure drop
    DP_tp = f*(L_v*N_cp/Dh)*G**2*rho
    
    # Port pressure drop
    G_p = 4*(m_dot/(np.pi*D_p**2))

    DP_port = 1.4*G_p**2/(2*rho)

    # Static head loss
    DP_stat = -rho*g*L_v # negative because downward flow <-> Condenser

    # The acceleration pressure drop for condensation is expressed as : ??? 
    
    DP_tot = DP_tp + DP_port + DP_stat
    
    return DP_tot

def _fanning_tube(Re):
    """
    Smooth circular tube Fanning friction factor as used in ACHP's
    Lockhart-Martinelli implementation (laminar 16/Re, turbulent Blasius-type
    0.046 Re^-0.2, linear blend for 1000 < Re < 2000).
    """
    if Re < 1000:
        return 16 / Re
    elif Re > 2000:
        return 0.046 / Re ** 0.2
    else:
        w = (Re - 1000) / (2000 - 1000)
        return (1 - w) * 16 / Re + w * 0.046 / Re ** 0.2

def _zeta_martin(Re, phi):
    """Martin (1996) / VDI Heat Atlas N6 Darcy-type friction factor.
    phi: corrugation angle from the main flow direction [rad]."""
    Re = max(Re, 1.0)
    # Straight-channel (phi = 0) and cross-corrugated (phi = 90 deg) factors
    if Re < 2000:
        zeta0, zeta1_0 = 64 / Re, 597 / Re + 3.85
    else:
        zeta0 = (1.8 * np.log10(Re) - 1.5) ** -2     # Konakov; check vs. the source
        zeta1_0 = 39 / Re ** 0.289

    a, b, c = 3.8, 0.18, 0.36
    zeta1 = a * zeta1_0
    rhs = (np.cos(phi) / np.sqrt(b * np.tan(phi) + c * np.sin(phi) + zeta0 / np.cos(phi))
           + (1 - np.cos(phi)) / np.sqrt(zeta1))
    return 1.0 / rhs**2

def dp_lockhart_martinelli(G, x, rho_l, rho_v, mu_l, mu_v, D_h, L,
                                  C=None, chevron_angle=None):
    """
    Two-phase FRICTIONAL pressure drop with the Lockhart-Martinelli separated
    flow model, in the analytical form of Chisholm (1967).

    Original validity: adiabatic two-phase flow in horizontal, smooth circular
    tubes (Lockhart & Martinelli, 1949). Applying it to a plate heat exchanger
    is an extrapolation.

    Inputs
    ------
    G             : total mass flux in ONE channel/tube [kg/m^2/s]
    x             : vapour quality [-] (use the zone-average quality for a zone)
    rho_l, rho_v  : saturated liquid / vapour density [kg/m^3]
    mu_l, mu_v    : saturated liquid / vapour dynamic viscosity [Pa.s]
    D_h           : hydraulic diameter [m] (tube ID, or 2*b/enlargement_factor for a PHE)
    L             : flow length [m] (L_p, or w_i * L_p for a moving boundary zone)
    C             : Chisholm constant [-], optional. If None, chosen from the
                    laminar/turbulent regime of each phase (5, 10, 12 or 20)
    chevron_angle : corrugation angle FROM THE MAIN FLOW DIRECTION [rad], optional.
                    If None  -> smooth tube friction factors (original method).
                    If given -> Martin's PHE friction factor for each phase
                                (adaptation, C not refitted for plates).

    Output
    -------
    dp : frictional two-phase pressure drop over the length L [Pa]
         (acceleration, gravity and port losses NOT included)
    """
    x = float(x)

    # 1. Reynolds number of each phase flowing ALONE at its actual flow rate
    Re_l = G * (1 - x) * D_h / mu_l
    Re_v = G * x * D_h / mu_v

    # 2. Single-phase frictional pressure gradients [Pa/m]
    def dpdz_single(Re, G_phase, rho):
        if G_phase <= 0:
            return 0.0
        if chevron_angle is None:
            f_F = _fanning_tube(Re)                      # Fanning
            return 2 * f_F * G_phase**2 / (rho * D_h)
        else:
            zeta = _zeta_martin(Re, chevron_angle)       # Darcy
            return zeta * G_phase**2 / (2 * rho * D_h)

    dpdz_l = dpdz_single(Re_l, G * (1 - x), rho_l)
    dpdz_v = dpdz_single(Re_v, G * x, rho_v)

    # 3. Chisholm constant from the flow regime of each phase
    #    (1500 used as transition Re to ensure continuity, as in ACHP)
    if C is None:
        if Re_l > 1500 and Re_v > 1500:
            C = 20.0    # turbulent liquid - turbulent vapour
        elif Re_l < 1500 and Re_v > 1500:
            C = 12.0    # laminar liquid - turbulent vapour
        elif Re_l > 1500 and Re_v < 1500:
            C = 10.0    # turbulent liquid - laminar vapour
        else:
            C = 5.0     # laminar liquid - laminar vapour

    # 4. Two-phase gradient (Chisholm form, identical to phi_l^2*dpdz_l = phi_v^2*dpdz_v)
    dpdz = dpdz_l + C * np.sqrt(dpdz_l * dpdz_v) + dpdz_v

    dp = dpdz * L
    return dp