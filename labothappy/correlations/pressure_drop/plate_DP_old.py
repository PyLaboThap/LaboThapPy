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

