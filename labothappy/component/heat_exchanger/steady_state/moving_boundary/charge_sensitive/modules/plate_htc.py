# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:47:52 2023

@author: Basile
"""

import numpy as np
from scipy.optimize import fsolve

#%%
def water_plate_HTC(mu, Pr, k, G, Dh):
    """
    Calibrated heat transfer coeffcient correlation for water side
    
    Inputs
    ----------
    mu : Viscosity [kg/(m*s)]
    
    Pr : Prandtl Number [/]
    
    k : thermal conductivity [W/(m*K)]
        
    G : Mass flux [kg/(m^2 * s)]
    
    Dh : Spacing between plates [m]

    Outputs
    -------
    h_conv : HTC in convection
    
    Reference
    -------
    Refrigerant R134a vaporisation heat transfer and pressure drop inside a small brazed plate heat exchanger
    G.A. Longo, A. Gasparella

    """
    
    # Bounds on Re (arbitrary) # !!! shall be checked
    Re_max = 1e6
    Re_min = 5
    
    # Reynolds number
    Re = G*Dh/mu
    
    if Re <= Re_max and Re >= Re_min:
        # Nusselt number
        Nu = 0.277*Re**(0.766)*Pr**(0.333)
    else: 
        print("Reynolds number for water out of bounds.")
        return 0
        
    # HTC
    h_conv = (k/Dh)*Nu

    return h_conv

def simple_plate_HTC(mu, Pr, k, G, Dh):
    # Reynolds number
    Re = G*Dh/mu
    
    if Re < 5*1e5:
        Nu = 0.3387*Re**(1/2)*Pr**(1/3)/(1+(0.0468/Pr)**(2/3))**(1/4)
    else:
        Nu = 0.0296*Re**(4/5)*Pr**(1/3)
        
    h_conv = (k/Dh)*Nu
    
    return h_conv

def muley_manglik_BPHEX_HTC(mu, mu_w, Pr, k, G, Dh, chevron_angle):
    # Reynolds number
    Re = G*Dh/mu
    
    beta = 180*chevron_angle/np.pi
    
    C = 0.2668 - 0.006967*beta + 7.244*1e-5*beta**2
    C_2 = Re**(0.728 + 0.0543*np.sin((2*np.pi*beta/90) + 3.7))
    
    Nu = C * C_2 * Pr**(1/3) * (mu/mu_w)**(0.14)
        
    h_conv = (k/Dh)*Nu
    
    return h_conv

def martin_BPHEX_HTC(mu, mu_w, Pr, k, G, Dh, chevron_angle):
    "Martin Holger: Correlation from the VDI Atlas"
 
    beta = chevron_angle
    Re = G*Dh/mu
    
    "Factor for correlations: provided by Focke et al."
    if Re >= 2000: # Regime : Turbulent
        xhi_0   = (1.8*np.log(Re)-1.5)**-2
        xhi_1_0 = 39/Re**0.289
    elif Re < 2000: # Regime: Laminar
        xhi_0   = 64/Re
        xhi_1_0 = 597/Re +3.85
    
    "Constant given by Martin"
    a = 3.8
    b = 0.18 
    c = 0.36    
    
    "Factor xhi"
    xhi_1 = a*xhi_1_0
    
    "Friction factor"    
    f = (np.cos(beta)/np.sqrt(b*np.tan(beta) + c*np.sin(beta) + xhi_0/np.cos(beta)) +(1 - np.cos(beta))/np.sqrt(xhi_1))**(-2)
    
    "Hagen number"
    Hg = f*Re**2/2
    
    "Extracted from the comparison with Heavear et al. [10]"
    c_q = 0.122
    q = 0.374
    
    "Nusslet number:"
    Nu = c_q*Pr**(1/3)*(mu/mu_w)**(1/6)*(2*Hg*np.sin(2*beta))**q
    
    "Heat Transffer Coefficient [W m^-2]:"
    hcv = Nu*k/Dh
    
    return hcv 


#%%
def han_boiling_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, i_fg, G, DT_log, Qdot, hconv_h, Dh, theta, pitch_co):
    """
    Inputs
    ------
    x        : Vapor quality [-]
    mu_l     : Liquid viscosity [Pa*s]
    k_l      : Liquid thermal conductivity [W/(m*K)]
    Pr_l     : Liquid Prandtl Number [-]
    rho_l    : Liquid density [kg/m^3]
    rho_v    : Vapor density [kg/m^3]
    i_fg     : Vaporisation heat [J/(kg)]
    G        : Mass flux [kg/(m^2 * s)]
    DT_log   : Logarithmic temperature difference [K]    
    Qdot     : Heat Transfer rate [W]
    hconv_h  : Convection heat transfer coefficient [W/(m^2*K)]
    Dh       : Plate spacing [m]
    theta    : Chevron angle [°]
    pitch_co : Corrugated pitch [m]
    
    Outputs
    -------
    h_cond : Condensation HTC [W/(m^2*K)]
    Nu     : Nusselt Number [-]
    DP_tot : Total Pressure Drop [Pa]
    
    Reference
    ---------
    HFC-410A vaporisation inside a commercial brazed plate heat exchanger
    Giovanni A. Longo, Andrea Gasparella
    """
    
    def iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        Bo_g = max(Bo_g, 1e-8)
        Nu = Ge1*Re_eq**Ge2*Bo_g**0.3*Pr_l**0.4
        h = Nu*k_l/Dh
        U = (1/h +  1/hconv_h)**-1
        A_tp = AU_tp/U
        q = Qdot/A_tp
        Bo = q/(G_eq*i_fg)
        res_Bo = (Bo-Bo_g)/Bo_g
        return res_Bo, Nu, h, U, A_tp, q, Bo
    
    def res_iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        res_Bo,_,_,_,_,_,_ = iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
        
        return res_Bo
    
    G_eq = G * ( (1 - x) + x * (rho_l/rho_v)**0.5)
    Re_eq = G_eq*Dh/mu_l
    AU_tp = Qdot/DT_log
    Ge1 = 2.81*(pitch_co/Dh)**(-0.041)*(theta)**(-2.83)
    Ge2 = 0.746*(pitch_co/Dh)**(-0.082)*(theta)**(0.61)
    Bo_0 = 0.5
    f_Bo = lambda xx: res_iter_Han_boiling(xx, Ge1, Ge2,Re_eq, Pr_l, k_l, hconv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    sol = fsolve(f_Bo, Bo_0,)
    Bo_sol = sol[0]
    _, Nu, h, _, _, _, _ = iter_Han_boiling(Bo_sol, Ge1, Ge2,Re_eq, Pr_l, k_l, hconv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    
    h_boiling = h
    
    return h_boiling, Nu

def han_cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, G, Dh, pitch_co, beta, L_v, N_cp, m_dot, D_p):
    """
    Inputs
    ------
    x        : Vapor quality [-]
    mu_l     : Liquid viscosity [Pa*s]
    k_l      : Liquid thermal conductivity [W/(m*K)]
    Pr_l     : Liquid Prandtl Number [-]
    rho_l    : Liquid density [kg/m^3]
    rho_v    : Vapor density [kg/m^3]
    G        : Mass flux [kg/(m^2 * s)]
    Dh       : Plate spacing [m]
    pitch_co : Corrugated pitch [m]
    beta     : Chevron angle [°]
    L_v      : Vertical length between fluid ports [m] 
    N_cp     : Number of canals [-]
    m_dot    : Flowrate [kg/s]
    
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
    
    G_eq = G * ( (1-x) + x * (rho_l/rho_v)**0.5)
    Re_eq = G_eq*Dh/mu_l
    
    # Heat Transfer
    Ge1 = 11.22*(pitch_co/Dh)**(-2.83)*(theta)**(-4.5)
    Ge2 = 0.35*(pitch_co/Dh)**(0.23)*(theta)**(1.48)
    
    Nu = Ge1*Re_eq**Ge2*Pr_l**(1/3)
    h_cond = Nu*k_l/Dh
    
    # Pressure drop
    Ge3 = 3521.1*(pitch_co/Dh)**(4.17)*(theta)**(-7.75)
    Ge4 = -1.024*(pitch_co/Dh)**(0.0925)*(theta)**(-1.3)
    
    f = Ge3*Re_eq**Ge4
    
    # Two phase related pressure drop
    DP_tp = f*(L_v*N_cp/Dh)*G_eq**2*rho_l
    
    # Port pressure drop
    m_dot_eq = m_dot*(1 - x + x*(rho_l/rho_v)**0.5)
    G_p = 4*(m_dot_eq/(np.pi*D_p**2))
    rho_m = 1/( (x/rho_v) + (1 - x)/rho_l )

    DP_port = 1.4*G_p**2/(2*rho_m)

    # Static head loss
    DP_stat = -rho_m*g*L_v # negative because downward flow <-> Condenser

    # The acceleration pressure drop for condensation is expressed as : ??? 
    
    DP_tot = DP_tp + DP_port + DP_stat
    
    return h_cond, Nu, DP_tot

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