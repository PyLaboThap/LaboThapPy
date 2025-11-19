# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 09:04:20 2025

@author: Basile
"""

import numpy as np

def nozzle_losses(v3, Re_3, alpha3, S3, c, b3):
    """
    v3 : Stator outlet speed [m/s]
    Re_3 : Chord stator Reynolds number [-]
    alpha3 : Stator outlet angle [rad]
    S3 : Stator pitch [m]
    c : Stator chord length [m]
    b3 : Stator passage width [m]
    
    Rodgers (1967)
    """    
    Dh_n = (v3**2 / 2) * (0.05/Re_3**0.2) * ((3*np.tan(alpha3))/(S3/c) + (S3*np.cos(alpha3))/b3) 
    
    return Dh_n

def rotor_incidence_losses(alpha4, n_blade_R, w4, beta4):
    """
    alpha4 : Rotor absolute inlet flow angle [rad]
    n_blade_R : Rotor blade count [-]
    w4 : Rotor relative inlet flow velocity [m/s]
    beta4 : Rotor relative inlet flow angle [rad]
    
    Wasserbauer & Glassman (1975)
    """    
    
    beta4_opt = np.arctan(np.tan(alpha4)/(1-(n_blade_R/1.98)))
    Dh_i = (w4**2 /2)*(np.sin(beta4 - beta4_opt))**2
    
    return Dh_i

def rotor_passage_losses(beta4, beta5, L_z, r4, r5h, r5s, r5, b4, b5, n_blade_R, w4, w5):
    """
    beta4 : Rotor relative inlet flow angle [rad]
    beta5 : Rotor relative outlet flow angle [rad]
    L_z : Rotor axial length [m]
    r4 : Rotor inlet radius [m]
    r5h : Rotor outlet hub radius [m]
    r5s : Rotor outlet shroud radius [m]
    r5 : Rotor outlet mean radius [m]
    b4 : Rotor inlet passage width [m]
    b5 : Rotor outlet passage width [m]
    n_blade_R : Rotor blade count [-]
    w4 : Rotor relative inlet flow velocity [m/s]    
    w5 : Rotor relative outlet flow velocity [m/s]    
        
    Takes into account both secondary flow and friction losses 
    Moustapha et al. (2003)
    """    
    
    K_p = 0.11 # Based on experimental data
    
    # Mean flow turning and c factor
    beta_mean = np.arctan((np.tan(beta4) + np.tan(beta5))/2)
    c = L_z/np.cos(beta_mean)
    
    # Hydraulic length in the rotor
    L_h = (np.pi/4)*((L_z-b4/2) + (r4 - r5s - b5/2))
    
    # Mean hydraulic diameter 
    Dh4 = 4*np.pi*r4**b4/(2*np.pi*r4 + n_blade_R*b4)
    Dh5 = 2*np.pi*(r5s**2 - r5h**2)/(np.pi*(r5s-r5h) + n_blade_R*b5)
    
    Dhmean = (Dh4+Dh5)/2

    # Correction factor to account for high secondary losses in high specific speed turbines with sharp meridional curvature
    if (r4-r5s)/b5 > 0.2:
        mf = 1
    elif (r4-r5s)/b5 <= 0.2:
        mf = 2
    
    factor = mf*K_p*(w4**2+w5**2)/2
    
    Dh_p_friction = factor*(L_h/Dhmean)
    Dh_p_secondary = factor*(0.68*(1-(r5/r4)**2))*np.cos(beta5)/(b5/c)
    
    return Dh_p_friction + Dh_p_secondary

def rotor_clearance_losses(cl_a, cl_r, L_z, r4, r5h, r5s, r5, b4, b5, n_blade_R, u4, vm4, vm5):
    """
    cl_a : Axial clearance [m]
    cl_r : Radial clearance [m]
    L_z : Rotor axial length [m]
    r4 : Rotor inlet radius [m]
    r5h : Rotor outlet hub radius [m]
    r5s : Rotor outlet shroud radius [m]
    r5 : Rotor outlet mean radius [m]
    b4 : Rotor inlet passage width [m]
    b5 : Rotor outlet passage width [m]
    n_blade_R : Rotor blade count [-]
    u4 : Rotor inlet blade speed [-]
    vm4 : Rotor inlet meridional speed [-]
    vm5 : Rotor outlet meridional speed [-]
    
    Moustapha et al. (2003)
    """
    
    K_a = 0.4 # Axial related coef
    K_r = 0.75 # Radial related coef
    K_ar = -0.3 # Axial-Radial interaction related coef
    
    C_a = (1-r5s/r4)/(vm4*b4)
    C_r = (r5s/r4)*(L_z-b4)/(vm5*r5*b5)
    
    Dh_cl = (n_blade_R*u4**3/(8*np.pi))*(K_a*cl_a*C_a + K_r*cl_r*C_r + K_ar*np.sqrt(cl_a*cl_r*C_a*C_r))
    
    if np.isnan(Dh_cl):
        return 0
    
    return Dh_cl

def rotor_trailing_edge_losses(w5, n_blade_R, t5, r5h, r5s, beta5, gamma5, Ma5_pr):
    """
    w5 : Rotor relative outlet flow velocity [m/s]    
    n_blade_R : Rotor blade count [-]
    t5 : Rotor outlet blade thickness [m]
    r5h : Rotor outlet hub radius [m]
    r5s : Rotor outlet shroud radius [m]
    beta5 : Rotor outlet relative flow angle [rad]
    gamma5 : Rotor outlet adiabatic coefficient [-]
    Ma5_pr : Rotor outlet relative Mach number [-]   
    
    Meroni (2018)
    """
    
    Dh_TE = (w5**2/2)*((n_blade_R*t5)/(np.pi*(r5s + r5h)*np.cos(beta5)))*(1+((gamma5-1)*Ma5_pr**2)/2)**(gamma5/(1-gamma5))
    
    return Dh_TE

def rotor_losses(alpha4, beta4, beta5, b4, b5, cl_a, cl_r, gamma5, L_z, Ma5_pr, n_blade_R, r4, r5, r5h, r5s, t5, u4, vm4, vm5, w4, w5):
    """
    alpha4 : Rotor absolute inlet flow angle [rad]
    beta4 : Rotor relative inlet flow angle [rad]
    beta5 : Rotor relative outlet flow angle [rad]
    b4 : Rotor inlet passage width [m]
    b5 : Rotor outlet passage width [m]
    cl_a : Axial clearance [m]
    cl_r : Radial clearance [m]
    gamma5 : Rotor outlet adiabatic coefficient [-]
    L_z : Rotor axial length [m]
    Ma5_pr : Rotor outlet relative Mach number [-]   
    n_blade_R : Rotor blade count [-]
    r4 : Rotor inlet radius [m]
    r5 : Rotor outlet mean radius [m]
    r5h : Rotor outlet hub radius [m]
    r5s : Rotor outlet shroud radius [m]
    t5 : Rotor outlet blade thickness [m]
    u4 : Rotor inlet blade speed [-]
    vm4 : Rotor inlet meridional speed [-]
    vm5 : Rotor outlet meridional speed [-]
    w4 : Rotor relative inlet flow velocity [m/s]    
    w5 : Rotor relative outlet flow velocity [m/s]  
    
    Aqel Thesis (2023)
    """
    
    # Incidence losses
    Dh_inc = rotor_incidence_losses(alpha4, n_blade_R, w4, beta4)
    
    # Friction and secondary losses
    Dh_p = rotor_passage_losses(beta4, beta5, L_z, r4, r5h, r5s, r5, b4, b5, n_blade_R, w4, w5)
    
    # Clearance losses
    Dh_cl = rotor_clearance_losses(cl_a, cl_r, L_z, r4, r5h, r5s, r5, b4, b5, n_blade_R, u4, vm4, vm5)
    
    # Trailing edge losses
    Dh_TE = rotor_trailing_edge_losses(w5, n_blade_R, t5, r5h, r5s, beta5, gamma5, Ma5_pr)
    
    # Total rotor losses
    Dh_tot = Dh_inc + Dh_p + Dh_cl + Dh_TE
    
    Dh_dict = {
        'Dh_tot' : Dh_tot,
        'Dh_inc' : Dh_inc,
        'Dh_p' : Dh_p,
        'Dh_cl' : Dh_cl,
        'Dh_TE' : Dh_TE,        
        }
    
    return Dh_dict


def windage_losses(omega, rho4, r4):
    """
    Or disk friction losses : It does not increase entropy in the flow it is more of a mechanical loss that reduces the produced power
    
    Daily & Nece (1960)
    """
    
    
    return

