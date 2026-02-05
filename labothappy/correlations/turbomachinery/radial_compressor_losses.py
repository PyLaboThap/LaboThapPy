# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 09:04:20 2025

@author: Basile
"""

import numpy as np

def rotor_incidence_losses(A1, A1_th, beta1, w1, xhi1):
    """
    A1 : Rotor inlet area [m2] (from blade pitch * height)
    A1_th : Rotor inlet throat area [m2] (from inlet throat opening * height)
    beta1 : Rotor Relative flow angle [rad]
    w1 : Rotor inlet relative velocity [m/s]
    xhi1 : Rotor Inlet Blade Angle [rad]

    Reference
    ---------
    Fortran program for predicting off-design performance of centrifugal
    compressors. Tech Rep (NASA); 1973. - Galvas M. 
    """    

    # Optimal incidence angle
    beta1_opt = np.arctan((A1/A1_th)*xhi1)
    
    # Losses
    Dh_inc = 0.5* w1**2 * np.sin(beta1-beta1_opt)**2
    
    return Dh_inc

def rotor_friciton_losses(beta1h, beta1s, b2, C_fi, L_z, n_bl_r, r1h, r1s, r2, w1, w1_th, w2, xhi2):
    """
    beta1h : Rotor inlet hub angle [rad]
    beta1s : Rotor inlet shroud angle [rad]
    b2     : Rotor outlet passage height [m]
    C_fi   : Rotor friction coefficient [-] : 0.004 ? 
    L_z    : Rotor axial length [m]
    n_bl_r : Rotor blade number [-]
    r1h    : Rotor inlet hub radius [m]
    r1s    : Rotor inlet shroud radius [m]
    r2     : Rotor outlet radius [m]
    w1     : Rotor inlet relative velocity [m/s]
    w1_th  : Rotor inlet throat relative velocity [m/s]
    w2     : Rotor outlet relative velocity [m/s]
    xhi2   : Rotor Outlet Blade Angle [rad]
    
    References
    ---------
    Aungier RH. Mean streamline aerodynamic performance analysis of centrifugal
    compressors. J Turbomach 1995;117(3):360.
    
    Fortran program for predicting off-design performance of centrifugal
    compressors. Tech Rep (NASA); 1973. - Galvas M. 
    
    Jansen W. A method for calculating the flow in a centrifugal impeller when entropy
    gradients are present. In: Royal Society conference on internal aerodynamics 
    (turbomachinery), 19-21 July, Cambridge, UK; 1967. p. 133–46.
    """    
    
    # Max mean velocity
    w_bar = np.max([np.sqrt( (w1**2  + w2**2)/2 ), np.sqrt( (w1_th**2  + w2**2)/2 )])    
    
    # Blade length
    fact1 = 2*r2 - (r1s+r1h) - b2 + 2*L_z
    fact2 = 2/((np.cos(beta1s)+np.cos(beta1h))/2 + np.cos(xhi2))
    
    L_b = (np.pi/8)*fact1*(fact2)
    
    # Hydraulic diameter
    term1 = 2*r2 / (n_bl_r/(np.pi*np.cos(xhi2)) + 2*r2/b2)
    term2 = 2*r1s / (2/(1 - r1h/r1s) + 2*n_bl_r / (np.pi*(1+r1h/r1s)) * np.sqrt(1+np.tan(beta1s)**2*(1 + (r1h/r1s)**2 / 2)))
    
    D_hb = term1 + term2
    
    # Losses
    Dh_f = 4* C_fi * (L_b/D_hb) * w_bar**2
    
    return Dh_f

def rotor_blade_loading_losses(C_df, Dh0, n_bl_r, r1s, r2, u2, w1, w1s, w2):
    """
    C_df   : Rotor friction coefficient [-] : 0.004 ? 
    Dh0    : Stage specific enthalpy difference [J/kg]
    n_bl_r : Rotor blade number [-]
    r1s    : Rotor inlet shroud radius [m]
    r2     : Rotor outlet radius [m]
    u2     : Rotor outlet tangential velocity [m/s]
    w1     : Rotor inlet relative velocity [m/s]
    w1s    : Rotor inlet shroud relative velocity [m/s]
    w2     : Rotor outlet relative velocity [m/s]
    
    References
    ---------
    Coppage J, Dallenbach F. Study of supersonic radial compressors for refrigeration
    and pressurization systems. Tech Rep, Garret Corp; 1956.
    """    
    
    # Blade Diffusion Factor 
    D_f = 1- w2/w1 + C_df*(Dh0/u2**2)*(w2/w1s)/(n_bl_r/np.pi * (1-r1s/r2) + 2*r1s/r2)
    
    # Losses
    Dh_bl = 0.05*D_f**2*u2**2
    
    return Dh_bl

def rotor_clearance_losses(b2, eps_a, eps_r, n_bl_r, rho1, rho2, r1h, r1s, r2, vu2, v1m):
    """
    b2     : Rotor outlet passage height [m]
    eps_a  : Rotor axial clearance [m] 
    eps_r  : Rotor radial clearance [m] 
    n_bl_r : Rotor blade number [-]
    rho1   : Rotor inlet density [kg/m^3]
    rho2   : Rotor outlet density [kg/m^3]
    r1h    : Rotor inlet hub radius [m]
    r1s    : Rotor inlet shroud radius [m]
    r2     : Rotor outlet radius [m]
    vu2    : Rotor outlet absolute tangential velocity [m/s]
    v1m    : Rotor inlet absolute meridional velocity [m/s]
    
    References
    ---------
    Jansen W. A method for calculating the flow in a centrifugal impeller when entropy
    gradients are present. In: Royal Society conference on internal aerodynamics 
    (turbomachinery), 19-21 July, Cambridge, UK; 1967. p. 133–46.
    """    
    
    eps_mean = (eps_a + eps_r)/2
    fact1 = 0.6*eps_mean/b2 * np.abs(vu2)
    fact2 = (4*np.pi/(b2*n_bl_r)) * ((r1s**2 - r1h**2)/((r2-r1s) * (1+rho2/rho1))) * np.abs(vu2) * v1m
    
    # Losses
    Dh_cl = fact1 * np.sqrt(fact2)
    
    return Dh_cl

def rotor_mixing_losses(alpha2, v2):
    """
    alpha2 : Rotor outlet absolute angle [rad]
    v2     : Rotor outlet absolute velocity [m/s]
    
    References
    ---------
    Jansen W. A method for calculating the flow in a centrifugal impeller when entropy
    gradients are present. In: Royal Society conference on internal aerodynamics 
    (turbomachinery), 19-21 July, Cambridge, UK; 1967. p. 133–46.
    
    Japikse D. Centrifugal compressor design and performance. Wilder, VT, USA:
    Concepts ETI Inc.; 1996.
    """    
    
    # Wake fraction at impeller exit
    # xhi = 0.93 * eps_w**2 + 0.07*eps_w where xhi = 0.15 and b_star = 1 : solving 2nd degree equation for eps_w
    b_star = 1
    xhi = 0.15
    eps_w = (-0.07 + np.sqrt(0.07**2 + 4*0.93*xhi))/2*0.93
    
    fact1 = 1/(1+np.tan(alpha2)**2)
    fact2 = (1-eps_w-b_star/(1-eps_w))
    
    # Losses
    Dh_ml = 0.5*v2**2 * fact1 * fact2**2
    
    return Dh_ml

def rotor_disc_friction_losses(b2, eps_b, mdot, mu2, rho1, rho2, r2, u2):
    """
    b2    : Rotor outlet passage height [m]
    eps_b : Back face clearance gap [m]
    mdot  : Flow rate [kg/s]
    mu2   : Rotor outlet viscosity [Pa*s]
    rho1  : Rotor inlet density [kg/m^3]
    rho2  : Rotor outlet density [kg/m^3]
    r2    : Rotor outlet radius [m]
    u2    : Rotor outlet tangential velocity [m/s]
    
    References
    ---------
    Jansen W. A method for calculating the flow in a centrifugal impeller when entropy
    gradients are present. In: Royal Society conference on internal aerodynamics 
    (turbomachinery), 19-21 July, Cambridge, UK; 1967. p. 133–46.
    """
    # Reynolds number
    Re2 = rho2*u2*r2/mu2
    
    # Total pressure loss coefficient
    if Re2 < 3*1e5:
        Kf = 3.7*(eps_b/b2)**0.1 / Re2**0.5
    else:
        Kf = 0.102*(eps_b/b2)**0.1 / Re2**0.2  
    
    # Mean density
    rho_mean = (rho1+rho2)/2
    
    # Losses
    Dh_df = 0.25*(rho_mean*u2**3*r2**2*Kf)/mdot
    
    return Dh_df

def rotor_recirculation_losses(alpha2, C_df, Dh0, n_bl_r, r1s, r2, u2, w1, w1s, w2):
    """
    alpha2 : Rotor outlet absolute angle [rad]
    C_df   : Rotor friction coefficient [-] : 0.004 ? 
    Dh0    : Stage specific enthalpy difference [J/kg]
    n_bl_r : Rotor blade number [-]
    r1s    : Rotor inlet shroud radius [m]
    r2     : Rotor outlet radius [m]
    u2     : Rotor outlet tangential velocity [m/s]
    w1     : Rotor inlet relative velocity [m/s]
    w1s    : Rotor inlet shroud relative velocity [m/s]
    w2     : Rotor outlet relative velocity [m/s]
    
    References
    ---------
    Oh HW, Yoon ES, Chung MK. An optimum set of loss models for performance
    prediction of centrifugal compressors. Proc Inst Mech Eng, Part A: J Power Energy
    1997;211(4):331–8.
    """    
    
    # Blade Diffusion Factor 
    D_f = 1- w2/w1 + C_df*(Dh0/u2**2)*(w2/w1s)/(n_bl_r/np.pi * (1-r1s/r2) + 2*r1s/r2)
    
    # Losses
    Dh_rc = 8*1e-5*np.sinh(3.5*alpha2**3)*D_f**2*u2**2
    
    return Dh_rc

#%%

def stator_incidence_losses(A4, A4_th, beta4, w4, xhi4):
    """
    A4    : Stator inlet area [m2] (from blade pitch * height)
    A4_th : Stator inlet throat area [m2] (from inlet throat opening * height)
    beta4 : Stator Relative flow angle [rad]
    w4    : Stator inlet relative velocity [m/s]
    xhi4  : Stator Inlet Blade Angle [rad]

     
    Reference
    ---------
    Conrad O, Raif K, Of MW. The calculation of performance maps for centrifugal
    compressors with vane-island diffusers. Proceedings of the twenty-fifth annual in
    ternational gas turbine conference and exhibit and twenty-second annual fluids
    engineering conference, March 9-13. New Orleans, Louisiana, USA: American
    Society of Mechanical Engineers; 1980. p. 135–47.
    """    

    # Optimal incidence angle
    beta4_opt = np.atan((A4/A4_th)*xhi4)
    
    # Losses
    Dh_inc = 0.5 * w4**2 * np.sin(beta4_opt-beta4)
    
    return Dh_inc

def stator_friction_losses(C_f, r3, r5, vm, xhi3, xhi5):
    """
    C_f  : Stator friction coefficient [-] (=0.004 ?) : C_f =  k(1.8*1e5/ Re)**0.2 with k = 0.005
    r3   : Stator inlet diameter [m]
    r5   : Stator outlet diameter [m]
    vm   : Stator meridional velocity [m/s]
    xhi3 : Stator Inlet Blade Angle [rad]
    xhi5 : Stator Outlet Blade Angle [rad]

    Reference
    ---------
    Conrad O, Raif K, Of MW. The calculation of performance maps for centrifugal
    compressors with vane-island diffusers. Proceedings of the twenty-fifth annual in
    ternational gas turbine conference and exhibit and twenty-second annual fluids
    engineering conference, March 9-13. New Orleans, Louisiana, USA: American
    Society of Mechanical Engineers; 1980. p. 135–47.
    """    

    L_b_over_D_h = (r5-r3)/np.cos( (xhi3 + xhi5)/2 )
    
    # Losses
    Dh_f = 2 * C_f * vm**2 * L_b_over_D_h
    
    return Dh_f

#%%

def radial_compressor_rotor_losses(A1, A1_th, alpha2, beta1, beta1h, beta1s, b2, C_df, C_fi, Dh0, eps_a, eps_b, eps_r, L_z, mdot, mu2, n_bl_r,
                                   rho1, rho2, r1h, r1s, r2, u2, vu2, v1m, v2, w1, w1_th, w1s, w2, xhi1, xhi2):
    """
    A1 : Rotor inlet area [m2] (from blade pitch * height)
    A1_th : Rotor inlet throat area [m2] (from inlet throat opening * height)
    alpha2 : Rotor outlet absolute angle [rad]
    beta1 : Rotor Relative flow angle [rad]
    beta1h : Rotor inlet hub angle [rad]
    beta1s : Rotor inlet shroud angle [rad]
    b2     : Rotor outlet passage height [m]
    C_df   : Rotor friction coefficient [-] : 0.004 ? 
    C_fi   : Rotor friction coefficient [-] : 0.004 ? 
    Dh0    : Stage specific enthalpy difference [J/kg]
    eps_a  : Rotor axial clearance [m] 
    eps_b : Back face clearance gap [m]
    eps_r  : Rotor radial clearance [m] 
    L_z    : Rotor axial length [m]
    mdot  : Flow rate [kg/s]
    mu2   : Rotor outlet viscosity [Pa*s]
    n_bl_r : Rotor blade number [-]
    rho1   : Rotor inlet density [kg/m^3]
    rho2   : Rotor outlet density [kg/m^3]
    r1h    : Rotor inlet hub radius [m]
    r1s    : Rotor inlet shroud radius [m]
    r2     : Rotor outlet radius [m]
    u2     : Rotor outlet tangential velocity [m/s]
    vu2    : Rotor outlet absolute tangential velocity [m/s]
    v1m    : Rotor inlet absolute meridional velocity [m/s]
    v2     : Rotor outlet absolute velocity [m/s]
    w1     : Rotor inlet relative velocity [m/s]
    w1_th  : Rotor inlet throat relative velocity [m/s]
    w1s    : Rotor inlet shroud relative velocity [m/s]
    w2     : Rotor outlet relative velocity [m/s]
    xhi1 : Rotor Inlet Blade Angle [rad]
    xhi2   : Rotor Outlet Blade Angle [rad]
    """
    
    Dh_rot = {}
    
    Dh_rot['inc'] = rotor_incidence_losses(A1, A1_th, beta1, w1, xhi1)
    
    Dh_rot['f'] = rotor_friciton_losses(beta1h, beta1s, b2, C_fi, L_z, n_bl_r, r1h, r1s, r2, w1, w1_th, w2, xhi2)
    
    Dh_rot['bl'] = rotor_blade_loading_losses(C_df, Dh0, n_bl_r, r1s, r2, u2, w1, w1s, w2)
    
    Dh_rot['cl'] = rotor_clearance_losses(b2, eps_a, eps_r, n_bl_r, rho1, rho2, r1h, r1s, r2, vu2, v1m)
    
    Dh_rot['ml'] = rotor_mixing_losses(alpha2, v2)
    
    Dh_rot['df'] = rotor_disc_friction_losses(b2, eps_b, mdot, mu2, rho1, rho2, r2, u2)
    
    Dh_rot['rc'] = rotor_recirculation_losses(alpha2, C_df, Dh0, n_bl_r, r1s, r2, u2, w1, w1s, w2)
    
    Dh_rot['tot'] = Dh_rot['inc'] + Dh_rot['f'] + Dh_rot['bl'] + Dh_rot['cl'] + Dh_rot['ml'] + Dh_rot['df'] + Dh_rot['rc']
    
    return Dh_rot

def radial_compressor_stator_losses():
    
    
    return

