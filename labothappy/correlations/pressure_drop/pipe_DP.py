# -*- coding: utf-8 -*-
"""
@author: Basile Chaudoir
"""

from math import log10, inf
import numpy as np
import warnings
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

#%%
def gnielinski_pipe_DP(mu, rho, G, Dh, L):
    """
    Inputs
    ------
    
    mu   : Dynamic Viscosity [Pa*s]
    Pr   : Prantl number [-]
    Pr_w : Prandtl number at wall conditions [-]
    k    : Thermal conductivity [W/(m*K)]
    G    : Flow rate per cross section area [kg/(m^2 * s)]
    Dh   : Hydraulic diameter [m]
    L    : Flow length [m]
    
    Outputs
    -------
    
    DP : Pressure drop [Pa]
    
    Reference
    ---------
    Validation of the Gnielinski correlation for evaluation of heat transfer coefficient of enhanced tubes by non-linear 
    regression model: An experimental study of absorption refrigeration system
    
    Syed Muhammad Ammar, Chan Woo Park
    """
    #-------------------------------------------------------------------------
    Re_min = 0
    Re_max = 1e06
    Re = G*Dh/mu

    #-------------------------------------------------------------------------
    f = (1.8*log10(Re) - 1.5)**(-2)
    v_flow = G/rho
    
    DP = f*L*rho*v_flow**2/(2*Dh)    
    #-------------------------------------------------------------------------
    if Re >= Re_max or Re <=Re_min:
        # warnings.warn('Gnielinski singe-phase: Out of validity range --> Re = ', Re, ' is out of [', Re_min, ' - ', Re_max, '] !!!')
        warnings.warn('Gnielinski singe-phase: Reynolds Out of validity range !!!')
    #-------------------------------------------------------------------------
    # print(f"Re : {Re}")
    # print(f"f : {f}")
    # print(f"v_flow : {v_flow}")
    # print(f"DP : {DP}")
    # print(f"rho : {rho}")
    
    return DP

def Muller_Steinhagen_Heck_DP(fluid, G, P_sat, x, L, D_in):
    """
    For 1 phase flow
    
    Inputs
    ------
    
    mu   : Dynamic Viscosity [Pa*s]
    rho  : Density [kg/m^3]
    G    : Flow rate per cross section area [kg/(m^2 * s)]
    Dh   : Hydraulic diameter [m]
    L    : Flow length [m]
    
    Outputs
    -------
    
    DP : Pressure drop [Pa]
    
    Reference
    ---------
    GENERALIZED PRESSURE DROP CORRELATIONFOR EVAPORATION AND CONDENSATIONIN SMOOTH AND MICRO-FIN TUBES
    
    Choi, Kedzierski, Domanski
    """    
    mu_l, h_l, rho_l = PropsSI(('V','H','D'), 'P', P_sat, 'Q', 0, fluid)
    mu_v, h_v, rho_v = PropsSI(('V','H','D'), 'P', P_sat, 'Q', 1, fluid)
        
    Re_l = G*D_in/mu_l
    Re_v = G*D_in/mu_v
    
    f_l = 0.079/Re_l**0.25
    f_v = 0.079/Re_v**0.25
    
    A = f_l*2*G**2 / (D_in*rho_l)
    B = f_v*2*G**2 / (D_in*rho_v)
    
    Y = A + 2*(B-A)*x
    
    DP_dz = Y*(1-x)**(1/3) + B*x**3

    DP = DP_dz*L

    return DP

def Choi_DP(fluid, G, rho_out, rho_in, P_sat, x_o, x_i, L, D_in):
    """
    For 1 phase flow
    
    Inputs
    ------
    
    mu   : Dynamic Viscosity [Pa*s]
    rho  : Density [kg/m^3]
    G    : Flow rate per cross section area [kg/(m^2 * s)]
    Dh   : Hydraulic diameter [m]
    L    : Flow length [m]
    
    Outputs
    -------
    
    DP : Pressure drop [Pa]
    
    Reference
    ---------
    GENERALIZED PRESSURE DROP CORRELATIONFOR EVAPORATION AND CONDENSATIONIN SMOOTH AND MICRO-FIN TUBES
    
    Choi, Kedzierski, Domanski
    """
    g = 9.81
    
    mu_l, h_l = PropsSI(('V','H'), 'P', P_sat, 'Q', 0, fluid)
    h_v = PropsSI('H', 'P', P_sat, 'Q', 1, fluid)
    
    Dh_lv = h_v - h_l
    Re_f = G*D_in/mu_l
    
    K_f = (x_o - x_i)*Dh_lv/(g*L)
    
    f_N = 0.00506*Re_f**(-0.0951) * K_f**(0.1554)

    DP = (f_N*L*(1/rho_out + 1/rho_in)/D_in + (1/rho_out + 1/rho_in))*G**2 

    return DP


def Churchill_DP(mu, rho, G, Dh, L):
    """
    For 1 phase flow
    
    Inputs
    ------
    
    mu   : Dynamic Viscosity [Pa*s]
    rho  : Density [kg/m^3]
    G    : Flow rate per cross section area [kg/(m^2 * s)]
    Dh   : Hydraulic diameter [m]
    L    : Flow length [m]
    
    Outputs
    -------
    
    DP : Pressure drop [Pa]
    
    Reference
    ---------
    Pressure Drop Model For Condensation From Superheated Vapor
    
    Jiange Xiao
    """
    
    e = 0.045*1e-3 # m : Roughness of Steel Pipes
    Re = G*Dh/mu # Reynolds Numbers

    J_1 = (-2.457*np.log((7/Re)**0.9 + 0.27*(e/Dh)))**16
    J_2 = (37530/Re)**16

    f = 8*((8/Re)**12 + (1/J_1 + 1/J_2)**1.5)**(1/2)

    DP = f*L*G**2/(2*rho*Dh)    

    return DP

def pipe_internal_DP(mu, rho, m_dot, params):
    """
    Inputs
    ------
    
    mu   : Dynamic Viscosity [Pa*s]
    Np   : Number of tube passes [-]
    rho  : Density [kg/m^3]
    G    : Flow rate per cross section area [kg/(m^2 * s)]
    Dh   : Hydraulic diameter [m]
    L    : Flow length [m]
    
    Outputs
    -------
    
    DP : Pipe pressure drop [Pa]
    
    Reference
    ---------
    ?
    """
    
    A_in_1_tube = np.pi*(params["Tube_OD"] - 2*params["Tube_t"])**2 / 4
    A_in_tubes = A_in_1_tube*params["n_tubes"]

    G = m_dot/A_in_tubes

    # Reynolds number
    Re = G*(params["Tube_OD"] - 2*params["Tube_t"])/mu
    
    # Flow speed
    u = G/rho
    
    # Friction coefficient
    f = (1.8*log10(Re) - 1.5)**(-2)
    
    # Pressure drop (Np : number of passes)
    DP = (4*f*params["Tube_L"]*params["Tube_pass"]/params["Tube_OD"] + 4*params["Tube_pass"])*rho*(u**2)/2
    
    return DP
