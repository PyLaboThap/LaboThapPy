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
    For  1 phase flow
    
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

from CoolProp.CoolProp import PropsSI

def Muller_Steinhagen_Heck_DP(fluid, G, P_sat, x, L, D_in):
    """
    Pressure drop correlation for two-phase flow (condensation & evaporation)
    in smooth tubes, based on the Müller-Steinhagen and Heck model.

    Inputs:
    -------
    fluid : Fluid name (e.g., 'R134a')
    G : Mass flux [kg/(m^2*s)]
    P_sat : Saturation pressure [Pa]
    x : Vapor quality (mass fraction, 0 <= x <= 1)
    L : Tube length [m]
    D_in : Hydraulic diameter [m]

    Returns:
    --------
    DP : Pressure drop [Pa]

    Reference:
    ----------
    Müller-Steinhagen, H., & Heck, K. (1986). A general correlation for pressure drops in two-phase flow in pipes. 
    Int. J. Heat Mass Transfer, 29(3), 381–388.
    """    
    # Get liquid and vapor properties at saturation
    mu_l = PropsSI('V', 'P', P_sat, 'Q', 0, fluid)  # Dynamic viscosity [Pa·s]
    rho_l = PropsSI('D', 'P', P_sat, 'Q', 0, fluid)  # Density [kg/m³]
    mu_v = PropsSI('V', 'P', P_sat, 'Q', 1, fluid)
    rho_v = PropsSI('D', 'P', P_sat, 'Q', 1, fluid)

    # Reynolds numbers
    Re_l = G * D_in / mu_l
    Re_v = G * D_in / mu_v

    # Friction factors (Blasius correlation for turbulent flow)
    f_l = 0.079 / Re_l**0.25
    f_v = 0.079 / Re_v**0.25

    # Pressure gradient terms
    A = f_l * 2 * G**2 / (D_in * rho_l)
    B = f_v * 2 * G**2 / (D_in * rho_v)

    # Müller-Steinhagen-Heck pressure gradient model
    Y = A + 2 * (B - A) * x
    DP_dz = Y * (1 - x)**(1/3) + B * x**3

    # Total pressure drop over length L
    DP = DP_dz * L

    return DP


def Choi_DP(fluid, G, rho_out, rho_in, P_sat, x_o, x_i, L, D_in):
    """
    For 2 phase flow (both condensation and evaporation)
    
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
    For 2 phase flow (only condensation)
    
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

def Cheng_CO2_DP(G, D_in, L, P, h_in, mu, fluid):
    """
    Assumptions : 
    -------------
    D_in = 10 mm
    G = 496.7–1346.2 kg/m2
    heat flux (qw): 97.4 ~ 400.3 kW/m2
    Pressure (P): 7.53–23.51 MPa.
    fluid : CO2

    Froude number (Fr) : 7.58 × 10 5~1834
    Reynolds number (Re) : 6.18 × 10^4 ~ 5.35 × 10^5
    
    Reference :
    -----------
    Supercritical carbon dioxide heat transfer in horizontal tube based on the Froude number analysis (2024)
    
    Liangyuan Cheng a, Jinliang Xu a,b,*, Wenxuan Cao a, Kaiping Zhou a, Guanglin Liu a,b
    """
    
    # Determine T_plus and T_minus
    def find_pseudo_crit_T(P):
        """
        Assumptions
        -----------
        P in [7.5,14] MPa
        
        Reference
        ---------
        Investigation on the Properties of Supercritical CO2 Fluid and its Heat
        Transfer Characteristics (2012)
        
        Z. Yang & J. Yang
        """
        
        p = P*1e-6 # MPa
        T_pc = - 31.40 + 12.15*p - 0.6927*p**2 + 0.03160 * p**3 - 0.0007521 * p**4
        
        return T_pc + 273.15

    def find_h_cp_L(fluid):
        """
        Reference
        ---------
        The Latent Heat of Supercritical Fluids (2019)
        
        Daniel T. Banuti
        """
        
        P_L = PropsSI('PCRIT', fluid)*0.1
        
        if P_L < PropsSI('PTRIPLE', fluid):
            P_L = PropsSI('PTRIPLE', fluid)*1.05
        
        T_L = PropsSI('T', 'P', P_L,'Q',0,fluid)
        
        SC = 0.01
        
        cp_L = PropsSI('CPMASS', 'P', P_L,'T', T_L-SC,fluid)
        h_L = PropsSI('H', 'P', P_L,'T', T_L-SC,fluid)
        
        h_L_0 = h_L - cp_L*T_L
        
        return cp_L, h_L_0

    def T_plus_minus(P, fluid):
        """
        Calculates pseudo-boiling T_pc, T_minus, and T_plus for a supercritical fluid,
        based on Banuti (2015).

        Parameters
        ----------
        P : float
            Pressure in Pa
        fluid : str
            Fluid name (e.g., 'CO2')

        Returns
        -------
        T_pc : float
            Pseudo-boiling temperature [K]
        T_minus : float
            Lower transition limit [K]
        T_plus : float
            Upper transition limit [K]
        """
        T_pc = find_pseudo_crit_T(P)
        cp_pc = PropsSI('CPMASS', 'T', T_pc, 'P', P, fluid)
        h_0_pc = PropsSI('H', 'T', T_pc, 'P', P, fluid)

        # Find Liquid Reference properties
        cp_L, h_L_0 = find_h_cp_L(fluid)

        # Ideal gas cp approximation
        T_max = PropsSI('TMAX', fluid)
        cp_IG = PropsSI('CPMASS', 'T', T_max, 'P', P, fluid)

        # Calculate transition bounds
        T_plus = (h_0_pc - cp_pc * T_pc) / (cp_IG - cp_pc)
        T_minus = (h_L_0 - h_0_pc + cp_pc * T_pc) / (cp_pc - cp_L)
        
        return T_pc, T_minus, T_plus
    
    # Pseudo-critical temperature and definition of pseudo-boiling region
    # T_minus is liquid like temperature and T_plus is vapor like T
    T_pc, T_minus, T_plus = T_plus_minus(P, fluid)

    # Density
    rho = PropsSI('D', 'H', h_in, 'P', P, fluid)

    # Liquid-Like and Vapor-Like Enthalpies
    h_LL = PropsSI('H', 'T', T_minus, 'P', P, fluid)
    h_LV = PropsSI('H', 'T', T_plus, 'P', P, fluid)
    
    # Pseudo-vapor mass quality
    x_pb_in = (h_in - h_LL)/(h_LV - h_LL)
        
    Re = G*D_in/mu
    
    # 1) VL (Vapor Like) and LL (Liquid Like) regimes
            
    if x_pb_in < 0 or x_pb_in > 1:
        C = 1.5
        f = C*(1.82*np.log10(Re) - 1.64)**(-2)
    
    # 2) TPL (Two Phase Like) regime
    else:
        g = 9.81 # gravtity constant : m/s^2
        Fr_VL = G**2 * x_pb_in/(rho**2 *g*D_in)
    
        f = 4.379*Re**(-0.388) * Fr_VL**(-0.0167)
        
    DP = f*(L*G**2)/(2*D_in*rho)
    
    return DP
