# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 13:34:29 2023

@author: Basile
"""

import numpy as np
import CoolProp
import math
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
from scipy.interpolate import interp2d

"""
Functions in this library are used for preliminary computations related to the internal behaviours of the heat pipe. 
    - figures_of_merit : Figures of merit of the fluid inside the heat pipe used in thermal_res_esdu
    - thermal_res_esdu : Internal thermal resistances of the heat pipe
    - Delta_P_v : Pressure loss inside the heat pipe
"""

def figures_of_merit(fluid, T_sat):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    F_r : Filling ratio of the thermosyphon [/]
    T_sat : Saturation temperature of the fluid [K]
    
    ---- Outputs : --------
    
    phi_cond : Figure of merit for condensation
    phi_boil : Figure of merit for boiling

    ---- Reference(s) : --------
    
    Investigation of low Global Warming Potential working fluids for a closed two-phase thermosyphon
    Robert W. MacGregor, Peter A. Kew, David A. Reay
    
    """
    
    if T_sat > 273.15 + 373.6:
        T_sat = 273.15 + 373.6
    
    P_atm = 101325 # [Pa] : atmospheric pressure
    
    "1) Get fluid data"
    
    # Liquid properties
    (P_sat, k_l, rho_l, mu_l, cp_l) = PropsSI(('P','L','D','V','C'), 'T', T_sat, 'Q', 0, fluid) # Sat Pressure, thermal conductivity, density, viscosity, spec. heat capacity
    # Gaseous properties
    rho_v = PropsSI('D', 'T', T_sat, 'Q', 1, fluid) # density
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg
    
    "2) Compute figures of merit"
    
    # Condensation
    phi_cond = ((Dh_evap*k_l**3*rho_l**2)/mu_l)**(1/4)

    # Boiling
    phi_boil = (0.32*(rho_l**0.65)*(k_l**0.3)*(cp_l**0.7))/(rho_v**(0.25)*Dh_evap**(0.4)*mu_l**(0.1))*(P_sat/P_atm)**(0.23)
    
    return (phi_cond, phi_boil)

def thermal_res_esdu(fluid, F_r, D_i, L_c, L_e, Q_dot_radial, T_sat):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    F_r : Filling ratio of the thermosyphon [/]
    D_i : Themosyphon internal diameter [m]
    L_c : Length of condensation zone [m]
    L_e : Length of evaporation zone [m]
    Q_dot_radial : Radial heat flux [W]
    T_sat : Saturation temperature of the fluid [K]
    
    ---- Outputs : --------
    
    R_cond : Internal Thermal Resistance for condenser [K/W]
    R_evap : Internal Thermal Resistance for evaporator [K/W] (composed of R_pool : pool boiling resistance 
                                                               and R_film : Liquid film resistance) 
    Re_f : Reynolds number of the liquid flow in the thermosyphon
    
    ---- Reference(s) : --------
    
    Investigation of low Global Warming Potential working fluids for a closed two-phase thermosyphon
    Robert W. MacGregor, Peter A. Kew, David A. Reay
    
    """
  
    g = 9.81 # [m/s^2] : gravity acceleration constant   
    
    "1) Fluid properties"

    mu_l = PropsSI('V', 'T', T_sat, 'Q', 1, fluid) # liquid viscosity [Pa*s]
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg

    # Reynolds number
    Re_f = 4*Q_dot_radial/(Dh_evap*mu_l*np.pi*D_i)
    
    "2) Condenser resistance"
    
    # Figures of merit
    (phi_cond, phi_boil) = figures_of_merit(fluid, T_sat)
    
    # Intermediate resistance
    R_cond_EDSU = 0.235*Q_dot_radial**(1/3)/(D_i**(4/3)*g**(1/3)*L_c*phi_cond**(4/3))
    
    if Re_f > 1300:
        R_cond = R_cond_EDSU*191*Re_f**(-0.733)
    elif Re_f < 50: 
        R_cond = R_cond_EDSU		
    else:
        R_cond = R_cond_EDSU

    "3) Evaporator resistance"
    
    # Film resistance
    R_film = R_cond*(L_c/L_e)
    
    # Pool resistance
    R_pool = 1/(phi_boil*g**(0.2)*Q_dot_radial**(0.4)*(np.pi*D_i*L_e)**(0.6))

    # Evap resistance
    if R_film > R_pool:
        R_evap = R_pool
    else:
        R_evap = R_pool*F_r + R_film*(1-F_r)
    
    return (R_cond, R_evap, R_pool, R_film, Re_f)

def Delta_P_v(fluid, D_v, L_eff, Q_dot_radial, T_sat):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    D_v : Vapor diameter (assumed to be the thermosyphon internal diameter) [m]
    L_eff : Effective length of the thermosyphon [m]
    Q_dot_radial : Radial heat flux [W]
    T_sat : Saturation temperature of the fluid [K]
    
    ---- Outputs : --------
    
    DP_v : Pressure drop of fluid inside the thermosyphon [Pa]
    
    ---- Reference(s) : --------
    
    /
    
    """
    
    if T_sat > 273.15 + 370:
        T_sat = 273.15 + 370
    
    "1) Get fluid data"
    
    (P_sat, rho_v, mu_v) = PropsSI(('P','D','V'),'T',T_sat,'Q',0,fluid) # Gas : (Pressure : Pa, density : kg/m^3, viscosity : Pa*s) 
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg

    A_v = (np.pi*D_v**2)/4 # vapour cross section area
    # u_v = Q_dot_radial/(A_v*Dh_evap*rho_v)
    
    "2) Compute Pressure Difference"
    
    # Reynolds number
    Re_a = (Q_dot_radial*D_v)/(mu_v*A_v*Dh_evap)
    
    if Re_a < 2300:
        DeltaP_v = (32*mu_v*Q_dot_radial*L_eff)/(rho_v*A_v*Dh_evap*D_v**2)
    else:
        DeltaP_v = (0.3164/Re_a**(1/4))*(Q_dot_radial**2 * L_eff)/(2*rho_v*A_v**2*Dh_evap**2*D_v)
        
    return DeltaP_v

