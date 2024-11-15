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
def gnielinski_pipe_htc(mu, Pr, Pr_w, k, G, Dh, L):
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
    
    hConv : Convection heat transfer coefficient [W/(m^2 * K)]
    
    Reference
    ---------
    Validation of the Gnielinski correlation for evaluation of heat transfer coefficient of enhanced tubes by non-linear 
    regression model: An experimental study of absorption refrigeration system
    
    Syed Muhammad Ammar, Chan Woo Park
    """
    #-------------------------------------------------------------------------
    def gnielinski_laminar(Re, Pr, Dh, L):
        Nu_1 = 4.364
        Nu_2 = 1.953*(Re*Pr*Dh/L)**(1/3)
        Nu = (Nu_1**3 + 0.6**3 + (Nu_2 - 0.6)**3)**(1/3)
        return Nu
    def gnielinski_turbulent(Re, Pr):
        f = (1.8*log10(Re) - 1.5)**(-2)
        Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L)**(2/3))*(Pr/Pr_w)**(0.11)
        return Nu
    #-------------------------------------------------------------------------
    Re_min = 0
    Re_max = 1e06
    Re = G*Dh/mu
    #-------------------------------------------------------------------------
    if Re > 1e4: #fully turbulent
        Pr_min = 0.1
        Pr_max = 1000
        Nu = gnielinski_turbulent(Re, Pr)
    elif Re < 2300: #fully laminar
        Pr_min = 0.6
        Pr_max = inf
        Nu = gnielinski_laminar(Re, Pr, Dh, L)
    else: #transition zone
        Pr_min = 0.1
        Pr_max = 1000
        gamma = (Re - 2300)/(1e4 - 2300)
        Nu_lam2300 = gnielinski_laminar(2300, Pr, Dh, L)
        Nu_turb10000 = gnielinski_turbulent(1e4, Pr)
        Nu = (1-gamma)*Nu_lam2300 + gamma*Nu_turb10000
    #-------------------------------------------------------------------------
    hConv = Nu*k/Dh
    
    #-------------------------------------------------------------------------
    if Re >= Re_max or Re <=Re_min:
        # warnings.warn('Gnielinski singe-phase: Out of validity range --> Re = ', Re, ' is out of [', Re_min, ' - ', Re_max, '] !!!')
        warnings.warn('Gnielinski singe-phase: Reynolds Out of validity range !!!')
    if Pr >= Pr_max or Pr <= Pr_min:
        # warnings.warn('Gnielinski singe-phase: Out of validity range --> Re = ', Pr, ' is out of [', Pr_min, ' - ', Pr_max, '] !!!')
        warnings.warn('Gnielinski singe-phase: Prandtl Out of validity range  !!!')
    #-------------------------------------------------------------------------
    return hConv

def pipe_internal_DP(mu, Pr, Np, rho, G, Dh, L):
    """
    Inputs
    ------
    
    mu   : Dynamic Viscosity [Pa*s]
    Pr   : Prantl number [-]
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
    
    # Reynolds number
    Re = G*Dh/mu
    
    # Flow speed
    u = G/rho
    
    # Friction coefficient
    f = (1.8*log10(Re) - 1.5)**(-2)
    
    # Pressure drop (Np : number of passes)
    DP = (4*f*L*Np/Dh + 4*Np)*rho*(u**2)/2
    
    return DP

def horizontal_tube_internal_condensation(fluid,m_dot,P_sat,h_in,T_w,D_in):
    """
    Inputs
    ------
    
    fluid : fluid name [-]
    m_dot : Flowrate [kg/s]
    P_sat : Saturation pressure [Pa]
    h_in  : Inlet enthalpy [J/kg]
    T_w   : Wall temperature [K]
    D_in  : Pipe internal diameter [m]
    
    Outputs
    -------
    
    h : Condensation heat transfer coeffieicnt [W/(m^2 * K)]
    
    Reference
    ---------
    Incropera's principle of heat transfer
    """

    # Geometrical parameters
    A_in = np.pi*(D_in/2)**2
    g = 9.81 # gravity acceleration constant m/s^2

    # 2 phase properties
    x = PropsSI('Q','H',h_in,'P',P_sat,fluid)
    h_fg = PropsSI('H','Q',1,'P',P_sat,fluid) - PropsSI('H','Q',0,'P',P_sat,fluid)
    T_sat = PropsSI('T','Q',0.5,'P',P_sat,fluid)
    
    # Vapor properties
    rho_v = PropsSI('D','Q',1,'P',P_sat,fluid)
    mu_v = PropsSI('V','Q',1,'P',P_sat,fluid)
    
    u_v = m_dot/(rho_v*A_in)

    Re_v = (rho_v*u_v*D_in)/mu_v

    # Liquid properties
    mu_l = PropsSI('V','Q',0,'P',P_sat,fluid)
    k_l = PropsSI('L','Q',0,'P',P_sat,fluid)
    rho_l = PropsSI('D','Q',0,'P',P_sat,fluid)
    cp_l = PropsSI('C','Q',0,'P',P_sat,fluid)
    Pr_l = PropsSI('PRANDTL', 'Q',0,'P',P_sat,fluid)

    if Re_v <= 35000: # Low speed vapor flow
        # Dobson and Chato
        C = 0.555    
        h_2_fg = h_fg + 0.375*cp_l*(T_sat - T_w)
        Nu = C*((rho_l*g*(rho_l - rho_v)*h_2_fg*D_in**3)/(mu_l*k_l*(T_sat - T_w)))**(0.25)
    else:
        # Liquid and Reynolds and Prandtl numbers
        Re_Dl = 4*m_dot*(1-x)/(np.pi*D_in*mu_l)
        
        # Martinelli parameter
        x_tt = ((1-x)/x)**(0.9)*(rho_v/rho_l)**(0.5)*(mu_l/mu_v)**(0.1)
        
        # Nusselt
        Nu = 0.023*Re_Dl**0.8 * Pr_l **0.4 * (1 + (2.22/x_tt**(0.89)))
    
    h = Nu*k_l/D_in
        
    return h

def horizontal_tube_internal_boiling(fluid,Q_act,A,m_dot,P_sat,h_in,T_w,D_in,L):
    """
    Inputs
    ------
    
    fluid : fluid name [-]
    Q_act : Actual heat transfer [W]
    A     : Heat exchange area [m^2]
    m_dot : Flowrate [kg/s]
    P_sat : Saturation pressure [Pa]
    h_in  : Inlet enthalpy [J/kg]
    T_w   : Wall temperature [K]
    D_in  : Pipe internal diameter [m]
    L     : Flow length [m]
    
    Outputs
    -------
    
    h_tp : Boiling heat transfer coeffieicnt [W/(m^2 * K)]

    Reference
    ---------
    Boiling Heat Transfer Inside Plain Tubes - Wolverine Tube, INC
    
    Steiner-Taborek Asymptotic Model
    """
        
    # Geometrical parameters
    A_in = np.pi*(D_in/2)**2
    
    # 2 phase properties
    x = PropsSI('Q','H',h_in,'P',P_sat,fluid)
    G = m_dot/A_in
    sigma = PropsSI('I','P',P_sat,'H',h_in,fluid) 
    T_sat = PropsSI('T','P',P_sat,'Q',0,fluid) 
    h_lv = PropsSI('H','Q',1,'P',P_sat,fluid) - PropsSI('H','Q',0,'P',P_sat,fluid)
    P_r = P_sat/PropsSI('PCRIT','Cyclopentane')

    # Vapor properties
    rho_v = PropsSI('D','Q',1,'P',P_sat,fluid)
    mu_v = PropsSI('V','Q',1,'P',P_sat,fluid)
    k_v = PropsSI('L','Q',1,'P',P_sat,fluid)
    rho_v = PropsSI('D','Q',1,'P',P_sat,fluid)
    Pr_v = PropsSI('PRANDTL', 'Q',1,'P',P_sat,fluid)
    
    # Liquid properties
    mu_l = PropsSI('V','Q',0,'P',P_sat,fluid)
    k_l = PropsSI('L','Q',0,'P',P_sat,fluid)
    rho_l = PropsSI('D','Q',0,'P',P_sat,fluid)
    Pr_l = PropsSI('PRANDTL', 'Q',0,'P',P_sat,fluid)
    
    # Wall property
    Pr_w = PropsSI('PRANDTL','T',T_w,'P',P_sat,fluid)
    
    # Liquid convective HTC
    h_L = gnielinski_pipe_htc(mu_l, Pr_l, Pr_w, k_l, G, D_in, L)
    
    # Heat flux
    q_act = Q_act/A
    
    # Onset of nucleate boiling
    r_o = 0.3*1e-6
    q_ONB = (2*sigma*T_sat*h_L)/(r_o*rho_v*h_lv)
    
    if q_act < q_ONB:
        # Vapor convective HTC
        h_v = gnielinski_pipe_htc(mu_v, Pr_v, Pr_w, k_v, G, D_in, L)
        
        a = (1-x)**1.5 + 1.9*x**0.6 * (1-x)**0.1 * (rho_l/rho_v)**0.35
        b = (h_v/h_L)*x**0.01*(1+8*(1-x)**0.7)*(rho_l/rho_v)**0.67
        
        if a == 0:
            F_tp = 1
        elif b == 0:
            F_tp = 1
        else:
            F_tp = (a**(-2.2) + b**(-2))**(-0.5)
        
        h_tp = F_tp*h_v
        
    else:   
        # Parameters
        M = 72.15
        q_o = 20000
        h_nb_o = 3070*(2420/2840) # Assume h_nb,o,Cyclopentane = h_nb,o,Cyclo-hexane * (h_nb,o,n-pentane/h_nb,o,n-hexane)
        D_in_o = 0.01 # m
        
        # Two-phase multiplier F_tp
        F_tp = ((1-x)**1.5 + 1.9*x**0.6 * (rho_l/rho_v)**0.35)**1.1
        
        # Nucleate boiling pressure correction factor
        F_pf = 2.816*P_r**0.45 + (3.4 + (1.7/(1 - P_r**7)))*P_r**3.7
        
        # Nucleate boiling exponent
        nf = 0.8 - 0.1*np.e**(1.75*P_r)
        
        # Molecular function
        F_M = 0.377 + 0.199*np.log(M) + 0.000028427*M**2
        
        # R_p ratio
        R_p = 1e-6 # m : Assumption
        R_p_o = 1e-6 # m
        
        # Nucleate boilin multiplier F_nb
        F_nb = F_pf*(q_act/q_o)**nf * (R_p/R_p_o)**(0.133) * (D_in/D_in_o)**(-0.4) * F_M
                
        # Two-phase htc
        h_tp = ((h_nb_o*F_nb)**3 + (h_L*F_tp)**3)**(1/3)

    return h_tp

def pool_boiling(fluid, T_sat, T_tube):
    """
    ---- Inputs : -------- 
    
    fluid  : fluid name [-]
    T_sat  : External fluid saturation temperature [K]
    T_tube : Tube temperature [K]
    
    ---- Outputs : --------
    
    h_pool : Heat transfer coefficient for pool boiling [W/(m^2*K)]
    
    ---- Reference(s) : --------
    
    Van Long Le - Heat Pipe Code
    
    """
    
    g = 9.81 # m/s^2 : gravity acceleration constant
    
    # Fluid Data # !!! Replace this with the oil data !!!
    rho_l = PropsSI('D','T',T_sat,'Q',0,fluid) # density at saturated liquid in kg/m^3
    rho_v = PropsSI('D','T',T_sat,'Q',1,fluid) # density at saturated vapor in kg/m^3
    mu_l = PropsSI('V','T',T_sat,'Q',0,fluid) # viscosity at saturated liquid in Pa*s
    cp_l = PropsSI('C','T',T_sat,'Q',0,fluid) # specific heat at saturated liquid
    sigma = PropsSI('I','T',T_sat,'Q',0.5,fluid) # surface tension of liquid vapor equilibrium
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific heat
    
    # Heat Atlas related parameters
    C_nb = 0.013 # Pentane on polished nickel
    n = 1.7 # Prandtl exponent for other fluids than water
    DT_e = T_tube - T_sat # Temperature difference between fluid and tube
    Pr_l = PropsSI('Prandtl','T',T_sat,'Q',0,fluid)
    
    coef_1 = mu_l*Dh_evap
    coef_2 = g*(rho_l-rho_v)/sigma
    coef_3 = (cp_l*DT_e)/(C_nb*Dh_evap*Pr_l**n)
    q_pool = coef_1*coef_2**(1/2)*coef_3**3
    
    return q_pool 

def film_boiling(D_out,fluid,T_sat,T_tube,P_sat):
    """
    ---- Inputs : -------- 
    
    D_out  : Outer Tube diameter [m]
    fluid  : fluid name [-]
    T_sat  : External fluid saturation temperature [K]
    T_tube : Tube wall temperature [K]
    P_sat  : Saturation Pressure [Pa]
    
    ---- Outputs : --------
    
    h_film : Heat transfer coefficient for film boiling [W/(m^2*K)]
    
    ---- Reference(s) : --------
    
    Foundations of heat transfer : Frank P. Incropera
    
    """
    g = 9.81 # gravity constant # [m/s^2]
    e = 0.9 # steel emissivity
    sigma = 5.67*1e-8 # Stefan-Boltzmann constant [W/(m^2*K^4)]
    
    T_f = (T_sat + T_tube)/2 # film mean temperature [K]
    
    (k_v,rho_v,nu_v,cp_v) = PropsSI(('L','D','V','C'), 'T',T_f,'P',P_sat,fluid) # Film thermal conductivity, density, kinematic viscosity, enthalpy and capacity
    h_v = PropsSI('H', 'T',T_sat,'Q',1,fluid) # Saturated vapor enthalpy
    (rho_l,h_l) = PropsSI(('D','H'), 'T',T_sat,'Q',0,fluid) # Saturated liquid density and enthalpy
    C = 0.62 # correlation constant for horizontal cylinders
        
    Dh_evap = h_v - h_l
    Dh_evap_corrected = Dh_evap + 0.8*cp_v*(T_tube - T_sat) # Accounts for sensible heat to maintain temperature of the film above saturation
    
    coef = (g*(rho_l - rho_v)*Dh_evap_corrected*D_out**3)/(nu_v*k_v*(T_tube - T_sat))
    
    h_conv = (D_out/k_v)*C*coef**(1/4) # convection coefficient
    h_rad = e*sigma*(T_tube**4 - T_sat**4)/(T_tube - T_sat)
    
    # Define the equation as a function
    def equation_to_solve(x, y, z):
        return x**(4/3) - (y**(4/3) + z*x**(1/3))
    
    # Define a function with parameters y and z
    def equation_to_solve_with_parameters(x):
        return equation_to_solve(x, h_conv, h_rad)
    
    # Use fsolve to find the numerical solution
    (h_film) = fsolve(equation_to_solve_with_parameters,x0 = [100])
    
    return h_film


def boiling_curve(D_out, fluid, T_sat, P_sat):
    """
    ---- Inputs : -------- 
    
    D_out : Outer Tube diameter [m]
    fluid (char) : fluid name (Cyclopentane)
    T_sat : External fluid saturation temperature [K]
    P_sat : External fluid saturation Pressure [Pa]
    
    ---- Outputs : --------
    
    h_final : Heat transfer coefficient vector for 0.01-1000 [K] of surface temperature difference [W/(m^2)]
                
    Surface temperature difference = T_wall - T_sat_fluid
    
    ---- Reference(s) : --------
    
    Foundations of heat transfer : Frank P. Incropera
    
    """
    
    # Vector creation
    DT = np.linspace(1,1000,1000) # Temperature differences
    q_pool = np.zeros(len(DT)) # pool boiling values vector
    q_film = np.zeros(len(DT)) # film boiling values vector
    q_final = np.zeros(len(DT)) # final boiling curve values vector

    g = 9.81 # gravity accelerattion constant
    
    # q_crit related parameters 
    C_crit = 0.15
    
    # q_min related parameters
    C_mhf = 0.09
    
    rho_l = PropsSI('D','T',T_sat,'Q',0,fluid) # density at saturated liquid in kg/m^3
    rho_v = PropsSI('D','T',T_sat,'Q',1,fluid) # density at saturated vapor in kg/m^3
    sigma = PropsSI('I','T',T_sat,'Q',0,fluid) # surface tension of liquid vapor equilibrium
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific heat
    
    # q value ending nucleate boiling phase
    q_crit = C_crit*Dh_evap*rho_v*((sigma*g*(rho_l-rho_v))/rho_v**2)**(1/4)
    
    # q value starting film boiling phase
    q_min = C_mhf*rho_v*Dh_evap*((g*sigma*(rho_l-rho_v))/(rho_l+rho_v)**2)**(1/4)
    
    for i in range(len(DT)):
        T_surf = DT[i] + T_sat
        q_pool[i] = pool_boiling(fluid, T_sat, T_surf)
        h_film = film_boiling(D_out, fluid, T_sat, T_surf, P_sat)
        q_film[i] = h_film*(DT[i])
    
    # Critical point : end of nucleate boiling
    DT_crit = np.argmin(np.abs(q_pool - q_crit))
    
    # Leidenfrost point : start of film boiling
    DT_min = np.argmin(np.abs(q_film - q_min))
    
    # Transition boiling phase slope (considered as linear)
    Delta = (q_min-q_crit)/(DT_min - DT_crit)

    for i in range(len(DT)):
        if DT[i] <= DT_crit:
            q_final[i] = q_pool[i]
        elif DT[i] > DT_min:
            q_final[i] = q_film[i]
        else:
            q_final[i] = q_crit + Delta*(DT[i] - DT_crit)
    
    h_final = q_final/DT
    
    DT = np.concatenate((np.array([0]),DT),axis = 0)
    h_final = np.concatenate((np.array([0]),h_final),axis = 0)
    
    return h_final, DT
