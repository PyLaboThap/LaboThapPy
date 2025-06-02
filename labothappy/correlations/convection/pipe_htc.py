# -*- coding: utf-8 -*-
"""
@author: Basile Chaudoir
"""

from math import log10, inf
import numpy as np
import warnings
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP

#%%

def gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L):
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
    
    Heat Exchanger design based on economic optiisation
    Caputo et al.
    """
    #-------------------------------------------------------------------------
    def Stephan_Preusser(Re, Pr, Dh, L):
        Nu_1 = 0.0677*(Re*Pr*(Dh/L)**1.33)
        Nu_2 = 1 + 0.1*Pr*(Re*Dh/L)**(0.3)
        Nu = 3.657 + Nu_1/Nu_2
        return Nu
    def gnielinski_transition(Re, Pr):
        f = (1.82*np.log10(Re) - 1.64)**(-2) # (1.8*log10(Re) - 1.5)**(-2)
        Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L)**(2/3)) #*(Pr/Pr_w)**(0.11)
        return Nu
    def sieder_tate(Re, Pr):
        Nu = 0.027*Re**0.8*Pr**(1/3)*(mu/mu_w)**0.14
        return Nu
    #-------------------------------------------------------------------------
    Re_min = 0
    Re_max = 1e06
    Re = G*Dh/mu

    #-------------------------------------------------------------------------
    if Re > 10000: #fully turbulent
        Pr_min = 0.1
        Pr_max = 1000
        Nu = sieder_tate(Re, Pr)
    if Re < 2300: #fully laminar
        Pr_min = 0.6
        Pr_max = inf
        Nu = Stephan_Preusser(Re, Pr, Dh, L)
    else: #transition zone    
        Pr_min = 0.1
        Pr_max = 1000
        Nu = gnielinski_transition(Re, Pr)
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
    return hConv, Re, Pr

def dittus_boetler_heating(mu, Pr, k, G, Dh):
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
    
    h : Convection [W/(m^2 * K)]
    
    Reference
    ---------
    Incropera's Foundation of heat transfer
    """

    # Reynolds number
    Re = G*Dh/mu
    
    # Nusselt number
    Nu = 0.023*Re**0.8 * Pr**(0.4)
    
    # HTC 
    h = Nu*k/Dh
    
    return h

def dittus_boetler_cooling(mu, Pr, k, G, Dh):
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
    
    h : Convection [W/(m^2 * K)]
    
    Reference
    ---------
    Incropera's Foundation of heat transfer
    """

    # Reynolds number
    Re = G*Dh/mu
    
    # Nusselt number
    Nu = 0.023*Re**0.8 * Pr**(0.3)
    
    # HTC 
    h = Nu*k/Dh
    
    return h

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
    Pr_l = PropsSI('PRANDTL','Q',0,'P',P_sat,fluid)

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

def steiner_taborek_internal_boiling(fluid,Q_act,A,m_dot,P_sat,h_in,T_w,D_in,L):
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

def horizontal_flow_boiling(fluid, G, P_sat, x, D_in, q):
    """
    Inputs
    ------
    
    fluid : fluid name [-]
    G     : Flowrate per area unit [kg/(m**2 * s)]
    P_sat : Saturation pressure [Pa]
    x     : Vapor mass fraction [-]
    D_in  : Pipe internal diameter [m]
    q     : Heat flux [W/(m**2 * K)]

    Outputs
    -------
    
    alpha : Evaporation heat transfer coefficient [W/(m^2 * K)]
    
    Reference
    ---------
    VDI Heat Atlas 2010 
    """
    
    if fluid == 'Helium' or fluid == "He":
        q_0 = 1000
    elif fluid == 'H2' or fluid == 'N2' or fluid == 'Neon' or fluid == 'Air' or fluid == 'O2' or fluid == 'Ar':
        q_0 = 10000
    elif fluid == 'Hydrogen' or fluid == 'Nitrogen' or fluid == 'Neon' or fluid == 'Air' or fluid == 'Oxygen' or fluid == 'Argon':
        q_0 = 10000
    elif fluid == 'H2O' or fluid == 'NH3' or fluid == 'CO2' or fluid == 'SF6':
        q_0 = 150000
    elif fluid == 'Water' or fluid == 'Ammonia':
        q_0 = 150000
    else:
        q_0 = 20000
    
    a_0_dict = {
        # Substances with their corresponding formulas
        "Methane": 8060, "CH4": 8060,
        "Ethane": 5210, "C2H6": 5210,
        "Propane": 4000, "C3H8": 4000,
        "n-Butane": 3300, "C4H10": 3300,
        "n-Pentane": 3070, "C5H12": 3070,
        "Cyclopentane": 3070,
        "Isopentane": 2940,
        "n-Hexane": 2840, "C6H14": 2840,
        "n-Heptane": 2420, "C7H16": 2420,
        "Cyclohexane": 2420, "C6H12": 2420,
        "Benzene": 2730, "C6H6": 2730,
        "Toluene": 2910, "C7H8": 2910,
        "Diphenyl": 2030, "C12H10": 2030,
        "Methanol": 2770, "CH4O": 2770,
        "Ethanol": 3690, "C2H6O": 3690,
        "n-Propanol": 3170, "C3H8O": 3170,
        "Isopropanol": 2920, "C3H8O": 2920,
        "n-Butanol": 2750, "C4H10O": 2750,
        "Isobutanol": 2940, "C4H10O": 2940,
        "Acetone": 3270, "C3H6O": 3270,
        "R11": 2690, "CCl3F": 2690,
        "R12": 3290, "CCl2F2": 3290,
        "R13": 3910, "CClF3": 3910,
        "CBrF3": 3380,
        "R22": 3930, "CHClF2": 3930,
        "R23": 4870, "CHF3": 4870,
        "R113": 2180, "C2Cl3F3": 2180,
        "R114": 2460, "C2Cl2F4": 2460,
        "R115": 2890, "C2ClF5": 2890,
        "R123": 2600, "C2H2Cl2F3": 2600,
        "R134a": 3500, "C2H2F4": 3500,
        "R152a": 4000, "C2H4F2": 4000,
        "R226": 3700, "C3HClF6": 3700,
        "R227": 3800, "C3HF7": 3800,
        "RC318": 2710, "C4F8": 2710,
        "R502": 2900, "CHClF2" : 2900, "C2F5Cl": 2900,
        "Chloromethane": 4790, "CH3Cl": 4790,
        "Tetrachloromethane": 2320, "CCl4": 2320,
        "Tetrafluoromethane": 4500, "CF4": 4500,
        "Helium": 1990, "He": 1990,
        "Hydrogen": 12220, "H2": 12220,
        "Neon": 8920, "Ne": 8920,
        "Nitrogen": 4380, "N2": 4380,
        "Argon": 3870, "Ar": 3870,
        "Oxygen": 4120, "O2": 4120,
    } # W/(m**2 * K)
    
    a_0 = a_0_dict[fluid]
        
    AS = CP.AbstractState("BICUBIC&HEOS", fluid)

    P_crit = AS.p_critical()    
    P_01 = 0.1*P_crit

    # Clamp x early
    x = min(max(x, 0), 1)
    
    # Liquid state
    AS.update(CP.PQ_INPUTS, P_01, 0)
    MM_fluid = AS.molar_mass()
    rho_l_01 = AS.rhomass()
    h_l_01 = AS.hmass()
    
    # Vapor state
    AS.update(CP.PQ_INPUTS, P_01, 1)
    h_v_01 = AS.hmass()
    rho_v_01 = AS.rhomass()
    
    if fluid == 'R134a':
        sigma_01 = 0.010449499360652493
    else:
        sigma_01 = PropsSI('I', 'P', P_01, 'Q', 0.5, fluid)
    
    # n exponent
    P_star = P_sat/P_crit
    n = 0.9 - 0.36*P_star**0.13
    
    # Factor depending on the fluid molar mass
    MM_H2 = 0.00201588 # H2 Molar mass
    C_F = 0.789*(MM_fluid/MM_H2)**0.11
    
    # F_p : Contribution of pressure
    F_p = 2.692*P_star**0.43 + (1.65*P_star**6.5 / (1-P_star**4.4))
    
    # F_d : Contribution of tube diameter
    F_d = (1e-2/D_in)**0.5

    # F_W : Contribution of tube wall roughness    
    Ra = 2.3*1e-6 # between 1.6 and 6.3 *1e-6 for carbon steel pipes
    Ra_o = 1e-6
    F_W = (Ra/Ra_o)**0.133
    
    # F_G : Contribution of mass flowrate
    G_0 = 100 # kg/(m**2 * s)
    F_G = (G/G_0)**0.25
    
    # F_x : Contribution of vapor quality
    g = 9.81 # m/s**2
    Dh_evap_01 = h_v_01 - h_l_01 # Evaporation specific heat
    
    q_cr_01 = 0.13*Dh_evap_01*rho_v_01**0.5 * (sigma_01*g*(rho_l_01-rho_v_01))**0.25
    q_cr_PB = 2.79*q_cr_01*P_star**0.4*(1-P_star)
    F_x = 1 - P_star**0.1 * (q/q_cr_PB)**0.3 * x
    
    # Heat transfer coefficient
    alpha = a_0*C_F*(q/q_0)**n * F_p * F_d * F_W * F_G * F_x # W/(m*2 * K)
    
    # print(f"a_0 : {a_0}")
    # print(f"C_F : {C_F}")
    # print(f"q : {q}")
    # print(f"q_0 : {q_0}")
    # print(f"n : {n}")
    # print(f"F_p : {F_p}")
    # print(f"F_d : {F_d}")
    # print(f"F_W : {F_W}")
    # print(f"F_G : {F_G}")
    # print(f"F_x : {F_x}")
    # print("-----------------")
    
    return alpha

# def flow_boiling_gungor_winterton(fluid, G, P_sat, x, D_in, q, mu_l, Pr_l, k_l):
#     """
#     Inputs
#     ------
    
#     fluid : fluid name [-]
#     G     : Flowrate per area unit [kg/(m**2 * s)]
#     P_sat : Saturation pressure [Pa]
#     x     : Vapor mass fraction [-]
#     D_in  : Pipe internal diameter [m]
#     q     : Heat flux [W/(m**2 * K)]

#     Outputs
#     -------
    
#     h_tp : Evaporation heat transfer coeffieicnt [W/(m^2 * K)]
    
#     Reference
#     ---------
    
#     """
    
#     # Pipe roughness    
#     h_f = dittus_boetler_heating(mu_l, Pr_l, k_l, G, D_in)
#     h_lv = PropsSI('H', 'P', P_sat, 'Q', 1, fluid) - PropsSI('H', 'P', P_sat, 'Q', 0, fluid)
    
#     rho_v, mu_v = PropsSI(('D','V'), 'P', P_sat, 'Q', 1, fluid)
#     rho_l, mu_l = PropsSI(('D','V'), 'P', P_sat, 'Q', 0, fluid)

#     mu_tp = PropsSI('V', 'P', P_sat, 'Q', x, fluid)

    
#     P_crit = PropsSI('PCRIT', fluid)
#     MM = PropsSI('M', fluid)
    
#     # Froude Number
#     g = 9.81 # m/s^2
#     v_l = G/rho_l
#     Fr = (v_l/(g*D_in)**0.5)
    
#     # print(Fr)
    
#     # Reduced pressure
#     P_r = P_sat/P_crit
    
#     # Lockhart-Martinelli parameter
#     X_tt = ((1-x)/x)**0.9 * (rho_v/rho_l)**0.5 * (mu_l/mu_v)**0.1
    
#     # Boiling number
#     Bo = q/(G*h_lv)
        
#     if Fr <= 0.05:
#         E_2 = Fr**(0.1 - 2*Fr)
#         S_2 = Fr**0.5
#     else:
#         E_2 = 1
#         S_2 = 1
        
#     # Convection enhancement factor
#     E = E_2*(1 + 24000*Bo**1.16 + 1.23 * (1/X_tt)**0.86)
    
#     # Boiling suppression factor
#     A_in = (np.pi/4)*D_in**2
#     P_wet = np.pi*D_in
#     D_e = 4*A_in/P_wet
    
#     # Re_l = G*(1-x)*D_e/mu_l
#     # S = S_2*((1 + 0.00000253 * Re_l**1.17)**(-1))

#     Re_tp = G*D_e/mu_tp
#     S = S_2*((1 + 0.00000253 * Re_tp**1.17)**(-1))
    
#     # Nucleate Bpiling : Cooper Correlation
#     h_nb = 55*P_r**(0.12)  * (-np.log10(P_r))**(-0.55) * (MM*1e3)**(-0.5) * q**0.67 # - 0.2*np.log(Ra)
    
#     h_tp = E*h_f + S*h_nb
    
#     # print(f"E : {E}")
#     # print(f"h_f : {h_f}")
#     # print(f"S : {S}")
#     # print(f"h_nb : {h_nb}")
    
#     return h_tp

def flow_boiling_gungor_winterton(fluid, G, P_sat, x, D_in, q, mu_l, Pr_l, k_l):
    """
    Inputs
    ------
    
    fluid : fluid name [-]
    G     : Flowrate per area unit [kg/(m**2 * s)]
    P_sat : Saturation pressure [Pa]
    x     : Vapor mass fraction [-]
    D_in  : Pipe internal diameter [m]
    q     : Heat flux [W/(m**2 * K)]

    Outputs
    -------
    
    h_tp : Evaporation heat transfer coeffieicnt [W/(m^2 * K)]
    
    Reference
    ---------
    
    """
    
    AS = CP.AbstractState("BICUBIC&HEOS", fluid)
    
    AS.update(CP.PQ_INPUTS, P_sat, 1)
    
    h_v = AS.hmass()
    rho_v = AS.rhomass()
    mu_v = AS.viscosity()

    AS.update(CP.PQ_INPUTS, P_sat, 0)

    h_l = AS.hmass()
    rho_l = AS.rhomass()
    mu_l = AS.viscosity()
    
    AS.update(CP.PQ_INPUTS, P_sat, x)    
    
    mu_tp = AS.viscosity()
    P_crit = AS.p_critical()
    MM = AS.molar_mass()

    # Pipe roughness    
    h_f = dittus_boetler_heating(mu_l, Pr_l, k_l, G, D_in)
    h_lv = h_v - h_l
    
    # Froude Number
    g = 9.81 # m/s^2
    v_l = G/rho_l
    Fr = (v_l/(g*D_in)**0.5)
    
    # print(Fr)
    
    # Reduced pressure
    P_r = P_sat/P_crit
    
    # Lockhart-Martinelli parameter
    X_tt = ((1-x)/x)**0.9 * (rho_v/rho_l)**0.5 * (mu_l/mu_v)**0.1
    
    # Boiling number
    Bo = q/(G*h_lv)
        
    if Fr <= 0.05:
        E_2 = Fr**(0.1 - 2*Fr)
        S_2 = Fr**0.5
    else:
        E_2 = 1
        S_2 = 1
        
    # Convection enhancement factor
    E = E_2*(1 + 24000*Bo**1.16 + 1.23 * (1/X_tt)**0.86)
    
    # Boiling suppression factor
    A_in = (np.pi/4)*D_in**2
    P_wet = np.pi*D_in
    D_e = 4*A_in/P_wet
    
    # Re_l = G*(1-x)*D_e/mu_l
    # S = S_2*((1 + 0.00000253 * Re_l**1.17)**(-1))

    Re_tp = G*D_e/mu_tp
    S = S_2*((1 + 0.00000115*E**2 * Re_tp**1.17)**(-1))
    
    # Nucleate Bpiling : Cooper Correlation
    h_nb = 55*P_r**(0.12)  * (-np.log10(P_r))**(-0.55) * (MM*1e3)**(-0.5) * q**0.67 # - 0.2*np.log(Ra)
    
    h_tp = E*h_f + S*h_nb
    
    # print(f"E : {E}")
    # print(f"h_f : {h_f}")
    # print(f"S : {S}")
    # print(f"h_nb : {h_nb}")
    
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


def Cheng_sCO2(G, q, T_w, P, h_in, h_out, mu, k, D_in, fluid):
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
        cv_pc = PropsSI('CVMASS', 'T', T_pc, 'P', P, fluid)
        h_0_pc = PropsSI('H', 'T', T_pc, 'P', P, fluid)
        MM = PropsSI('M', fluid)

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

    # Wall Enthalpy and density
    h_w = PropsSI('H', 'T', T_w, 'P', P, fluid)
    h_pc = PropsSI('H', 'T', T_pc, 'P', P, fluid)

    rho_w = PropsSI('D', 'T', T_w, 'P', P, fluid)

    # Supercritical Boiling Number (wall conditions)
    SBO = q/(G*h_pc)
    SBO_w = q/(G*h_w)

    # Liquid-Like and Vapor-Like Enthalpies
    h_LL = PropsSI('H', 'T', T_minus, 'P', P, fluid)
    h_LV = PropsSI('H', 'T', T_plus, 'P', P, fluid)
    
    # Pseudo-vapor mass quality
    x_pb_in = (h_in - h_LL)/(h_LV - h_LL)
    x_pb_out = (h_out - h_LL)/(h_LV - h_LL)
    
    x_ave = 0.5*(x_pb_in + x_pb_out)
    
    Re = G*D_in/mu
    
    # 1) VL (Vapor Like) and LL (Liquid Like) regimes
    
    T_in = PropsSI('T', 'H', h_in, 'P', P, fluid)
    cp_ave = (h_w - h_in)/(T_w - T_in)
    
    Pr_ave = mu*cp_ave/k
    
    if x_ave < 0 or x_ave > 1:
        Nu = 0.019*(Re**0.79)*(Pr_ave**0.51)
        h_conv = Nu*k/D_in
        
        return h_conv
    
    # 2) TPL (Two Phase Like) regime
    else:
        g = 9.81 # gravtity constant : m/s^2
        Fr_VL_w = G**2 * x_pb_in/(rho_w**2 *g*D_in)
        
        Nu_top = 0.127*(Re**1.065) * (Pr_ave**1.053) * (SBO_w**0.66) * (Fr_VL_w**(-0.012))
        Nu_bot = 1.512*(Re**0.878) * (Pr_ave**0.685) * (SBO_w**0.69) * (Fr_VL_w**(-0.0213))
        
        h_top = Nu_top*k/D_in
        h_bot = Nu_bot*k/D_in
                
        # print(f"Re : {Re}")
        # print(f"Pr_ave : {Pr_ave}")
        
        # print(f"SBO_w : {SBO_w}")
        # print(f"SBO : {SBO}")
        
        # print(f"Fr_VL_w : {Fr_VL_w}")
        
        # print(f"Nu_top : {Nu_top}")
        # print(f"Nu_bot : {Nu_bot}")
        
        # print(f"h_top : {h_top}")
        # print(f"h_bot : {h_bot}")     
        
        return (h_top + h_bot)/2
        
def Liu_sCO2(G, P, T_w, k, rho, mu, cp, D_in, fluid):
    """
    Assumptions 
    -----------
    D_in = 4-10.7 [mm]
    T_in = 25-67 [°C]
    P_in = 7.5-8.5 [MPa]
    G = 74-200 [kg/(m^2 * s)]
    fluid : CO2
    
    Reference
    ---------
    Experimental study on heat transfer and pressure drop of supercritical CO2 cooled in a large tube (2014)
    Zhan-Bin Liu, Ya-Ling He*, Yi-Fan Yang, Ji-You Fei
    """
    AS = CP.AbstractState("BICUBIC&HEOS", fluid)  
    AS.update(CP.PT_INPUTS, P, T_w)
    
    rho_w = AS.rhomass()
    cp_w = AS.cpmass()
    mu_w = AS.viscosity()
    Pr_w = AS.Prandtl()
        
    Re_w = G*D_in/mu_w
    
    Nu = 0.01*Re_w**0.9 * Pr_w**0.5 * (rho_w/rho)**0.906 * (cp_w/cp)** (-0.585)
    
    h_conv = Nu*k/D_in
    
    return h_conv