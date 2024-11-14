# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 13:54:18 2023

@author: Basile
"""

import __init__

import numpy as np
import CoolProp
import math
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize
from scipy.interpolate import interp2d, interp1d
from library.component.sizing.heat_exchanger.heat_pipe_HTX.modules.HP_internal import Delta_P_v, thermal_res_esdu
from library.component.sizing.heat_exchanger.heat_pipe_HTX.modules.HP_h_coeffs import radiative_coeff, external_flow_inline_bank, pool_boiling, external_flow_staggered_bank, ext_conv_boil #, h_cond_Th66, external_flow_finnedtubebank

"""
Functions in this library are used for computations related to a single heat pipe tube.  
    - operating_limits : Defines limits on power based on several limiting behaviors susceptible to happen in a heat pipe
    - equ_syst : Contains a system solved iteratively in order to compute temperatures inside the heat pipe
    - thermosyphon_model : Allows to compute final tube outputs by solving equ_syst's system
"""

def operating_limits(fluid, F_r, D_i, L_e, beta, geo, T_sat):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    F_r : Filling ratio of the thermosyphon [/]
    D_i : Themosyphon internal diameter [m]
    L_e : Length of evaporation zone [m]
    beta : Thermosyphon inclination angle [rad]
    geo (char) : geometry of the thermosyphon (circular or annular)
    T_sat : Saturation temperature of the fluid [K]
    
    ---- Outputs : --------
    
    Q_dot_ent (Entrainment limit) [W] : Even when there is sufficient liquid present in the thermosyphon to prevent dry-out occuring, 
                                        the overall rate of heat transfer is subject to another limit; this occurs when the rate of 
                                        entrainment of liquid by the vapour prevents the downward flow of liquid.
                                        
    Q_dot_boil (Boiling limit) [W] : Boiling limit occurs when a stable film of vapour is formed between the liquid and the heated 
                                     wall of the evaporator.
        
    Q_dot_son (Sonic limit) [W] : At low operating pressures, the vapour velocity may be appreciate compared with the sonic velocity 
                                  in the vapour.
        
    Q_dot_dry (Dryout limit) [W] : As applied to thermosyphons, the term "dryout" implies that the volume of the liquid fill is not 
                                   sufficient to cover all of the pipe above the pool with a film of liquid. Thus with a vertical pipe,
                                   most of the falling film of liquid would have evaporated before reaching the pool, leaving dry 
                                   patches, with a few rivulets of liquid returning to the pool; with an inclined pipe, dry patches 
                                   would appear at the top of the evaporator. the available evidence suggest that dryout in a vertical 
                                   thermosyphon is avoided if the volume of liquid fill meets the contitions called for in 
                                   ESDU 81038 Section 2.3 .
    
    Q_dot_vap (Vapor pressure limit) [W] : When operating a thermosyphon at a pressure substantially below atmospheric, the pressure 
                                           drop of the vapor may be significant compared to the pressure in the evaporator. Vapor 
                                           pressures are low but necessarily exceed zero at temperature close to the bottom of the 
                                           operational range of a heat pipe. The mimimum vapour pressure, which occurs at the closed 
                                           end of the condenser, can be very small. The pressure drop in the vapour duct, Deltap_v, is 
                                           then constrained by this effectively zero pressure and by the low vapour pressure existing 
                                           at the closed end of the evaporator. Because Deltap_v increases with the overall rate of heat 
                                           transfer, the constraint on Deltap_v requires Q_dot (if not otherwise limited) to be 
                                           approximately limited to a value, called the vapour limit in this text.
    
    ---- Reference(s) : --------
    
    ESDU. Heat Pipes - Performance of Two phase Closed Thermosyphons. London, U.K.: Engineering Sciences Data Unit; 1981."
    Enhancement of heat transport in thermosyphon airpreheater at high temperature with binary working fluid : A case study of TEG–water (F_1 with respect to Bo)
    
    """
    
    g = 9.81 # m/s^2 : gravity acceleration constant
    
    "1) Get fluid data"
    
    (P_sat, rho_v, mu_v, sigma) = PropsSI(('P','D','V','I'),'T',T_sat,'Q',0,fluid) # Gas : (Pressure : Pa, density : kg/m^3, viscosity : Pa*s, surface tension : N/m) 
    (rho_l, mu_l) = PropsSI(('D','V'),'T',T_sat,'Q',1,fluid) # Liquid : (Density : kg/m^3, viscosity : Pa*s)
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg
        
    "2) Parameter F_1"
    
    # Bond number (or Eötvös number, dimensionless number measuring the importance of gravitational forces compared to surface tension forces for the movement of liquid front)
    if geo == 'circular':
        Bo = D_i*(g*abs(rho_l-rho_v)/sigma)**(1/2)
    elif geo == 'annular':
        Bo = (D_i/2)*(g*abs(rho_l-rho_v)/sigma)**(1/2)
    else:
        print("Error. Geometry shall be defined to either 'annular' or 'circular'. \n")
        return
    
    # F_1 parameter
    
    if Bo > 11:
        F_1 = 8.2
    else:
        F_1 = -0.0331*Bo**2 + 0.8161*Bo + 3.2134

    "3) Parameter F_2"
    
    # Dimensionless pressure parameter
    
    K_p = P_sat/(g*abs(rho_l-rho_v)*sigma)**(1/2)
    
    # F_2 parameter
    
    if K_p <= 4e4:
        F_2 = K_p**(-0.17)
    else:
        F_2 = 0.165    
        
    "4) Parameter F_3 : function of inclination of the pipe (beta in rad)"
    
    if beta == np.pi/2:
        F_3 = 1
    else:
        if Bo>=4:
            F_3 = 0.306075777 + 3.88126365*beta - 8.71941343*beta**2 + 13.4517469*beta**3 - 11.9947367*beta**4 + 5.32966132*beta**5 - 0.93010743*beta**6
        elif Bo >=2:
            F_3 = 0.222344657 + 2.71186727*beta - 6.78785522*beta**2 + 13.1807234*beta**3 - 13.6843299*beta**4 + 6.73969856*beta**5 - 1.26256636*beta**6
        else:
            F_3 = 0.213964282 + 0.922627057*beta + 0.285274233*beta**2 - 0.721665592*beta**3 + 0.235369565*beta**4
    
    "5) Compute limits"
    
    # Entrainment limit (F_1*F_2*F_3 group can be called Kutateladze number)
    Q_dot_ent = F_1*F_2*F_3*Dh_evap*(rho_v**(1/2))*(sigma*g*abs(rho_l-rho_v))**(1/4)
    
    # Boiling limit
    Q_dot_boil = 0.12*rho_v*Dh_evap*((sigma*abs(rho_l-rho_v)*g)/rho_v**2)**(1/4)
        
    # Sonic limit
    Q_dot_son = 0.474*Dh_evap*(rho_v*P_sat)**(1/2)
    
    # Dryout limit
    Q_dot_dry = (2*(Dh_evap**2)*sigma*rho_v)**(3/7)*(rho_l*abs(rho_l-rho_v)*g*Dh_evap/(3*mu_l*L_e))**(1/7)*(F_r/(447*(1-F_r)))**(4/7)

    # Vapor pressure limit
    Q_dot_vap = (D_i**2)*Dh_evap*P_sat*rho_v/(64*mu_v*L_e)
    
    return (Q_dot_ent, Q_dot_boil, Q_dot_son, Q_dot_dry, Q_dot_vap)

def equ_syst(p, *param):
    
    "System to be solved using fsolve in the thermosyphon_model function, look there for variable descriptions"

    # Unknowns and parameters
    (P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_su, T_e_o_ex, P_evap_ex, R_tube_tot) = p # Unknowns
    (fluid, fluid_cd, F_r, D_i, L_e, L_c, L_eff, T_e_o_su, T_c_o_ex, R_axial, A_c_o, A_e_o, R_wc, R_we, S_T, S_L, P_cond_ex, u_wf_in, N_col, D_o, M_dot_wf, L_M, p_CO2, p_H2O, P_evap_su, M_dot_g, u_gas_fume, h_wf, foul, arrang) = param # Parameters

    if np.isnan(T_e_o_ex) or np.isnan(T_we_o): # or np.isnan(h_e_o):
        print("NAN_in")
        print(T_e_o_ex)
        print(T_we_o)
    
    # heat coefficients and resistances for the condenser  

    try: 
        h_wf = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,fluid_cd)    
        flag_1_phase = 1
    except:
        flag_1_phase = 0

    if flag_1_phase == 1:        
        if arrang == 'Inline':
            (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_inline_bank(fluid_cd, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
            h_wf_ex = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,fluid_cd)
        else:
            (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_staggered_bank(fluid_cd, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
    else: 
        h_c_o = ext_conv_boil(D_o, fluid_cd, T_c_o_ex, T_wc_o, u_wf_in)        
    
    R_c_o = (1+foul)/(A_c_o*h_c_o) # condenser external thermal resistance : takes into account the fouling factor
    
    # heat coefficients and resistances for the evaporator         
    h_r_e_o = radiative_coeff(p_CO2, p_H2O, (T_e_o_su + T_e_o_ex)/2, T_we_o, L_M)
    
    if arrang == 'Inline':
        (h_c_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_inline_bank('Air', T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[0:5]
    else: 
        (h_c_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_staggered_bank('Air', T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[0:5]

    h_e_o = h_r_e_o + h_c_e_o
    R_e_o = (1+foul)/(A_e_o*h_e_o) # condenser external thermal resistance : takes into account the fouling factor
    
    # Pressure drop

    f1 = PropsSI('P', 'T', T_v_e, 'Q', 1, fluid) - P_v_e # Saturation temperature at evaporator
    f2 = Delta_P_v(fluid, D_i, L_eff, Q_dot_radial, T_v_e) - DP_v
    
    f3 = P_v_e - DP_v - P_v_c # Pressure at condenser side => Impacts condenser saturation temperaturee
    f4 = PropsSI('T', 'P', P_v_c, 'Q', 1, fluid) - T_v_c
    
    # Thermosiphon internal sides thermal resistances
    
    (R_cond_i, R_evap_i, R_pool, R_film, Re_f) = thermal_res_esdu(fluid, F_r, D_i, L_c, L_e, Q_dot_radial, T_v_e)
    
    # Heat balances
    
    # Evaporator
    f5 = (T_e_o_su - T_we_o)/R_e_o - Q_dot_axial - Q_dot_radial # Heat transfer from outside fluid to evap wall
    f6 = (T_we_o - T_we_i)/R_we - Q_dot_radial # Heat transfer through evap wall
    f7 = (T_we_i - T_v_e)/R_evap_i - Q_dot_radial # Heat transfer from evap wall to thermosyphon fluid
    
    # Condenser
    f8 = (T_v_c - T_wc_i)/R_cond_i - Q_dot_radial # Heat transfer from thermosyphon fluid to cond wall 
    f9 = (T_wc_i - T_wc_o)/R_wc - Q_dot_radial # Heat transfer through cond wall
    f10 = (T_wc_o - T_c_o_su)/R_c_o - Q_dot_axial - Q_dot_radial # Heat transfer from cond wall to outside fluid 

    # Heat transfer through thermosyphon wall in the axial direction
    f11 = (T_we_i+T_we_o)/2 - T_we # average wall temperature of evaporator section
    f12 = (T_wc_i+T_wc_o)/2 - T_wc # average wall temperature of condenser section
    
    f13 = (T_we - T_wc)/R_axial - Q_dot_axial
    
    if flag_1_phase == 1:
        cp_wf = PropsSI('C', 'P', P_cond_ex, 'T', T_c_o_su, fluid_cd) # [J/(kg*K)]
        f14 = T_c_o_su + (Q_dot_axial + Q_dot_radial)/(M_dot_wf*cp_wf) - T_c_o_ex
    else:
        f14 = T_c_o_su - T_c_o_ex        
    
    cp_g = PropsSI('C', 'P', P_evap_su,'T',(T_e_o_ex + T_e_o_su)/2,'Air')
    
    f15 = T_e_o_su - (Q_dot_axial + Q_dot_radial)/(M_dot_g*cp_g) - T_e_o_ex
    f16 = P_evap_su - DP_evap - P_evap_ex  
    
    R_DP = abs((T_v_e - T_v_c)/Q_dot_radial)
    A_in = (np.pi/4)*D_i**2
    
    f17 = ((R_we + R_evap_i+ R_DP + R_cond_i + R_wc)**(-1) + R_axial**(-1))**(-1) - R_tube_tot

    return (f1, f2, f3, f4, f5[0], f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17)

def thermosyphon_model(beta, D_i, D_o, F_r, fluid_cd, fluid_ev, fluid, geo, k_pipe, L_a, L_c, L_e, L_M, M_dot_g, M_dot_wf, N_col, p_CO2, p_H2O, P_cond_ex, P_evap_su, S_L, S_T, T_c_o_ex, T_e_o_su, u_gas_fume, u_wf_in, A_casing, W_casing, arrang, h_wf, foul, vect_init):

    """
    ---- Inputs : -------- 

    beta : Thermosyphon inclination angle [rad]
    D_i : Themosyphon internal diameter [m]
    D_o : Themosyphon external diameter [m]        
    F_r : Filling ratio of the thermosyphon [/]
    fluid (char) : fluid name
    geo : tube geometry
    k_pipe : Thermal conductivity of the pipe [W/(m*K)]
    L_a : Length of adiabatic zone [m]
    L_c : Length of condensation zone [m]
    L_e : Length of evaporation zone [m]
    L_M : Mean beam length of a tube in the configuration [m]
    M_dot_g : Mass flowrate of fume gases in the evaporator [kg/s]
    M_dot_wf : Mass flowrate of working fluid in the condenser [kg/s]
    N_col : Number of tube columns [/]
    p_CO2 : Partial pressure of CO2 in the fume gases [/]
    p_H2O : Partial pressure of H2O in the fume gases [/]
    P_cond_su : Pressure of workingg fluid at the supply of the tube in the condenser side [/]
    P_evap_su : Pressure of fume gases at the supply of the tube in the evaporator side [/]
    S_L : Longitudinal pitch [m]
    S_T : Tranversal pitch [m]
    T_c_o : Temperature of the fluid supply outside the condenser [K]
    T_e_o : Temperature of the fluid supply outside the evaporator (fume gases) [K]
    u_gas_fume : Flow speed of fume gases in the evaporator [m/s]
    u_wf_in : Flow speed of oil in the condenser [m/s]
    A_casing : Casing area [m^2]
    W_casing : Casing width [m]
    arrang : Tube arrangement [-]
    h_wf : Inlet working fluid enthalpy [J/kg]
    foul : Fouling factor [-]
    vect_init : Vector with initial conditions of the equation system
    
    ---- Outputs : --------
    
    Q_dot : Heat transfer rate through the thermpsyphon [W]
    T_we_o : Surface temperature at the outside of the evaporator [K]
    P_v_e : Pressure in the evaporation zone [Pa] 
    T_v_e : Temperature in the evaporation zone [K]
    m_dot_v : Gas mass flowrate in the thermosyphon [kg/s]
    P_cond_ex : Pressure of oil after the tube [Pa]
    P_evap_ex : Pressure of fume gases after the tube [Pa]
    
    (Q_dot : Heat transfer rate through the thermpsyphon [W]
    
    Q_dot_ent : Entrainment limit [W]
    Q_dot_boil : Boiling limit [W]
    Q_dot_son : Sonic limit [W]
    Q_dot_vap : Vapor pressure limit [W]
    (Look in operating_limits function description for more info)
        
    T_we_o : Surface temperature at the outside of the evaporator [K]
    T_we_i : Surface temperature at the inside of the evaporator [K]
    T_v_e : Temperature in the evaporation zone [K]
    T_wc_o : Surface temperature at the outside of the condenser [K]
    T_wc_i : Surface temperature at the outside of the condenser [K]
    T_v_c : Temperature in the condensing zone [K]
    Q_dot_axial : Heat transfer rate in the thermosyphon wall in the axial direction [K]
    R_axial : Heat transfer resistance of the thermosyphon wall in the axial direction [W/K]
    P_v_e : Pressure in the evaporation zone [Pa]
    P_v_c : Pressure in the condensing zone [Pa]
    m_dot_v : Gas mass flowrate in the thermosyphon [kg/s]
    V_dot_v : Gas volume flowrate in the thermosyphon [m^3/s]
    V_v : Gas velocity in the thermosyphon [m/s]
    V_dot_l : Liquid volume flowrate in the thermosyphon [m^3/s])
    
    ---- Reference(s) : --------
    
    - Investigation of low Global Warming Potential working fluids for a closed two-phase thermosyphon
    Robert W. MacGregor, Peter A. Kew, David A. Reay
    
    - ESDU. Heat Pipes - Performance of Two phase Closed Thermosyphons. London, U.K.: Engineering Sciences Data Unit; 1981."
    
    """    
    
    "1) Design parameters"
    
    L_eff = L_a + (L_e + L_c)/2 # Thermosyphon effective length
    
    A_e_o = np.pi*L_e*D_o # evaporator ext area
    A_e_i = np.pi*L_e*D_i # evaporator int area
    A_c_o = np.pi*L_c*D_o 
        
    A_c_i = np.pi*L_c*D_i # condenser int area
    A_axial = (np.pi/4)*(D_o**2 - D_i**2) # tube annuli cross section area
    
    "2) Determine thermal resistances"
    
    R_we = np.log(D_o/D_i)/(2*np.pi*L_e*k_pipe) # evaporator wall thermal resistance
    R_wc = np.log(D_o/D_i)/(2*np.pi*L_c*k_pipe) # condenser wall thermal resistance
    R_axial = L_eff/(k_pipe*A_axial) # tube axial thermal resistance
    
    "3) System of equations : determine temperatures and heat rates"
    
    x_init = vect_init
    syst_param = (fluid, fluid_cd, F_r, D_i, L_e, L_c, L_eff, T_e_o_su, T_c_o_ex, R_axial, A_c_o, A_e_o, R_wc, R_we, S_T, S_L, P_cond_ex, u_wf_in, N_col, D_o, M_dot_wf, L_M, p_CO2, p_H2O, P_evap_su, M_dot_g, u_gas_fume, h_wf, foul, arrang)

    (sol) = fsolve(equ_syst, x0 = x_init, args = syst_param)
    
    (P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_su, T_e_o_ex, P_evap_ex, R_tube_tot) = sol

    # Compute air velocities
    V_max_air, a_gas_fume = external_flow_inline_bank('Air', T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[4:6]
    
    "4) Results computation"
    
    # Phase change heats at both thermosyphon ends
    Dh_evap_e = PropsSI('H', 'P', P_v_e, 'Q', 1, fluid)-PropsSI('H', 'P', P_v_e, 'Q', 0, fluid)
    Dh_evap_c = PropsSI('H', 'P', P_v_c, 'Q', 1, fluid)-PropsSI('H', 'P', P_v_c, 'Q', 0, fluid)
    
    # heat transfer coefficients
    try: 
        h_wf = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,fluid_cd)    
        flag_1_phase = 1
    except:
        flag_1_phase = 0
    
    if flag_1_phase == 1:        
        if arrang == 'Inline':
            (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_inline_bank(fluid_cd, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
            h_wf_ex = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,fluid_cd)
        else:
            (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_staggered_bank(fluid_cd, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
    else: 
        h_c_o = ext_conv_boil(D_o, fluid_cd, T_c_o_ex, T_wc_o, u_wf_in)        
        
    if arrang == 'Inline':
        (h_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_inline_bank(fluid_ev, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[0:5]
    else: 
        (h_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_staggered_bank(fluid_ev, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[0:5]
    
    R_c_o = 1/(A_c_o*h_c_o)
    R_e_o = 1/(A_e_o*h_e_o)

    # Results
    
    Q_dot = Q_dot_radial + Q_dot_axial
    
    m_dot_v = Q_dot_radial/Dh_evap_e # gas mass flow rate
    V_dot_v = m_dot_v/PropsSI('D', 'P', P_v_e, 'Q', 1, fluid) # gas volumic flow rate
    V_v = V_dot_v/((np.pi*D_i**2)/4) # gas velocity in thermosyphon
    
    m_dot_l = Q_dot_radial/Dh_evap_c # gas mass flow rate
    V_dot_l = m_dot_l/PropsSI('D', 'P', P_v_c, 'Q', 0, fluid) # gas volumic flow rate
    V_l = V_dot_l/((np.pi*D_i**2)/4) # gas velocity in thermosyphon
    
    "5) Operating limits"
    (Q_dot_ent, Q_dot_boil, Q_dot_son, Q_dot_dry, Q_dot_vap) = operating_limits(fluid, F_r, D_i, L_e, beta, geo, T_v_e)

    return (Q_dot, T_wc_o, T_we_o, P_v_e, T_v_e, m_dot_v, P_evap_ex, V_max_air, a_gas_fume, R_tube_tot, R_e_o, R_c_o)  
