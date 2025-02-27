# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 09:52:51 2024

@author: Basile
"""

import numpy as np
import CoolProp
import math
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
from scipy.interpolate import interp1d, interp2d

"""
Functions in this library are used for computations realted to the external heat transfer coefficients between the heat pipe and its external fluids. 
    - radiative_coeff : Radiative heat transfer coefficient between fume gases and the heat pipe
    - euler_coeff : Euler coefficient to compute the pressure losses (used in external_flow_inline_bank)
    - external_flow_inline_bank : Convective heat transfer coefficient between fume gases and the heat pipe + Pressure loss of the fume gas when passing through the tube
    - h_cond_Th66 : Convective coefficient between heat pipe and Therminol 66 without phase change (Shell-and-tube)
    - pool_boiling : Convective coefficient between heat pipe and a fluid boiling in a pool 
"""

def euler_coeff(a,b,Re):
    """
    ---- Inputs : -------- 
        
    a : Transversal pitch to diameter ratio of tube bank [/]
    b : Longitudinal pitch to diameter ratio of tube bank [/]
    Re : Reynolds number of flow across tube bank [/]
    
    ---- Outputs : --------
    
    Eu : Euler number [/]
    
    ---- Reference(s) : --------
    
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/
    
    """
    
    # Coefficients determination
    
    if a != b:
        print("The studied tube bank is not square, Euler number can't be computed.")
        return 0
    else:
        if a != 1.25 and a != 1.5 and a != 2:
            a_arr = np.array([1.25,1.5,2])
            a = min(a_arr, key=lambda x:abs(x-a))
            
        if a == 1.25:
            if (Re >= 3 and Re <= 2e3):
                c = np.array([0.272, 0.207*1e3, 0.102*1e3, -0.286*1e4, 0])
                
            elif (Re > 2e3 and Re <= 2e6):
                c = np.array([0.267, 0.249*1e4, -0.927*1e7, 0.1*1e11, 0])
                            
            else:
                # print("Re value is not valid for Euler number computation")
                return 0
        elif a == 1.5:
            if (Re >= 3 and Re <= 2e3):
                c = np.array([0.263, 0.867*1e2, -0.202, 0, 0])
                
            elif (Re > 2e3 and Re <= 2e6):
                c = np.array([0.235, 0.197*1e4, -0.124*1e8, 0.312*1e11, -0.274*1e14])
                         
            else:
                # print("Re value is not valid for Euler number computation",Re)
                return 0
        elif a == 2:
            if (Re >= 7 and Re <= 800):
                c = np.array([0.188, 0.566*1e2, -0.646*1e3, 0.601*1e4, -0.183*1e5])
                
            elif (Re > 800 and Re <= 2e6):
                c = np.array([0.247, -0.595, 0.15, -0.137, 0.396])
                        
            else:
                # print("Re value is not valid for Euler number computation")
                return 0
        else:
            print("S_T/D_out ratio is different from 1.25, 1.5 and 2. No data is thus available to compute Euler number")
            return 0
        
        # Euler number computation
        Eu = 0
        
        for i in range(5):
            Eu = Eu + c[i]/Re**i
        
        # print('Eu' ,Eu)
    return Eu

def euler_coeff_stag(a,b,Re):
    """
    ---- Inputs : -------- 
        
    a : Transversal pitch to diameter ratio of tube bank [/]
    b : Longitudinal pitch to diameter ratio of tube bank [/]
    Re : Reynolds number of flow across tube bank [/]
    
    ---- Outputs : --------
    
    Eu : Euler number [/]
    
    ---- Reference(s) : --------
    
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/
    
    """

    # Coefficients determination
    if b != 1.25 and b != 1.5 and b != 2 and b != 2.5:
        arr_b = np.array([1.25,1.5,2,2.5])
        b = min(arr_b, key=lambda x:abs(x-b))
    if b == 1.25:
        if (Re >= 3 and Re <= 1e3):
            c = np.array([0.795, 0.247*1e3, 0.335*1e3, -0.155*1e4, 0.241*1e4])
            
        elif (Re > 1e3 and Re <= 2e6):
            c = np.array([0.245, 0.339*1e4, -0.984*1e7, 0.132*1e11, -0.599*1e13])
                        
        else:
            # print("Re value is not valid for Euler number computation")
            return 0
        
    elif b == 1.5:
        if (Re >= 3 and Re <= 1e3):
            c = np.array([0.683, 0.111*1e3, -0.973*1e2, -0.426*1e3, 0.574*1e3])
            
        elif (Re > 1e3 and Re <= 2e6):
            c = np.array([0.203, 0.248*1e4, -0.758*1e7, 0.104*1e11, -0.482*1e13])
                     
        else:
            # print("Re value is not valid for Euler number computation")
            return 0
        
    elif b == 2:
        if (Re >= 7 and Re <= 1e2):
            c = np.array([0.713, 0.448*1e2, -0.126*1e3, -0.582*1e3, 0])
            
        elif (Re > 1e2 and Re <= 1e4):
            c = np.array([0.343, 0.303*1e3, -0.717*1e5, 0.88*1e7, -0.38*1e9])
            
        elif (Re > 1e4 and Re <= 1e6):
            c = np.array([0.162, 0.181*1e4, 0.792*1e8, -0.165*1e13, 0.872*1e16])  
            
        else:
            # print("Re value is not valid for Euler number computation")
            return 0
        
    elif b == 2.5:
        if (Re >= 1e2 and Re <= 5e3):
            c = np.array([0.33, 0.989*1e2, -0.148*1e5, 0.192*1e7, -0.862*1e8])
            
        elif (Re > 5e3 and Re <= 2e6):
            c = np.array([0.119, 0.498*1e4, -0.507*1e8, 0.251*1e12, -0.463*1e15])
                     
        else:
            # print("Re value is not valid for Euler number computation")
            return 0
        
    else:
        print("S_L/D_out ratio is different from 1.25, 1.5 and 2. No data is thus available to compute Euler number")
        return 0
    
    # Euler number computation
    Eu = 0
    
    for i in range(5):
        Eu = Eu + c[i]/Re**i
        
    return Eu

def DP_external_flow_inline_bank(fluid, T_in, P_mean, u_flow, D_out, S_T, S_L):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    P_mean : Supply fluid pressure [Pa]
    u_flow : Supply fluid velocity [m/s]
    D_out : Outer diameter of a tube [m]
    S_T : Tube pitch : transversal [m] (écartement)
    S_L : Tube pitch : longitudinal [m]
    
    ---- Outputs : --------
    
    DP : Pressure loss caused by the evaporator [Pa]
    
    ---- Reference(s) : --------
    
    Nusselt number computation : Frank P Incropera : - Foundations of heat transfer 
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/ (Pressure drop and V_max computations)
    
    """
    
    # Pitch to diameter ratio
    a = S_T/D_out
    b = S_L/D_out
    
    (rho, mu) = PropsSI((('D','V')), 'T', T_in, 'P', P_mean, fluid) # flow prop (density, viscosity, thermal conductivity)

    # Reynolds number
    V_max = S_T/(S_T - D_out) * u_flow
    Re_max = rho*V_max*D_out/mu

    # pressure drop 
    Eu = euler_coeff(a,b,Re_max) # Euler number	
    DP_row = (1/2)*Eu*rho*V_max**2 # average pressure drop across 1 row
        
    return DP_row

def DP_external_flow_staggered_bank(fluid, T_in, P_mean, u_flow, D_out, S_T, S_L):

    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    P_mean : Supply fluid pressure [Pa]
    u_flow : Supply fluid velocity [m/s]
    D_out : Outer diameter of a tube [m]
    S_T : Tube pitch : transversal [m] (écartement)
    S_L : Tube pitch : longitudinal [m]
    
    ---- Outputs : --------
    
    DP : Pressure loss caused by the evaporator [Pa]
    
    ---- Reference(s) : --------
    
    Nusselt number computation : Frank P Incropera : - Foundations of heat transfer 
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/ (Pressure drop and V_max computations)
    
    """
    
    # pitch to diameter ratios
    a=S_T/D_out
    b=S_L/D_out

    (rho, mu) = PropsSI((('D','V')), 'T', T_in, 'P', P_mean, fluid) # flow prop (density, viscosity, thermal conductivity)

    # Reynolds number
    V_max = S_T/(S_T - D_out) * u_flow
    S_D	= (S_L**2+(S_T/2)**2)**(1/2)
    
    if(S_D< (S_T + D_out)/2):
        V_max = S_T/(2*(S_D-D_out))*u_flow
    else:
        V_max = S_T/(S_T-D_out)*u_flow
        
    Re_max = rho*V_max*D_out/mu
  
    # pressure drop 
    Eu = euler_coeff_stag(a,b,Re_max) # Euler number	
    DP_row = (1/2)*Eu*rho*V_max**2 # average pressure drop across 1 row
        
    return DP_row 

def tube_bank_DP(fluid, T_in, P_mean, u_flow, D_out, S_T, S_L, arrang):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    P_in : Supply fluid pressure [Pa]
    u_flow : Supply fluid velocity [m/s]
    D_out : Outer diameter of a tube [m]
    S_T : Tube pitch : transversal [m] (écartement)
    S_L : Tube pitch : longitudinal [m]
    arrangement : Inline (I) or Staggered (S)
        
    ---- Outputs : --------
    
    DP : Pressure loss caused by the evaporator [Pa]
    
    ---- Reference(s) : --------
    
    Nusselt number computation : Frank P Incropera : - Foundations of heat transfer 
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/ (Pressure drop and V_max computations)
    
    """
    if arrang == "Inline":
        DP_row = DP_external_flow_inline_bank(fluid, T_in, P_mean, u_flow, D_out, S_T, S_L)
    elif arrang == "Staggered":
        DP_row = DP_external_flow_staggered_bank(fluid, T_in, P_mean, u_flow, D_out, S_T, S_L)
    else:
        print("Tube arrangement not specified by either I (Inline) or S (Staggered)")
    
    return DP_row


def external_flow_inline_bank(fluid, T_in, T_out, T_w, P_mean, u_flow, N_row, D_out, S_T, S_L):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    T_out : Exhaust fluid temperature [K]
    T_w : Fluid wall temperature at the thermosiphon tube [K]
    P_mean : Supply fluid pressure [Pa]
    u_flow : Supply fluid velocity [m/s]
    N_row : Number of tube rows [/]
    D_out : Outer diameter of a tube [m]
    S_T : Tube pitch : transversal [m] (écartement)
    S_L : Tube pitch : longitudinal [m]
    
    ---- Outputs : --------
    
    h : Heat transfer coefficient at the evaporator (convective) [W/(m^2*K)]
    DP : Pressure loss caused by the evaporator [Pa]
    Nu : Nusselt number [/]
    Re : Reynolds number [/]
    
    ---- Reference(s) : --------
    
    Nusselt number computation : Frank P Incropera : - Foundations of heat transfer 
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/ (Pressure drop and V_max computations)
    
    """
    
    gamma = 1.4 # [/]
    R = 285.71429 # J/(kg*K)
    
    # Fluid properties
    T_m = (T_in + T_out)/2 # mean temperature of fluid, used to calculate fluid properties
    
    (rho, mu, k, Pr) = PropsSI((('D','V','L','PRANDTL')), 'T', T_m, 'P', P_mean, fluid) # flow prop (density, viscosity, thermal conductivity)
    (rho_w, mu_w, k_w, Pr_w) = PropsSI((('D','V','L','PRANDTL')), 'T', T_w, 'P', P_mean, fluid) # near wall prop

    # Reynolds number
    V_max = S_T/(S_T - D_out) * u_flow
    Re_max = rho*V_max*D_out/mu
        
    a_gas = (gamma*R*T_in)**(1/2)
    
    if V_max >= a_gas:
        print('Supersonic regime in the tubes')
        print(V_max,'>',a_gas)
          
    # C and m parameters determination
    if(10<Re_max) and (Re_max<100): 
        C = 0.80
        m = 0.40
 
    if(100<Re_max) and (Re_max<1000):
        Nu = 0.3 + (0.62*Re_max**(1/2)*Pr**(1/3))/(1+(0.4/Pr)**(2/3))**(1/4)*(1+(Re_max/282000)**(5/8))**(4/5) # "Churchill and Bernstein- valid for all Re tested experimentally and 0.2 < Pr ---FLUID properties are evaluated at film temperature"
        
        Pe = Pr*Re_max
        
        if (Pe<=0.2):  
            print('The Peclet number (Pr*Re) in External_Flow_Cylinder should be greater than 0.2. The value is %.2f', Pe)
                    
    if(1000<Re_max) and (Re_max<2e5):
        if(S_T/S_L>0.7):
            C = 0.27
            m = 0.63
        else:
            print('S_T/S_L<0.7. For values of S_T/S_L< 0.7, experimental results show that heat transfer is inefficient and aligned tubes should not be used.')

    if(2e5<Re_max) and (Re_max<2e6):
        C = 0.021
        m = 0.84
        
    if (N_row<20) and ((Re_max<1e2) or (Re_max>1e3)):	# Apply a correction factor when there is less than 20 rows of tube as thee first rows influence significantly the flow on the following and thus influencing the Nusselt number
        C_2_val = [0.7,0.8,0.86,0.9,0.92,0.95,0.97,0.98,0.99]
        N_L_val = [1,2,3,4,5,7,10,13,16]
        
        f_C2 = interp1d(N_L_val,C_2_val)
        C_2 = f_C2(N_row)
        
        Nu = C_2 * C * Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)
     
    if (N_row>=20) and ((Re_max<1e2) or (Re_max>1e3)):
        Nu=C*Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)
    
    h = Nu*k/D_out 
    
    return h

def external_flow_staggered_bank(fluid, T_in, T_out, T_w, P_mean, u_flow, N_row, D_out, S_T, S_L):

    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    T_out : Exhaust fluid temperature [K]
    T_w : Fluid wall temperature at the thermosiphon tube [K]
    P_mean : Supply fluid pressure [Pa]
    u_flow : Supply fluid velocity [m/s]
    N_row : Number of tube rows [/]
    D_out : Outer diameter of a tube [m]
    S_T : Tube pitch : transversal [m] (écartement)
    S_L : Tube pitch : longitudinal [m]
    
    ---- Outputs : --------
    
    h : Heat transfer coefficient at the evaporator (convective) [W/(m^2*K)]
    DP : Pressure loss caused by the evaporator [Pa]
    Nu : Nusselt number [/]
    Re : Reynolds number [/]
    
    ---- Reference(s) : --------
    
    Nusselt number computation : Frank P Incropera : - Foundations of heat transfer 
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/ (Pressure drop and V_max computations)
    
    """
    
    # pitch to diameter ratios
    a=S_T/D_out
    
    # mean temperature of fluid, used to calculate fluid properties
    T_m	=(T_in+T_out)/2	

    (rho, mu, k, Pr) = PropsSI((('D','V','L','PRANDTL')), 'T', T_m, 'P', P_mean, fluid) # flow prop (density, viscosity, thermal conductivity)
    (rho_w, mu_w, k_w, Pr_w) = PropsSI((('D','V','L','PRANDTL')), 'T', T_w, 'P', P_mean, fluid) # near wall prop

    # Reynolds number
    V_max = S_T/(S_T - D_out) * u_flow
    Re_max = rho*V_max*D_out/mu

    S_D	= (S_L**2+(S_T/2)**2)**(1/2)
    
    if(S_D< (S_T + D_out)/2):
        V_max = S_T/(2*(S_D-D_out))*u_flow
    else:
        V_max = S_T/(S_T-D_out)*u_flow
  
    Re = rho*V_max*D_out/mu

    # Various warnings
    if (Re_max>2E6):
        print('Re is out of range for ExternalFlow_Inline_Bank.  The maximum value is 2E6 while the value supplied is %.2f',Re_max)
    if (Re_max<30):
        print('Re is out of range for ExternalFlow_Inline_Bank.  The minimum value is 30 while the value supplied is %.2f',Re_max)
    if (Pr<0.5):
        print('The range of Prandtl number in ExternalFlow_Inline_Bank should be greater than 0.5. The value is %.2f',Pr)
    if (Pr>500):
        print('The range of Prandtl number in ExternalFlow_Inline_Bank should be less than 500. The value is %.2f',Pr)
    if(a<1.25):
        print('The transversal pitch, S_T/D is out of range for ExternalFlow_Inline_Bank. The minimum value is 1.25 while the value supplied is %.2f',a)
    if(a>2.5):
        print('The transversal pitch, S_T/D is out of range for ExternalFlow_Inline_Bank. The maximum value is 2.5 while the value supplied is %.2f',a)
   
    if(10<Re) and (Re<100): 
        C = 0.90
        m = 0.40

    if(100<Re) and (Re<1000):
        Nu = 0.3 + (0.62*Re_max**(1/2)*Pr**(1/3))/(1+(0.4/Pr)**(2/3))**(1/4)*(1+(Re_max/282000)**(5/8))**(4/5) # "Churchill and Bernstein- valid for all Re tested experimentally and 0.2 < Pr ---FLUID properties are evaluated at film temperature"
        Pe = Pr*Re_max
        
        if (Pe<=0.2):  
            print('The Peclet number (Pr*Re) in External_Flow_Cylinder should be greater than 0.2. The value is %.2f', Pe)
                    
    if(1000<Re) and (Re<2e5):
        if(S_T/S_L<2):
            C=0.35*(S_T/S_L)**(1/5)
            m=0.60
        else:
            C=0.40	
            m=0.60
            
    if(2e5<Re) and (Re<2e6):
        C = 0.021
        m = 0.84
            
    if (N_row<20) and ((Re_max<1e2) or (Re_max>1e3)):	# Apply a correction factor when there is less than 20 rows of tube as thee first rows influence significantly the flow on the following and thus influencing the Nusselt number
        C_2_val = [0.64,0.76,0.84,0.89,0.92,0.95,0.97,0.98,0.99]
        N_L_val = [1,2,3,4,5,7,10,13,16]
        
        f_C2 = interp1d(N_L_val,C_2_val)
        C_2 = f_C2(N_row)
        
        Nu = C_2 * C * Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)
        
    if (N_row>=20) and ((Re_max<1e2) or (Re_max>1e3)):
        Nu=C*Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)

    # pressure drop 
    h = Nu*k/D_out 
        
    return h

def tube_bank_htc_1P(fluid, T_in, T_out, T_w, P_mean, u_flow, N_row, D_out, S_T, S_L, arrang):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    T_out : Exhaust fluid temperature [K]
    T_w : Fluid wall temperature at the thermosiphon tube [K]
    P_in : Supply fluid pressure [Pa]
    u_flow : Supply fluid velocity [m/s]
    N_row : Number of tube rows [/]
    D_out : Outer diameter of a tube [m]
    S_T : Tube pitch : transversal [m] (écartement)
    S_L : Tube pitch : longitudinal [m]
    arrangement : Inline (I) or Staggered (S)
        
    ---- Outputs : --------
    
    h : Heat transfer coefficient at the evaporator (convective) [W/(m^2*K)]
    DP : Pressure loss caused by the evaporator [Pa]
    Nu : Nusselt number [/]
    Re : Reynolds number [/]
    
    ---- Reference(s) : --------
    
    Nusselt number computation : Frank P Incropera : - Foundations of heat transfer 
    Tube Banks, Crossflow over - Beale, Steven : https://www.thermopedia.com/cn/content/1211/ (Pressure drop and V_max computations)
    
    """
    if arrang == "Inline":
        h = external_flow_inline_bank(fluid, T_in, T_out, T_w, P_mean, u_flow, N_row, D_out, S_T, S_L)
    elif arrang == "Staggered":
        h = external_flow_staggered_bank(fluid, T_in, T_out, T_w, P_mean, u_flow, N_row, D_out, S_T, S_L)
    else:
        print("Tube arrangement not specified by either I (Inline) or S (Staggered)")
    
    return h

def ext_tube_conv_boil(D_out, fluid, T_sat, T_tube, V_flow):
    """
    ---- Inputs : -------- 
    
    D_out : Outer Tube diameter [m]
    fluid (char) : fluid name (Cyclopentane)
    T_sat : External fluid saturation temperature [K]
    T_tube : In-flow Tube temperature [K]    
    V_flow : External fluid flow speed [m/s]
    
    ---- Outputs : --------
    
    h_conv_boil : Heat transfer coefficient related to external forced convection boiling for a crossflow over a tube [W/(m^2)]
    
    ---- Reference(s) : --------
    
    Foundations of heat transfer : Frank P. Incropera
    
    """
    rho_l = PropsSI('D','T',T_sat,'Q',0,fluid) # density at saturated liquid in kg/m^3
    rho_v = PropsSI('D','T',T_sat,'Q',1,fluid) # density at saturated vapor in kg/m^3
    sigma = PropsSI('I','T',T_sat,'Q',0.5,fluid) # surface tension of liquid vapor equilibrium
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific heat
    
    # Weber number
    We_D = (rho_v*V_flow**2*D_out)/sigma
    coef_low = (1/np.pi)*(1 + (4/We_D)**(1/3))
    
    # Lienhard and Eichhorn correlations
    
    if coef_low >= (0.275/np.pi)*(rho_l/rho_v)**(1/2) + 1: # low velocity region
        q_max = rho_v*V_flow*Dh_evap*coef_low
    else: # high velocity region
        term_1 = ((rho_l/rho_v)**(3/4))/(169*np.pi)
        term_2 = ((rho_l/rho_v)**(1/2))/(19.2*np.pi*We_D**(1/3))
        q_max = rho_v*V_flow*Dh_evap*(term_1 + term_2)
    
    h_conv_boil = q_max/(T_tube - T_sat)
    
    return h_conv_boil

def ext_tube_film_condens(D_out, fluid, T_sat, T_w, V_flow):
    """
    ---- Inputs : -------- 
    
    D_out : Outer Tube diameter [m]
    fluid (char) : fluid name (Cyclopentane)
    T_sat : External fluid saturation temperature [K]
    T_w : Tube wall temperature [K]    
    V_flow : External fluid flow speed [m/s]
    
    ---- Outputs : --------
    
    h_conv_boil : Heat transfer coefficient related to external forced convection boiling for a crossflow over a tube [W/(m^2)]
    
    ---- Reference(s) : --------
    
    Foundations of heat transfer : Frank P. Incropera
    
    """
    
    C = 0.729 # [-] For tubes / 0.826 for spheres
    g = 9.81 # [m/s^2] : gravity acceleration constant    
    
    (rho_l, mu_l, k_l) = PropsSI(('D','V','L'),'T',T_sat,'Q',0,fluid) # density at saturated liquid in kg/m^3, liquid viscosity in Pa*s, liquid conductivity in W/(m*K)
    rho_v = PropsSI('D','T',T_sat,'Q',1,fluid) # density at saturated vapor in kg/m^3
    Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific heat
    
    # Weber number
    fact_1 = rho_l*g*(rho_l-rho_v)*Dh_evap*D_out**3
    fact_2 = mu_l*k_l*(T_sat - T_w)
    
    h_film_cond = C*(fact_1/fact_2)**(1/4)
    
    return h_film_cond
