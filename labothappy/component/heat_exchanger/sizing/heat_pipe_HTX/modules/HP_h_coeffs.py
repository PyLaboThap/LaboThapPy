# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 13:31:50 2023

@author: Basile
"""

import numpy as np
import CoolProp
import math
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
from scipy.interpolate import interp1d, RectBivariateSpline

"""
Functions in this library are used for computations realted to the external heat transfer coefficients between the heat pipe and its external fluids. 
    - radiative_coeff : Radiative heat transfer coefficient between fume gases and the heat pipe
    - euler_coeff : Euler coefficient to compute the pressure losses (used in external_flow_inline_bank)
    - euler_coeff_stag : Euler coefficient to compute the pressure losses (used in external_flow_staggered_bank)
    - external_flow_inline_bank : Convective heat transfer coefficient between fume gases and the heat pipe + Pressure loss of the fume gas when passing through the tube
    - external_flow_staggered_bank : Convective heat transfer coefficient between fume gases and the heat pipe + Pressure loss of the fume gas when passing through the tube
    - film_boiling : Used in boiling_curve
    - pool_boiling : Used in boiling_curve
    - boiling_curve : Convective coefficient curve between heat pipe and a fluid boiling in a pool
    - ext_conv_boil : Convective coefficient between heat pipe and a fluid boiling with a forced convective 
"""

def radiative_coeff(p_CO2,p_H2O,T_g,T_w,L):
    """
    ---- Inputs : -------- 
        
    p_CO2 : CO2 partial pressure in fume gases [/]
    p_H2O : H2O partial pressure in fume gases [/]
    T_g : Mean fume gas temperature [K]
    T_w : Mean thermosyphon evaporator wall temperature [K]
    L : Mean beam length [m]
        
    ---- Outputs : --------
    
    h_r : Radiative heat transfer coefficient between fume gases and evaporator wall [W/(m^2 * K)]
    
    ---- Reference(s) : --------
    
    ?
    
    """ 
    
    sigma = 5.67*1e-8 # W/(m^2 * K^4) Stefan-Boltzmann constant
    
    # Emissivity
    pL = (p_CO2+p_H2O)*L
    p_ratio = p_CO2 / p_H2O
    
    # Interpolation of a_i parameters : values in Long's code
    p_ratio_val = np.array([0,0.5,1,2,3,1e9])
    T_g_val = np.array([0, 500, 1000, 1500, 2000])
    
    a_0_val = np.array([[2.007, 2.137, 2.266, 2.395, 2.41],
                    [2.436, 2.506, 2.575, 2.645, 2.65],
                    [2.455, 2.532, 2.609, 2.686, 2.703],
                    [2.474, 2.556, 2.637, 2.718, 2.748],
                    [2.478, 2.561, 2.643, 2.726, 2.758],
                    [2.382, 2.491, 2.599, 2.708, 2.771]])
    
    a_1_val = np.array([[0.082 , 0.1281, 0.1742, 0.2203, 0.2602],
                    [0.154 , 0.2166, 0.2792, 0.3418, 0.4279],
                    [0.1497, 0.2148, 0.2799, 0.345 , 0.444],
                    [0.1397, 0.206 , 0.2723, 0.3386, 0.4464],
                    [0.1435, 0.2075, 0.2715, 0.3355, 0.4372],
                    [0.1107, 0.2061, 0.3015, 0.3969, 0.5099]])
    
    a_2_val = np.array([[-0.0304, -0.0347, -0.039 , -0.0433, -0.0651],
                    [-0.0574, -0.0611, -0.0648, -0.0685, -0.0674],
                    [-0.0603, -0.0674, -0.0745, -0.0816, -0.0859],
                    [-0.0432, -0.0618, -0.0804, -0.099 , -0.1086],
                    [-0.0486, -0.0651, -0.0816, -0.0981, -0.1122],
                    [-0.0265, -0.0613, -0.0961, -0.1309, -0.1646]])
    
    a_3_val = np.array([[0.00076, 0.00238, 0.004  , 0.00562, -0.00155],
                    [0.0137 , 0.0077 , 0.0017 ,	-0.0043, -0.012],
                    [0.006  , 0.0027 , -0.0006,	-0.0039, -0.0135],
                    [0.015  , 0.009  , 0.003  ,	-0.003 , -0.0139],
                    [0.0066 , 0.0059 , 0.0052 ,	0.0045 , -0.0065],
                    [0.03324, 0.02257, 0.0119 ,	0.00123, -0.0165]])

    f_0 = RectBivariateSpline(p_ratio_val, T_g_val, a_0_val)
    f_1 = RectBivariateSpline(p_ratio_val, T_g_val, a_1_val)
    f_2 = RectBivariateSpline(p_ratio_val, T_g_val, a_2_val)
    f_3 = RectBivariateSpline(p_ratio_val, T_g_val, a_3_val)

    # f_0 = RectBivariateSpline(T_g_val, p_ratio_val,a_0_val)
    # f_1 = RectBivariateSpline(T_g_val, p_ratio_val,a_1_val)
    # f_2 = RectBivariateSpline(T_g_val, p_ratio_val,a_2_val)
    # f_3 = RectBivariateSpline(T_g_val, p_ratio_val,a_3_val)

    a_0 = f_0(p_ratio,T_g)
    a_1 = f_1(p_ratio,T_g)
    a_2 = f_2(p_ratio,T_g)
    a_3 = f_3(p_ratio,T_g)

    # a_0 = f_0(T_g,p_ratio)
    # a_1 = f_1(T_g,p_ratio)
    # a_2 = f_2(T_g,p_ratio)
    # a_3 = f_3(T_g,p_ratio)
    
    epsilon = 10**(a_0+a_1*np.log10(pL)+a_2*(np.log10(pL))**2+a_3*(np.log10(pL))**3)/T_g
    
    # Absorptivity
    alpha = epsilon*(T_g/T_w)**0.5  
    
    # Heat transfer coefficient
    h_r = sigma*(epsilon*T_g**4-alpha*T_w**4)/(T_g-T_w)
    
    if h_r < 0 :
        h_r = [0]

    if np.isnan(h_r):
        h_r = [0]
    
    if np.isinf(h_r):
        h_r = [0]
    
    if T_w < (epsilon**4 * T_g)/(alpha**4):
        h_r = [0]
    
    return h_r[0]

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

def external_flow_inline_bank(fluid, T_in, T_out, T_w, P_in, u_in, N_col, D_out, S_T, S_L):
    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    T_out : Exhaust fluid temperature [K]
    T_w : Fluid wall temperature at the thermosiphon tube [K]
    P_in : Supply fluid pressure [Pa]
    u_in : Supply fluid velocity [m/s]
    N_col : Number of tube rows [/]
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
    
    # Pitch to diameter ratio
    a = S_T/D_out
    b = S_L/D_out
    
    # Fluid properties
    T_m = (T_in + T_out)/2 # mean temperature of fluid, used to calculate fluid properties
    
    if fluid == 'Air':
        (rho, mu, k, Pr) = PropsSI((('D','V','L','PRANDTL')), 'T', T_m, 'P', P_in, fluid) # flow prop (density, viscosity, thermal conductivity, Prandtl Number)
        (rho_w, mu_w, k_w, Pr_w) = PropsSI((('D','V','L','PRANDTL')), 'T', T_w, 'P', P_in, fluid) # near wall prop  
    elif fluid == 'Oil':
        T_f_val =   np.array([100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370]) # °C
        rho_f_val = np.array([955,948,941,934,928,921,914,907,899,892,885,878,870,863,856,848,840,832,825,817,809,800,792,783,775,766,757,748]) # kg/m^3
        cp_f_val =  np.array([1.84,1.87,1.91,1.94,1.98,2.01,2.05,2.09,2.12,2.16,2.19,2.23,2.27,2.30,2.34,2.38,2.42,2.45,2.49,2.53,2.57,2.61,2.65,2.69,2.73,2.77,2.81,2.85])*1e3 # J/(kg*K)
        k_f_val =   np.array([0.1135,0.1128,0.1121,0.1114,0.1107,0.1099,0.1091,0.1083,0.1074,0.1065,0.1056,0.1046,0.1036,0.1026,0.1015,0.1004,0.0993,0.0982,0.0970,0.0958,0.0946,0.0933,0.0920,0.0906,0.0893,0.0879,0.0865,0.0850]) # W/(m*K)
        mu_f_val =  np.array([3.60,2.92,2.42,2.05,1.75,1.52,1.33,1.18,1.06,0.95,0.86,0.784,0.718,0.661,0.611,0.567,0.529,0.495,0.464,0.437,0.413,0.391,0.371, 0.353,0.336,0.321,0.308,0.295])*1e-3 # Pa*s
                
        f_rho = interp1d(T_f_val, rho_f_val)
        f_cp = interp1d(T_f_val, cp_f_val)
        f_k = interp1d(T_f_val, k_f_val)
        f_mu = interp1d(T_f_val, mu_f_val)    
        
        # Oil supply properties
        rho = f_rho(T_in-273.15)
        cp = f_cp(T_in-273.15)
        k = f_k(T_in-273.15)
        mu = f_mu(T_in-273.15)
        Pr = (mu*cp)/k

        # Wall mean properties        
        if T_w-273.15 > 370:
            rho_w = f_rho(370)
            cp_w = f_cp(370)
            k_w = f_k(370)
            mu_w = f_mu(370)
            Pr_w = (mu_w*cp_w)/k_w
        else:
            rho_w = f_rho(T_w-273.15)
            cp_w = f_cp(T_w-273.15)
            k_w = f_k(T_w-273.15)
            mu_w = f_mu(T_w-273.15)
            Pr_w = (mu_w*cp_w)/k_w
    else:   
        (rho, mu, k, Pr) = PropsSI((('D','V','L','PRANDTL')), 'T', T_m, 'P', P_in, fluid) # flow prop (density, viscosity, thermal conductivity, Prandtl Number)
        try:
            (rho_w, mu_w, k_w, Pr_w) = PropsSI((('D','V','L','PRANDTL')), 'T', T_w, 'P', P_in, fluid) # near wall prop
        except:
            (rho_w, mu_w, k_w, Pr_w) = PropsSI((('D','V','L','PRANDTL')), 'T', T_w-1, 'P', P_in, fluid) # near wall prop

    # Reynolds number
    V_max = S_T/(S_T - D_out) * u_in
    Re_max = rho*V_max*D_out/mu
        
    a_gas = (gamma*R*T_in)**(1/2)
    
    if V_max >= a_gas:
        print('Supersonic regime in the tubes')
        print(V_max,'>',a_gas)
        
    # Various warnings
    if (Re_max>2E6):
        print('Re is out of range for ExternalFlow_Inline_Bank.  The maximum value is 2E6 while the value supplied is %.2f',Re_max)
    if (Re_max<30):
        print('Re is out of range for ExternalFlow_Inline_Bank.  The minimum value is 30 while the value supplied is %.2f',Re_max)
    if (Pr<0.5):
        print('The range of Prandtl number in ExternalFlow_Inline_Bank should be greater than 0.5. The value is %.2f',Pr)
    if (Pr>500):
        print('The range of Prandtl number in ExternalFlow_Inline_Bank should be less than 500. The value is %.2f',Pr)
    if(b<1.25):
        print('The longitudinal pitch, S_L/D is out of range for ExternalFlow_Inline_Bank. The minimum value is 1.25 while the value supplied is %.2f',b)
    if(b>2.5):
        print('The longitudinal pitch, S_L/D is out of range for ExternalFlow_Inline_Bank. The maximum value is 2.5 while the value supplied is %.2f',b)
        
    # C and m parameters determination
    if(10<Re_max) and (Re_max<100): 
        C = 0.80
        m = 0.40
 
    if(100<Re_max) and (Re_max<1000):
        Nu = 0.3 + (0.62*Re_max**(1/2)*Pr**(1/3))/(1+(0.4/Pr)**(2/3))**(1/4)*(1+(Re_max/282000)**(5/8))**(4/5) # "Churchill and Bernstein- valid for all Re tested experimentally and 0.2 < Pr ---FLUID properties are evaluated at film temperature"
        
        Pe = Pr*Re_max
        
        if (Pe<=0.2):  
            print('The Peclet number (Pr*Re) in External_Flow_Cylinder should be greater than 0.2. The value is %.2f', Pe)
            
        C_d = 10/Re_max**(2/3) + 1.0 # Drag coefficient
        
    if(1000<Re_max) and (Re_max<2e5):
        if(S_T/S_L>0.7):
            C = 0.27
            m = 0.63
        else:
            print('S_T/S_L<0.7. For values of S_T/S_L< 0.7, experimental results show that heat transfer is inefficient and aligned tubes should not be used.')

    if(2e5<Re_max) and (Re_max<2e6):
        C = 0.021
        m = 0.84
        
    if (N_col<20) and ((Re_max<1e2) or (Re_max>1e3)):	# Apply a correction factor when there is less than 20 rows of tube as thee first rows influence significantly the flow on the following and thus influencing the Nusselt number
        C_2_val = [0.7,0.8,0.86,0.9,0.92,0.95,0.97,0.98,0.99]
        N_L_val = [1,2,3,4,5,7,10,13,16]
        
        f_C2 = interp1d(N_L_val,C_2_val)
        C_2 = f_C2(N_col)
        
        Nu = C_2 * C * Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)
     
    if (N_col>=20) and ((Re_max<1e2) or (Re_max>1e3)):
        Nu=C*Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)

    # pressure drop 
    Eu = euler_coeff(a,b,Re_max) # Euler number	
    DP_row = (1/2)*Eu*rho*V_max**2 # average pressure drop across 1 row
    h = Nu*k/D_out 
    
    # print('V_max',V_max)
    # print('DP',DP_row)
    
    return (h, DP_row, Nu, Re_max, V_max, a_gas)

def external_flow_staggered_bank(fluid, T_in, T_out, T_w, P_in, u_in, N_col, D_out, S_T, S_L):

    """
    ---- Inputs : -------- 
        
    fluid (char) : fluid name
    T_in : Supply Fluid temperature [K]
    T_out : Exhaust fluid temperature [K]
    T_w : Fluid wall temperature at the thermosiphon tube [K]
    P_in : Supply fluid pressure [Pa]
    u_in : Supply fluid velocity [m/s]
    N_col : Number of tube rows [/]
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
    
    # pitch to diameter ratios
    a=S_T/D_out
    b=S_L/D_out
    
    # mean temperature of fluid, used to calculate fluid properties
    T_m	=(T_in+T_out)/2	

    if fluid == 'Air':
        (rho, mu, k, Pr) = PropsSI((('D','V','L','PRANDTL')), 'T', T_m, 'P', P_in, fluid) # flow prop (density, viscosity, thermal conductivity)
        (rho_w, mu_w, k_w, Pr_w) = PropsSI((('D','V','L','PRANDTL')), 'T', T_w, 'P', P_in, fluid) # near wall prop
    elif fluid == 'Oil':
        T_f_val =   np.array([100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370]) # °C
        rho_f_val = np.array([955,948,941,934,928,921,914,907,899,892,885,878,870,863,856,848,840,832,825,817,809,800,792,783,775,766,757,748]) # kg/m^3
        cp_f_val =  np.array([1.84,1.87,1.91,1.94,1.98,2.01,2.05,2.09,2.12,2.16,2.19,2.23,2.27,2.30,2.34,2.38,2.42,2.45,2.49,2.53,2.57,2.61,2.65,2.69,2.73,2.77,2.81,2.85])*1e3 # J/(kg*K)
        k_f_val =   np.array([0.1135,0.1128,0.1121,0.1114,0.1107,0.1099,0.1091,0.1083,0.1074,0.1065,0.1056,0.1046,0.1036,0.1026,0.1015,0.1004,0.0993,0.0982,0.0970,0.0958,0.0946,0.0933,0.0920,0.0906,0.0893,0.0879,0.0865,0.0850]) # W/(m*K)
        mu_f_val =  np.array([3.60,2.92,2.42,2.05,1.75,1.52,1.33,1.18,1.06,0.95,0.86,0.784,0.718,0.661,0.611,0.567,0.529,0.495,0.464,0.437,0.413,0.391,0.371, 0.353,0.336,0.321,0.308,0.295])*1e-3 # Pa*s
                
        f_rho = interp1d(T_f_val, rho_f_val)
        f_cp = interp1d(T_f_val, cp_f_val)
        f_k = interp1d(T_f_val, k_f_val)
        f_mu = interp1d(T_f_val, mu_f_val)    
        
        # Oil supply properties
        rho = f_rho(T_in-273.15)
        cp = f_cp(T_in-273.15)
        k = f_k(T_in-273.15)
        mu = f_mu(T_in-273.15)
        Pr = (mu*cp)/k

        # Wall mean properties        
        if T_w-273.15 > 370:
            rho_w = f_rho(370)
            cp_w = f_cp(370)
            k_w = f_k(370)
            mu_w = f_mu(370)
            Pr_w = (mu_w*cp_w)/k_w
        else:
            rho_w = f_rho(T_w-273.15)
            cp_w = f_cp(T_w-273.15)
            k_w = f_k(T_w-273.15)
            mu_w = f_mu(T_w-273.15)
            Pr_w = (mu_w*cp_w)/k_w
        
    # Reynolds number
    V_max = S_T/(S_T - D_out) * u_in
    Re_max = rho*V_max*D_out/mu

    a_gas = (gamma*R*T_in)**(1/2)

    S_D	= (S_L**2+(S_T/2)**2)**(1/2)
    
    if(S_D< (S_T + D_out)/2):
        V_max = S_T/(2*(S_D-D_out))*u_in
    else:
        V_max = S_T/(S_T-D_out)*u_in
  
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
            
        C_d = 10/Re_max**(2/3) + 1.0 # Drag coefficient
        
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
            
    if (N_col<20) and ((Re_max<1e2) or (Re_max>1e3)):	# Apply a correction factor when there is less than 20 rows of tube as thee first rows influence significantly the flow on the following and thus influencing the Nusselt number
        C_2_val = [0.64,0.76,0.84,0.89,0.92,0.95,0.97,0.98,0.99]
        N_L_val = [1,2,3,4,5,7,10,13,16]
        
        f_C2 = interp1d(N_L_val,C_2_val)
        C_2 = f_C2(N_col)
        
        Nu = C_2 * C * Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)
        
    if (N_col>=20) and ((Re_max<1e2) or (Re_max>1e3)):
        Nu=C*Re_max**m*Pr**0.36*(Pr/Pr_w)**(1/4)

    # pressure drop 
    Eu = euler_coeff_stag(a,b,Re_max) # Euler number	
    DP_row = (1/2)*Eu*rho*V_max**2 # average pressure drop across 1 row
    h = Nu*k/D_out 
        
    return (h, DP_row, Nu, Re, V_max, a_gas)

def pool_boiling(fluid, T_sat, T_tube):
    """
    ---- Inputs : -------- 
    
    Corr : Correlation used for the heat transfer coefficient (name)
    fluid (char) : fluid name (oil)
    T_sat : External fluid saturation temperature [K]
    T_tube : Tube temperature [K]
    
    ---- Outputs : --------
    
    h_pool : Heat transfer coefficient at the condenser [W/(m^2*K)]
    
    ---- Reference(s) : --------
    
    Foundations of heat transfer : Frank P. Incropera
    
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
    
    D_out : Outer Tube diameter [m]
    fluid (char) : fluid name (oil)
    T_sat : External fluid saturation temperature [K]
    P_sat = Saturation Pressure [Pa]
    q_flux : Condenser design heat rate [W]
    
    ---- Outputs : --------
    
    h_film : Heat transfer coefficient at the condenser [W/(m^2*K)]
    
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
    
    # Initial guess for x
    initial_guess = 1.0
    
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
    sigma = PropsSI('I','T',T_sat,'Q',0.5,fluid) # surface tension of liquid vapor equilibrium
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
        else:
            q_final[i] = q_crit
    
    h_final = q_final/DT
    
    DT = np.concatenate((np.array([0]),DT),axis = 0)
    h_final = np.concatenate((np.array([0]),h_final),axis = 0)
    
    return h_final, DT

def ext_conv_boil(D_out, fluid, T_sat, T_tube, V_flow):
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
