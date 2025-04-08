# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:36:53 2024

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP

def s_max_kern(Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout): # Maximum flow section (m**2)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    Shell_ID        : Shell inside diameter [m]
    central_spacing : Spacing between two baffles [m]
    tube_layout     : Angle of tube layout [°]
    
    Output
    ------
    
    s_max : Maximum flow area in the tube bank [m^2]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio # Tube Pitch
    s_max = (Shell_ID / p_T)*(p_T -  Tube_OD) * central_spacing

    return s_max


def s_max(Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout): # Maximum flow section (m**2)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    Shell_ID        : Shell inside diameter [m]
    central_spacing : Spacing between two baffles [m]
    tube_layout     : Angle of tube layout [°]
    
    Output
    ------
    
    s_max : Maximum flow area in the tube bank [m^2]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio # Tube Pitch
    s_max_coef = (Shell_ID / p_T)*(p_T -  Tube_OD) * central_spacing

    if tube_layout == 0: # Squared Inline tube arrangement
    
        s_max = s_max_coef
    
    elif tube_layout == 45: # Squared Staggered tube arrangement
    
        s_max = (2/3**(1/2))*s_max_coef

    elif tube_layout == 60: # Triangular Staggered tube arrangement
    
        s_max = (2**(1/2))*s_max_coef
    
    else:
        
        print("Tube arrangement shall be either square Inline (tube_layout = 0), square Staggered (tube_layout = 45) or triangular staggered (tube_layout = 60)")
        return 0
        
    return s_max

def s_L(Baffle_cut, Shell_ID, n_tubes, Tube_OD):
    """
    Input
    -----
    
    Baffle_cut      : Proportions of shell area not blocked by the baffles [-]
    Shell_ID        : Shell inside diameter [m]
    n_tubes         : Number of tubes [-]
    Tube_OD         : Tube outside diameter [m]
    
    Output
    ------
    
    s_L : Maximum flow area in the tube bank in longitudinal shell direction [m^2]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    S_L = (Baffle_cut/100)*(np.pi * Shell_ID/4 - n_tubes*np.pi*Tube_OD**2/4)
    
    return S_L

def d_h(Tube_OD, pitch_ratio, tube_layout): # Hydraulic diameter (m)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    tube_layout     : Angle of tube layout [°]
    
    Output
    ------
    
    D_h : Flow Hydraulic diameter in tube banks [m]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio # Tube Pitch

    if tube_layout == 0 or tube_layout == 45: # Squared tube arrangement
    
        D_h = 4*( (p_T**2 - np.pi*(Tube_OD**2 / 4)) / (np.pi*Tube_OD) )

    else:
    
        D_h = 2 * (3**(1/2) / np.pi) * (p_T**2 / Tube_OD) - Tube_OD
        
    return D_h

#%%

def shell_htc_kern(m_dot, T_wall, T_in, P_in, fluid, params):
    """
    Inputs
    ----------
    
    m_dot  : Flow rate [kg/s]
    T_wall : Wall Temperature [K]
    T_in   : Inlet Temperature [K]
    P_in   : Inlet pressure [Pa]
    fluid  : fluid name [-]
    params : HTX parameters [-]

    Outputs
    -------
    
    h = shell side heat transfer coefficient [W/(m^2 * K)]    
    DP = shell side pressure drop [Pa]
    
    References
    -------
    Process Heat Transfer - D. Q. Kern
    """
    
    "1) HTC"
    
    AS = CP.AbstractState("BICUBIC&HEOS", fluid)
    
    AS.update(CP.PT_INPUTS, P_in, T_in)
    
    rho = AS.rhomass()
    mu = AS.viscosity()
    Pr = AS.Prandtl()
    k = AS.conductivity()

    AS.update(CP.PT_INPUTS, P_in, T_wall)
    
    mu_w = AS.viscosity()

    S_T = s_max_kern(params['Tube_OD'], params['pitch_ratio'], params['Shell_ID'], params['central_spacing'], params['tube_layout']) # m^2
    D_hydro = d_h(params['Tube_OD'], params['pitch_ratio'], params['tube_layout'])

    V_t = m_dot/(S_T*rho)
    
    Re = rho*V_t*(D_hydro/mu)
    
    # print(Re)    
    
    if Re < 2e3:
        Nu = 0.427 * Re**0.528 * Pr**(1/3) 
    else:
        # McAdams, if no Baffles -> HYP : mu = mu_w
        Nu = 0.36*Pr**(1/3) * Re**(0.55) * (mu/mu_w)**(0.14)  

    h = Nu*k/D_hydro
    
    # print(f"h : {h}")
    # print(f"Re : {Re}")
    # print(f"Pr : {Pr}")
        
    return h, Re, Pr

#%%

def shell_htc_donohue(m_dot, T_in, P_in, fluid, params):
    """
    Inputs
    ----------
    
    m_dot  : Flow rate [kg/s]
    T_in   : Inlet Temperature [K]
    P_in   : Inlet pressure [Pa]
    fluid  : fluid name [-]
    geom   : HTX geometry [-]

    Outputs
    -------
    
    h = shell side heat transfer coefficient [W/(m^2 * K)]    
    
    References
    -------
    Process Heat Transfer - D. Q. Kern
    """
    
    rho = PropsSI('D','T',T_in, 'P',P_in,fluid)
    mu = PropsSI('V','T',T_in, 'P',P_in,fluid)
    Pr = PropsSI('PRANDTL','T',T_in, 'P',P_in,fluid)
    k = PropsSI('L','T',T_in, 'P',P_in,fluid)
    
    S_T = s_max(params['Tube_OD'], params['pitch_ratio'], params['Shell_ID'], params['central_spacing'], params['tube_layout']) # m^2
    D_hydro = d_h(params['Tube_OD'], params['pitch_ratio'], params['tube_layout'])
    
    V_t = m_dot/(S_T*rho)
        
    S_L = s_L(params['Baffle_cut'], params['Shell_ID'], params['n_tubes'], params['Tube_OD']) # m^2
    V_L = m_dot/(S_L*rho)
    
    V_R = (V_L*V_t)**(1/2)
    
    Re = rho*V_R*(params['Tube_OD']/mu)

    # HYP : mu_in_wf = mu_w_in_wf
    Nu = 0.22*Re**(0.6)*Pr**(0.33) # * (mu_in_wf/mu_w_in_wf)**(0.14)  
    
    h = Nu*k/params['Tube_OD']   
    
    return h

#%%

def bell_delaware_coefs(layout_angle, Re):
    """

    Parameters
    ----------
    layout_angle : Tube layout angle
    Re : Reynolds number

    Returns
    -------
    a : coefficients
    b : coefficients
        
    Reference
    -------   
    https://cheguide.com/heat_exchanger_rating.html

    """
    # Define the breakpoints and corresponding 'a' and 'b' values
    breakpoints = [1, 10, 1e2, 1e3, 1e4, 1e5]
    
    if layout_angle == 0: # Squared Inline tube arrangement
        a_values = [
            [0.97, -0.667, 1.187, 0.37],
            [0.97, -0.667, 1.187, 0.37],
            [0.9, -0.631, 1.187, 0.37],
            [0.408, -0.46, 1.187, 0.37],
            [0.107, -0.266, 1.187, 0.37],
            [0.37, -0.395, 1.187, 0.37],
            #[0.37, -0.395, 1.187, 0.37]  # Extend the last value for extrapolation
        ]
        b_values = [
            [35, -1, 6.3, 0.378],
            [35, -1, 6.3, 0.378],
            [32.1, -0.963, 6.3, 0.378],
            [6.09, -0.602, 6.3, 0.378],
            [0.0815, 0.022, 6.3, 0.378],
            [0.391, -0.148, 6.3, 0.378],
            #[0.391, -0.148, 6.3, 0.378]  # Extend the last value for extrapolation
        ]
        
    elif layout_angle == 45: # Squared Staggered tube arrangement: 
        a_values = [
            [1.55, -0.667, 1.93, 0.5],
            [1.55, -0.667, 1.93, 0.5],
            [0.498, -0.656, 1.93, 0.5],
            [0.73, -0.5, 1.93, 0.5],
            [0.37, -0.396, 1.93, 0.5],
            [0.37, -0.396, 1.93, 0.5],
        ]   
        b_values = [
            [32, -1, 6.59, 0.52],
            [32, -1, 6.59, 0.52],
            [26.2, -0.913, 6.59, 0.52],
            [3.5, -0.476, 6.59, 0.52],
            [0.333, -0.136, 6.59, 0.52],
            [0.303, -0.126, 6.59, 0.52],
        ]
    
    elif layout_angle == 60: # Triangular Staggered tube arrangement
        a_values = [
            [1.4, -0.667, 1.45, 0.519],
            [1.4, -0.667, 1.45, 0.519],
            [1.36, -0.657, 1.45, 0.519],
            [0.593, -0.477, 1.45, 0.519],
            [0.321, -0.388, 1.45, 0.519],
            [0.321, -0.388, 1.45, 0.519],
        ]
        b_values = [
            [48, -1, 7, 0.5],
            [48, -1, 7, 0.5],
            [45.1, -0.973, 7, 0.5],
            [4.57, -0.476, 7, 0.5],
            [0.486, -0.152, 7, 0.5],
            [0.372, -0.123, 7, 0.5],
        ]
    
    else:
        print("Tube arrangement shall be either square Inline (tube_layout = 0), square Staggered (tube_layout = 45) or triangular staggered (tube_layout = 60)")

    # Interpolate on the log scale
    log_Re = np.log10(Re)
    log_breakpoints = np.log10(breakpoints)
    
    a = np.zeros_like(a_values[0])
    b = np.zeros_like(b_values[0])
    
    for i in range(len(a)):
        a[i] = np.interp(log_Re, log_breakpoints, [a_values[j][i] for j in range(len(a_values))])
        b[i] = np.interp(log_Re, log_breakpoints, [b_values[j][i] for j in range(len(b_values))])
    
    return a, b

def shell_bell_delaware_htc(m_dot_shell, T_shell, T_shell_w, P_shell, shell_fluid, params):
    """
    Inputs
    ----------
    m_dot : shell fluid flowrate [kg/s]
    T_in : 
    T_w : 
    P_in : 
    fluid : 
    geom : 
    A_in_one_tube : 

    Outputs
    -------
    h : Shell-Side heat transfer coefficient [W/(m^2 * K)]
    DP_shell : Shell Pressure Drop [bar]
    
    References
    -------
    https://cheguide.com/heat_exchanger_rating.html
    """
    
    "1) HTC"
    
    "1.1) Ideal HTC"
    
    # Bulk fluid thermodynamical properties
    
    mu = PropsSI('V','T',T_shell, 'P',P_shell,shell_fluid)
    cp = PropsSI('C','T',T_shell, 'P',P_shell,shell_fluid)
    Pr = PropsSI('PRANDTL','T',T_shell, 'P',P_shell,shell_fluid)

    try:
        # Wall fluid thermodynamical properties
        mu_w = PropsSI('V','T',T_shell_w, 'P',P_shell,shell_fluid)
    except: 
        mu_w = PropsSI('V','T',T_shell, 'P',P_shell,shell_fluid)

    # Effective sections and hydraulic diameter
    S_T = s_max(params['Tube_OD'], params['pitch_ratio'], params['Shell_ID'], params['central_spacing'], params['tube_layout']) # m^2
    S_L = s_L(params['Baffle_cut'], params['Shell_ID'], params['n_tubes'], params['Tube_OD']) # m^2

    # Transversal to tube Flow speed and Reynolds number
    G_s = m_dot_shell/S_T
    
    Re = G_s*(params['Tube_OD']/mu)
        
    # Coefficients required for modelling tube bank flow
    a_vec, b_vec = bell_delaware_coefs(params['tube_layout'], Re)
    a = a_vec[2]/(1+0.14*Re**a_vec[3])
    j_i = a_vec[0]*(1.33/params['pitch_ratio'])**a * Re**a_vec[1]
        
    # Ideal Tube bank heat transfer coefficient
    h_id = j_i * cp * (m_dot_shell/S_T) * (Pr)**(-2/3) * (mu/mu_w)**0.14 
        
    "1.2) Correction factor for Baffle Window Flow : J_c"
    
    # Diameter effectively covered by the tube
    D_CTL = params['D_OTL'] - params['Tube_OD']
    
    theta_CTL = 2*np.arccos((params['Shell_ID'] * (1 - 2*params['Baffle_cut']/100))/D_CTL)
    F_w = (theta_CTL - np.sin(theta_CTL))/(2*np.pi)
    F_c = 1 - 2*F_w
    
    J_c = 0.55 + 0.72*F_c
    
    "1.3) Correction factor for Baffle Leakage : J_l"

    # Clearances (tube to baffle and shell to baffle)
    S_tb = params['clear_TB'] # (np.pi/4)*((geom.Tube_OD + geom.clear_TB)**2 - geom.Tube_OD**2)*geom.n_tubes*(1-F_w)
    S_sb = params['clear_BS']

    r_L = (S_sb + S_tb)/S_L
    r_S = S_sb /(S_sb + S_tb)

    J_l = 0.44*(1-r_S) + (1 - 0.44*(1-r_S))*np.e**(-2.2*r_L)
        
    "1.4) Correction factor for Bundle Bypass effects : J_b"
    
    # Tube pitch
    P_T = params['Tube_OD']*params['pitch_ratio']

    # Tube pitch relative to the flow    
    if params['tube_layout'] == 0: # Square inline tube arrangement
        P_p = P_T # for 90° layout
    elif params['tube_layout'] == 60: # Triangular staggered tube arrangement
        P_p = P_T * (3**(0.5)/2) 
    elif params['tube_layout'] == 45: # Square staggered tube arrangement
        P_p = P_T / (2**(0.5)) 
    
    # Number of effective tube rows crossed between two baffles
    N_tcc = (params['Shell_ID']/P_p)*(1 - 2*params['Baffle_cut']/100)
    
    # Number of sealing strips
    N_ss = params['N_strips']
    
    # Sealing strips ratio
    r_ss = N_ss/N_tcc

    if r_ss > 0.5:
        J_b = 1
    else:
        # Bypass flow area
        S_b = params['central_spacing'] * (params['Shell_ID'] - params['D_OTL'] - params['Tube_OD']/2)
        
        if Re < 100:
            C_j = 1.35
        else:
            C_j = 1.25
            
        J_b = np.e**(-C_j*(S_b / S_T)*(1 - (2*r_ss)**(1/3)))
    
    "1.5) Correction factor for adverse temperature gradient : J_r"
    
    # Effective number of baffles
    N_B = 1 + int((params['Tube_L'] - params['Tubesheet_t'] - params['inlet_spacing'] - params['outlet_spacing']) / params['central_spacing'])
    
    # Number of effective tube rows crossed in the main flow region ?
    N_tcw = (0.8/P_p)*(params['Shell_ID']*(params['Baffle_cut']/100) - ((params['Shell_ID']-(params['D_OTL']-params['Tube_OD']))/2))
            
    if Re >= 100:
        J_r = 1
    else:
        N_C = (N_tcw + N_tcc)*(1 + N_B)
        J_rl = (10/N_C)*0.18   

        if Re <= 20:
            J_r = J_rl
        else:
            J_r = J_rl + (20-Re)*(J_rl - 1)/80


    "1.6) Correction factor for unequal baffle spacing : J_s"
    
    if Re >= 100:
        n1 = 0.6
    else:
        n1 = 1/3
        
    J_s = ((N_B-1)+(params['inlet_spacing']/ params['central_spacing'])**(1-n1) + (params['outlet_spacing']/params['central_spacing'])**(1-n1))/((N_B-1)+(params['inlet_spacing']/params['central_spacing']) + (params['outlet_spacing']/params['central_spacing']))
    
    "1.7) Final Heat Transfer Coefficient"
    
    h = h_id*J_c*J_l*J_b*J_r*J_s
    
    return h
