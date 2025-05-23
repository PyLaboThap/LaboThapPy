# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:36:53 2024

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

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
        
        D_h = (4*(0.43*p_T**2 - (0.5*np.pi*Tube_OD**2 /4)))/(0.5*np.pi*Tube_OD)
        # D_h = 2 * (3**(1/2) / np.pi) * (p_T**2 / Tube_OD) - Tube_OD
        
    return D_h

#%%

def shell_DP_kern(m_dot, T_wall, h_in, P_in, fluid, params):
    """
    Inputs
    ----------
    
    m_dot  : Flow rate [kg/s]
    T_wall : Wall Temperature [K]
    h_in   : Inlet Enthalpy [K]
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
    
    (rho, mu) = PropsSI(('D','V'),'H',h_in, 'P',P_in,fluid)
    
    mu_w = PropsSI('V','T',T_wall,'P',P_in,fluid)

    S_T = s_max(params['Tube_OD'], params['pitch_ratio'], params['Shell_ID'], params['central_spacing'], params['tube_layout']) # m^2
    D_hydro = d_h(params['Tube_OD'], params['pitch_ratio'], params['tube_layout'])

    V_t = m_dot/(S_T*rho)
    
    Re = rho*V_t*(D_hydro/mu)
    # print(Re)
    
    "2) DP"
    
    f = np.e**(0.576 - 0.19*np.log(Re)) # [-] : Friction coefficient
    G = m_dot/S_T # kg/(m^2 * s)
    phi_s = (mu/mu_w)**(0.14)  

    DP = (f*G**2 * params['Shell_ID'] * (params['cross_passes'] + 1))/(2*rho*D_hydro*phi_s)
    
    return DP

# def shell_DP_kern(m_dot, T_wall, h_in, P_in, fluid, params):
#     """
#     Inputs
#     ----------
    
#     m_dot  : Flow rate [kg/s]
#     T_wall : Wall Temperature [K]
#     h_in   : Inlet Enthalpy [K]
#     P_in   : Inlet pressure [Pa]
#     fluid  : fluid name [-]
#     params : HTX parameters [-]

#     Outputs
#     -------
    
#     h = shell side heat transfer coefficient [W/(m^2 * K)]    
#     DP = shell side pressure drop [Pa]
    
#     References
#     -------
#     Process Heat Transfer - D. Q. Kern
#     """
    
#     "1) HTC"
    
#     (rho, mu) = PropsSI(('D','V'),'H',h_in, 'P',P_in,fluid)
    
#     mu_w = PropsSI('V','T',T_wall,'P',P_in,fluid)

#     S_T = s_max_kern(params['Tube_OD'], params['pitch_ratio'], params['Shell_ID'], params['central_spacing'], params['tube_layout']) # m^2
#     D_hydro = d_h(params['Tube_OD'], params['pitch_ratio'], params['tube_layout'])
    
#     V_s = m_dot/(S_T*rho)
    
#     Re = rho*V_s*(D_hydro/mu)
    
#     "2) DP"
#     bo = 0.72
    
#     f = 2*bo*Re**(-0.15) # [-] : Friction coefficient

#     DP = f*(rho*V_s**2 / 2)*(params['Tube_L']/params['central_spacing'])*(params['Shell_ID']/D_hydro) # (f*G**2 * params['Shell_ID'] * (params['cross_passes'] + 1))/(2*rho*D_hydro*phi_s)
    
#     return DP

#%%

def shell_DP_donohue(m_dot, T_in, P_in, fluid, params):
    """
    Inputs
    ----------
    
    m_dot  : Flow rate [kg/s]
    T_in   : Inlet Temperature [K]
    P_in   : Inlet pressure [Pa]
    fluid  : fluid name [-]
    params : HTX geometrical parameters [-]

    Outputs
    -------
    
    DP = shell side pressure drop [Pa]
    
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
    
    "2) DP"
    
    f = np.e**(0.576 - 0.19*np.log(Re)) # [-] : Friction coefficient
    G = m_dot/S_T # kg/(m^2 * s)
    phi_s = 1 # (mu/mu_w)**(0.14)  

    DP = (f*G**2 * params['Shell_ID'] * (params['cross_passes'] + 1))/(2*rho*D_hydro*phi_s)
    
    return DP


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

#%% 

def shell_bell_delaware_DP(m_dot_shell, h_shell, P_shell, shell_fluid, params):
    """
    Inputs
    ----------
    m_dot : shell fluid flowrate
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
    
    "1) DP"

    "1.1) DP for ideal Tube bank"
    
    # Bulk fluid thermodynamical properties
    
    rho = PropsSI('D','H',h_shell, 'P',P_shell,shell_fluid)
    mu = PropsSI('V','H',h_shell, 'P',P_shell,shell_fluid)
    cp = PropsSI('C','H',h_shell, 'P',P_shell,shell_fluid)

    # Wall fluid thermodynamical properties
    mu_w = PropsSI('V','H',h_shell, 'P',P_shell,shell_fluid)
    
    # Effective sections and hydraulic diameter  
    S_T = s_max(params['Tube_OD'], params['pitch_ratio'], params['Shell_ID'], params['central_spacing'], params['tube_layout']) # m^2
    S_L = s_L(params['Baffle_cut'], params['Shell_ID'], params['n_tubes'], params['Tube_OD']) # m^2
    D_hydro = d_h(params['Tube_OD'], params['pitch_ratio'], params['tube_layout']) # m
    
    # Transversal to tube Flow speed and Reynolds number
    V_t = m_dot_shell/(S_T*rho)
    Re = rho*V_t*(D_hydro/mu)
    
    a_vec, b_vec = bell_delaware_coefs(params["tube_layout"], Re)
    
    # Coefficients required for modelling tube bank flow
    b = b_vec[2]/(1+0.14*Re**b_vec[3])
    f_i = b_vec[0]*(1.33/params["pitch_ratio"])**b * Re**b_vec[1]
    
    # Shell transverse mass Flux
    G_shell = m_dot_shell / (S_T)
    
    # Tube pitch
    P_T = params["Tube_OD"]*params["pitch_ratio"]
    
    # Tube pitch relative to the flow    
    if params["tube_layout"] == 0: # Square inline tube arrangement
        P_p = P_T # for 90° layout
    elif params["tube_layout"] == 60: # Triangular staggered tube arrangement
        P_p = P_T * (3**(0.5)/2) 
    elif params["tube_layout"] == 45: # Square staggered tube arrangement
        P_p = P_T / (2**(0.5))     

    # Number of effective tube rows crossed between two baffles
    N_tcc = (params["Shell_ID"]/P_p)*(1 - 2*params["Baffle_cut"]/100)
    
    # Ideal Pressure drop    
    DP_ideal = 2*f_i*(G_shell**2 / rho) * (mu/mu_w)**(0.14) * N_tcc
    
    "1.2) Correction factor for Baffle Leakage : R_l"

    # Diameter effectively covered by the tube
    D_CTL = params["D_OTL"] - params["Tube_OD"]
    
    theta_CTL = 2*np.arccos((params["Shell_ID"] * (1 - 2*params["Baffle_cut"]/100))/D_CTL)
    F_w = (theta_CTL - np.sin(theta_CTL))/(2*np.pi)

    # Clearances (tube to baffle and shell to baffle)
    S_tb = params["clear_TB"] # (np.pi/4)*((geom.Tube_OD + geom.clear_TB)**2 - geom.Tube_OD**2)*geom.n_tubes*(1-F_w)
    S_sb = params["clear_BS"]

    r_L = (S_sb + S_tb)/S_L
    r_S = S_sb /(S_sb + S_tb)

    p = 0.8 - 0.15*(1 + r_S)    
    R_l = np.e**(-1.33*(1+r_S)*r_L**p)
    
    "1.3) Pressure drop for window section : DP_w"
    
    theta_DS = 2*np.arccos(1 - 2*params["Baffle_cut"]/100)
    
    SWG = (params["Shell_ID"]**2/8)*(theta_DS - np.sin(theta_DS))
    SWT = params["n_tubes"]*F_w*(np.pi*params["Tube_OD"]**2/4)
    SW = SWG - SWT
    GW = m_dot_shell/(S_T*SW)**0.5
    DW = 4*SW /(np.pi*params["Tube_OD"]*params["n_tubes"]*F_w + theta_DS*params["Shell_ID"])
    
    # Number of effective tube rows crossed in the main flow region ?
    N_tcw = (0.8/P_p)*(params["Shell_ID"]*(params["Baffle_cut"]/100) - ((params["Shell_ID"]-(params["D_OTL"]-params["Tube_OD"]))/2))
    
    if Re >= 100: # turbulent flow
        DP_w = params["cross_passes"]*R_l*(2+0.6*N_tcw)*GW**2/(2*rho)
        
    else: # laminar flow
        a_coef = (26*GW*cp)/rho
        b_coef = N_tcw/(P_T-params["Tube_OD"])
        c_coef = params["central_spacing"]/DW**2
        d_coef = GW**2/rho
        
        DP_w = params["cross_passes"]*R_l*(a_coef*(b_coef + c_coef) + d_coef)
        
    "1.4) Correction factor for Bundle Bypass effect : R_b"

    # Number of sealing strips
    N_ss = params["N_strips"]
    
    # Sealing strips ratio
    r_ss = N_ss/N_tcc

    if r_ss < 0.5:
        if Re >= 100:
            C_r = 3.7
        else:
            C_r = 4.5
        
        # Bypass flow area
        S_b = params["central_spacing"] * (params["Shell_ID"] - params["D_OTL"] - params["Tube_OD"]/2)
        
        R_b = np.e**(-C_r*(S_b/S_T)*(1-(2*r_ss)**(1/3)))
        
    else:
        R_b = 1
        
    "1.5) Correction factor for unequal baffle spacing : R_s"

    if Re >= 100:
        n = 0.2
    else:
        n = 1
    
    R_s = 0.5*((params["central_spacing"]/params["inlet_spacing"])**(2-n) + (params["central_spacing"]/params["outlet_spacing"])**(2-n))
    
    "1.6) Final DP"
    
    # Pressure drop in Central Baffle spaces
    DP_PC = (params["cross_passes"] - 1)*DP_ideal*R_l*R_b
    
    # Pressure drop in entrance & exit baffle spaces
    DP_PE = DP_ideal*(1 + N_tcw/N_tcc)*R_b*R_s

    # Final DP
    DP_shell = DP_PE + DP_PC + DP_w

    return DP_shell
