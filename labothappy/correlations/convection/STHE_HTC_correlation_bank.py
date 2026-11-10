# -*- coding: utf-8 -*-
"""
Shell-side heat transfer coefficient of shell & tube heat exchangers, in the
(AS, geom, cond) formalism of the BPHX correlations.

cond : dictionary of the cell conditions built by the model (_cond_dict). The shell correlations read
       m_dot (flow of the whole side), mu, mu_w (optional: without it mu_w = mu), Pr, k, cp.
geom : geometry of the side: all the heat exchanger parameters (Tube_OD, pitch_ratio, Shell_ID,
       central_spacing, tube_layout, ...). n_parallel (optional, default 1) divides the flow rate.
"""

import numpy as np

#%% GEOMETRY HELPERS

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


def s_max(Tube_OD, pitch_ratio, Shell_ID, D_OTL, central_spacing, tube_layout): # Maximum flow section (m**2)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    Shell_ID        : Shell inside diameter [m]
    D_OTL           : Outer tube limit diameter [m]
    central_spacing : Spacing between two baffles [m]
    tube_layout     : Angle of tube layout [°]: 0 square inline, 45 square staggered, 60 triangular staggered
    
    Output
    ------
    
    s_max : Maximum flow area in the tube bank [m^2]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio
    
    # TODO: the three effective pitch factors below are kept from the previous version and are not
    #       verified; the labels of 45 and 60 were swapped in the previous version of this function.
    if tube_layout == 0: # Square inline tube arrangement
    
        p_T_eff = p_T
    
    elif tube_layout == 45: # Square staggered tube arrangement
    
        p_T_eff = (2**(1/2))*p_T
    
    elif tube_layout == 60: # Triangular staggered tube arrangement
    
        p_T_eff = (2/3**(1/2))*p_T
    
    else:
        raise ValueError("Tube arrangement shall be either square inline (tube_layout = 0), "
                         "square staggered (tube_layout = 45) or triangular staggered (tube_layout = 60)")

    # (the bundle-to-shell gap term is not divided by p_T_eff: both terms are lengths)
    s_max = central_spacing * ((Shell_ID - D_OTL) + (D_OTL - Tube_OD)*(p_T - Tube_OD)/p_T_eff)
    
    return s_max


def s_L(Baffle_cut, Shell_ID, n_tubes, Tube_OD):
    """
    Input
    -----
    
    Baffle_cut      : Baffle cut [% of the shell diameter]
    Shell_ID        : Shell inside diameter [m]
    n_tubes         : Number of tubes [-]
    Tube_OD         : Tube outside diameter [m]
    
    Output
    ------
    
    s_L : Flow area in the longitudinal shell direction (baffle window) [m^2]:
          the fraction Baffle_cut of the shell section, minus the tubes in it (uniform tube density)
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    return (Baffle_cut/100)*(np.pi/4)*(Shell_ID**2 - n_tubes*Tube_OD**2)


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

#%% KERN

def STHE_kern_shell_HTC(AS, geom, cond):
    """
    Kern shell-side heat transfer coefficient [W/(m^2 K)].

    cond : m_dot (flow of the whole side) [kg/s], mu, mu_w (optional) [Pa s], Pr, k [W/(m K)]
    geom : Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout; n_parallel (optional)

    Reference
    ---------
    Process Heat Transfer - D. Q. Kern
    """
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)        # flow rate per shell
    mu, Pr, k = cond['mu'], cond['Pr'], cond['k']
    mu_w = cond.get('mu_w', mu)

    S_T = s_max_kern(geom['Tube_OD'], geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])          # m^2
    D_hydro = d_h(geom['Tube_OD'], geom['pitch_ratio'], geom['tube_layout'])

    # Reynolds number on the equivalent diameter and the cross-flow mass velocity
    Re = (m_dot / S_T) * D_hydro / mu

    # McAdams (with baffles)
    Nu = 0.36 * Pr**(1/3) * Re**0.55 * (mu / mu_w)**0.14

    return Nu * k / D_hydro

#%% DONOHUE

def STHE_donohue_shell_HTC(AS, geom, cond):
    """
    Donohue shell-side heat transfer coefficient [W/(m^2 K)]. No viscosity correction.

    cond : m_dot (flow of the whole side) [kg/s], mu [Pa s], Pr, k [W/(m K)]
    geom : Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout, Baffle_cut [%], n_tubes;
           n_parallel (optional)

    Reference
    ---------
    Process Heat Transfer - D. Q. Kern
    """
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)        # flow rate per shell
    mu, Pr, k = cond['mu'], cond['Pr'], cond['k']
    d_o = geom['Tube_OD']

    # Cross-flow section (Kern) and longitudinal section (baffle window)
    S_T = s_max_kern(d_o, geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])          # m^2
    S_L = s_L(geom['Baffle_cut'], geom['Shell_ID'], geom['n_tubes'], d_o)   # m^2

    # Re = rho V_R d_o / mu with V_R = sqrt(V_L V_T) and V = m_dot / (rho S): the density cancels
    Re = m_dot * d_o / (mu * np.sqrt(S_T * S_L))

    Nu = 0.22 * Re**0.6 * Pr**0.33

    return Nu * k / d_o

#%% BELL-DELAWARE

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
        ]
        
        b_values = [
            [35, -1, 6.3, 0.378],
            [35, -1, 6.3, 0.378],
            [32.1, -0.963, 6.3, 0.378],
            [6.09, -0.602, 6.3, 0.378],
            [0.0815, 0.022, 6.3, 0.378],
            [0.391, -0.148, 6.3, 0.378],
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
        raise ValueError("Tube arrangement shall be either square inline (tube_layout = 0), "
                         "square staggered (tube_layout = 45) or triangular staggered (tube_layout = 60)")

    # Interpolate on the log scale
    log_Re = np.log10(Re)
    log_breakpoints = np.log10(breakpoints)
    
    a = np.zeros_like(a_values[0])
    b = np.zeros_like(b_values[0])
    
    for i in range(len(a)):
        a[i] = np.interp(log_Re, log_breakpoints, [a_values[j][i] for j in range(len(a_values))])
        b[i] = np.interp(log_Re, log_breakpoints, [b_values[j][i] for j in range(len(b_values))])
    
    return a, b


def STHE_bell_delaware_shell_HTC(AS, geom, cond):
    """
    Bell-Delaware shell-side heat transfer coefficient [W/(m^2 K)].

    cond : m_dot (flow of the whole side) [kg/s], mu, mu_w (optional) [Pa s], cp [J/(kg K)], Pr
    geom : Tube_OD, pitch_ratio, Shell_ID, D_OTL, central_spacing, tube_layout, Baffle_cut [%],
           clear_TB, clear_BS, n_tubes, N_strips, Tube_L, Tubesheet_t, inlet_spacing, outlet_spacing;
           n_parallel (optional)

    References
    -------
    https://cheguide.com/heat_exchanger_rating.html

    + 

    https://fr.scribd.com/document/351507985/Shell-and-Tube-HX

    tube_layout : 
        30° or 60° (depends on definition) - Clean fluid (laminar + turbulent) (best for compactness)
        45° - Fouled + Laminar
        0 or 90° - Fouled + Turbulent (best for maintenance)
    """
    m_dot_shell = cond['m_dot'] / geom.get('n_parallel', 1)  # flow rate per shell
    mu, cp, Pr = cond['mu'], cond['cp'], cond['Pr']
    mu_w = cond.get('mu_w', mu)
    p = geom                                                  # geometry (all the heat exchanger parameters)

    "1) HTC"
    
    "1.1) Ideal HTC"

    # Effective sections and hydraulic diameter
    S_m = s_max(p['Tube_OD'], p['pitch_ratio'], p['Shell_ID'], p['D_OTL'], p['central_spacing'], p['tube_layout']) # m^2
    
    # Transversal to tube Flow speed and Reynolds number
    G_s = m_dot_shell/S_m
    
    Re = G_s*(p['Tube_OD']/mu)
        
    # Coefficients required for modelling tube bank flow
    a_vec, b_vec = bell_delaware_coefs(p['tube_layout'], Re)
    a = a_vec[2]/(1+0.14*Re**a_vec[3])
    j_i = a_vec[0]*(1.33/p['pitch_ratio'])**a * Re**a_vec[1]
        
    # Ideal Tube bank heat transfer coefficient
    h_id = j_i * cp * (m_dot_shell/S_m) * (Pr)**(-2/3) * (mu/mu_w)**0.14 
        
    "1.2) Correction factor for Baffle Window Flow : J_c"
    
    # Diameter effectively covered by the tube
    D_CTL = p['D_OTL'] - p['Tube_OD']
    
    theta_CTL = 2*np.arccos((p['Shell_ID'] * (1 - 2*p['Baffle_cut']/100))/D_CTL)
    F_w = (theta_CTL - np.sin(theta_CTL))/(2*np.pi)
    F_c = 1 - 2*F_w
    
    J_c = 0.55 + 0.72*F_c
    
    "1.3) Correction factor for Baffle Leakage : J_l"
    theta_ds = 2*np.arccos(1-2*(p['Baffle_cut']/100))
    
    # Clearances (tube to baffle and shell to baffle)
    S_tb = (np.pi/4)*((p['Tube_OD'] + p['clear_TB'])**2 - p['Tube_OD']**2)*p['n_tubes']*(1-F_w) # clear_TB - max - 0.794 mm in TEMA stdrads and between 0.397 mm 
    S_sb = 0.00436 * p['Shell_ID'] * p['clear_BS'] * (180/np.pi) * (2*np.pi - theta_ds) # - L_sb = 2 mm if D_s < 0.4m and = 1.6 + 0.004*D_s
    
    r_L = (S_sb + S_tb)/S_m
    r_S = S_sb /(S_sb + S_tb)

    J_l = 0.44*(1-r_S) + (1 - 0.44*(1-r_S))*np.e**(-2.2*r_L)
        
    "1.4) Correction factor for Bundle Bypass effects : J_b"
    
    # Tube pitch
    P_T = p['Tube_OD']*p['pitch_ratio']

    # Tube pitch relative to the flow    
    if p['tube_layout'] == 0: # Square inline tube arrangement
        P_p = P_T # for 90° layout
    elif p['tube_layout'] == 60: # Triangular staggered tube arrangement
        P_p = P_T * (3**(0.5)/2) 
    elif p['tube_layout'] == 45: # Square staggered tube arrangement
        P_p = P_T / (2**(0.5)) 
    else:
        raise ValueError(f"Unknown tube_layout {p['tube_layout']!r}: 0, 45 or 60.")
    
    # Number of effective tube rows crossed between two baffles
    N_tcc = (p['Shell_ID']/P_p)*(1 - 2*p['Baffle_cut']/100)
    
    # Number of sealing strips
    N_ss = p['N_strips']
    
    # Sealing strips ratio
    r_ss = N_ss/N_tcc

    if r_ss > 0.5:
        J_b = 1
    else:
        # Bypass flow area
        S_b = p['central_spacing'] * (p['Shell_ID'] - p['D_OTL'] - p['Tube_OD']/2) # TODO: check the - Tube_OD/2
        
        if Re < 100:
            C_j = 1.35
        else:
            C_j = 1.25
            
        J_b = np.e**(-C_j*(S_b / S_m)*(1 - (2*r_ss)**(1/3)))
    
    "1.5) Correction factor for adverse temperature gradient : J_r"
    
    # Effective number of baffles
    N_B = 1 + int((p['Tube_L'] - p['Tubesheet_t'] - p['inlet_spacing'] - p['outlet_spacing']) / p['central_spacing'])
    
    # Number of effective tube rows crossed in the main flow region ?
    N_tcw = (0.8/P_p)*(p['Shell_ID']*(p['Baffle_cut']/100) - ((p['Shell_ID']-(p['D_OTL']-p['Tube_OD']))/2))
            
    if Re >= 100:
        J_r = 1
    else:
        N_C = (N_tcw + N_tcc)*(1 + N_B)
        J_rl = (10/N_C)**0.18   

        if Re <= 20:
            J_r = J_rl
        else:
            J_r = J_rl + (20-Re)*(J_rl - 1)/80


    "1.6) Correction factor for unequal baffle spacing : J_s"
    
    if Re >= 100:
        n1 = 0.6
    else:
        n1 = 1/3
        
    J_s = ((N_B-1)+(p['inlet_spacing']/ p['central_spacing'])**(1-n1) + (p['outlet_spacing']/p['central_spacing'])**(1-n1))/((N_B-1)+(p['inlet_spacing']/p['central_spacing']) + (p['outlet_spacing']/p['central_spacing']))
    
    "1.7) Final Heat Transfer Coefficient"
    
    h = h_id*J_c*J_l*J_b*J_r*J_s

    return h