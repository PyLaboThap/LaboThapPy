# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 17:34:23 2024

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

def htc_tube_and_fins_annular(fluid, params, P_in, T_in, m_dot_in):
    """
    Parameters
    ----------
    fluid : Fluid Name
    
    params : HTX parameters dictionnary
    
    P_in : Inlet fluid pressure [Pa]
    
    T_in : Inlet fluid temperature [K]
    
    m_dot_in : Inlet fluid flow rate [kg/s]

    Returns
    -------
    A_tot : Total outer area
        
    A_in : Tube internal area
    
    R_cond : Conduction resistance
    
    h_rdc : Finned side heat transfer coefficient
    
    Source
    -------
    HANDBOOK FOR TRANSVERSELY FINNED TUBE HEAT EXCHANGER DESIGN 
    EUGENE PIS’MENNYI, GEORGIY POLUPAN, IGNACIO CARVAJAL-MARISCAL, FLORENCIO SANCHEZ-SILVA, IGOR PIORO
    
    """
    
    # We consider here a staggered bank
    
    "Air data (Pure air considered)"
    
    rho_in = PropsSI('D','P',P_in,'T',T_in,fluid)
    V_dot_in = m_dot_in/rho_in # m^3/s
    
    "Geom data"

    n_tubes = params['n_tubes']
    Tube_ID = params['Tube_OD'] - 2*params['Tube_t']
    n_tpr = n_tubes/(params['n_rows']*params['n_passes'])
    
    HTX_L = params['L']
    HTX_W = params['w']
    
    Fin_L = (params['Fin_OD'] - params['Tube_OD'])/2
    
    "HT Area computations"
    
    # Fin HT area
    N_fins = params['Tube_L']*params['Fin_per_m'] - 1                                # Number of Fins
    Fin_spacing = (params['Tube_L'] - N_fins*params['Fin_t'])/(N_fins+1)

    A_out_fin = 2*np.pi*(params['Fin_OD']/2)*params['Fin_t']*N_fins*n_tubes*params['n_passes']
    A_out_tube = 2*np.pi*(params['Tube_OD']/2)*n_tubes*params['n_passes']*(params['Tube_L'] - N_fins*params['Fin_t'])
    A_out_plate_fin = 2*np.pi*((params['Fin_OD']/2)**2 - (params['Tube_OD']/2)**2)*N_fins*n_tubes*params['n_passes']
    
    A_out_tot = A_out_fin + A_out_tube + A_out_plate_fin
    
    # A_tot = params['A_finned']
    
    Af_A = (A_out_tot - A_out_tube) / A_out_tot # A_fin/A_tot
    At_A = A_out_tube / A_out_tot # A_tube/A_tot
 
    "Internal Tube Area"
    
    A_in = np.pi*Tube_ID*params['Tube_L']*n_tubes*params['n_passes']
    
    "R_cond"

    Tube_nom_OD = Tube_ID + (params['Fin_OD']*(N_fins*params['Fin_t']) + params['Tube_OD']*(params['Tube_L']-N_fins*params['Fin_t']))/params['Tube_L']

    fact_cond_1 = np.log(Tube_nom_OD/Tube_ID)
    fact_cond_2 = 2*np.pi*params['k_fin']*params['Tube_L']*n_tubes*params['n_passes']
    R_cond = fact_cond_1/fact_cond_2  

    "Gas velocity - Flow area"
    
    # Fin conventional length
    D_fin_c = params['Tube_OD'] + (2*Fin_L*params['Fin_t'])/Fin_spacing
    
    Tube_diag_pitch = np.sqrt(2)*params['pitch'] # Square staggered bank
    psi_c = (params['pitch'] - D_fin_c)/(Tube_diag_pitch - D_fin_c)

    if psi_c > 2:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*params['Tube_L'])*(2/psi_c)
    else:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*params['Tube_L'])

    u_in_air = V_dot_in/S_flow
    
    # print("HTX_W", HTX_W)
    # print("n_tpr*D_fin_c", n_tpr*D_fin_c)
    
    "h_c computation"
    
    k_g, Pr_g, nu_g = PropsSI(("L","PRANDTL","V"),'T',T_in,'P',P_in,fluid)
    
    "Fin Coefficient - Psi_f"
    
    Psi_f = 1/(2*params['Tube_OD']*Fin_spacing)*(params['Fin_OD']**2 - params['Tube_OD']**2 + 2*params['Fin_OD']*params['Fin_t']) + 1 - (params['Fin_t']/Fin_spacing)
                                           
    "Bundle shape param X"
    
    sigma_1 = params['pitch_ratio']
    sigma_2 = params['pitch_ratio']
    
    X = sigma_1/sigma_2 - 1.26/Psi_f - 2
        
    "n"
    n = 0.7 + 0.08*np.tanh(X) + 0.005*Psi_f
        
    "C_q"
    C_q = (1.36 - np.tanh(X))*(1.1/(Psi_f + 8) - 0.014)

    "C_z"
    if params['n_rows']*params['n_passes'] > 8:
        C_z = 1
    
    elif sigma_1/sigma_2  > 2:
        C_z = 3.5*(params['n_rows']*params['n_passes'])**0.03 - 2.72
    else:
        C_z = 3.15*(params['n_rows']*params['n_passes'])**0.05 - 2.5
    
    h_c = 1.13*C_z*C_q*(k_g/params['Tube_OD'])*(u_in_air*params['Tube_OD']/nu_g)**n * Pr_g**(0.33)    
        
    "Fin Efficiency"
    
    eta_f, psi_eta = eta_fin_straight_rect(params['Fin_t'], Fin_L, h_c, params['k_fin'], params['Fin_OD'], params['Tube_OD'])
    
    "Reduced HTC"
    
    h_rdc = (Af_A*eta_f*psi_eta + At_A)*h_c

    return h_rdc, A_out_tot

def htc_tube_and_fins_square(fluid, params, P_in, T_in, m_dot_in):
    """
    Parameters
    ----------
    fluid : Fluid Name
    
    params : HTX parameters dictionnary
    
    P_in : Inlet fluid pressure [Pa]
    
    T_in : Inlet fluid temperature [K]
    
    m_dot_in : Inlet fluid flow rate [kg/s]

    Returns
    -------
    A_tot : Total outer area
        
    A_in : Tube internal area
    
    R_cond : Conduction resistance
    
    h_rdc : Finned side heat transfer coefficient
    
    Source
    -------
    HANDBOOK FOR TRANSVERSELY FINNED TUBE HEAT EXCHANGER DESIGN 
    EUGENE PIS’MENNYI, GEORGIY POLUPAN, IGNACIO CARVAJAL-MARISCAL, FLORENCIO SANCHEZ-SILVA, IGOR PIORO
    
    """
    
    # We consider here a staggered bank
    
    "Gas data"
    
    rho_in = PropsSI('D','P',P_in,'T',T_in,fluid)
    V_dot_in = m_dot_in/rho_in # m^3/s
    
    "Geom data"

    n_tubes = params['n_tubes']
    Tube_ID = params['Tube_OD'] - 2*params['Tube_t']
    n_tpr = n_tubes/(params['n_rows'])
    
    # print("n_tpr", n_tpr)
    
    HTX_L = params['L']
    HTX_W = params['w']
    
    Fin_L = (params['Fin_OD'] - params['Tube_OD'])/2
    
    "HT Area computations"
    
    # Fin HT area
    N_fins = params['Tube_L']*params['Fin_per_m'] - 1                                # Number of Fins
    Fin_spacing = (params['Tube_L'] - N_fins*params['Fin_t'])/(N_fins+1)
   
    A_r = 2*(params['Fin_OD']**2 - 0.785*params['Tube_OD']**2 + 2*params['Fin_OD']*params['Fin_t'])*(params['Tube_L']/Fin_spacing)*params['n_tubes']*params['n_passes']  # 
    
    L_t = params['Tube_L'] - N_fins*params['Fin_t']
    A_t = np.pi*params['Tube_OD']*(params['Tube_L']*(1 - params['Fin_t']/Fin_spacing)*params['n_tubes']*params['n_passes'] + L_t)
    
    A_tot = A_r + A_t    
    
    Af_A = A_r/A_tot # A_fin/A_tot
    At_A = A_t / A_tot # A_tube/A_tot
 
    "Internal Tube Area"
    
    A_in = np.pi*Tube_ID*params['Tube_L']*n_tubes
    
    "R_cond"

    Tube_nom_OD = Tube_ID + (params['Fin_OD']*(N_fins*params['Fin_t']) + params['Tube_OD']*(params['Tube_L']-N_fins*params['Fin_t']))/params['Tube_L']

    fact_cond_1 = np.log(Tube_nom_OD/Tube_ID)
    fact_cond_2 = 2*np.pi*params['k_fin']*params['Tube_L']*n_tubes*params['n_passes']
    R_cond = fact_cond_1/fact_cond_2  
    
    "Gas velocity - Flow area"
    
    # Fin conventional length
    D_fin_c = params['Tube_OD'] + (2*Fin_L*params['Fin_t'])/Fin_spacing
    
    Tube_diag_pitch = np.sqrt(2)*params['pitch'] # Square staggered bank
    psi_c = (params['pitch'] - D_fin_c)/(Tube_diag_pitch - D_fin_c)

    if psi_c > 2:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*params['Tube_L'])*(2/psi_c)
    else:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*params['Tube_L'])

    u_in = V_dot_in/S_flow
    
    # print("HTX_W", HTX_W)
    # print("n_tpr*D_fin_c", n_tpr*D_fin_c)

    "h_c computation"
    
    k_g, Pr_g, nu_g = PropsSI(("L","PRANDTL","V"),'T',T_in,'P',P_in,fluid)
    
    "Fin Coefficient - Psi_f"
    
    # Psi_f = 1/(2*geom.Tube_OD*Fin_spacing)*(geom.Fin_OD**2 - geom.Tube_OD**2 + 2*geom.Fin_OD*geom.Fin_t) + 1 - (geom.Fin_t/Fin_spacing)

    Psi_f = 2*(params['Fin_OD']**2 - 0.785*params['Tube_OD']**2 + 2*params['Fin_OD']*params['Fin_t'])/(np.pi*params['Tube_OD']*Fin_spacing) + 1 - (params['Fin_t']/Fin_spacing)

    "Bundle shape param X"
    
    sigma_1 = params['pitch_ratio']
    sigma_2 = params['pitch_ratio']
    
    X = sigma_1/sigma_2 - 1.26/Psi_f - 2
        
    "n"
    n = 0.7 + 0.08*np.tanh(X) + 0.005*Psi_f
        
    "C_q"
    C_q = (1.36 - np.tanh(X))*(1.1/(Psi_f + 8) - 0.014)

    "C_z"
    
    if params['n_rows']*params['n_passes'] > 8:
        C_z = 1
    elif sigma_1/sigma_2  > 2:
        C_z = 3.5*(params['n_rows']*params['n_passes'])**0.03 - 2.72
    else:
        C_z = 3.15*(params['n_rows']*params['n_passes'])**0.05 - 2.5
    
    h_c = 1.13*C_z*C_q*(k_g/params['Tube_OD'])*(u_in*params['Tube_OD']/nu_g)**n * Pr_g**(0.33)    

    "Fin Efficiency"
    
    eta_f, psi_eta = eta_fin_straight_rect(params['Fin_t'], Fin_L, h_c, params['k_fin'], params['Fin_OD'], params['Tube_OD'])
    
    "Reduced HTC"
    
    h_rdc = (Af_A*eta_f*psi_eta + At_A)*h_c
    
    return h_rdc, A_tot

def htc_tube_and_fins(fluid, params, P_in, T_in, m_dot_in, Fin_type):
    """
    Parameters
    ----------
    fluid : Fluid Name
    
    params : HTX parameters dictionnary
    
    P_in : Inlet fluid pressure [Pa]
    
    T_in : Inlet fluid temperature [K]
    
    m_dot_in : Inlet fluid flow rate [kg/s]

    Fin_type : Annular or Square [-]

    Returns
    -------
    A_tot : Total outer area
        
    A_in : Tube internal area
    
    R_cond : Conduction resistance
    
    h_rdc : Finned side heat transfer coefficient
    
    Source
    -------
    HANDBOOK FOR TRANSVERSELY FINNED TUBE HEAT EXCHANGER DESIGN 
    EUGENE PIS’MENNYI, GEORGIY POLUPAN, IGNACIO CARVAJAL-MARISCAL, FLORENCIO SANCHEZ-SILVA, IGOR PIORO
    
    """
    
    if Fin_type == "Annular":
        h_rdc, A_tot = htc_tube_and_fins_annular(fluid, params, P_in, T_in, m_dot_in)
        
    elif Fin_type == "Square":
        h_rdc, A_tot = htc_tube_and_fins_square(fluid, params, P_in, T_in, m_dot_in)
    
    return h_rdc, A_tot

def eta_fin_straight_rect(t_fin, L_fin, h_c, k_fin, D_fin, D_tube):
    """
    This function determines the efficiency of a straight fin with a rectangular profile given the dimensions, convection coefficient, and fin conductivity
    +
    This function returns the efficiency of a straight fin with a rectangular profile given the dimensionless parameter mL
    
    Inputs:
    h - heat transfer coefficient
    k - fin conductivity
    L - fin length or depth
    t - base thickness of fin
    
    Output
    e_f - Fin efficiency
    """

    if(h_c < 0):
        print('The heat transfer coefficient given is less than zero. The value for h is', h_c)
        return 
	
    # Computation of mL dimensionless parameter
    m = np.sqrt(2*h_c/(k_fin*t_fin)) 
    L_fin_prime = L_fin*(1+(0.191+0.054*(D_fin/D_tube))*np.log(D_fin/D_tube))
    mL=m*L_fin_prime

    if mL == 0:
        eta_f = 1.0
    else:
        eta_f = np.tanh(mL)/mL		
        
    # Corr factor to fin efficiency
    psi_eta = 1 - 0.016*(D_fin/D_tube - 1)*(1+ np.tanh(2*mL - 1))
    
    return eta_f, psi_eta

def DP_tube_and_fins(fluid, geom, P_in, T_in, m_dot_in):
    """

    Parameters
    ----------
    fluid : Fluid Name
    
    geom : ACC geometry class
    
    P_in : Inlet fluid pressure [Pa]
    
    T_in : Inlet fluid temperature [K]
    
    m_dot_in : Inlet fluid flow rate [kg/s]

    Returns
    -------
    DP : Pressure drop across the tube bank [Pa]

    Source
    -------
    HANDBOOK FOR TRANSVERSELY FINNED TUBE HEAT EXCHANGER DESIGN 
    EUGENE PIS’MENNYI, GEORGIY POLUPAN, IGNACIO CARVAJAL-MARISCAL, FLORENCIO SANCHEZ-SILVA, IGOR PIORO
    """

    # We consider here a staggered bank
    
    "Air data (Pure air considered)"
    
    rho_g, nu_g = PropsSI(('D','V'),'P',P_in,'T',T_in,fluid)
    V_dot_in = m_dot_in/rho_g # m^3/s
    
    "Geom data"

    n_tubes = geom.n_tubes*geom.n_passes
    Tube_ID = geom.Tube_OD - 2*geom.Tube_t
    n_tpr = n_tubes/geom.n_rows
    
    HTX_L = geom.L
    HTX_W = geom.w
    
    Fin_L = (geom.Fin_OD - geom.Tube_OD)/2
    
    "Gas velocity - Flow area"

    # Fin HT area
    N_fins = geom.Tube_L*geom.Fin_per_m - 1 # Number of Fins
    Fin_spacing = (geom.Tube_L - N_fins*geom.Fin_t)/(N_fins+1)
    
    # Fin conventional length
    D_fin_c = geom.Tube_OD + (2*Fin_L*geom.Fin_t)/Fin_spacing
    
    Tube_diag_pitch = np.sqrt(2)*geom.pitch # Square staggered bank
    psi_c = (geom.pitch - D_fin_c)/(Tube_diag_pitch - D_fin_c)

    if psi_c > 2:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*geom.Tube_L)*(2/psi_c)
    else:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*geom.Tube_L)

    u_in_air = V_dot_in/S_flow

    "Resistance coefficient for a tube row - xi,o"
    
    # n and C_r
    
    S_1 = geom.pitch
    S_2 = geom.pitch
    
    fact_1 = np.pi*(geom.Tube_OD*Fin_spacing + 2*Fin_L*geom.Fin_t + 2*Fin_L*(Fin_L + geom.Tube_OD))
    fact_2 = geom.pitch*Fin_spacing - (geom.Tube_OD*Fin_spacing + 2*Fin_L*geom.Fin_t)
    
    A_totF = fact_1/fact_2 # Reduced length of developped surface
        
    n = 0.17*(A_totF)**0.25 * (S_1/S_2)**0.57 * np.exp(-0.36*(S_1/S_2))
    C_r = 2.8*(A_totF)**0.53 * (S_1/S_2)**1.3 * np.exp(-0.9*(S_1/S_2))

    # C_z
    
    if geom.n_rows >= 6:
        C_z = 1
    else:
        C_z = np.exp(0.1*(6/geom.n_rows - 1))
        
    # D_eq
    
    D_eq = 2*(Fin_spacing*(geom.pitch-geom.Tube_OD) - 2*Fin_L*geom.Fin_t)/(2*Fin_L + Fin_spacing)

    # xi_o
    
    xi_o = C_z*C_r*(u_in_air*D_eq/nu_g)**(-n)
    
    "Pressure drop computation"
    
    C_op = 1.1 # Value from the book : correction factor taking account of real operating conditions of the heat-transfer surface
    DP = C_op*xi_o*geom.n_rows*(rho_g*u_in_air**2)/2
    
    return DP
