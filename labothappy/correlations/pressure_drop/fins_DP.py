# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 17:34:23 2024

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI


def DP_tube_and_fins(fluid, params, P_in, T_in, m_dot_in):
    
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

    n_tubes = params['n_tubes']
    Tube_ID = params['Tube_OD'] - 2*params['Tube_t']
    n_tpr = n_tubes/params['n_rows']
    
    HTX_L = params['Tube_L']
    HTX_W = params['w']
    
    Fin_L = (params['Fin_OD'] - params['Tube_OD'])/2
    
    "Gas velocity - Flow area"

    # Fin HT area
    N_fins = params['Tube_L']*params['Fin_per_m'] - 1 # Number of Fins
    Fin_spacing = (params['Tube_L'] - N_fins*params['Fin_t'])/(N_fins+1)
    
    # Fin conventional length
    D_fin_c = params['Tube_OD'] + (2*Fin_L*params['Fin_t'])/Fin_spacing
    
    Tube_diag_pitch = np.sqrt(2)*np.sqrt(params['pitch_V']*params['pitch_H']) # Square staggered bank
    psi_c = (np.sqrt(params['pitch_V']*params['pitch_H']) - D_fin_c)/(Tube_diag_pitch - D_fin_c)

    if psi_c > 2:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*params['Tube_L'])*(2/psi_c)
    else:
        S_flow = (HTX_L*HTX_W - n_tpr*D_fin_c*params['Tube_L'])

    u_in_air = V_dot_in/S_flow

    "Resistance coefficient for a tube row - xi,o"
    
    # n and C_r
    
    S_1 = params['pitch_V']
    S_2 = params['pitch_H']
    
    fact_1 = np.pi*(params['Tube_OD']*Fin_spacing + 2*Fin_L*params['Fin_t'] + 2*Fin_L*(Fin_L + params['Tube_OD']))
    fact_2 = np.sqrt(params['pitch_V']*params['pitch_H'])*Fin_spacing - (params['Tube_OD']*Fin_spacing + 2*Fin_L*params['Fin_t'])
    
    A_totF = fact_1/fact_2 # Reduced length of developped surface
        
    n = 0.17*(A_totF)**0.25 * (S_1/S_2)**0.57 * np.exp(-0.36*(S_1/S_2))
    C_r = 2.8*(A_totF)**0.53 * (S_1/S_2)**1.3 * np.exp(-0.9*(S_1/S_2))

    # C_z
    
    if params['n_rows'] >= 6:
        C_z = 1
    else:
        C_z = np.exp(0.1*(6/params['n_rows'] - 1))
        
    # D_eq
    
    D_eq = 2*(Fin_spacing*(np.sqrt(params['pitch_V']*params['pitch_H'])-params['Tube_OD']) - 2*Fin_L*params['Fin_t'])/(2*Fin_L + Fin_spacing)

    # xi_o
    
    xi_o = C_z*C_r*(u_in_air*D_eq/nu_g)**(-n)
    
    "Pressure drop computation"
    
    C_op = 1.1 # Value from the book : correction factor taking account of real operating conditions of the heat-transfer surface
    DP = C_op*xi_o*params['n_rows']*(rho_g*u_in_air**2)/2
    
    return DP
