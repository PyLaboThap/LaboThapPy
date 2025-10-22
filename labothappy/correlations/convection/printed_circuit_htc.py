# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 17:13:48 2025

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP

# def PCHE_Lee(alpha, D_c, k, L_c, mdot, mu, N_c, N_p, Pr, rho, R_p, side, t_2, t_3):
#     """
#     Reference
#     ---------
#     1) Design and Dynamic Modeling of Printed Circuit Heat Exchangers for Supercritical Carbon Dioxide Brayton Power Cycles (2018) 
    
#     Yuan Jiang, Eric Liese, Stephen E. Zitney, Debangsu Bhattacharyya
    
#     2) Thermal-hydraulic performance analysis of zigzag channel PCHE according to bending angle (2024)
    
#     Yoomyeong Lee, Seongmin Lee, Hong Beom Park, Kyoung Woo Seo, Donghwi Lee
#     """

#     # Flow speed
#     if side == "Cold":
#         v = mdot/(rho*(1/(1+R_p))*N_c*N_p*(np.pi*D_c**2)/8)
#     else:
#         v = mdot/(rho*(R_p/(1+R_p))*N_c*N_p*(np.pi*D_c**2)/8)
    
#     D_h = np.pi*D_c/(2+np.pi)
#     Re = (rho*v*D_h)/mu # Re in 0.5 to 2000
    
#     if alpha == 10:
#         p = 12.93*1e-3 # m
#         h = 1.14*1e-3 # m
#     elif alpha == 20:
#         p = 12.09*1e-3 # m
#         h = 2.2*1e-3 # m
#     elif alpha == 30:
#         p = 10.95*1e-3 # m
#         h = 3.16*1e-3 # m
#     elif alpha == 40:
#         p = 9.53*1e-3 # m
#         h = 4*1e-3 # m
    
#     Nu = 0.0479*Re*0.7089*Pr**(1/3)*(h/p)**0.3
#     f = 29.1955*Re**(-0.5223)*(h/p)**0.864
    
#     h_conv = (Nu*k)/D_h
#     DP = f*L_c*rho*v**2/(2*D_h)
    
#     # Cold side HX area (used to evaluate heat transfer)
#     A_HX = 1/(1+R_p)*N_c*N_p*(1 + np.pi/2)*D_c*L_c
    
#     # t_e : conduction thickness
#     t_e = ((D_c + t_3)*(D_c/2 + t_2) - (np.pi*D_c**2))/(D_c + t_3)
    
#     return h_conv, DP


def PCHE_Lee(alpha, D_c, G, k, L_c, mu, Pr, rho):
    """
    Reference
    ---------
    1) Design and Dynamic Modeling of Printed Circuit Heat Exchangers for Supercritical Carbon Dioxide Brayton Power Cycles (2018) 
    
    Yuan Jiang, Eric Liese, Stephen E. Zitney, Debangsu Bhattacharyya
    
    2) Thermal-hydraulic performance analysis of zigzag channel PCHE according to bending angle (2024)
    
    Yoomyeong Lee, Seongmin Lee, Hong Beom Park, Kyoung Woo Seo, Donghwi Lee
    """
    
    D_h = np.pi*D_c/(2+np.pi)
    Re = (G*D_h)/mu # Re in 0.5 to 2000
     
    if Re <= 2000:
        if alpha <= 10:
            p = 12.93*1e-3 # m
            h = 1.14*1e-3 # m
        elif alpha <= 20:
            p = 12.09*1e-3 # m
            h = 2.2*1e-3 # m
        elif alpha <= 30:
            p = 10.95*1e-3 # m
            h = 3.16*1e-3 # m
        elif alpha <= 40:
            p = 9.53*1e-3 # m
            h = 4*1e-3 # m
        
        Nu = 0.0479*Re**0.7089*Pr**(1/3)*(h/p)**0.3
        f = 29.1955*Re**(-0.5223)*(h/p)**0.864
            
    else:
        # Nu = 0.0176*Re**0.809*Pr**0.33
        
        Nu = 0.0845*Re**0.721*Pr**0.33
        
        # Nu = 0.0292*Re**0.8138        
    
    h_conv = (Nu*k)/D_h

    # DP = f*L_c*rho*v**2/(2*D_h)
    
    return h_conv

