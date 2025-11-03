# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 17:13:48 2025

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP


def PCHE_conv(alpha, D_c, G, k, L_c, mu, mu_w, Pr, T, type_channel):
    """
    
    type_channel : Straight / Zigzag / Wavy / Airfoil
    
    
    Reference
    ---------
    Evaluation of thermal-hydraulic performance and economics of Printed Circuit Heat Exchanger (PCHE) for 
    recuperators of Sodium cooled Fast Reactors (SFRs) using CO2 and N2 as working fluids (2022)

    Su Won Lee a, Seong Min Shin a, SungKun Chung b, HangJin Jo a,b,
    """

    Dh = np.pi*D_c/(2+np.pi)
    Re = (G*Dh)/mu 
    
    m_dot = G*(np.pi*(D_c**2)/8)
    
    if type_channel == "Straight":
         def Stephan_Preusser(Re, Pr, Dh, L):
             Nu_1 = 0.0677*(Re*Pr*(Dh/L)**1.33)
             Nu_2 = 1 + 0.1*Pr*(Re*Dh/L)**(0.3)
             Nu = 3.657 + Nu_1/Nu_2
             return Nu
         def gnielinski_transition(Re, Pr):
             f = (1.82*np.log10(Re) - 1.64)**(-2) # (1.8*log10(Re) - 1.5)**(-2)
             Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L_c)**(2/3)) #*(Pr/Pr_w)**(0.11)
             return Nu
         def sieder_tate(Re, Pr):
             Nu = 0.027*Re**0.8*Pr**(1/3)*(mu/mu_w)**0.14
             return Nu
         
         #-------------------------------------------------------------------------
         Re_min = 0
         Re_max = 1e06
         Re = G*Dh/mu

         #-------------------------------------------------------------------------
         # if Re > 10000: #fully turbulent
         #     Pr_min = 0.1
         #     Pr_max = 1000
         #     Nu = sieder_tate(Re, Pr)
         if Re < 2300: #fully laminar
             Pr_min = 0.6
             Pr_max = np.inf
             Nu = Stephan_Preusser(Re, Pr, Dh, L_c)
         else: #transition zone    
             Pr_min = 0.1
             Pr_max = 1000
             Nu = gnielinski_transition(Re, Pr)
             #-------------------------------------------------------------------------
         
    elif type_channel == "Zigzag":
        # Angle akpha in degrees
        Nu = 2.2631*1e6*Re**(1.0515**1e-3)*m_dot**(6.5405*1e-1)*alpha**(-2.0665*1e-4)*(L_c/D_c)**(8.0885*1e-5)*2.2967*1e6*m_dot**(6.5429*1e-1)*(T-273.15)**(-1.3319*1e-3)
        
    hConv = Nu*k/Dh
    
    
    return hConv


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

