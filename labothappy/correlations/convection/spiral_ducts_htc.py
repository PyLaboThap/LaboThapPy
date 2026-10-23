# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:19:26 2026

@author: Basile
"""

import numpy as np
from labothappy.correlations.convection.pipe_htc import gnielinski_pipe_htc

def curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    """
    Fluid Flow in Helically Coiled Pipes (2023)
    
    Leonardo Di G. Sigalotti, Carlos E. Alvarado-Rodríguez and Otto Rendón 
    """
    
    curv_ratio = Dh/(2*R_curv)
    Re = G*Dh/mu
    Re_cr = 12730*curv_ratio**0.2 # Ito - Reynolds for laminar - turbulent transition depending on curvature
    
    if Re < Re_cr: 
        h = laminar_curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv)
        
    elif Re >= 10**5: # Schmidt correction from Nellis and Klein
        h = turbulent_curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv)
    
    else: # Xhan correlation from Sigalotti et. al. : between Re_cr and 10e5
        h = transition_curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv)
    
    return h

def laminar_curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    
    curv_ratio = Dh/(2*R_curv)
    Re = G*Dh/mu
    De = Re*np.sqrt(curv_ratio) # Dean number
    
    if curv_ratio <= 0.01:
        h_str, _, _ = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
        h = h_str
    elif 0.01 <= curv_ratio <= 0.1: # Xhan correlation from Sigalotti et. al.
        Nu = (2.153 + 0.318*De**0.643)*Pr**0.177
        h = Nu*k/Dh
    else: 
        # Dean number only can not determine characteristics of fluid flow and heat transfer in helically coiled tubes with large curvature
        # Zhang et. al.
        Nu = 0.8532*De**0.4822 * curv_ratio ** (-0.03451)
        h = Nu*k/Dh

    return h

def transition_curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    
    curv_ratio = Dh/(2*R_curv)
    Re = G*Dh/mu
    De = Re*np.sqrt(curv_ratio) # Dean number
    
    if curv_ratio <= 0.01:
        h_str, _ , _ = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
        h = h_str 
    elif 0.01 <= curv_ratio <= 0.1: # Xhan correlation from Sigalotti et. al.
        num = 0.1033*np.sqrt(De)
        den = ((1+1.729/De)**0.5 - (1.729/De)**0.5)**3
        f = 16*num/(Re*den) # Friction factor from Ito correlation
        De = Re*np.sqrt(curv_ratio) # Dean number
        Nu = 0.04 * f * Re**1.13 * curv_ratio**0.077 * Pr**0.4
        h = Nu*k/Dh
    elif 0.1 <= curv_ratio <= 0.3:
        # Numerical prediction of turbulent flow and heat transfer in helically coiled pipes - Ivan Di Piazza, Michele Ciofalo
        num = 0.1033*np.sqrt(De)
        den = ((1+1.729/De)**0.5 - (1.729/De)**0.5)**3
        f = 16*num/(Re*den) # Friction factor from Ito correlation
        Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L)**(2/3)) #*(Pr/Pr_w)**(0.11)
        h = Nu*k/Dh
    else: # No other correlations for this
        # Numerical prediction of turbulent flow and heat transfer in helically coiled pipes - Ivan Di Piazza, Michele Ciofalo
        num = 0.1033*np.sqrt(De)
        den = ((1+1.729/De)**0.5 - (1.729/De)**0.5)**3
        f = 16*num/(Re*den) # Friction factor from Ito correlation
        Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L)**(2/3)) #*(Pr/Pr_w)**(0.11)
        h = Nu*k/Dh
            
    return h

def turbulent_curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    
    curv_ratio = Dh/(2*R_curv)
    Re = G*Dh/mu
    De = Re*np.sqrt(curv_ratio) # Dean number
     
    if curv_ratio <= 0.01:
        h_str, _, _ = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
        h = h_str 
    
    elif 0.01 <= curv_ratio <= 0.2: # Schmidt correction from Nellis and Klein
        h_str, _, _ = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
        corr = 1+3.6*(1-curv_ratio)*(curv_ratio**0.8)
        h = h_str*corr 

    elif 0.2 <= curv_ratio <= 0.3:
        # Numerical prediction of turbulent flow and heat transfer in helically coiled pipes - Ivan Di Piazza, Michele Ciofalo
        num = 0.1033*np.sqrt(De)
        den = ((1+1.729/De)**0.5 - (1.729/De)**0.5)**3
        f = 16*num/(Re*den) # Friction factor from Ito correlation
        Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L)**(2/3)) #*(Pr/Pr_w)**(0.11)
        h = Nu*k/Dh
        
    else: # No other correlations for this
        # Numerical prediction of turbulent flow and heat transfer in helically coiled pipes - Ivan Di Piazza, Michele Ciofalo
        num = 0.1033*np.sqrt(De)
        den = ((1+1.729/De)**0.5 - (1.729/De)**0.5)**3
        f = 16*num/(Re*den) # Friction factor from Ito correlation
        Nu = (((f/8)*(Re-1000)*Pr) / (1+12.7*(f/8)**(1/2) * (Pr**(2/3)-1)) )*(1 + (Dh/L)**(2/3)) #*(Pr/Pr_w)**(0.11)
        h = Nu*k/Dh
        
    return h

def curved_rect_duct_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    
    curv_ratio = Dh/2*R_curv
    Re = G*Dh/mu
    De = Re*np.sqrt(curv_ratio) # Dean number

    if Re < 2300: # MODELING LAMINAR HEAT TRANSFER IN A CURVED RECTANGULAR DUCT WITH A COMPUTATIONAL FLUID DYNAMICS CODE - Facao and Oliveira
        if De < 1000:
            Nu = 0.8508 * De**0.5207
        elif De > 1200:
            Nu = 38.5
        else:
            Nu_1000 = 0.8508 * De**0.5207
            Nu_1200 = 38.5
            
            slope = (Nu_1200 - Nu_1000)/(1200-1000)
            Nu = Nu_1000 + slope*(De - 1000)
        
        h = Nu*k/Dh
        
    else:
        # If turbulent, it does not depend on the 
        h = curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv)

    return h

