# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:19:26 2026

@author: Basile
"""

import numpy as np
from labothappy.correlations.convection.pipe_htc import gnielinski_pipe_htc

def curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    """
    Performance Evaluation of a Helical Coil Heat
    Exchanger Working under Supercritical Conditions in
    a Solar Organic Rankine Cycle Installation
    +
    Fluid Flow in Helically Coiled Pipes
    """
    
    curv_ratio = Dh/2*R_curv
    Re = G*Dh/mu
    
    
    if Re < 20000: # Xin and Ebadian
        De = Re*np.sqrt(curv_ratio) # Dean number
        Nu = (2.153 + 0.318*De**0.643)*Pr**0.177
        h = Nu*k/Dh
        
    if Re >= 20000:
        h_str = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
        corr = 1+3.6*(1-curv_ratio)*(curv_ratio**0.8) # Nellis and Klein
        h = h_str*corr 
        
    return h

def curved_rect_duct_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv):
    
    curv_ratio = Dh/2*R_curv
    Re = G*Dh/mu

    if Re < 2300:
        Nu = 0
        h_str = 0
    else:
        # If turbulent, it does not depend on the 
        h = curved_pipe_1P_htc(mu, Pr, mu_w, k, G, Dh, L, R_curv)

    return h

