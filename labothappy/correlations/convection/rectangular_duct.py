# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:36:57 2026

@author: basil
"""

from labothappy.correlations.convection.pipe_htc import gnielinski_pipe_htc

def rect_duct_1P_htc(mu, Pr, mu_w, k, G, Dh, L, AR):
    
    Re = G*Dh/mu

    if Re < 2300: # Laminar
        # Shah and London (1978)
        Nu = 7.541*(1 - 2.61*AR + 4.97*AR**2 - 5.119*AR**3 + 2.702*AR**4 - 0.548*AR**5)
        h = Nu*k/Dh
        
    elif Re > 10000: # as in pipes
        h = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
    
    else: # Transition - interpolation with turbulent 
        Nu_lam = 7.541*(1 - 2.61*AR + 4.97*AR**2 - 5.119*AR**3 + 2.702*AR**4 - 0.548*AR**5)
        h_lam = Nu_lam*k/Dh

        h_turb = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, Dh, L)
        slope = (h_turb - h_lam)/(10000 - 2300) 
        h = h_lam + slope*(Re-2300)
    
    return h

