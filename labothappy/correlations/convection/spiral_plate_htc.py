# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:01:15 2026

@author: basil
"""

def single_phase_spiral_plate_htc(mu, Pr, mu_w, k, G, Dh, R_curv, ch_w):
    
    curv_ratio = Dh/(2*R_curv)
    Re = G*Dh/mu
    Re_cr = 12730*curv_ratio**0.2 # Ito - Reynolds for laminar - turbulent transition depending on curvature
    
    cp = Pr*k/mu
    
    if Re <= Re_cr: # Laminar
        h = cp*G*(1+3.5*(curv_ratio))* (0.023*Re**(-0.2) * Pr**(-2/3))
    else: # Turbulent
        h = 1.86*cp*G*Re**(-2/3) * Pr**(-2/3) * (ch_w/(2*R_curv))**(1/6)

    return h

