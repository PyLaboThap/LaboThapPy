# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:01:37 2026

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI


#%% Circulators

def magna1_100_120_GF():
    
    # 1) Specific Speed
    Omega_rated = 3000/60 # Hz
    Omega_pp_rads = Omega_rated*(2*np.pi) # rad/s
    
    fluid = "Water"
    g = 9.81 # m/s2
    Q = 38.25/3600 # m3/s
    eff_tot = 0.649 # [-]
    eff_mot = 0.95 # [-]
    eff_VFD = 1 # [-]    
    
    eff_hyd = eff_tot/(eff_VFD*eff_mot)
    
    rho_1 = PropsSI("D", "P", 1e5, "T", 273.15+15, fluid)
    
    # Height difference
    h_2 = h_1 = 0
    v_1 = v_2 = 0
    p_2 = 0.912*1e5
    p_1 = 0
    
    DH = (p_2 - p_1)/(rho_1*g) + (h_2 - h_1) + (v_2**2 - v_1**2)/(2*g)
    
    omega_s = (Omega_pp_rads*np.sqrt(Q))/((g*DH)**(3/4))
    omega_s_imp = 2733*omega_s
        
    # Data curve 
    
    
    
    return omega_s_imp


#%%

if __name__ == "__main__":
    Omega = 3000 # RPM
    
    omega_s = magna1_100_120_GF(Omega)




