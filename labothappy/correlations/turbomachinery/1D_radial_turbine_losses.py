# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 09:04:20 2025

@author: Basile
"""

import numpy as np

def nozzle_losses(v4, Re_N, alpha4, s, c, b4):
    """
    v4 : Rotor inlet speed [m/s]
    Re_N : 
    alpha4 : 
    s : 
    c : 
    b4 : 
    """    
    Dh_n = (v4**2 / 2) * (0.05/Re_N**0.2) * ((3*np.tan(alpha4))/(s/c) + (s*np.cos(alpha4))/b4) 
    
    return Dh_n



