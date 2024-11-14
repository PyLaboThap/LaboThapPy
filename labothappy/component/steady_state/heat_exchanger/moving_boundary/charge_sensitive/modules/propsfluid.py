# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 11:38:04 

Modification w/r to previous version:
    -T_mean is arbitrarely set to the critical temperature in order to avoid errors.

@author: jvega
"""

import CoolProp.CoolProp as CP
import numpy as np
from scipy.interpolate import interp1d

def propsfluid(T_mean, P_mean, T_wall, fluid, incompr_flag):
    
    #Force T_wall to be under TCRIT
    # if not incompr_flag:
    #     T_crit = CP.PropsSI("TCRIT", fluid)
    # T_mean = min(T_mean, T_crit)
    # T_wall = min(T_wall, T_crit)
    #-----------------------------------------------------------------------
    mu = CP.PropsSI('V',        'T', T_mean, 'P', P_mean, fluid)
    Pr = CP.PropsSI('Prandtl',  'T', T_mean, 'P', P_mean, fluid)
    k = CP.PropsSI('L',        'T', T_mean, 'P', P_mean, fluid)
    mu_wall = CP.PropsSI("V", "T", T_wall, "P", P_mean, fluid)
    mu_rat = mu/mu_wall
    
    return mu, Pr, k, mu_wall, mu_rat, 0, 0
        