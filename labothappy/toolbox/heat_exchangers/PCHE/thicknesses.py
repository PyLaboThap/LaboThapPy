# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 17:42:12 2025

@author: basil
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

def PCHE_thicknesses(D_c, P_des, T_des):
    """
    Inputs
    ------
    
    D_c : Channel diameter [m]
    P_des : Design Pressure [Pa]
    T_des : Design fluid temperature [K]

    Outputs
    -------
    t_1 : Inter-channel horizontal thickness [m]
    t_2 : Plate thickness [m]
    
    Reference
    ---------
    Structural Assessment of Printed Circuit Heat Exchangers in Supercritical CO2 
    Waste Heat Recovery Systems for Ship Applications (2022)
    
    WANG Jian, YAN Xinping, LU Mingjian, SUN Yuwei, WANG Jiawei
    """

    E = 0.7 # Diffusion welding factor : 0.7 is a recommended value from source
    
    "Maximum allowable stress - SS 316L : 100000hr lifetime"
        
    T_S_interp = np.array([40, 67, 100, 125, 150, 175, 200, 225, 250, 275,
                          300, 325, 350, 375, 400, 425, 450]) + 273.15  # [K] : Temperature vector
    # [MPa] : Max stress with respect to temperature vector
    S_interp = np.array([116.42, 117.1, 117, 116.6, 115.75, 114.1, 111.1, 107.5, 103.8,
                        100.9, 98, 95.7, 94.2, 92.7, 91.3, 89, 88.3])
    
    S_fun = interp1d(T_S_interp, S_interp, kind='linear')
    
    S_calc = S_fun(T_des)*1e6  # [Pa] : Young Modulus
    
    SE = S_calc* E
    
    "t1 calculation"
    h = D_c/2
    H = D_c
    
    t1_1 = (P_des*h)/(2*SE)
    t1_2 = ((P_des*H)+np.sqrt(P_des**2 * h**2 + 12*SE*P_des*H**2))/(6*SE)
    t1_3 = ((P_des*H)-np.sqrt(P_des**2 * h**2 + 12*SE*P_des*H**2))/(6*SE)
    
    t1 = max(t1_1, t1_2, t1_3)
    
    "t2 calculation"
    t2 = P_des*H/SE
        
    return np.array([t1, t2])*1.5

#%% Test

if __name__ == "__main__":
    # Test case    
    D_c = 2*1e-3
    P_des = 140*1e5
    T_des = 250 + 273.15
    
    t1, t2 = PCHE_thicknesses(D_c, P_des, T_des)
