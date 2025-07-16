# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 11:33:04 2025

@author: basil
"""

import numpy as np
from scipy.interpolate import interp1d

def tube_sheet_thickness(D_o, pitch, T_shell, P_des, G):
    """
    Inputs
    ----------
        - D_o : Input outer tube diameter [m]
        - pitch : Tube Pitch [m]
        - T_shell : Input shell temperature [K]
        - P_des : Design shell Pressure [Pa]
        - G : Gasket Diameter [m]
    
    Outputs
    -------
        - t_tube_sheet : Tubesheet thickness preventing bending [m]
        
    Reference
    ---------
    Mechanical Design of Shell and Tube Type Heat Exchanger as per ASME Section VIII Div.1 and TEMA Codes for Two Tubes
    
    """
    T_S_interp = np.array([0, 93.33, 204.444, 315.556, 371.111,
                          398.889, 426.667]) + 273.15  # [K] : Temperature vector
    # [MPa] : Max stress with respect to temperature vector
    S_interp = np.array([158.57942, 158.57942, 158.57942,
                        134.44777, 131.00039, 103.42136, 82.737088])
    
    S_fun = interp1d(T_S_interp, S_interp, kind='linear')
    
    """
    
    Max allowable pressure depending on pipe outside diameter and thickness
    If under critical pressure, associated saturation temperature
    
    """

    "Compute P_max for inputs"
    
    S_tube_calc = S_fun(T_shell)*1e6  # [Pa]
    
    F = 1 # floating tubesheet
    
    eta = 1 - 0.785/(pitch/D_o)**2 # Square tube pattern
    
    t_tube_sheet = (F*G/3)*np.sqrt(P_des/(eta*S_tube_calc))
    
    return t_tube_sheet