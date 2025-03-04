# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 14:39:18 2025

@author: Basile
"""

import numpy as np
from scipy.interpolate import interp1d

def shell_thickness(D_i, T_shell, P_des):
    """
    Inputs
    ----------
        - D_i : Inner diameter [m]
        - T_shell : Shell temperature [K]
        - P_des : Design Pressure [Pa]
    
    Outputs
    -------
        - t : Minimum allowable thickness [m]
        
    Reference
    ---------
    2007 ASME BPV Code 
    
    """
    T_S_interp = np.array([0, 93.33, 204.444, 315.556, 371.111,
                          398.889, 426.667]) + 273.15  # [K] : Temperature vector
    # [MPa] : Max stress with respect to temperature vector
    S_interp = np.array([158.57942, 158.57942, 158.57942,
                        134.44777, 131.00039, 103.42136, 82.737088])
    
    S_fun = interp1d(T_S_interp, S_interp, kind='linear')
    
    """
    
    Max allowable internal pressure depending on pipe outside diameter and thickness
    
    """

    "Compute P_max for inputs"
    
    S_tube_calc = S_fun(T_shell)*1e6  # [Pa]

    from scipy.optimize import fsolve
    
    # Define the equation as a function
    def equation(t, P_des, S_tube_calc, D_i):
        return P_des - S_tube_calc * ((2*t - 0.01*D_i) / (D_i - (t - 0.005*D_i)))
    
    # Solve for t
    t_initial_guess = 0.001
    t_solution = fsolve(equation, t_initial_guess, args=(P_des, S_tube_calc, D_i))
    return t_solution[0]
