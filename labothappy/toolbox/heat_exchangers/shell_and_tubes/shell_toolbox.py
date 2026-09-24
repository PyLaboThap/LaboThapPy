# -*- coding: utf-8 -*-
# LaboThApPy: An open-source tool for advanced thermodynamic cycle simulation, designed for researchers and engineers.
#
# Copyright (C) <2025> <Université catholique de Louvain (UCLouvain), Belgique
#                       Université de Liège (ULiège), Belgique
#                       Université de Mons (UMONS), Belgique>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# List of the contributors to the development of LaboThApPy: see AUTHORS file.
# Description and complete License: see NOTICE & LICENSE files.

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

    P_des = max(P_des, 7*1e5)

    "Compute P_max for inputs"
    
    S_tube_calc = S_fun(T_shell)*1e6  # [Pa]

    from scipy.optimize import fsolve
    
    # Define the equation as a function
    def equation(t, P_des, S_tube_calc, D_i):
        return P_des - S_tube_calc * ((2*t - 0.01*D_i) / (D_i - (t - 0.005*D_i)))
    
    # Solve for t
    t_initial_guess = 0.001
    t_solution = fsolve(equation, t_initial_guess, args=(P_des, S_tube_calc, D_i))
    
    corrosion_allowance = 1.6*1e-3 # mm as per TEMA code
    
    return t_solution[0] + corrosion_allowance


