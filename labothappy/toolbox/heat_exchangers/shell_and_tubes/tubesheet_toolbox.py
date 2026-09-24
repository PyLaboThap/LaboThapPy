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
    
    corrosion_allowance = 1.6*1e-3 # m
    
    return t_tube_sheet + corrosion_allowance