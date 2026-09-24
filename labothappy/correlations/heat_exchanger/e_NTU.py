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
Created on Fri Feb 28 15:17:12 2025

@author: marie
"""
import numpy as np

def e_NTU(NTU, C_r, params):
    
    if params['Flow_Type'] == "CounterFlow":
        eps = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))    
    
    
    elif params['Flow_Type'] == "ParallelFlow":
        eps = (1 - np.exp(-NTU * (1 + C_r))) / (1 + C_r) 
        
        
    elif params['Flow_Type'] == "CrossFlow_Unmixed":        
        eps = 1 - np.exp((1 / C_r) * (NTU ** 0.22) * (np.exp(-C_r * (NTU ** 0.78)) - 1))

        
    elif params['Flow_Type'] == "CrossFlow_Mixed":
        eps = (1 / C_r) * (1 - np.exp(-C_r * (1 - np.exp(-NTU))))
        
        
    elif params['Flow_Type'] == "ShellAndTube_1_2":
        eps = 2 * (1 + C_r + np.sqrt(1 + C_r**2) * (1 + np.exp(-NTU * np.sqrt(1 + C_r**2))) / (1 - np.exp(-NTU * np.sqrt(1 + C_r**2))))**-1

    elif params['Flow_Type'] == "ShellAndTube_n_passes":
        eps_1 = 2 / (1 + C_r + np.sqrt(1 + C_r**2)) * (1 - np.exp(-NTU * np.sqrt(1 + C_r**2)))
        n = params["n_shell_pass"]
        eps = (( (1 - eps_1 * C_r) / (1 - eps_1) )**n - 1) / ( ( (1 - eps_1 * C_r) / (1 - eps_1) )**n - C_r)

    else:
        raise ValueError(f"Flow_Type '{params['Flow_Type']}' not recognized or not implemented")
    
    return eps



#Correlation from Fundamentals Of Heat And Mass Transfer Frank P Incropera
