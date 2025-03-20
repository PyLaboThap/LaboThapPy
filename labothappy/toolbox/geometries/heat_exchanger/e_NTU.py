# -*- coding: utf-8 -*-
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
