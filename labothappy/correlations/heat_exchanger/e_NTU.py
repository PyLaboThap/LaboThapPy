# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:17:12 2025

@author: marie
"""
import numpy as np

def e_NTU(NTU, C_r, params):
    """
    Sources:
    S. Kakaç, H. Liu, and A. Pramuanjaroenkij, Heat Exchangers, 0 edn. CRC Press, 2012. doi: 10.1201/b11784.
    Page 60
    
    G. F. Nellis and S. A. Klein, Heat transfer, 1. paperback ed. Cambridge: Cambridge University Press, 2012.
    Page 856
    """
    
    # ADD if C_r = 1 ??
    
    
    if params['Flow_Type'] == "CounterFlow":
        if C_r > 0.999: # better C_r == 1.00?
            eps = NTU / (1+ NTU)
        else:
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
        
    elif params['Flow_Type'] == "one_fluid" or C_r < 1e-5 : # better write C_r == 0 ?
        eps = 1 - np.exp(-NTU)

    
            
    else:
        raise ValueError(f"Flow_Type '{params['Flow_Type']}' not recognized or not implemented")
    
    return eps



#Correlation from Fundamentals Of Heat And Mass Transfer Frank P Incropera
