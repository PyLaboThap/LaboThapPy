# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:17:12 2025

@author: marie
"""
import numpy as np



# Side functions first
def Crossflow_Cmax_mixed_Cmin_unmixed(NTU, C_r):
    eps = (1 - np.exp(C_r* (np.exp(-NTU)-1))) / C_r
    return eps

def Crossflow_Cmin_mixed_Cmax_unmixed(NTU, C_r):
    eps = 1 - np.exp(-(1 - np.exp(-C_r * NTU)) / C_r)
    return eps
    
def Crossflow_both_mixed(NTU, C_r):
    eps = (1 / C_r) * (1 - np.exp(-C_r * (1 - np.exp(-NTU))))
    return eps
    
def Crossflow_both_unmixed(NTU, C_r):
    eps = 1 - np.exp((1 / C_r) * (NTU ** 0.22) * (np.exp(-C_r * (NTU ** 0.78)) - 1))
    return eps
    

# Main function
def e_NTU(NTU, C_r, params,  C_c = None, C_h = None):
    """
    Sources:
    S. Kakaç, H. Liu, and A. Pramuanjaroenkij, Heat Exchangers, 0 edn. CRC Press, 2012. doi: 10.1201/b11784.
    Page 60
    
    G. F. Nellis and S. A. Klein, Heat transfer, 1. paperback ed. Cambridge: Cambridge University Press, 2012.
    Page 856
    """
    
    
    
    if params['Flow_Type'] == "CounterFlow":
        if C_r > 0.999: # if C_r = 1
            eps = NTU / (1+ NTU)
        else:
            eps = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))    
    
    
    elif params['Flow_Type'] == "ParallelFlow":
        eps = (1 - np.exp(-NTU * (1 + C_r))) / (1 + C_r) 
        
        
    elif params['Flow_Type'] == "CrossFlow_Unmixed":        
        eps = Crossflow_both_unmixed(NTU, C_r)

        
    elif params['Flow_Type'] == "CrossFlow_Mixed":
        eps = Crossflow_both_mixed(NTU, C_r)
        
        
        
    elif params['Flow_Type'] == "CrossFlow":
        
        
        # If both are mixed. Equivalent to setting 'CrossFlow_Mixed'
        if params['H_mixing'] == "Mixed" and params['C_mixing'] == "Mixed":
            eps = Crossflow_both_mixed(NTU, C_r)
        
        # If both are unmixed. Equivalent to setting 'CrossFlow_Unmixed'
        elif params['H_mixing'] == "Unmixed" and params['C_mixing'] == "Unmixed":
            eps = Crossflow_both_unmixed(NTU, C_r)
        
        elif C_c ==  None or C_h == None:
            raise ValueError("When using 'CrossFlow' with the hot and cold flow having different mixing, please call the function eNTU with the cold and hot heat capacity rate of each flow C_dot_cold ('C_c') and C_dot_hot ('C_h').")
        
        elif params['H_mixing'] == "Unmixed" and params['C_mixing'] == "Mixed":
            if C_c > C_h:
                eps = Crossflow_Cmax_mixed_Cmin_unmixed(NTU, C_r)
            elif C_c < C_h:
                eps = Crossflow_Cmin_mixed_Cmax_unmixed(NTU, C_r)  
            else:
                raise ValueError("C_h and C_c are equal. Impossible to find a proper correlation")
              
        elif params['H_mixing'] == "Mixed" and params['C_mixing'] == "Unmixed":
            if C_c < C_h:
                eps = Crossflow_Cmax_mixed_Cmin_unmixed(NTU, C_r)
            elif C_c > C_h:
                eps = Crossflow_Cmin_mixed_Cmax_unmixed(NTU, C_r)
            else:
                raise ValueError("C_h and C_c are equal. Impossible to find a proper correlation")    
            
        else:
            raise ValueError("If 'Flow_Type' = 'CrossFlow' selected, please provide the mixing of the fluids via the parameters 'C_mixing' and 'H_mixing'. Their values shall be equal to 'Unmixed' or 'Mixed'")
            


    elif params['Flow_Type'] == "ShellAndTube_1_2":
        eps = 2 * (1 + C_r + np.sqrt(1 + C_r**2) * (1 + np.exp(-NTU * np.sqrt(1 + C_r**2))) / (1 - np.exp(-NTU * np.sqrt(1 + C_r**2))))**-1

    elif params['Flow_Type'] == "ShellAndTube_n_passes":
        eps_1 = 2 / (1 + C_r + np.sqrt(1 + C_r**2)) * (1 - np.exp(-NTU * np.sqrt(1 + C_r**2)))
        n = params["n_shell_pass"]
        eps = (( (1 - eps_1 * C_r) / (1 - eps_1) )**n - 1) / ( ( (1 - eps_1 * C_r) / (1 - eps_1) )**n - C_r)
        
    elif params['Flow_Type'] == "one_fluid" or C_r < 1e-5 : # if C_r == 0 
        eps = 1 - np.exp(-NTU)

    
            
    else:
        raise ValueError(f"Flow_Type '{params['Flow_Type']}' not recognized or not implemented. The implemented flow types are:\n - CounterFlow\n - ParallelFlow\n - CrossFlow\n - CrossFlow_Unmixed\n - CrossFlow_Mixed\n - ShellAndTube_1_2\n - ShellAndTube_n_passes\n - one_fluid")
    
    return eps



#Correlation from Fundamentals Of Heat And Mass Transfer Frank P Incropera
