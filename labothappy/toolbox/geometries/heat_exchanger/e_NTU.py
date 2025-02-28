# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:17:12 2025

@author: marie
"""
import numpy as np

def e_NTU(NTU, C_r, Flow_Type):
    
    if Flow_Type == "CounterFlow":
        eps = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))    
        return eps
    
    
#TODO : add some configurations (MP)  