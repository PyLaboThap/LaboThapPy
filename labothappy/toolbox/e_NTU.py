# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:25:34 2025

@author: marie
"""
import numpy as np

def e_NTU(self, NTU, C_r):
    
    if self.params['Flow_Type'] == "CounterFlow":
        eps = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))    
        return eps
    
    