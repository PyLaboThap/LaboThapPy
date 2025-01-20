# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:46:33 2024

@author: Basile
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d

class Spliter_Geom(object):
    
    def __init__(self, **kwargs):
        
        self.n_inlet = None
        self.V = None # Volume [m^3]

    def set_parameters(self, name, **kwargs):
        
        if name == "DECAGONE_ACC_Spliter":
            self.V = 2 # [m^3]
            self.flow_coef = [0.5,0.5]
           