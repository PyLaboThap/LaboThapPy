# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""

from Port.Mass_connector import Mass_connector
from CoolProp.CoolProp import PropsSI
import numpy as np

class Spliter(object):
    
    class geom():
            pass 

    def __init__(self, geom):
        
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.defined = False
        
        "Geometry sub-object"
        self.geom = geom # parameters 
        
        "Input"
        self.point_su = [None]

        "Outputs"
        self.point_ex = np.full(len(geom.flow_coef), None)
                
#%%    
    def check_calculable(self):
        if self.point_su[0] != None:
#            if (self.point_su[0].completely_known):
             self.calculable = True
        
    def check_parametrized(self):
        if self.geom != None:
            self.parametrized = True
           
#%%
    def set_parameters(self, **kwargs):
            """
            Set parameters of the heat exchanger.
    
            Parameters
            ----------
            **kwargs : dict
                Key-value pairs representing parameters and their values.
                
                Example of call : heat_exchanger.set_parameters(D_o=0.05, t=0.002, H_core=2.0)
            """
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger.")
            
            self.check_parametrized()
            
#%%

    def solve(self):
        
        "1) Compute output"
        
        for i in range(len(self.point_ex)):
            self.point_ex[i] = Mass_connector()
            self.point_ex[i].set_fluid(self.point_su[0].fluid)
            self.point_ex[i].set_p(self.point_su[0].p)
            self.point_ex[i].set_h(self.point_su[0].h)
            self.point_ex[i].set_m_dot(self.point_su[0].m_dot*self.geom.flow_coef[i])
        
        self.defined = True
        
        return
        
        
            