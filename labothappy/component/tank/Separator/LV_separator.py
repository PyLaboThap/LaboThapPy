# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""

from Port.Mass_connector import Mass_connector
from CoolProp.CoolProp import PropsSI
import numpy as np

class LV_separator(object):
    
    class geom():
            pass 

    def __init__(self):
        
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.defined = False
        
        "Geometry sub-object"
        self.geom = None # parameters 
        
        "Input"
        self.point_su = [None]

        "Outputs"
        self.point_ex = [None, None]
                
#%%    
    def check_calculable(self):
        if self.point_su[0] != None:
            if (self.point_su[0].completely_known):
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
        
        "1) Create output mass connectors from the input"
                
        for i in range(len(self.point_ex)):
            self.point_ex[i] = Mass_connector()  # Replace the element with a Mass_connector object
            
            self.point_ex[i].set_fluid(self.point_su[0].fluid)
            self.point_ex[i].set_p(self.point_su[0].p)
            
        self.point_ex[0].set_x(1)
        self.point_ex[1].set_x(0)
        
        "2) Mass Flowrates"

        if self.point_su[0].x >= 0:  # LV Saturation
            m_dot_l = (1-self.point_su[0].x) * self.point_su[0].m_dot
            m_dot_v = self.point_su[0].x * self.point_su[0].m_dot

            self.point_ex[0].set_m_dot(m_dot_v)
            self.point_ex[1].set_m_dot(m_dot_l)            
        
        else:  # 1 phase
            T_sat = PropsSI('T','P',self.point_su[0].p,'Q',0.5,self.point_su[0].fluid)
            
            if self.point_su[0].T > T_sat: # Vapor phase
                self.point_ex[0].set_m_dot(self.point_su[0].m_dot)
                self.point_ex[1].set_m_dot(0)                
            else:
                self.point_ex[1].set_m_dot(self.point_su[0].m_dot)
                self.point_ex[0].set_m_dot(0)                  

        # print("su :", self.point_su[0].m_dot)
        # print("ex1 :", self.point_ex[0].m_dot)
        # print("ex2 :", self.point_ex[1].m_dot)

        self.defined = True    
        
        return
        
            
            
            