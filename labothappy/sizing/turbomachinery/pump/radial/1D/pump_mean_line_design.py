# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 12:44:38 2025

@author: basil
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, root_scalar

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Omega = 10000 # [RPM]

class RadialPumpMeanLineDesign():
    
    def __init__(self, fluid):
        
        self.fluid = fluid
        
        self.inputs = {}
        self.params = {}
        
        self.total_states  = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3,4,5])
        self.static_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3,4,5])
        self.AS = CP.AbstractState('HEOS', fluid)
    
    def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.total_states['H'][position] = self.AS.hmass()            
        self.total_states['S'][position] = self.AS.smass()            
        self.total_states['P'][position] = self.AS.p()            
        self.total_states['D'][position] = self.AS.rhomass()            

        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1  
            
        self.total_states['V'][position] = self.AS.viscosity()            
        
        return
    
    def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.static_states['H'][position] = self.AS.hmass()            
        self.static_states['S'][position] = self.AS.smass()            
        self.static_states['P'][position] = self.AS.p()            
        self.static_states['D'][position] = self.AS.rhomass()    
        
        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1            
            
        self.static_states['V'][position] = self.AS.viscosity()            

        return

#%% 
    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value

#%%
    # ---------------- Design related metods --------------------------------------------------------------

    def design(self):
        self.update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        self.update_static_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        
        "------------- 1) Isentropic Pumping Calculation -----------------------------------------------" 
        s_in = self.total_states['S'][1]
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s_in)
        
        h_is_ex = self.AS.hmass()
        Dh0s = h_is_ex - self.total_states['H'][1]
                
        self.Dh0 = self.inputs['W_dot']/self.inputs['mdot']
        self.eta_is = self.Dh0/Dh0s

        hout = self.total_states['H'][1] - Dh0s * self.eta_is
        
        self.AS.update(CP.HmassP_INPUTS, hout, self.inputs['p_ex'])
        rho_out = self.AS.rhomass()
        
        return 0 
        
    
if __name__ == "__main__":
    
    fluid = 'CO2'
    
    Pump = RadialPumpMeanLineDesign(fluid)
    
    Pump.set_inputs(
        mdot = 100, # kg/s
        W_dot = 4.69*1e6, # W
        p_su = 40*1e5, # Pa
        T0_su = 273.15 + 15, # K
        p_ex = 140*1e5, # Pa
        )
    
    Pump.set_parameters(
        psi = 1, # [-]
        phi = 0.4, # [-]
        xhi = 0.4, # [-]
        r5_r4_ratio = 0.5, # [-] : Maybe impose r5/r4 ratio instead - Iterate on it
        )
    
    Pump.design()
    
    