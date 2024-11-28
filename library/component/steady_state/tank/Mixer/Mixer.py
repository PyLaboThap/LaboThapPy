# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
import numpy as np
from component.base_component import BaseComponent

class Mixer(BaseComponent): 
    """
    Component: Mixer
    -----------------------------------------------------------
    Connectors:
        su (array of MassConnector): Mass connectors for the supply side.
        ex (MassConnector): Mass connector for the exhaust side.
    -----------------------------------------------------------
    Parameters:
        n_inlet : Number of supply inlets
    -----------------------------------------------------------
    Inputs:
        su_p: Supply side pressures (shall all be the same).
        su_h: Supply side enthalpies.
        su_m_dot: Supply side flowrates.
        su_fluid: Suction side fluids (no mixture for now, shall be the same).
    -----------------------------------------------------------
    Ouputs:
        ex_p: Exhaust side pressure (same as inlet pressures).
        ex_h: Exhaust side enthalpy (from conservation of internal energy).
        ex_m_dot: Exhaust side flowrate (sum of all inlet flowrates).
        ex_fluid: Exhaust side fluid (no mixture for now).
    """
    def __init__(self, n_inlets):
                
        "Input"
        self.su = np.full(n_inlets, MassConnector())

        self.params['n_inlets'] = n_inlets

        "Outputs"
        self.ex = MassConnector()

    def get_required_inputs(self):
        if self.inputs == {}:
            for i in range(self.params['n_inlets']):
                if self.su[i].fluid is not None:
                    self.inputs['su_'+ str(i) +'_fluid'] = self.su[i].fluid
                if self.su[i].T is not None:
                    self.inputs['su_'+ str(i) +'_T'] = self.su[i].T
                elif self.su[i].h is not None:
                    self.inputs['su_'+ str(i) +'_h'] = self.su[i].h
                if self.su[i].p is not None:
                    self.inputs['su_'+ str(i) +'_p'] = self.su[i].p
                if self.su[i].m_dot is not None:
                    self.inputs['su_'+ str(i) +'_m_dot'] = self.su[i].m_dot
        
        if self.inputs != {}:
            for i in range(self.params['n_inlets']):
                self.su[i].set_fluid(self.inputs['su_'+ str(i) +'_fluid'])
                if ('su_'+ str(i) +'_T') in self.inputs:
                    self.su[i].set_T(self.inputs['su_'+ str(i) +'_T'])
                elif ('su_'+ str(i) +'_h') in self.inputs:
                    self.su[i].set_h(self.inputs['su_'+ str(i) +'_h'])
                if ('su_'+ str(i) +'_p') in self.inputs:
                    self.su[i].set_p(self.inputs['su_'+ str(i) +'_p'])
                if ('su_'+ str(i) +'_m_dot') in self.inputs:
                    self.su[i].set_m_dot(self.inputs['su_'+ str(i) +'_m_dot'])

        req_inputs = []
        for i in range(self.params['n_inlets']):
            req_inputs.append('su_'+str(i)+'_p')
            req_inputs.append('su_'+str(i)+'_h')
            req_inputs.append('su_'+str(i)+'_fluid')
            req_inputs.append('su_'+str(i)+'_m_dot')

        return req_inputs
    
    def get_required_parameters(self):
        return [
            'n_inlets',
        ]
    
    def print_setup(self):
        print("=== Mixer Setup ===")
        print("Connectors:")

        for i in range(self.params['n_inlets']):
            print(f"  - su ({i}): fluid={self.su[i].fluid}, T={self.su[i].T}, p={self.su[i].p}, m_dot={self.su[i].m_dot}")

        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nInputs:")
        for input in self.get_required_inputs():
            if input in self.inputs:
                print(f"  - {input}: {self.inputs[input]}")
            else:
                print(f"  - {input}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")

        print("======================")

#%%

    def solve(self):
        
        "1) Check if pressures are the same"
        
        for i in range(len(self.point_su)):
            if self.point_su[0].p != self.point_su[i].p:
                print("Not the same pressure at all mixer inputs")
                return
        
        "2) Compute output"
        
        self.point_ex[0] = MassConnector()
        self.point_ex[0].set_fluid(self.point_su[0].fluid)
        self.point_ex[0].set_p(self.point_su[0].p)
        
        m_dot_out = 0
        
        for i in range(len(self.point_su)):
            m_dot_out = m_dot_out + self.point_su[i].m_dot
            
        self.point_ex[0].set_m_dot(m_dot_out)
        
        h_out_mean = 0
        
        for i in range(len(self.point_su)):
            h_out_mean = h_out_mean + self.point_su[i].m_dot*self.point_su[i].h
        
        h_out_mean = h_out_mean/m_dot_out
        
        self.point_ex[0].set_h(h_out_mean)
        
        self.defined = True
        
        return
        
        
            