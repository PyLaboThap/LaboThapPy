# -*- coding: utf-8 -*-
"""
Created on Fri May 31 14:52:28 2024

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

class PumpPolynEff(BaseComponent):
    def __init__(self):

        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_pp = WorkConnector()
        self.N_pp = None

#%%
    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return ['su_p', 'su_T', 'ex_p', 'su_fluid', 'N_pp']

    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        elif self.su.h is not None:
            self.inputs['su_h'] = self.su.h
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p
        if self.N_pp is not None:
            self.inputs['N_pp'] = self.N_pp

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        elif 'su_h' in self.inputs:
            self.su.set_h(self.inputs['su_h'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])
        if 'N_pp' in self.inputs:
            self.W_pp.set_N(self.inputs['N_pp'])
            self.N_pp = self.inputs['N_pp']

    def get_required_parameters(self):
        return ['min_flowrate', 'rated_flowrate', 'max_flowrate', 'N_pp_rated', 'pump_voltage', 'pump_phases', 'eta_v', 'V_swept',
                'eta_max_motor', 'W_dot_el_rated', 'coefs_pump', 'eta_tot', 'eta_m']

    def print_setup(self):
        print("=== Pump Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nInputs:")
        for input_name in self.get_required_inputs():
            if input_name in self.input_values:
                print(f"  - {input_name}: {self.input_values[input_name]}")
            else:
                print(f"  - {input_name}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")
        print("======================")
    
    def eta_el(self,P_ratio):
        coefs_20p = [78.74503721229180,1.54402269709448,-0.04662069008665,0.00069559243591,-0.00000499382422,0.00000001349770] # !!! Aribtrary 
        coefs_m20 = [0.82025554862776,4.78234707015054,0.73842411551209,-0.06398392686793,0.00134594665523]
        
        eta_ratio = 0
        
        if P_ratio >= 20:
            for i in range(len(coefs_20p)):
                eta_ratio = eta_ratio + coefs_20p[i]*P_ratio**i
        else:
            for i in range(len(coefs_m20)):
                eta_ratio = eta_ratio + coefs_m20[i]*P_ratio**i
    
        return eta_ratio*self.params['eta_max_motor']
    
    
    def solve(self):
        
        self.check_calculable()
        self.check_parametrized()
        
        if self.calculable and self.parametrized:

            #MODELLING PART
            if self.su.p < self.ex.p and self.N_pp > 0:      
                
                "Flowrate"
                self.V_dot = self.N_pp*self.params['V_swept']*self.params['eta_v'] # m^3/s
                self.m_dot = self.V_dot*self.su.D # kg/s

                self.su.set_m_dot(self.m_dot) 
                self.ex.set_m_dot(self.m_dot) 
                
                "Power and Outlet State Computation"

                self.eta_tot = self.params['eta_tot'](self.V_dot*3600)/100
                p_ex = self.ex.p
                
                def res_compute_pump_power(eta_el,p_ex):
                    self.eta_is = self.eta_tot / (self.params['eta_m'] * self.params['eta_v'] * eta_el)
                                        
                    h_ex_s = PropsSI('H', 'P',p_ex, 'S', self.su.s, self.su.fluid)
                    h_ex = self.su.h + (h_ex_s - self.su.h) / self.eta_is
                    
                    self.ex.set_h(h_ex) 
                    
                    self.W_dot_f = (h_ex - self.su.h) * self.su.m_dot
                    self.W_dot_el = self.W_dot_f / self.eta_tot
                    
                    eta_el_g = self.eta_el(100*self.W_dot_el / self.params['W_dot_el_rated'])/100

                    return eta_el - eta_el_g
                
                # Define the lambda function
                f = lambda eta_el: res_compute_pump_power(eta_el,p_ex)
                
                # Use fsolve to find the root
                initial_guess = self.params['eta_max_motor']
                out_fsolve, info, ier, msg = fsolve(f, initial_guess, full_output=True)

                self.W_pp.set_W_dot(self.W_dot_f)

                # # Outlet State
                # h_ex_s = PropsSI('H', 'P', self.point_ex[0].p,'S',self.point_su[0].s,self.point_su[0].fluid)                
                # h_ex = self.point_su[0].h + (h_ex_s - self.point_su[0].h)/self.geom.eta_is
                
                # self.W_dot_f = (h_ex - self.point_su[0].h)*self.point_su[0].m_dot
                
                
                
                # self.W_dot_el = self.W_dot_f/self.eta_tot
                                    
            # self.point_ex.set_p(P_ex)
            self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")
        



