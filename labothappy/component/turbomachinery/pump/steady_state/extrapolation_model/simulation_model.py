# -*- coding: utf-8 -*-
"""
Created on Fri May 31 14:52:28 2024

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

class PumpExtrapolationModel(BaseComponent):
    def __init__(self):

        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_pp = WorkConnector()

#%%
    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return ['su_p', 'su_T', 'ex_p', 'su_fluid', 'Omega_pp']

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
        if self.W_pp.N is not None:
            self.inputs['Omega_pp'] = self.W_pp.N*60

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
        if 'Omega_pp' in self.inputs:
            self.W_pp.set_N(self.inputs['Omega_pp']/60)
            self.Omega_pp = self.inputs['Omega_pp']

    def get_required_parameters(self):
        return ['Omega_rated','min_flowrate', 'rated_flowrate', 'max_flowrate', 'PI_rated', 'D_p', 'V_dot_curve', 'Delta_H_curve', 'eta_is_curve',
                'NPSH_r_curve', 'eta_m', 'eta_max_motor', 'W_dot_el_rated']

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

#%%

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
        
        g = 9.81 # m/s^2
        
        if not self.calculable:
            print("Component is not calculable")
            return
        
        if not self.parametrized:
            print("Component is not parametrized")
            return

        #MODELLING PART
        if self.su.p > self.ex.p or self.Omega_pp < 0: 
            print("Supply pressure is higher than exhaust pressure or rotation speed is negative")
        else:     

            self.V_dot_flag = 0
            
            "Speed extrapolation of curves"

            self.V_dot_curve = self.params['V_dot_curve']*(self.Omega_pp/self.params['Omega_rated'])               
            self.DH_curve = self.params['Delta_H_curve']*(self.Omega_pp/self.params['Omega_rated'])**2
            self.eta_is_curve = self.params['eta_is_curve']
            self.NPSH_r_curve = self.params['NPSH_r_curve']*(self.Omega_pp/self.params['Omega_rated'])**2
            
            # Interpolation of extrapolated curves
            DH_V = interp1d(self.DH_curve, self.V_dot_curve, kind='linear', fill_value='extrapolate')
            V_eta = interp1d(self.V_dot_curve, self.eta_is_curve, kind='linear', fill_value='extrapolate')
            V_NPSH_r = interp1d(self.V_dot_curve, self.NPSH_r_curve, kind='linear', fill_value='extrapolate')
            
            "Height Difference"
            DP = self.ex.p - self.su.p
            self.DH = DP/(g*self.su.D)
            
            "Flowrate"
            self.V_dot = DH_V(self.DH)
            
            if self.V_dot > self.V_dot_curve[-1] or self.V_dot < self.V_dot_curve[0]:
                self.V_dot_flag = 1
            
            self.m_dot = self.su.D*self.V_dot/3600
            self.su.set_m_dot(self.m_dot)
            self.ex.set_m_dot(self.m_dot)
            
            "Isentropic Efficiency"
            self.eta_is = V_eta(self.V_dot)
            
            "NPSH_r"
            self.NPSH_r = V_NPSH_r(self.V_dot)
            
            "Outlet State Computation"
            
            h_ex_s = PropsSI('H', 'P', self.ex.p,'S',self.su.s,self.su.fluid)                
            h_ex = self.su.h + (h_ex_s - self.su.h)/self.eta_is
            self.ex.set_h(h_ex)
            
            "Power and Outlet State Computation"
            
            self.W_dot_hyd = self.su.m_dot*(h_ex_s - self.su.h)
            
            Delta_h = self.ex.h - self.su.h
            self.W_dot_wf = self.su.m_dot*Delta_h # W
            
            self.eta_el = self.eta_el(100*(self.W_dot_wf/(self.params['eta_m']*self.eta_is)) / self.params['W_dot_el_rated'])/100
            
            self.W_dot_el = self.W_dot_wf/(self.params['eta_m']*self.eta_el)

            self.W_pp.set_W_dot(self.W_dot_wf)

            if self.V_dot_flag:
                print("Flowrate outside possible range. Impossible operation.")
            
        self.defined = True
