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

class Turb_polyn_eff(BaseComponent):
    def __init__(self):

        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector() 
        self.W_turb = WorkConnector()

    def get_required_inputs(self):
        
        if self.inputs == {}:
            if self.su.T is not None:
                self.inputs['su_T'] = self.su.T
            elif self.su.h is not None:
                self.inputs['su_h'] = self.su.h
            if self.su.p is not None:
                self.inputs['su_p'] = self.su.p
            if self.ex.p is not None:
                self.inputs['ex_p'] = self.ex.p
            if self.work_exp.speed is not None:
                self.inputs['N_rot'] = self.work_exp.speed
            if self.su.fluid is not None:
                self.inputs['su_fluid'] = self.su.fluid
        
        if self.inputs != {}:
            self.su.set_fluid(self.inputs['su_fluid'])
            if 'su_T' in self.inputs:
                self.su.set_T(self.inputs['su_T'])
            elif 'su_h' in self.inputs:
                self.su.set_h(self.inputs['su_h'])
            if 'su_p' in self.inputs:
                self.su.set_p(self.inputs['su_p'])
            if 'ex_p' in self.inputs:
                self.ex.set_p(self.inputs['ex_p'])
            if 'N_rot' in self.inputs:
                self.W_turb.set_N(self.inputs['N_rot'])
            if 'su_fluid' in self.inputs:
                self.su.set_fluid(self.inputs['su_fluid'])

        return ['su_p', 'su_T', 'ex_p', 'N_rot', 'su_fluid']

    def get_required_parameters(self):
        return [
            'D_inlet', 'N_turb_rated', 'turb_voltage', 'turb_phases', 'eta_max_motor', 
            'W_dot_el_rated', 'eta_m', 'eta_is_coefs', 'eta_is_coefs_red', 'A_th'
        ]
    
    def print_setup(self):
        print("=== Expander Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - W_dot: speed={self.work_exp.speed}")
        print(f"  - Q_dot_amb: temperature_in={self.heat_amb.temperature_in}")

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
    
    def eta_el_fun(self,P_ratio):

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
    
                
    def eta_is_turb(self):
        
        a = self.params['eta_is_coefs']
        p_cond = self.ex.p*1e-5
        m_dot = self.su.m_dot
        
        eta = (
        a[0]
        + a[1] * p_cond
        + a[2] * m_dot
        + a[3] * p_cond ** 2
        + a[4] * p_cond * m_dot
        + a[5] * m_dot ** 2
        + a[6] * p_cond ** 3
        + a[7] * p_cond ** 2 * m_dot
        + a[8] * p_cond * m_dot ** 2
        + a[9] * m_dot ** 3
        )
        
        return eta
    
    def solve(self):
        
        self.check_calculable()
        self.check_parametrized()
        
        if self.calculable and self.parametrized:
            
            "Speed of sound"
            R = 8.3144 # Ideal gas constant
            R_M = R/PropsSI('M',self.su.fluid)
            
            (cp_in,cv_in) = PropsSI(('C','CVMASS'),'T', self.su.T, 'P', self.su.p, self.su.fluid)
            gamma = cp_in/cv_in # Heat capacity ratio      
            
            self.a = np.sqrt(gamma*R_M*self.su.T) # m/s
            
            "Flowrate : Turbine considered as always choked"
            self.m_dot = self.a*self.su.D*self.params['A_th']
            self.su.set_m_dot(self.m_dot)
            self.su.set_m_dot(self.m_dot)
            
            "Computation of isentropic efficiency from map"
            
            self.eta_is = self.eta_is_turb()
            
            "Outlet state"
            
            h_ex_is = PropsSI('H', 'P', self.ex.p, 'S', self.su.s, self.su.fluid)
            h_ex = self.su.h - (self.eta_is/100)*(self.su.h - h_ex_is)
            
            self.ex.set_h(h_ex)
            self.ex.set_m_dot(self.su.m_dot)
            
            "Work"
            
            self.W_dot = (self.su.h - h_ex)*self.su.m_dot
            self.eta_el = self.eta_el_fun(100*self.W_dot/self.params['W_dot_el_rated'])/100
            self.W_el = self.W_dot*self.eta_el
            
            self.W_turb.set_W_dot(self.W_dot)

            self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")
        
    
    
    