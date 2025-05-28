# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:32:39 2024

@author: Basile

source : 
Semi-empirical correlation to model heat losses 
along solar parabolic trough collectors 
Rémi Dickes, Vincent Lemort and Sylvain Quoilin 

"""

import __init__
#import component.solar.parabolictroughcollector as parabolictroughcollector
from CoolProp.CoolProp import PropsSI
import numpy as np
from geometries.solar.parabolictrough_geometry import PT_Collector_Geom
from component.base_component import BaseComponent

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

class PT_collector(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.Q_amb = HeatConnector()
        self.v_wind = None
        self.Theta = None
        self.DNI = None

    def get_required_inputs(self):
            self.sync_inputs()
            # Return a list of required inputs
            # DNI : Direct Natural Irradiation
            # Incidence angle
            return ['su_p', 'su_m_dot', 'su_T', 'T_amb', 'DNI', 'Theta', 'v_wind', 'su_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.m_dot is not None:
            self.inputs['su_m_dot'] = self.su.m_dot
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.v_wind is not None:
            self.inputs['v_wind'] = self.v_wind
        if self.Q_amb.T_cold is not None:
            self.inputs['T_amb'] = self.Q_amb.T_cold
        if self.Theta is not None:
            self.inputs['Theta'] = self.Theta
        if self.DNI is not None:
            self.inputs['DNI'] = self.DNI


    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_m_dot' in self.inputs:
            self.su.set_m_dot(self.inputs['su_m_dot'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'v_wind' in self.inputs:
            self.v_wind = self.inputs['v_wind']
        if 'T_amb' in self.inputs:
            self.Q_amb.set_T_cold(self.inputs['T_amb'])
        if 'Theta' in self.inputs:
            self.Theta = self.inputs['Theta']
        if 'DNI' in self.inputs:
            self.DNI = self.inputs['DNI']

    def get_required_parameters(self):
        return [
                'coll_eff', 'L', 'W', 'A', 'A_r', 'm', 'L_f',

                'Tube_OD', 'Tube_V',

                'alpha_r', 'refl_m', 'epsilon_r', 'envel_tau', 'eta_opt',
                    
                'V_w_max_tr', 'V_w_max_st', 'Vdot_min', 'Vdot_max', 'T_f_min', 'T_f_max',

                'a', 'n_disc'
        ]
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - H_su: fluid={self.su['H'].fluid}, T={self.su['H'].T}, p={self.su['H'].p}, m_dot={self.su['H'].m_dot}")
        print(f"  - C_su: fluid={self.su['C'].fluid}, T={self.su['C'].T}, p={self.su['C'].p}, m_dot={self.su['C'].m_dot}")

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
    def heat_losses(self, k):
        "Calibrated for soponova_microcsp collector"        

        T_amb = self.Q_amb.T_cold - 273.15 # °C
        T_htf = self.T[k] - 273.15 # °C

        Terms = []
        Terms.append(self.params['a'][0]) 
        Terms.append(self.params['a'][1]*(T_htf - T_amb))
        Terms.append(self.params['a'][2]*(T_htf - T_amb)**2) 
        Terms.append(self.DNI*np.cos(self.Theta)*(self.params['a'][3]*T_htf**2))
        Terms.append(self.DNI*np.cos(self.Theta)*(self.params['a'][4]*np.sqrt(self.v_wind)))
        Terms.append(self.params['a'][5]*T_htf**3)
        Terms.append(self.v_wind*self.params['a'][6])
        Terms.append(self.v_wind*self.params['a'][7]*(T_htf - T_amb))
        Terms.append(np.sqrt(self.v_wind)*self.params['a'][8])
        Terms.append(np.sqrt(self.v_wind)*self.params['a'][9]*(T_htf - T_amb))

        Terms_np = np.array(Terms) # W/m (of line collector)

        return sum(Terms_np)

    def Q_dot_abs(self, k):

            # Sun absorbed power
            eta_opt = self.params['eta_opt']*self.params['alpha_r']*self.params['refl_m']*self.params['envel_tau']
            Q_dot_sun_raw = self.DNI*np.cos(self.Theta)*self.params['W']# *eta_opt # self.L_disc
            Q_dot_sun = Q_dot_sun_raw*eta_opt
            Opt_losses = Q_dot_sun_raw*(1-eta_opt)

            # Heat loss power
            HL = self.heat_losses(k)
            Q_dot_loss = HL #*self.L_disc
            Q_dot_abs = Q_dot_sun - Q_dot_loss

            return Q_dot_abs


    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:
            
            self.T = np.zeros(self.params['n_disc']+1)
            self.p = np.zeros(self.params['n_disc']+1)
            self.h = np.zeros(self.params['n_disc']+1)
            self.Q_dot_disc = np.zeros(self.params['n_disc'])
            self.Q_dot_abs_disc = np.zeros(self.params['n_disc'])
            self.eta_coll = np.zeros(self.params['n_disc'])

            self.T[0] = self.su.T
            self.p[0] = self.su.p
            self.h[0] = self.su.h

            self.L_disc = self.params['L']/self.params['n_disc']

            for i in range(len(self.Q_dot_disc)):
                self.Q_dot_disc[i] = self.DNI*np.cos(self.Theta)*self.L_disc*self.params['W']
                self.Q_dot_abs_disc[i] = self.Q_dot_disc[i]*self.params['coll_eff'](self.DNI, (self.T[i] - self.Q_amb.T_cold))
                self.eta_coll[i] = self.params['coll_eff'](self.DNI, (self.T[i] - self.Q_amb.T_cold))

                self.h[i+1] = self.h[i] + self.Q_dot_abs_disc[i]/self.su.m_dot
                self.p[i+1] = self.p[i]
                self.T[i+1] = PropsSI('T','H',self.h[i+1],'P',self.p[i+1],self.su.fluid)

            self.Q_dot = sum(self.Q_dot_abs_disc)

            self.ex.set_fluid(self.su.fluid)
            self.ex.set_m_dot(self.su.m_dot)
            self.ex.set_h(self.h[-1])
            self.ex.set_p(self.p[-1])

            return
        
        else:
            if not self.calculable:
                print("Input of the component not completely known. Required inputs:")
                for input in self.get_required_inputs():
                    if input not in self.inputs:
                        print(f"  - {input}")
                        
            if not self.parametrized:
                print("Parameters of the component not completely known. Required parameters:")
                for param in self.get_required_parameters():
                    if param not in self.params:
                        print(f"  - {param}")
    
    def print_results(self):
        if self.defined:
            print("=== Heat Exchanger Results ===")
            print(f"  - H_ex: fluid={self.ex['H'].fluid}, T={self.ex['H'].T}, p={self.ex['H'].p}, m_dot={self.ex['H'].m_dot}")
            print(f"  - C_ex: fluid={self.ex['C'].fluid}, T={self.ex['C'].T}, p={self.ex['C'].p}, m_dot={self.ex['C'].m_dot}")
            print(f"  - Q_dot: temperature_in={self.Q_dot}")

        else:
            print("Heat Exchanger component is not defined. Ensure it is solved first.")

#-------------------------------------------------------------