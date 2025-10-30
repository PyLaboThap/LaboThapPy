# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:32:39 2024

@author: Basile

source : 
Semi-empirical correlation to model heat losses 
along solar parabolic trough collectors 
Rémi Dickes, Vincent Lemort and Sylvain Quoilin 
https://hdl.handle.net/2268/182680

"""

import __init__
#import component.solar.parabolictroughcollector as parabolictroughcollector
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import numpy as np
from component.base_component import BaseComponent

from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector


class PT_collector(BaseComponent):
    """
    Component: Parabolic Trough Collector
    
    Model: Semi empirical
    
    **Description**
    
    This model determines the thermal power absorbed by a parabolic trough collector and its heat transfer fluid exhaust specific enthalpy .

    
    **Assumptions**:
        - Heat losses to the ambiant are considered.

    **Connectors**:
    
        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

        Q_amb (HeatConnector): Heat connector to the ambiant.
        
    **Parameters**:
        
        coll_eff: Collector efficieny [-]
        
        a [1..9]: semi empirical parameters 
            
        n_disc: number of discretisations in the lenght of the heat collection element. [-]
        
        envel_tau: transmittance of the heat collection element envelop. [-]
        
        alpha_r: Receiver absorptivity. [-]
        
        refl_m: Mirror reflectivity. [-]
        
        epsilon_r: Receiver emittance @ 400°C. [-]
        
        eta_opt: Overall Optical Efficiency. [-]

        V_w_max_tr: maximum wind speed (stowed). [m/s]
        
        V_w_max_st: minimum wind speed (stowed). [m/s]
        
        Vdot_min: minimum recommended flowrate. [m^3/s]
        
        Vdot_max: maximum recommended flowrate. [m^3/s]
        
        T_f_min: Minimum operating fluid temperature. [K]
        
        T_f_max: Maximum operating fluid temperature. [K]
        
        Geometrical parameters:
            
            L: length of the parabolic trough collector. [m]
            
            W: Width of the parabolic trough collector. [m]
            
            A: Area of the parabolic trough collector. [m^2]
            
            A_r: Reflective area. [m^2]
            
            m: Collector weight. [kg]
            
            L_f: Focal Length. [m]
            
            Tube_OD: Tube external diameter. [m]
            
            Tube_V: Tube volume. [m^3]
            
    **Inputs**:
        
        fluid: Fluid at the suction of the PT collector. [-]
        
        m_dot: Mass flow rate in the PT collector (at suction). [kg/s]
        
        T_su: Temperature at the suction. [K]
        
        p_su: Pressure at the suction. [Pa]
        
        v_wind: Wind speed on the PT collector's heat collection element. [m/s]
        
        T_amb: Ambiant temperature. [K]
        
        Theta: Incidence angle. [rad]
        
        DNI: Direct Normal Irradiance. [W/m^2]
        
    **Ouputs**:

        h_ex: Exhaust side specific enthalpy. [J/kg] 

        Q_dot: Thermal Power absorbed by the PT collector. [W]

    """
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
            return ['P_su', 'm_dot', 'T_su', 'T_amb', 'DNI', 'Theta', 'v_wind', 'fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['fluid'] = self.su.fluid
        if self.su.m_dot is not None:
            self.inputs['m_dot'] = self.su.m_dot
        if self.su.T is not None:
            self.inputs['T_su'] = self.su.T
        if self.su.p is not None:
            self.inputs['P_su'] = self.su.p
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
        if 'fluid' in self.inputs:
            self.su.set_fluid(self.inputs['fluid'])
        if 'm_dot' in self.inputs:
            self.su.set_m_dot(self.inputs['m_dot'])
        if 'T_su' in self.inputs:
            self.su.set_T(self.inputs['T_su'])
        if 'P_su' in self.inputs:
            self.su.set_p(self.inputs['P_su'])
        if 'v_wind' in self.inputs:
            self.v_wind = self.inputs['v_wind']
        if 'T_amb' in self.inputs:
            self.Q_amb.set_T_cold(self.inputs['T_amb'])
        if 'Theta' in self.inputs:
            self.Theta = self.inputs['Theta']
        if 'DNI' in self.inputs:
            self.DNI = self.inputs['DNI']



    # def set_inputs(self, **kwargs):
    #     """Set inputs directly through a dictionary and update connector properties."""
    #     self.inputs.update(kwargs)

    #     # Define mappings from input keys to methods
    #     input_methods = {
    #         # su connector inputs
    #         'fluid':     lambda val: self.su.set_fluid(val),
    #         'x_su':      lambda val: self.su.set_x(val),
    #         'T_su':      lambda val: self.su.set_T(val),
    #         'h_su':      lambda val: self.su.set_h(val),
    #         'P_su':      lambda val: self.su.set_p(val),
    #         'm_dot':     lambda val: self.su.set_m_dot(val),

    #         # su_1 connector inputs
    #         'fluid_su_1': lambda val: self.su_1.set_fluid(val),
    #         'T_su_1':    lambda val: self.su_1.set_T(val),
    #         'h_su_1':    lambda val: self.su_1.set_h(val),
    #         'P_su_1':    lambda val: self.su_1.set_p(val),
    #         'm_dot_su_1': lambda val: self.su_1.set_m_dot(val),

    #         # su_2 connector inputs
    #         'fluid_su_2': lambda val: self.su_2.set_fluid(val),
    #         'T_su_2':    lambda val: self.su_2.set_T(val),
    #         'h_su_2':    lambda val: self.su_2.set_h(val),
    #         'P_su_2':    lambda val: self.su_2.set_p(val),
    #         'm_dot_su_2': lambda val: self.su_2.set_m_dot(val),

    #         # su_H connector inputs
    #         'fluid_H': lambda val: self.su_H.set_fluid(val),
    #         'T_su_H':    lambda val: self.su_H.set_T(val),
    #         'h_su_H':    lambda val: self.su_H.set_h(val),
    #         'P_su_H':    lambda val: self.su_H.set_p(val),
    #         'm_dot_H':   lambda val: self.su_H.set_m_dot(val),

    #         # su_C connector inputs
    #         'fluid_C': lambda val: self.su_C.set_fluid(val),
    #         'T_su_C':    lambda val: self.su_C.set_T(val),
    #         'h_su_C':    lambda val: self.su_C.set_h(val),
    #         'P_su_C':    lambda val: self.su_C.set_p(val),
    #         'm_dot_C':   lambda val: self.su_C.set_m_dot(val),

    #         # ex connector inputs
    #         'P_ex':      lambda val: self.ex.set_p(val),
    #         'T_ex':      lambda val: self.ex.set_T(val),
    #         'h_ex':      lambda val: self.ex.set_h(val),
    #         'x_ex':      lambda val: self.ex.set_x(val),

    #         # ex_1 connector inputs
    #         'P_ex_1':    lambda val: self.ex_1.set_p(val),
    #         'T_ex_1':    lambda val: self.ex_1.set_T(val),
    #         'h_ex_1':    lambda val: self.ex_1.set_h(val),

    #         # ex_2 connector inputs
    #         'P_ex_2':    lambda val: self.ex_2.set_p(val),
    #         'T_ex_2':    lambda val: self.ex_2.set_T(val),
    #         'h_ex_2':    lambda val: self.ex_2.set_h(val),         

    #         # ex_H connector inputs
    #         'P_ex_H':    lambda val: self.ex_H.set_p(val),
    #         'T_ex_H':    lambda val: self.ex_H.set_T(val),
    #         'h_ex_H':    lambda val: self.ex_H.set_h(val),

    #         # ex_C connector inputs
    #         'P_ex_C':    lambda val: self.ex_C.set_p(val),
    #         'T_ex_C':    lambda val: self.ex_C.set_T(val),
    #         'h_ex_C':    lambda val: self.ex_C.set_h(val),

    #         # W connector inputs
    #         'N_rot':     lambda val: self.W.set_N(val),

    #         # Q_amb connector inputs
    #         'T_amb':     lambda val: self.Q_amb.set_T_cold(val),
    #         'DNI':       lambda val: self.Q_amb.set_DNI(val),
    #         'Theta':     lambda val: self.Q_amb.set_Theta(val), "Irradiance angle"
    #         'v_wind':    lambda val: self.solar.set_v_wind(val),

    #     }

    #     unknown_keys = []  # To collect any keys that do not match the input methods

    #     for key, value in self.inputs.items():
    #         method = input_methods.get(key)
    #         if method:
    #             try:
    #                 method(value)
    #             except Exception as e:
    #                 # Optionally log the exception or raise with more context
    #                 pass  # Replace with logging if desired
    #         else:
    #             unknown_keys.append(key)

    #     if unknown_keys:
    #         raise ValueError(f"Unrecognized input keys: {', '.join(unknown_keys)}")
    #     return


    # def sync_inputs(self):
    #     """Synchronize the inputs dictionary with the connector states."""

    #     # Lazy getters: only access if the connector exists
    #     attribute_map = {
    #         # su connectors
    #         'fluid':     lambda: self.su.fluid if hasattr(self,'su') else None,
    #         'T_su':      lambda: self.su.T if hasattr(self,'su') else None,
    #         'h_su':      lambda: self.su.h if hasattr(self,'su') else None,
    #         'P_su':      lambda: self.su.p if hasattr(self,'su') else None,
    #         'm_dot':     lambda: self.su.m_dot if hasattr(self, 'su') else None,

    #         # su_H connector
    #         'fluid_H':   lambda: self.su_H.fluid if hasattr(self,'su_H') else None,
    #         'T_su_H':    lambda: self.su_H.T if hasattr(self,'su_H') else None,
    #         'h_su_H':    lambda: self.su_H.h if hasattr(self,'su_H') else None,
    #         'P_su_H':    lambda: self.su_H.p if hasattr(self,'su_H') else None,
    #         'm_dot_H':   lambda: self.su_H.m_dot if hasattr(self,'su_H') else None,

    #         # su_C connector
    #         'fluid_C':   lambda: self.su_C.fluid if hasattr(self,'su_C') else None,
    #         'T_su_C':    lambda: self.su_C.T if hasattr(self,'su_C') else None,
    #         'h_su_C':    lambda: self.su_C.h if hasattr(self,'su_C') else None,
    #         'P_su_C':    lambda: self.su_C.p if hasattr(self,'su_C') else None,
    #         'm_dot_C':   lambda: self.su_C.m_dot if hasattr(self,'su_C') else None,

    #         # ex connector
    #         'P_ex':      lambda: self.ex.p if hasattr(self,'ex') else None,
    #         'T_ex':      lambda: self.ex.T if hasattr(self,'ex') else None,
    #         'h_ex':      lambda: self.ex.h if hasattr(self,'ex') else None,

    #         # ex_C connector
    #         'P_ex_C':    lambda: self.ex_C.p if hasattr(self,'ex_C') else None,
    #         'T_ex_C':    lambda: self.ex_C.T if hasattr(self,'ex_C') else None,
    #         'h_ex_C':    lambda: self.ex_C.h if hasattr(self,'ex_C') else None,

    #         # ex_H connector
    #         'P_ex_H':    lambda: self.ex_H.p if hasattr(self,'ex_H') else None,
    #         'T_ex_H':    lambda: self.ex_H.T if hasattr(self,'ex_H') else None,
    #         'h_ex_H':    lambda: self.ex_H.h if hasattr(self,'ex_H') else None,

    #         # W connector
    #         'N_rot':     lambda: self.W.N if hasattr(self,'W') else None,

    #         # Q_amb connector
    #         'T_amb':     lambda: self.Q_amb.T_cold if hasattr(self,'Q_amb') else None,
    #         'DNI':       lambda: self.Q_amb.DNI if hasattr(self,'Q_amb') else None,
    #         'Theta':       lambda: self.Q_amb.Theta if hasattr(self,'Q_amb') else None,
    #         'v_wind':       lambda: self.solar.v_wind if hasattr(self,'Q_amb') else None,
    #     }

    #     self.inputs = getattr(self,'inputs',{})

    #     for key, getter in attribute_map.items():
    #         try:
    #             value = getter()
    #             if value is not None:
    #                 self.inputs[key] = value
    #         except Exception:
    #             pass  # Optional: add logging for debugging

    #     return
    
    
    
    
    
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

        self.AS = CP.AbstractState('HEOS', self.su.fluid)
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
                self.AS.update(CP.HmassP_INPUTS, self.h[i+1], self.p[i+1])               
                self.T[i+1] = self.AS.T()

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