import sys
import os

# Get the absolute path of the directory that contains the script (simulation_model.py)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Determine the project root directory (which contains both 'connector' and 'component')
project_root = os.path.abspath(os.path.join(current_dir, '..', '..', '..'))

# Add the project root to sys.path if it's not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)

#%%

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

# from component.heat_exchanger.moving_boundary.simple_model.modules.U import U_Gnielinski_calibrated, U_DittusBoelter, U_Cooper_calibrater, U_Thonon

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np
import math

class HXPinchCst(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su_C = MassConnector() # Working fluid supply
        self.su_H = MassConnector() # Secondary fluid supply
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.Q_dot = HeatConnector()
        self.guesses = {}

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_H_fluid', 'su_H_T', 'su_H_cp', 'su_H_m_dot']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su_C.fluid is not None:
            self.inputs['su_C_fluid'] = self.su_C.fluid
        if self.su_C.h is not None:
            self.inputs['su_C_h'] = self.su_C.h
        if self.su_C.m_dot is not None:
            self.inputs['su_C_m_dot'] = self.su_C.m_dot
        if self.su_H.fluid is not None:
            self.inputs['su_H_fluid'] = self.su_H.fluid
        if self.su_H.T is not None:
            self.inputs['su_H_T'] = self.su_H.T
        if self.su_H.cp is not None:
            self.inputs['su_H_cp'] = self.su_H.cp
        if self.su_H.m_dot is not None:
            self.inputs['su_H_m_dot'] = self.su_H.m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'su_C_fluid' in self.inputs:
            self.su_C.set_fluid(self.inputs['su_C_fluid'])
        if 'su_C_h' in self.inputs:
            self.su_C.set_h(self.inputs['su_C_h'])
        if 'su_C_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['su_C_m_dot'])
        if 'su_H_fluid' in self.inputs:
            self.su_H.set_fluid(self.inputs['su_H_fluid'])
        if 'su_H_T' in self.inputs:
            self.su_H.set_T(self.inputs['su_H_T'])
        if 'su_H_cp' in self.inputs:
            self.su_H.set_cp(self.inputs['su_H_cp'])
        if 'su_H_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['su_H_m_dot'])

        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_H_fluid', 'su_H_T', 'su_H_cp', 'su_H_m_dot']
    
    def get_required_parameters(self):
        return [
            'Pinch', # pinch point
            'Delta_T_sh_sc', # Superheating or subcooling
            'type_HX' # Evaporator or condenser
        ]
    
    def get_required_guesses(self):
        return [ 'P_sat' ]
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")

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

    def system_evap(self, x):
        P_ev = x[0]
        # Ensure the pressure is non-negative
        if P_ev < 0:
            P_ev = 10000
        
        # Get the temperature of the evaporator based on the pressure and quality
        T_ev = PropsSI('T', 'P', P_ev, 'Q', 0.5, self.su_C.fluid)
        
        # Refrigerant side calculations
        # Liquid zone
        h_ev_su = self.su_C.h
        h_ev_l = PropsSI('H', 'P', P_ev, 'Q', 0, self.su_C.fluid)
        Q_dot_ev_l = self.su_C.m_dot * (h_ev_l - h_ev_su)

        if Q_dot_ev_l < 0:
            Q_dot_ev_l = 0
        
        # Two-phase zone
        h_ev_v = PropsSI('H', 'P', P_ev, 'Q', 1, self.su_C.fluid)
        Q_dot_ev_tp = self.su_C.m_dot * (h_ev_v - h_ev_l)
        
        # Vapor zone
        self.T_C_ex = T_ev + self.params['Delta_T_sh_sc']
        h_ev_ex = PropsSI('H', 'P', P_ev, 'T', self.T_C_ex, self.su_C.fluid)
        Q_dot_ev_v = self.su_C.m_dot * (h_ev_ex - h_ev_v)
        
        # Total heat transfer
        Q_dot_ev = Q_dot_ev_l + Q_dot_ev_tp + Q_dot_ev_v
        
        # Secondary fluid side calculations
        self.T_H_ex = self.su_H.T - Q_dot_ev / (self.su_H.m_dot * self.su_H.cp)
        self.T_H_v = self.su_H.T - Q_dot_ev_v / (self.su_H.m_dot * self.su_H.cp)
        self.T_H_l = self.T_H_v - Q_dot_ev_tp / (self.su_H.m_dot * self.su_H.cp)
        
        # Calculate pinch point and residual
        PP = min(self.su_H.T - self.T_C_ex, self.T_H_l - T_ev)
        res = abs(PP - self.params['Pinch']) / self.params['Pinch']
        
        # Update the state of the working fluid
        self.Q = Q_dot_ev
        self.P_sat = P_ev
        
        return res
    
    def system_cond(self, x):
        P_cd = x[0]
        
        # Ensure the pressure is non-negative
        if P_cd < 0:
            P_cd = 10000
        
        # Get the temperature of the condenser based on pressure and quality
        T_cd = PropsSI('T', 'P', P_cd, 'Q', 0.5, self.su_H.fluid)
        
        # Refrigerant side calculations
        # Initial state enthalpy
        h_wf_su_cd = self.su_H.h
        
        # Vapor zone
        h_wf_cd_v = PropsSI('H', 'P', P_cd, 'Q', 1, self.su_H.fluid)
        Q_dot_cd_v = self.su_H.m_dot * (h_wf_su_cd - h_wf_cd_v)

        if Q_dot_cd_v < 0:
            Q_dot_cd_v = 0
        
        # Two-phase zone
        h_wf_cd_l = PropsSI('H', 'P', P_cd, 'Q', 0, self.su_H.fluid)
        Q_dot_cd_tp = self.su_H.m_dot * (h_wf_cd_v - h_wf_cd_l)
        
        # Liquid zone
        self.T_H_ex = T_cd - self.params['Delta_T_sh_sc']
        h_H_ex_cd = PropsSI('H', 'P', P_cd, 'T', self.T_H_ex, self.su_H.fluid)
        Q_dot_cd_l = self.su_H.m_dot * (h_wf_cd_l - h_H_ex_cd)
        
        # Total heat transfer
        Q_dot_cd = Q_dot_cd_v + Q_dot_cd_tp + Q_dot_cd_l
        
        # Water side calculations
        self.T_C_ex = self.su_C.T + Q_dot_cd / (self.su_C.m_dot * self.su_C.cp)
        T_C_cd_l = self.su_C.T + Q_dot_cd_l / (self.su_C.m_dot * self.su_C.cp)
        T_C_cd_v = T_C_cd_l + Q_dot_cd_tp / (self.su_C.m_dot * self.su_C.cp)
        
        # Pinch point position
        Pinch_cd = min(self.T_H_ex - self.su_C.T, T_cd - T_C_cd_v)

        # Calculate residual
        res = abs(Pinch_cd - self.params['Pinch']) / self.params['Pinch']
        
        # Update the state of the working fluid
        self.Q = Q_dot_cd
        self.P_sat = P_cd
        
        return res

    def solve(self):
        # Ensure all required checks are performed

        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            print("HTX IS NOT CALCULABLE")
            return

        # Determine the type of heat exchanger and set the initial guess for pressure
        if self.params['type_HX'] == 'evaporator':
            P_ev_guess = self.guesses.get('P_sat', PropsSI('P', 'T', 120 + 273.15, 'Q', 0.5, self.su_C.fluid)) # Guess the saturation pressure, first checks if P_sat is in the guesses dictionary, if not it calculates it
            x = [P_ev_guess]

            try:
                """EVAPORATOR MODEL"""
                fsolve(self.system_evap, x)

                """Update connectors after the calculations"""
                self.update_connectors()

                # Mark the model as solved if successful
                self.solved = True
            except Exception as e:
                # Handle any errors that occur during solving
                self.solved = False
                print(f"Convergence problem in evaporator model: {e}")

        elif self.params['type_HX'] == 'condenser':
            P_cd_guess = self.guesses.get('P_sat', PropsSI('P', 'T', 25 + 273.15, 'Q', 0.5, self.su_H.fluid)) # Guess the saturation pressure, first checks if P_sat is in the guesses dictionary, if not it calculates it
            x = [P_cd_guess]

            try:
                """CONDENSER MODEL"""
                fsolve(self.system_cond, x)

                """Update connectors after the calculations"""
                self.update_connectors()

                # Mark the model as solved if successful
                self.solved = True
            except Exception as e:
                # Handle any errors that occur during solving
                self.solved = False
                print(f"Convergence problem in condenser model: {e}")


    def update_connectors(self):
        
        "Mass Connectors"

        if self.params['type_HX'] == 'evaporator':

            self.su_C.set_p(self.P_sat)

            self.ex_C.set_fluid(self.su_C.fluid)
            self.ex_C.set_T(self.T_C_ex)
            self.ex_C.set_p(self.P_sat)
            self.ex_C.set_m_dot(self.su_C.m_dot)

            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_m_dot(self.su_H.m_dot)
            self.ex_H.set_T(self.T_H_ex)
            
            "Heat conector"
            self.Q_dot.set_Q_dot(self.Q)

        else: 

            self.su_H.set_p(self.P_sat)

            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_T(self.T_H_ex)
            self.ex_H.set_p(self.P_sat)
            self.ex_H.set_m_dot(self.su_H.m_dot)

            self.ex_C.set_fluid(self.su_C.fluid)
            self.ex_C.set_m_dot(self.su_C.m_dot)
            self.ex_C.set_T(self.T_C_ex)
            
            "Heat conector"
            self.Q_dot.set_Q_dot(self.Q)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q_dot.Q_dot}")

        if self.params['type_HX'] == 'evaporator':
            print(f"P_sat: {self.su_C.p}")
        else:
            print(f"P_sat: {self.su_H.p}")

        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")
        print("======================")





