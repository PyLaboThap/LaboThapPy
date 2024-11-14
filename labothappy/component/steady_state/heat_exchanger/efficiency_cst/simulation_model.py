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

class HXEffCst(BaseComponent):
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
        if self.su_C.T is not None:
            self.inputs['su_C_T'] = self.su_C.T
        if self.su_C.m_dot is not None:
            self.inputs['su_C_m_dot'] = self.su_C.m_dot
        if self.su_C.p is not None:
            self.inputs['su_C_p'] = self.su_C.p

        if self.su_H.fluid is not None:
            self.inputs['su_H_fluid'] = self.su_H.fluid
        if self.su_H.T is not None:
            self.inputs['su_H_T'] = self.su_H.T
        if self.su_H.h is not None:
            self.inputs['su_H_h'] = self.su_H.h
        if self.su_H.cp is not None:
            self.inputs['su_H_cp'] = self.su_H.cp
        if self.su_H.m_dot is not None:
            self.inputs['su_H_m_dot'] = self.su_H.m_dot
        if self.su_H.p is not None:
            self.inputs['su_H_p'] = self.su_H.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'su_C_fluid' in self.inputs:
            self.su_C.set_fluid(self.inputs['su_C_fluid'])
        if 'su_C_T' in self.inputs:
            self.su_C.set_T(self.inputs['su_C_T'])
        if 'su_C_h' in self.inputs:
            self.su_C.set_h(self.inputs['su_C_h'])
        if 'su_C_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['su_C_m_dot'])
        if 'su_C_p' in self.inputs:
            self.su_C.set_p(self.inputs['su_C_p'])

        if 'su_H_fluid' in self.inputs:
            self.su_H.set_fluid(self.inputs['su_H_fluid'])
        if 'su_H_T' in self.inputs:
            self.su_H.set_T(self.inputs['su_H_T'])
        if 'su_H_h' in self.inputs:
            self.su_H.set_h(self.inputs['su_H_h'])
        if 'su_H_cp' in self.inputs:
            self.su_H.set_cp(self.inputs['su_H_cp'])
        if 'su_H_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['su_H_m_dot'])
        if 'su_H_p' in self.inputs:
            self.su_H.set_p(self.inputs['su_H_p'])

        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_H_fluid', 'su_H_T', 'su_H_cp', 'su_H_m_dot']
    
    def get_required_parameters(self):
        return [
            'eta', # Efficiency
        ]
    
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

    def solve(self):
        # Ensure all required checks are performed

        if self.su_H.m_dot is None: 
            self.su_H.m_dot = self.su_C.m_dot
        
        if self.su_C.m_dot is None:
            self.su_C.m_dot = self.su_H.m_dot            

        self.check_calculable()
        self.check_parametrized()

        if not self.calculable:
            print("HTX IS NOT CALCULABLE")
            return

        if not self.parametrized:
            print("HTX IS NOT PARAMETRIZED")
            return

        "Heat Capacity Fluxes"

        C_cp = PropsSI('C','P',self.su_C.p,'H',self.su_C.h,self.su_C.fluid)
        H_cp = PropsSI('C','P',self.su_H.p,'H',self.su_H.h,self.su_H.fluid)

        C_Cdot = C_cp*self.su_C.m_dot
        H_Cdot = H_cp*self.su_H.m_dot
        
        Cdot_min = min(C_Cdot,H_Cdot)

        "Heat Transfer Rate"
        Q_dot = Cdot_min*self.params['eta']*(self.su_H.T - self.su_C.T)

        "Outlet state"

        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_p(self.su_C.p)
        self.ex_C.set_m_dot(self.su_C.m_dot)
        self.ex_C.set_h(self.su_C.h + Q_dot/self.ex_C.m_dot)

        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_p(self.su_H.p)
        self.ex_H.set_m_dot(self.su_H.m_dot)
        self.ex_H.set_h(self.su_H.h - Q_dot/self.ex_H.m_dot)
        
        self.Q_dot.set_Q_dot(Q_dot)

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





