#%%

import __init__

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
    """
    Component: Counterflow Heat Exchanger with Constant Effectiveness (HXEffCst)
    
    Model: Simplified Heat Exchanger Model with Constant Effectiveness
    
    **Description**:
    
        This model simulates a counterflow heat exchanger with a fixed effectiveness (η). It estimates the maximum possible heat transfer based on inlet conditions and applies the given effectiveness to compute the actual heat transfer rate. The model assumes a simplified configuration without discretization or detailed internal states. It is suitable for quick, steady-state evaluations in system-level models.
    
    **Assumptions**:
    
        - Steady-state operation.
        - Constant effectiveness (η) is applied across the entire heat exchanger.
        - No phase change or pressure drop effects are modeled.
        - Uniform flow distribution on both hot and cold sides.
        - Fluid properties are retrieved from CoolProp.
    
    **Connectors**:
    
        su_C (MassConnector): Mass connector for the cold-side supply.
        
        su_H (MassConnector): Mass connector for the hot-side supply.
        
        ex_C (MassConnector): Mass connector for the cold-side exhaust.
        
        ex_H (MassConnector): Mass connector for the hot-side exhaust.
        
        Q_dot (HeatConnector): Heat transfer connector representing the total exchanged heat.
    
    **Parameters**:
    
        eta (float): Effectiveness of the heat exchanger [-].
    
    **Inputs**:
    
        su_H_fluid (str): Hot-side fluid.
        
        su_H_h (float): Hot-side inlet specific enthalpy [J/kg].
        
        su_H_p (float): Hot-side inlet pressure [Pa].
        
        su_H_m_dot (float): Hot-side mass flow rate [kg/s].
    
        su_C_fluid (str): Cold-side fluid.
        
        su_C_h (float): Cold-side inlet specific enthalpy [J/kg].
        
        su_C_p (float): Cold-side inlet pressure [Pa].
        
        su_C_m_dot (float): Cold-side mass flow rate [kg/s].
    
    **Outputs**:
    
        ex_C_h: Cold-Side Exhaust specific enthalpy at outlet [J/kg].
        
        ex_C_p: Cold-Side Exhaust pressure at outlet [Pa].
                    
        ex_C_h: Hot-Side specific enthalpy at outlet [J/kg].
        
        ex_C_p: Hot-Side pressure at outlet [Pa].
        
        Q_dot: Total heat transfer rate across the exchanger [W].

    """

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
        return['Csu_fluid', 'Csu_h', 'Csu_p', 'Csu_m_dot', 'Hsu_fluid', 'Hsu_h', 'Hsu_p', 'Hsu_m_dot']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su_C.fluid is not None:
            self.inputs['Csu_fluid'] = self.su_C.fluid
        if self.su_C.h is not None:
            self.inputs['Csu_h'] = self.su_C.h
        if self.su_C.T is not None:
            self.inputs['Csu_T'] = self.su_C.T
        if self.su_C.m_dot is not None:
            self.inputs['Csu_m_dot'] = self.su_C.m_dot
        if self.su_C.p is not None:
            self.inputs['Csu_p'] = self.su_C.p

        if self.su_H.fluid is not None:
            self.inputs['Hsu_fluid'] = self.su_H.fluid
        if self.su_H.T is not None:
            self.inputs['Hsu_T'] = self.su_H.T
        if self.su_H.h is not None:
            self.inputs['Hsu_h'] = self.su_H.h
        if self.su_H.cp is not None:
            self.inputs['Hsu_cp'] = self.su_H.cp
        if self.su_H.m_dot is not None:
            self.inputs['Hsu_m_dot'] = self.su_H.m_dot
        if self.su_H.p is not None:
            self.inputs['Hsu_p'] = self.su_H.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'Csu_fluid' in self.inputs:
            self.su_C.set_fluid(self.inputs['Csu_fluid'])
        if 'Csu_T' in self.inputs:
            self.su_C.set_T(self.inputs['Csu_T'])
        if 'Csu_h' in self.inputs:
            self.su_C.set_h(self.inputs['Csu_h'])
        if 'Csu_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['Csu_m_dot'])
        if 'Csu_p' in self.inputs:
            self.su_C.set_p(self.inputs['Csu_p'])

        if 'Hsu_fluid' in self.inputs:
            self.su_H.set_fluid(self.inputs['Hsu_fluid'])
        if 'Hsu_T' in self.inputs:
            self.su_H.set_T(self.inputs['Hsu_T'])
        if 'Hsu_h' in self.inputs:
            self.su_H.set_h(self.inputs['Hsu_h'])
        if 'Hsu_cp' in self.inputs:
            self.su_H.set_cp(self.inputs['Hsu_cp'])
        if 'Hsu_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['Hsu_m_dot'])
        if 'Hsu_p' in self.inputs:
            self.su_H.set_p(self.inputs['Hsu_p'])

        return['Csu_fluid', 'Csu_h', 'Csu_p', 'Csu_m_dot', 'Hsu_fluid', 'Hsu_h', 'Hsu_p', 'Hsu_m_dot']
    
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

        "Define Q_dot_max through enthalpies"
        
        H_h_id = PropsSI('H', 'P', self.su_H.p, 'T', self.su_C.T, self.su_H.fluid)
        H_c_id = PropsSI('H', 'P', self.su_C.p, 'T', self.su_H.T, self.su_C.fluid)
        
        Q_dot_maxh = self.su_H.m_dot*abs(H_h_id-self.su_H.h)
        Q_dot_maxc = self.su_C.m_dot*abs(H_c_id-self.su_C.h)
        
        Q_dot_max = min(Q_dot_maxh,Q_dot_maxc)
        
        "Heat Transfer Rate"
        Q_dot = self.params['eta']*Q_dot_max

        "Outlet states"   
        self.update_connectors(Q_dot)
        self.solved = True
        return

    def update_connectors(self, Q_dot):
        
        "Mass Connectors"
        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)
        self.ex_C.set_h(self.su_C.h + Q_dot/self.su_C.m_dot)
        self.ex_C.set_p(self.su_C.p)

        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_m_dot(self.su_H.m_dot)
        self.ex_H.set_h(self.su_H.h - Q_dot/self.su_H.m_dot)
        self.ex_H.set_p(self.su_H.p)

        "Heat conector"
        self.Q_dot.set_Q_dot(Q_dot)

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





