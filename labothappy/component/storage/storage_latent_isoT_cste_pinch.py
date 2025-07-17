import __init__

"""
import sys
import os

# Get the absolute path of the directory that contains the script (simulation_model.py)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Determine the project root directory (which contains both 'connector' and 'component')
project_root = os.path.abspath(os.path.join(current_dir, '..', '..')) 

# Add the project root to sys.path if it's not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)
"""

#%%

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent


from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, root, minimize
import numpy as np
import math

class StorageLatentIsothermalCstePinch(BaseComponent):
    """
    Component: Heat Exchanger with constant pinch point.

    **Description**:
    
        Simulates a Heat Exchanger with a constant pinch point.
        The Pinch Point is understood as being the location on the heat exchanger length where
        the temperature difference between the hot and the cold fluid is minimal.

    **Assumptions**:

        - Steady-state operation
        - No pressure drops considered
        - No loss to the ambient considered.

    **Connectors**:

        su_H (MassConnector): Mass connector for the hot suction side.
        su_C (MassConnector): Mass connector for the cold suction side.

        ex_H (MassConnector): Mass connector for the hot exhaust side.
        ex_C (MassConnector): Mass connector for the cold exhaust side.

        Q_dot (HeatConnector): Heat connector for the heat transfer between the fluids

    **Parameters**:

        Pinch: Pinch point temperature difference [K] or [°C]
        
        Delta_T_sh_sc: Superheating or subcooling, depending if the HEX is an evaporator (superheating) or a condenser (subcooling)
            
        type_HX: HX type, i.e. evaporator or condenser

    **Inputs**:
        
        su_C_fluid: Cold suction side fluid. [-]

        su_C_h: Cold suction side enthalpy. [J/kg]

        su_C_p: Cold suction side pressure. [Pa]

        su_C_m_dot: Cold suction side mass flow rate. [kg/s]

        su_H_fluid: Hot suction side fluid. [-]

        su_H_h: Hot suction side enthalpy. [J/kg]

        su_H_p: Hot suction side pressure. [Pa]

        su_H_m_dot: Hot suction side mass flow rate. [kg/s]

    **Outputs**:

        ex_C_h: Cold exhaust side enthalpy. [J/kg]

        ex_C_p: Cold exhaust side pressure. [Pa]

        ex_H_h: Hot exhaust side enthalpy. [J/kg]

        ex_H_p: Hot exhaust side pressure. [Pa]

        Q_dot: Heat Exchanger's heat duty. [W]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Working fluid supply
        self.ex = MassConnector()

        self.sto_fluid = MassConnector()

        self.Q_dot = HeatConnector()

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return ['fluid', 'sto_fluid', 'h_su', 'T_su', 'm_dot']
    
    def get_required_parameters(self):
        return [
            'Pinch', # pinch point
            'Delta_T_sh_sc', # Superheating or subcooling
        ]
    
    def get_required_guesses(self):
        return []
    
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

        self.check_calculable()
        self.check_parametrized()

        if not self.calculable:
            print("StorageLatentIsothermalCstePinch is not calculable, check inputs.")
            return
        
        if not self.parametrized:
            print("StorageLatentIsothermalCstePinch is not parametrized, check parameters.")
            return

        P_sto = 101325 # Pa : Latent storage at ambient pressure
        
        # Get temperature at the triple point (freezing/melting point)
        self.sto_fluid.set_T(PropsSI("T_triple", "P", P_sto, "Q", 0, "Water"))
        self.sto_fluid.set_p(P_sto)

        if self.su.T <= self.sto_fluid.T:
            if self.su.T <= self.sto_fluid.T - self.params['Pinch'] - self.params['Delta_T_sh_sc']:
                self.T_ex = self.sto_fluid.T - self.params['Pinch']
                self.T_sat = self.T_ex - self.params['Delta_T_sh_sc']
                
                self.P_sat = PropsSI('P', 'Q', 0, 'T', self.T_sat, self.su.fluid)
                self.h_ex = PropsSI('H', 'P', self.P_sat, 'T', self.T_ex, self.su.fluid)

                self.h_sat_v = PropsSI('H', 'Q', 1, 'T', self.T_sat, self.su.fluid)
                self.h_sat_l = PropsSI('H', 'Q', 0, 'T', self.T_sat, self.su.fluid)
                
                self.Q_dot_3 = self.su.m_dot*(self.h_ex - self.h_sat_v)
                self.Q_dot_2 = self.su.m_dot*(self.h_sat_v - self.h_sat_l)

                h_su = PropsSI('H', 'P', self.P_sat, 'T', self.su.T, self.su.fluid)

                self.Q_dot_1 = self.su.m_dot*(self.h_sat_l - h_su)
                
                self.Q = self.Q_dot_1 + self.Q_dot_2 + self.Q_dot_3
                
            else: # Pinch and SC_SH cannot be satisfied
                # self.Q = 0

                # self.T_ex = self.su.T
                # self.P_sat = self.su.p
                # self.h_ex = self.su.h

                return
        else:
            if self.su.T >= self.sto_fluid.T + self.params['Pinch'] + self.params['Delta_T_sh_sc']:
                self.T_ex = self.sto_fluid.T + self.params['Pinch']
                self.T_sat = self.T_ex + self.params['Delta_T_sh_sc']
                
                self.P_sat = PropsSI('P', 'Q', 0, 'T', self.T_sat, self.su.fluid)
                self.h_ex = PropsSI('H', 'P', self.P_sat, 'T', self.T_ex, self.su.fluid)

                self.h_sat_v = PropsSI('H', 'Q', 1, 'T', self.T_sat, self.su.fluid)
                self.h_sat_l = PropsSI('H', 'Q', 0, 'T', self.T_sat, self.su.fluid)
                
                self.Q_dot_3 = self.su.m_dot*(self.h_sat_l - self.h_ex)
                self.Q_dot_2 = self.su.m_dot*(self.h_sat_v - self.h_sat_l)

                h_su = PropsSI('H', 'P', self.P_sat, 'T', self.su.T, self.su.fluid)
                
                self.Q_dot_1 = self.su.m_dot*(h_su - self.h_sat_v)

                self.Q = self.Q_dot_1 + self.Q_dot_2 + self.Q_dot_3
            else: # Pinch and SC_SH cannot be satisfied
                raise ValueError("Pinch and SC_SH are not satisfied")
                return
        
        self.solved = True
        self.update_connectors()


    def update_connectors(self):
        
        "Mass Connectors"

        self.su.set_p(self.P_sat)

        self.ex.set_fluid(self.su.fluid)
        self.ex.set_T(self.T_ex)
        self.ex.set_p(self.P_sat)
        self.ex.set_m_dot(self.su.m_dot)
        
        "Heat conector"
        self.Q_dot.set_Q_dot(self.Q)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q_dot.Q_dot}")

        print(f"P_sat: {self.su.p}")

        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - sto: fluid={self.sto_fluid.fluid}, T={self.sto_fluid.T}, p={self.sto_fluid.p}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")
        print("======================")

    def plot_disc(self):
        import matplotlib.pyplot as plt
        
        plt.figure()
        
        plt.plot([0, self.Q_dot_1 ]                          , [self.su.T, self.T_sat]  , 'g', label='Fluid')
        plt.plot([self.Q_dot_1, self.Q_dot_1+self.Q_dot_2]   , [self.T_sat, self.T_sat], 'g')
        plt.plot([self.Q_dot_1+self.Q_dot_2, self.Q]         , [self.T_sat, self.ex.T]  , 'g')

        plt.plot([0, self.Q]                                 , [self.sto_fluid.T, self.sto_fluid.T], 'r', label='Storage')

        plt.xlabel("Q_dot [W]")
        plt.ylabel("T [K]")
        plt.grid()
        plt.legend()
        plt.show()