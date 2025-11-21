import __init__

#%%

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, root, minimize
import numpy as np
import math

class StorageLatentIsothermalCstePinch(BaseComponent):
    """
    **Component**: Heat Exchanger with constant pinch point.

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
        
        fluid: Working fluid. [-]

        h_su: Suction side enthalpy. [J/kg]

        T_su: Suction side temperature. [K]

        m_dot: Mass flow rate of working fluid. [kg/s]

    **Outputs**:

        h_ex: Exhaust side enthalpy. [J/kg]

        T_ex: Exhaust side temperature. [K]

        P_sat: Saturation pressure at T_ex. [Pa]

        Q_dot: Heat transfer rate. [W]

    """ 
        
    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Working fluid supply
        self.ex = MassConnector()

        self.sto_fluid = MassConnector()
        self.DP = 0

        self.Q_dot = HeatConnector()

    def get_required_inputs(self): 
        # Return a list of required inputs
        return ['fluid', 'sto_fluid', 'h_su', 'T_su', 'm_dot']
    
    def get_required_parameters(self):
        return [
            'Pinch', # pinch point
            'Delta_T_sh_sc', # Superheating or subcooling
        ]
    
    def get_required_guesses(self):
        return []

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
        
        self.AS_sto = CP.AbstractState('HEOS', self.inputs['sto_fluid'])
        self.AS = CP.AbstractState('HEOS', self.su.fluid)
        
        # Get temperature at the triple point (freezing/melting point)
        try:    
            self.sto_fluid.set_T(self.params['T_sto'])
        except:
            self.sto_fluid.set_T(self.AS_sto.Ttriple())
        
        self.sto_fluid.set_p(P_sto)

        if 'DP' in self.params:
            self.DP = self.params['DP']

        if self.su.T <= self.sto_fluid.T:
            if  self.su.T - (self.sto_fluid.T - self.params['Pinch'] - self.params['Delta_T_sh_sc']) <= 1e-2:
                self.T_ex = self.sto_fluid.T - self.params['Pinch']
                self.T_sat = self.T_ex - self.params['Delta_T_sh_sc']
                
                self.AS.update(CP.QT_INPUTS, 0, self.T_sat)
                self.P_sat = self.AS.p()

                self.AS.update(CP.PT_INPUTS, self.P_sat, self.T_ex)                
                self.h_ex = self.AS.hmass()

                self.AS.update(CP.QT_INPUTS, 1, self.T_sat)
                self.h_sat_v = self.AS.hmass()
                                
                self.AS.update(CP.QT_INPUTS, 0, self.T_sat)
                self.h_sat_l = self.AS.hmass()
                
                self.Q_dot_3 = self.su.m_dot*(self.h_ex - self.h_sat_v)
                self.Q_dot_2 = self.su.m_dot*(self.h_sat_v - self.h_sat_l)

                try:
                    self.AS.update(CP.PT_INPUTS, self.P_sat,self.su.T)
                    h_su = self.AS.hmass()
                except:
                    h_su = self.su.h

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
                
                self.AS.update(CP.QT_INPUTS, 0, self.T_sat)
                self.P_sat = self.AS.p()                
                
                self.AS.update(CP.PT_INPUTS, self.P_sat, self.T_ex)                
                self.h_ex = self.AS.hmass()
                
                self.AS.update(CP.QT_INPUTS, 1, self.T_sat)
                self.h_sat_v = self.AS.hmass()
                                
                self.AS.update(CP.QT_INPUTS, 0, self.T_sat)
                self.h_sat_l = self.AS.hmass()
                
                self.Q_dot_3 = self.su.m_dot*(self.h_sat_l - self.h_ex)
                self.Q_dot_2 = self.su.m_dot*(self.h_sat_v - self.h_sat_l)

                try:
                    self.AS.update(CP.PT_INPUTS, self.P_sat,self.su.T)
                    h_su = self.AS.hmass()
                except:
                    h_su = self.su.h
                    
                self.Q_dot_1 = self.su.m_dot*(h_su - self.h_sat_v)

                self.Q = self.Q_dot_1 + self.Q_dot_2 + self.Q_dot_3
                
            else: # Pinch and SC_SH cannot be satisfied
                # raise ValueError("Pinch and SC_SH are not satisfied")
                return
        
        self.solved = True
        self.update_connectors()


    def update_connectors(self):
        
        "Mass Connectors"

        self.su.set_p(self.P_sat)

        self.ex.set_fluid(self.su.fluid)
        self.ex.set_T(self.T_ex)
        self.ex.set_p(self.P_sat - self.DP)
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