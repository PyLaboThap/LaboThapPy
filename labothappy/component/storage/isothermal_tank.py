import __init__

#%%

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

from scipy.optimize import fsolve, root, minimize
from CoolProp.CoolProp import PropsSI

import CoolProp.CoolProp as CP
import numpy as np
import math

class IsothermalTank(BaseComponent):
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

    def __init__(self, fluid):
        super().__init__()
        self.su = MassConnector() # Working fluid supply
        self.ex = MassConnector()
        self.sto_fluid = MassConnector()
        self.Q_dot = HeatConnector()

        self.AS = CP.AbstractState('HEOS', fluid)

        self.su.set_fluid(fluid)
        self.ex.set_fluid(fluid)
        self.sto_fluid.set_fluid(fluid)

        self.mode = 'Charge'
        self.FR = 0 # Filling Rate
        self.SOC = 0 # State of Charge

        self.init = {}
        self.initialized = False

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required 
        
        if self.mode == 'Charge':
            return ['fluid', 'T_su', 'P_su', 'm_dot']

        elif self.mode == 'Discharge':
            return ['fluid', 'T_ex', 'P_ex', 'm_dot']
            
        else:
            raise ValueError("Mode is not 'Charge' nor 'Discharge'")
            
        return

    def get_required_parameters(self):
        return [
            'Vol', # Volume
            'HeatLoss', # Heat loss
            'DT' # Time interval in seconds
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

    def set_initial_conditions(self, T_init, P_init, FR_init):
        """Set the initial conditions for the thermocline."""
        
        # Get temperature at the triple point (freezing/melting point)
        self.sto_fluid.set_T(T_init)
        self.sto_fluid.set_p(P_init)            
        
        self.FR = FR_init
        self.SOC = self.FR
        
        self.init['P'] = P_init
        self.init['T'] = T_init
        self.init['FR'] = FR_init
        
        self.initialized = True
        
        return

    def set_mode(self, mode):
        """Set the initial conditions for the thermocline."""
        
        self.mode = mode
        return

    def system(self):
        
        if self.mode == 'Charge':
            self.AS.update(CP.PT_INPUTS, self.su.p, self.su.T)
            
            rho_in = self.AS.rhomass()
            Vdot_in = self.su.m_dot/rho_in
            vol_init = self.FR*self.params['Vol']
            vol_end = min(vol_init + Vdot_in*self.params['DT'], self.params['Vol'])
            self.FR = min(1, vol_end/self.params['Vol'])
            
            h_in = self.AS.hmass()
            self.AS.update(CP.PT_INPUTS, self.sto_fluid.p, self.sto_fluid.T)
            h_init = self.AS.hmass()
            rho_init = self.AS.rhomass()
            h_end = (vol_init*rho_init*h_init + h_in*self.su.m_dot*self.params['DT'])/(vol_init*rho_init + self.su.m_dot*self.params['DT'])
            self.sto_fluid.set_h(h_end)
            self.sto_fluid.set_p(self.sto_fluid.p)
            self.Q_dot.set_Q_dot(h_in*self.su.m_dot)
        
        elif self.mode == 'Discharge':
            
            self.AS.update(CP.PT_INPUTS, self.ex.p, self.ex.T)
            
            rho_out = self.AS.rhomass()
            Vdot_out = self.ex.m_dot*rho_out
            vol_init = self.FR*self.params['Vol']
            vol_end = vol_init - Vdot_out*self.params['DT']
            self.FR = min(0, vol_end/self.params['Vol'])
            
            h_out = self.AS.hmass()
            self.AS.update(CP.PT_INPUTS, self.sto_fluid.p, self.sto_fluid.T)
            h_init = self.AS.hmass()
            rho_init = self.AS.rhomass()
            h_end = (vol_init*rho_init*h_init - h_out*self.ex.mdot*self.params['DT'])/(vol_init*rho_init - self.ex.mdot*self.params['DT'])
            self.sto_fluid.set_h(h_end)
            self.sto_fluid.set_p(self.sto_fluid.p)
            self.Q_dot.set_Q_dot(h_out*self.ex.mdot)
            
        else:
            raise ValueError("Mode is not 'Charge' nor 'Discharge'")
            
        T_max = max(self.su.T, self.sto_fluid.T)
        self.AS.update(CP.PT_INPUTS, self.sto_fluid.p, T_max)
        h_max = self.AS.hmass()
        rho_T_max = self.AS.rhomass()
        
        self.SOC = (self.sto_fluid.h*self.FR*self.params['Vol']*self.sto_fluid.D) / (h_max*self.params['Vol']*rho_T_max)
              
        return

    def solve(self):
        # Ensure all required checks are performed

        self.check_calculable()
        self.check_parametrized()

        if not self.calculable:
            print("IsothermalTank is not calculable, check inputs.")
            return
        
        if not self.parametrized:
            print("IsothermalTank is not parametrized, check parameters.")
            return

        if not self.initialized:
            print("IsothermalTank is not initialized, check initial conditions.")
            return   
        
        self.system()
        self.solved = True
        
        # self.update_connectors()
        return

    # def update_connectors(self):
        
    #     "Mass Connectors"

    #     self.sto_fluid.set_h(self.su.fluid)
    #     self.ex.set_T(self.T_ex)
    #     self.ex.set_p(self.P_sat)
    #     self.ex.set_m_dot(self.su.m_dot)
        
    #     "Heat conector"
    #     self.Q_dot.set_Q_dot(self.Q)

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

test = 'Discharge'

if test == 'Charge':        
    tank = IsothermalTank('Water')
    
    tank.set_inputs(
        fluid = 'Water',
        T_su = 95 + 273.15,
        P_su = 10*1e5,
        m_dot = 1
        )
    
    tank.set_parameters(
        Vol = 3,
        DT = 60,
        HeatLoss = 0
        )
    
    tank.set_mode("Discharge")
    
    tank.set_initial_conditions(95 + 273.15, 10*1e5, 0.5)
    
    tank.solve()

elif test == 'Discharge':
    tank = IsothermalTank('Water')
    
    tank.set_inputs(
        fluid = 'Water',
        T_ex = 95 + 273.15,
        P_ex = 10*1e5,
        m_dot = 1
        )
    
    tank.set_parameters(
        Vol = 3,
        DT = 60,
        HeatLoss = 0
        )
    
    tank.set_mode("Discharge")
    
    tank.set_initial_conditions(95 + 273.15, 10*1e5, 0.5)
    
    tank.solve()

