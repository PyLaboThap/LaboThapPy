# -*- coding: utf-8 -*-
import numpy as np

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

import numpy as np

class PumpPolyEff(BaseComponent):
    """
        Component: Peripherical pump with polynomial efficiency model

        Model: The model is based on the characteristic curves of the pump. The efficiency is calculated based on the volume flow rate and the head.

        **Descritpion**:

            This model is used to simulate the performance of a peripherical pump. 

        **Assumptions**:

            - Steady-state operation.

        **Connectors**:

            su (MassConnector): Mass connector for the suction side.

            ex (MassConnector): Mass connector for the exhaust side.

            W_pp (WorkConnector): Work connector.

        **Parameters**:

            /    
        
        **Inputs**:

            su_p: Suction side pressure. [Pa]

            su_T: Suction side temperature. [K]

            ex_p: Exhaust side pressure. [Pa]

            m_dot: Mass flow rate. [kg/s]

            su_fluid: Suction side fluid. [-]

        **Ouputs**:

            ex_h: Exhaust side specific enthalpy. [J/kg]

            ex_T: Exhaust side temperature. [K]

            W_dot_pp: Pump power. [W]

            N_pp: Pump rotational speed. [Hz]

    """
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_pp = WorkConnector()

    def get_required_inputs(self):
            self.sync_inputs()
            # Return a list of required inputs
            return ['su_p', 'su_T', 'ex_p', 'm_dot', 'su_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.su.m_dot is not None:
            self.inputs['m_dot'] = self.su.m_dot
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])
        if 'm_dot' in self.inputs:
            self.su.set_m_dot(self.inputs['m_dot'])


    def get_required_parameters(self):
        return [
        ]
    
    def print_setup(self):
        print("=== Pump Setup ===")
       
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

    def frequency(self, Head_value, V_dot_value):
        # For the frequency w=60Hz, we can determine the Head corresponding to the volume flow rate based on the graph
        w1 = 60
        Head_w1 = -23.8 * V_dot_value + 94.76

        # Using the similitudes, we can now determine the real frequency
        w = np.sqrt(Head_value/Head_w1) * w1
        return w
    
    def pump_consumption(self, w_value, V_dot_value):      
        # For the frequency w=60Hz, we can determine the Power corresponding to the volume flow rate
        w1 = 60
        W_dot_pp_w1 = -0.2 * V_dot_value + 1.8

        # Using the similitudes, we can now determine the real Power
        W_dot_pp_kW = W_dot_pp_w1 * (w_value/w1)**3
        W_dot_pp = W_dot_pp_kW*1000
        return W_dot_pp

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("PumpPolyEff could not be solved. It is not calculable and/or not parametrized")
            return

        try:
            # Calculate the necessary parameters for the pump model
            m_dot = self.su.m_dot  # Mass flow rate at the suction point
            P_su = self.su.p  # Pressure at the suction point
            h_su = self.su.h  # Enthalpy at the suction point
            P_ex = self.ex.p  # Pressure at the discharge point
            rho_su = self.su.D  # Density at the suction point
            
            "Transform for water"
            rho_water = 1000  # Density of water [kg/m^3]
            m_dot_water = m_dot * rho_water / rho_su  # Mass flow rate of water [kg/s]
            V_dot_water = (m_dot_water/1000)*3600 # Volume flow rate of water [m^3/h]

            DeltaP_ref = P_ex - P_su  # Pressure difference across the pump
            DeltaP_water = DeltaP_ref * rho_water / rho_su # Pressure difference across the pump for water [Pa]
            g = 9.81  # Acceleration due to gravity

            Head_water = DeltaP_water / (1000 * g)  # Calculate the head based on pressure difference and density

            # Map the frequency and pump consumption based on the calculated head and mass flow rate
            self.w = self.frequency(Head_water, V_dot_water)  # Frequency based on head and volume flow rate [Hz]

            self.W_dot_pp = self.pump_consumption(self.w, V_dot_water)  # Pump consumption based on frequency and volume flow rate [kW]

            h_ex = h_su + self.W_dot_pp / m_dot  # Calculate the enthalpy at the discharge point
            self.h_ex = h_ex

            if h_ex == None:
                self.convergence = False
            else:
                self.convergence = True

            if self.convergence:
                self.update_connectors()
                self.solved = True

        except Exception as e:
            print(f"PumpPolyEff could not be solved. Error: {e}")
            self.solved = False

    def update_connectors(self):
        """Update the connectors with the calculated values."""
        
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_h(self.h_ex)

        self.W_pp.set_W_dot(self.W_dot_pp)
        self.W_pp.set_N(self.w*60)

    def print_results(self):
        print("=== Pump Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_pp: {self.W_pp.W_dot} [W]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Pump Results ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")
        print("Work connector:")
        print(f"  - W_dot_pp: {self.W_pp.W_dot} [W]")
        print(f"  - N_pp: {self.W_pp.N} [RPM]")
        print("=========================")



