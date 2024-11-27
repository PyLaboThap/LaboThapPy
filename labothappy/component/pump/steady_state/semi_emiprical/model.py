"""
Author: Elise Neven
email: elise.neven@uliege.be
Date: 18/11/2024
"""

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

from CoolProp.CoolProp import PropsSI
import math


class PumpSE(BaseComponent):
    """
    Component: Pump

    Model: Based on the Thesis of Rémi Dickes. Semi-empirical model.

    **Description**:
        This model is used to simulate the performances and the behavior of pumps.
        The parameters need to be calibrated with experimental data.

    **Assumption**:
        Steady-state operation.

    **Connectors**:
        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

        W_pp (WorkConnector): Work connector to connect to a potential motor model.

    **Parameters**:
        V_s: machine displacement volume [m^3]
        V: volume of the pump [m^3]
        A_leak: leakage surface area [m2]
        W_dot_loss: constant power losses [W]
        K_0_loss: term for the proportional losses [-]

    **Inputs**:
        su_p: inlet pressure of the WF [Pa]
        su_h: inlet temperature of the WF [J/kg]
        ex_p: outlet pressure of the WF [Pa]
        fluid: nature of the fluid (string) [-]
        N_pp: the pump speed [rpm]

    **Outputs**:
        ex_T =  exhaust temperature [K]
        ex_h =  exhaust enthalpy [J/kg]
        m_dot = fluid mass flow rate [kg/s]
        W_dot = mechanical power [W]
        epsilon_is = isentropic efficiency [-]
        epsilon_vol = volumetric efficiency [-]
        M = mass of fluid inside the pump [kg]
    
    """
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_pp = WorkConnector()

    def get_required_inputs(self):
            self.sync_inputs()
            # Return a list of required inputs
            return ['su_p', 'su_h', 'ex_p', 'N_pp', 'fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p
        if self.W_pp.N is not None:
            self.inputs['N_pp'] = self.W_pp.N

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'fluid' in self.inputs:
            self.su.set_fluid(self.inputs['fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])
        if 'N_pp' in self.inputs:
            self.W_pp.set_N(self.inputs['N_pp'])

    def get_required_parameters(self):
        return [
            'V_s', 'V', 'A_leak', 'W_dot_loss', 'K_0_loss'
        ]
    
    def print_setup(self):
        print("=== Expander Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - W_pp: N={self.W_pp.N}, W_dot={self.W_pp.W_dot}")

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
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("PumpSE could not be solved. It is not calculable and/or not parametrized")
            return
        
        # PUMP MODEL
        try:
            T_su = PropsSI('T','P',self.su.p,'H',self.su.h,self.su.fluid)
            s_su = PropsSI('S','P',self.su.p,'H',self.su.h,self.su.fluid)
            rho_su = PropsSI('D','P',self.su.p,'H', self.su.h, self.su.fluid)
            h_ex_s = PropsSI('H','P', self.ex.p,'S',self.su.s,self.su.fluid)
            self.m_dot = max(1e-3, (self.W_pp.N/60*self.V_s*rho_su)-(self.A_leak*math.sqrt(2*rho_su*(self.ex.p-self.su.p))));
            self.epsilon_vol = self.m_dot/(self.W_pp.N/60*self.V_s*rho_su)
            self.W_dot = self.W_dot_loss + self.K_0_loss*self.m_dot/rho_su*(self.ex.p-self.su.p);
            self.epsilon_is = (self.m_dot*(h_ex_s-self.su.h))/self.W_dot
            self.h_ex = self.su.h+self.W_dot/self.m_dot

            h_min = PropsSI('H','P',5e4,'T',253.15,self.su.fluid)
            h_max =  PropsSI('H','P',4e6,'T',500,self.su.fluid)
            if self.h_ex > h_min and self.h_ex < h_max:
                self.solved = True
                self.update_connectors()
            else:
                self.solved = False

        except Exception as e:
            print(f"PumpSE could not be solved. Error: {e}")
            self.solved = False

    def update_connectors(self):
        """Update the connectors with the calculated values."""
        self.su.set_m_dot(self.m_dot)
        self.ex.set_m_dot(self.m_dot)
        
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.m_dot)
        self.ex.set_h(self.h_ex)

        self.W_pp.set_W_dot(self.W_dot)

    def print_results(self):
        print("=== Expander Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_exp: {self.W_exp.W_dot} [W]")
        print(f"  - epsilon_is: {self.epsilon_is} [-]")
        print(f"  - m_dot: {self.m_dot} [kg/s]")
        print(f"  - epsilon_v: {self.epsilon_v} [-]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Expander Results ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")
        print("Work connector:")
        print(f"  - W_dot_exp: {self.W_exp.W_dot} [W]")
        print("=========================")
        print("Heat connector:")
        print(f"  - Q_dot_amb: {self.Q_amb.Q_dot} [W]")
        print(f"  - T_hot: {self.Q_amb.T_hot} [K]")
        print(f"  - T_cold: {self.Q_amb.T_cold} [K]")
        print("=========================")
