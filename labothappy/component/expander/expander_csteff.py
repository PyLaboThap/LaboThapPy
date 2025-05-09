
import numpy as np
from component.base_component import BaseComponent
from connector.heat_connector import HeatConnector
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from CoolProp.CoolProp import PropsSI


class ExpanderCstEff(BaseComponent):
    """
    Component: Expander

    Model: Constant isentropic efficiency model

    Reference: /

    **Description**:

        This model determines the exhaust specific enthalpy and the exhaust temperature of an expander. This model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.
        - Isentropic efficiency stays constant for all the conditions.

    **Connectors**:

        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

    **Parameters**:

        eta_is: Isentropic efficiency. [-]

    **Inputs**:

        P_su: Suction side pressure. [Pa]

        T_su: Suction side temperature. [K]
        
        m_dot: Mass flow flow rate. [kg/s]

        P_ex: Exhaust side pressure. [Pa]

        fluid: fluid. [-]

    **Ouputs**:

        h_ex: Exhaust side specific enthalpy. [J/kg]

        T_ex: Exhaust side temperature. [K]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()  # Mass_connector
        self.W_exp = WorkConnector()

    def get_required_inputs(self):  # Used in check_calculablle to see if all of the required inputs are set
        # Return a list of required inputs
        return ["P_su", "T_su", "P_ex", "m_dot", "fluid"]

    def get_required_parameters(self):
        return ["eta_is"]

    def solve(self):
        # Perform checks to ensure the model can be calculated and has parameters
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print(
                "ExpanderCstEff could not be solved. It is not calculable and/or not parametrized"
            )
            return
        try:
            """EXPANDER MODEL"""

            # Calculate the outlet enthalpy based on isentropic efficiency
            h_ex_is = PropsSI("H", "P", self.ex.p, "S", self.su.s, self.su.fluid)
            h_ex = self.su.h - (self.su.h - h_ex_is) / self.params["eta_is"]
            w_exp = self.su.h - h_ex
            self.ex.set_m_dot(self.su.m_dot)
            W_dot_exp = self.su.m_dot * w_exp

            # Update connectors after the calculations
            self.update_connectors(h_ex, w_exp, self.ex.p, W_dot_exp)

            # Mark the model as solved if successful
            self.solved = True

        except Exception as e:
            # Handle any errors that occur during solving
            self.solved = False
            print(f"Solving problem in expander model: {e}")
            return

    def update_connectors(self, h_ex, w_exp, p_ex, W_dot_exp):
        """Update the connectors with the calculated values."""
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_h(h_ex)
        self.ex.set_p(p_ex)
        self.ex.set_m_dot(self.su.m_dot)
        self.W_exp.set_W_dot(W_dot_exp)

    def print_results(self):
        print("=== Expander Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_exp: {self.W_exp.W_dot} [W]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Expander Results ===")
        print("Mass connectors:")
        print(
            f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]"
        )
        print(
            f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]"
        )
        print("=========================")
        print("Work connector:")
        print(f"  - W_dot_exp: {self.W_exp.W_dot} [W]")
        print("=========================")


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
