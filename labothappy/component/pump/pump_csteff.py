
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

import CoolProp.CoolProp as CP

class PumpCstEff(BaseComponent):
    """
    **Component**: Pump

    **Model**: Constant efficiency

    **Descritpion**:

        This model determines the exhaust specific enthalpy and the exhaust temperature of a pump. This model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.
        - Efficiency stays constant for all the conditions.

    **Connectors**:

        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

        W (WorkConnector): Work connector for the pump power consumption.

    **Parameters**:

        eta_is: Isentropic efficiency. [-]

    **Inputs**:

        P_su: Suction side pressure. [Pa]

        T_su: Suction side temperature. [K]
        
        m_dot : Suction side flowrate [kg/s]

        P_ex: Exhaust side pressure. [Pa]

        fluid: Suction side fluid. [-]

    **Ouputs**:

        h_ex: Exhaust side specific enthalpy. [J/kg] 

        T_ex: Exhaust side temperature. [K]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W = WorkConnector()

    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su', 'T_su', 'P_ex', 'fluid']

    def get_required_parameters(self):
        return ['eta_is']

    def solve(self):
        # Perform checks to ensure the model can be calculated and has parameters
        self.check_calculable()
        self.check_parametrized()
            
        if not (self.calculable and self.parametrized):
            print("Pump model is either not calculable or not parametrized.")
            self.solved = False
            return

        try:
            self.AS = CP.AbstractState('HEOS', self.su.fluid)
            """PUMP MODEL"""
            # Calculate the outlet enthalpy based on isentropic efficiency
            self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
            h_ex_is = self.AS.hmass()
            h_ex = self.su.h + (h_ex_is - self.su.h) / self.params['eta_is']
            w_pp = h_ex - self.su.h
            W_dot_pp = self.su.m_dot*w_pp

            """Update connectors after the calculations"""
            self.update_connectors(h_ex, w_pp, W_dot_pp)

            # Mark the model as solved if successful
            self.solved = True
        except Exception as e:
            # Handle any errors that occur during solving
            self.solved = False
            print(f"Convergence problem in pump model: {e}")

        
    def update_connectors(self, h_ex, w_pp, W_dot_pp):
        """Update the connectors with the calculated values."""
        self.ex.set_h(h_ex)
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.W.set_w(w_pp)
        self.W.set_W_dot(W_dot_pp)

    def print_results(self):
        print("=== Pump Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_pp: {self.W.W_dot} [W]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Pump States ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")
        print("Work connector:")
        print(f"  - W_dot_pp: {self.W.W_dot} [W]")
        print("=========================")
