from component.base_component import BaseComponent
from connector.mass_connector import MassConnector


class ValveIsenthalpic(BaseComponent):
    """
    **Component**: Valve

    **Model**: Isenthalpic

    **Descritpion**:

        This model determines the outlet conditions for an isenthalpic valve. This model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.

    **Connectors**:

        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

    **Inputs**:

        P_su: Suction side pressure. [Pa]

        T_su: Suction side temperature. [K]

        P_ex: Exhaust side pressure. [Pa]

        fluid: Suction side fluid. [-]

    **Ouputs**:

        h_ex: Exhaust side specific enthalpy. [J/kg] 

        T_ex: Exhaust side temperature. [K]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Mass_connector for the suction side
        self.ex = MassConnector() # Mass_connector for the exhaust side

    def get_required_inputs(self):
            # Return a list of required inputs
            return ['P_su', 'T_su', 'P_ex', 'fluid']

    def get_required_parameters(self):
        return []

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("IsenthalpicValve could not be solved. It is not calculable and/or not parametrized")
            return
        try:
            
            h_ex = self.su.h
            self.update_connectors(h_ex)

            self.solved = True
        except Exception as e:
            print(f"Error: {e}")
            self.solved = False
            return
    
    def update_connectors(self, h_ex):
        self.ex.set_h(h_ex)
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)

    def print_results(self):
        print("=== Expansion Valve Results ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nResults:")
        print(f"  - h_ex: {self.ex.h}")
        print(f"  - T_ex: {self.ex.T}")
        print("=========================")

    def print_states_connectors(self):
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")

