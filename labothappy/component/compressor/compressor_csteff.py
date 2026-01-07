import __init__

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP

class CompressorCstEff(BaseComponent):
    """
    **Component**: Compressor

    **Model**: Constant isentropic efficiency

    **Descritpion**:

        This model determines the exhaust specific enthalpy and the exhaust temperature of a compressor. This model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.
        - Isentropic efficiency stays constant for all the conditions.

    **Connectors**:

        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

        W (WorkConnector): Work connector for the mechanical work.

    **Parameters**:

        eta_is: Isentropic efficiency. [-]

    **Inputs**:

        P_su: Suction side pressure. [Pa]

        T_su: Suction side temperature. [K]

        P_ex: Exhaust side pressure. [Pa]

        fluid: Working fluid. [-]

        m_dot: Mass flow rate of working fluid. [kg/s]

    **Ouputs**:

        h_ex: Exhaust side specific enthalpy. [J/kg] 

        T_ex: Exhaust side temperature. [K]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Mass_connector for the suction side
        self.ex = MassConnector() # Mass_connector for the exhaust side
        self.W = WorkConnector()

    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su', 'T_su', 'P_ex', 'fluid', 'm_dot']

    def get_required_parameters(self):
        return [
            'eta_is',
        ]

    def solve(self):
        print('coucou')
        self.check_calculable()
        self.check_parametrized()
        self.print_setup()
        self.AS = CP.AbstractState('HEOS', self.su.fluid)

        try:
            self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
            h_ex_is = self.AS.hmass()
            
            self.AS.T()
            
            h_ex = self.su.h + (h_ex_is - self.su.h) / self.params['eta_is']
            w = h_ex - self.su.h
            W_dot = self.su.m_dot*w
            self.update_connectors(h_ex, w, W_dot)

            self.solved = True
        except Exception as e:
            print(f"Error: {e}")
            self.solved = False
            return
    
    def update_connectors(self, h_ex, w, W_dot):
        self.ex.set_h(h_ex)
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.W.set_w(w)
        self.W.set_W_dot(W_dot)

    def print_results(self):
        print("=== Compressor Results ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, h={self.su.h}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, h={self.ex.h}")

        print("\nResults:")
        print(f"  - h_ex: {self.ex.h}")
        print(f"  - T_ex: {self.ex.T}")
        print("=========================")

    def print_states_connectors(self):
        print("=== Compressor Results ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")






