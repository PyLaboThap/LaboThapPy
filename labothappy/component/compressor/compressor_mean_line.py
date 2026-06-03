import __init__

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import numpy as np

class RadialCompressorMeanLine(BaseComponent):
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

        _STATE_KEYS = ('H', 'S', 'P', 'D', 'A', 'V')

        # State dicts – replaces pd.DataFrame(columns=[…], index=[1,2,3,4,5])
        self.total_states  = {k: {i: np.nan for i in range(1, 6)} for k in _STATE_KEYS}
        self.static_states = {k: {i: np.nan for i in range(1, 6)} for k in _STATE_KEYS}
        
    
    def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.total_states['H'][position] = self.AS.hmass()            
        self.total_states['S'][position] = self.AS.smass()            
        self.total_states['P'][position] = self.AS.p()            
        self.total_states['D'][position] = self.AS.rhomass()            

        try:        
            self.total_states['A'][position] = self.AS.speed_sound()
        except:
            self.total_states['A'][position] = -1  
            
        self.total_states['V'][position] = self.AS.viscosity()
        
        return
    
    def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.static_states['H'][position] = self.AS.hmass()
        self.static_states['S'][position] = self.AS.smass()
        self.static_states['P'][position] = self.AS.p()
        self.static_states['D'][position] = self.AS.rhomass()
        
        try:        
            self.static_states['A'][position] = self.AS.speed_sound()
        except:
            self.static_states['A'][position] = -1
            
        self.static_states['V'][position] = self.AS.viscosity()

        return
    
    # def get_required_inputs(self):
    #     # Return a list of required inputs
    #     return ['P_su', 'T_su', 'P_ex', 'fluid', 'm_dot']

    # def get_required_parameters(self):
    #     return [
    #         'eta_is',
    #     ]

    def solve_rotor(self):
        
        return

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if not self.calculable:
            if self.print_flag:
                print("RadialCompressorMeanLine could not be solved. It is not calculable.")
            self.solved = False
            return
        
        if not self.parametrized:
            if self.print_flag:
                print("RadialCompressorMeanLine could not be solved. It is not parametrized.")
            self.solved = False
            return
        
        self.AS = CP.AbstractState('HEOS', self.su.fluid)

        self.update_total_AS(CP.PT_INPUTS, self.inputs['P_su'], self.inputs['T_su'], 1)
        self.update_static_AS(CP.PT_INPUTS, self.inputs['P_su'], self.inputs['T_su'], 1)

        self.solve_rotor()

        return
    
    def update_connectors(self, h_ex, w, W_dot):
        
        self.ex.reset()
        
        self.ex.set_h(h_ex)
        self.ex.set_p(self.inputs['P_ex'])
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

if __name__ == "__main__":
    
    Comp = RadialCompressorMeanLine()
    
    Comp.set_inputs(
        fluid = 'CO2',
        P_su = 76.9*1e5, # Total inlet P : Pa
        T_su = 305.97, # Total inlet T : K
        m_dot  = 2.15, # kg/s
        N_rot = 50000, # RPM
        )

    Comp.set_parameters()
    
    Comp.solve()
    
    
    


