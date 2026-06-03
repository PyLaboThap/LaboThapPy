import __init__

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import minimize, brentq

class CompressorRadialMeanLine(BaseComponent):
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

        self._STATE_KEYS = ('H', 'S', 'P', 'D', 'A', 'V')

        self.total_states  = {k: {i: np.nan for i in range(1, 6)} for k in self._STATE_KEYS}
        self.static_states  = {k: {i: np.nan for i in range(1, 6)} for k in self._STATE_KEYS}

    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su', 'T_su', 'N_rot', 'fluid', 'm_dot']

    def get_required_parameters(self):
        return []

#%%

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

#%%

    def solve_Rotor(self):

        self.impeller_inlet_blockage = (
            1 - self.params['n_blade_R'] * self.params['t_b']
            / (np.pi * self.params['r1'] * 2 * np.sin(np.pi * abs(self.inputs['xhi1']) / 180))
        )
        self.A1_th = np.pi * (self.params['r1s']**2 - self.params['r1h']**2) * self.impeller_inlet_blockage

        self.Vel_Tri_R['alpha1'] = alpha1 = self.params['xhi_igv']
        self.Vel_Tri_R['u1'] = u1 = self.params['N_rot'] * self.params['r1']
        self.u1h = self.params['N_rot'] * self.params['r1h']
        self.u1s = self.params['N_rot'] * self.params['r1s']

        def compute_h1_new(h1): 
            self.update_static_AS(CP.HmassSmass_INPUTS, h1, self.total_states['S'][1], 1)
            
            self.Vel_Tri_R['vm1'] = vm1 = self.inputs['mdot'] / (self.static_states['D'][1] * self.A1_th)
            self.Vel_Tri_R['vu1'] = vu1 = vm1*np.tan(self.Vel_Tri_R['alpha1'])
            self.Vel_Tri_R['v1']  = v1  = np.sqrt(vm1**2 + vu1**2)

            self.Vel_Tri_R['wu1'] = wu1 = vu1 - self.Vel_Tri_R['u1']
            self.Vel_Tri_R['w1']  = w1  = np.sqrt(vm1**2 + wu1**2)
            self.Vel_Tri_R['beta1'] = np.arccos(vm1/w1)
            
            self.beta1h = np.arctan(self.u1h / vm1)
            self.beta1s = np.arctan(self.u1s / vm1)
            
            self.w1s = vm1 / np.cos(self.beta1s)
            self.w1h = vm1 / np.cos(self.beta1h)
                        
            h1_new = self.total_states['H'][1] - v1**2 / 2
                        
            return h1_new

        def residual_h1(h1):
            return h1 - compute_h1_new(h1)
        
        h1_min = self.total_states['H'][1] * 0.95
        h1_max = self.total_states['H'][1]
        
        h1_solution = brentq(residual_h1, h1_min, h1_max, xtol=1e-6)

#%%

    def solve(self):
        
        # 0) Setup the solving
        
        self.check_calculable()
        self.check_parametrized()
        # self.print_setup()
        
        if not self.calculable:
            self.solved = False
            print("CompressorRadialMeanLine could not be solved. It is not calculable.")
            return
            
        if not self.parametrized:
            self.solved = False
            print("CompressorRadialMeanLine could not be solved. It is not parametrized.")
            return
        
        # 1) Setup the inlet

        self.AS = CP.AbstractState('HEOS', self.su.fluid)
        self.update_total_AS(CP.HmassP_INPUTS, self.su.h, self.su.p, 1)
        
        # 2) Solve 
        self.solve_Rotor()
        
        # try:
        #     self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        #     self.h_ex_is = self.AS.hmass()
            
        #     self.AS.T()
            
        #     h_ex = self.su.h + (self.h_ex_is - self.su.h) / self.params['eta_is']
        #     w = h_ex - self.su.h
        #     W_dot = self.su.m_dot*w
        #     self.update_connectors(h_ex, w, W_dot)

        #     self.solved = True
        #     # self.print_states_connectors()
        # except Exception as e:
        #     print(f"Error: {e}")
        #     self.solved = False
        #     return
    
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
    
    Comp = CompressorRadialMeanLine()

    Comp.set_inputs(
        fluid = 'CO2',
        m_dot  = 2.15,
        P_su = 76.9*1e5,
        T_su = 305.97,
        N_rot  = 55000,
    )

    Comp.set_parameters(
        xhi_igv = 1.2095095940256784
        )

    Comp.solve()
