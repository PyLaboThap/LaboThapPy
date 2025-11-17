import sys
import os

# Get the absolute path of the directory that contains the script (simulation_model.py)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Determine the project root directory (which contains both 'connector' and 'component')
project_root = os.path.abspath(os.path.join(current_dir, '..', '..', '..'))

# Add the project root to sys.path if it's not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)

#%%

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

# from component.heat_exchanger.moving_boundary.simple_model.modules.U import U_Gnielinski_calibrated, U_DittusBoelter, U_Cooper_calibrater, U_Thonon

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

import CoolProp.CoolProp as CP
import numpy as np
import math

class HXEffCstDisc(BaseComponent):
    """
    Component: Counterflow Heat Exchanger with Constant Effectiveness (HXEffCstDisc)
    
    Model: Discretized Counterflow Heat Exchanger with Fixed Effectiveness and Pinch Check
    
    **Description**:
    
        This model simulates a counterflow heat exchanger with constant effectiveness, discretized into segments along the flow direction. The model uses energy balances and basic thermodynamics (via CoolProp) to determine outlet conditions for both hot and cold streams. It iteratively adjusts the effectiveness to satisfy a minimum pinch point temperature difference (Pinch_min). The model is suitable for on-design, steady-state simulations of heat exchangers where detailed flow dynamics are simplified.
    
    **Assumptions**:
    
        - Steady-state operation.
        - Constant effectiveness per iteration, adjusted to meet Pinch_min.
        - Uniform discretization of the heat exchanger along its length.
        - Uniform mass flow rates at inlets and outlets.
        - Fluid properties are obtained via CoolProp.
    
    **Connectors**:
    
        su_C (MassConnector): Mass connector for the cold-side supply.
        su_H (MassConnector): Mass connector for the hot-side supply.
        ex_C (MassConnector): Mass connector for the cold-side exhaust.
        ex_H (MassConnector): Mass connector for the hot-side exhaust.
        Q_dot (HeatConnector): Heat transfer connector for the total exchanged heat.
    
    **Parameters**:
    
        eta (float): Initial effectiveness of the heat exchanger [-].
        
        n_disc (int): Number of discretization segments along the exchanger length.
        
        Pinch_min (float): Minimum allowable pinch temperature difference [K].
    
    **Inputs**:
    
        su_H_fluid (str): Hot-side fluid.
        
        su_H_h (float): Hot-side inlet specific enthalpy [J/kg].
        
        su_H_p (float): Hot-side inlet pressure [Pa].
        
        su_H_m_dot (float): Hot-side mass flow rate [kg/s].
    
        su_C_fluid (str): Cold-side fluid.
        
        su_C_h (float): Cold-side inlet specific enthalpy [J/kg].
        
        su_C_p (float): Cold-side inlet pressure [Pa].
        
        su_C__m_dot (float): Cold-side mass flow rate [kg/s].
    
    **Outputs**:
    
        ex_C_h: Cold-Side Exhaust specific enthalpy at outlet [J/kg].
        
        ex_C_p: Cold-Side Exhaust pressure at outlet [Pa].
                    
        ex_C_h: Hot-Side specific enthalpy at outlet [J/kg].
        
        ex_C_p: Hot-Side pressure at outlet [Pa].
        
        Q_dot: Total heat transfer rate across the exchanger [W].
        
        DT_pinch: Minimum temperature difference between hot and cold streams across all segments [K].
        
        h_hot, h_cold: Arrays of enthalpies across discretization points [J/kg].
        T_hot, T_cold: Arrays of temperatures across discretization points [K].
    
    """
    
    def __init__(self):
        super().__init__()
        self.su_C = MassConnector() # Working fluid supply
        self.su_H = MassConnector() # Secondary fluid supply
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.Q_dot = HeatConnector()
        self.guesses = {}
        self.DT_pinch = -2
        
        self.DP_h = 0
        self.DP_c = 0
        
        self.h_hot = None
        self.h_cold = None
        self.T_hot = None
        self.T_cold = None

        self.eta_pinch = 1
        self.effectiveness = 1

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']
    
    def get_required_parameters(self):
        return [
            'eta_max', 'n_disc', 'Pinch_min' # Efficiency
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

    def pinch_residual(self, eta):
        """Return DT_pinch(eta) so we can root-find or bisect on it."""
        self.counterflow_discretized(eta)
        # we want DT_pinch = target (e.g. 0 or Pinch_min)
        return self.DT_pinch

    def find_Q_dot_max(self):
        # --- External enthalpy-based Q_dot_max ---
        self.AS_H.update(CP.PT_INPUTS, self.su_H.p, self.su_C.T)
        H_h_id = self.AS_H.hmass()
        
        self.AS_C.update(CP.PT_INPUTS, self.su_C.p, self.su_H.T)
        H_c_id = self.AS_C.hmass()
        
        Q_dot_maxh = self.su_H.m_dot * (self.su_H.h - H_h_id)
        Q_dot_maxc = self.su_C.m_dot * (H_c_id - self.su_C.h)
        
        self.Q_dot_max_ext = np.min([Q_dot_maxh, Q_dot_maxc])
        self.Q_dot_max = self.Q_dot_max_ext
        
        # Quick check: if at eta=1 the pinch is already > 0,
        # then internal limit = external limit and we can skip search.
        self.pinch_residual(1.0)
        if self.DT_pinch > 1e-3:
            self.eta_pinch = 1.0
            self.Q_dot_max_int = self.Q  # at eta=1
            self.Q_dot_max = min(self.Q_dot_max_ext, self.Q_dot_max_int)
            return
    
        # Otherwise, find eta_pinch where DT_pinch ~ 0 using bisection
        eta_low, eta_high = 0.0, 1.0
        tol_int = 1e-3
        max_it = 30
    
        for _ in range(max_it):
            eta_mid = 0.5*(eta_low + eta_high)
            self.pinch_residual(eta_mid)
    
            if self.DT_pinch > tol_int:
                # pinch OK → can increase eta
                eta_low = eta_mid
            else:
                # pinch too small → decrease eta
                eta_high = eta_mid
    
            if abs(eta_high - eta_low) < 1e-3:
                break
    
        self.eta_pinch = eta_low
        # one last evaluation at eta_pinch for consistency
        self.pinch_residual(self.eta_pinch)
        self.Q_dot_max_int = self.Q
        self.Q_dot_max = min(self.Q_dot_max_ext, self.Q_dot_max_int)
        
        return

    def counterflow_discretized(self, eta):
                
        "Heat Transfer Rate"
        n = self.params['n_disc']
        
        self.Q = eta*self.Q_dot_max
        Q_dot_seg = self.Q / n
        
        # Set inlet enthalpies
        self.h_hot[0] = self.su_H.h
        self.T_hot[0] = self.su_H.T
        
        h_cold_out = self.su_C.h + self.Q/self.su_C.m_dot
        self.h_cold[0] = h_cold_out
        
        self.AS_C.update(CP.HmassP_INPUTS, h_cold_out, self.su_C.p)
        self.T_cold[0] = self.AS_C.T()
                
        for i in range(n):
            # Hot side: forward direction
            self.h_hot[i+1] = self.h_hot[i] - Q_dot_seg / self.su_H.m_dot
            
            self.AS_H.update(CP.HmassP_INPUTS, self.h_hot[i+1], self.p_hot[i+1])
            self.T_hot[i+1] = self.AS_H.T()
                        
            # Cold side: reverse direction
            self.h_cold[i+1] = self.h_cold[i] - Q_dot_seg / self.su_C.m_dot
            
            self.AS_C.update(CP.HmassP_INPUTS, self.h_cold[i+1], self.p_cold[i+1])
            self.T_cold[i+1] = self.AS_C.T()
                    
        self.DT_pinch = np.min(self.T_hot - self.T_cold)
        
        return 

    def solve(self):
        # Ensure all required checks are performed

        if self.su_H.m_dot is None: 
            self.su_H.m_dot = self.su_C.m_dot
        
        if self.su_C.m_dot is None:
            self.su_C.m_dot = self.su_H.m_dot            

        self.check_calculable()
        self.check_parametrized()

        if self.su_H.T < self.su_C.T:
            save = self.su_C
            self.su_C = self.su_H
            self.su_H = save

        if 'DP_h' in self.params:
            self.DP_h = self.params['DP_h']
            
        if 'DP_c' in self.params:
            self.DP_c = self.params['DP_c']

        if not self.calculable:
            print("HTX IS NOT CALCULABLE")
            return

        if not self.parametrized:
            print("HTX IS NOT PARAMETRIZED")
            return

        self.AS_H = CP.AbstractState('BICUBIC&HEOS', self.su_H.fluid)
        self.AS_C = CP.AbstractState('BICUBIC&HEOS', self.su_C.fluid)

        if self.T_hot is None:
            # Allocate arrays
            self.h_hot = np.zeros(self.params['n_disc']+1)
            self.h_cold = np.zeros(self.params['n_disc']+1)
            self.T_hot = np.zeros(self.params['n_disc']+1)
            self.T_cold = np.zeros(self.params['n_disc']+1)

            self.p_hot = np.zeros(self.params['n_disc']+1)
            self.p_cold = np.zeros(self.params['n_disc']+1)
                
            # Establish pressure drops 
            DP_h_disc = self.DP_h/self.params['n_disc']
            DP_c_disc = self.DP_c/self.params['n_disc']
            
            self.p_hot[0] = self.su_H.p
            self.p_cold[0] = self.su_C.p
            
            for i in range(self.params['n_disc']):
                self.p_hot[i+1] = self.p_hot[i] - DP_h_disc
                self.p_cold[i+1] = self.p_cold[i] - DP_c_disc
            
        "Find Q_dot_max"

        self.DT_pinch = -1

        self.find_Q_dot_max()

        self.DT_pinch = -1

        self.epsilon = self.params['eta_max']

        while self.DT_pinch <= self.params['Pinch_min']:
            self.counterflow_discretized(self.epsilon)

            if self.DT_pinch <= self.params['Pinch_min']:
                self.epsilon -= 0.01
                                
                if self.epsilon <= 0:
                    self.solved = False
                    
                    if self.print_flag:
                        print("No eta satisfies Pinch_min in HXEffCstDisc")
                        
                    return

        "Outlet states"   
        
        self.update_connectors()
        self.solved = True
        return

    def update_connectors(self):
        
        "Mass Connectors"
        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)
        self.ex_C.set_h(self.su_C.h + self.Q/self.su_C.m_dot)
        self.ex_C.set_p(self.p_cold[-1])
        
        self.AS_C.update(CP.HmassP_INPUTS, self.ex_C.h, self.ex_C.p)
        self.ex_C.set_T(self.AS_C.T())

        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_m_dot(self.su_H.m_dot)
        self.ex_H.set_h(self.su_H.h - self.Q/self.su_H.m_dot)
        self.ex_H.set_p(self.p_hot[-1])
        
        self.AS_H.update(CP.HmassP_INPUTS, self.ex_H.h, self.ex_H.p)
        self.ex_H.set_T(self.AS_H.T())
        
        "Heat conector"
        self.Q_dot.set_Q_dot(self.Q)

        return

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q_dot.Q_dot}")
        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")
        print("======================")

    def plot_disc(self):
        import matplotlib.pyplot as plt

        x = np.arange(len(self.T_hot))
        plt.plot(x, self.T_hot, label='Hot Fluid [°C]', marker='o')
        plt.plot(x, self.T_cold, label='Cold Fluid [°C]', marker='s')
        plt.xlabel('Discretization')
        plt.ylabel('Temperature [°C]')
        plt.title('Counterflow Heat Exchanger Temperature Profiles')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        return



