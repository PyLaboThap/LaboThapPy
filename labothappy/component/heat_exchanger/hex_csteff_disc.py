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
import numpy as np
import math

class HXEffCstDisc(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su_C = MassConnector() # Working fluid supply
        self.su_H = MassConnector() # Secondary fluid supply
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.Q_dot = HeatConnector()
        self.guesses = {}
        self.DT_pinch = -1

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return['Csu_fluid', 'Csu_h', 'Csu_p', 'Csu_m_dot', 'Hsu_fluid', 'Hsu_h', 'Hsu_p', 'Hsu_m_dot']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su_C.fluid is not None:
            self.inputs['Csu_fluid'] = self.su_C.fluid
        if self.su_C.h is not None:
            self.inputs['Csu_h'] = self.su_C.h
        if self.su_C.T is not None:
            self.inputs['Csu_T'] = self.su_C.T
        if self.su_C.m_dot is not None:
            self.inputs['Csu_m_dot'] = self.su_C.m_dot
        if self.su_C.p is not None:
            self.inputs['Csu_p'] = self.su_C.p

        if self.su_H.fluid is not None:
            self.inputs['Hsu_fluid'] = self.su_H.fluid
        if self.su_H.T is not None:
            self.inputs['Hsu_T'] = self.su_H.T
        if self.su_H.h is not None:
            self.inputs['Hsu_h'] = self.su_H.h
        if self.su_H.cp is not None:
            self.inputs['Hsu_cp'] = self.su_H.cp
        if self.su_H.m_dot is not None:
            self.inputs['Hsu_m_dot'] = self.su_H.m_dot
        if self.su_H.p is not None:
            self.inputs['Hsu_p'] = self.su_H.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'Csu_fluid' in self.inputs:
            self.su_C.set_fluid(self.inputs['Csu_fluid'])
        if 'Csu_T' in self.inputs:
            self.su_C.set_T(self.inputs['Csu_T'])
        if 'Csu_h' in self.inputs:
            self.su_C.set_h(self.inputs['Csu_h'])
        if 'Csu_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['Csu_m_dot'])
        if 'Csu_p' in self.inputs:
            self.su_C.set_p(self.inputs['Csu_p'])

        if 'Hsu_fluid' in self.inputs:
            self.su_H.set_fluid(self.inputs['Hsu_fluid'])
        if 'Hsu_T' in self.inputs:
            self.su_H.set_T(self.inputs['Hsu_T'])
        if 'Hsu_h' in self.inputs:
            self.su_H.set_h(self.inputs['Hsu_h'])
        if 'Hsu_cp' in self.inputs:
            self.su_H.set_cp(self.inputs['Hsu_cp'])
        if 'Hsu_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['Hsu_m_dot'])
        if 'Hsu_p' in self.inputs:
            self.su_H.set_p(self.inputs['Hsu_p'])

        return['Csu_fluid', 'Csu_h', 'Csu_p', 'Csu_m_dot', 'Hsu_fluid', 'Hsu_h', 'Hsu_p', 'Hsu_m_dot']
    
    def get_required_parameters(self):
        return [
            'eta', 'n_disc', 'Pinch_min' # Efficiency
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

    def counterflow_discretized(self):
        
        "Define Q_dot_max through enthalpies"
        
        H_h_id = PropsSI('H', 'P', self.su_H.p, 'T', self.su_C.T, self.su_H.fluid)
        H_c_id = PropsSI('H', 'P', self.su_C.p, 'T', self.su_H.T, self.su_C.fluid)
        
        Q_dot_maxh = self.su_H.m_dot*(self.su_H.h- H_h_id)
        Q_dot_maxc = self.su_C.m_dot*(H_c_id-self.su_C.h)
        
        # print(f"self.su_C.T : {self.su_C.T}")
        # print(f"self.su_H.T : {self.su_H.T}")
        
        # print(f"Q_dot_maxh : {Q_dot_maxh}")
        # print(f"Q_dot_maxc : {Q_dot_maxc}")
        
        # print("------------------------------------------------")
        
        self.Q_dot_max = min(Q_dot_maxh,Q_dot_maxc)
                    
        "Heat Transfer Rate"
        self.Q = self.params['eta']*self.Q_dot_max    
        Q_dot_seg = self.Q / self.params['n_disc']
    
        # Set inlet enthalpies
        self.h_hot[0] = self.su_H.h
        self.T_hot[0] = self.su_H.T
        
        h_cold_out = self.su_C.h + self.Q/self.su_C.m_dot
        self.h_cold[0] = h_cold_out
        
        self.T_cold[0] = PropsSI('T', 'H', h_cold_out, 'P', self.su_C.p, self.su_C.fluid)
    
        for i in range(self.params['n_disc']):
            # Hot side: forward direction
            self.h_hot[i+1] = self.h_hot[i] - Q_dot_seg / self.su_H.m_dot
            self.T_hot[i+1] = PropsSI('T', 'H', self.h_hot[i+1], 'P', self.su_H.p, self.su_H.fluid)
    
            # Cold side: reverse direction
            self.h_cold[i+1] = self.h_cold[i] - Q_dot_seg / self.su_C.m_dot
            self.T_cold[i+1] = PropsSI('T', 'H', self.h_cold[i+1], 'P', self.su_C.p, self.su_C.fluid)
        
        self.DT_pinch = min(self.T_hot - self.T_cold)
        
        return 

    def solve(self):
        # Ensure all required checks are performed

        if self.su_H.m_dot is None: 
            self.su_H.m_dot = self.su_C.m_dot
        
        if self.su_C.m_dot is None:
            self.su_C.m_dot = self.su_H.m_dot            

        self.check_calculable()
        self.check_parametrized()

        if not self.calculable:
            print("HTX IS NOT CALCULABLE")
            return

        if not self.parametrized:
            print("HTX IS NOT PARAMETRIZED")
            return

        # Allocate arrays
        self.h_hot = np.zeros(self.params['n_disc']+1)
        self.h_cold = np.zeros(self.params['n_disc']+1)
        self.T_hot = np.zeros(self.params['n_disc']+1)
        self.T_cold = np.zeros(self.params['n_disc']+1)

        self.DT_pinch = -1 

        while self.DT_pinch <= self.params['Pinch_min']:
            self.counterflow_discretized()

            if self.DT_pinch <= self.params['Pinch_min']:
                self.params['eta'] -= 0.01
                
                if self.params['eta'] <= 0:
                    self.solved = False
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
        self.ex_C.set_p(self.su_C.p)
        self.ex_C.set_T(PropsSI("T", "H", self.ex_C.h, "P", self.ex_C.p, self.su_C.fluid))

        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_m_dot(self.su_H.m_dot)
        self.ex_H.set_h(self.su_H.h - self.Q/self.su_H.m_dot)
        self.ex_H.set_p(self.su_H.p)
        self.ex_H.set_T(PropsSI("T", "H", self.ex_H.h, "P", self.ex_H.p, self.su_H.fluid))

        "Heat conector"
        self.Q_dot.set_Q_dot(self.Q)

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



