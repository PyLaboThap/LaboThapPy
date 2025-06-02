"""
Created on Aug 03 21:31:37 2023

@author: 
"""
import os
import sys

def find_project_root(starting_dir):
    markers = ['connector', 'component']
    current_dir = starting_dir
    while True:
        if all(os.path.isdir(os.path.join(current_dir, marker)) for marker in markers):
            return current_dir
        parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
        if parent_dir == current_dir:
            return None
        current_dir = parent_dir

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = find_project_root(current_dir)

if project_root and project_root not in sys.path:
    sys.path.insert(0, project_root)
    
###########################

from toolbox.pump.pump_csteff import PumpCstEff
import numpy as np

# Example usage
PP = PumpCstEff()

# Set initial conditions
PP.su.set_properties(P=319296.56, T=331.03, fluid='R1233ZDE')
PP.su.set_m_dot(1.0)  
PP.ex.set_properties(P=606240.14, fluid='R1233ZDE')
PP.set_parameters(eta_is=0.9)
PP.solve()
PP.print_results()

# Second case
PP.su.set_properties(P=400000, fluid='R1233ZDE')
PP.su.set_m_dot(1.0)  
PP.ex.set_properties(P=800000, fluid='R1233ZDE')
PP.solve()
PP.print_results()

# Third case (set T via inputs)
PP.set_inputs(su_T=400)
PP.su.set_m_dot(1.0)  
PP.solve()
PP.print_states_connectors()


