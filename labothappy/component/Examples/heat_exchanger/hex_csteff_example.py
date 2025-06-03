# import __init__
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

###

from toolbox.heat_exchanger.hex_csteff import HXEffCst

import numpy as np

"Simple test - Model for recuperators in simple thermodynamic studies "

#Exo ORC M&S
HTX = HXEffCst()

# Set input conditions
HTX.set_inputs(
    Csu_fluid='Water',
    Csu_T=273.15 + 60,
    Csu_m_dot=1,
    Csu_p=3e5,

    Hsu_fluid='CO2',
    Hsu_T=273.15 + 160,
    Hsu_m_dot=0.08,
    Hsu_p=130e5,
)

# Set heat exchanger parameters
HTX.set_parameters(**{
    'eta': 0.95,
})

# Solve the model
HTX.solve()

# Output results
HTX.print_results()
HTX.print_states_connectors()