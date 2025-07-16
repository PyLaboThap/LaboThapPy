import __init__
import os
import sys

from component.heat_exchanger.hex_csteff import HXEffCst

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

    Hsu_fluid='Water',
    Hsu_T=273.15 + 110,
    Hsu_m_dot=0.3,
    Hsu_p=10e5,
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
