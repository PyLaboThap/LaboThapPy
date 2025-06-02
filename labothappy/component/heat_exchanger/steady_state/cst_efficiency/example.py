import __init__

from component.heat_exchanger.steady_state.cst_efficiency.simulation_model import HXEffCst

import numpy as np

"Simple test - Model for recuperators in simple thermodynamic studies "

#Exo ORC M&S
HX = HXEffCst()

HX.set_inputs(
    Csu_fluid = 'Water',
    Csu_T = 273.155 + 15,
    Csu_m_dot = 0.1,
    Csu_p = 5e5,

    Hsu_fluid = 'CO2',
    Hsu_T = 450,
    Hsu_m_dot = 0.16,
    Hsu_p = 150*1e5,
)

HX.set_parameters(**{
    'eta': 0.92,
})

HX.solve()
HX.print_results()
HX.print_states_connectors()
