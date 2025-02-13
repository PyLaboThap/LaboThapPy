import __init__

from component.heat_exchanger.steady_state.cst_efficiency.simulation_model import HXEffCst

import numpy as np

"Simple test - Model for recuperators in simple thermodynamic studies "

#Exo ORC M&S
HTX = HXEffCst()

HTX.set_inputs(
    su_C_fluid = 'Water',
    su_C_T = 273.155 + 60,
    su_C_m_dot = 1,
    su_C_p = 3e5,

    su_H_fluid = 'CO2',
    su_H_T = 273.15 + 160,
    su_H_m_dot = 0.08,
    su_H_p = 130*1e5,
)

HTX.set_parameters(**{
    'eta': 0.95,
})

HTX.solve()
HTX.print_results()
HTX.print_states_connectors()
