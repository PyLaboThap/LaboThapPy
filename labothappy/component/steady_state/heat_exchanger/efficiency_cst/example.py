import __init__
import component.steadystate.heat_exchanger.pinch_cst
from component.steadystate.heat_exchanger.pinch_cst import simulation_model

from simulation_model import HXEffCst
import numpy as np

"Simple test - Model for recuperators in simple thermodynamic studies "

#Exo ORC M&S
HTX = HXEffCst()

HTX.set_inputs(
    su_C_fluid = 'Cyclopentane',
    su_C_T = 273.155 + 25.85,
    su_C_m_dot = 26.61,
    su_C_p = 1048 * 1e3,

    su_H_fluid = 'Cyclopentane',
    su_H_T = 273.15 + 79.86,
    su_H_m_dot = 26.61,
    su_H_p = 48.13 * 1e3,
)

HTX.set_parameters(**{
    'eta': 0.8,
})

HTX.solve()
HTX.print_results()
HTX.print_states_connectors()
