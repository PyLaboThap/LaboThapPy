
import __init__
from component.heat_exchanger.steady_state.pinch_cst.hex_cstpinch import HXPinchCst

# from simulation_model import HXPinchCst
import numpy as np

"Evaporator test"

#Exo ORC M&S
EVAP = HXPinchCst()

EVAP.set_inputs(
    su_C_fluid = 'CO2',
    su_C_h = 200000,
    su_C_m_dot = 0.08,
    su_H_fluid = 'Water', #Oil
    su_H_T = 15+273.15,
    su_H_cp = 4188,
    su_H_m_dot = 1,
)

EVAP.set_parameters(**{
    'Pinch': 3,
    'Delta_T_sh_sc': 3,
    'type_HX': 'evaporator'
})

EVAP.solve()
EVAP.print_results()
EVAP.print_states_connectors()

# # "Condenser test"

# COND = HXPinchCst()

# COND.set_inputs(
#     fluid_wf = 'R245fa',
#     su_wf_h = 512299,
#     su_wf_m_dot = 0.06,
#     fluid_sf = 'Water',
#     su_sf_T = 30+273.15,
#     su_sf_cp = 4187,
#     su_sf_m_dot = 0.4
# )

# COND.set_parameters(**{
#     'Pinch': 5,
#     'Delta_T_sh_sc': 5,
#     'type_HX': 'condenser'
# })

# COND.solve()
# COND.print_results()
# COND.print_states_connectors()
