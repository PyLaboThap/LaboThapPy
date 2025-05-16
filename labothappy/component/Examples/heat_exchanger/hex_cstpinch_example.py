
import __init__
from component.heat_exchanger.hex_cstpinch import HXPinchCst

# from simulation_model import HXPinchCst
import numpy as np

"Evaporator test"

# # Exo ORC M&S
EVAP = HXPinchCst()

EVAP.set_inputs(
    su_C_fluid = 'Cyclopentane',
    su_C_T = 130+273.15,
    su_C_p = 831.8*1e3,
    su_C_m_dot = 51.03,

    su_H_fluid = 'Water', #Oil
    su_H_T = 145+273.15,
    su_H_p = 5*1e5,
    su_H_m_dot = 400,
)

EVAP.set_parameters(**{
    'Pinch': 4,
    'Delta_T_sh_sc': 10,
    'type_HX': 'evaporator'
})

EVAP.solve()
EVAP.print_results()
EVAP.print_states_connectors()
EVAP.plot_disc()

# "Condenser test"

# COND = HXPinchCst()

# COND.set_inputs(
#     su_H_fluid = 'Cyclopentane',
#     su_H_T = 41.2+273.15,
#     su_H_p = 68.3*1e3,
#     su_H_m_dot = 46.18,
    
#     su_C_fluid = 'Air',
#     su_C_T = 20+273.15,
#     su_C_p = 1e5,
#     su_C_m_dot = 1911
# )

# COND.set_parameters(**{
#     'Pinch': 5,
#     'Delta_T_sh_sc': 5,
#     'type_HX': 'condenser'
# })

# COND = HXPinchCst()

# COND.set_inputs(
#     su_H_fluid = 'Water',
#     su_H_T = 130+273.15,
#     su_H_p = 5*1e5,
#     su_H_m_dot = 0.16,
    
#     su_C_fluid = 'Water',
#     su_C_T = 90+273.15,
#     su_C_p = 100e5,
#     su_C_m_dot = 1000
    
# )

# COND.set_parameters(**{
#     'Pinch': 10,
#     'Delta_T_sh_sc': 10,
#     'type_HX': 'condenser'
# })

COND.solve()
COND.print_results()
COND.print_states_connectors()
COND.plot_disc()
