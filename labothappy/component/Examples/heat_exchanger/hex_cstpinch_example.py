
import __init__
from component.heat_exchanger.hex_cstpinch import HXPinchCst

# from simulation_model import HXPinchCst
import numpy as np

"Evaporator test"

# -------------------------------------------------------------------------------------------------------------

# Exo ORC M&S
# EVAP = HXPinchCst()

# EVAP.set_inputs(
#     fluid_C = 'Cyclopentane',
#     T_su_C = 130+273.15,
#     P_su_C = 831.8*1e3,
#     m_dot_C = 51.03,

#     fluid_H = 'Water', #Oil
#     T_su_H = 145+273.15,
#     P_su_H = 5*1e5,
#     m_dot_H = 400,
# )

# EVAP.set_parameters(**{
#     'Pinch': 4,
#     'Delta_T_sh_sc': 10,
#     'type_HX': 'evaporator'
# })

# EVAP.solve()
# EVAP.print_results()
# EVAP.print_states_connectors()
# EVAP.plot_disc()

# -------------------------------------------------------------------------------------------------------------

"Condenser test"

# COND = HXPinchCst()

# COND.set_inputs(
#     fluid_H = 'Cyclopentane',
#     T_su_H = 41.2+273.15,
#     P_su_H = 68.3*1e3,
#     m_dot_H = 46.18,
    
#     fluid_C = 'Air',
#     T_su_C = 20+273.15,
#     P_su_C = 1e5,
#     m_dot_C = 1911
# )

# COND.set_parameters(**{
#     'Pinch': 5,
#     'Delta_T_sh_sc': 5,
#     'type_HX': 'condenser'
# })

# -------------------------------------------------------------------------------------------------------------

COND = HXPinchCst()

COND.set_inputs(
    fluid_H = 'Water',
    T_su_H = 130+273.15,
    P_su_H = 5*1e5,
    m_dot_H = 0.16,
    
    fluid_C = 'Water',
    T_su_C = 90+273.15,
    P_su_C = 100e5,
    m_dot_C = 1000
    
)

COND.set_parameters(**{
    'Pinch': 10,
    'Delta_T_sh_sc': 10,
    'type_HX': 'condenser',
    'DP_c': 10*1e3,
    'DP_h': 10*1e3
})

# -------------------------------------------------------------------------------------------------------------

COND.solve()
COND.print_results()
COND.print_states_connectors()
COND.plot_disc()
