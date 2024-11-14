
import __init__
from component.steadystate.heat_exchanger.pinch_cst import simulation_model

from simulation_model import HXPinchCst
import numpy as np

"Evaporator test"

#Exo ORC M&S
EVAP = HXPinchCst()

EVAP.set_inputs(
    fluid_wf = 'R245fa',
    su_wf_h = 249154,
    su_wf_m_dot = 0.0638,
    fluid_sf = 'INCOMP::PNF2', #Oil
    su_sf_T = 150+273.15,
    su_sf_cp = 2000,
    su_sf_m_dot = 0.363,
)

EVAP.set_parameters(**{
    'Pinch': 5,
    'Delta_T_sh_sc': 5,
    'type_HX': 'evaporator'
})

EVAP.solve()
EVAP.print_results()
EVAP.print_states_connectors()

# "Condenser test"

COND = HXPinchCst()

COND.set_inputs(
    fluid_wf = 'R245fa',
    su_wf_h = 512299,
    su_wf_m_dot = 0.06,
    fluid_sf = 'Water',
    su_sf_T = 30+273.15,
    su_sf_cp = 4187,
    su_sf_m_dot = 0.4
)

COND.set_parameters(**{
    'Pinch': 5,
    'Delta_T_sh_sc': 5,
    'type_HX': 'condenser'
})

COND.solve()
COND.print_results()
COND.print_states_connectors()
