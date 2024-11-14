from component.heat_exchanger.moving_boundary.simple_model.simulation_model import HXSimpleMB
import numpy as np


"Evaporator test"

#Exo ORC M&S
EVAP = HXSimpleMB()

EVAP.set_inputs(
    fluid_wf = 'R245fa',
    su_wf_T = 60+273.15,
    su_wf_m_dot = 0.4621,
    su_wf_x = 0,
    fluid_sf = 'Water',
    su_sf_T = 70+273.15,
    su_sf_p = 2e5,
    su_sf_m_dot = 2.425,
    ex_sf_T = 62+273.15
)

EVAP.set_parameters(**{
    'HX_type': 'evaporator',
    'HX_D': 0.06, #Diamètre de port d'entré
    'HX_A': 17.8, #Surface d'échange
    'min_pinch': 2,
    'Delta_T_sup_or_sub': 5
})

EVAP.solve()
EVAP.print_results()