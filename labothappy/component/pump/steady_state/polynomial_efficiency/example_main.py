# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:07:43 2024

@author: Basile
"""
import __init__

from connector.mass_connector import MassConnector
from component.steady_state.pump.polynomial_efficiency.simulation_model import PumpPolynEff
from component.steady_state.pump.polynomial_efficiency.modules.c_pump_polyn_geom import GeometryPolynPump
from CoolProp.CoolProp import PropsSI    

"--------- 1) EXAMPLE : DECAGONE PUMP ------------------------------------------------------------------------------------------"

"Create pump object"
Pump = PumpPolynEff()

"Set Inputs"
Pump.set_inputs(
                su_fluid = 'Cyclopentane',
                su_T = 38 + 273.15, # K
                su_p = 1*1e5, # Pa
                ex_p = 31.5*1e5 # Pa
                )

"Geometry"
Pump_geom = GeometryPolynPump()
Pump_geom.set_parameters("DECAGONE_PUMP") 

"Set Parameters"
Pump.set_parameters(
                    min_flowrate = Pump_geom.min_flowrate, rated_flowrate = Pump_geom.rated_flowrate, max_flowrate = Pump_geom.max_flowrate,
                    N_pp_rated = Pump_geom.N_pp_rated, pump_voltage = Pump_geom.pump_voltage, pump_phases = Pump_geom.pump_phases, 
                    eta_v = Pump_geom.eta_v, V_swept = Pump_geom.V_swept, eta_max_motor = Pump_geom.eta_max_motor, W_dot_el_rated = Pump_geom.W_dot_el_rated,
                    coefs_pump = Pump_geom.coefs_pump, eta_tot = Pump_geom.eta_tot, eta_m = Pump_geom.eta_m, 

                    N_pp = 50
                    )


Pump.solve()
