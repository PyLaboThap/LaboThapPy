# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:07:43 2024

@author: Basile
"""
import __init__

from connector.mass_connector import MassConnector
from component.steady_state.turbomachinery.pump.extrapolation_model.simulation_model import PumpExtrapolationModel
from component.steady_state.turbomachinery.pump.extrapolation_model.modules.geometry_extrapolation_pump import Geometry_extrapol_pump
from CoolProp.CoolProp import PropsSI    

"--------- 1) DECAGONE PUMP ------------------------------------------------------------------------------------------"
# Instantiate Pump
Pump = PumpExtrapolationModel()

# Set Input
Pump.set_inputs(
                su_fluid = 'Cyclopentane',
                su_T = 32+273.15,
                su_p = 1*1e5,
                ex_p = 31.5*1e5,
                Omega_pp = 3000
                )

"Geometry"

Pump_geom = Geometry_extrapol_pump()
Pump_geom.set_parameters("DECAGONE_PUMP")

"Parameters"
Pump.set_parameters(
                    Omega_rated = Pump_geom.Omega_rated, min_flowrate = Pump_geom.min_flowrate, rated_flowrate = Pump_geom.rated_flowrate, max_flowrate = Pump_geom.max_flowrate,
                    PI_rated = Pump_geom.PI_rated, D_p = Pump_geom.D_p, V_dot_curve = Pump_geom.V_dot_curve, Delta_H_curve = Pump_geom.Delta_H_curve, eta_is_curve = Pump_geom.eta_is_curve,
                    NPSH_r_curve = Pump_geom.NPSH_r_curve, eta_m = Pump_geom.eta_m, eta_max_motor = Pump_geom.eta_max_motor, W_dot_el_rated = Pump_geom.W_dot_el_rated
)

Pump.solve()
print(Pump.W_dot_wf)
