# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 08:39:27 2026

@author: titouanjanod
"""

from connector.heat_connector import HeatConnector
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from machine.base_circuit import BaseCircuit
from machine.circuit_rec import RecursiveCircuit

from toolbox.geometries.heat_exchanger.c_geometry_HXs_Zorlu import Zorlu_HXs
from component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from component.heat_exchanger.hex_eNTU import HexeNTU

from component.base_component import BaseComponent
from component.heat_exchanger.hex_cstpinch import HexCstPinch
from component.heat_exchanger.hex_csteff import HexCstEff

from component.pump.pump_csteff import PumpCstEff
from component.expander.expander_csteff import ExpanderCstEff

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import AbstractState
import CoolProp.CoolProp as CP



"""
Main inputs that could change ; for initialisation
"""
N_rot_pp = 3000 # RPM
T_brine = 104.2 # °C
p_brine = 1.59e5 # Pa
fluid_wh = "Water" # to be changed by implementing thing of Polimi/Turboden
m_dot_brine_max = 857.9 # kg/s
WH_use_pre = 0.2 # 20%

T_amb = 15 # °C
p_amb = 99.3e3 # Pa
Glide_air = 10 # K
fluid_air = 'Air'
m_dot_air = 2557 # kg/s

m_dot_eqwater = 1000 # kg/s
T_PCM = 134.8 # °C
fluid_eq = 'Water'
p_eqwater = 8e5 # Pa 



ORC = RecursiveCircuit('Cyclopentane')


# importing components geometries
geom_obj = Zorlu_HXs() # note that the geometries are in an assynchronised file, since it contains sensitive data


"PREHEATER instanciation"
EV = HexCstPinch()
# CD = hex_cstpinch()  
EV.set_parameters(**{
    'Pinch': 10,
    'Delta_T_sh_sc': 3 })


PH = HexMBChargeSensitive('Shell&Tube')
# PH.set_inputs(
#     fluid_H = 'Water',
#     T_su_H = 104.2 + 273.15, # K
#     P_su_H = 2.1*1e5, # Pa
#     m_dot_H = 171.6, # kg/s

#     # Second fluid
#     fluid_C = 'Cyclopentane',
#     T_su_C = 50.1 + 273.15, # K
#     P_su_C = 7.75*1e5, # Pa
#     m_dot_C = 62.74, # kg/s  
#     )
Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
geom_obj.set_parameters("ORC_preheater")
PH.set_parameters(**geom_obj.geom)
PH.set_parameters(n_disc=30)
PH.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31
PH.set_DP(**geom_obj.DP) 

# PH.solve()
# PH.print_states_connectors()



"""
Recuperator instanciation
"""
REC = HexeNTU("Tube&Fins")

# REC.set_inputs(
#     fluid_H = 'Cyclopentane',
#     T_su_H = 64.4 + 273.15, # K
#     P_su_H = 0.8*1e5, # Pa
#     m_dot_H = 62.71, # kg/s

#     # Second fluid
#     fluid_C = 'Cyclopentane',
#     T_su_C = 31.27 + 273.15, # K
#     P_su_C = 8.258*1e5, # Pa
#     m_dot_C = 62.74, # kg/s  
#     )
Corr_C = "Gnielinski"
Corr_H = "Tube_And_Fins"
REC.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'

geom_obj = Zorlu_HXs()
geom_obj.set_parameters("ORC_recuperator")
REC.set_parameters(**geom_obj.geom)
geom_obj.set_DP("ORC_recuperator")
REC.set_DP(**geom_obj.DP)

# REC.solve()
# REC.print_results()
# print(f"\n  - Q_dot = {REC.Q_hex.Q_dot/1000} kW")


"""
ACC instanciation
"""
ACC = HexeNTU("Tube&Fins")

# ACC.set_inputs(
#     fluid_H = 'Cyclopentane',
#     T_su_H = 38.46 + 273.15, # K
#     P_su_H = 0.7003*1e5, # Pa
#     m_dot_H = 62.71, # kg/s

#     # Second fluid
#     fluid_C = 'Air',
#     T_su_C = 15 + 273.15, # K
#     P_su_C = 1.01325*1e5, # Pa
#     m_dot_C = 2557, # kg/s  
#     )

Corr_C = "Tube_And_Fins"
Corr_H = "Gnielinski"
ACC.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'

geom_obj = Zorlu_HXs()
geom_obj.set_parameters("ORC_ACC")
ACC.set_parameters(**geom_obj.geom)
geom_obj.set_DP("ORC_ACC")
# ACC.set_DP()
ACC.set_DP(**geom_obj.DP)

# ACC.solve()
# ACC.print_results()
# print(f"\n  - Q_dot = {ACC.Q_hex.Q_dot/1000} kW")


PP = PumpCstEff()
PP.set_parameters(eta_is=0.75)
PP.set_inputs(N_rot = N_rot_pp)

TURB = ExpanderCstEff()
TURB.set_parameters(eta_is=0.86)

ORC.add_component(EV, "EV")
ORC.add_component(PH, "PH")
ORC.add_component(REC, "REC")
ORC.add_component(ACC, "ACC")
ORC.add_component(PP, "PP")
ORC.add_component(TURB, "TURB")


# ORC.link_components("PP", "", component2_name, input_connector)
ORC.link_components("PP", "m-ex", "REC", "m-su_C")
ORC.link_components("REC", "m-ex_C", "PH", "m-su_C")
ORC.link_components("PH", "m-ex_C", "EV", "m-su_C")
ORC.link_components("EV", "m-ex_C", "TURB", "m-su")
ORC.link_components("TURB", "m-ex", "REC", "m-su_H")
ORC.link_components("REC", "m-ex_H", "ACC", "m-su_H")
ORC.link_components("ACC", "m-ex_H", "PP", "m-su")


m_dot_brine = m_dot_brine_max * WH_use_pre
Geo_brine = MassConnector()
ORC.add_source("geo_brine", Geo_brine, ORC.components["PH"], "m-su_H")
ORC.set_source_properties(T= T_brine, fluid= fluid_wh, m_dot= m_dot_brine, target='geo_brine', P = p_brine)


Cooling_air = MassConnector()
ORC.add_source("Cooling_air", Cooling_air, ORC.components["ACC"], "m-su_C")
ORC.set_source_properties(T = T_amb, fluid = fluid_air, m_dot = m_dot_air, target='Cooling_air', P = p_amb)
# ORC.set_source_properties(T= , fluid= , m_dot= , target='', P = )

Eq_water = MassConnector()
ORC.add_source("Eq_water", Eq_water, ORC.components["EV"], "m-su_H")
ORC.set_source_properties(T = T_PCM, fluid= fluid_eq , m_dot= m_dot_eqwater, target='Eq_water', P = p_eqwater)


p_low = 0.7*1e5 
p_high = 7.66*1e5

ORC.set_cycle_guess(target='PP:su', SC = 1, p = p_low, m_dot = 62.74)
ORC.set_cycle_guess(target='PP:ex', p = p_high)

ORC.set_fixed_properties(target='PP:su', SC = 1)

ORC.set_iteration_variable(target=['ACC:ex_H'], variable='p', objective = 'PP:su-SC', tol = 1e-2, rel = 1, damping_factor = 0.2)


ORC.set_residual_variable(target = "EV:ex_C" , variable='h', tol = 1e-2)
ORC.set_residual_variable(target = "CD:ex_H" , variable='h', tol = 1e-2)
ORC.set_residual_variable(target = "EV:ex_C" , variable='m_dot', tol = 1e-1)
ORC.set_residual_variable(target = "CD:ex_H" , variable='m_dot', tol = 1e-1)



if __name__ == "__main__":
    d=1
    # ORC.solve()


