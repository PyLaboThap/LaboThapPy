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
Main inputs that could change
"""
N_rot_pp = 3000


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

Geo_brine = MassConnector()
ORC.add_source("geo_brine", Geo_brine, ORC.components["PH"], "m-su_H")

Cooling_air = MassConnector()
ORC.add_source("Cooling_air", Cooling_air, ORC.components["ACC"], "m-su_C")

Eq_water = MassConnector()
ORC.add_source("Eq_water", Eq_water, ORC.components["EV"], "m-su_H")



ORC.link_components("PH", "", "Spliter", "m-su")

ORC.link_components("TURB", "m-ex", "Mixer", "m-su_1")
ORC.link_components("Expander_2", "m-ex", "Mixer", "m-su_2")
ORC.link_components("Expander_3", "m-ex", "Mixer", "m-su_3")

ORC.link_components("Mixer", "m-ex", "Condenser", "m-su_H")
ORC.link_components("Condenser", "m-ex_H", "Pump", "m-su")


if __name__ == "__main__":
    
    d=1


