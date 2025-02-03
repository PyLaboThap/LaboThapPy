# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

from machine.circuit import Circuit
from CoolProp.CoolProp import PropsSI

from connector.mass_connector import MassConnector

from component.heat_exchanger.steady_state.pinch_cst.simulation_model import HXPinchCst
from component.heat_exchanger.steady_state.cst_efficiency.simulation_model import HXEffCst
from component.volumetric_machine.compressor.steady_state.constant_isentropic_efficiency.simulation_model import CompressorCstEff 
from component.valve.isenthalpic_valve import Isenthalpic_Valve


CO2_HP = Circuit('CO2')

# Create components
Compressor = CompressorCstEff()
GasCooler = HXEffCst()
Valve = Isenthalpic_Valve()
Evaporator = HXPinchCst()

#%% COMPRESSOR PARAMETERS

Compressor.set_parameters(eta_is=0.8)

#%% GASCOOLER PARAMETERS

GasCooler.set_parameters(**{
    'eta': 0.95,
})

#%% EVAPORATOR PARAMETERS

Evaporator.set_parameters(**{
    'Pinch': 3,
    'Delta_T_sh_sc': 3,
    'type_HX': 'evaporator'
})

#%% ADD AND LINK COMPONENTS

CO2_HP.add_component(Compressor, "Compressor")
CO2_HP.add_component(GasCooler, "GasCooler")
CO2_HP.add_component(Valve, "Valve")
CO2_HP.add_component(Evaporator, "Evaporator")

# # Link components+
# ORC.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
# ORC.link_components("Evaporator", "m-ex_C", "Spliter", "m-su")

# ORC.link_components("Spliter", "m-ex_1", "Expander_1", "m-su")
# ORC.link_components("Spliter", "m-ex_2", "Expander_2", "m-su")
# ORC.link_components("Spliter", "m-ex_3", "Expander_3", "m-su")

# ORC.link_components("Expander_1", "m-ex", "Mixer", "m-su_1")
# ORC.link_components("Expander_2", "m-ex", "Mixer", "m-su_2")
# ORC.link_components("Expander_3", "m-ex", "Mixer", "m-su_3")

# ORC.link_components("Mixer", "m-ex", "Condenser", "m-su_H")
# ORC.link_components("Condenser", "m-ex_H", "Pump", "m-su")