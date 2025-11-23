# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 15:00:41 2025

WorkShop of Basile

@author: titouanjanod
"""

# import __init__
# import numpy as np
from connector.mass_connector import MassConnector
# from connector.work_connector import WorkConnector
from machine.circuit_rec import RecursiveCircuit
from component.compressor.compressor_csteff import CompressorCstEff
from component.heat_exchanger.hex_cstpinch import HexCstPinch
from component.valve.valve_isenthalpic import ValveIsenthalpic

# from CoolProp.CoolProp import PropsSI

fluid = 'propane'
sf = 'water'
HP = RecursiveCircuit(fluid)

# Guesses
mdot = 0.2
SC_cd = 3
SH_ev = 3
p_low = 7e5
p_high = 15e5

CP = CompressorCstEff()
CD = HexCstPinch()
EV = HexCstPinch()
VV = ValveIsenthalpic()

CD_sink = MassConnector()
EV_source = MassConnector()


CP.set_parameters(eta_is=0.7)
CD.set_parameters(Pinch=3 , Delta_T_sh_sc=3, HX_type='condenser')
EV.set_parameters(Pinch=3 , Delta_T_sh_sc=3, HX_type='evaporator')
VV.set_parameters()


HP.add_component(CP, "CP")
HP.add_component(CD, "CD")
HP.add_component(EV, "EV")
HP.add_component(VV, "VV")

HP.link_components("CP", "m-ex", "CD", "m-su_H")
HP.link_components("CD", "m-ex_H", "VV", "m-su")
HP.link_components("VV", "m-ex", "EV", "m-su_C")
HP.link_components("EV", "m-ex_C", "CP", "m-su")

HP.add_source("CD_water", CD_sink, HP.components["CD"], "m-su_C")
HP.set_source_properties(T=40+273.15, fluid=sf, m_dot=0.5, P=3e5, target="CD_water")


HP.add_source("EV_water", EV_source, HP.components["EV"], "m-su_C")
HP.set_source_properties(T=20+273.15, fluid=sf, m_dot=5, P=3e5, target="EV_water")

HP.set_cycle_guess(target='CP:su', m_dot=mdot, SH=SH_ev, p=p_low)
HP.set_cycle_guess(target='CP:ex', p=p_high)

# HP.set_cycle_guess(target='CD:ex_H', p=p_high, m_dot=mdot, SC=SC_cd)
# HP.set_cycle_guess(target='VV:ex', p=p_low)


HP.set_cycle_guess(target='VV:su', p=p_high, m_dot=mdot, SC=SC_cd)
HP.set_cycle_guess(target='VV:ex', p=p_low)



HP.set_residual_variable(target="VV:ex", variable='h', tolerance=1e-3)
HP.set_residual_variable(target="VV:ex", variable='p', tolerance=1e-3)


HP.set_residual_variable(target="EV:ex_C", variable='h', tolerance=1e-3)
HP.set_residual_variable(target="EV:ex_C", variable='p', tolerance=1e-3)


HP.set_residual_variable(target="CD:ex_H", variable='h', tolerance=1e-3)
HP.set_residual_variable(target="CD:ex_H", variable='p', tolerance=1e-3)

HP.set_residual_variable(target="CP:ex", variable='h', tolerance=1e-3)
HP.set_residual_variable(target="CP:ex", variable='p', tolerance=1e-3)



HP.solve()




