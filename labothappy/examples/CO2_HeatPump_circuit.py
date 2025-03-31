# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

from CoolProp.CoolProp import PropsSI

from component.heat_exchanger.steady_state.cst_efficiency.simulation_model import (
    HXEffCst,
)
from component.heat_exchanger.steady_state.pinch_cst.simulation_model import HXPinchCst
from component.valve.isenthalpic_valve import Isenthalpic_Valve
from component.volumetric_machine.compressor.steady_state.constant_isentropic_efficiency.simulation_model import (
    CompressorCstEff,
)
from connector.mass_connector import MassConnector
from machine.circuit import Circuit

case_study = "IHX"

if case_study == "Simple":
    CO2_HP = Circuit("CO2")

    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCst()
    Valve = Isenthalpic_Valve()
    Evaporator = HXPinchCst()

    #  COMPRESSOR PARAMETERS

    Compressor.set_parameters(eta_is=0.8)

    # GASCOOLER PARAMETERS

    GasCooler.set_parameters(
        **{
            "eta": 0.95,
        }
    )

    #  EVAPORATOR PARAMETERS

    Evaporator.set_parameters(
        **{"Pinch": 3, "Delta_T_sh_sc": 3, "type_HX": "evaporator"}
    )

    #  ADD AND LINK COMPONENTS

    # Add components
    CO2_HP.add_component(Compressor, "Compressor")
    CO2_HP.add_component(GasCooler, "GasCooler")
    CO2_HP.add_component(Valve, "Valve")
    CO2_HP.add_component(Evaporator, "Evaporator")

    # Link components
    CO2_HP.link_components("Compressor", "m-ex", "GasCooler", "m-su_H")
    CO2_HP.link_components("GasCooler", "m-ex_H", "Valve", "m-su")
    CO2_HP.link_components("Valve", "m-ex", "Evaporator", "m-su_C")
    CO2_HP.link_components("Evaporator", "m-ex_C", "Compressor", "m-su")

    # SOURCES AND SINKS

    Gas_cooler_source = MassConnector()
    CO2_HP.add_source(
        "GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C"
    )
    CO2_HP.set_source_properties(
        T=60 + 273.15, fluid="Water", m_dot=1, target="GC_Water", P=3e5
    )

    EV_source = MassConnector()
    CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
    CO2_HP.set_source_properties(
        T=15 + 273.15, fluid="Water", m_dot=1, target="EV_Water", P=2e5
    )

    # CYCLE GUESSES

    P_low = 40 * 1e5
    P_high = 160 * 1e5

    CO2_HP.set_cycle_guess(target="Compressor:su", m_dot=0.08, SH=5, p=P_low)
    CO2_HP.set_cycle_guess(target="Compressor:ex", p=P_high)

    CO2_HP.set_cycle_guess(target="Valve:ex", p=P_low)
    # CO2_HP.set_cycle_guess(target='Evaporator:ex_C', p = P_low, T = 273.15 + 10, m_dot = 0.08)

    # CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target="Evaporator:ex_C", variable="h", tolerance=1e-3)
    CO2_HP.set_residual_variable(target="GasCooler:ex_H", variable="h", tolerance=1e-3)

    CO2_HP.solve()

if case_study == "IHX":

    CO2_HP = Circuit("CO2")

    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCst()
    IHX = HXEffCst()
    Valve = Isenthalpic_Valve()
    Evaporator = HXPinchCst()

    # COMPRESSOR PARAMETERS

    Compressor.set_parameters(eta_is=0.8)

    # GASCOOLER PARAMETERS

    GasCooler.set_parameters(
        **{
            "eta": 0.95,
        }
    )

    # IHX PARAMETERS

    IHX.set_parameters(
        **{
            "eta": 0,
        }
    )

    # EVAPORATOR PARAMETERS

    Evaporator.set_parameters(
        **{"Pinch": 3, "Delta_T_sh_sc": 5, "type_HX": "evaporator"}
    )

    # ADD AND LINK COMPONENTS

    # Add components
    CO2_HP.add_component(Compressor, "Compressor")
    CO2_HP.add_component(GasCooler, "GasCooler")
    CO2_HP.add_component(IHX, "IHX")
    CO2_HP.add_component(Valve, "Valve")
    CO2_HP.add_component(Evaporator, "Evaporator")

    # Link components
    CO2_HP.link_components("Compressor", "m-ex", "GasCooler", "m-su_H")
    CO2_HP.link_components("GasCooler", "m-ex_H", "IHX", "m-su_H")
    CO2_HP.link_components("IHX", "m-ex_H", "Valve", "m-su")
    CO2_HP.link_components("Valve", "m-ex", "Evaporator", "m-su_C")
    CO2_HP.link_components("Evaporator", "m-ex_C", "IHX", "m-su_C")
    CO2_HP.link_components("IHX", "m-ex_C", "Compressor", "m-su")

    # SOURCES AND SINKS

    Gas_cooler_source = MassConnector()
    CO2_HP.add_source(
        "GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C"
    )
    CO2_HP.set_source_properties(
        T=60 + 273.15, fluid="Water", m_dot=2, target="GC_Water", P=3e5
    )

    EV_source = MassConnector()
    CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
    CO2_HP.set_source_properties(
        T=15 + 273.15, fluid="Water", m_dot=2, target="EV_Water", P=2e5
    )

    # CYCLE GUESSES

    P_low = 40 * 1e5
    P_high = 160 * 1e5

    CO2_HP.set_cycle_guess(target="Compressor:su", m_dot=0.08, SH=20, p=P_low)
    CO2_HP.set_cycle_guess(target="Compressor:ex", p=P_high)

    CO2_HP.set_cycle_guess(target="Valve:su", p=P_high, T=30 + 273.15, m_dot=0.08)
    CO2_HP.set_cycle_guess(target="Valve:ex", p=P_low)

    # CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target="IHX:su_C", variable="h", tolerance=1e-3)
    CO2_HP.set_residual_variable(target="IHX:su_H", variable="h", tolerance=1e-3)
    CO2_HP.set_residual_variable(target="IHX:ex_C", variable="h", tolerance=1e-3)
    CO2_HP.set_residual_variable(target="IHX:ex_H", variable="h", tolerance=1e-3)

    CO2_HP.solve()
