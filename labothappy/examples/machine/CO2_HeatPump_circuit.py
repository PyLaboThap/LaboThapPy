# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

import __init__

from machine.circuit import Circuit
from CoolProp.CoolProp import PropsSI

from connector.mass_connector import MassConnector

from component.heat_exchanger.steady_state.pinch_cst.simulation_model import HXPinchCst
from component.heat_exchanger.steady_state.cst_efficiency.simulation_model import HXEffCst
from component.volumetric_machine.compressor.steady_state.constant_isentropic_efficiency.simulation_model import CompressorCstEff 
from component.valve.isenthalpic_valve import Isenthalpic_Valve

def basic_CO2_HP(HSource, CSource, eta_cp, eta_gc, PP_ev, SH_ev, P_low, P_high):
    CO2_HP = Circuit('CO2')
    
    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCst()
    Valve = Isenthalpic_Valve()
    Evaporator = HXPinchCst()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)
    
    #%% GASCOOLER PARAMETERS
    
    GasCooler.set_parameters(**{
        'eta': eta_gc,
    })
    
    #%% EVAPORATOR PARAMETERS
    
    Evaporator.set_parameters(**{
        'Pinch': PP_ev,
        'Delta_T_sh_sc': SH_ev,
        'type_HX': 'evaporator'
    })
    
    #%% ADD AND LINK COMPONENTS
    
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
    
    #%% SOURCES AND SINKS
    
    Gas_cooler_source = MassConnector()
    CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
    CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
    EV_source = MassConnector()
    CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
    CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
    #%% CYCLE GUESSES
    
    CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = 0.08, SH = 5, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)
        
    CO2_HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='h', tolerance= 1e-3)
    
    return CO2_HP

def IHX_CO2_HP(HSource, CSource, eta_cp, eta_gc, eta_IHX, PP_ev, SH_ev, P_low, P_high):
    CO2_HP = Circuit('CO2')
    
    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCst()
    IHX = HXEffCst()
    Valve = Isenthalpic_Valve()
    Evaporator = HXPinchCst()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)
    
    #%% GASCOOLER PARAMETERS
    
    GasCooler.set_parameters(**{
        'eta': eta_gc,
    })
    
    #%% IHX PARAMETERS
    
    IHX.set_parameters(**{
        'eta': eta_IHX,
    })
        
    #%% EVAPORATOR PARAMETERS
    
    Evaporator.set_parameters(**{
        'Pinch': PP_ev,
        'Delta_T_sh_sc': SH_ev,
        'type_HX': 'evaporator'
    })
    
    #%% ADD AND LINK COMPONENTS
    
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
    
    #%% SOURCES AND SINKS
    
    Gas_cooler_source = MassConnector()
    CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
    CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
    EV_source = MassConnector()
    CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
    CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
    #%% CYCLE GUESSES
    
    CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = 0.08, SH = 20, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)

    CO2_HP.set_cycle_guess(target='Valve:su', p = P_high, T = 30+273.15, m_dot = 0.08)    
    CO2_HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='IHX:su_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='IHX:su_H', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='IHX:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='IHX:ex_H', variable='h', tolerance= 1e-3)
    
    return CO2_HP

T_cold_source = -10 + 273.15
P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
P_crit_CO2 = PropsSI('PCRIT','CO2')

P_low_guess = min(1.2*P_sat_T_CSource,0.99*P_crit_CO2)

P_high = 160*1e5

HSource = MassConnector()
HSource.set_properties(fluid = 'Water', T = 150 + 273.15, p = 3e5, m_dot = 8)

CSource = MassConnector()
CSource.set_properties(fluid = 'R22', T = T_cold_source, p = 5e5, m_dot = 8)

CO2_HP = IHX_CO2_HP(HSource, CSource, 0.7, 0.95, 0.8, 3, 3, P_low_guess, P_high)
    
CO2_HP.solve()


# if __name__ == "__main__":
#     P_low_guess = 40*1e5
#     P_high = 160*1e5
    
#     HSource = MassConnector()
#     HSource.set_properties(fluid = 'Water', T = 60 + 273.15, p = 3e5, m_dot = 2)
    
#     CSource = MassConnector()
#     CSource.set_properties(fluid = 'Water', T = 15 + 273.15, p = 3e5, m_dot = 2)
    
#     CO2_HP = IHX_CO2_HP(HSource, CSource, 0.7, 0.95, 0.8, 3, 3, P_low_guess, P_high)
    
#     CO2_HP.solve()
