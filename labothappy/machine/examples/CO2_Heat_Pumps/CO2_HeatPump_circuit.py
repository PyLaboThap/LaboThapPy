# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

import __init__

from machine.circuit_rec import RecursiveCircuit
from CoolProp.CoolProp import PropsSI

from connector.mass_connector import MassConnector

from component.heat_exchanger.hex_cstpinch import HXPinchCst
from component.heat_exchanger.hex_csteff import HXEffCst
from component.heat_exchanger.hex_csteff_disc import HXEffCstDisc
from component.compressor.compressor_csteff import CompressorCstEff 
from component.valve.isenthalpic_valve_P_ex import IsenthalpicValve_P_ex
from component.expander.expander_csteff import ExpanderCstEff
from component.tank.tank_mixer import Mixer
from component.tank.tank_LV_separator import LV_Separator

def basic_CO2_HP(HSource, CSource, eta_cp, eta_gc, PP_ev, SH_ev, P_low, P_high):
    
    CO2_HP = RecursiveCircuit('CO2')
    
    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCst()
    Valve = IsenthalpicValve_P_ex()
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
    
    CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = 0.16, SH = 5, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)
        
    CO2_HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Valve:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Valve:ex', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='p', tolerance= 1e-3)
    
    return CO2_HP

def Exp_CO2_HP(HSource, CSource, eta_cp, eta_exp, eta_gc, PP_ev, SH_ev, P_low, P_high):
    CO2_HP = RecursiveCircuit('CO2')
    
    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCst()
    Expander = ExpanderCstEff()
    Evaporator = HXPinchCst()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)

    #%% EXPANDER PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)
    
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
    CO2_HP.add_component(Expander, "Expander")
    CO2_HP.add_component(Evaporator, "Evaporator")
    
    # Link components
    CO2_HP.link_components("Compressor", "m-ex", "GasCooler", "m-su_H")
    CO2_HP.link_components("GasCooler", "m-ex_H", "Expander", "m-su")
    CO2_HP.link_components("Expander", "m-ex", "Evaporator", "m-su_C")
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
        
    CO2_HP.set_cycle_guess(target='Expander:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='h', tolerance= 1e-3)
    
    return CO2_HP

def IHX_CO2_HP(HSource, CSource, eta_cp, eta_gc, eta_IHX, PP_ev, SH_ev, P_low, P_high, m_dot, print_flag):
    CO2_HP = RecursiveCircuit('CO2')
    
    n_disc_HX = 50
    PP_min_HX = 5
    
    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCstDisc()
    IHX = HXEffCstDisc()
    Valve = IsenthalpicValve_P_ex()
    Evaporator = HXPinchCst()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)
    
    #%% GASCOOLER PARAMETERS
    
    GasCooler.set_parameters(**{
        'eta': eta_gc, 'n_disc' : n_disc_HX, 'Pinch_min' : PP_min_HX
    })
    
    #%% IHX PARAMETERS
    
    IHX.set_parameters(**{
        'eta': eta_IHX, 'n_disc' : n_disc_HX, 'Pinch_min' : PP_min_HX
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
    
    if not print_flag:
        CO2_HP.mute_print()
    
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
    
    CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = m_dot, SH = 20, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)

    CO2_HP.set_cycle_guess(target='Valve:su', p = P_high, T = 30+273.15, m_dot = m_dot)    
    CO2_HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Valve:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Valve:ex', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='p', tolerance= 1e-3)
    
    return CO2_HP

def IHX_EXP_CO2_HP(HSource, CSource, eta_cp, eta_gc, eta_IHX, eta_exp, PP_ev, SH_ev, P_low, P_high, m_dot):
    CO2_HP = RecursiveCircuit('CO2')
    
    n_disc_HX = 50
    PP_min_HX = 5
    
    # Create components
    Compressor = CompressorCstEff()
    GasCooler = HXEffCstDisc()
    IHX = HXEffCstDisc()
    Expander = ExpanderCstEff()
    Evaporator = HXPinchCst()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)
    
    #%% GASCOOLER PARAMETERS
    
    GasCooler.set_parameters(**{
        'eta': eta_gc, 'n_disc' : n_disc_HX, 'Pinch_min' : PP_min_HX
    })
    
    #%% IHX PARAMETERS
    
    IHX.set_parameters(**{
        'eta': eta_IHX, 'n_disc' : n_disc_HX, 'Pinch_min' : PP_min_HX
    })
        
    #%% EXPANDER PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)
    
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
    CO2_HP.add_component(Expander, "Expander")
    CO2_HP.add_component(Evaporator, "Evaporator")
    
    # Link components
    CO2_HP.link_components("Compressor", "m-ex", "GasCooler", "m-su_H")
    CO2_HP.link_components("GasCooler", "m-ex_H", "IHX", "m-su_H")
    CO2_HP.link_components("IHX", "m-ex_H", "Expander", "m-su")
    CO2_HP.link_components("Expander", "m-ex", "Evaporator", "m-su_C")
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
    
    CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = m_dot, SH = 20, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)

    CO2_HP.set_cycle_guess(target='Expander:su', p = P_high, T = 30+273.15, m_dot = m_dot)    
    CO2_HP.set_cycle_guess(target='Expander:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Expander:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Expander:ex', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)

    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='p', tolerance= 1e-3)
    
    return CO2_HP

def Flash_CO2_HP_Parallel_CP(HSource, CSource, eta_cp, eta_gc, PP_ev, SH_ev, P_low, P_mid, P_high):
    CO2_HP = RecursiveCircuit('CO2')

    # Create components
    Compressor_1 = CompressorCstEff()
    Compressor_2 = CompressorCstEff()
    GasCooler = HXEffCst()
    Valve_1 = IsenthalpicValve_P_ex()
    Valve_2 = IsenthalpicValve_P_ex()
    Evaporator = HXPinchCst()
    Mixer_cp = Mixer(n_inlets = 2)
    Separator = LV_Separator()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor_1.set_parameters(eta_is=eta_cp)
    Compressor_2.set_parameters(eta_is=eta_cp)
    
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
    CO2_HP.add_component(Compressor_1, "Compressor_1")
    CO2_HP.add_component(Compressor_2, "Compressor_2")
    
    CO2_HP.add_component(Mixer_cp, "Mixer")
    CO2_HP.add_component(GasCooler, "GasCooler")
    CO2_HP.add_component(Valve_1, "Valve_HP")

    CO2_HP.add_component(Separator, "Separator")
    
    CO2_HP.add_component(Valve_2, "Valve_LP")
    CO2_HP.add_component(Evaporator, "Evaporator")
    
    # Link components
    CO2_HP.link_components("Compressor_1", "m-ex", "Mixer", "m-su_1")
    CO2_HP.link_components("Compressor_2", "m-ex", "Mixer", "m-su_2")
    
    CO2_HP.link_components("Mixer", "m-ex", "GasCooler", "m-su_H")
    CO2_HP.link_components("GasCooler", "m-ex_H", "Valve_HP", "m-su")
    CO2_HP.link_components("Valve_HP", "m-ex", "Separator", "m-su")

    CO2_HP.link_components("Separator", "m-ex_l", "Valve_LP", "m-su")
    
    CO2_HP.link_components("Separator", "m-ex_v", "Compressor_1", "m-su")
    
    CO2_HP.link_components("Valve_LP", "m-ex", "Evaporator", "m-su_C")
    CO2_HP.link_components("Evaporator", "m-ex_C", "Compressor_2", "m-su")
    
    #%% SOURCES AND SINKS
    
    Gas_cooler_source = MassConnector()
    CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
    CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
    EV_source = MassConnector()
    CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
    CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
    #%% CYCLE GUESSES
    
    m_dot_tot = 0.16
    
    CO2_HP.set_cycle_guess(target='Compressor_1:su', m_dot = 0.5*m_dot_tot, SH = 3, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor_1:ex', p = P_high)

    CO2_HP.set_cycle_guess(target='Compressor_2:su', m_dot = 0.5*m_dot_tot, SH = 3, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor_2:ex', p = P_high)

    CO2_HP.set_cycle_guess(target='Valve_HP:ex', p = P_mid)
    CO2_HP.set_cycle_guess(target='Valve_LP:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Mixer:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Mixer:ex', variable='m_dot', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='m_dot', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Valve_HP:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Valve_HP:ex', variable='m_dot', tolerance= 1e-3)
    
    return CO2_HP

def Flash_CO2_HP_Series_CP(HSource, CSource, eta_cp, eta_gc, PP_ev, SH_ev, P_low, P_mid, P_high):
    CO2_HP = RecursiveCircuit('CO2')

    # Create components
    Compressor_1 = CompressorCstEff()
    Compressor_2 = CompressorCstEff()
    GasCooler = HXEffCst()
    Valve_1 = IsenthalpicValve_P_ex()
    Valve_2 = IsenthalpicValve_P_ex()
    Evaporator = HXPinchCst()
    Mixer_cp = Mixer(n_inlets = 2)
    Separator = LV_Separator()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor_1.set_parameters(eta_is=eta_cp)
    Compressor_2.set_parameters(eta_is=eta_cp)
    
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
    CO2_HP.add_component(Compressor_1, "Compressor_1")
    CO2_HP.add_component(Compressor_2, "Compressor_2")
    
    CO2_HP.add_component(Mixer_cp, "Mixer")
    CO2_HP.add_component(GasCooler, "GasCooler")
    CO2_HP.add_component(Valve_1, "Valve_HP")

    CO2_HP.add_component(Separator, "Separator")
    
    CO2_HP.add_component(Valve_2, "Valve_LP")
    CO2_HP.add_component(Evaporator, "Evaporator")
    
    # Link components
    CO2_HP.link_components("Evaporator", "m-ex_C", "Compressor_1", "m-su")
    CO2_HP.link_components("Compressor_1", "m-ex", "Mixer", "m-su_1")
    CO2_HP.link_components("Mixer", "m-ex", "Compressor_2", "m-su")
    
    CO2_HP.link_components("Compressor_2", "m-ex", "GasCooler", "m-su_H")
    CO2_HP.link_components("GasCooler", "m-ex_H", "Valve_HP", "m-su")
    CO2_HP.link_components("Valve_HP", "m-ex", "Separator", "m-su")

    CO2_HP.link_components("Separator", "m-ex_l", "Valve_LP", "m-su")
    CO2_HP.link_components("Separator", "m-ex_v", "Mixer", "m-su_2")
    
    CO2_HP.link_components("Valve_LP", "m-ex", "Evaporator", "m-su_C")
    
    #%% SOURCES AND SINKS
    
    Gas_cooler_source = MassConnector()
    CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
    CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
    EV_source = MassConnector()
    CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
    CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
    #%% CYCLE GUESSES
    
    m_dot_tot = 0.1
    
    CO2_HP.set_cycle_guess(target='Compressor_1:su', m_dot = 0.5*m_dot_tot, SH = 3, p = P_low)
    CO2_HP.set_cycle_guess(target='Compressor_1:ex', p = P_mid)

    CO2_HP.set_cycle_guess(target='Compressor_2:su', m_dot = m_dot_tot, SH = 3, p = P_mid)
    CO2_HP.set_cycle_guess(target='Compressor_2:ex', p = P_high)

    CO2_HP.set_cycle_guess(target='Valve_HP:ex', p = P_mid)
    CO2_HP.set_cycle_guess(target='Valve_LP:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Mixer:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Mixer:ex', variable='m_dot', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='m_dot', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Valve_HP:ex', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='Valve_HP:ex', variable='m_dot', tolerance= 1e-3)
    
    return CO2_HP

if __name__ == "__main__":

    study_case = "Flash"    

    # Pressure levels
    P_low_guess = 40*1e5
    P_high = 180*1e5
    
    # Compressor param
    eta_compressor = 0.7
    
    # GasCooler param
    eta_GC = 0.9

    # Evap Param
    DT_pp_ev = 3
    SH_ev = 3
    
    # Hot Source
    T_high = 100 + 273.15
    p_high = 3e5
    fluid_high = 'Water'
    m_dot_high = 2

    # Cold Source
    T_low = 15 + 273.15
    fluid_low = 'Water'
    p_low = 3e5
    m_dot_low = 2

    if study_case == "IHX":

        eta_IHX = 0.8

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_high, p = p_high, m_dot = 2)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'Water', T = T_low, p = p_low, m_dot = 2)
        
        CO2_HP = IHX_CO2_HP(HSource, CSource, eta_compressor, eta_GC, eta_IHX, DT_pp_ev, SH_ev, P_low_guess, P_high)
        
        CO2_HP.solve()        

    elif study_case == "Simple":

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_high, p = p_high, m_dot = 2)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'Water', T = T_low, p = p_low, m_dot = 2)
        
        CO2_HP = basic_CO2_HP(HSource, CSource, eta_compressor, eta_GC, DT_pp_ev, SH_ev, P_low_guess, P_high)
        
        CO2_HP.solve()      

    elif study_case == "Expander":

        # Compressor param
        eta_exp = 0.7

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_high, p = p_high, m_dot = 2)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'Water', T = T_low, p = p_low, m_dot = 2)
        
        CO2_HP = Exp_CO2_HP(HSource, CSource, eta_compressor, eta_exp, eta_GC, DT_pp_ev, SH_ev, P_low_guess, P_high)
        
        CO2_HP.solve()      
    
        COP = CO2_HP.components["GasCooler"].model.Q_dot.Q_dot/(CO2_HP.components["Compressor"].model.W_dot.W_dot - CO2_HP.components["Expander"].model.W_exp.W_dot)
    
    # print(f"COP : {COP}")

# """
# T_cold_source = -10 + 273.15
# P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
# P_crit_CO2 = PropsSI('PCRIT','CO2')

# P_low_guess = min(1.2*P_sat_T_CSource,0.99*P_crit_CO2)

# P_high = 160*1e5

# HSource = MassConnector()
# HSource.set_properties(fluid = 'Water', T = 150 + 273.15, p = 3e5, m_dot = 8)

# CSource = MassConnector()
# CSource.set_properties(fluid = 'R22', T = T_cold_source, p = 5e5, m_dot = 8)

# CO2_HP = IHX_CO2_HP(HSource, CSource, 0.7, 0.95, 0.8, 3, 3, P_low_guess, P_high)
    
# CO2_HP.solve()

# """