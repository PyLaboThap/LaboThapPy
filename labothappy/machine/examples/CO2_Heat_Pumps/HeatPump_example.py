# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

from machine.circuit_rec import RecursiveCircuit
from connector.mass_connector import MassConnector
from component.heat_exchanger.hex_csteff import HexCstEff
from component.heat_exchanger.hex_cstpinch import HexCstPinch
from component.compressor.compressor_csteff import CompressorCstEff 
from component.valve.valve_isenthalpic import ValveIsenthalpic

from CoolProp.CoolProp import PropsSI

def basic_HP(fluid, HSource, CSource, eta_cp, PP_cd, SC_cd, PP_ev, SH_ev, P_low, P_high, mdot):
    
    HP = RecursiveCircuit(fluid)
    
    # Create components
    Compressor = CompressorCstEff()
    Condenser = HexCstPinch()
    Valve = ValveIsenthalpic()
    Evaporator = HexCstPinch()
    
    #%% COMPRESSOR PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)
    
    #%% GASCOOLER PARAMETERS
    
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'HX_type': 'condenser'
    })
    
    #%% EVAPORATOR PARAMETERS

    Evaporator.set_parameters(**{
        'Pinch': PP_ev,
        'Delta_T_sh_sc': SH_ev,
        'HX_type': 'evaporator'
    })
    
    #%% ADD AND LINK COMPONENTS
    
    # Add components
    HP.add_component(Compressor, "Compressor")
    HP.add_component(Condenser, "Condenser")
    HP.add_component(Valve, "Valve")
    HP.add_component(Evaporator, "Evaporator")
    
    # Link components
    HP.link_components("Compressor", "m-ex", "Condenser", "m-su_H")
    HP.link_components("Condenser", "m-ex_H", "Valve", "m-su")
    HP.link_components("Valve", "m-ex", "Evaporator", "m-su_C")
    HP.link_components("Evaporator", "m-ex_C", "Compressor", "m-su")
    
    #%% SOURCES AND SINKS
    
    CD_source = MassConnector()
    HP.add_source("CD_Water", CD_source, HP.components["Condenser"], "m-su_C")
    HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='CD_Water', P = HSource.p)
    
    EV_source = MassConnector()
    HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
    HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
    #%% CYCLE GUESSES
    
    HP.set_cycle_guess(target='Compressor:su', m_dot = mdot, SH = SH_ev, p = P_low)
    HP.set_cycle_guess(target='Compressor:ex', p = P_high)

    print(P_low)        
    HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    HP.set_residual_variable(target='Valve:ex', variable='h', tolerance= 1e-3)
    HP.set_residual_variable(target='Valve:ex', variable='p', tolerance= 1e-3)

    HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    HP.set_residual_variable(target='Evaporator:ex_C', variable='p', tolerance= 1e-3)

    HP.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
    HP.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)

    HP.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)
    HP.set_residual_variable(target='Condenser:ex_H', variable='p', tolerance= 1e-3)
    
    return HP

# def IHX_CO2_HP(HSource, T_cold_source, eta_cp, eta_gc, eta_IHX, PP_ev, SH_ev, P_low, P_high, m_dot, mute_print_flag):
#     CO2_HP = RecursiveCircuit('CO2')
    
#     n_disc_HX = 50
#     PP_min_HX = 5
    
#     # Create components
#     Compressor = CompressorCstEff()
#     GasCooler = HXEffCstDisc()
#     IHX = HXEffCstDisc()
#     Valve = IsenthalpicValve_P_ex()
#     Evaporator = StorageLatentIsothermalCstePinch()
    
#     #%% COMPRESSOR PARAMETERS
    
#     Compressor.set_parameters(eta_is=eta_cp)
    
#     #%% GASCOOLER PARAMETERS
    
#     GasCooler.set_parameters(**{
#         'eta_max': eta_gc, 'n_disc' : n_disc_HX, 'Pinch_min' : PP_min_HX
#     })
    
#     #%% IHX PARAMETERS
    
#     IHX.set_parameters(**{
#         'eta_max': eta_IHX, 'n_disc' : n_disc_HX, 'Pinch_min' : PP_min_HX
#     })
        
#     #%% EVAPORATOR PARAMETERS

#     Evaporator.set_inputs(**{
#         'sto_fluid': 'Water',
#     })
    
#     Evaporator.set_parameters(**{
#         'Pinch': PP_ev,
#         'Delta_T_sh_sc': SH_ev,
#         'T_sto' : T_cold_source,
#     })
    
    
#     #%% ADD AND LINK COMPONENTS
    
#     # Add components
#     CO2_HP.add_component(Compressor, "Compressor")
#     CO2_HP.add_component(GasCooler, "GasCooler")
#     CO2_HP.add_component(IHX, "IHX")
#     CO2_HP.add_component(Valve, "Valve")
#     CO2_HP.add_component(Evaporator, "Evaporator")
    
#     if mute_print_flag:
#         CO2_HP.mute_print()
        
#     # Link components
#     CO2_HP.link_components("Compressor", "m-ex", "GasCooler", "m-su_H")
#     CO2_HP.link_components("GasCooler", "m-ex_H", "IHX", "m-su_H")
#     CO2_HP.link_components("IHX", "m-ex_H", "Valve", "m-su")
#     CO2_HP.link_components("Valve", "m-ex", "Evaporator", "m-su")
#     CO2_HP.link_components("Evaporator", "m-ex", "IHX", "m-su_C")
#     CO2_HP.link_components("IHX", "m-ex_C", "Compressor", "m-su")
    
#     #%% SOURCES AND SINKS
    
#     Gas_cooler_source = MassConnector()
#     CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
#     CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
#     # EV_source = MassConnector()
#     # CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
#     # CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
#     #%% CYCLE GUESSES
    
#     T_ex_valve_guess = T_cold_source - PP_ev - SH_ev - 5
    
#     CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = m_dot, SH = 20, p = P_low)
#     CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)

#     CO2_HP.set_cycle_guess(target='Valve:su', p = P_high, T = T_ex_valve_guess, m_dot = m_dot)    
#     CO2_HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
#     #%% CYCLE RESIDUAL VARIABLES
#     CO2_HP.set_residual_variable(target='Valve:ex', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Valve:ex', variable='p', tolerance= 1e-3)

#     CO2_HP.set_residual_variable(target='Evaporator:ex', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Evaporator:ex', variable='p', tolerance= 1e-3)

#     CO2_HP.set_residual_variable(target='Compressor:ex', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Compressor:ex', variable='p', tolerance= 1e-3)

#     CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='GasCooler:ex_H', variable='p', tolerance= 1e-3)
    
#     return CO2_HP


if __name__ == "__main__":
    
    study_case = "Simple"    

    # HP Fluid
    fluid = 'Propane'
    
    # Compressor param
    eta_cp = 0.7
    
    # Condenser param
    PP_cd = 3
    SC_cd = 3

    # Evap Param
    PP_ev = 3
    SH_ev = 3
    
    # Hot Source
    T_HS = 40 + 273.15
    p_HS = 3e5
    fluid_HS = 'Water'
    m_dot_HS = 2

    # Cold Source
    T_CS = 20 + 273.15
    fluid_CS = 'Water'
    p_CS = 3e5
    m_dot_CS = 100

    # Pressure Guesses
    P_high_guess = PropsSI('P', 'T', T_HS, 'Q', 0.5, fluid)
    P_low_guess  = PropsSI('P', 'T', T_CS-5, 'Q', 0.5, fluid)
    
    mdot = 0.1
    
    if study_case == "Simple":

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_HS, p = p_HS, m_dot = m_dot_HS)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = fluid_CS, T = T_CS, p = p_CS, m_dot = m_dot_CS)
        
        HP_example = basic_HP(fluid, HSource, CSource, eta_cp, PP_cd, SC_cd, PP_ev, SH_ev, P_low_guess, P_high_guess, mdot)
        HP_example.solve()        

    # elif study_case == "Simple":
        
    #     eta_IHX = 0.8
    #     HSource = MassConnector()
    #     HSource.set_properties(fluid = 'Water', T = T_high, p = p_high, m_dot = 2)
        
    #     CSource = MassConnector()
    #     CSource.set_properties(fluid = 'Water', T = T_low, p = p_low, m_dot = 2)
        
    #     CO2_HP = basic_CO2_HP(HSource, CSource, eta_compressor, eta_GC, DT_pp_ev, SH_ev, P_low_guess, P_high)
        
    #     CO2_HP.solve()      

    