# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

from labothappy.machine.circuit_fpi import CircuitFPI
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.compressor.compressor_csteff import CompressorCstEff 
from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic

from CoolProp.CoolProp import PropsSI

def basic_HP(fluid, HSource, CSource, eta_cp, PP_cd, SC_cd, PP_ev, SH_ev, P_low, P_high, mdot):
    
    HP = CircuitFPI(fluid)
    
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

def basic_IHX_HP(fluid, HSource, CSource, eta_cp, eff_rec, PP_cd, SC_cd, PP_ev, SH_ev, P_low, P_high, mdot):
    
    HP = CircuitFPI(fluid)
    
    # Create components
    Compressor = CompressorCstEff()
    Condenser = HexCstPinch()
    Valve = ValveIsenthalpic()
    Evaporator = HexCstPinch()
    Recuperator = HexCstEff()

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

    #%% RECUPERATOR PARAMETERS
    Recuperator.set_parameters(
        eta=eff_rec)
    
    #%% ADD AND LINK COMPONENTS
    
    # Add components
    HP.add_component(Compressor, "Compressor")
    HP.add_component(Condenser, "Condenser")
    HP.add_component(Valve, "Valve")
    HP.add_component(Evaporator, "Evaporator")
    HP.add_component(Recuperator, "Recuperator")
    
    # Link components with mass connectors
    HP.link_components("Compressor", "m-ex", "Condenser", "m-su_H")
    HP.link_components("Condenser", "m-ex_H", "Recuperator", "m-su_H")
    HP.link_components("Recuperator", "m-ex_H", "Valve", "m-su")
    HP.link_components("Valve", "m-ex", "Evaporator", "m-su_C")
    HP.link_components("Evaporator", "m-ex_C", "Recuperator", "m-su_C")
    HP.link_components("Recuperator", "m-ex_C", "Compressor", "m-su")
    
    #%% SOURCES AND SINKS
    
    CD_source = MassConnector()
    HP.add_source("CD_Water", CD_source, HP.components["Condenser"], "m-su_C")
    HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='CD_Water', P = HSource.p)
    
    EV_source = MassConnector()
    HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
    HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
    #%% CYCLE GUESSES

    h_su_vlv_guess = PropsSI('H', 'P', P_high, 'Q', 0, fluid) - SC_cd #- 10000
    T_sat_LP_guess = PropsSI("T", "P", P_low, "Q", 1, fluid)
    
    HP.set_cycle_guess(target='Compressor:su', m_dot = mdot, T=T_sat_LP_guess+SH_ev, p = P_low)
    HP.set_cycle_guess(target='Compressor:ex', p = P_high)

    HP.set_cycle_guess(target='Valve:ex', m_dot = mdot, h=h_su_vlv_guess, p = P_low)
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

if __name__ == "__main__":
    
    study_case = "Zorlu"    
    
    if study_case == "Example":
        
        fluid = 'Propane'
        
        # Hot Source
        T_HS = 60 + 273.15
        p_HS = 3e5
        fluid_HS = 'Water'
        m_dot_HS = 2
    
        # Cold Source
        T_CS = 20 + 273.15
        fluid_CS = 'Water'
        p_CS = 3e5
        m_dot_CS = 10
    
        # Pressure Guesses
        P_high_guess = PropsSI('P', 'T', T_HS, 'Q', 0.5, fluid)
        P_low_guess  = PropsSI('P', 'T', T_CS-5, 'Q', 0.5, fluid)
        
        mdot = 0.1        

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_HS, p = p_HS, m_dot = m_dot_HS)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = fluid_CS, T = T_CS, p = p_CS, m_dot = m_dot_CS)
        
        HP_example = basic_HP(fluid=fluid, HSource=HSource, CSource=CSource, eta_cp=0.7, PP_cd=3, SC_cd=3, PP_ev=3, SH_ev=3, P_low=P_low_guess, P_high=P_high_guess, mdot=mdot)
        HP_example.solve(method='wegstein')          

        Compressor = HP_example.components['Compressor'].model

        print(f"Converged at P_HP = {Compressor.ex.p}, P_LP = {Compressor.su.p}")
    
    elif study_case == "Zorlu":

        fluid = "Cyclopentane"
        
        # Hot Source
        T_HS = 141+273.15
        p_HS = 5e5
        fluid_HS = 'Water'
        m_dot_HS = 10000
    
        # Cold Source
        T_CS = 113.1+273.15
        fluid_CS = 'Water'
        p_CS = 2e5
        m_dot_CS = 500  # kg/s
    
        # Pressure Guesses
        P_high_guess = PropsSI('P', 'T', T_HS, 'Q', 0.5, fluid)
        P_low_guess  = PropsSI('P', 'T', T_CS-5, 'Q', 0.5, fluid)
        
        mdot = 20       

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_HS, p = p_HS, m_dot = m_dot_HS)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = fluid_CS, T = T_CS, p = p_CS, m_dot = m_dot_CS)
        
        HP_zorlu = basic_IHX_HP(fluid=fluid, HSource=HSource, CSource=CSource, eta_cp=0.7, eff_rec=0.8, PP_cd=3, SC_cd=3, PP_ev=3, SH_ev=3, P_low=P_low_guess, P_high=P_high_guess, mdot=mdot)
        HP_zorlu.solve(method='wegstein')          

        Compressor = HP_zorlu.components['Compressor'].model

        print(f"Converged at P_HP = {Compressor.ex.p}, P_LP = {Compressor.su.p}")
        
