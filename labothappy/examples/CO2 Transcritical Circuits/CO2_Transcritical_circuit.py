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
from component.heat_exchanger.hex_csteff_disc import HXEffCstDisc
from component.expander.expander_csteff import ExpanderCstEff 
from component.turbomachinery.pump.steady_state.constant_efficiency.simulation_model import PumpCstEff 
# from component.valve.isenthalpic_valve_P_ex import Isenthalpic_Valve_P_ex
# from component.valve.isenthalpic_valve_x_ex import Isenthalpic_Valve_x_ex
# from component.tank.mixer.simulation_model import Mixer
# from component.tank.Separator.LV_separator import LV_Separator

def basic_CO2_TC(HSource, CSource, eta_pp, eta_exp, eta_gh, PP_cd, SC_cd, P_low, P_high):
    CO2_HP = Circuit('CO2')
    
    # Create components
    Expander = ExpanderCstEff()
    GasHeater = HXEffCstDisc()
    Pump = PumpCstEff()
    Condenser = HXPinchCst()
    
    #%% Pump PARAMETERS
    
    Pump.set_parameters(eta_is=eta_pp)

    #%% Expander PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)
    
    #%% GASCOOLER PARAMETERS
    
    GasHeater.set_parameters(**{
        'eta': eta_gh, 'n_disc' : 20, 'Pinch_min' : 5
    })
    
    #%% EVAPORATOR PARAMETERS
    
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'type_HX': 'condenser'
    })
    
    #%% ADD AND LINK COMPONENTS
    
    # Add components
    CO2_HP.add_component(Expander, "Expander")
    CO2_HP.add_component(GasHeater, "GasHeater")
    CO2_HP.add_component(Pump, "Pump")
    CO2_HP.add_component(Condenser, "Condenser")
            
    # Link components
    CO2_HP.link_components("Pump", "m-ex", "GasHeater", "m-su_C")
    CO2_HP.link_components("GasHeater", "m-ex_C", "Expander", "m-su")
    CO2_HP.link_components("Expander", "m-ex", "Condenser", "m-su_H")
    CO2_HP.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    #%% SOURCES AND SINKS
    
    Gas_heater_source = MassConnector()
    CO2_HP.add_source("GH_Water", Gas_heater_source, CO2_HP.components["GasHeater"], "m-su_H")
    CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GH_Water', P = HSource.p)
    
    CD_source = MassConnector()
    CO2_HP.add_source("CD_Water", CD_source, CO2_HP.components["Condenser"], "m-su_C")
    CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='CD_Water', P = CSource.p)
    
    #%% CYCLE GUESSES
    
    CO2_HP.set_cycle_guess(target='Pump:su', m_dot = 0.08, SC = 5, p = P_low)
    CO2_HP.set_cycle_guess(target='Pump:ex', p = P_high)
        
    CO2_HP.set_cycle_guess(target='Expander:ex', p = P_low)
    
    #%% CYCLE RESIDUAL VARIABLES
    CO2_HP.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)
    CO2_HP.set_residual_variable(target='GasHeater:ex_C', variable='h', tolerance= 1e-3)
    
    return CO2_HP

# def IHX_CO2_HP(HSource, CSource, eta_cp, eta_gc, eta_IHX, PP_ev, SH_ev, P_low, P_high):
#     CO2_HP = Circuit('CO2')
    
#     # Create components
#     Compressor = CompressorCstEff()
#     GasCooler = HXEffCst()
#     IHX = HXEffCst()
#     Valve = Isenthalpic_Valve_P_ex()
#     Evaporator = HXPinchCst()
    
#     #%% COMPRESSOR PARAMETERS
    
#     Compressor.set_parameters(eta_is=eta_cp)
    
#     #%% GASCOOLER PARAMETERS
    
#     GasCooler.set_parameters(**{
#         'eta': eta_gc,
#     })
    
#     #%% IHX PARAMETERS
    
#     IHX.set_parameters(**{
#         'eta': eta_IHX,
#     })
        
#     #%% EVAPORATOR PARAMETERS
    
#     Evaporator.set_parameters(**{
#         'Pinch': PP_ev,
#         'Delta_T_sh_sc': SH_ev,
#         'type_HX': 'evaporator'
#     })
    
#     #%% ADD AND LINK COMPONENTS
    
#     # Add components
#     CO2_HP.add_component(Compressor, "Compressor")
#     CO2_HP.add_component(GasCooler, "GasCooler")
#     CO2_HP.add_component(IHX, "IHX")
#     CO2_HP.add_component(Valve, "Valve")
#     CO2_HP.add_component(Evaporator, "Evaporator")
    
#     # Link components
#     CO2_HP.link_components("Compressor", "m-ex", "GasCooler", "m-su_H")
#     CO2_HP.link_components("GasCooler", "m-ex_H", "IHX", "m-su_H")
#     CO2_HP.link_components("IHX", "m-ex_H", "Valve", "m-su")
#     CO2_HP.link_components("Valve", "m-ex", "Evaporator", "m-su_C")
#     CO2_HP.link_components("Evaporator", "m-ex_C", "IHX", "m-su_C")
#     CO2_HP.link_components("IHX", "m-ex_C", "Compressor", "m-su")
    
#     #%% SOURCES AND SINKS
    
#     Gas_cooler_source = MassConnector()
#     CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
#     CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
#     EV_source = MassConnector()
#     CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
#     CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
#     #%% CYCLE GUESSES
    
#     CO2_HP.set_cycle_guess(target='Compressor:su', m_dot = 0.08, SH = 20, p = P_low)
#     CO2_HP.set_cycle_guess(target='Compressor:ex', p = P_high)

#     CO2_HP.set_cycle_guess(target='Valve:su', p = P_high, T = 30+273.15, m_dot = 0.08)    
#     CO2_HP.set_cycle_guess(target='Valve:ex', p = P_low)
    
#     #%% CYCLE RESIDUAL VARIABLES
#     CO2_HP.set_residual_variable(target='IHX:su_C', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='IHX:su_H', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='IHX:ex_C', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='IHX:ex_H', variable='h', tolerance= 1e-3)
    
#     return CO2_HP

# def Flash_CO2_HP(HSource, CSource, eta_cp, eta_gc, PP_ev, SH_ev, P_low, P_mid, P_high):
#     CO2_HP = Circuit('CO2')

#     # Create components
#     Compressor_1 = CompressorCstEff()
#     Compressor_2 = CompressorCstEff()
#     GasCooler = HXEffCst()
#     Valve_1 = Isenthalpic_Valve_P_ex()
#     Valve_2 = Isenthalpic_Valve_P_ex()
#     Evaporator = HXPinchCst()
#     Mixer_cp = Mixer(n_inlets = 2)
#     Separator = LV_Separator()
    
#     #%% COMPRESSOR PARAMETERS
    
#     Compressor_1.set_parameters(eta_is=eta_cp)
#     Compressor_2.set_parameters(eta_is=eta_cp)
    
#     #%% GASCOOLER PARAMETERS
    
#     GasCooler.set_parameters(**{
#         'eta': eta_gc,
#     })
    
#     #%% EVAPORATOR PARAMETERS
    
#     Evaporator.set_parameters(**{
#         'Pinch': PP_ev,
#         'Delta_T_sh_sc': SH_ev,
#         'type_HX': 'evaporator'
#     })
    
#     #%% ADD AND LINK COMPONENTS
    
#     # Add components
#     CO2_HP.add_component(Compressor_1, "Compressor_1")
#     CO2_HP.add_component(Compressor_2, "Compressor_2")
    
#     CO2_HP.add_component(Mixer_cp, "Mixer")
#     CO2_HP.add_component(GasCooler, "GasCooler")
#     CO2_HP.add_component(Valve_1, "Valve_HP")

#     CO2_HP.add_component(Separator, "Separator")
    
#     CO2_HP.add_component(Valve_2, "Valve_LP")
#     CO2_HP.add_component(Evaporator, "Evaporator")
    
#     # Link components
#     CO2_HP.link_components("Compressor_1", "m-ex", "Mixer", "m-su_1")
#     CO2_HP.link_components("Compressor_2", "m-ex", "Mixer", "m-su_2")
    
#     CO2_HP.link_components("Mixer", "m-ex", "GasCooler", "m-su_H")
#     CO2_HP.link_components("GasCooler", "m-ex_H", "Valve_HP", "m-su")
#     CO2_HP.link_components("Valve_HP", "m-ex", "Separator", "m-su")

#     CO2_HP.link_components("Separator", "m-ex_l", "Valve_LP", "m-su")
    
#     CO2_HP.link_components("Separator", "m-ex_v", "Compressor_1", "m-su")
    
#     CO2_HP.link_components("Valve_LP", "m-ex", "Evaporator", "m-su_C")
#     CO2_HP.link_components("Evaporator", "m-ex_C", "Compressor_2", "m-su")
    
#     #%% SOURCES AND SINKS
    
#     Gas_cooler_source = MassConnector()
#     CO2_HP.add_source("GC_Water", Gas_cooler_source, CO2_HP.components["GasCooler"], "m-su_C")
#     CO2_HP.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GC_Water', P = HSource.p)
    
#     EV_source = MassConnector()
#     CO2_HP.add_source("EV_Water", EV_source, CO2_HP.components["Evaporator"], "m-su_H")
#     CO2_HP.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='EV_Water', P = CSource.p)
    
#     #%% CYCLE GUESSES
    
#     m_dot_tot = 0.08
    
#     CO2_HP.set_cycle_guess(target='Compressor_1:su', m_dot = 0.5*m_dot_tot, SH = 3, p = P_low)
#     CO2_HP.set_cycle_guess(target='Compressor_1:ex', p = P_high)

#     CO2_HP.set_cycle_guess(target='Compressor_2:su', m_dot = 0.5*m_dot_tot, SH = 1, p = P_low)
#     CO2_HP.set_cycle_guess(target='Compressor_2:ex', p = P_high)

#     CO2_HP.set_cycle_guess(target='Valve_HP:ex', p = P_mid)
#     CO2_HP.set_cycle_guess(target='Valve_LP:ex', p = P_low)
    
#     #%% CYCLE RESIDUAL VARIABLES
#     CO2_HP.set_residual_variable(target='Mixer:ex', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Mixer:ex', variable='m_dot', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Evaporator:ex_C', variable='m_dot', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Valve_HP:ex', variable='h', tolerance= 1e-3)
#     CO2_HP.set_residual_variable(target='Valve_HP:ex', variable='m_dot', tolerance= 1e-3)
    
#     return CO2_HP

if __name__ == "__main__":

    T_cold_source = 0+273.15
    T_hot_source = 130+273.15

    eta_is_pp = 0.7
    eta_is_exp = 0.8
    eta_gh = 0.9
    # eta_IHX = 0.3

    PPTD_cd = 5
    SC_cd = 5

    P_high = 140*1e5
    P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
    P_crit_CO2 = PropsSI('PCRIT','CO2')

    P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)    

    study_case = "Simple"    

    HSource = MassConnector()
    HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = 0.1)
    
    CSource = MassConnector()
    CSource.set_properties(fluid = 'R22', T = T_cold_source, p = 5e5, m_dot = 1)
    
    CO2_TC = basic_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, PPTD_cd, SC_cd, P_low_guess, P_high)
    
    CO2_TC.solve()


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