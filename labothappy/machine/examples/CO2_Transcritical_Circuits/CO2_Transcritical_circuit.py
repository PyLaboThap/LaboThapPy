# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:31:53 2025

@author: Basile
"""

import __init__

from machine.circuit_rec import RecursiveCircuit
from CoolProp.CoolProp import PropsSI

from labothappy.connector.mass_connector import MassConnector

from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc
from labothappy.component.expander.expander_csteff import ExpanderCstEff 
from labothappy.component.pump.pump_csteff import PumpCstEff 
from labothappy.component.storage.storage_latent_isoT_cste_pinch import StorageLatentIsothermalCstePinch
from labothappy.component.compressor.compressor_csteff import CompressorCstEff

from labothappy.component.tank.tank_mixer import TankMixer
from labothappy.component.tank.tank_spliter import TankSpliter

#%%

def basic_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_pp, eta_exp, eta_gh, 
               PP_cd, SC_cd, P_low, P_high, m_dot, DP_h_gh = 0, DP_c_gh = 0, DP_h_cond = 0,
               DP_c_cond = 0,mute_print_flag=1):
    
    CO2_TC = RecursiveCircuit('CO2')
    
    # Create components
    Expander = ExpanderCstEff()
    GasHeater = HexCstEffDisc()
    Pump = PumpCstEff()
    Condenser = HexCstPinch()
    
    # Pump PARAMETERS
    
    Pump.set_parameters(eta_is=eta_pp)

    # Expander PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)
    
    # GASCOOLER PARAMETERS
    
    GasHeater.set_parameters(**{
        'eta_max': eta_gh, 'n_disc' : 20, 'Pinch_min' : 5
    })
    
    # EVAPORATOR PARAMETERS
    
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'HX_type': 'condenser'
    })
    
    # ADD AND LINK COMPONENTS
    
    # Add components
    CO2_TC.add_component(Expander, "Expander")
    CO2_TC.add_component(GasHeater, "GasHeater")
    CO2_TC.add_component(Pump, "Pump")
    CO2_TC.add_component(Condenser, "Condenser")
            
    if mute_print_flag:
        CO2_TC.mute_print()
    
    # Link components
    CO2_TC.link_components("Pump", "m-ex", "GasHeater", "m-su_C")
    CO2_TC.link_components("GasHeater", "m-ex_C", "Expander", "m-su")
    CO2_TC.link_components("Expander", "m-ex", "Condenser", "m-su_H")
    CO2_TC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    # SOURCES AND SINKS
    
    Gas_heater_source = MassConnector()
    CO2_TC.add_source("GH_Water", Gas_heater_source, CO2_TC.components["GasHeater"], "m-su_H")
    CO2_TC.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GH_Water', P = HSource.p)
    
    CD_source = MassConnector()
    CO2_TC.add_source("CD_Water", CD_source, CO2_TC.components["Condenser"], "m-su_C")
    CO2_TC.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='CD_Water', P = CSource.p)
    
    # CYCLE GUESSES
    
    CO2_TC.set_cycle_guess(target='Pump:su', m_dot = m_dot, SC = 5, p = P_low)
    CO2_TC.set_cycle_guess(target='Pump:ex', p = P_high)
        
    CO2_TC.set_cycle_guess(target='Expander:ex', p = P_low)
    
    # CYCLE RESIDUAL VARIABLES
    CO2_TC.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='GasHeater:ex_C', variable='h', tolerance= 1e-3)
    
    return CO2_TC

#%%

def REC_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_pp, eta_exp, eta_gh, 
               eta_rec, PP_cd, SC_cd, P_low, P_high, m_dot, DP_h_rec = 0, DP_c_rec = 0, 
               DP_h_gh = 0, DP_c_gh = 0, DP_h_cond = 0, DP_c_cond = 0, mute_print_flag=1):

    CO2_TC = RecursiveCircuit('CO2')
    
    # Create components
    Expander = ExpanderCstEff()
    GasHeater = HexCstEffDisc()
    Rec = HexCstEffDisc()
    Pump = PumpCstEff()
    Condenser = HexCstPinch()
    
    # Pump PARAMETERS
    
    Pump.set_parameters(eta_is=eta_pp)

    # Expander PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)

    # Recuperator PARAMETERS
    
    Rec.set_parameters(**{
        'eta_max': eta_rec, 'n_disc' : 20, 'Pinch_min' : Pinch_min_REC, 'DP_h' : DP_h_rec, 'DP_c' : DP_c_rec,
    })    
    
    # GASCOOLER PARAMETERS
    
    GasHeater.set_parameters(**{
        'eta_max': eta_gh, 'n_disc' : 50, 'Pinch_min' : Pinch_min_GH, 'DP_h' : DP_h_gh, 'DP_c' : DP_c_gh,
    })
    
    # EVAPORATOR PARAMETERS
    
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'HX_type': 'condenser',
        'DP_h' : DP_h_cond, 
        'DP_c' : DP_c_cond, 
    })
    
    # ADD AND LINK COMPONENTS
    
    # Add components
    CO2_TC.add_component(Expander, "Expander")
    CO2_TC.add_component(GasHeater, "GasHeater")
    CO2_TC.add_component(Pump, "Pump")
    CO2_TC.add_component(Condenser, "Condenser")
    CO2_TC.add_component(Rec, "Recuperator")
            
    if mute_print_flag:
        CO2_TC.mute_print()
    
    # Link components
    CO2_TC.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
    CO2_TC.link_components("Recuperator", "m-ex_C", "GasHeater", "m-su_C")
    CO2_TC.link_components("GasHeater", "m-ex_C", "Expander", "m-su")
    CO2_TC.link_components("Expander", "m-ex", "Recuperator", "m-su_H")
    CO2_TC.link_components("Recuperator", "m-ex_H", "Condenser", "m-su_H")
    CO2_TC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    # SOURCES AND SINKS
    
    Gas_heater_source = MassConnector()
    CO2_TC.add_source("GH_Water", Gas_heater_source, CO2_TC.components["GasHeater"], "m-su_H")
    CO2_TC.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GH_Water', P = HSource.p)
    
    CD_source = MassConnector()
    CO2_TC.add_source("CD_Water", CD_source, CO2_TC.components["Condenser"], "m-su_C")
    CO2_TC.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='CD_Water', P = CSource.p)
    
    # CYCLE GUESSES
    
    CO2_TC.set_cycle_guess(target='Pump:su', m_dot = m_dot, SC = 5, p = P_low)
    CO2_TC.set_cycle_guess(target='Pump:ex', p = P_high)

    CO2_TC.set_cycle_guess(target='Expander:su', p = P_high, T = HSource.T, m_dot = m_dot)    
    CO2_TC.set_cycle_guess(target='Expander:ex', p = P_low)
    
    # ITERATION VARIABLES
    
    CO2_TC.set_iteration_variable(target=['Expander:ex'], variable='p', objective = 'Link:Condenser:su_H-p', tol = 1e-2, rel = 1, damping_factor = 0.2, cycle = CO2_TC)
    
    # CYCLE RESIDUAL VARIABLES
    
    CO2_TC.set_residual_variable(target='Recuperator:su_C', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:su_C', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Expander:ex', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Expander:ex', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:ex_H', variable='p', tolerance= 1e-5)
    
    return CO2_TC

#%% 

def REC_CO2_TC_sto(HSource, T_cold_source, Pinch_min_GH, Pinch_min_REC, eta_pp, eta_exp, eta_gh, 
               eta_rec, PP_cd, SC_cd, P_low, P_high, m_dot, DP_h_rec = 0, DP_c_rec = 0, 
               DP_h_gh = 0, DP_c_gh = 0, DP_cond = 0,mute_print_flag=1):
    
    CO2_TC = RecursiveCircuit('CO2')
    
    # Create components
    Expander = ExpanderCstEff()
    GasHeater = HexCstEffDisc()
    Rec = HexCstEffDisc()
    Pump = PumpCstEff()
    Condenser = StorageLatentIsothermalCstePinch()
    
    # Pump PARAMETERS
    
    Pump.set_parameters(eta_is=eta_pp)

    # Expander PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)

    # Recuperator PARAMETERS
    
    Rec.set_parameters(**{
        'eta_max': eta_rec, 'n_disc' : 20, 'Pinch_min' : Pinch_min_REC, 'DP_h' : DP_h_rec, 'DP_c' : DP_c_rec,
    })    
    
    # GASCOOLER PARAMETERS
    
    GasHeater.set_parameters(**{
        'eta_max': eta_gh, 'n_disc' : 20, 'Pinch_min' : Pinch_min_GH, 'DP_h' : DP_h_gh, 'DP_c' : DP_c_gh,
    })
    
    # EVAPORATOR PARAMETERS
    
    Condenser.set_inputs(**{
        'sto_fluid': 'Water',
    })
    
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'T_sto' : T_cold_source,
        'DP' : DP_cond, 
    })
    
    # ADD AND LINK COMPONENTS
    
    # Add components
    CO2_TC.add_component(Expander, "Expander")
    CO2_TC.add_component(GasHeater, "GasHeater")
    CO2_TC.add_component(Pump, "Pump")
    CO2_TC.add_component(Condenser, "Condenser")
    CO2_TC.add_component(Rec, "Recuperator")
            
    if mute_print_flag:
        CO2_TC.mute_print()
    
    # Link components
    CO2_TC.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
    CO2_TC.link_components("Recuperator", "m-ex_C", "GasHeater", "m-su_C")
    CO2_TC.link_components("GasHeater", "m-ex_C", "Expander", "m-su")
    CO2_TC.link_components("Expander", "m-ex", "Recuperator", "m-su_H")
    CO2_TC.link_components("Recuperator", "m-ex_H", "Condenser", "m-su")
    CO2_TC.link_components("Condenser", "m-ex", "Pump", "m-su")
    
    # SOURCES AND SINKS
    
    Gas_heater_source = MassConnector()
    CO2_TC.add_source("GH_Water", Gas_heater_source, CO2_TC.components["GasHeater"], "m-su_H")
    CO2_TC.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GH_Water', P = HSource.p)
       
    # CYCLE GUESSES
    
    CO2_TC.set_cycle_guess(target='Pump:su', m_dot = m_dot, SC = 5, p = P_low)
    CO2_TC.set_cycle_guess(target='Pump:ex', p = P_high)

    CO2_TC.set_cycle_guess(target='Expander:su', p = P_high, T = HSource.T, m_dot = m_dot)    
    CO2_TC.set_cycle_guess(target='Expander:ex', p = P_low)
    
    # ITERATION VARIABLES
    
    CO2_TC.set_iteration_variable(target=['Expander:ex'], variable='p', objective = 'Link:Condenser:su-p', tol = 1e-2, rel = 1, damping_factor = 0.2, cycle = CO2_TC)
    
    # CYCLE RESIDUAL VARIABLES
    
    CO2_TC.set_residual_variable(target='Recuperator:su_C', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:su_C', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:su_H', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:su_H', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Expander:ex', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Expander:ex', variable='h', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:ex_C', variable='p', tolerance= 1e-5)
    CO2_TC.set_residual_variable(target='Recuperator:ex_H', variable='p', tolerance= 1e-5)
    
    return CO2_TC

def Recomp_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_pp, eta_exp, eta_cp, eta_rec, eta_gh, 
               PP_cd, SC_cd, P_low, P_high, m_dot, spliter_frac = 0.5, DP_h_gh = 0, DP_c_gh = 0, DP_h_cond = 0,
               DP_c_cond = 0 ,mute_print_flag=1):
    
    CO2_TC = RecursiveCircuit('CO2')
    
    rep_spliter = [spliter_frac, 1-spliter_frac]
    
    # Create components
    Expander = ExpanderCstEff()
    GasHeater = HexCstEffDisc()
    Compressor = CompressorCstEff()
    RecupHT = HexCstEffDisc()
    RecupLT = HexCstEffDisc()
    Pump = PumpCstEff()
    Condenser = HexCstPinch()
    Spliter = TankSpliter(outlet_repartition=rep_spliter)
    Mixer = TankMixer(n_inlets = 2)
    
    # Pump PARAMETERS
    
    Pump.set_parameters(eta_is=eta_pp)

    # Expander PARAMETERS
    
    Expander.set_parameters(eta_is=eta_exp)
    
    # GASCOOLER PARAMETERS
    
    GasHeater.set_parameters(**{
        'eta_max': eta_gh, 'n_disc' : 20, 'Pinch_min' : 5
    })
    
    # EVAPORATOR PARAMETERS
    
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'HX_type': 'condenser'
    })
    
    # Recup LT PARAMETERS
    
    RecupLT.set_parameters(**{
        'eta_max': eta_rec, 
        'n_disc' : 20, 
        'Pinch_min' : Pinch_min_REC, 
        })

    # Recup HT PARAMETERS
    
    RecupHT.set_parameters(**{
        'eta_max': eta_rec, 
        'n_disc' : 20, 
        'Pinch_min' : Pinch_min_REC, 
        })    
    
    # GASCOOLER PARAMETERS
    
    Compressor.set_parameters(eta_is=eta_cp)
    
    # ADD AND LINK COMPONENTS
    
    # Add components
    CO2_TC.add_component(Expander, "Expander")
    CO2_TC.add_component(GasHeater, "GasHeater")
    CO2_TC.add_component(Pump, "Pump")
    CO2_TC.add_component(Condenser, "Condenser")
    CO2_TC.add_component(Compressor, "Compressor")
    CO2_TC.add_component(RecupLT, "RecupLT")
    CO2_TC.add_component(RecupHT, "RecupHT")
    CO2_TC.add_component(Spliter, "Spliter")
    CO2_TC.add_component(Mixer, "Mixer")
          
    if mute_print_flag:
        CO2_TC.mute_print()
    
    # Link components
    CO2_TC.link_components("Pump", "m-ex", "RecupLT", "m-su_C")
    CO2_TC.link_components("RecupLT", "m-ex_C", "Mixer", "m-su_1")
    CO2_TC.link_components("Mixer", "m-ex", "RecupHT", "m-su_C")
    CO2_TC.link_components("RecupHT", "m-ex_C", "GasHeater", "m-su_C")
    CO2_TC.link_components("GasHeater", "m-ex_C", "Expander", "m-su")
    
    CO2_TC.link_components("Expander", "m-ex", "RecupHT", "m-su_H")
    CO2_TC.link_components("RecupHT", "m-ex_H", "RecupLT", "m-su_H")
    CO2_TC.link_components("RecupLT", "m-ex_H", "Spliter", "m-su")

    CO2_TC.link_components("Spliter", "m-ex_1", "Condenser", "m-su_H")
    CO2_TC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    CO2_TC.link_components("Spliter", "m-ex_2", "Compressor", "m-su")
    CO2_TC.link_components("Compressor", "m-ex", "Mixer", "m-su_2")
    
    # SOURCES AND SINKS
    
    Gas_heater_source = MassConnector()
    CO2_TC.add_source("GH_Water", Gas_heater_source, CO2_TC.components["GasHeater"], "m-su_H")
    CO2_TC.set_source_properties(T=HSource.T, fluid=HSource.fluid, m_dot=HSource.m_dot, target='GH_Water', P = HSource.p)
    
    CD_source = MassConnector()
    CO2_TC.add_source("CD_Water", CD_source, CO2_TC.components["Condenser"], "m-su_C")
    CO2_TC.set_source_properties(T=CSource.T, fluid=CSource.fluid, m_dot=CSource.m_dot, target='CD_Water', P = CSource.p)
    
    # CYCLE GUESSES
    
    CO2_TC.set_cycle_guess(target='Pump:su', m_dot = m_dot*rep_spliter[0], SC = 5, p = P_low)
    CO2_TC.set_cycle_guess(target='Pump:ex', p = P_high)
        
    CO2_TC.set_cycle_guess(target='Expander:ex', p = P_low)
    
    CO2_TC.set_cycle_guess(target='Compressor:su', m_dot = m_dot*rep_spliter[1], T = CSource.T + 10, p = P_low)
    CO2_TC.set_cycle_guess(target='Compressor:ex', p = P_high)

    CO2_TC.set_cycle_guess(target='RecupHT:su_C', m_dot = m_dot, T = CSource.T+40, p = P_high)
    CO2_TC.set_cycle_guess(target='RecupHT:su_H', m_dot = m_dot, T = CSource.T+70, p = P_low)
        
    CO2_TC.set_fixed_properties(target='Pump:su', SC = SC_cd)    
    CO2_TC.set_iteration_variable(target=['Expander:ex'], variable='p', objective = 'Link:Condenser:su_H-p', tol = 1e-2, rel = 1, damping_factor = 0.2, cycle = CO2_TC)
    
    # CYCLE RESIDUAL VARIABLES
    CO2_TC.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='GasHeater:ex_C', variable='h', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='Expander:ex', variable='p', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='Expander:ex', variable='h', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='RecupLT:ex_H', variable='p', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='RecupLT:ex_H', variable='h', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='Spliter:ex_2', variable='h', tolerance= 1e-3)
    CO2_TC.set_residual_variable(target='Spliter:ex_2', variable='p', tolerance= 1e-3)
    
    return CO2_TC

#%%

if __name__ == "__main__": 

    study_case = "Recomp"

    if study_case == "Simple":
        T_cold_source = 0.1+273.15
        T_hot_source = 130+273.15

        m_dot = 0.08

        eta_is_pp = 0.7
        eta_is_exp = 0.8
        eta_gh = 0.9
        eta_rec = 0.3

        PPTD_cd = 5
        SC_cd = 5

        Pinch_min_GH = 3
        Pinch_min_REC = 3

        P_high = 140*1e5
        P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
        P_crit_CO2 = PropsSI('PCRIT','CO2')

        P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)   
        
        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = 0.1)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'Water', T = T_cold_source, p = 5e5, m_dot = 10)
        
        CO2_TC = basic_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_is_pp, eta_is_exp, eta_gh, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot,
                            DP_h_gh = 50*1e3, DP_c_gh = 2*1e5, DP_h_cond = 1*1e5, DP_c_cond = 50*1e3, mute_print_flag=0)
                
        CO2_TC.solve()

    elif study_case == "Recup_sto":
        T_cold_source = 0.1+273.15
        T_hot_source = 130+273.15

        m_dot = 0.08

        eta_is_pp = 0.7
        eta_is_exp = 0.8
        eta_gh = 0.9
        eta_rec = 0.3

        PPTD_cd = 5
        SC_cd = 5

        Pinch_min_GH = 3
        Pinch_min_REC = 3

        P_high = 140*1e5
        P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
        P_crit_CO2 = PropsSI('PCRIT','CO2')

        P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)   
        
        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = 0.1)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'R22', T = T_cold_source, p = 5e5, m_dot = 1)
        
        CO2_TC = REC_CO2_TC_sto(HSource, T_cold_source, Pinch_min_GH, Pinch_min_REC, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot,
                            DP_h_rec = 1*1e5, DP_c_rec = 2*1e5, DP_h_gh = 50*1e3, DP_c_gh = 2*1e5, DP_cond = 1*1e5, mute_print_flag=0)
                
        CO2_TC.solve()
     
    elif study_case == "Recup":
        
        T_cold_source = 0.1+273.15
        T_hot_source = 130+273.15

        m_dot = 3.00123861e+01

        eta_is_pp = 0.8
        eta_is_exp = 0.9
        eta_gh = 0.95
        eta_rec = 0.8

        PPTD_cd = 5
        SC_cd = 0.1

        Pinch_min_GH = 5
        Pinch_min_REC = 0

        P_high = 140*1e5
        P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
        P_crit_CO2 = PropsSI('PCRIT','CO2')

        P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)   
        
        T_cold_source = 0.1+273.15
        
        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = m_dot*1)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'Water', T = T_cold_source, p = 5e5, m_dot = 10000)
        
        DP_h_rec = 1*1e5
        DP_c_rec = 2*1e5
        
        DP_h_gh = 50*1e3
        DP_c_gh = 2*1e5
        
        DP_h_cond = 1*1e5
        DP_c_cond = 50*1e3
        
        # for i in range(100):
        CO2_TC = REC_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot,
                            DP_h_rec = DP_h_rec, DP_c_rec = DP_c_rec, DP_h_gh = DP_h_gh, DP_c_gh = DP_c_gh, DP_h_cond = DP_h_cond, DP_c_cond = DP_c_cond, mute_print_flag=0)
                
        CO2_TC.solve()
                
    elif study_case == "Recomp":
        
        T_cold_source = 0.1+273.15
        T_hot_source = 130+273.15

        m_dot = 3.00123861e+01

        eta_is_pp = 0.8
        eta_is_cp = 0.8
        eta_is_exp = 0.9
        eta_gh = 0.95
        eta_rec = 0.95

        PPTD_cd = 5
        SC_cd = 0.1

        Pinch_min_GH = 5
        Pinch_min_REC = 0

        P_high = 120*1e5
        P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
        P_crit_CO2 = PropsSI('PCRIT','CO2')

        P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)   
        
        T_cold_source = 0.1+273.15
        
        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = m_dot*1)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = 'Water', T = T_cold_source, p = 5e5, m_dot = 10000)
        
        DP_h_rec = 1*1e5
        DP_c_rec = 2*1e5
        
        DP_h_gh = 50*1e3
        DP_c_gh = 2*1e5
        
        DP_h_cond = 1*1e5
        DP_c_cond = 50*1e3
        
        CO2_TC = Recomp_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_is_pp, eta_is_exp, 
                               eta_is_cp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, 
                               m_dot, spliter_frac = 1, mute_print_flag=1)
        
        # CO2_TC = Recomp_CO2_TC(HSource, CSource, Pinch_min_GH, Pinch_min_REC, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot,
        #                     DP_h_rec = DP_h_rec, DP_c_rec = DP_c_rec, DP_h_gh = DP_h_gh, DP_c_gh = DP_c_gh, DP_h_cond = DP_h_cond, DP_c_cond = DP_c_cond, mute_print_flag=0)
         
        CO2_TC.solve()
                
        eta = (CO2_TC.components['Expander'].model.W.W_dot - CO2_TC.components['Pump'].model.W.W_dot - CO2_TC.components['Compressor'].model.W.W_dot)/(CO2_TC.components['GasHeater'].model.Q)
        
        print(f"eta_th : {eta}")
        print(f"Q_dot Recup HT : {CO2_TC.components['RecupHT'].model.Q}")
        print(f"Q_dot Recup LT : {CO2_TC.components['RecupLT'].model.Q}")
        
    CO2_TC.plot_cycle_Ts()
    CO2_TC.Ts_gif()