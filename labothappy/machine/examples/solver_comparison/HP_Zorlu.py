from labothappy.machine.circuit_fpi import CircuitFPI

from labothappy.connector.mass_connector import MassConnector

from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

from CoolProp.CoolProp import PropsSI

import time
import numpy as np

# SH_vec = np.linspace(1, 10, 10)
# SC_vec = np.linspace(1, 10, 10)

# PP_ev_vec = np.linspace(1,10,10)
# PP_cd_vec = np.linspace(1,10,10)

SH_vec = np.linspace(3,3,1)
SC_vec = np.linspace(3,3,1)

PP_ev_vec = np.linspace(3,3,1)
PP_cd_vec = np.linspace(3,3,1)

eff_rec_vec = np.linspace(0.8, 0.8, 1)

avg_success_time = 0
successes = 0
tries = 0

for PP_ev in PP_ev_vec:
    for PP_cd in PP_cd_vec:
        for eff_rec_val in eff_rec_vec:
            for SH in SH_vec: 
                for SC in SC_vec:
                    tries +=1
                    
                    fluid = "Cyclopentane"
                    HP = CircuitFPI(fluid)
                    
                    "Ignore debug printing"
                    HP.mute_print()
                    HP.mute_plot()
                    
                    "Create components"
                    Compressor = CompressorCstEff()
                    Condenser = HexCstPinch()
                    ExpansionValve = ValveIsenthalpic()
                    Evaporator = HexCstPinch()
                    Recuperator = HexCstEff()
                    
                    "Set component parameters"
                    # Compressor
                    eta_is_cp = 0.8 # -
                    
                    # Condenser
                    Pinch_cd = PP_cd  # K
                    SC_cd = SC # 3  # K
                    DP_h_cd = 0
                    
                    # Evaporator
                    Pinch_ev = PP_ev  # K
                    SH_ev = SH # 5 # K
                    DP_c_ev = 0
        
                    # Recuperator
                    eff_rec = eff_rec_val # 0.8 -
                    DP_h_rec = 0
                    DP_c_rec = 0
                    
                    # Add fluid sources
                    CD_source = MassConnector('Water')
                    T_su_w_cd = 141+273.15
                    P_su_w_cd = 5e5
                    m_dot_w_cd = 10000  # kg/s
                    EV_source = MassConnector('Water')
                    T_su_w_ev = 113.1+273.15
                    P_su_w_ev = 2*1e5
                    m_dot_w_ev = 500  # kg/s
                    
                    #%% Inputs
                    m_dot_ref = 20 # kg/s
                    
                    Compressor.set_parameters(
                        eta_is=eta_is_cp)
                    
                    Condenser.set_parameters(
                        Pinch=Pinch_cd, 
                        Delta_T_sh_sc=SC_cd, 
                        HX_type="condenser",
                        DP_c = 0*1e3,
                        DP_h = DP_h_cd)
                    
                    Evaporator.set_parameters(
                        Pinch=Pinch_ev,
                        Delta_T_sh_sc=SH_ev, 
                        HX_type="evaporator",
                        DP_c = DP_c_ev,
                        DP_h = 0*1e3)
                    
                    Recuperator.set_parameters(
                        eta=eff_rec,
                        DP_c = DP_c_rec,
                        DP_h = DP_h_rec)
                    
                    # Add components to circuit
                    HP.add_component(Compressor, "Compressor")
                    HP.add_component(Condenser, "Condenser")
                    HP.add_component(ExpansionValve, "ExpansionValve")
                    HP.add_component(Evaporator, "Evaporator")
                    HP.add_component(Recuperator, "Recuperator")
                    
                    # Link components with mass connectors
                    HP.link_components("Compressor", "m-ex", "Condenser", "m-su_H")
                    HP.link_components("Condenser", "m-ex_H", "Recuperator", "m-su_H")
                    HP.link_components("Recuperator", "m-ex_H", "ExpansionValve", "m-su")
                    HP.link_components("ExpansionValve", "m-ex", "Evaporator", "m-su_C")
                    HP.link_components("Evaporator", "m-ex_C", "Recuperator", "m-su_C")
                    HP.link_components("Recuperator", "m-ex_C", "Compressor", "m-su")
                    
                    HP.add_source("CD_Water", CD_source, HP.components["Condenser"], "m-su_C")
                    HP.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
                    HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
                    HP.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="EV_Water")
                                    
                    #%% Cycle guess values
                    
                    P_LP_guess = PropsSI("P", "T", T_su_w_ev+50, "Q", 1, fluid)
                    T_sat_LP_guess = PropsSI("T", "P", P_LP_guess, "Q", 1, fluid)
                    
                    P_HP_guess = PropsSI("P", "T", T_su_w_cd+50, "Q", 0, fluid)
                    h_su_vlv_guess = PropsSI('H', 'P', P_HP_guess, 'Q', 0, fluid) - SC_cd - 10000
                    
                    HP.set_cycle_guess(target="Compressor:su", m_dot = m_dot_ref)
                    HP.set_cycle_guess(target="ExpansionValve:su", m_dot = m_dot_ref, h=h_su_vlv_guess, p = P_HP_guess)
                    HP.set_cycle_guess(target="ExpansionValve:ex", p=P_LP_guess)
                    
                    HP.set_cycle_guess(target="Compressor:su", p=P_LP_guess, T=T_sat_LP_guess+10)
                    HP.set_cycle_guess(target="Compressor:ex", p=P_HP_guess)
                    
                    
                    start = time.perf_counter()
                    HP.solve(max_iter=100, method='wegstein')
                    end = time.perf_counter()
                    
                    elapsed = end - start
                    
                    if HP.converged:
                        successes += 1
                        avg_success_time += elapsed
                        # print(f"Converged in {HP._iteration_count} iterations !")
                    else:
                        print(f"PP_cd : {SC}")
                        print(f"PP_ev : {SH}")
                        print(f"------------------")
                    
                    Q_cd = HP.components['Condenser'].model.Q.Q_dot
                    Q_ev = HP.components['Evaporator'].model.Q.Q_dot
                    W_cp = HP.components['Compressor'].model.W.W_dot
                    Q_rec = HP.components['Recuperator'].model.Q.Q_dot
                    COP = Q_cd/W_cp
            
                    if np.mod(tries,100) == 0:
                        print(f"Tries : {tries}")
            
avg_success_time = avg_success_time/successes
conv_prop = successes/tries

print(f"avg_success_time : {avg_success_time}")
print(f"conv_prop : {conv_prop}")

# HP.plot_cycle_Ts()
