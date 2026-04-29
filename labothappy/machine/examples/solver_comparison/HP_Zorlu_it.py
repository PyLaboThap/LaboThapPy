from labothappy.machine.circuit_it import IterativeCircuit

from labothappy.connector.mass_connector import MassConnector

from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

from CoolProp.CoolProp import PropsSI

import time
import numpy as np

# SC_cd_vec = np.linspace(1,10,10)
# SH_ev_vec = np.linspace(1,10,10)
# # eff_rec_vec = np.linspace(0, 1, 11)

# PP_ev_vec = np.linspace(1,10,10)
# PP_cd_vec = np.linspace(1,10,10)

SC_cd_vec = np.linspace(3,3,1)
SH_ev_vec = np.linspace(3,3,1)
# eff_rec_vec = np.linspace(0, 1, 11)

PP_ev_vec = np.linspace(3,3,1)
PP_cd_vec = np.linspace(3,3,1)

# Instanciate Circuit
fluid = "Cyclopentane"

tries = 0
successes = 0
failures = 0

success_time = 0

# for eff_rec_val in eff_rec_vec:
for PP_ev in PP_ev_vec:
    for PP_cd in PP_cd_vec:
        for SC_cd in SC_cd_vec:
            for SH_ev in SH_ev_vec:
                tries += 1
        
                HP = IterativeCircuit(fluid)
                 
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
                Pinch_cd = 10  # K
                SC_cd = SC_cd # 3  # K
                DP_h_cd = 0
                
                # Evaporator
                Pinch_ev = 3  # K
                SH_ev = SH_ev # 5 # K
                DP_c_ev = 0
        
                # Recuperator
                eff_rec = 0.8 # 0.8 -
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
                HP.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="Condenser:su_C")
                HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
                HP.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="Evaporator:su_H")
                                
                #%% Inputs
                m_dot_ref = 20  # kg/s
                
                #%% Cycle guess values
                
                # Set iteration variables
                P_HP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
                P_LP_guess = PropsSI("P", "T", T_su_w_ev-10, "Q", 1, fluid)
                h_su_vlv_guess = PropsSI('H', 'P', P_HP_guess, 'Q', 0, fluid)
                # h_su_vlv_guess = PropsSI('H', 'P', P_HP_guess, 'T', T_su_vlv_guess, fluid)
                
                HP.set_cycle_input(target="ExpansionValve:su", m_dot = m_dot_ref)
                HP.set_cycle_input(target="Compressor:su", m_dot = m_dot_ref)
                HP.set_cycle_guess(target="Compressor:su", SH=10)     
                
                HP.set_iteration_variable(
                    target=["ExpansionValve:su", "Compressor:ex"],
                    variable="p",
                    guess=P_HP_guess,
                    tolerance=1e-6
                )
                
                HP.set_iteration_variable(
                    target=["ExpansionValve:ex", "Compressor:su"],
                    variable="p",
                    guess=P_LP_guess,
                    tolerance=1e-6
                )
                
                HP.set_iteration_variable(
                    target="ExpansionValve:su",
                    variable="h",
                    guess=h_su_vlv_guess,
                    tolerance=1e-6
                )
                
                #%%
                
                # Fermeture croisée — p HP : calculée sur Condenser:ex_H doit rejoindre le guess sur ExpansionValve:su
                HP.set_residual_variable(
                    target="ExpansionValve:su-p",
                    tolerance=1e-3
                )
                
                HP.set_residual_variable(
                    target="Compressor:su-p",
                    tolerance=1e-3
                )
                
                HP.set_residual_variable(
                    target="ExpansionValve:su-h",
                    tolerance=1e-3
                )
                
                start = time.perf_counter()
                HP.solve(method='lm')
                end = time.perf_counter()
        
                elapsed = end - start
        
                if HP.converged:
                    # print(f"Converged in {HP.n_it} iterations !") # res_energy : {HP.res_energy}")
                    successes += 1
                    success_time += elapsed
                else:
                    print(f"Failure... (conv)")
                    
                    failures += 1
        
                if np.mod(tries,100) == 0:
                    print(f"Tries : {tries}")

avg_success_time = success_time/successes
conv_prop = successes/tries

print(f"avg_success_time : {avg_success_time}")
print(f"conv_prop : {conv_prop}")    
        
# HP.plot_cycle_Ts() 

# print(f"Success : {successes}/{tries} {round(successes/tries * 100,2)} %")
# print(f"Success Avg Time : {round(success_time/(successes+1e-14),5)} s")
   