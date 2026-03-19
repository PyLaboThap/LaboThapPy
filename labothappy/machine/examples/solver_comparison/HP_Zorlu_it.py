from labothappy.machine.circuit_it import IterativeCircuit

from labothappy.connector.mass_connector import MassConnector

from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

from CoolProp.CoolProp import PropsSI

import time
import numpy as np

# T_guess_cd = np.linspace(120,160,21) + 273.15
# T_guess_ev = np.linspace(100,130,16) + 273.15

T_guess_cd = [141   + 2 + 273.15]
T_guess_ev = [113.1 - 2 + 273.15]

# SC_cd_vec = np.linspace(1,10,10)
# SH_ev_vec = np.linspace(1,10,10)

SC_cd_vec = np.linspace(3,3,1)
SH_ev_vec = np.linspace(1,1,1)

# Instanciate Circuit
fluid = "Cyclopentane"

tries = len(T_guess_cd)*len(T_guess_ev)*len(SC_cd_vec)*len(SH_ev_vec)
successes = 0
failures = 0

success_time = 0

for T_cd in T_guess_cd:
    for T_ev in T_guess_ev:
        for SC_cd in SC_cd_vec:
            for SH_ev in SH_ev_vec:
                HP = IterativeCircuit(fluid)
                 
                # Ignore debug printing
                # HP.mute_print()
                HP.mute_plot()
                
                # Create components
                Compressor = CompressorCstEff()
                Condenser = HexCstPinch()
                ExpansionValve = ValveIsenthalpic()
                Evaporator = HexCstPinch()
                Recuperator = HexCstEff()
                
                # Set component parameters
                eta_is_cp = 0.8 # -
                
                Pinch_cd = 10  # K
                SC_cd = SC_cd # 3  # K
                
                Pinch_ev = 3  # K
                SH_ev = SH_ev # 5 # K
                
                eff_rec = 0.8 # -
                
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
                    DP_h = 0*1e3)
                
                Evaporator.set_parameters(
                    Pinch=Pinch_ev,
                    Delta_T_sh_sc=SH_ev, 
                    HX_type="evaporator",
                    DP_c = 0*1e3,
                    DP_h = 0*1e3)
                
                Recuperator.set_parameters(
                    eta=eff_rec,
                    DP_c = 0*1e3,
                    DP_h = 0*1e3)
                
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
                m_dot_ref = 19.31* 1.0352594687116445  # kg/s
                
                #%% Cycle guess values
                
                # Set iteration variables
                P_HP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
                P_LP_guess = PropsSI("P", "T", T_su_w_ev-10, "Q", 1, fluid)
                h_su_vlv_guess = PropsSI('H', 'P', P_HP_guess, 'Q', 0, fluid) - SC_cd - 0.1
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
                HP.solve()
                end = time.perf_counter()

                elapsed = end - start

                if HP.converged:
                    print(f"Success in {HP.n_it} iteration(s) !")
                    successes += 1
                    success_time += elapsed
                else:
                    print(f"Failure... (conv)")
                    failures += 1
                    
                # except:
                #     print(f"Failure...")
                
    HP.plot_cycle_Ts() 
    
    # HP.Ts_gif()
    
print(f"Success : {successes}/{tries} {round(successes/tries * 100,2)} %")
print(f"Success Avg Time : {round(success_time/(successes+1e-14),5)} s")
   
# from labothappy.machine.circuit_it import IterativeCircuit
# from labothappy.connector.mass_connector import MassConnector
# from labothappy.component.compressor.compressor_csteff import CompressorCstEff
# from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
# from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic
# from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

# from CoolProp.CoolProp import PropsSI

# # Instanciate Circuit
# fluid = "Cyclopentane"
# HP = IterativeCircuit(fluid)

# # Ignore debug printing
# # HP.mute_print()

# # Create components
# Compressor = CompressorCstEff()
# Condenser = HexCstPinch()
# ExpansionValve = ValveIsenthalpic()
# Evaporator = HexCstPinch()
# Recuperator = HexCstEff()

# # Set component parameters
# eta_is_cp = 0.8 # -

# Pinch_cd = 10  # K
# SC_cd = 3  # K

# Pinch_ev = 3  # K
# SH_ev = 5 # K

# eff_rec = 0.8 # -

# Compressor.set_parameters(
#     eta_is=eta_is_cp)

# Condenser.set_parameters(
#     Pinch=Pinch_cd, 
#     Delta_T_sh_sc=SC_cd, 
#     HX_type="condenser",
#     DP_c = 0*1e3,
#     DP_h = 0*1e3)

# Evaporator.set_parameters(
#     Pinch=Pinch_ev, 
#     Delta_T_sh_sc=SH_ev, 
#     HX_type="evaporator",
#     DP_c = 0*1e3,
#     DP_h = 0*1e3)

# Recuperator.set_parameters(
#     eta=eff_rec,
#     DP_c = 0*1e3,
#     DP_h = 0*1e3)

# # Add components to circuit
# HP.add_component(Compressor, "Compressor")
# HP.add_component(Condenser, "Condenser")
# HP.add_component(ExpansionValve, "ExpansionValve")
# HP.add_component(Evaporator, "Evaporator")
# HP.add_component(Recuperator, "Recuperator")

# # Link components with mass connectors
# HP.link_components("Compressor", "m-ex", "Condenser", "m-su_H")
# HP.link_components("Condenser", "m-ex_H", "Recuperator", "m-su_H")
# HP.link_components("Recuperator", "m-ex_H", "ExpansionValve", "m-su")
# HP.link_components("ExpansionValve", "m-ex", "Evaporator", "m-su_C")
# HP.link_components("Evaporator", "m-ex_C", "Recuperator", "m-su_C")
# HP.link_components("Recuperator", "m-ex_C", "Compressor", "m-su")

# # Add fluid sources
# CD_source = MassConnector('Water')
# T_su_w_cd = 141+273.15
# P_su_w_cd = 5e5
# m_dot_w_cd = 10000  # kg/s
# EV_source = MassConnector('Water')
# T_su_w_ev = 113.1+273.15
# P_su_w_ev = 189*1e3
# m_dot_w_ev = 497.6  # kg/s

# HP.add_source("CD_Water", CD_source, HP.components["Condenser"], "m-su_C")
# HP.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
# HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
# HP.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="EV_Water")

# #%% Inputs
# m_dot_ref = 19.31* 1.0352594687116445  # kg/s

# #%% Cycle guess values

# HP.set_cycle_input(target="ExpansionValve:su", m_dot = m_dot_ref)
# HP.set_cycle_input(target="Compressor:su", m_dot = m_dot_ref)
# HP.set_cycle_input(target="Evaporator:ex_C", SH=SH_ev)
# HP.set_cycle_input(target="Condenser:ex_H", m_dot = m_dot_ref, SC=SC_cd)

# # Set iteration variables
# P_HP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
# P_LP_guess = PropsSI("P", "T", T_su_w_ev-10, "Q", 1, fluid)
# h_su_vlv_guess = PropsSI('H', 'P', P_HP_guess, 'Q', 0, fluid) - SC_cd - 0.1
# # h_su_vlv_guess = PropsSI('H', 'P', P_HP_guess, 'T', T_su_vlv_guess, fluid)

# HP.set_iteration_variable(
#     target=["ExpansionValve:su", "Compressor:ex", "Condenser:ex_H"],
#     variable="p",
#     guess=P_HP_guess,
#     tolerance=1e-6
# )

# HP.set_iteration_variable(
#     target="ExpansionValve:ex",
#     variable="p",
#     guess=P_LP_guess,
#     tolerance=1e-6
# )

# HP.set_iteration_variable(
#     target="ExpansionValve:su",
#     variable="h",
#     guess=h_su_vlv_guess,
#     tolerance=1e-6
# )

# # Set residual variables
# HP.set_residual_variable(
#     pre_target="ExpansionValve:su",
#     post_target="Condenser:ex_H",
#     variable="p",
#     tolerance=1e-3
# )

# HP.set_residual_variable(
#     pre_target="ExpansionValve:ex",
#     post_target="Evaporator:su_C",
#     variable="p",
#     tolerance=1e-3
# )

# HP.set_residual_variable(
#     pre_target="ExpansionValve:su",
#     post_target="ExpansionValve:su",
#     variable="h",
#     tolerance=1e-3
# )

# # Solve circuit
# HP.solve()

# HP.plot_cycle_Ts()

# print(f"Converged at P_HP = {Compressor.ex.p}, P_LP = {Compressor.su.p}")
