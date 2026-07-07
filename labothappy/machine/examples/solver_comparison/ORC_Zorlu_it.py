from labothappy.machine.circuit_it import IterativeCircuit

from labothappy.connector.mass_connector import MassConnector

from labothappy.component.pump.pump_csteff import PumpCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

from CoolProp.CoolProp import PropsSI

import time
import numpy as np

#%% Instanciate Circuit

# T_guess_cd = np.linspace(-5,5,11) + 24+273.15
# T_guess_ev = np.linspace(-5,5,11) + 141+273.15

T_guess_cd = [24+273.15]
T_guess_ev = [141+273.15] 

# T_guess_cd = [24+273.15]
# T_guess_ev = [141+273.15] 

# SC_cd_vec = np.linspace(1,1,1)
# SH_ev_vec = np.linspace(1,10,10)

# PP_ev_vec = np.linspace(1,10,10)
# PP_cd_vec = np.linspace(1,10,10)

SC_cd_vec = np.linspace(1,1,1)
SH_ev_vec = np.linspace(3,3,1)

PP_ev_vec = np.linspace(4,4,1)
PP_cd_vec = np.linspace(10,10,1)

# PP_ev : 2.0
# PP_cd : 6.0
# SH_ev : 2.0
# SC_cd : 3.0

# Instanciate Circuit
fluid = "Cyclopentane"

successes = 0
failures = 0

success_time = 0
tries = 0

for PP_ev in PP_ev_vec:
    for PP_cd in PP_cd_vec:
        for T_cd in T_guess_cd:
            for T_ev in T_guess_ev:
                for SC_cd in SC_cd_vec:
                    for SH_ev in SH_ev_vec:
                        
                        # print(f"="*20)

                        # print(f"PP_ev : {PP_ev}")
                        # print(f"PP_cd : {PP_cd}")
                        # print(f"SC_cd : {SC_cd}")
                        # print(f"SH_ev : {SH_ev}")
                        
                        tries += 1
                        
                        ORC = IterativeCircuit(fluid)
        
                        # Ignore debug printing
                        ORC.mute_print()
                        
                        # Create components
                        Pump = PumpCstEff()
                        Condenser = HexCstPinch()
                        Expander = ExpanderCstEff()
                        Evaporator = HexCstPinch()
                        Recuperator = HexCstEff()
                        Preheater = HexCstEff()
                        
                        #%% Set component parameters
                        eta_is_pp = 0.7
                        
                        eta_is_exp = 0.8
                        
                        Pinch_cd = PP_cd  # K
                        SC_cd = SC_cd # 1 # K
                        
                        Pinch_ev = PP_ev  # K
                        SH_ev = SH_ev # 3 # K
                        
                        eff_rec = 0.8 
                        eff_pre = 0.8
                        
                        Pump.set_parameters(eta_is=eta_is_pp)
                        
                        Condenser.set_parameters(
                            Pinch=Pinch_cd, 
                            Delta_T_sh_sc=SC_cd, 
                            HX_type="condenser", 
                            DP_c = 0,
                            DP_h = 0*1e3)
                        
                        Expander.set_parameters(eta_is=eta_is_exp)
                        
                        Evaporator.set_parameters(
                            Pinch=Pinch_ev, 
                            Delta_T_sh_sc=SH_ev, 
                            HX_type="evaporator",
                            DP_c = 0*1e3,
                            DP_h = 0)
                        
                        Recuperator.set_parameters(
                            eta=eff_rec,
                            DP_c = 0*1e3,
                            DP_h = 0*1e3)
                        
                        Preheater.set_parameters(
                            eta=eff_pre,
                            DP_c = 0*1e3,
                            DP_h = 0)
                        
                        #%% Add components to circuit
                        ORC.add_component(Pump, "Pump")
                        ORC.add_component(Condenser, "Condenser")
                        ORC.add_component(Expander, "Expander")
                        ORC.add_component(Evaporator, "Evaporator")
                        ORC.add_component(Recuperator, "Recuperator")
                        ORC.add_component(Preheater, "Preheater")
                        
                        # Link components with mass connectors
                        ORC.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
                        ORC.link_components("Recuperator", "m-ex_C", "Preheater", "m-su_C")
                        ORC.link_components("Preheater", "m-ex_C", "Evaporator", "m-su_C")
                        ORC.link_components("Evaporator", "m-ex_C", "Expander", "m-su")
                        ORC.link_components("Expander", "m-ex", "Recuperator", "m-su_H")
                        ORC.link_components("Recuperator", "m-ex_H", "Condenser", "m-su_H")
                        ORC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
                        
                        #%% Add fluid sources
                        CD_source = MassConnector('Water')
                        T_su_w_cd = 24+273.15
                        P_su_w_cd = 2e5
                        m_dot_w_cd = 900 # kg/s
                        
                        EV_source = MassConnector('Water')
                        T_su_w_ev = 141+273.15
                        P_su_w_ev = 10e5
                        m_dot_w_ev = 10000  # kg/s
                        
                        PRE_source = MassConnector('Water')
                        T_su_w_pre = 113.1+273.15
                        P_su_w_pre = 2e5
                        m_dot_w_pre = 60  # kg/s
                        
                        ORC.add_source("CD_Water", CD_source, ORC.components["Condenser"], "m-su_C")
                        ORC.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="Condenser:su_C")
                        
                        ORC.add_source("EV_Water", EV_source, ORC.components["Evaporator"], "m-su_H")
                        ORC.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="Evaporator:su_H")
                        
                        ORC.add_source("PRE_Water", PRE_source, ORC.components["Preheater"], "m-su_H")
                        ORC.set_source_properties(T=T_su_w_pre, fluid='Water', P=P_su_w_pre, m_dot = m_dot_w_pre, target="Preheater:su_H")
                        
                        #%% Inputs
                        m_dot_ref = 34.51 # kg/s
                        
                        #%% Cycle guess values
                        T_sat_guess_cd = T_cd + 3 + 5
                        T_sat_guess_ev = T_ev - 3 - 10
                        
                        P_low = PropsSI("P", "T", T_sat_guess_cd, "Q", 1, fluid)
                        P_high = PropsSI("P", "T", T_sat_guess_ev, "Q", 0, fluid)
                        
                        h_pp_guess = PropsSI("H", "T", T_sat_guess_cd-SC_cd, "P", P_low, fluid)
                        h_exp_guess = PropsSI("H", "T", T_sat_guess_ev+SH_ev, "P", P_high, fluid)
        
                        #%%
                                        
                        ORC.set_cycle_guess(target="Pump:su", m_dot = m_dot_ref, h=h_pp_guess, p=P_low)
                        ORC.set_cycle_guess(target="Pump:ex", p=P_high)
                                        
                        # ORC.set_cycle_guess(target="Recuperator:su_C", p=P_high, m_dot=m_dot_ref, T=T_su_w_cd+Pinch_cd)
                        
                        ORC.set_cycle_guess(target="Expander:su", p=P_high, m_dot = m_dot_ref, SH=SH_ev)
                        ORC.set_cycle_guess(target="Expander:ex", p=P_low)
                        
                        #%% CYCLE FIXED VARIABLES AND ITERATION VARIABLE
                        
                        ORC.set_iteration_variable(
                            target  = 'Expander:ex',
                            variable = 'p',
                            guess = P_low,
                            tolerance=1e-6
                            )
                                        
                        ORC.set_iteration_variable(
                            target  = ['Pump:ex', 'Expander:su'],
                            variable = 'p',
                            guess=P_high,
                            tolerance=1e-6
                            )
                        
                        ORC.set_iteration_variable(
                            target="Expander:su",
                            variable="h",
                            guess=h_exp_guess,
                            tolerance=1e-6
                            )
                        
                        #%% CYCLE FIXED VARIABLES AND ITERATION VARIABLE
        
                        ORC.set_residual_variable(
                            target="Expander:su-p",
                            tolerance=1e-3
                        )
                        
                        ORC.set_residual_variable(
                            target="Pump:su-p",
                            tolerance=1e-3
                        )
                        
                        ORC.set_residual_variable(
                            target="Expander:su-h",
                            tolerance=1e-3
                        )
                        
                        start = time.perf_counter()
                        ORC.solve(method='newton')
                        end = time.perf_counter()
        
                        elapsed = end - start
        
                        if ORC.converged:
                            # print(f"Success !")
                            successes += 1
                            success_time += elapsed
                        else:
                            # print(f"Failure... (conv)")
                            print(f"="*20)
                            print(f"PP_ev : {PP_ev}")
                            print(f"PP_cd : {PP_cd}")
                            print(f"SH_ev : {SH_ev}")
                            print(f"SC_cd : {SC_cd}")     
                            
                        if np.mod(tries,100) == 0:
                            print(f"Tries : {tries}")
                
                
avg_success_time = success_time/(successes+1e-9)
conv_prop = successes/tries

print(f"avg_success_time : {avg_success_time}")
print(f"conv_prop : {conv_prop}")    

# # print(f"Success : {successes}/{tries} {round(successes/tries * 100,2)} %")
# # print(f"Success Avg Time : {round(success_time/(successes+1e-8),5)} s")                
# # ORC.plot_cycle_Ts()

# # import matplotlib.pyplot as plt
# # plt.show()

# # # print(f"Converged at P_HP = {Pump.ex.p}, P_LP = {Pump.su.p}")

# # print("Pump")
# # print(Pump.W.W_dot)
# # print((Pump.ex.h - Pump.su.h)*Pump.su.m_dot)

# from labothappy.machine.circuit_it import IterativeCircuit
# from labothappy.connector.mass_connector import MassConnector
# from labothappy.component.pump.pump_csteff import PumpCstEff
# from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
# from labothappy.component.expander.expander_csteff import ExpanderCstEff
# from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
# from CoolProp.CoolProp import PropsSI

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
# import time

# # ── Parameter sweep ────────────────────────────────────────────────────────────
# T_guess_cd = [24  + 273.15]
# T_guess_ev = [141 + 273.15]

# PP_ev_vec  = np.linspace(1, 10, 10)
# PP_cd_vec  = np.linspace(1, 10, 10)
# SH_ev_vec  = np.linspace(1, 10, 10)
# SC_cd_vec  = np.linspace(1, 10, 10)

# fluid      = "Cyclopentane"

# # Result array shape: (PP_ev, PP_cd, SH_ev, SC_cd)
# n = [len(PP_ev_vec), len(PP_cd_vec), len(SH_ev_vec), len(SC_cd_vec)]
# converged_arr = np.zeros(n, dtype=bool)

# labels     = ["PP_ev (K)", "PP_cd (K)", "SH_ev (K)", "SC_cd (K)"]
# vecs       = [PP_ev_vec, PP_cd_vec, SH_ev_vec, SC_cd_vec]

# total = np.prod(n)
# done  = 0

# # ── Sweep ──────────────────────────────────────────────────────────────────────
# for i0, PP_ev in enumerate(PP_ev_vec):
#     for i1, PP_cd in enumerate(PP_cd_vec):
#         for i2, SH_ev in enumerate(SH_ev_vec):
#             for i3, SC_cd in enumerate(SC_cd_vec):

#                 done += 1
#                 if done % 500 == 0:
#                     print(f"Progress: {done}/{total}")

#                 for T_cd in T_guess_cd:
#                     for T_ev in T_guess_ev:

#                         ORC = IterativeCircuit(fluid)
#                         ORC.mute_print()

#                         Pump        = PumpCstEff()
#                         Condenser   = HexCstPinch()
#                         Expander    = ExpanderCstEff()
#                         Evaporator  = HexCstPinch()
#                         Recuperator = HexCstEff()
#                         Preheater   = HexCstEff()

#                         Pump.set_parameters(eta_is=0.7)
#                         Condenser.set_parameters(Pinch=PP_cd, Delta_T_sh_sc=SC_cd,
#                                                  HX_type="condenser", DP_c=0, DP_h=0)
#                         Expander.set_parameters(eta_is=0.8)
#                         Evaporator.set_parameters(Pinch=PP_ev, Delta_T_sh_sc=SH_ev,
#                                                   HX_type="evaporator", DP_c=0, DP_h=0)
#                         Recuperator.set_parameters(eta=0.8, DP_c=0, DP_h=0)
#                         Preheater.set_parameters(eta=0.8, DP_c=0, DP_h=0)

#                         ORC.add_component(Pump,        "Pump")
#                         ORC.add_component(Condenser,   "Condenser")
#                         ORC.add_component(Expander,    "Expander")
#                         ORC.add_component(Evaporator,  "Evaporator")
#                         ORC.add_component(Recuperator, "Recuperator")
#                         ORC.add_component(Preheater,   "Preheater")

#                         ORC.link_components("Pump",        "m-ex",    "Recuperator", "m-su_C")
#                         ORC.link_components("Recuperator", "m-ex_C",  "Preheater",   "m-su_C")
#                         ORC.link_components("Preheater",   "m-ex_C",  "Evaporator",  "m-su_C")
#                         ORC.link_components("Evaporator",  "m-ex_C",  "Expander",    "m-su")
#                         ORC.link_components("Expander",    "m-ex",    "Recuperator", "m-su_H")
#                         ORC.link_components("Recuperator", "m-ex_H",  "Condenser",   "m-su_H")
#                         ORC.link_components("Condenser",   "m-ex_H",  "Pump",        "m-su")

#                         CD_source = MassConnector('Water')
#                         EV_source = MassConnector('Water')
#                         PRE_source = MassConnector('Water')

#                         ORC.add_source("CD_Water",  CD_source,  ORC.components["Condenser"],  "m-su_C")
#                         ORC.set_source_properties(T=24+273.15,    fluid='Water', P=2e5,   m_dot=900,   target="Condenser:su_C")
#                         ORC.add_source("EV_Water",  EV_source,  ORC.components["Evaporator"], "m-su_H")
#                         ORC.set_source_properties(T=141+273.15,   fluid='Water', P=10e5,  m_dot=10000, target="Evaporator:su_H")
#                         ORC.add_source("PRE_Water", PRE_source, ORC.components["Preheater"],  "m-su_H")
#                         ORC.set_source_properties(T=113.1+273.15, fluid='Water', P=2e5,   m_dot=60,    target="Preheater:su_H")

#                         m_dot_ref = 34.51
#                         P_low  = PropsSI("P", "T", T_cd,         "Q", 1, fluid)
#                         P_high = PropsSI("P", "T", T_ev,         "Q", 0, fluid)
#                         h_pp   = PropsSI("H", "T", T_cd - SC_cd, "P", P_low,  fluid)
#                         h_exp  = PropsSI("H", "T", T_ev + 2*SH_ev, "P", P_high, fluid)

#                         ORC.set_cycle_guess(target="Pump:su",     m_dot=m_dot_ref, h=h_pp,  p=P_low)
#                         ORC.set_cycle_guess(target="Pump:ex",     p=P_high)
#                         ORC.set_cycle_guess(target="Expander:su", p=P_high, m_dot=m_dot_ref, SH=SH_ev*2)
#                         ORC.set_cycle_guess(target="Expander:ex", p=P_low)

#                         ORC.set_iteration_variable(target='Expander:ex',              variable='p', guess=P_low,  tolerance=1e-6)
#                         ORC.set_iteration_variable(target=['Pump:ex', 'Expander:su'], variable='p', guess=P_high, tolerance=1e-6)
#                         ORC.set_iteration_variable(target="Expander:su",              variable='h', guess=h_exp,  tolerance=1e-6)

#                         ORC.set_residual_variable(target="Expander:su-p", tolerance=1e-3)
#                         ORC.set_residual_variable(target="Pump:su-p",     tolerance=1e-3)
#                         ORC.set_residual_variable(target="Expander:su-h", tolerance=1e-3)

#                         ORC.solve(method='fsolve')
#                         converged_arr[i0, i1, i2, i3] = ORC.converged

# # ── Save raw results ───────────────────────────────────────────────────────────
# np.save("convergence_results.npy", converged_arr)
# print("Results saved to convergence_results.npy")

# # ── Plot: 6 convergence-rate maps ─────────────────────────────────────────────
# from itertools import combinations

# axis_pairs = list(combinations(range(4), 2))   # (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)

# fig, axes = plt.subplots(2, 3, figsize=(14, 9))
# axes = axes.flatten()

# cmap = plt.cm.RdYlGn   # red = never converges, green = always converges

# for ax, (a, b) in zip(axes, axis_pairs):
#     # Average convergence rate over the two remaining axes
#     other_axes = tuple(i for i in range(4) if i not in (a, b))
#     rate = converged_arr.mean(axis=other_axes)   # shape (n_a, n_b)

#     # Orient so X=b, Y=a (first index on Y axis)
#     im = ax.imshow(
#         rate,
#         origin='lower',
#         aspect='auto',
#         cmap=cmap,
#         vmin=0, vmax=1,
#         extent=[vecs[b][0], vecs[b][-1], vecs[a][0], vecs[a][-1]]
#     )

#     ax.set_xlabel(labels[b], fontsize=11)
#     ax.set_ylabel(labels[a], fontsize=11)
#     ax.set_title(f"{labels[a].split()[0]} vs {labels[b].split()[0]}", fontsize=12, fontweight='bold')

#     cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
#     cb.set_label("Convergence rate", fontsize=9)
#     cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
#     cb.set_ticklabels(["0%", "25%", "50%", "75%", "100%"])

# fig.suptitle("ORC convergence rate maps (fsolve)\naveraged over remaining parameters",
#              fontsize=14, fontweight='bold', y=1.01)

# plt.tight_layout()
# plt.savefig("convergence_maps.pdf", bbox_inches='tight', dpi=150)
# plt.savefig("convergence_maps.png", bbox_inches='tight', dpi=150)
# plt.show()
# print("Figures saved.")

import matplotlib.pyplot as plt

P_evap = []

for x in ORC.x_log:
    P_evap.append(x[1]*1e-5)
    
plt.plot(P_evap)
plt.ylabel("$P_{evap}$ [bar]", fontsize=12)
plt.xlabel("Iteration [-]", fontsize=12)
plt.plot([0,len(P_evap)], [6.938, 6.938])

plt.grid()

plt.show()
    
