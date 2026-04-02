from labothappy.machine.circuit import Circuit

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

SC_cd_vec = np.linspace(1,10,10)
SH_ev_vec = np.linspace(1,10,10)

PP_ev_vec = np.linspace(1,10,10)
PP_cd_vec = np.linspace(1,10,10)

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
                        
                        tries += 1
                        
                        ORC = Circuit(fluid)
        
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
                        ORC.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
                        
                        ORC.add_source("EV_Water", EV_source, ORC.components["Evaporator"], "m-su_H")
                        ORC.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="EV_Water")
                        
                        ORC.add_source("PRE_Water", PRE_source, ORC.components["Preheater"], "m-su_H")
                        ORC.set_source_properties(T=T_su_w_pre, fluid='Water', P=P_su_w_pre, m_dot = m_dot_w_pre, target="PRE_Water")
                        
                        #%% Inputs
                        m_dot_ref = 34.51 # kg/s
                        
                        #%% Cycle guess values
                        T_cd_guess = T_cd
                        T_ev_guess = T_ev
                        
                        P_low = PropsSI("P", "T", T_cd_guess, "Q", 1, fluid)
                        P_high = PropsSI("P", "T", T_ev_guess, "Q", 0, fluid)
                        
                        h_pp_guess = PropsSI("H", "T", T_cd_guess-SC_cd, "P", P_low, fluid)
        
                        #%%
                        
                        ORC.set_cycle_guess(target="Pump:su", m_dot = m_dot_ref, h=h_pp_guess, p=P_low)
                        ORC.set_cycle_guess(target="Pump:ex", p=P_high)
                        
                        # ORC.set_cycle_guess(target="Recuperator:su_C", p=P_high, m_dot=m_dot_ref, T=T_su_w_cd+Pinch_cd)
                        
                        ORC.set_cycle_guess(target="Expander:su", p=P_high, m_dot = m_dot_ref, SH=SH_ev*2)
                        ORC.set_cycle_guess(target="Expander:ex", p=P_low)
                        
                        #%% CYCLE FIXED VARIABLES AND ITERATION VARIABLE
                                        
                        ORC.set_iteration_variable(
                            it_var  = 'Expander:ex-p',
                            objective = 'Condenser:ex_H-p',
                            obj_type = "Link"
                        )
                                        
                        ORC.set_iteration_variable(
                            it_var  = 'Pump:ex-p',
                            objective = 'Evaporator:ex_C-p',
                            obj_type = "Link"
                        )
                        
                        start = time.perf_counter()
                        ORC.solve(max_iter=100, method='wegstein')
                        end = time.perf_counter()
        
                        elapsed = end - start
        
                        if ORC.converged:
                            # print(f"Success in {ORC._iteration_count} iterations !")
                            successes += 1
                            success_time += elapsed
                        else:
                            print(f"Failure... (conv)")
                        
                        if np.mod(tries,100) == 0:
                            print(f"Tries : {tries}")
          

avg_success_time = success_time/successes
conv_prop = successes/tries

print(f"avg_success_time : {avg_success_time}")
print(f"conv_prop : {conv_prop}")          
          
# print(f"Success : {successes}/{tries} {round(successes/tries * 100,2)} %")
# print(f"Success Avg Time : {round(success_time/(successes+1e-8),5)} s")                
# ORC.plot_cycle_Ts()

# import matplotlib.pyplot as plt
# plt.show()

# # print(f"Converged at P_HP = {Pump.ex.p}, P_LP = {Pump.su.p}")

# print("Pump")
# print(Pump.W.W_dot)
# print((Pump.ex.h - Pump.su.h)*Pump.su.m_dot)

#%%

# import matplotlib.pyplot as plt
# import numpy as np

# methods = ["FPI", "Wegstein", "Newton", "Powell", "lm", "CoolSolve"]

# case_study = "HP"

# if case_study == "HP":
#     convergence = [100, 100, 100, 100, 100, 100]
#     avg_time = [35, 11, 70, 35, 57, 40]

# elif case_study == "ORC":
#     convergence = [100, 100, 99, 86, 95, 51]
#     avg_time = [17, 15, 56, 43, 62, 40]

# x = np.arange(len(methods))
# width = 0.35

# fontsize = 14

# fig, ax1 = plt.subplots(figsize=(7, 3))

# bars1 = ax1.bar(x - width / 2, convergence, width, label="Convergence [%]", color="#378ADD")
# ax1.set_ylabel("Convergence [%]", fontsize=12, color="#378ADD")
# ax1.tick_params(axis="y", labelcolor="#378ADD")
# ax1.set_ylim(0, 110)
# ax1.set_xticks(x)
# ax1.set_xticklabels(methods, fontsize=11)

# ax2 = ax1.twinx()
# bars2 = ax2.bar(x + width / 2, avg_time, width, label="Average Time [ms]", color="#E8832A")
# ax2.set_ylabel("Average Time [ms]", fontsize=12, color="#E8832A")
# ax2.tick_params(axis="y", labelcolor="#E8832A")
# ax2.set_ylim(0, 110)

# for bar in bars1:
#     ax1.text(
#         bar.get_x() + bar.get_width() / 2,
#         bar.get_height() + 1.5,
#         f"{int(bar.get_height())}",
#         ha="center", va="bottom", fontsize=9, color="#378ADD"
#     )

# for bar in bars2:
#     ax2.text(
#         bar.get_x() + bar.get_width() / 2,
#         bar.get_height() + 1,
#         f"{int(bar.get_height())}",
#         ha="center", va="bottom", fontsize=9, color="#E8832A"
#     )

# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.grid(axis="y", linestyle="--", alpha=0.3)

# plt.tight_layout()
# plt.savefig("/home/claude/table8_bar_chart.png", dpi=150)
# plt.show()
