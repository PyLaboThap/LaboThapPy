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

# T_guess_ev = np.linspace(120,160,21) + 273.15
# T_guess_cd = np.linspace(100,130,16) + 273.15

T_guess_cd = np.linspace(-5,5,11) + 24+273.15
T_guess_ev = np.linspace(-5,5,11) + 141+273.15

SC_cd_vec = np.linspace(1,10,10)
SH_ev_vec = np.linspace(1,10,10)

# Instanciate Circuit
fluid = "Cyclopentane"

successes = 0
failures = 0

success_time = 0
tries = 0

for T_cd in T_guess_cd:
    for T_ev in T_guess_ev:
        for SC_cd in SC_cd_vec:
            for SH_ev in SH_ev_vec:
                
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
                
                Pinch_cd = 3  # K
                SC_cd = SC_cd # 1 # K
                
                Pinch_ev = 10  # K
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
                T_sat_guess_cd = T_cd
                T_sat_guess_ev = T_ev
                
                P_low = PropsSI("P", "T", T_sat_guess_cd, "Q", 1, fluid)
                P_high = PropsSI("P", "T", T_sat_guess_ev, "Q", 0, fluid)
                
                h_pp_guess = PropsSI("H", "T", T_sat_guess_cd-SC_cd, "P", P_low, fluid)
                h_exp_guess = PropsSI("H", "T", T_sat_guess_ev+2*SH_ev, "P", P_high, fluid)

                #%%
                                
                ORC.set_cycle_guess(target="Pump:su", m_dot = m_dot_ref, h=h_pp_guess, p=P_low)
                ORC.set_cycle_guess(target="Pump:ex", p=P_high)
                                
                # ORC.set_cycle_guess(target="Recuperator:su_C", p=P_high, m_dot=m_dot_ref, T=T_su_w_cd+Pinch_cd)
                
                ORC.set_cycle_guess(target="Expander:su", p=P_high, m_dot = m_dot_ref, SH=SH_ev*2)
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
                    pass
                    # print(f"Failure... (conv)")
                
                if np.mod(tries,100) == 0:
                    print(f"Tries : {tries}")
                
                
# print(f"Success : {successes}/{tries} {round(successes/tries * 100,2)} %")
# print(f"Success Avg Time : {round(success_time/(successes+1e-8),5)} s")                
# ORC.plot_cycle_Ts()

# import matplotlib.pyplot as plt
# plt.show()

# # print(f"Converged at P_HP = {Pump.ex.p}, P_LP = {Pump.su.p}")

# print("Pump")
# print(Pump.W.W_dot)
# print((Pump.ex.h - Pump.su.h)*Pump.su.m_dot)

