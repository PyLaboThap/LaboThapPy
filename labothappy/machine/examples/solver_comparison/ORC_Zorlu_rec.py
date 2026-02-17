from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.pump.pump_csteff import PumpCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

from CoolProp.CoolProp import PropsSI

import time
import numpy as np

#%% Instanciate Circuit

T_guess_ev = np.linspace(120,160,21) + 273.15
T_guess_cd = np.linspace(100,130,16) + 273.15

SC_cd_vec = np.linspace(1,1,1)
SH_ev_vec = np.linspace(3,3,1)

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

                ORC = RecursiveCircuit(fluid)

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
                
                Pinch_cd = 5  # K
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
                    DP_h = 15*1e3)
                
                Expander.set_parameters(eta_is=eta_is_exp)
                
                Evaporator.set_parameters(
                    Pinch=Pinch_ev, 
                    Delta_T_sh_sc=SH_ev, 
                    HX_type="evaporator",
                    DP_c = 50*1e3,
                    DP_h = 0)
                
                Recuperator.set_parameters(
                    eta=eff_rec,
                    DP_c = 50*1e3,
                    DP_h = 10*1e3)
                
                Preheater.set_parameters(
                    eta=eff_pre,
                    DP_c = 50*1e3,
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
                P_low = PropsSI("P", "T", T_cd+Pinch_cd+SC_cd, "Q", 1, fluid)
                P_high = PropsSI("P", "T", T_ev-Pinch_ev-SH_ev, "Q", 0, fluid)
                
                #%%
                
                ORC.set_cycle_guess(target="Pump:su", m_dot = m_dot_ref, SC=SC_cd*2, p=P_low)
                ORC.set_cycle_guess(target="Pump:ex", p=P_high)
                
                # ORC.set_cycle_guess(target="Recuperator:su_C", p=P_high, m_dot=m_dot_ref, T=T_su_w_cd+Pinch_cd)
                
                ORC.set_cycle_guess(target="Expander:su", p=P_high, m_dot = m_dot_ref, SH=SH_ev*2)
                ORC.set_cycle_guess(target="Expander:ex", p=P_low)
                
                #%% Set residual variables
                ORC.set_residual_variable(target="Pump:ex", variable="h", tolerance=1e-3)
                ORC.set_residual_variable(target="Pump:ex", variable="p", tolerance=1e-3)
                
                ORC.set_residual_variable(target="Condenser:ex_H", variable="h", tolerance=1e-3)
                ORC.set_residual_variable(target="Condenser:ex_H", variable="p", tolerance=1e-3)
                
                ORC.set_residual_variable(target="Evaporator:ex_C", variable="h", tolerance=1e-3)
                ORC.set_residual_variable(target="Evaporator:ex_C", variable="p", tolerance=1e-3)
                
                ORC.set_residual_variable(target="Expander:ex", variable="h", tolerance=1e-3)
                ORC.set_residual_variable(target="Expander:ex", variable="p", tolerance=1e-3)
                
                #%% CYCLE FIXED VARIABLES AND ITERATION VARIABLE
                damping = 0.2 # 0.2
                
                ORC.set_iteration_variable(target=['Expander:ex'], variable='p', objective = 'Link:Condenser:su_H-p', tol = 1e-2, rel = 1, damping_factor = damping, cycle = ORC)
                ORC.set_iteration_variable(target=['Pump:ex'], variable='p', objective = 'Link:Evaporator:su_C-p', tol = 1e-2, rel = 1, damping_factor = damping, cycle = ORC)
                
                try: 
                    start = time.perf_counter()
                    ORC.solve()                    
                    end = time.perf_counter()

                    elapsed = end - start

                    if ORC.converged:
                        print(f"Success !")
                        successes += 1
                        success_time += elapsed
                    else:
                        print(f"Failure... (conv)")
                        failures += 1
                except:
                    print(f"Failure...")
                
    # HP.plot_cycle_Ts() 
    
print(f"Success : {successes}/{tries} {round(successes/tries * 100,2)} %")
print(f"Success Avg Time : {round(success_time/successes,5)} s")                
# ORC.plot_cycle_Ts()

# print(f"Converged at P_HP = {Pump.ex.p}, P_LP = {Pump.su.p}")
