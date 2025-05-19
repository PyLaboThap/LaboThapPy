# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:05:01 2025

@author: Basile
"""

from machine.circuit import Circuit
from CoolProp.CoolProp import PropsSI

from connector.mass_connector import MassConnector

from component.heat_exchanger.hex_MB_charge_sensitive import HeatExchangerMB
from component.heat_exchanger.steady_state.finite_volumes.cross_flow_tube_and_fins.simulation_model import CrossFlowTubeAndFinsHTX

from component.heat_exchanger.steady_state.finite_volumes.cross_flow_tube_and_fins.modules.geometry_tube_and_fins_hx import CFTubeAndFinsGeom
from toolbox.geometries.heat_exchanger.geometry_shell_and_tube_hx import ShellAndTubeGeom
from toolbox.geometries.heat_exchanger.geometry_tube_and_fins_hx import TubeAndFinsGeom

from component.turbomachinery.turbine.steady_state.polyn_isentropic_eff.simulation_model import Turb_polyn_eff
from component.turbomachinery.turbine.steady_state.polyn_isentropic_eff.modules.geometry_turb_polyn import Geometry_polyn_turb

from component.turbomachinery.pump.steady_state.extrapolation_model.simulation_model_P_out import PumpExtrapolationModelP_out
from component.turbomachinery.pump.steady_state.extrapolation_model.modules.geometry_extrapolation_pump import Geometry_extrapol_pump

from component.tank.spliter.simulation_model import Spliter
from component.tank.tank_mixer import Mixer 

#%% Define Circuit
ORC = Circuit('Cyclopentane')

#%% Create components
Turbine = Turb_polyn_eff()

Condenser_1 = CrossFlowTubeAndFinsHTX()
Condenser_2 = CrossFlowTubeAndFinsHTX()

Pump = PumpExtrapolationModelP_out()

Spliter_ACC = Spliter(outlet_repartition=[0.5,0.5])
Mixer_ACC = Mixer(n_inlets=2)

#%% CONDENSER PARAMETERS

# TubeAndFinsGeom

# Condenser 1
ACC_geom = CFTubeAndFinsGeom()
ACC_geom.set_parameters("DECAGONE_ACC")

Condenser_1.set_parameters(
    A_finned = ACC_geom.A_finned, A_flow = ACC_geom.A_flow, A_in_tot = ACC_geom.A_in_tot, A_out_tot = ACC_geom.A_out_tot, A_unfinned = ACC_geom.A_unfinned, # 5
    B_V_tot = ACC_geom.B_V_tot, Fin_OD = ACC_geom.Fin_OD, Fin_per_m = ACC_geom.Fin_per_m, Fin_t = ACC_geom.Fin_t, Fin_type = ACC_geom.Fin_type, # 10
    Finned_tube_flag = ACC_geom.Tube_t, L = ACC_geom.L, T_V_tot = ACC_geom.T_V_tot, Tube_L = ACC_geom.Tube_L, Tube_OD = ACC_geom.Tube_OD, # 15
    Tube_cond = ACC_geom.Tube_cond, Tube_t = ACC_geom.Tube_t, fouling = ACC_geom.fouling, h = ACC_geom.h, k_fin = ACC_geom.k_fin, # 20
    Tube_pass = ACC_geom.n_passes, n_rows = ACC_geom.n_rows, n_tubes = ACC_geom.n_tubes, pitch = ACC_geom.pitch, pitch_ratio = ACC_geom.pitch_ratio, # 25
    tube_arrang = ACC_geom.tube_arrang, w = ACC_geom.w, # 27

    Fin_Side = 'C', H_DP_ON = True, C_DP_ON = True, n_disc = 30)          

Condenser_2.set_parameters(
    A_finned = ACC_geom.A_finned, A_flow = ACC_geom.A_flow, A_in_tot = ACC_geom.A_in_tot, A_out_tot = ACC_geom.A_out_tot, A_unfinned = ACC_geom.A_unfinned, # 5
    B_V_tot = ACC_geom.B_V_tot, Fin_OD = ACC_geom.Fin_OD, Fin_per_m = ACC_geom.Fin_per_m, Fin_t = ACC_geom.Fin_t, Fin_type = ACC_geom.Fin_type, # 10
    Finned_tube_flag = ACC_geom.Tube_t, L = ACC_geom.L, T_V_tot = ACC_geom.T_V_tot, Tube_L = ACC_geom.Tube_L, Tube_OD = ACC_geom.Tube_OD, # 15
    Tube_cond = ACC_geom.Tube_cond, Tube_t = ACC_geom.Tube_t, fouling = ACC_geom.fouling, h = ACC_geom.h, k_fin = ACC_geom.k_fin, # 20
    Tube_pass = ACC_geom.n_passes, n_rows = ACC_geom.n_rows, n_tubes = ACC_geom.n_tubes, pitch = ACC_geom.pitch, pitch_ratio = ACC_geom.pitch_ratio, # 25
    tube_arrang = ACC_geom.tube_arrang, w = ACC_geom.w, # 27

    Fin_Side = 'C', H_DP_ON = True, C_DP_ON = True, n_disc = 30)   

#%% EVAPORATOR PARAMETERS

Evaporator = HeatExchangerMB('Shell&Tube')

# Evaporator
EVAP_geom = ShellAndTubeGeom()
EVAP_geom.set_parameters("DECAGONE_EVAP_Equ")

Evaporator.set_parameters(
    A_eff = EVAP_geom.A_eff, Baffle_cut = EVAP_geom.Baffle_cut, D_OTL = EVAP_geom.D_OTL, N_strips = EVAP_geom.N_strips, S_V_tot = EVAP_geom.S_V_tot, # 5
    Shell_ID = EVAP_geom.Shell_ID, T_V_tot = EVAP_geom.T_V_tot, Tube_L = EVAP_geom.Tube_L, Tube_OD = EVAP_geom.Tube_OD, Tube_pass = EVAP_geom.n_tube_passes, # 10
    Tube_t = EVAP_geom.Tube_t, Tubesheet_t = EVAP_geom.Tubesheet_t, central_spacing = EVAP_geom.central_spacing, clear_BS = EVAP_geom.clear_BS, clear_TB = EVAP_geom.clear_TB, # 15
    cross_passes = EVAP_geom.cross_passes, foul_s = EVAP_geom.foul_s, foul_t = EVAP_geom.foul_t, inlet_spacing = EVAP_geom.inlet_spacing, n_series = EVAP_geom.n_series, # 20
    n_tubes = EVAP_geom.n_tubes, outlet_spacing = EVAP_geom.outlet_spacing, pitch_ratio = EVAP_geom.pitch_ratio, tube_cond = EVAP_geom.tube_cond, tube_layout = EVAP_geom.tube_layout, # 25

    Shell_Side = 'H', # 26
    Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True, n_disc = 30) # 30


Corr_H_ev = {"1P" : "Shell_Bell_Delaware_HTC", "2P" : "Shell_Bell_Delaware_HTC"}
Corr_C_ev = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

Corr_H_DP = "Shell_Kern_DP"
Corr_C_DP = "Choi_DP"

# Set the heat transfer coefficients correlations of the evaporator           
Evaporator.set_htc(htc_type = 'Correlation', Corr_H = Corr_H_ev, Corr_C = Corr_C_ev)

# Set the pressure drop correlations of the condenser
# Evaporator.set_DP()
Evaporator.set_DP(DP_type="Correlation", Corr_H=Corr_H_DP, Corr_C=Corr_C_DP)


#%% RECUPERATOR PARAMETERS

Recuperator = HeatExchangerMB('Tube&Fins')

# Recuperator
rec_geom = TubeAndFinsGeom()
rec_geom.set_parameters("DECAGONE_RECUP")

Recuperator.set_parameters(
    A_finned = rec_geom.A_finned, A_flow = rec_geom.A_flow, A_in_tot = rec_geom.A_in_tot, A_out_tot = rec_geom.A_out_tot, A_unfinned = rec_geom.A_unfinned, # 5
    B_V_tot = rec_geom.B_V_tot, Fin_OD = rec_geom.Fin_OD, Fin_per_m = rec_geom.Fin_per_m, Fin_t = rec_geom.Fin_t, Fin_type = rec_geom.Fin_type, # 10
    Finned_tube_flag = rec_geom.Tube_t, L = rec_geom.L, T_V_tot = rec_geom.T_V_tot, Tube_L = rec_geom.Tube_L, Tube_OD = rec_geom.Tube_OD, # 15
    Tube_cond = rec_geom.Tube_cond, Tube_t = rec_geom.Tube_t, fouling = rec_geom.fouling, h = rec_geom.h, k_fin = rec_geom.k_fin, # 20
    Tube_pass = rec_geom.n_passes, n_rows = rec_geom.n_rows, n_tubes = rec_geom.n_tubes, pitch = rec_geom.pitch, pitch_ratio = rec_geom.pitch_ratio, # 25
    tube_arrang = rec_geom.tube_arrang, w = rec_geom.w, n_series = 1, # 28

    Fin_Side = 'H', # 28

    Flow_Type = 'CrossFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 30) # 32

Corr_H_rec = {"1P" : "Tube_And_Fins", "2P" : "ext_tube_film_condens"}
Corr_C_rec = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

# Set the pressure drop correlations of the condenser
Recuperator.set_DP()

# Set the heat transfer coefficients correlations of the condenser           
Recuperator.set_htc(htc_type = 'Correlation', Corr_H = Corr_H_rec, Corr_C = Corr_C_rec)

#%% EXPANDER PARAMETERS AND INPUTS

"Geometry"

Turb_geom = Geometry_polyn_turb()
Turb_geom.set_parameters("DECAGONE_TURB") 

"Set params"
Turbine.set_parameters(
                    D_inlet = Turb_geom.D_inlet, N_turb_rated = Turb_geom.N_turb_rated, turb_voltage = Turb_geom.turb_voltage, turb_phases = Turb_geom.turb_phases,
                    eta_max_motor = Turb_geom.eta_max_motor, W_dot_el_rated = Turb_geom.W_dot_el_rated, eta_m = Turb_geom.eta_m, eta_is_coefs = Turb_geom.eta_is_coefs,
                    eta_is_coefs_red = Turb_geom.eta_is_coefs_red, A_th = Turb_geom.A_th 
                    )

Turbine.set_inputs(N_rot = Turb_geom.N_turb_rated)

#%% PUMP INPUTS

"Geometry"

Pump_geom = Geometry_extrapol_pump()
Pump_geom.set_parameters("DECAGONE_PUMP")

"Parameters"
Pump.set_parameters(
                    Omega_rated = Pump_geom.Omega_rated, min_flowrate = Pump_geom.min_flowrate, rated_flowrate = Pump_geom.rated_flowrate, max_flowrate = Pump_geom.max_flowrate,
                    PI_rated = Pump_geom.PI_rated, D_p = Pump_geom.D_p, V_dot_curve = Pump_geom.V_dot_curve, Delta_H_curve = Pump_geom.Delta_H_curve, eta_is_curve = Pump_geom.eta_is_curve,
                    NPSH_r_curve = Pump_geom.NPSH_r_curve, eta_m = Pump_geom.eta_m, eta_max_motor = Pump_geom.eta_max_motor, W_dot_el_rated = Pump_geom.W_dot_el_rated
)

Pump.set_inputs(Omega_pp = 2950)

#%% ADD AND LINK COMPONENTS
ORC.add_component(Turbine, "Turbine")

ORC.add_component(Mixer_ACC, "Mixer")
ORC.add_component(Spliter_ACC, "Spliter")

ORC.add_component(Condenser_1, "ACC_1")
ORC.add_component(Condenser_2, "ACC_2")
ORC.add_component(Recuperator, "Recuperator")
ORC.add_component(Pump, "Pump") # :!\ Est-ce que les noms des composants sont importants?
ORC.add_component(Evaporator, "Evaporator")

# Link components+
ORC.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
ORC.link_components("Recuperator", "m-ex_C", "Evaporator", "m-su_C")
ORC.link_components("Evaporator", "m-ex_C", "Turbine", "m-su")

ORC.link_components("Turbine", "m-ex", "Recuperator", "m-su_H")
ORC.link_components("Recuperator", "m-ex_H", "Spliter", "m-su")

ORC.link_components("Spliter", "m-ex_1", "ACC_1", "m-su_H")
ORC.link_components("Spliter", "m-ex_2", "ACC_2", "m-su_H")

ORC.link_components("ACC_1", "m-ex_H", "Mixer", "m-su_1")
ORC.link_components("ACC_2", "m-ex_H", "Mixer", "m-su_2")

ORC.link_components("Mixer", "m-ex", "Pump", "m-su")

#%% SOURCES AND SINKS

Hot_oil_source = MassConnector()
ORC.add_source("EV_Oil", Hot_oil_source, ORC.components["Evaporator"], "m-su_H")
ORC.set_source_properties(T=310 + 273.15, fluid='INCOMP::T66', m_dot= 19.42, target='EV_Oil', P = 3.25*1e5)  # 19.42

ACC_Air_source_1 = MassConnector()
ORC.add_source("CD_Water", ACC_Air_source_1, ORC.components["ACC_1"], "m-su_C")
ORC.set_source_properties(T=12 + 273.15, fluid='Air', m_dot=158.5, target='CD_Water', P = 1.05e5)

ACC_Air_source_2 = MassConnector()
ORC.add_source("CD_Water", ACC_Air_source_2, ORC.components["ACC_2"], "m-su_C")
ORC.set_source_properties(T=12 + 273.15, fluid='Air', m_dot=158.5, target='CD_Water', P = 1.05e5)

#%% CYCLE GUESSES

T_high = 230+273.15

P_low = 0.9*1e5
P_high = 30*1e5

ORC.set_cycle_guess(target='Turbine:su', p = P_high, T = T_high)
ORC.set_cycle_guess(target='Turbine:ex', p = P_low)

ORC.set_cycle_guess(target='Pump:su', m_dot = 14, SC = 5, p = P_low)

#%% CYCLE FIXED VARIABLES AND ITERATION VARIABLE

ORC.set_fixed_properties(target='Pump:su', SC = 3)    
ORC.set_iteration_variable(target=['Turbine:ex'], variable='p', objective = 'Pump:su-SC', tol = 1e-2, rel = 1, damping_factor = 0.2)

#%% CYCLE RESIDUAL VARIABLES
ORC.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
ORC.set_residual_variable(target='Mixer:ex', variable='h', tolerance= 1e-3)
ORC.set_residual_variable(target='Evaporator:ex_C', variable='m_dot', tolerance= 1e-2)
ORC.set_residual_variable(target='Mixer:ex', variable='m_dot', tolerance= 1e-2)

ORC.solve()

