# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:43:31 2026

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

from labothappy.machine.circuit import Circuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.expander.turbine_mean_line_Aungier import AxialTurbineMeanLine
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity

# -------- 1) Instanciate Circuit --------
fluid = 'CO2'
CO2_TC = Circuit(fluid)

# -------- 2) Create components --------
Turbine = AxialTurbineMeanLine(fluid)
GasHeater = HexMBChargeSensitive('Shell&Tube')
Recuperator = HexMBChargeSensitive('PCHE')
Pump = PumpCurveSimilarity()
Condenser = HexMBChargeSensitive('Shell&Tube')

# -------- 3) Set component parameters --------
# 3.1) Turbine

Turb_params = {
    "type": "Axial Turbine",
    "mdot_rated": 311.29856776173153,
    "Wdot_rated": 20735708.67089974,
    "N_rot_rated": 1798.1320562528444,
    "total_to_static_efficiency": 0.927789124530079,
    "DP_rated": 2.84,
    "n_stages": 7,
    "p0_su": 14654008.433412187,
    "T0_su": 406.1822040393341,
    "p_ex": 5159550.900413497,
    "r_m": 0.19067839204343143,
    "delta_tip": 0.0004,
    "N_lw": 0,
    "D_lw": 0,
    "e_blade": 2e-06,
    "stator": {
        "h_blade_S": [0.024298971423751134, 0.028007674498303256, 0.03261359682274655, 0.038339585990648024, 0.045469718325074224, 0.054372703804122034, 0.06553303324259455, 0.07209225754646112],
        "chord_S": [0.007018333827634321, 0.007490359631955781, 0.00808167675881387, 0.008810188506889813, 0.009696931119931251, 0.010767749805129607, 0.012054864777445069, 0.012790316070895836],
        "xhi_S1": [-0.6274864722618412, -0.6274864722618412, -0.6274864722618412, -0.6274864722618412, -0.6274864722618412, -0.6274864722618412, -0.6274864722618412, -0.6274864722618412],
        "xhi_S2": [1.1314694189847898, 1.1314694189847898, 1.1314694189847898, 1.1314694189847898, 1.1314694189847898, 1.1314694189847898, 1.1314694189847898, 1.1314694189847898],
        "pitch_S": [0.005504391139311712, 0.005874595053608749, 0.006338357655037761, 0.00690972027606392, 0.007605181378648812, 0.00844501103442938, 0.0094544791536278, 0.010031284372998545],
        "o_S": [0.0013860567332180695, 0.0014792775118815668, 0.001596057235570641, 0.001739931642010707, 0.0019150551969202856, 0.0021265320923097414, 0.0023807255259106295, 0.0025259704290851734],
        "t_TE_S": [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
        "t_blade_S": [0.0024564168396720124, 0.002621625871184523, 0.0028285868655848543, 0.0030835659774114343, 0.003393925891975938, 0.003768712431795362, 0.004219202672105774, 0.004476610624813542],
        "n_blade_S": [218, 204, 189, 173, 158, 142, 127, 119],
        "R_c_S": [0.010836790096933203, 0.01156563040111291, 0.012478664711229602, 0.013603536951698384, 0.014972728529576805, 0.016626146222203153, 0.01861354025733473,0.019749127633110847]},
    "rotor": {
        "h_blade_R": [0.026263942408487173, 0.030438887206607458, 0.03562526346906793, 0.042076712413706994, 0.05011970128497465, 0.06018139437431045, 0.07282683153995556, None],
        "chord_R": [0.007269199906964899, 0.007804550313630929, 0.008468446881519582, 0.009280327566310602, 0.0102636000429607, 0.011447445439118742, 0.012868315398534612, None],
        "xhi_R1": [0.5366146308074933, 0.5366146308074933, 0.5366146308074933, 0.5366146308074933, 0.5366146308074933, 0.5366146308074933, 0.5366146308074933, None],
        "xhi_R2": [-1.1673268351691564, -1.1673268351691564, -1.1673268351691564, -1.1673268351691564, -1.1673268351691564, -1.1673268351691564, -1.1673268351691564, None],
        "pitch_R": [0.00628592795237799, 0.006748863919010937, 0.007322958186191646, 0.008025019424820382, 0.008875289059014475, 0.009899001016622598, 0.011127676291604209, None],
        "o_R": [0.0019234337389752406, 0.002065087710187251,0.002240755056554132, 0.0024555790703569644, 0.002715753433478063, 0.00302899948611711, 0.003404962350478334, None],
        "t_TE_R": [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
        "t_blade_R": [0.0025442199674377145, 0.002731592609770825, 0.0029639564085318535,  0.0032481146482087105, 0.0035922600150362447, 0.00400660590369156, 0.004503910389487113, None],
        "n_blade_R": [191, 178, 164, 149, 135, 121, 108, None],
        "R_c_R": [0.011224144576060682, 0.012050762421239916, 0.013075864379763989, 0.01432946399199564, 0.015847704317868035, 0.017675643025324723, 0.019869563959179458, None]},
    }
   
Turbine.set_parameters(
    r_m = Turb_params['r_m'],
    nStages = Turb_params['n_stages'],
    mdot_rated = Turb_params['mdot_rated'],
    DP_rated = Turb_params['DP_rated'],
    damping = 0.3,
    delta_tip = Turb_params['delta_tip'],
    N_lw = Turb_params['N_lw'],
    D_lw = Turb_params['D_lw'],
    e_blade = Turb_params['e_blade']
    )
       
Turbine.set_stage_parameters(
    h_blade_S = Turb_params['stator']["h_blade_S"],
    chord_S = Turb_params['stator']["chord_S"],
    xhi_S1 = Turb_params['stator']["xhi_S1"],
    xhi_S2 = Turb_params['stator']["xhi_S2"],
    pitch_S = Turb_params['stator']["pitch_S"], 
    o_S = Turb_params['stator']["o_S"], 
    t_TE_S = Turb_params['stator']["t_TE_S"], 
    t_blade_S = Turb_params['stator']["t_blade_S"], 
    n_blade_S = Turb_params['stator']["n_blade_S"], 
    R_c_S = Turb_params['stator']["R_c_S"], 
    
    h_blade_R = Turb_params['rotor']["h_blade_R"],
    chord_R = Turb_params['rotor']["chord_R"],
    xhi_R1 = Turb_params['rotor']["xhi_R1"],
    xhi_R2 = Turb_params['rotor']["xhi_R2"],
    pitch_R = Turb_params['rotor']["pitch_R"],
    o_R = Turb_params['rotor']["o_R"],
    t_TE_R = Turb_params['rotor']["t_TE_R"],
    t_blade_R = Turb_params['rotor']["t_blade_R"],
    n_blade_R = Turb_params['rotor']["n_blade_R"],
    R_c_R = Turb_params['rotor']["R_c_R"],
    )

# 3.2) GasHeater
"Geometry Loading"
params_GH = {'n_series': 1,
          'n_parallel': 3,
          'foul_t': 0.000176,
          'foul_s': 0.000176,
          'tube_cond': 20,
          'Shell_Side': 'H',
          'Flow_Type': 'Shell&Tube',
          'n_disc': 20,
          'Shell_ID': 1.8288,
          'Tube_L': 5.832,
          'Tube_OD': 0.009525,
          'Tube_t': 0.0005588,
          'central_spacing': 0.37,
          'Tube_pass': 1,
          'n_tubes': 16164,
          'pitch_ratio': 1.33,
          'tube_layout': 0,
          'Baffle_cut': 34.973}

GasHeater.set_parameters(
    Baffle_cut = params_GH['Baffle_cut'], Shell_ID = params_GH['Shell_ID'], Tube_L = params_GH['Tube_L'], 
    Tube_OD = params_GH['Tube_OD'], Tube_pass = params_GH['Tube_pass'], Tube_t = params_GH['Tube_t'],
    central_spacing = params_GH['central_spacing'], foul_s = params_GH['foul_s'],
    foul_t = params_GH['foul_t'], n_series = params_GH['n_series'], n_parallel = params_GH['n_parallel'], n_tubes = params_GH['n_tubes'], 
    pitch_ratio = params_GH['pitch_ratio'], tube_cond = params_GH['tube_cond'], tube_layout = params_GH['tube_layout'],

    Shell_Side = params_GH['Shell_Side'],

    Flow_Type = params_GH['Flow_Type'], n_disc = params_GH['n_disc'])

GH_H_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
GH_C_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Flow_boiling"}

GH_H_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
GH_C_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Gnielinski_DP"}

GasHeater.set_htc(htc_type="Correlation_Disc", Corr_H=GH_H_Corr, Corr_C=GH_C_Corr) # 'User-Defined' or 'Correlation' # 28
GasHeater.set_DP(DP_type="Correlation_Disc", Corr_C=GH_C_DP, Corr_H=GH_H_DP)

# 3.3) Condenser
params_CD = {"n_series": 1,
            "n_parallel": 10.0,
            "foul_t": 0.000176,
            "foul_s": 0.000176,
            "tube_cond": 20,
            "Shell_Side": "C",
            'Flow_Type': 'Shell&Tube',
            'n_disc': 20,
            "Shell_ID": 1.6764,
            "Tube_L": 6.8,
            "Tube_OD": 0.0254,
            "Tube_t": 0.0012446,
            "central_spacing": 1.7,
            "Tube_pass": 1.0,
            "cross_passes": 3,
            "n_tubes": 2106,
            "pitch_ratio": 1.25,
            "tube_layout": 0.0,
            "Baffle_cut": 18.0,
          }

Condenser.set_parameters(
    Baffle_cut = params_CD['Baffle_cut'], Shell_ID = params_CD['Shell_ID'], Tube_L = params_CD['Tube_L'], 
    Tube_OD = params_CD['Tube_OD'], Tube_pass = params_CD['Tube_pass'], Tube_t = params_CD['Tube_t'],
    central_spacing = params_CD['central_spacing'], foul_s = params_CD['foul_s'],
    foul_t = params_CD['foul_t'], n_series = params_CD['n_series'], n_parallel = params_CD['n_parallel'], n_tubes = params_CD['n_tubes'], 
    pitch_ratio = params_CD['pitch_ratio'], tube_cond = params_CD['tube_cond'], tube_layout = params_CD['tube_layout'],

    Shell_Side = params_CD['Shell_Side'],

    Flow_Type = params_CD['Flow_Type'], n_disc = params_CD['n_disc'])

CD_H_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Thome_Condensation"}
CD_C_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

CD_H_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Choi_DP"}
CD_C_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}

Condenser.set_htc(htc_type="Correlation_Disc", Corr_H=CD_H_Corr, Corr_C=CD_C_Corr) # 'User-Defined' or 'Correlation' # 28
Condenser.set_DP(DP_type="Correlation_Disc", Corr_C=CD_C_DP, Corr_H=CD_H_DP)

# 3.4) Pump
D_H_rated = 1211 # m
V_dot_rated = 1298 # m3/h
eta_rated = 0.867 # -

# Example characteristic curves and parameters
V_dot_curve =  np.array([76.2, 91.0, 113.7, 136.3, 158.9, 181.8, 204.7, 227.6, 250.2, 272.9, 295.5])* V_dot_rated / 250.2
eta_is_curve = np.array([40.4, 47.2,  56.0,  62.6,  68.4,  72.8,  76.1,  78.5,  79.2,  78.7,  76.8])* eta_rated / 79.2
Delta_H_curve = np.array([1406.5, 1394.7, 1370.9, 1342.4, 1309.0, 1270.8, 1232.6, 1177.6, 1105.7, 1014.6, 918.7])* D_H_rated / 1105.7
NPSH_r_curve = np.array([1.1, 1.1, 1.25, 1.4, 1.6, 1.8, 1.9, 2, 3, 3.85, 4.7])  # m, increases again near max flow
N_rated = 2900 # RPM
# Reference point: water at 20°C and 1 atm for a rated speed of 2900 RPM

Pump = PumpCurveSimilarity()

# Set Parameters
Pump.set_parameters(
    V_dot_curve = V_dot_curve,
    Delta_H_curve = Delta_H_curve,
    eta_is_curve = eta_is_curve,
    NPSH_r_curve = NPSH_r_curve,
    N_rot_rated = N_rated,
    mode = "P_N",  # Mode can be "M_N", "P_M", or "P_N"
)

# 3.5) Recuperator

params_REC = {'alpha': 32.62, # Channel zigzag angle
          'D_c': 2.42*1e-3, # Channel diameter
          'C_V_tot' : 1, 
          'H_V_tot' : 1, 
          'k_cond': 60, # plate conductivity
          'L_c': 0.7432303013776589, # channel length
          'N_c': 736, # n channels per plate
          'N_p': 563, # n plates
          'R_p': 1, # n_hot_channel_row / n_cold_channel_row
          't_2': 0.0012282802564224898, # Horizontal pitch
          't_3': 0.0009428803890487963, # Plate_thickness
          'type_channel' : 'Zigzag',
          "AS_Type" : "HEOS"} 

Recuperator.set_parameters(
    alpha = params_REC['alpha'], C_V_tot = params_REC['C_V_tot'], H_V_tot = params_REC['H_V_tot'], D_c = params_REC['D_c'], k_cond = params_REC['k_cond'], L_c = params_REC['L_c'], 
    N_c = params_REC['N_c'], N_p = params_REC['N_p'], R_p = params_REC['R_p'], t_2 = params_REC['t_2'], t_3 = params_REC['t_3'],
    
    Flow_Type = 'CounterFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 20)

REC_Corr_H = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
REC_Corr_C = {"1P" : "Gnielinski", "SC" : "Gnielinski"}   

REC_Corr_H_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach"}
REC_Corr_C_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach"}  

Recuperator.set_htc(htc_type = 'Correlation_Disc', Corr_H = REC_Corr_H, Corr_C = REC_Corr_C) # 'User-Defined' or 'Correlation' # 28
Recuperator.set_DP(DP_type="Correlation_Disc", Corr_C=REC_Corr_C_DP, Corr_H=REC_Corr_H_DP)

#%%

# -------- 4) Add components to circuit --------
CO2_TC.add_component(Turbine, "Turbine")
CO2_TC.add_component(Condenser, "Condenser")
CO2_TC.add_component(Pump, "Pump")
CO2_TC.add_component(GasHeater, "GasHeater")
CO2_TC.add_component(Recuperator, "Recuperator")

# -------- 5) Link components with mass connectors --------
CO2_TC.link_components("Turbine", "m-ex", "Recuperator", "m-su_H")
CO2_TC.link_components("Recuperator", "m-ex_H", "Condenser", "m-su_H")
CO2_TC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
CO2_TC.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
CO2_TC.link_components("Recuperator", "m-ex_C", "GasHeater", "m-su_C")
CO2_TC.link_components("GasHeater", "m-ex_C", "Turbine", "m-su")

# -------- 6) Add fluid sources --------
CD_source = MassConnector('Water')
T_su_w_cd = 0.1 + 273.15
P_su_w_cd = 5e5
m_dot_w_cd = 4542 # kg/s
GH_source = MassConnector('Water')
T_su_w_gh = 150+273.15
P_su_w_gh = 5e5
m_dot_w_gh = 192.24 # kg/s

CO2_TC.add_source("CD_Water", CD_source, CO2_TC.components["Condenser"], "m-su_C")
CO2_TC.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
CO2_TC.add_source("GH_Water", GH_source, CO2_TC.components["GasHeater"], "m-su_H")
CO2_TC.set_source_properties(T=T_su_w_gh, fluid='Water', P=P_su_w_gh, m_dot = m_dot_w_gh, target="GH_Water")

# -------- 7) Set cycle inputs --------

N_pp = 2900 # 

# Set Inputs
Pump.set_inputs(
    N_rot=N_pp,  # Rotational speed in RPM
    fluid="CO2",  # Actual fluid type
)

SC_cd = 5  # K
N_turb = 1800 # RPM

Turbine.set_inputs(
      N_rot = N_turb, 
    )

# -------- 8) Set cycle guesses --------

mdot_CO2_guess = 250 # kg/s
P_LP_guess = PropsSI("P", "T", T_su_w_cd+5, "Q", 0, fluid)
P_HP_guess = 160*1e5

T_sat_LP_guess = PropsSI("T", "P", P_LP_guess, "Q", 0.5, fluid)
h_SC_guess = PropsSI("H", "P", P_LP_guess, "T", T_sat_LP_guess - SC_cd, fluid)

h_ex_exp_guess = 470000 # J/kg

CO2_TC.set_cycle_guess(target="Pump:su", m_dot = mdot_CO2_guess, h=h_SC_guess, p=P_LP_guess)
CO2_TC.set_cycle_guess(target="Pump:ex", p=P_HP_guess)

CO2_TC.set_cycle_guess(target="Turbine:ex", m_dot = mdot_CO2_guess, p=P_LP_guess, h = h_ex_exp_guess)

# orc.set_iteration_variable(
#     it_var  = 'Expander:ex-p',
#     objective = 'Pump:su-SC',
#     target_value = SC_cd,
#     obj_type = "Target_val",
#     damping_factor = 0.3,
# )

# orc.set_iteration_variable(
#     it_var  = 'Pump:ex-p',
#     objective = 'Expander:W-N_rot',
#     target_value = N_exp,
#     obj_type = "Target_val",
#     damping_factor = 0.1,
# )

# # -------- 8) Solve — swap method here for comparison --------
METHOD = 'wegstein'   # <-- change to compare: 'successive_substitution',
# #                       #     'wegstein', 'fsolve', 'lm', 'broyden1', 'anderson'

# # METHOD = 'anderson'   # <-- change to compare: 'successive_substitution',
#                       #     'wegstein', 'fsolve', 'lm', 'broyden1', 'anderson'

CO2_TC.solve(method=METHOD, max_iter=100)

# print(f"\n[{METHOD}] Converged: {orc.converged}")
# print(f"  P_HP = {Expander.su.p:.2f} Pa")
# print(f"  P_LP = {Expander.ex.p:.2f} Pa")
# print(f"  T_su_expander = {Expander.su.T - 273.15:.2f} °C")
# print(f"  m_dot = {Pump.su.m_dot:.4f} kg/s")
# print(f"  W_exp = {Expander.W.W_dot:.2f} W")
# print(f"  W_pump = {Pump.W.W_dot:.2f} W")
# print(f"  Q_ev = {Evaporator.Q.Q_dot:.2f} W")
# print(f"  Q_cd = {Condenser.Q.Q_dot:.2f} W")

# """
# 1) Fixer des variables constantes 
    
# 2) Réfléchir à ce que ça implique niveau variables d'itération
    
# 3) Quid DeltaP ? 
    
# """
