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
from labothappy.component.expander.turbine_mean_line_Aungier_new import AxialTurbineMeanLine
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

Turb_params = {'type': 'Axial Turbine',
    'mdot_rated': 318.437021666738,
    'Wdot_rated': 15244151.612636594,
    'N_rot_rated': 2864.775003723627,
    'total_to_static_efficiency': 0.8851473369105862,
    'DP_rated': 2.93,
    'n_stages': 7,
    'p0_su': 15309670.5,
    'T0_su': 406.4,
    'p_ex': 5220928,
    'r_m': 0.20956432565203412,
    'delta_tip': 0.0004,
    'N_lw': 0,
    'D_lw': 0,
    'e_blade': 2e-06,
    'stator': {'h_blade_S': [0.024581131784973807,
      0.02731115991390685,
      0.030535241526643685,
      0.0343428949996286,
      0.03884110173212556,
      0.04415656370209263,
      0.05044235012970171,
      0.05396446226956401],
     'chord_S': [0.017387241421509592,
      0.018189646227485024,
      0.019157515021682635,
      0.020309513815586707,
      0.021666371275761282,
      0.023250883041521515,
      0.025089330758218745,
      0.026110623182803244],
     'xhi_S1': [-0.5493774236565352,
      -0.5493774236565352,
      -0.5493774236565352,
      -0.5493774236565352,
      -0.5493774236565352,
      -0.5493774236565352,
      -0.5493774236565352,
      -0.5493774236565352],
     'xhi_S2': [1.1268112059950612,
      1.1268112059950612,
      1.1268112059950612,
      1.1268112059950612,
      1.1268112059950612,
      1.1268112059950612,
      1.1268112059950612,
      1.1268112059950612],
     'pitch_S': [0.014118294382457852,
      0.014769840363216227,
      0.015555741716329778,
      0.01649115508675976,
      0.01759291197811346,
      0.01887952225855718,
      0.020372326403956538,
      0.021201607297417914],
     'o_S': [0.00407007938170583,
      0.004257909709554279,
      0.00448447255112079,
      0.004754137325719167,
      0.005071756287732223,
      0.0054426655373106475,
      0.005873017193716411,
      0.006112085665777279],
     't_TE_S': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
     't_blade_S': [0.0017387241421509593,
      0.0018189646227485025,
      0.0019157515021682636,
      0.002030951381558671,
      0.0021666371275761284,
      0.0023250883041521517,
      0.0025089330758218745,
      0.0026110623182803248],
     'n_blade_S': [93, 89, 85, 80, 75, 70, 65, 62],
     'R_c_S': [0.02593658923532629,
      0.02713353838605762,
      0.02857731056016104,
      0.03029574989120514,
      0.03231977738022937,
      0.034683397336424615,
      0.03742581415244761,
      0.03894927848261039]},
    'rotor': {'h_blade_R': [0.02599623192810525,
      0.028975662485097392,
      0.03249353115560146,
      0.03664835013019067,
      0.041556335411602884,
      0.04735747516384924,
      0.054218932168004906,
      None],
     'chord_R': [0.017804442744401783,
      0.018691926981905842,
      0.019754111574023726,
      0.0210108023910845,
      0.02248355410347207,
      0.02419698360942641,
      0.02617885482072475,
      None],
     'xhi_R1': [0.5033968740474793,
      0.5033968740474793,
      0.5033968740474793,
      0.5033968740474793,
      0.5033968740474793,
      0.5033968740474793,
      0.5033968740474793,
      None],
     'xhi_R2': [-1.1544956444218015,
      -1.1544956444218015,
      -1.1544956444218015,
      -1.1544956444218015,
      -1.1544956444218015,
      -1.1544956444218015,
      -1.1544956444218015,
      None],
     'pitch_R': [0.015148494029818542,
      0.015903589253319174,
      0.01680732419089942,
      0.017876550204477332,
      0.01912960658162921,
      0.02058743803494024,
      0.022273666839919944,
      None],
     'o_R': [0.004776675516045738,
      0.005014774752793548,
      0.005299743578132708,
      0.005636895621769747,
      0.006032013691274512,
      0.006491702145839364,
      0.007023409643056108,
      None],
     't_TE_R': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
     't_blade_R': [0.0019275069515591894,
      0.0020395300073626607,
      0.002174650407702253,
      0.0022904090712424882,
      0.0024358345376652743,
      0.002590844225248941,
      0.002758858045024202,
      None],
     'n_blade_R': [87, 83, 78, 74, 69, 64, 59, None],
     'R_c_R': [0.02655892943743007,
      0.02788279177780649,
      0.029467255051191132,
      0.03134186372129275,
      0.03353877094100647,
      0.036094697795777626,
      0.03905105977860335,
      None]}}
    
Turbine.set_parameters(
    r_m = Turb_params['r_m'],
    nStages = Turb_params['n_stages'],
    mdot_rated = Turb_params['mdot_rated'],
    DP_rated = Turb_params['DP_rated'],
    damping = 0.3,
    delta_tip = Turb_params['delta_tip'],
    N_lw = Turb_params['N_lw'],
    D_lw = Turb_params['D_lw'],
    e_blade = Turb_params['e_blade'],
    solve_type = 'map'
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

# MAP_SAVE_PATH = r"C:\Users\Basile\Desktop\Travail\Thèse\Travail\WP1\Turbomachines\Save Maps\turb_map.parquet"   # <-- set your path here

# df_map = Turbine.load_map_df(MAP_SAVE_PATH)

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
GH_C_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Choi_DP"}

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

REC_Corr_H = {"1P" : "Gnielinski", "SC" : "Gnielinski", "2P" : "Thome_Condensation"}
REC_Corr_C = {"1P" : "Gnielinski", "SC" : "Gnielinski", "2P" : "Flow_boiling"}   

REC_Corr_H_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach", "2P" : "Choi_DP"}
REC_Corr_C_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach", "2P" : "Choi_DP"}  

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
m_dot_w_gh = 210 #192.24 # kg/s

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
N_turb = 2865 # RPM

Turbine.set_inputs(
      N_rot = N_turb, 
    )

# -------- 8) Set cycle guesses --------

mdot_CO2_guess = 320 # kg/s
P_LP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
P_HP_guess = 153*1e5

T_sat_LP_guess = PropsSI("T", "P", P_LP_guess, "Q", 0.5, fluid)
h_SC_guess = PropsSI("H", "P", P_LP_guess, "T", T_sat_LP_guess - SC_cd, fluid)

h_su_exp_guess = PropsSI("H", "P", P_HP_guess, "T", T_su_w_gh - 5, fluid) # J/kg
# h_ex_exp_guess = 470000 # J/kg

CO2_TC.set_cycle_guess(target="Pump:su", m_dot = mdot_CO2_guess, h=h_SC_guess, p=P_LP_guess)
CO2_TC.set_cycle_guess(target="Pump:ex", p=P_HP_guess)

CO2_TC.set_cycle_guess(target="Turbine:su", m_dot = mdot_CO2_guess, p=P_HP_guess, h = h_su_exp_guess)
# CO2_TC.set_cycle_guess(target="Turbine:ex", p=P_LP_guess, h = h_ex_exp_guess)

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
