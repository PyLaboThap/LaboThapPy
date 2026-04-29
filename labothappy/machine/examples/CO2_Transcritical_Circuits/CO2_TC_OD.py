# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:43:31 2026

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt

from labothappy.machine.circuit_fpi import CircuitFPI
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.expander.turbine_mean_line_Aungier import AxialTurbineMeanLine
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity

import os
os.environ["NUMEXPR_MAX_THREADS"] = "1"    # no threading inside workers
os.environ["OMP_NUM_THREADS"]     = "1"    # same for OpenMP (numpy, scipy)
os.environ["MKL_NUM_THREADS"]     = "1"    # same for MKL
os.environ["OPENBLAS_NUM_THREADS"]= "1"    # same for OpenBLAS

def CO2_OD_TC(mdot, P_high, mute_print=1):
    # -------- 1) Instanciate Circuit --------
    fluid = 'CO2'
    CO2_TC = CircuitFPI(fluid)
    
    if mute_print:
        CO2_TC.mute_print()
    
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
        # solve_type = 'map'
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
        mode = "P_M",  # Mode can be "M_N", "P_M", or "P_N"
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
    m_dot_w_gh = 192.24 # kg/s
    
    CO2_TC.add_source("CD_Water", CD_source, CO2_TC.components["Condenser"], "m-su_C")
    CO2_TC.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
    CO2_TC.add_source("GH_Water", GH_source, CO2_TC.components["GasHeater"], "m-su_H")
    CO2_TC.set_source_properties(T=T_su_w_gh, fluid='Water', P=P_su_w_gh, m_dot = m_dot_w_gh, target="GH_Water")

    # -------- 7) Set cycle inputs --------
    
    # Set Inputs
    Pump.set_inputs(
        fluid="CO2",  # Actual fluid type
    )
    
    N_turb = Turb_params['N_rot_rated'] # RPM
    
    Turbine.set_inputs(
          N_rot = N_turb, 
        )

    # -------- 8) Set cycle guesses --------
    
    mdot_CO2_guess = mdot # kg/s
    P_LP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
    P_HP_guess = P_high
    
    T_sat_LP_guess = PropsSI("T", "P", P_LP_guess, "Q", 0.5, fluid)
    h_SC_guess = PropsSI("H", "P", P_LP_guess, "T", T_sat_LP_guess - 5, fluid)
    
    h_su_exp_guess = PropsSI("H", "P", P_HP_guess, "T", T_su_w_gh - 10, fluid) # J/kg
    
    CO2_TC.set_cycle_guess(target="Pump:su", m_dot = mdot_CO2_guess, h=h_SC_guess, p=P_LP_guess)
    CO2_TC.set_cycle_guess(target="Pump:ex", p=P_HP_guess)
    
    CO2_TC.set_cycle_guess(target="Turbine:su", m_dot = mdot_CO2_guess, p=P_HP_guess, h = h_su_exp_guess)

    return CO2_TC

#%%

case_study = "Simulation"

if case_study == "Simulation":

    # # # -------- 8) Solve — swap method here for comparison --------
    METHOD = 'wegstein'   # <-- change to compare: 'successive_substitution',

    mdot = 320 # kg/s
    P_high = 153*1e5 # Pa
    CO2_TC = CO2_OD_TC(mdot, P_high, mute_print=0)
    CO2_TC.solve(method=METHOD, max_iter=100)

    print(f"SC : {CO2_TC.components['Pump'].model.su.SC}")
    print(f"N_pp : {CO2_TC.components['Pump'].model.W.N_rot}")

    W_pp = CO2_TC.components['Pump'].model.W.W_dot
    W_turb = CO2_TC.components['Turbine'].model.W.W_dot
    Q_gh = CO2_TC.components['GasHeater'].model.Q.Q_dot

    eta = (W_turb - W_pp)/Q_gh

    print(f"W_pp : {W_pp}")
    print(f"W_turb : {W_turb}")
    print(f"Q_gh : {Q_gh}")
    print(f"eta : {eta}")

elif case_study == "Sensitivity":
    
    # -------- TARGET SPECIFICATIONS --------
    N_pump_target = 2900   # RPM
    SC_cd_target  = 5.0    # K
    
    # -------- SCALING FACTORS --------
    P_scale    = 1e5
    mdot_scale = 1.0
    
    last_obj = [np.inf]
        
    import numpy as np
    import matplotlib.pyplot as plt
    
    # -------- SWEEP RANGES --------
    P_sweep    = np.linspace(140, 160, 40) * 1e5   # Pa
    mdot_sweep = np.linspace(290, 330, 40)          # kg/s
    
    from joblib import Parallel, delayed
    from tqdm import tqdm
    
    # -------- SINGLE POINT EVALUATION --------
    def evaluate_point(i, j, m, P):
        try:
            CO2_TC = CO2_OD_TC(m, P)
            CO2_TC.mute_print()
            CO2_TC.solve(method='wegstein', max_iter=100)
        except Exception as e:
            print(f"  [{i+1}/{len(mdot_sweep)}, {j+1}/{len(P_sweep)}]"
                  f" [EXCEPTION] m={m:.1f} kg/s | P={P/1e5:.1f} bar — {e}")
            return i, j, np.nan, np.nan, np.nan, np.nan
    
        if not CO2_TC.converged:
            print(f"  [{i+1}/{len(mdot_sweep)}, {j+1}/{len(P_sweep)}]"
                  f" [NOT CONVERGED] m={m:.1f} kg/s | P={P/1e5:.1f} bar")
            return i, j, np.nan, np.nan, np.nan, np.nan
    
        N_pump_actual = CO2_TC.components["Pump"].model.W.N_rot
        SC_actual     = CO2_TC.components["Condenser"].model.ex_H.SC
        
        if SC_actual is None:
            SC_actual = 0
    
        r0 = (N_pump_actual - N_pump_target) / N_pump_target
        r1 = (SC_actual     - SC_cd_target)  / SC_cd_target
    
        W_dot_turb = CO2_TC.components["Turbine"].model.W.W_dot
        W_dot_pp   = CO2_TC.components["Pump"].model.W.W_dot
        Q_dot_gh   = CO2_TC.components["GasHeater"].model.Q.Q_dot
        eta        = (W_dot_turb - W_dot_pp) / Q_dot_gh
    
        obj = r0**2 + r1**2
    
        print(f"  [{i+1}/{len(mdot_sweep)}, {j+1}/{len(P_sweep)}]"
              f" P={P/1e5:.1f} bar | m={m:.1f} kg/s"
              f" | N_pump={N_pump_actual:.1f} RPM | SC={SC_actual:.2f} K"
              f" | obj={obj:.2e}")
    
        return i, j, obj, N_pump_actual, SC_actual, eta
    
    
    # -------- BUILD TASK LIST --------
    tasks = [
        (i, j, m, P)
        for i, m in enumerate(mdot_sweep)
        for j, P in enumerate(P_sweep)
    ]
    
    # -------- RUN IN PARALLEL --------
    results = Parallel(n_jobs=-2, backend='loky', verbose=0)(
        delayed(evaluate_point)(i, j, m, P)
        for i, j, m, P in tqdm(tasks, desc="2D Sweep", unit="pt", total=len(tasks))
    )
    
    # -------- FILL ARRAYS --------
    obj_2D    = np.full((len(mdot_sweep), len(P_sweep)), np.nan)
    N_pump_2D = np.full((len(mdot_sweep), len(P_sweep)), np.nan)
    SC_2D     = np.full((len(mdot_sweep), len(P_sweep)), np.nan)
    eta_2D    = np.full((len(mdot_sweep), len(P_sweep)), np.nan)
    
    for i, j, obj, N_pump, SC, eta in results:
        obj_2D[i, j]    = obj
        N_pump_2D[i, j] = N_pump
        SC_2D[i, j]     = SC
        eta_2D[i, j]    = eta
    
    # -------- PLOT --------
    P_grid, mdot_grid = np.meshgrid(P_sweep / 1e5, mdot_sweep)
    
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fig.suptitle("2D Sensitivity Analysis", fontsize=14)
    
    datasets = [
        (obj_2D,    "Objective [-]",  "Objective"),
        (N_pump_2D, "N_pump [RPM]",   "N_pump"),
        (SC_2D,     "SC [K]",         "SC"),
        (eta_2D,    "eta [-]",       "Efficiency"),
    ]
    
    for ax, (data, label, title) in zip(axes, datasets):
        cf = ax.contourf(P_grid, mdot_grid, data, levels=20, cmap='viridis')
        fig.colorbar(cf, ax=ax, label=label)
    
        # Overlay contour lines for readability
        cs = ax.contour(P_grid, mdot_grid, data, levels=10, colors='white',
                        linewidths=0.5, alpha=0.5)
        ax.clabel(cs, inline=True, fontsize=7, fmt='%.2g')
    
        # Mark initial guess
        ax.plot(P_high / 1e5, mdot, 'r*', markersize=12, label='Initial guess')
    
        ax.set_xlabel("P_HP [bar]")
        ax.set_ylabel("m_dot [kg/s]")
        ax.set_title(title)
        ax.legend(fontsize=8)
    
    # Add target contours on N_pump and SC plots
    for ax, data, target in [
        (axes[1], N_pump_2D, N_pump_target),
        (axes[2], SC_2D,     SC_cd_target),
    ]:
        ax.contour(P_grid, mdot_grid, data, levels=[target],
                    colors='red', linewidths=2, linestyles='--')
        ax.plot([], [], 'r--', linewidth=2, label='Target')
        ax.legend(fontsize=8)
    
    plt.tight_layout()
    plt.show()


elif case_study == "Optimization_minimize":

    mdot = 320 # kg/s
    P_high = 153*1e5 # Pa
    
    # -------- TARGET SPECIFICATIONS --------
    N_pump_target = 2900   # RPM
    SC_cd_target  = 5.0    # K
    
    # -------- SCALING FACTORS --------
    P_scale    = 1e5
    mdot_scale = 1.0
    
    #%% Minimize

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import minimize, NonlinearConstraint
    
    # -------- CONSTRAINT TOLERANCES --------
    tol_N_pump = 0.01 * N_pump_target   # 1% = 29 RPM
    tol_SC     = 0.02 * SC_cd_target    # 2% = 0.1 K
    
    # -------- SCALING FACTORS --------
    P_scale    = 1e5
    mdot_scale = 1.0
    
    # -------- INITIAL GUESS --------
    x0 = [P_high / P_scale,
          mdot    / mdot_scale]
    
    # -------- BOUNDS --------
    bounds = [(130, 170),    # P_HP in bar
              (250, 400)]    # m_dot in kg/s
    
    # -------- CACHE --------
    _cache = {'x': None, 'result': None}
    
    # -------- CORE EVALUATION --------
    def evaluate_cycle(x):
        x = np.asarray(x)
    
        # Return cached result if x unchanged
        if _cache['x'] is not None and np.allclose(x, _cache['x'], rtol=0, atol=1e-12):
            return _cache['result']
    
        P_high    = x[0] * P_scale
        m_dot_CO2 = x[1] * mdot_scale
    
        if P_high < 1e6 or m_dot_CO2 < 1.0:
            _cache['x']      = x.copy()
            _cache['result'] = None
            return None
    
        try:
            CO2_TC = CO2_OD_TC(m_dot_CO2, P_high)
            CO2_TC.mute_print()
            CO2_TC.solve(method='wegstein', max_iter=100)
        except Exception as e:
            print(f"  [EXCEPTION] {e}")
            _cache['x']      = x.copy()
            _cache['result'] = None
            return None
    
        if not CO2_TC.converged:
            print(f"  [NOT CONVERGED] P_HP={P_high/1e5:.2f} bar | mdot={m_dot_CO2:.2f} kg/s")
            _cache['x']      = x.copy()
            _cache['result'] = None
            return None
    
        N_pump_actual = CO2_TC.components["Pump"].model.W.N_rot
        SC_actual     = CO2_TC.components["Condenser"].model.ex_H.SC
        W_dot_turb    = CO2_TC.components["Turbine"].model.W.W_dot
        W_dot_pp      = CO2_TC.components["Pump"].model.W.W_dot
        Q_dot_gh      = CO2_TC.components["GasHeater"].model.Q.Q_dot
        eta           = (W_dot_turb - W_dot_pp) / Q_dot_gh
    
        print(f"  P_HP={P_high/1e5:.2f} bar | mdot={m_dot_CO2:.2f} kg/s"
              f" | N_pump={N_pump_actual:.1f} RPM | SC={SC_actual:.2f} K"
              f" | eta={eta:.4f}")
    
        result = (eta, N_pump_actual, SC_actual)
        _cache['x']      = x.copy()
        _cache['result'] = result
        return result
    
    
    # -------- OBJECTIVE --------
    def efficiency_objective(x):
        out = evaluate_cycle(x)
        if out is None:
            return 1e6
        eta, _, _ = out
        return -eta
    
    
    # -------- CONSTRAINTS (vector-valued — single call for both) --------
    def all_outputs(x):
        out = evaluate_cycle(x)
        if out is None:
            return [N_pump_target * 10, SC_cd_target * 10]   # large violation
        _, N_pump, SC = out
        return [N_pump, SC]
    
    nc = NonlinearConstraint(
        fun = all_outputs,
        lb  = [N_pump_target - tol_N_pump, SC_cd_target - tol_SC],
        ub  = [N_pump_target + tol_N_pump, SC_cd_target + tol_SC],
    )
    
    # -------- SOLVE --------
    result = minimize(
        efficiency_objective,
        x0,
        method      = 'trust-constr',
        bounds      = bounds,
        constraints = nc,
        options     = {
            'maxiter'             : 100,
            'disp'                : True,
            'gtol'                : 1e-4,
            'xtol'                : 1e-4,
            'finite_diff_rel_step': [1.0 / 170, 1.0 / 400],
        }
    )
    
    # -------- RESULTS --------
    P_HP_sol  = result.x[0] * P_scale
    m_dot_sol = result.x[1] * mdot_scale
    
    # Re-evaluate solution point to get all outputs
    out_sol = evaluate_cycle(result.x)
    eta_sol, N_pump_sol, SC_sol = out_sol if out_sol is not None else (np.nan, np.nan, np.nan)
    
    print(f"\nOptimization converged : {result.success}")
    print(f"  Message              : {result.message}")
    print(f"  Iterations           : {result.nit}")
    print(f"  Function evals       : {result.nfev}")
    print(f"  P_HP                 = {P_HP_sol/1e5:.3f} bar")
    print(f"  m_dot                = {m_dot_sol:.3f} kg/s")
    print(f"  Efficiency           = {eta_sol:.4f}")
    print(f"  N_pump               = {N_pump_sol:.1f} RPM"
          f"  (target: {N_pump_target} ± {tol_N_pump:.1f})")
    print(f"  SC                   = {SC_sol:.2f} K"
          f"  (target: {SC_cd_target} ± {tol_SC:.2f})")
         
#%% PSO
    
elif case_study == "Optimization_PSO":

    mdot = 320 # kg/s
    P_high = 153*1e5 # Pa
    
    # -------- TARGET SPECIFICATIONS --------
    N_pump_target = 2900   # RPM
    SC_cd_target  = 5.0    # K
    
    # -------- SCALING FACTORS --------
    P_scale    = 1e5
    mdot_scale = 1.0

    import pyswarms as ps
    from pyswarms.single import GlobalBestPSO
    import numpy as np
    import matplotlib.pyplot as plt
    from joblib import Parallel, delayed
    from tqdm import tqdm
    
    # -------- CONSTRAINT TOLERANCES --------
    tol_N_pump = 0.02 * N_pump_target
    tol_SC     = 0.1
    
    # -------- BOUNDS --------
    # pyswarms expects (min_bound, max_bound) as arrays of shape (n_dims,)
    P_lo,    P_hi    = 130, 160    # bar
    mdot_lo, mdot_hi = 250, 350    # kg/s
    
    min_bounds = np.array([P_lo,    mdot_lo])
    max_bounds = np.array([P_hi,    mdot_hi])
    bounds     = (min_bounds, max_bounds)
    
    # -------- SINGLE PARTICLE EVALUATION --------
    def evaluate_particle(x):
        """
        x : array [P_HP_bar, mdot_kgs]
        Returns fitness (scalar) — pyswarms minimises this.
        """
        P_high    = x[0] * P_scale
        m_dot_CO2 = x[1] * mdot_scale
    
        try:
            CO2_TC = CO2_OD_TC(m_dot_CO2, P_high)
            CO2_TC.mute_print()
            CO2_TC.solve(method='wegstein', max_iter=100)
        except Exception as e:
            return np.inf
    
        if not CO2_TC.converged:
            return np.inf
    
        N_pump_actual = CO2_TC.components["Pump"].model.W.N_rot
        SC_actual     = CO2_TC.components["Condenser"].model.ex_H.SC
        W_dot_turb    = CO2_TC.components["Turbine"].model.W.W_dot
        W_dot_pp      = CO2_TC.components["Pump"].model.W.W_dot
        Q_dot_gh      = CO2_TC.components["GasHeater"].model.Q.Q_dot
        eta           = (W_dot_turb - W_dot_pp) / Q_dot_gh
    
        if SC_actual is None or N_pump_actual is None:
            return 1e3
    
        # -------- PENALTY --------
        viol_N  = max(0, abs(N_pump_actual - N_pump_target) - tol_N_pump) / N_pump_target
        viol_SC = max(0, abs(SC_actual     - SC_cd_target)  - tol_SC)     / SC_cd_target
        
        penalty = 1 * (viol_N + viol_SC)
    
        fitness = -eta + penalty
    
        return fitness
    
    
    # -------- BATCH OBJECTIVE (pyswarms interface) --------
    # pyswarms passes the entire swarm at once: particles shape (n_particles, n_dims)
    # Must return array of shape (n_particles,)
    def batch_objective(particles):
        """Parallel evaluation of all particles in the swarm."""
        results = Parallel(n_jobs=-1, backend='loky', verbose=0, return_as='generator')(
            delayed(evaluate_particle)(particles[i])
            for i in range(len(particles))
        )
        fitnesses = np.array(list(tqdm(
            results,
            desc="Evaluating swarm",
            unit="particle",
            total=len(particles),
            leave=False,
        )))
        return fitnesses
    
    
    # -------- PSO HYPERPARAMETERS --------
    options = {
        'c1': 1.5,   # cognitive weight
        'c2': 1.5,   # social weight
        'w' : 0.7,   # inertia
    }
    
    # -------- INITIALISE OPTIMIZER --------
    optimizer = GlobalBestPSO(
        n_particles = 50,
        dimensions  = 2,
        options     = options,
        bounds      = bounds,
        init_pos    = None,   # replace with np.array([...]) to seed particles
    )
    
    # Seed one particle at best guess from 2D map if available
    try:
        init_pos = np.random.uniform(
            low  = min_bounds,
            high = max_bounds,
            size = (50, 2)
        )
        init_pos[0] = [P_sweep[j_best] / P_scale, mdot_sweep[i_best] / mdot_scale]
    
        optimizer = GlobalBestPSO(
            n_particles = 50,
            dimensions  = 2,
            options     = options,
            bounds      = bounds,
            init_pos    = init_pos,
        )
        print(f"  Seeded particle 0 at P={init_pos[0,0]:.2f} bar | mdot={init_pos[0,1]:.2f} kg/s")
    except NameError:
        print("  No 2D map available — random initialisation")
    
    # -------- RUN OPTIMISATION --------
    best_fitness, best_pos = optimizer.optimize(
        batch_objective,
        iters = 30,
    )
    
    # -------- RESULTS --------
    P_HP_sol  = best_pos[0] * P_scale
    m_dot_sol = best_pos[1] * mdot_scale
    
    # Re-evaluate best position to extract all outputs
    CO2_TC_sol = CO2_OD_TC(m_dot_sol, P_HP_sol)
    CO2_TC_sol.mute_print()
    CO2_TC_sol.solve(method='wegstein', max_iter=100)
    
    N_pump_sol = CO2_TC_sol.components["Pump"].model.W.N_rot
    SC_sol     = CO2_TC_sol.components["Condenser"].model.ex_H.SC
    W_dot_turb = CO2_TC_sol.components["Turbine"].model.W.W_dot
    W_dot_pp   = CO2_TC_sol.components["Pump"].model.W.W_dot
    Q_dot_gh   = CO2_TC_sol.components["GasHeater"].model.Q.Q_dot
    eta_sol    = (W_dot_turb - W_dot_pp) / Q_dot_gh
    
    print(f"\nPSO Optimisation Results")
    print(f"  P_HP       = {P_HP_sol/1e5:.3f} bar")
    print(f"  m_dot      = {m_dot_sol:.3f} kg/s")
    print(f"  Efficiency = {eta_sol:.4f}")
    print(f"  N_pump     = {N_pump_sol:.1f} RPM  (target: {N_pump_target} ± {tol_N_pump:.1f})")
    print(f"  SC         = {SC_sol:.2f} K    (target: {SC_cd_target} ± {tol_SC:.2f})")
    print(f"  Fitness    = {best_fitness:.4e}")
    
    # -------- CONVERGENCE PLOT --------
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(optimizer.cost_history, 'o-')
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best fitness [-]")
    ax.set_title("PSO Convergence")
    plt.tight_layout()
    plt.show()