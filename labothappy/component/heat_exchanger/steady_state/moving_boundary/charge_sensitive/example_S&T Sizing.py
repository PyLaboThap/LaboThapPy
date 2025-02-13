"""
Supplemental code for paper:
I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of 
Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014
"""

"""
Modification w/r to previous version:
    - Putting some order in the Objective Function "for" loops. Sparing some
    lines of code.
    - x_di_c correct calculation.
"""

# from __future__ import division, print_function

from simulation_model_AS_DP import HeatExchangerMB
from modules.geometry_shell_and_tube_hx import ShellAndTubeGeom
from CoolProp.CoolProp import PropsSI    
import numpy as np

#%%

import time
start_time = time.time()   

"------------ Shell and Tube HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HeatExchangerMB('Shell&Tube')

HX.set_inputs(
    # First fluid
    Hsu_fluid = 'Cyclopentane',
    Hsu_T = 273.15 + 38.43, # K
    Hsu_p = 51.75*1e3, # Pa
    Hsu_m_dot = 32.85, # kg/s

    # Second fluid
    Csu_fluid = 'Water',
    Csu_T = 273.15 + 24, # K
    Csu_p = 101.3*1e3, # Pa
    Csu_m_dot = 1000, # kg/s  # Make sure to include fluid information
)

"Correlation Loading And Setting"

Corr_H = {"1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}
Corr_C = {"1P" : "shell_htc_kern", "2P" : "shell_htc_kern"}

HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

"Parameters Setting"

params = {'tube_layout': 45, 'Tube_pass': 2, 'n_series': 1, 'Baffle_cut': 25, 'foul_t': 0, 'foul_s': 0, 'tube_cond': 50, 'Shell_Side': 'C', 'Flow_Type': 'Shell&Tube', 'H_DP_ON': True, 'C_DP_ON': True, 'n_disc': 30, 'A_eff': 461.22350065882466, 'S_V_tot': 3.2783973897608107, 'Shell_ID': 0.889, 'T_V_tot': np.float64(1.7905077622091865), 'Tube_L': 20, 'Tube_OD': 0.0254, 'Tube_t': np.float64(0.00277), 'central_spacing': np.float64(1.033), 'cross_passes': 9, 'n_tubes': 289, 'pitch_ratio': 1.25}

HX.set_parameters(
    A_eff = params['A_eff'], Baffle_cut = params['Baffle_cut'], S_V_tot = params['S_V_tot'],
    Shell_ID = params['Shell_ID'], T_V_tot = params['T_V_tot'], Tube_L = params['Tube_L'], 
    Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
    central_spacing = params['central_spacing'], cross_passes = params['cross_passes'], foul_s = params['foul_s'],
    foul_t = params['foul_t'], n_series = params['n_series'], n_tubes = params['n_tubes'], 
    pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],

    Shell_Side = params['Shell_Side'],

    Flow_Type = params['Flow_Type'], H_DP_ON = params['H_DP_ON'], C_DP_ON = params['C_DP_ON'], n_disc = params['n_disc']) 

"Correlation Loading And Setting"

Corr_H_DP = "pipe_internal_DP"
Corr_C_DP = "shell_DP_kern"

HX.set_DP(DP_type="User-Defined", UD_H_DP=1e3, UD_C_DP=1e3)
# HX.set_DP(DP_type = "Correlation", Corr_H=Corr_H_DP, Corr_C=Corr_C_DP)

HX.solve()

