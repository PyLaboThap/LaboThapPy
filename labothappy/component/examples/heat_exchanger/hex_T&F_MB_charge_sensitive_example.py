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
import __init__
from component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive

#%%

import time
start_time = time.time()   

"--------- Tube and Fins HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HexMBChargeSensitive('Tube&Fins')

# "Setting inputs"

# -------------------------------------------------------------------------------------------------------------

# # DECAGONE Recuperator HTX case

HX.set_inputs(
    # First fluid
    fluid_H = 'Cyclopentane',
    T_su_H = 133.8 + 273.15, # K
    P_su_H = 0.8*1e5, # Pa
    m_dot_H = 13.76, # kg/s

    # Second fluid
    fluid_C = 'Cyclopentane',
    T_su_C = 35.1 + 273.15, # K
    P_su_C = 31.5*1e5, # Pa
    m_dot_C = 13.76, # kg/s  # Make sure to include fluid information
)

"Correlation Loading"

# Corr_H = {"1P" : "Tube_And_Fins", "2P" : "ext_tube_film_condens"}
Corr_H = {"1P" : "Tube_And_Fins", "2P" : "Tube_And_Fins"}
Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

Corr_H_DP = {"1P" : "Tube_And_Fins_DP", "2P" : "Tube_And_Fins_DP"}
Corr_C_DP = {"1P" : "Gnielinski_DP", "2P" : "Choi_DP"}

# -------------------------------------------------------------------------------------------------------------

# DECAGONE ACC HTX case

# HX.set_inputs(
#     # First fluid
#     fluid_H = 'Cyclopentane',
#     T_su_H = 53.6 + 273.15, # K
#     P_su_H = 0.7*1e5, # Pa
#     m_dot_H = 13.8/2, # kg/s

#     # Second fluid
#     Csu_fluid = 'Air',
#     Csu_T = 12 + 273.15, # K
#     Csu_p = 1.05*1e5, # Pa
#     Csu_m_dot = 158.5, # kg/s  # Make sure to include fluid information
# )

# "Geometry Loading"

# HX_geom = TubeAndFinsGeom()
# HX_geom.set_parameters("DECAGONE_ACC") 

# Fin_Side = 'C'

# "Correlation Loading"

# Corr_H = {"1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}
# Corr_C = {"1P" : "Tube_And_Fins", "2P" : "Boiling_curve"}

# -------------------------------------------------------------------------------------------------------------

"Parameters Setting"

params = {
        'A_flow' : 3.36, # [m2]    
        'fouling' : 0.9, # [(m^2*K)/W]
        'Fin_OD' : 0.0350625, # [m]
        'Fin_per_m' : 400, # [-]
        'Fin_t' : 0.00025, # [m]
        'k_fin' : 50, # [W/(m*K)] : Steel
        'Fin_type' : "Square", # [-]
        'Finned_tube_flag' : 1, # [-]
        'h' : 1.6, # [m]
        'Tube_pass' : 6, # [-]
        'n_rows' : 4, # [-]
        'n_series' : 1, # [-]
        'n_parallel' : 1, # [-]
        'pitch_ratio' : 2.125, # [-]
        'tube_arrang' : "Staggered", # [-]
        'Tube_cond' : 50, # [(m^2*K)/W]
        'Tube_L' : 4.5, # [m]
        'Tube_OD' : 0.0165, # [m]        
        'Tube_t' : 0.0015, # [m]        
        'n_tubes' : 480, # [-]
        'w' : 0.787, # [m]
        
        'Fin_Side' : 'H'
        }

HX.set_parameters(
    A_flow = params['A_flow'], Fin_OD = params['Fin_OD'], Fin_per_m = params['Fin_per_m'], Fin_t = params['Fin_t'], Fin_type = params['Fin_type'], # 10
    Finned_tube_flag = params['Tube_t'], Tube_L = params['Tube_L'], Tube_OD = params['Tube_OD'], # 15
    Tube_cond = params['Tube_cond'], Tube_t = params['Tube_t'], fouling = params['fouling'], h = params['h'], k_fin = params['k_fin'], # 20
    Tube_pass = params['Tube_pass'], n_rows = params['n_rows'], n_series = params['n_series'], n_tubes = params['n_tubes'], pitch_ratio = params['pitch_ratio'], tube_arrang = params['tube_arrang'], w = params['w'], # 30

    Fin_Side = params['Fin_Side'], # 28

    Flow_Type = 'CrossFlow', n_disc = 50) # 32

# User defined values

# UD_H_HTC = {'Liquid':29,
#             'Vapor' : 29,
#             'Two-Phase' : 29,
#             'Vapor-wet' : 29,
#             'Dryout' : 29,
#             'Transcritical' : 29}

# UD_C_HTC = {'Liquid':10000,
#             'Vapor' : 10000,
#             'Two-Phase' : 10000,
#             'Vapor-wet' : 10000,
#             'Dryout' : 10000,
#             'Transcritical' : 10000}

# HX.set_HTC(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'
HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 

# HX.set_DP() # equivalent to HX.set_DP(DP_type = None)
# HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
# HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

"Solve the component"
HX.solve()
HX.plot_cells()