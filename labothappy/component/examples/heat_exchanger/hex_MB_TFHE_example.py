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
# import __init__
from labothappy.component.heat_exchanger.hex_MB import HexMB

#%%

import time
start_time = time.time()   

"--------- Tube and Fins HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

case_study = "REC"

HX = HexMB('TFHE')

# "Setting inputs"

# -------------------------------------------------------------------------------------------------------------
if case_study == "REC":

    # DECAGONE Recuperator HTX case
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
    # Fin side (hot): Pismennyi on the gas; for a two-phase / vapor-wet cell the same correlation is reused
    # (phase change ignored, warning). Tube side (cold): pipe correlations, None = automatic choice
    # (cyclopentane = hydrocarbon: Steiner-Taborek / VDI / Gungor-Winterton for boiling).
    Corr_H = {"1P" : "Pismennyi", "2P" : "Pismennyi"}
    Corr_C = {"1P" : "Gnielinski", "2P" : None}
    
    # Pressure drop correlations do not exist yet for the finned side: use DP_type = "User-Defined" below
    Corr_H_DP = {"1P" : "Pismennyi", "2P" : "Pismennyi"}
    Corr_C_DP = {"1P" : "Churchill", "2P" : "Friedel"}

# -------------------------------------------------------------------------------------------------------------

# elif case_study == "REC":

#     # DECAGONE ACC HTX case
    
#     HX.set_inputs(
#         # First fluid
#         fluid_H = 'Cyclopentane',
#         T_su_H = 53.6 + 273.15, # K
#         P_su_H = 0.7*1e5, # Pa
#         m_dot_H = 13.8/2, # kg/s
    
#         # Second fluid
#         Csu_fluid = 'Air',
#         Csu_T = 12 + 273.15, # K
#         Csu_p = 1.05*1e5, # Pa
#         Csu_m_dot = 158.5, # kg/s  # Make sure to include fluid information
#     )
    
#     "Geometry Loading"
    
#     HX_geom = TubeAndFinsGeom()
#     HX_geom.set_parameters("DECAGONE_ACC") 
    
#     fin_side = 'C'
    
#     "Correlation Loading"
    
#     Corr_H = {"1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}
#     Corr_C = {"1P" : "Tube_And_Fins", "2P" : "Boiling_curve"}

# -------------------------------------------------------------------------------------------------------------

"Parameters Setting"

params = {
        'A_flow' : 3.36, # [m2]    
        'R_f' : 0.9, # [(m^2*K)/W]
        'D_f_o' : 0.0350625, # [m]
        'f_density' : 400, # [1/m]
        't_f' : 0.00025, # [m]
        'k_f' : 50, # [W/(m*K)] : Steel
        'fin_type' : "Square", # [-]
        'finned' : 1, # [-]
        'h' : 1.6, # [m]
        'n_passes' : 6, # [-]
        'n_rows' : 4, # [-]
        'n_series' : 1, # [-]
        'n_parallel' : 1, # [-]
        'PR' : 2.125, # [-]
        'theta' : 45, # [deg] tube layout: 30, 45 or 60 = staggered, 0 or 90 = inline
        'k_wall' : 50, # [W/(m*K)]
        'L_t' : 4.5, # [m]
        'D_t_o' : 0.0165, # [m]        
        't_t' : 0.0015, # [m]        
        'n_t' : 480, # [-]
        'w' : 0.787, # [m]
        
        'fin_side' : 'H'
        }

HX.set_parameters(
    A_flow = params['A_flow'], D_f_o = params['D_f_o'], f_density = params['f_density'], t_f = params['t_f'], fin_type = params['fin_type'], # 10
    finned = params['finned'], L_t = params['L_t'], D_t_o = params['D_t_o'], # 15
    k_wall = params['k_wall'], t_t = params['t_t'], R_f = params['R_f'], h = params['h'], k_f = params['k_f'], # 20
    n_passes = params['n_passes'], n_rows = params['n_rows'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], PR = params['PR'], theta = params['theta'], w = params['w'], # 30

    fin_side = params['fin_side'], # 28

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
# HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

"Solve the component"
HX.solve()
HX.plot_cells()