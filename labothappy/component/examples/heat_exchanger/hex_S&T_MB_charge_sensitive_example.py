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

from component.heat_exchanger.hex_MB_charge_sensitive import HeatExchangerMB
# from toolbox.geometries.heat_exchanger.geometry_shell_and_tube_hx import ShellAndTubeGeom
from CoolProp.CoolProp import PropsSI    


#%%

import time
start_time = time.time()   

"------------ Shell and Tube HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HeatExchangerMB('Shell&Tube')

# -------------------------------------------------------------------------------------------------------------

# Methanol Sensible HT Case
# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = 'Methanol',
#     Hsu_T = 95 + 273.15, # K
#     Hsu_p = 10e5, # Pa
#     Hsu_m_dot = 27.8, # kg/s

#     # Second fluid
#     Csu_fluid = 'Water',
#     Csu_T = 25 + 273.15, # K
#     Csu_p = 5*1e5, # Pa
#     Csu_m_dot = 68.9, # kg/s  # Make sure to include fluid information
# )

# "Geometry Loading"


# params = {'n_series': 1,
#           'foul_t': 0.0002,
#           'foul_s': 0.00033,
#           'tube_cond': 50,
#           'Overdesign': 0,
#           'Shell_Side': 'H',
#           'Flow_Type': 'Shell&Tube',
#           'H_DP_ON': True,
#           'C_DP_ON': True,
#           'n_disc': 50,
#           'A_eff': 262.8,
#           'S_V_tot': 3.081237509643974,
#           'Shell_ID': 0.83,
#           'T_V_tot': 2.1307042370289317,
#           'Tube_L': 3.379,
#           'Tube_OD': 0.016,
#           'Tube_t': 0.016*0.1,
#           'central_spacing': 0.5,
#           'Tube_pass': 2,
#           'cross_passes': 6,
#           'n_tubes': 1567,
#           'pitch_ratio': 1.25,
#           'tube_layout': 60,
#           'Baffle_cut': 25}

# HX.set_parameters(
#     A_eff = params['A_eff'], Baffle_cut = params['Baffle_cut'], S_V_tot = params['S_V_tot'],
#     Shell_ID = params['Shell_ID'], T_V_tot = params['T_V_tot'], Tube_L = params['Tube_L'], 
#     Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
#     central_spacing = params['central_spacing'], cross_passes = params['cross_passes'], foul_s = params['foul_s'],
#     foul_t = params['foul_t'], n_series = params['n_series'], n_tubes = params['n_tubes'], 
#     pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],

#     Shell_Side = params['Shell_Side'],

#     Flow_Type = params['Flow_Type'], H_DP_ON = params['H_DP_ON'], C_DP_ON = params['C_DP_ON'], n_disc = params['n_disc'])

# # "Correlation Loading"

# Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
# Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

# Corr_H_DP = "Shell_Kern_DP"
# Corr_C_DP = "Gnielinski_DP"

# -------------------------------------------------------------------------------------------------------------

# Sizing code example case
# R134a evaporator case

HX.set_inputs(
              # Hot Fluid
              T_su_H = 273.15 + 26, # K
              P_su_H = 10*1e5, # 51.75*1e3, # Pa
              m_dot_H = 5.35, # kg/s
              fluid_H = 'Water',
              
              # Cold Fluid
              h_su_C = PropsSI('H','T', 273.15+7,'Q',0,'R134a')+1, # K
              P_su_C = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # 51.75*1e3, # Pa
              m_dot_C = 1.62, # kg/s
              fluid_C = 'R134a'
              )

"Geometry Loading"

params = {'n_series': 1,
          'foul_t': 0.000176,
          'foul_s': 0.000176,
          'tube_cond': 50,
          'Overdesign': 0,
          'Shell_Side': 'H',
          'Flow_Type': 'Shell&Tube',
          'H_DP_ON': True,
          'C_DP_ON': True,
          'n_disc': 30,
          'A_eff': 30.2,
          'S_V_tot': 3.081237509643974,
          'Shell_ID': 0.173,
          'T_V_tot': 2.1307042370289317,
          'Tube_L': 4.18,
          'Tube_OD': 0.01,
          'Tube_t': 0.001,
          'central_spacing': 0.381,
          'Tube_pass': 2,
          'cross_passes': 9,
          'n_tubes': 230,
          'pitch_ratio': 1.25,
          'tube_layout': 60,
          'Baffle_cut': 25}

HX.set_parameters(
    A_eff = params['A_eff'], Baffle_cut = params['Baffle_cut'], S_V_tot = params['S_V_tot'],
    Shell_ID = params['Shell_ID'], T_V_tot = params['T_V_tot'], Tube_L = params['Tube_L'], 
    Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
    central_spacing = params['central_spacing'], cross_passes = params['cross_passes'], foul_s = params['foul_s'],
    foul_t = params['foul_t'], n_series = params['n_series'], n_tubes = params['n_tubes'], 
    pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],

    Shell_Side = params['Shell_Side'],

    Flow_Type = params['Flow_Type'], H_DP_ON = params['H_DP_ON'], C_DP_ON = params['C_DP_ON'], n_disc = params['n_disc'])

# "Correlation Loading"

Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

Corr_H_DP = "Shell_Kern_DP"
Corr_C_DP = "Muller_Steinhagen_Heck_DP"

# -------------------------------------------------------------------------------------------------------------

# # Sizing code example case
# CO2 GasCooler Case

# HX.set_inputs(
#               # Hot Fluid
#               Hsu_T = 273.15 + 100, # K
#               Hsu_p = 100*1e5, # 51.75*1e3, # Pa
#               Hsu_m_dot = 100, # kg/s
#               Hsu_fluid = 'CO2',
              
#               # Cold Fluid
#               Csu_T = 273.15 + 20, # K
#               Csu_p = 5*1e5, # 51.75*1e3, # Pa
#               Csu_m_dot = 47.82, # kg/s
#               Csu_fluid = 'Water'
#               )

# "Geometry Loading"

# params = {'n_series': 1,
#         'foul_t': 0,
#         'foul_s': 0,
#         'tube_cond': 50,
#         'Overdesign': 0,
#         'Shell_Side': 'C',
#         'Flow_Type': 'Shell&Tube',
#         'H_DP_ON': True,
#         'C_DP_ON': True,
#         'n_disc': 30,
#         'A_eff': 856,
#         'S_V_tot': 1000,
#         'Shell_ID': 1.43,
#         'T_V_tot': 100,
#         'Tube_L': 9.21,
#         'Tube_OD': 0.02,
#         'Tube_t': 0.003,
#         'central_spacing': 0.5*1.43,
#         'Tube_pass': 1,
#         'cross_passes': 12,
#         'n_tubes': 1480,
#         'pitch_ratio': 1.5,
#         'tube_layout': 45,
#         'Baffle_cut': 25}

# HX.set_parameters(
#     A_eff = params['A_eff'], Baffle_cut = params['Baffle_cut'], S_V_tot = params['S_V_tot'],
#     Shell_ID = params['Shell_ID'], T_V_tot = params['T_V_tot'], Tube_L = params['Tube_L'], 
#     Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
#     central_spacing = params['central_spacing'], cross_passes = params['cross_passes'], foul_s = params['foul_s'],
#     foul_t = params['foul_t'], n_series = params['n_series'], n_tubes = params['n_tubes'], 
#     pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],

#     Shell_Side = params['Shell_Side'],

#     Flow_Type = params['Flow_Type'], H_DP_ON = params['H_DP_ON'], C_DP_ON = params['C_DP_ON'], n_disc = params['n_disc'])

# # "Correlation Loading"

# Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
# Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

# Corr_C_DP = "Shell_Kern_DP"
# Corr_H_DP = "Cheng_CO2_DP"

# -------------------------------------------------------------------------------------------------------------

# HX.set_inputs(
#               # Hot Fluid
#               Hsu_T = 273.15 + 105.4, # K
#               Hsu_p = 90*1e5, # 51.75*1e3, # Pa
#               Hsu_m_dot = 1590.82/4, # kg/s
#               Hsu_fluid = 'CO2',
              
#               # Cold Fluid
#               Csu_T = 273.15 + 25, # K
#               Csu_p = 5*1e5, # 51.75*1e3, # Pa
#               Csu_m_dot = 980/4, # kg/s
#               Csu_fluid = 'Water'
#               )

# params = {'htc_type': 'Correlation',
#          'A_eff': 6944.3238615708515,
#          'Baffle_cut': 26.726000000000003,
#          'S_V_tot': 16.218940644233115,
#          'Shell_ID': 1.6764,
#          'T_V_tot': 12.88333832771875,
#          'Tube_L': 14.84,
#          'Tube_OD': 0.009524999999999999,
#          'Tube_pass': 1,
#          'Tube_t': 0.0005587999999999999,
#          'central_spacing': 0.742,
#          'cross_passes': 19,
#          'foul_s': 0,
#          'foul_t': 0,
#          'n_series': 1,
#          'n_tubes': 15638,
#          'pitch_ratio': 1.33,
#          'tube_cond': 50,
#          'tube_layout': 60,
#          'Shell_Side': 'C',
#          'Flow_Type': 'Shell&Tube',
#          'H_DP_ON': True,
#          'C_DP_ON': True,
#          'n_disc': 1,
#          'DP_type': 'Correlation',
#          'Overdesign': 0}

# HX.set_parameters(
#     A_eff = params['A_eff'], Baffle_cut = params['Baffle_cut'], S_V_tot = params['S_V_tot'],
#     Shell_ID = params['Shell_ID'], T_V_tot = params['T_V_tot'], Tube_L = params['Tube_L'], 
#     Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
#     central_spacing = params['central_spacing'], cross_passes = params['cross_passes'], foul_s = params['foul_s'],
#     foul_t = params['foul_t'], n_series = params['n_series'], n_tubes = params['n_tubes'], 
#     pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],

#     Shell_Side = params['Shell_Side'],

#     Flow_Type = params['Flow_Type'], H_DP_ON = params['H_DP_ON'], C_DP_ON = params['C_DP_ON'], n_disc = params['n_disc'])

# # "Correlation Loading"

# Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
# Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

# Corr_C_DP = "Shell_Kern_DP"
# Corr_H_DP = "Cheng_CO2_DP"

# -------------------------------------------------------------------------------------------------------------

"Parameters Setting"

HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

# UD_H_HTC = {'Liquid': 1500,
#             'Vapor' : 100,
#             'Two-Phase' : 1000,
#             'Vapor-wet' : 100,
#             'Dryout' : 1000,
#             'Transcritical' : 1500}

# UD_C_HTC = {'Liquid': 1350,
#             'Vapor' : 100,
#             'Two-Phase' : 1000,
#             'Vapor-wet' : 100,
#             'Dryout' : 10000,
#             'Transcritical' : 1350}

# HX.set_htc(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'

# HX.set_parameters(
#     A_eff = HX_geom.A_eff, Baffle_cut = HX_geom.Baffle_cut, D_OTL = HX_geom.D_OTL, N_strips = HX_geom.N_strips, S_V_tot = HX_geom.S_V_tot, # 5
#     Shell_ID = HX_geom.Shell_ID, T_V_tot = HX_geom.T_V_tot, Tube_L = HX_geom.Tube_L, Tube_OD = HX_geom.Tube_OD, Tube_pass = HX_geom.Tube_pass, # 10
#     Tube_t = HX_geom.Tube_t, Tubesheet_t = HX_geom.Tubesheet_t, central_spacing = HX_geom.central_spacing, clear_BS = HX_geom.clear_BS, clear_TB = HX_geom.clear_TB, # 15
#     cross_passes = HX_geom.cross_passes, foul_s = HX_geom.foul_s, foul_t = HX_geom.foul_t, inlet_spacing = HX_geom.inlet_spacing, n_series = HX_geom.n_series, # 20
#     n_tubes = HX_geom.n_tubes, outlet_spacing = HX_geom.outlet_spacing, pitch_ratio = HX_geom.pitch_ratio, tube_cond = HX_geom.tube_cond, tube_layout = HX_geom.tube_layout, # 25

#     Shell_Side = 'H', # 26
#     Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True, n_disc = 50) # 30

# Corr_H_DP = "Shell_Bell_Delaware_DP"
# Corr_C_DP = "Gnielinski_DP"

HX.set_DP(DP_type="Correlation", Corr_H=Corr_H_DP, Corr_C=Corr_C_DP)
# HX.set_DP()

"Solve the component"

HX.solve()
HX.plot_cells()
