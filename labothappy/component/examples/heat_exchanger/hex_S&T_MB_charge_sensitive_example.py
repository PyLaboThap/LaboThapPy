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
from toolbox.geometries.heat_exchanger.geometry_shell_and_tube_hx import ShellAndTubeGeom
from CoolProp.CoolProp import PropsSI    


#%%

import time
start_time = time.time()   

"------------ Shell and Tube HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HeatExchangerMB('Shell&Tube')

# # "Setting inputs"
# # -------------------------------------------------------------------------------------------------------------
# # DECAGONE Evaporator case
# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = 'INCOMP::T66',
#     Hsu_T = 310 + 273.15, # K
#     Hsu_p = 3.25*1e5, # Pa
#     Hsu_m_dot = 19.42, # kg/s

#     # Second fluid
#     Csu_fluid = 'Cyclopentane',
#     Csu_T = 95.1 + 273.15, # 95.1 + 273.15, # K
#     Csu_p = 31.5*1e5, # 31.5*1e5, # Pa
#     Csu_m_dot = 13.84, # 13.84, # kg/s  # Make sure to include fluid informastion
# )

# "Geometry Loading"

# HX_geom = ShellAndTubeGeom()
# HX_geom.set_parameters("DECAGONE_EVAP_Equ") 

# "Correlation Loading"

# # Corr_H = {"1P" : "Shell_Bell_Delaware_HTC", "2P" : "Shell_Bell_Delaware_HTC"}
# Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
# Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

# -------------------------------------------------------------------------------------------------------------

# Sam HTX case
# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = 'Water',
#     Hsu_T = 90 + 273.15, # K
#     Hsu_p = 3*1e5, # Pa
#     Hsu_m_dot = 5.5188, # kg/s

#     # Second fluid
#     Csu_fluid = 'Water',
#     Csu_T = 10 + 273.15, # K
#     Csu_p = 1e5, # Pa
#     Csu_m_dot = 80, # kg/s  # Make sure to include fluid information
# )

# "Geometry Loading"

# HX_geom = ShellAndTubeGeom()
# HX_geom.set_parameters("Valid_Sam") 

# "Correlation Loading"

# Corr_H = {"1P" : "Shell_Bell_Delaware_HTC", "2P" : "Shell_Bell_Delaware_HTC"}
# Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

# -------------------------------------------------------------------------------------------------------------

# # Kakac Condenser HTX case
# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = 'R22',
#     Hsu_T = 37.1 + 273.15, # K
#     Hsu_p = PropsSI("P","T",37 + 273.15,"Q",0,"R22"), # Pa
#     Hsu_m_dot = 0.737, # kg/s

#     # Second fluid
#     Csu_fluid = 'Water',
#     Csu_T = 18 + 273.15, # K
#     Csu_p = 4*1e5, # Pa
#     Csu_m_dot = 1.36, # kg/s  # Make sure to include fluid information
# )

# "Geometry Loading"

# HX_geom = ShellAndTubeGeom()
# HX_geom.set_parameters("kakac_Ex_9_3") 

# "Correlation Loading"

# Corr_H = {"1P" : "Shell_Bell_Delaware_HTC", "2P" : "ext_tube_film_condens"}
# Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

# -------------------------------------------------------------------------------------------------------------

# # Sizing code example case

HX.set_inputs(
              # Hot Fluid
              Hsu_T = 273.15 + 26, # K
              Hsu_p = 2*1e5, # 51.75*1e3, # Pa
              Hsu_m_dot = 5.35, # kg/s
              Hsu_fluid = 'Water',
              
              # Cold Fluid
              Csu_h = PropsSI('H','T', 273.15+7,'Q',0,'R134a'), # K
              Csu_p = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # 51.75*1e3, # Pa
              Csu_m_dot = 1.62, # kg/s
              Csu_fluid = 'R134a'
              )

# # # "Geometry Loading"

# # params =  {'tube_layout': 45, 'Tube_pass': 2, 'n_series': 1, 'Baffle_cut': 25, 'foul_t': 0, 'foul_s': 0, 'tube_cond': 50, 'Shell_Side': 'C', 'Flow_Type': 'Shell&Tube', 'H_DP_ON': True, 'C_DP_ON': True, 'n_disc': 30, 'A_eff': 697.3292368049083, 'S_V_tot': 10.268864099111516, 'Shell_ID': 1.524, 'T_V_tot': 2.42544652488181, 'Tube_L': 7.45, 'Tube_OD': 0.038099999999999995, 'Tube_t': 0.00277, 'central_spacing': 0.745, 'cross_passes': 9, 'n_tubes': 391, 'pitch_ratio': 1.25}
params = {'n_series': 1,
         'foul_t': 0.000176,
         'foul_s': 0.000176,
         'tube_cond': 50,
         'Overdesign': 0,
         'Shell_Side': 'H',
         'Flow_Type': 'Shell&Tube',
         'H_DP_ON': True,
         'C_DP_ON': True,
         'n_disc': 5,
         'A_eff': 548.8559657840013,
         'S_V_tot': 3.081237509643974,
         'Shell_ID': 0.8382,
         'T_V_tot': 2.1307042370289317,
         'Tube_L': 11.9,
         'Tube_OD': 0.0254,
         'Tube_t': 0.00277,
         'central_spacing': 1.106,
         'Tube_pass': 2,
         'cross_passes': 10,
         'n_tubes': 578,
         'pitch_ratio': 1.25,
         'tube_layout': 60,
         'Baffle_cut': 36.748}

HX.set_parameters(
    A_eff = params['A_eff'], Baffle_cut = params['Baffle_cut'], S_V_tot = params['S_V_tot'],
    Shell_ID = params['Shell_ID'], T_V_tot = params['T_V_tot'], Tube_L = params['Tube_L'], 
    Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
    central_spacing = params['central_spacing'], cross_passes = params['cross_passes'], foul_s = params['foul_s'],
    foul_t = params['foul_t'], n_series = params['n_series'], n_tubes = params['n_tubes'], 
    pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],

    Shell_Side = params['Shell_Side'],

    Flow_Type = params['Flow_Type'], H_DP_ON = params['H_DP_ON'], C_DP_ON = params['C_DP_ON'], n_disc = params['n_disc'])

# # "Correlation Loading"

Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}
Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

Corr_H_DP = "Shell_Kern_DP"
Corr_C_DP = "Gnielinski_DP"

# -------------------------------------------------------------------------------------------------------------

"Parameters Setting"

HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

# UD_H_HTC = {'Liquid':100,
#             'Vapor' : 100,
#             'Two-Phase' : 1000,
#             'Vapor-wet' : 100,
#             'Dryout' : 1000,
#             'Transcritical' : 200}

# UD_C_HTC = {'Liquid':100,
#             'Vapor' : 100,
#             'Two-Phase' : 1000,
#             'Vapor-wet' : 100,
#             'Dryout' : 10000,
#             'Transcritical' : 200}

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
