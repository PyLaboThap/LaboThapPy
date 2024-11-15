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

from simulation_model import HeatExchangerMB
from modules.geometry_shell_and_tube_hx import ShellAndTubeGeom
from CoolProp.CoolProp import PropsSI    


#%%

import time
start_time = time.time()   

"------------ Shell and Tube HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HeatExchangerMB('Shell&Tube')

# "Setting inputs"

# -------------------------------------------------------------------------------------------------------------

# DECAGONE Evaporator case
HX.set_inputs(
    # First fluid
    Hsu_fluid = 'INCOMP::S800',
    Hsu_T = 310 + 273.15, # K
    Hsu_p = 3.25*1e5, # Pa
    Hsu_m_dot = 19.42, # kg/s

    # Second fluid
    Csu_fluid = 'Cyclopentane',
    Csu_T = 95.1 + 273.15, # K
    Csu_p = 31.5*1e5, # Pa
    Csu_m_dot = 13.84, # kg/s  # Make sure to include fluid information
)

"Geometry Loading"

HX_geom = ShellAndTubeGeom()
HX_geom.set_parameters("DECAGONE_EVAP_Equ") 

"Correlation Loading"

Corr_H = {"1P" : "Shell_Bell_Delaware_HTC", "2P" : "Shell_Bell_Delaware_HTC"}
Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

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

"Parameters Setting"

HX.set_parameters(
    A_eff = HX_geom.A_eff, Baffle_cut = HX_geom.Baffle_cut, D_OTL = HX_geom.D_OTL, N_strips = HX_geom.N_strips, S_V_tot = HX_geom.S_V_tot, # 5
    Shell_ID = HX_geom.Shell_ID, T_V_tot = HX_geom.T_V_tot, Tube_L = HX_geom.Tube_L, Tube_OD = HX_geom.Tube_OD, Tube_pass = HX_geom.Tube_pass, # 10
    Tube_t = HX_geom.Tube_t, Tubesheet_t = HX_geom.Tubesheet_t, central_spacing = HX_geom.central_spacing, clear_BS = HX_geom.clear_BS, clear_TB = HX_geom.clear_TB, # 15
    cross_passes = HX_geom.cross_passes, foul_s = HX_geom.foul_s, foul_t = HX_geom.foul_t, inlet_spacing = HX_geom.inlet_spacing, n_series = HX_geom.n_series, # 20
    n_tubes = HX_geom.n_tubes, outlet_spacing = HX_geom.outlet_spacing, pitch_ratio = HX_geom.pitch_ratio, tube_cond = HX_geom.tube_cond, tube_layout = HX_geom.tube_layout, # 25

    Shell_Side = 'C', # 26

    Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True, n_disc = 100) # 30


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

HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31
# HX.set_HTC(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'

"Solve the component"
HX.solve()

HX.plot_cells()