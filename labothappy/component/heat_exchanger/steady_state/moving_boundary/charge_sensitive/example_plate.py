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
from modules.geometry_plate_hx_swep import PlateGeomSWEP

#%%

# import time
# start_time = time.time()   

# "------------ Plate HX -------------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HeatExchangerMB('Plate')

# "Setting inputs"

# # ---------------------------------------------------------------------------------------------------------

# # Condenser Case
# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = 'Cyclopentane',
#     Hsu_T = 139 + 273.15, # K
#     Hsu_p = 0.77*1e5, # Pa
#     Hsu_m_dot = 0.014, # kg/s

#     # Second fluid
#     Csu_fluid = 'Water',
#     Csu_T = 12 + 273.15, # K
#     Csu_p = 5*1e5, # Pa
#     Csu_m_dot = 0.2, # kg/s  # Make sure to include fluid information
# )

# "Geometry Loading"

# HX_geom = PlateGeomSWEP()
# HX_geom.set_parameters("B20Hx24/1P") 

# Corr_H = {"1P" : "Gnielinski", "2P" : "Han_cond_BPHEX"}
# Corr_C = {"1P" : "Gnielinski", "2P" : "Han_Boiling_BPHEX_HTC"}

# ---------------------------------------------------------------------------------------------------------

# Evaporator Case
HX.set_inputs(
    # First fluid
    Hsu_fluid = 'INCOMP::T66',
    Hsu_T = 243 + 273.15, # K
    Hsu_p = 5*1e5, # Pa
    Hsu_m_dot = 0.4, # kg/s

    # Second fluid
    Csu_fluid = 'Cyclopentane',
    Csu_T = 41 + 273.15, # K
    Csu_p = 31.5*1e5, # Pa
    Csu_m_dot = 0.014, # kg/s  # Make sure to include fluid information
)

"Geometry Loading"

HX_geom = PlateGeomSWEP()
HX_geom.set_parameters("B20Hx24/1P") 

Corr_H = {"1P" : "Gnielinski", "2P" : "Han_cond_BPHEX"}
Corr_C = {"1P" : "Gnielinski", "2P" : "Han_Boiling_BPHEX_HTC"}

# ---------------------------------------------------------------------------------------------------------

# # Desuperheater Case
# HX.set_inputs(
#     # First fluid
#     Hsu_fluid = 'Cyclopentane',
#     Hsu_T = 205 + 273.15, # K
#     Hsu_p = 1*1e5, # Pa
#     Hsu_m_dot = 0.014, # kg/s

#     # Second fluid
#     Csu_fluid = 'Water',
#     Csu_T = 12 + 273.15, # K
#     Csu_p = 4*1e5, # Pa
#     Csu_m_dot = 0.2, # kg/s  # Make sure to include fluid information
# )

# "Geometry Loading"

# HX_geom = PlateGeomSWEP()
# HX_geom.set_parameters("B35TM0x10/1P") 

# Corr_H = {"1P" : "Gnielinski", "2P" : "Han_cond_BPHEX"}
# Corr_C = {"1P" : "Gnielinski", "2P" : "Han_Boiling_BPHEX_HTC"}

# ---------------------------------------------------------------------------------------------------------

"Parameters Setting"

HX.set_parameters(
    A_c = HX_geom.A_c, A_h = HX_geom.A_h, h = HX_geom.h, l = HX_geom.l, l_v = HX_geom.l_v, # 5
    C_CS = HX_geom.C_CS, C_Dh = HX_geom.C_Dh, C_V_tot = HX_geom.C_V_tot, C_canal_t = HX_geom.C_canal_t, C_n_canals = HX_geom.C_n_canals, # 10
    H_CS = HX_geom.H_CS, H_Dh = HX_geom.H_Dh, H_V_tot = HX_geom.H_V_tot, H_canal_t = HX_geom.H_canal_t, H_n_canals = HX_geom.H_n_canals, # 15
    casing_t = HX_geom.casing_t, chevron_angle = HX_geom.chevron_angle, fooling = HX_geom.fooling, # 18
    n_plates = HX_geom.n_plates, plate_cond = HX_geom.plate_cond, plate_pitch_co = HX_geom.plate_pitch_co, t_plates = HX_geom.t_plates, w = HX_geom.w, # 23
    
    Flow_Type = 'CounterFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 50) # 27

UD_H_HTC = {'Liquid':100,
            'Vapor' : 100,
            'Two-Phase' : 1000,
            'Vapor-wet' : 100,
            'Dryout' : 1000,
            'Transcritical' : 200}

UD_C_HTC = {'Liquid':100,
            'Vapor' : 100,
            'Two-Phase' : 1000,
            'Vapor-wet' : 100,
            'Dryout' : 10000,
            'Transcritical' : 200}

# HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 28
HX.set_htc(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'
HX.set_DP()

"Solve the component"
HX.solve()
HX.plot_cells()