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
from toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP

#%%

# import time
# start_time = time.time()   

# "------------ Plate HX -------------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HeatExchangerMB('PCHE')

# # "Setting inputs"

# # # ---------------------------------------------------------------------------------------------------------

# Case from:  
# Numerical modelling and transient analysis of a printed circuit heat exchanger 
# used as recuperator for supercritical CO2 heat to power conversion systems

HX.set_inputs(
    # First fluid
    fluid_H = 'CO2',
    T_su_H = 249 + 273.15, # K
    P_su_H = 96.4*1e5, # Pa
    m_dot_H = 5.35, # kg/s

    # Second fluid
    fluid_C = 'CO2',
    T_su_C = 52.77 + 273.15, # K
    P_su_C = 165.4*1e5, # Pa
    m_dot_C = 5.35, # kg/s  # Make sure to include fluid information
)

"Geometry Loading"
params = {'alpha': 40, # Channel zigzag angle
          'D_c': 2*1e-3, # Channel diameter
          'C_V_tot' : 1, 
          'H_V_tot' : 1, 
          'k_cond': 60, # plate conductivity
          'L_c': 1.5, # channel length
          'N_c': 112, # n channels per plate
          'N_p': 256, # n plates
          'R_p': 1, # n_hot_channel_row / n_cold_channel_row
          't_2': 0.63*1e-3, # Horizontal pitch
          't_3': 1*1e-3} # Plate_thickness

Corr_H = {"SC" : "Gnielinski", "1P" : "Gnielinski"}
Corr_C = {"SC" : "Gnielinski", "1P" : "Gnielinski"}

H_DP = "Gnielinski_DP"
C_DP = "Gnielinski_DP"
    
# # ---------------------------------------------------------------------------------------------------------

# "Parameters Setting"

HX.set_parameters(
    alpha = params['alpha'], C_V_tot = params['C_V_tot'], H_V_tot = params['H_V_tot'], D_c = params['D_c'], k_cond = params['k_cond'], L_c = params['L_c'], 
    N_c = params['N_c'], N_p = params['N_p'], R_p = params['R_p'], t_2 = params['t_2'], t_3 = params['t_3'],
    
    Flow_Type = 'CounterFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 50) # 27

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

HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 28
# HX.set_htc(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'
HX.set_DP(DP_type= 'Correlation', Corr_H = H_DP, Corr_C = C_DP)

# "Solve the component"
HX.solve()
# HX.plot_cells()

