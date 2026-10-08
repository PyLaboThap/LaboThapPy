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
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.plot_SPHE import plot_spiral

#%%

import time
start_time = time.time()   

"--------- Tube and Fins HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HexMBChargeSensitive('SPHE')

# "Setting inputs"

# -------------------------------------------------------------------------------------------------------------

# # DECAGONE Recuperator HTX case

HX.set_inputs(
    # First fluid
    fluid_H = 'Propane',
    T_su_H = 72.4 + 273.15, # K
    P_su_H = 19.6*1e5, # Pa
    m_dot_H = 0.0192, # kg/s
    
    # Second fluid
    fluid_C = 'Water',
    T_su_C = 50 + 273.15, # K
    P_su_C = 5*1e5, # Pa
    m_dot_C = 0.286, # kg/s  # Make sure to include fluid information
)

"Correlation Loading"

# Corr_H = {"1P" : "Tube_And_Fins", "2P" : "ext_tube_film_condens"}
Corr_H = {"1P" : "Tube_And_Fins", "2P" : "Tube_And_Fins"}
Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}

Corr_H_DP = {"1P" : "Tube_And_Fins_DP", "2P" : "Tube_And_Fins_DP"}
Corr_C_DP = {"1P" : "Gnielinski_DP", "2P" : "Choi_DP"}

# -------------------------------------------------------------------------------------------------------------

"Parameters Setting"

params = {
        'C_canal_t' : 0.003, # [m]
        'D_ext' : 0.2, # [m]
        'H' : 0.1, # [m]
        'H_canal_t' : 0.003, # [m]
        'r0_in' : 1*1e-3, # [m]
        't' : 0.001, # [m]
        
        'inner_channel' : 'cold'
        }

HX.set_parameters(
    C_canal_t = params['C_canal_t'], D_ext = params['D_ext'], H = params['H'], 
    H_canal_t = params['H_canal_t'], r0_in = params['r0_in'], t = params['t'],

    inner_channel = params['inner_channel'], Flow_Type = "CounterFlow", n_disc = 30) # 32

# User defined values

UD_H_HTC = {'Liquid': 5000,
            'Vapor' : 1000,
            'Two-Phase' : 10000,
            'Vapor-wet' : 10000,
            'Dryout' : 10000,
            'Transcritical' : 5000}

UD_C_HTC = {'Liquid': 5000,
            'Vapor' : 1000,
            'Two-Phase' : 10000,
            'Vapor-wet' : 10000,
            'Dryout' : 10000,
            'Transcritical' : 5000}

HX.set_htc(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'
# HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 

HX.set_DP() # equivalent to HX.set_DP(DP_type = None)
# HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
# HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
# HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

"Solve the component"
HX.solve()
HX.plot_cells()

"Plot the spiral"

plot_spiral(HX)
