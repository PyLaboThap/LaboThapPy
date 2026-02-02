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

from component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive

#%%

# import time
# start_time = time.time()

# "------------ Plate HX -------------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HexMBChargeSensitive('PCHE')

test_case = "test_CO2"

# # "Setting inputs"

# # # ---------------------------------------------------------------------------------------------------------

# Case from:  
# Numerical modelling and transient analysis of a printed circuit heat exchanger 
# used as recuperator for supercritical CO2 heat to power conversion systems

if test_case == "test_CO2":
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
              't_3': 1*1e-3, # Plate thickness
              'type_channel' : 'Zigzag'} 
    
    Corr_H = {"SC" : "Gnielinski", "1P" : "Gnielinski"}
    Corr_C = {"SC" : "Gnielinski", "1P" : "Gnielinski"}

    Corr_H_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach"}
    Corr_C_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach"}  

    # ---------------------------------------------------------------------------------------------------------
    # "Parameters Setting"
    
    HX.set_parameters(
        alpha = params['alpha'], C_V_tot = params['C_V_tot'], H_V_tot = params['H_V_tot'], D_c = params['D_c'], k_cond = params['k_cond'], L_c = params['L_c'], 
        N_c = params['N_c'], N_p = params['N_p'], R_p = params['R_p'], t_2 = params['t_2'], t_3 = params['t_3'], type_channel = params['type_channel'],
        
        Flow_Type = 'CounterFlow', H_DP_ON = True, C_DP_ON = True, n_disc = 50) # 27

if test_case == "TCO2_recup":

    HX.set_inputs(
        # First fluid
        fluid_H = 'CO2',
        T_su_H = 321.88, # K
        P_su_H = 5742510, # Pa
        m_dot_H = 415.93, # kg/s
    
        # Second fluid
        fluid_C = 'CO2',
        T_su_C = 305.61, # K
        P_su_C = 14153425, # Pa
        m_dot_C = 415.93, # kg/s  # Make sure to include fluid information
    )
    
    "Geometry Loading"
    params = {'alpha': 32.62, # Channel zigzag angle
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
              'type_channel' : 'Zigzag'} 

    Corr_H = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
    Corr_C = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
    
    # H_DP = "Gnielinski_DP"
    # C_DP = "Gnielinski_DP"    
    
    H_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach"}
    C_DP = {"SC" : "Darcy_Weisbach", "1P" : "Darcy_Weisbach"}  
    
    # ---------------------------------------------------------------------------------------------------------
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

# HX.set_DP() # equivalent to HX.set_DP(DP_type = None)
# HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
# HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

# "Solve the component"
HX.solve()
# HX.plot_cells()

