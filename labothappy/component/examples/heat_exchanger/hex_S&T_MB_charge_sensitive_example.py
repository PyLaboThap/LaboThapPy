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

from labothappy.toolbox.plots.plot_MB_STHE import plot_MB_STHE

from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from CoolProp.CoolProp import PropsSI    


#%%

import time
start_time = time.time()

"------------ Shell and Tube HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

case_study = 'R134a'
n_disc = 30
        
# # -------------------------------------------------------------------------------------------------------------
# # Methanol Sensible HT Case

HX = HexMBChargeSensitive('Shell&Tube')

if case_study == 'Methanol':

    HX.set_inputs(
        # First fluid
        fluid_H = 'Methanol',
        T_su_H = 95 + 273.15, # K
        P_su_H = 10e5, # Pa
        m_dot_H = 27.8, # kg/s
    
        # Second fluid
        fluid_C = 'Water',
        T_su_C = 25 + 273.15, # K
        P_su_C = 5*1e5, # Pa
        m_dot_C = 68.9, # kg/s  # Make sure to include fluid information
    )
    
    "Geometry Loading"
    params = {'n_series': 1,
              'n_parallel': 1,
              'foul_t': 0.0002,
              'foul_s': 0.00033,
              'tube_cond': 50,
              'Shell_Side': 'H',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'Shell_ID': 0.83,
              'Tube_L': 3.379,
              'Tube_OD': 0.016,
              'Tube_t': 0.016*0.1,
              'central_spacing': 0.5,
              'Tube_pass': 2,
              'n_tubes': 1567,
              'pitch_ratio': 1.25,
              'tube_layout': 60,
              'Baffle_cut': 25}
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
    Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_H_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    Corr_C_DP = {"1P" : "Gnielinski_DP", "2P" : "Choi_DP"}
# -------------------------------------------------------------------------------------------------------------

# Sizing code example case
# R134a evaporator case

if case_study == 'R134a':
    
    HX.set_inputs(
                  # Hot Fluid
                  T_su_H = 273.15 + 26, # K
                  P_su_H = 2*1e5, # 51.75*1e3, # Pa
                  m_dot_H = 5.35, # kg/s
                  fluid_H = 'Water',
                  
                  # Cold Fluid
                  h_su_C = PropsSI('H','T', 273.15+7,'Q',0,'R134a')+1, # kJ/kg
                  P_su_C = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # 51.75*1e3, # Pa
                  m_dot_C = 1.62, # kg/s
                  fluid_C = 'R134a'
                  )
    
    "Geometry Loading"

    params = {'n_series': 1,
              'n_parallel': 1,
              'foul_t': 0.000176,
              'foul_s': 0.000176,
              'tube_cond': 50,
              'Shell_Side': 'H',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'Shell_ID': 0.173,
              'Tube_L': 4.18,
              'Tube_OD': 0.01,
              'Tube_t': 0.001,
              'central_spacing': 0.381,
              'Tube_pass': 2,
              'n_tubes': 230,
              'pitch_ratio': 1.25,
              'tube_layout': 60,
              'Baffle_cut': 25}        
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
    # Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling_gungor_winterton"}
    Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_H_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    Corr_C_DP = {"1P" : "Gnielinski_DP", "2P" : "Muller_Steinhagen_Heck_DP"}

if case_study == 'R134a_cond':
            
    HX.set_inputs(
                  # Hot Fluid
                  T_su_C = 273.15 + 12, # K
                  P_su_C = 2*1e5, # 51.75*1e3, # Pa
                  m_dot_C = 20, # kg/s
                  fluid_C = 'Water',
                  
                  # Cold Fluid
                  h_su_H = PropsSI('H','T', 273.15+30,'Q',1,'R134a')-1, # kJ/kg
                  P_su_H = PropsSI('P','T', 273.15+30,'Q',1,'R134a'), # 51.75*1e3, # Pa
                  m_dot_H = 5.78, # kg/s
                  fluid_H = 'R134a'
                  )
    
    "Geometry Loading"
    
    params = {'n_series': 1,
              'n_parallel': 1,
              'foul_t': 0.0,
              'foul_s': 0.0,
              'tube_cond': 50,
              'Shell_Side': 'C',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'Shell_ID': 0.4,
              'Tube_L': 5.348,
              'Tube_OD': 0.01958,
              'Tube_t': 0.002,
              'central_spacing': 0.59,
              'Tube_pass': 1,
              'n_tubes': 194,
              'pitch_ratio': 1.25,
              'tube_layout': 60,
              'Baffle_cut': 20}
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"1P" : "Gnielinski", "2P" : "Thome_Condensation"}
    # Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling_gungor_winterton"}
    Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_C_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    Corr_H_DP = {"1P" : "Gnielinski_DP", "2P" : "Choi_DP"}

if case_study == 'R32':
    
    HX.set_inputs(
                  # Hot Fluid
                  T_su_H = 273.15 + 30, # K
                  P_su_H = 2*1e5, # 51.75*1e3, # Pa
                  m_dot_H = 36.3, # kg/s
                  fluid_H = 'Water',
                  
                  # Cold Fluid
                  h_su_C = PropsSI('H','T', 273.15+6,'Q',0,'R32')+1, # kJ/kg
                  P_su_C = PropsSI('P','T', 273.15+6,'Q',0,'R32'), # 51.75*1e3, # Pa
                  m_dot_C = 3.2716, # kg/s
                  fluid_C = 'R32'
                  )
    
    "Geometry Loading"
    
    params = {'n_series': 1,
              'n_parallel': 1,
              'foul_t': 0.000176,
              'foul_s': 0.000176,
              'tube_cond': 50,
              'Shell_Side': 'H',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'Shell_ID': 0.41353,
              'Tube_L': 5.48864,
              'Tube_OD': 0.02128,
              'Tube_t': 0.002128,
              'central_spacing': 0.94,
              'Tube_pass': 2,
              'n_tubes': 140,
              'pitch_ratio': 1.25,
              'tube_layout': 0,
              'Baffle_cut': 25}
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_C = {"1P" : "Gnielinski", "2P" : "choi_boiling"}
    Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_H_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    Corr_C_DP = {"1P" : "Gnielinski_DP", "2P" : "Muller_Steinhagen_Heck_DP"}

# -------------------------------------------------------------------------------------------------------------

if case_study == 'CO2_GasCooler_1':

    # Sizing code example case
    # CO2 GasCooler Case
    
    HX.set_inputs(
                  # Hot Fluid
                  T_su_H = 273.15 + 100, # K
                  P_su_H = 100*1e5, # Pa
                  m_dot_H = 100, # kg/s
                  fluid_H = 'CO2',
                  
                  # Cold Fluid
                  T_su_C = 273.15 + 20, # K
                  P_su_C = 5*1e5, # Pa
                  m_dot_C = 47.82, # kg/s
                  fluid_C = 'Water'
                  )
    
    "Geometry Loading"
    
    params = {'n_series': 1,
              'n_parallel': 1,
            'foul_t': 0,
            'foul_s': 0,
            'tube_cond': 50,
            'Shell_Side': 'C',
            'Flow_Type': 'Shell&Tube',
            'n_disc': 30,
            'Shell_ID': 1.43,
            'Tube_L': 9.21,
            'Tube_OD': 0.02,
            'Tube_t': 0.003,
            'central_spacing': 0.5*1.43,
            'Tube_pass': 1,
            'n_tubes': 1480,
            'pitch_ratio': 1.5,
            'tube_layout': 45,
            'Baffle_cut': 25}
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
    Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_C_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    Corr_H_DP = {"1P" : "Gnielinski_DP", "2P" : "Muller_Steinhagen_Heck_DP", "SC" : "Cheng_CO2_DP"} 

# -------------------------------------------------------------------------------------------------------------

if case_study == 'CO2_GasCooler_2':

    HX.set_inputs(
                  # Hot Fluid
                  T_su_H = 273.15 + 105.4, # K
                  P_su_H = 90*1e5, # 51.75*1e3, # Pa
                  m_dot_H = 1590.82/4, # kg/s
                  fluid_H = 'CO2',
                  
                  # Cold Fluid
                  T_su_C = 273.15 + 25, # K
                  P_su_C = 5*1e5, # 51.75*1e3, # Pa
                  m_dot_C = 980/4, # kg/s
                  fluid_C = 'Water'
                  )
    
    params = {'htc_type': 'Correlation',
              'Baffle_cut': 26.726000000000003,
              'Shell_ID': 1.6764,
              'Tube_L': 14.84,
              'Tube_OD': 0.009524999999999999,
              'Tube_pass': 1,
              'Tube_t': 0.0005587999999999999,
              'central_spacing': 0.742,
              'foul_s': 0,
              'foul_t': 0,
              'n_series': 1,
              'n_parallel': 1,
              'n_tubes': 15638,
              'pitch_ratio': 1.33,
              'tube_cond': 50,
              'tube_layout': 60,
              'Shell_Side': 'C',
              'Flow_Type': 'Shell&Tube',
              'n_disc': 1}
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
    Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_C_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    Corr_H_DP = {"1P" : "Gnielinski_DP", "2P" : "Muller_Steinhagen_Heck_DP", "SC" : "Cheng_CO2_DP"} 

# # -------------------------------------------------------------------------------------------------------------

if case_study == 'CO2_CD':

    # Sizing code example case
    # CO2 Condenser Case
    
    HX.set_inputs(
                  # Hot Fluid
                  T_su_H = 273.15 + 20, # K
                  P_su_H = 5087147.357957976, # 51.75*1e3, # Pa
                  m_dot_H = 30, # kg/s
                  fluid_H = 'CO2',
                  
                  # Cold Fluid
                  T_su_C = 273.15 + 3, # K
                  P_su_C = 5*1e5, # 51.75*1e3, # Pa
                  m_dot_C = 100, # kg/s
                  fluid_C = 'Water'
                  )
    
    "Geometry Loading"
    
    params = {'n_series': 1,
              'n_parallel' : 2,
              'foul_t': 0,
              'foul_s': 0,
              'tube_cond': 50,
              'Shell_Side': 'C',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'Shell_ID': 1.43,
              'Tube_L': 9.21,
              'Tube_OD': 0.02,
              'Tube_t': 0.003,
              'central_spacing': 1.43,
              'Tube_pass': 2,
              'n_tubes': 1480,
              'pitch_ratio': 1.5,
              'tube_layout': 45,
              'Baffle_cut': 25}
    
    HX.set_parameters(
        Baffle_cut = params['Baffle_cut'], Shell_ID = params['Shell_ID'], Tube_L = params['Tube_L'], 
        Tube_OD = params['Tube_OD'], Tube_pass = params['Tube_pass'], Tube_t = params['Tube_t'],
        central_spacing = params['central_spacing'], foul_s = params['foul_s'],
        foul_t = params['foul_t'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_tubes = params['n_tubes'], 
        pitch_ratio = params['pitch_ratio'], tube_cond = params['tube_cond'], tube_layout = params['tube_layout'],
    
        Shell_Side = params['Shell_Side'],
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Thome_Condensation"}
    Corr_C = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    Corr_H_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Choi_DP"}
    Corr_C_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
        
# # -------------------------------------------------------------------------------------------------------------

# "User-Defined Heat Transfer Coefficients Setting"

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

# # -------------------------------------------------------------------------------------------------------------

# HEAT TRANSFER COEFFICIENT SETTING

HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

# PRESSURE DROP SETTING

# HX.set_DP() # equivalent to HX.set_DP(DP_type = None)
# HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
# HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

"Solve the component"

HX.solve()  # the function you want to profile
HX.plot_cells()

plot_MB_STHE(HX, props='htc')
