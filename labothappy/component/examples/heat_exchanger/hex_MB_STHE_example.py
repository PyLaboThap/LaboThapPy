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

from labothappy.component.heat_exchanger.hex_MB import HexMB
from CoolProp.CoolProp import PropsSI    


#%%

import time
start_time = time.time()

"------------ Shell and Tube HTX ------------------------------------------------------------------------------------------"

"HTX Instanciation"

case_study = 'Methanol'
n_disc = 30
        
# # -------------------------------------------------------------------------------------------------------------
# # Methanol Sensible HT Case

HX = HexMB('STHE')

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
              'R_f_T': 0.0002,
              'R_f_S': 0.00033,
              'k_wall': 50,
              'shell_side': 'H',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'D_s': 0.83,
              'L_t': 3.379,
              'D_t_o': 0.016,
              't_t': 0.016*0.1,
              'B': 0.5,
              'n_passes': 2,
              'n_t': 1567,
              'PR': 1.25,
              'theta': 60,
              'BC': 25}
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
    
        charge_compute = False,

        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    Corr_C = {"1P" : "Gnielinski", "2P" : "VDI"}
    Corr_H = {"1P" : "Kern", "2P" : "Kern"}
    
    Corr_H_DP = {"1P" : "Kern", "2P" : "Kern"}
    Corr_C_DP = {"1P" : "Churchill", "2P" : "MSH"}
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
              'R_f_T': 0.000176,
              'R_f_S': 0.000176,
              'k_wall': 50,
              'shell_side': 'H',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'D_s': 0.173,
              'L_t': 4.18,
              'D_t_o': 0.01,
              't_t': 0.001,
              'B': 0.381,
              'n_passes': 2,
              'n_t': 230,
              'PR': 1.25,
              'theta': 60,
              'BC': 25}        
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
        
        chare_compute = False,
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_C = {"1P" : "Gnielinski", "2P" : "VDI"}
    Corr_H = {"1P" : "Kern", "2P" : "Kern"}
    
    Corr_H_DP = {"1P" : "Kern", "2P" : "Kern"}
    Corr_C_DP = {"1P" : "Churchill", "2P" : "MSH"}

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
              'R_f_T': 0.0,
              'R_f_S': 0.0,
              'k_wall': 50,
              'shell_side': 'C',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'D_s': 0.4,
              'L_t': 5.348,
              'D_t_o': 0.01958,
              't_t': 0.002,
              'B': 0.59,
              'n_passes': 1,
              'n_t': 194,
              'PR': 1.25,
              'theta': 60,
              'BC': 20}
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
    
        charge_compute = False,
    
        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"1P" : "Gnielinski", "2P" : "Thome_Condensation"}
    # Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling_gungor_winterton"}
    Corr_C = {"1P" : "Kern", "2P" : "Kern"}
    
    Corr_C_DP = {"1P" : "Churchill", "2P" : "Churchill"}
    Corr_H_DP = {"1P" : "Churchill", "2P" : "MSH"}

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
              'R_f_T': 0.000176,
              'R_f_S': 0.000176,
              'k_wall': 50,
              'shell_side': 'H',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'D_s': 0.41353,
              'L_t': 5.48864,
              'D_t_o': 0.02128,
              't_t': 0.002128,
              'B': 0.94,
              'n_passes': 2,
              'n_t': 140,
              'PR': 1.25,
              'theta': 0,
              'BC': 25}
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
    
        charge_compute = False,

        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_C = {"1P" : "Gnielinski", "2P" : "Choi"}
    Corr_H = {"1P" : "Kern", "2P" : "Kern"}
    
    Corr_H_DP = {"1P" : "Kern", "2P" : "Kern"}
    Corr_C_DP = {"1P" : "Konakov", "2P" : "MSH"}

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
            'R_f_T': 0,
            'R_f_S': 0,
            'k_wall': 50,
            'shell_side': 'C',
            'Flow_Type': 'Shell&Tube',
            'n_disc': 30,
            'D_s': 1.43,
            'L_t': 9.21,
            'D_t_o': 0.02,
            't_t': 0.003,
            'B': 0.5*1.43,
            'n_passes': 1,
            'n_t': 1480,
            'PR': 1.5,
            'theta': 45,
            'BC': 25}
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
    
        charge_compute = False,

        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
    Corr_C = {"1P" : "Kern", "2P" : "Kern"}
    
    Corr_C_DP = {"1P" : "Kern", "2P" : "Kern"}
    Corr_H_DP = {"1P" : "Churchill", "2P" : "MSH", "SC" : "Cheng-CO2"} 

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
              'BC': 26.726000000000003,
              'D_s': 1.6764,
              'L_t': 14.84,
              'D_t_o': 0.009524999999999999,
              'n_passes': 1,
              't_t': 0.0005587999999999999,
              'B': 0.742,
              'R_f_S': 0,
              'R_f_T': 0,
              'n_series': 1,
              'n_parallel': 1,
              'n_t': 15638,
              'PR': 1.33,
              'k_wall': 50,
              'theta': 60,
              'shell_side': 'C',
              'Flow_Type': 'Shell&Tube',
              'n_disc': 1}
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
    
        charge_compute = False,

        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
    Corr_C = {"1P" : "Kern", "2P" : "Kern"}
    
    Corr_C_DP = {"1P" : "Kern", "2P" : "Kern"}
    Corr_H_DP = {"1P" : "Konakov", "2P" : "MSH", "SC" : "Cheng-CO2"} 

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
              'R_f_T': 0,
              'R_f_S': 0,
              'k_wall': 50,
              'shell_side': 'C',
              'Flow_Type': 'Shell&Tube',
              'n_disc': int(n_disc),
              'D_s': 1.43,
              'L_t': 9.21,
              'D_t_o': 0.02,
              't_t': 0.003,
              'B': 1.43,
              'n_passes': 2,
              'n_t': 1480,
              'PR': 1.5,
              'theta': 45,
              'BC': 25}
    
    HX.set_parameters(
        BC = params['BC'], D_s = params['D_s'], L_t = params['L_t'], 
        D_t_o = params['D_t_o'], n_passes = params['n_passes'], t_t = params['t_t'],
        B = params['B'], R_f_S = params['R_f_S'],
        R_f_T = params['R_f_T'], n_series = params['n_series'], n_parallel = params['n_parallel'], n_t = params['n_t'], 
        PR = params['PR'], k_wall = params['k_wall'], theta = params['theta'],
    
        shell_side = params['shell_side'],
    
        charge_compute = False,

        Flow_Type = params['Flow_Type'], n_disc = params['n_disc'])
    
    # "Correlation Loading"
    
    Corr_H = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Thome_Condensation"}
    Corr_C = {"SC" : "Kern", "1P" : "Kern", "2P" : "Kern"}
    
    Corr_H_DP = {"SC" : "Konakov", "1P" : "Konakov", "2P" : "MSH"}
    Corr_C_DP = {"SC" : "Kern", "1P" : "Kern", "2P" : "Kern"}
        
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
