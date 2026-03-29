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
from component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive

# from toolbox.geometries.heat_exchanger.geometry_shell_and_tube_hx import ShellAndTubeGeom
from CoolProp.CoolProp import PropsSI
from toolbox.geometries.heat_exchanger.c_geometry_HXs_Zorlu import Zorlu_HXs

import numpy as np

HX_name = 'PREHEATER'

n_disc =30
    
if HX_name == "PREHEATER":
    HX = HexMBChargeSensitive('Shell&Tube')
    
    
    HX.set_inputs(
        fluid_H = 'Water',
        T_su_H = 104.2 + 273.15, # K
        P_su_H = 2.1*1e5, # Pa
        m_dot_H = 171.6, # kg/s
    
        # Second fluid
        fluid_C = 'Cyclopentane',
        T_su_C = 50.1 + 273.15, # K
        P_su_C = 7.75*1e5, # Pa
        m_dot_C = 62.74, # kg/s  
        )
    
    
    geom_obj = Zorlu_HXs()
    geom_obj.set_parameters("ORC_preheater")
    HX.set_parameters(**geom_obj.geom)
    HX.set_parameters(n_disc=30)
    
    # "Correlation Loading"
    
    Corr_C = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
    Corr_H = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    # Corr_H_DP = {"1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"}
    # Corr_C_DP = {"1P" : "Gnielinski_DP", "2P" : "Choi_DP"}
    # -------------------------------------------------------------------------------------------------------------
     
    # HEAT TRANSFER COEFFICIENT SETTING
    
    HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31
    
    # PRESSURE DROP SETTING
    
    geom_obj.set_DP("ORC_preheater")
    
    HX.set_DP(**geom_obj.DP) # equivalent to HX.set_DP(DP_type = None)
    # HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
    # HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
    # HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
    
    "Solve the component"
    
    HX.solve()  # the function you want to profile
    HX.print_states_connectors()
    
    
    epsilon_pre = HX.Q_dot / (HX.su_C.m_dot * HX.su_C.cp * (HX.su_H.T - HX.su_C.T))
    print(f'  - epsilon_pre = {epsilon_pre}')

elif HX_name == "RECUPERATOR":
    HX = HexMBChargeSensitive('Tube&Fins')
    
    HX.set_inputs(
        fluid_H = 'Cyclopentane',
        T_su_H = 64.4 + 273.15, # K
        P_su_H = 0.8*1e5, # Pa
        m_dot_H = 62.71, # kg/s
    
        # Second fluid
        fluid_C = 'Cyclopentane',
        T_su_C = 31.27 + 273.15, # K
        P_su_C = 8.258*1e5, # Pa
        m_dot_C = 62.74, # kg/s  
        )
    "Set HTC"
    Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}
    Corr_H = {"1P" : "Tube_And_Fins", "2P" : "Tube_And_Fins"}
    # Corr_H = {"1P" : "Tube&Fins", "2P" : "Tube&Fins"}
    HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
    
    
    "Set geometry"
    geom_obj = Zorlu_HXs()
    geom_obj.set_parameters("ORC_recuperator")
    HX.set_parameters(**geom_obj.geom)
    HX.set_parameters(n_disc=n_disc)
    
    "Set pressure drops"
    # geom_obj.set_DP("ORC_recuperator")
    # HX.set_DP(**geom_obj.DP) 
    HX.set_DP() # no pressure drops normally
        
    "Solve the component"
    
    HX.solve() 
    HX.print_states_connectors()
    
    epsilon_rec = HX.Q_dot / (HX.su_C.m_dot * HX.su_C.cp * (HX.su_H.T - HX.su_C.T))
    print(f'  - epsilon_rec = {epsilon_rec}')
    # "Verification of epsilon"
    # C_dot_C = HX.su_C.m_dot * HX.su_C.cp
    # C_dot_H = HX.su_H.m_dot * HX.su_H.cp
    # C_dot_min = min(C_dot_C , C_dot_H)
    # epsilon_rec = HX.Q_dot / (C_dot_min  * (HX.su_H.T - HX.su_C.T))
    # print(f'  - epsilon = {epsilon_rec}')
    
elif HX_name == "ACC":
    HX = HexMBChargeSensitive('Tube&Fins')
    
    HX.set_inputs(
        fluid_H = 'Cyclopentane',
        T_su_H = 64.4 + 273.15, # K
        P_su_H = 0.8*1e5, # Pa
        m_dot_H = 62.71, # kg/s
    
        # Second fluid
        fluid_C = 'Air',
        T_su_C = 31.27 + 273.15, # K
        P_su_C = 8.258*1e5, # Pa
        m_dot_C = 62.74, # kg/s  
        )
    "Set HTC"
    Corr_C = {"1P" : "Gnielinski", "2P" : "Boiling_curve"}
    Corr_H = {"1P" : "Tube_And_Fins", "2P" : "Tube_And_Fins"}
    # Corr_H = {"1P" : "Tube&Fins", "2P" : "Tube&Fins"}
    HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
    
    
    "Set geometry"
    geom_obj = Zorlu_HXs()
    geom_obj.set_parameters("ORC_ACC")
    HX.set_parameters(**geom_obj.geom)
    HX.set_parameters(n_disc=n_disc)
    
    "Set pressure drops"
    # geom_obj.set_DP("ORC_recuperator")
    # HX.set_DP(**geom_obj.DP) 
    HX.set_DP() # no pressure drops normally
        
    "Solve the component"
    
    HX.solve() 
    HX.print_states_connectors()
    
    epsilon_rec = HX.Q_dot / (HX.su_C.m_dot * HX.su_C.cp * (HX.su_H.T - HX.su_C.T))
    print(f'  - epsilon_ACC = {epsilon_ACC}')
    