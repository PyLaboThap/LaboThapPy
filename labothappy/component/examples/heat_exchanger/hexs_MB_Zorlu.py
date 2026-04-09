# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 13:06:16 2026

@author: titouanjanod
"""

# from __future__ import division, print_function

import __init__
from component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from component.heat_exchanger.hex_eNTU import HexeNTU

# from toolbox.geometries.heat_exchanger.geometry_shell_and_tube_hx import ShellAndTubeGeom
from CoolProp.CoolProp import PropsSI
from toolbox.geometries.heat_exchanger.c_geometry_HXs_Zorlu import Zorlu_HXs

import numpy as np

HX_name = 'ACC'

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
    HX.set_DP()

    #geom_obj.set_DP("ORC_preheater")
    
    #HX.set_DP(**geom_obj.DP)
    # HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
    # HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
    # HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
    
    "Solve the component"
    
    HX.solve()  # the function you want to profile
    HX.print_states_connectors()
    print(f"  - epsilon = {HX.epsilon_th}")
    print(f'  - Q_dot = {HX.Q_dot/1000} kW')
    

elif HX_name == "RECUPERATOR":
    Type_HX = 'epsNTU'   # 'epsNTU'  or 'MB'
    
    if Type_HX == 'MB':
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
        HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
        
        "Set geometry"
        geom_obj = Zorlu_HXs()
        geom_obj.set_parameters("ORC_recuperator")
        HX.set_parameters(**geom_obj.geom)
        HX.set_parameters(n_disc=n_disc)
        
        "Set pressure drops"
        # geom_obj.set_DP("ORC_recuperator")
        # HX.set_DP(**geom_obj.DP) 
        HX.set_DP() # no pressure drops 
            
        "Solve the component"
        
        HX.solve() 
        HX.print_states_connectors()
        print(f'  - epsilon = {HX.epsilon_th}')
        print(f'  - Q_dot = {HX.Q_dot/1000} kW')
    
    elif Type_HX =='epsNTU' :
        HX = HexeNTU("Tube&Fins")
        
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
        Corr_C = "Gnielinski"
        Corr_H = "Tube_And_Fins"
        HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
        
        
        
        geom_obj = Zorlu_HXs()
        geom_obj.set_parameters("ORC_recuperator")
        HX.set_parameters(**geom_obj.geom)
        geom_obj.set_DP("ORC_recuperator")
        HX.set_DP(**geom_obj.DP)
        
        
        HX.solve()
        HX.print_results()
        print(f"\n  - Q_dot = {HX.Q_hex.Q_dot/1000} kW")
        # HX.print_setup()
        
        
    else:
        print(f"Type_HX shall be equal to 'epsNTU' or 'Tube&Fins'")

    
    
    
    
elif HX_name == "ACC":
    
    Type_HX = 'epsNTU'   # 'epsNTU'  or 'MB'
    
    if Type_HX == 'MB':
        HX = HexMBChargeSensitive('Tube&Fins')
    
        
        HX.set_inputs(
            fluid_H = 'Cyclopentane',
            T_su_H = 38.46 + 273.15, # K
            P_su_H = 0.7003*1e5, # Pa
            m_dot_H = 62.71, # kg/s
        
            # Second fluid
            fluid_C = 'Air',
            T_su_C = 15 + 273.15, # K
            P_su_C = 1.01325*1e5, # Pa
            m_dot_C = 2557, # kg/s  
            )
        "Set HTC"
        Corr_H = {"1P" : "Gnielinski", "2P" : "Thome_Condensation"}
        Corr_C = {"1P" : "Tube_And_Fins", "2P" : "Tube_And_Fins"}
        # Corr_H = {"1P" : "Tube&Fins", "2P" : "Tube&Fins"}
        HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
        
        
        "Set geometry"
        geom_obj = Zorlu_HXs()
        geom_obj.set_parameters("ORC_ACC")
        HX.set_parameters(**geom_obj.geom)
        HX.set_parameters(n_disc=n_disc)
        
        "Set pressure drops"
        # geom_obj.set_DP("ORC_recuperator")
        HX.set_DP(**geom_obj.DP) 
        HX.set_DP() # no pressure drops normally
            
        "Solve the component"
        
        HX.solve() 
        HX.print_states_connectors()
        print(f'  - epsilon = {HX.epsilon_th}')
        
        
        
    elif Type_HX =='epsNTU' :
        HX = HexeNTU("Tube&Fins")

        HX.set_inputs(
            fluid_H = 'Cyclopentane',
            T_su_H = 38.46 + 273.15, # K
            P_su_H = 0.7003*1e5, # Pa
            m_dot_H = 62.71, # kg/s
        
            # Second fluid
            fluid_C = 'Air',
            T_su_C = 15 + 273.15, # K
            P_su_C = 1.01325*1e5, # Pa
            m_dot_C = 2557, # kg/s  
            )
        
        Corr_C = "Gnielinski"
        Corr_H = "Tube_And_Fins"
        HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
        
        
        
        geom_obj = Zorlu_HXs()
        geom_obj.set_parameters("ORC_ACC")
        HX.set_parameters(**geom_obj.geom)
        geom_obj.set_DP("ORC_ACC")
        HX.set_DP(**geom_obj.DP)
        
        
        HX.solve()
        # HX.print_results()
        # print(f"\n  - Q_dot = {HX.Q_hex.Q_dot/1000} kW")
        
        
    
    

elif HX_name == "IHX":
    HX = HexMBChargeSensitive('Shell&Tube')
    
    

    
    HX.set_inputs(
        fluid_H = 'Cyclopentane',
        T_su_H = 139.8 + 273.15, # K
        P_su_H = 9.744*1e5, # Pa
        m_dot_H = 35.3, # kg/s
    
        # Second fluid
        fluid_C = 'Cyclopentane',
        T_su_C = 97.6 + 273.15, # K
        P_su_C = 3.84*1e5, # Pa
        m_dot_C = 35.3, # kg/s  
        )
    
    "Set HTC"
    Corr_H = {"1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}
    Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    # Corr_H = {"1P" : "Tube&Fins", "2P" : "Tube&Fins"}
    HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
    
    
    "Set geometry"
    geom_obj = Zorlu_HXs()
    geom_obj.set_parameters("HTHP_IHX")
    HX.set_parameters(**geom_obj.geom)
    HX.set_parameters(n_disc=n_disc)
    
    "Set pressure drops"
    # geom_obj.set_DP("ORC_recuperator")
    # HX.set_DP(**geom_obj.DP) 
    HX.set_DP() # if no pressure drops 
        
    "Solve the component"
    
    HX.solve() 
    HX.print_states_connectors()
    print(f'  - Q_dot = {HX.Q_dot/1000} kW')
    
    "Verification of epsilon"
    C_dot_C = HX.su_C.m_dot * HX.su_C.cp
    C_dot_H = HX.su_H.m_dot * HX.su_H.cp
    C_dot_min = min(C_dot_C , C_dot_H)
    epsilon_ihx = HX.Q_dot / (C_dot_min  * (HX.su_H.T - HX.su_C.T))
    print(f'  - epsilon_ihx = {epsilon_ihx}')
    
    
elif HX_name == "HTHP_EVAP":
    HX = HexMBChargeSensitive('Shell&Tube')
    
    
    T_ex_IHX = 115.9 + 273.15 #K
    p_ex_IHX = 924.4 *1e3 # Pa
    h_ex_IHX = PropsSI("H", "P", p_ex_IHX, "T", T_ex_IHX, "Cyclopentane")
    h_su_ev = h_ex_IHX
    
    HX.set_inputs(
        fluid_H = 'Water',
        T_su_H = 104.2 + 273.15, # K
        P_su_H = 1.89*1e5, # Pa
        m_dot_H = 857.9, # kg/s
    
        # Second fluid
        fluid_C = 'Cyclopentane',
        h_su_C = h_su_ev,
        # T_su_C = 97.42 + 273.15 -0.01, # K  -0.01 to make it like a liquid 
        P_su_C = 3.99*1e5, # Pa
        m_dot_C = 35.3, # kg/s  
        )
    
    "Set HTC"
    Corr_H = {"1P" : "Gnielinski", "2P" : "Gnielinski"}
    Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    # Corr_H = {"1P" : "Tube&Fins", "2P" : "Tube&Fins"}
    HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation'
    
    
    "Set geometry"
    geom_obj = Zorlu_HXs()
    geom_obj.set_parameters("HTHP_evap")
    HX.set_parameters(**geom_obj.geom)
    HX.set_parameters(n_disc=n_disc)
    
    "Set pressure drops"
    # geom_obj.set_DP("HTHP_evap")
    # HX.set_DP(**geom_obj.DP) 
    HX.set_DP() # if no pressure drops 
        
    "Solve the component"
    
    HX.solve() 
    HX.print_states_connectors()
    print(f'  - Q_dot = {HX.Q_dot/1000} kW')
    
    "Verification of epsilon"
    # C_dot_C = HX.su_C.m_dot * HX.su_C.cp
    # C_dot_H = HX.su_H.m_dot * HX.su_H.cp
    # C_dot_min = min(C_dot_C , C_dot_H)
    # epsilon_ihx = HX.Q_dot / (C_dot_min  * (HX.su_H.T - HX.su_C.T))
    # print(f'  - epsilon_ihx = {epsilon_ihx}')