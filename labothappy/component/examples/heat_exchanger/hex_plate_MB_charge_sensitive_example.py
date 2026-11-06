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

# # from __future__ import division, print_function
# import __init__

import numpy as np

from labothappy.toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive_bis import HexMBChargeSensitive

#%%

# import time
# start_time = time.time()   

# "------------ Plate HX -------------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HexMBChargeSensitive('BPHE')

case_study = "C5_EVAP"

# BPHE correlation names (Corr_H / Corr_C): '1P' = single phase, '2P' = condensation on the
# hot side and evaporation on the cold side; None or 'auto' = automatic selection.
#   single phase : Sinnott, Kumar, Thonon, Muley_Manglik, Martin
#   condensation : Shah, Han, Longo
#   evaporation  : Cooper, Amalfi

# "Setting inputs"
if case_study == "C5_COND":
    # ---------------------------------------------------------------------------------------------------------
    # Condenser Case
    HX.set_inputs(
        # First fluid
        fluid_H = 'Cyclopentane',
        T_su_H = 139 + 273.15, # K
        P_su_H = 0.77*1e5, # Pa
        m_dot_H = 0.014, # kg/s
    
        # Second fluid
        fluid_C = 'Water',
        T_su_C = 12 + 273.15, # K
        P_su_C = 5*1e5, # Pa
        m_dot_C = 0.2, # kg/s  # Make sure to include fluid information
    )
    
    "Geometry Loading"
    
    HX_geom = PlateGeomSWEP()
    HX_geom.set_parameters("B20Hx24/1P") 
    
    Corr_H = {"1P" : "Martin", "2P" : "Han"}
    Corr_C = {"1P" : "Martin", "2P" : "Longo"} 
    
    HX.set_parameters(
        A_c = 0.616, 
        A_h = 0.616, 
        h = 0.041, 
        l = 0.324, 
        l_v = 0.2682, 
        
        n_channel_c = 12, 
        n_channel_h = 11, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 10*np.pi/180, # !!! arbitrary, 
        fooling = 0.705/100, 
        n_plates = 24, 
        plate_cond = 45, 
        w = 0.094, 
        
        Flow_Type = 'CounterFlow', n_disc = 50)

    HX.set_DP(DP_type="User-Defined", UD_H_DP = 6*1e3, UD_C_DP = 11*1e3)

if case_study == "C5_EVAP":
    # ---------------------------------------------------------------------------------------------------------
    # Evaporator Case
    HX.set_inputs(
        # First fluid
        fluid_H = 'INCOMP::T66',
        T_su_H = 243 + 273.15, # K
        P_su_H = 5*1e5, # Pa
        m_dot_H = 0.4, # kg/s
    
        # Second fluid
        fluid_C = 'Cyclopentane',
        T_su_C = 41 + 273.15, # K
        P_su_C = 31.5*1e5, # Pa
        m_dot_C = 0.014, # kg/s  # Make sure to include fluid information
    )
    
    "Geometry Loading"
    
    HX_geom = PlateGeomSWEP()
    HX_geom.set_parameters("B20Hx24/1P")
    
    Corr_H = {"1P" : "auto", "2P" : "Han"}
    Corr_C = {"1P" : "auto", "2P" : "Cooper"}

    HX.set_parameters(
        A_c = 0.616, 
        A_h = 0.616, 
        h = 0.041, 
        l = 0.324, 
        l_v = 0.2682, 
        
        n_channel_c = 12, 
        n_channel_h = 11, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 30*np.pi/180, # !!! arbitrary, 
        fooling = 0.705/100, 
        n_plates = 24, 
        plate_cond = 45, 
        w = 0.094, 
        
        Flow_Type = 'CounterFlow', n_disc = 50)

    HX.set_DP(DP_type="User-Defined", UD_H_DP = 6*1e3, UD_C_DP = 11*1e3)

if case_study == "C5_DSH":
    # ---------------------------------------------------------------------------------------------------------
    # Desuperheater Case
    HX.set_inputs(
        # First fluid
        fluid_H = 'Cyclopentane',
        T_su_H = 205 + 273.15, # K
        P_su_H = 1*1e5, # Pa
        m_dot_H = 0.014, # kg/s
    
        # Second fluid
        fluid_C = 'Water',
        T_su_C = 12 + 273.15, # K
        P_su_C = 5*1e5, # Pa
        m_dot_C = 0.08, # kg/s  # Make sure to include fluid information
    )
    
    "Geometry Loading"
    HX_geom = PlateGeomSWEP()
    HX_geom.set_parameters("B35TM0x10/1P") 
    
    Corr_H = {"1P" : "Martin", "2P" : "Longo"}
    Corr_C = {"1P" : "Martin", "2P" : "Cooper"}
    
    HX.set_parameters(
        A_c = 0.752, 
        A_h = 0.752, 
        h = 0.0446, 
        l = 0.393, 
        l_v = 0.324, # 5
        
        n_channel_c = 5, 
        n_channel_h = 4, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 10*np.pi/180, # !!! arbitrary, 
        fooling = 98.51/100, 
        n_plates = 10, 
        plate_cond = 45, 
        w = 0.243, # 23
        
        Flow_Type = 'CounterFlow', n_disc = 50)

    HX.set_DP(DP_type="User-Defined", UD_H_DP = 2*1e3, UD_C_DP = 0*1e3)

if case_study == "R1233zd(E)_COND":

    # ---------------------------------------------------------------------------------------------------------
    # Condenser Case
    HX.set_inputs(
        # First fluid
        fluid_H = 'R1233zd(E)',
        T_su_H = 334.3838266524477, # K
        P_su_H = 2*1e5, # Pa
        m_dot_H = 0.41973076344816307, # kg/s
    
        # Second fluid
        fluid_C = 'Water',
        T_su_C = 15 + 273.15, # K
        P_su_C = 2*1e5, # Pa
        m_dot_C = 2.6, # kg/s  # Make sure to include fluid information
    )

    "Geometry Loading"
    HX_geom = PlateGeomSWEP()
    HX_geom.set_parameters("P200THx140/1P_Condenser")
    
    Corr_H = {"1P" : "Martin", "2P" : "Longo"}
    Corr_C = {"1P" : "Martin", "2P" : "Cooper"}

    "Parameters Setting"

    HX.set_parameters(
        A_c = 17.8, 
        A_h = 17.8, 
        h = 0.346, 
        l = 0.525, 
        l_v = 0.450, # 5
        
        n_channel_c = 70, 
        n_channel_h = 69, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 30*np.pi/180, # !!! arbitrary, 
        fooling = 0.258/100, 
        n_plates = 140, 
        plate_cond = 45, 
        w = 0.243, # 23
        
        Flow_Type = 'CounterFlow', n_disc = 50)

    HX.set_DP(DP_type="User-Defined", UD_H_DP = 20*1e3, UD_C_DP = 50*1e3)

# ---------------------------------------------------------------------------------------------------------

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

# Set the heat transfer coefficients correlations
HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 28

# HX.set_DP() # equivalent to HX.set_DP(DP_type = None)
# HX.set_DP(DP_type="User-Defined", UD_C_DP = 10000, UD_H_DP = 10000) # Fixed User-Defined values, equally distributed over discretizations
# HX.set_DP(DP_type="Correlation_Global", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)
# HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

"Solve the component"
HX.solve()
HX.plot_cells()