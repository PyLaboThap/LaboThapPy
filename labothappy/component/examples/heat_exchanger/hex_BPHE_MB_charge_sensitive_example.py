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

import CoolProp.CoolProp as CP
from labothappy.toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive_bis import HexMBChargeSensitive

#%%

# import time
# start_time = time.time()   

# "------------ Plate HX -------------------------------------------------------------------------------------------------"

"HTX Instanciation"

HX = HexMBChargeSensitive('BPHE')

case_study = "C5_COND"

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
        
        n_channels_C = 12, 
        n_channels_H = 11, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 10*np.pi/180, # !!! arbitrary, 
        fooling = 0.705/100, 
        n_plates = 24, 
        plate_cond = 45, 
        w = 0.094, 
        
        charge_compute = False,

        Flow_Type = 'CounterFlow', n_disc = 50)

    # HX.set_DP(DP_type="User-Defined", UD_H_DP = 6*1e3, UD_C_DP = 11*1e3)

    Corr_H_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    Corr_C_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    
    
elif case_study == "C5_EVAP":
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
        
        n_channels_C = 12, 
        n_channels_H = 11, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 30*np.pi/180, # !!! arbitrary, 
        fooling = 0.705/100, 
        n_plates = 24, 
        plate_cond = 45, 
        w = 0.094, 
        
        charge_compute = False,

        Flow_Type = 'CounterFlow', n_disc = 50)

    # HX.set_DP(DP_type="User-Defined", UD_H_DP = 6*1e3, UD_C_DP = 11*1e3)

    Corr_H_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    Corr_C_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    
elif case_study == "C5_DSH":
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
        
        n_channels_C = 5, 
        n_channels_H = 4, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 10*np.pi/180, # !!! arbitrary, 
        fooling = 98.51/100, 
        n_plates = 10, 
        plate_cond = 45, 
        w = 0.243, # 23
        
        charge_compute = False,

        Flow_Type = 'CounterFlow', n_disc = 50)

    # HX.set_DP(DP_type="User-Defined", UD_H_DP = 2*1e3, UD_C_DP = 0*1e3)

    Corr_H_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    Corr_C_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    
elif case_study == "R1233zd(E)_COND":

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
        
        n_channels_C = 70, 
        n_channels_H = 69, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 30*np.pi/180, # !!! arbitrary, 
        fooling = 0.258/100, 
        n_plates = 140, 
        plate_cond = 45, 
        w = 0.243, # 23
        
        charge_compute = False,

        Flow_Type = 'CounterFlow', n_disc = 50)

    # HX.set_DP(DP_type="User-Defined", UD_H_DP = 20*1e3, UD_C_DP = 50*1e3)

    Corr_H_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    Corr_C_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    
    

elif case_study == 'Condenser_R1233ZDE':
    T_su_H=39.3+273.15
    P_su_H=CP.PropsSI('P', 'T', 22.1+273.15, 'Q', 1, 'R1233ZDE')
    h_su_H = CP.PropsSI('H', 'T', T_su_H, 'P', P_su_H, 'R1233ZDE')
    m_dot_H = 0.2

    T_su_C = 5+273.15
    h_su_C = CP.PropsSI('H', 'T', T_su_C, 'P', 1e5, 'Water')
    m_dot_C = 0.65

    HX.set_inputs(
        fluid_H='R1233ZDE', h_su_H=h_su_H, P_su_H=P_su_H, m_dot_H=m_dot_H,
        fluid_C='Water', h_su_C=h_su_C,  P_su_C=1e5, m_dot_C=m_dot_C
    )

    corrugation_amplitude = 0.002
    phi = 1.26                      # Surface enlargement factor (calibration parameter)
    c = 6 * phi - 1
    X = np.sqrt(((np.sqrt(8 * c**2 - 56) - c) / 7)**2 - 1)
    corrugation_pitch = np.pi * corrugation_amplitude / X

    HX.set_parameters(
        A_c = 17.8, 
        A_h = 17.8, 
        h = 0.346, 
        l = 0.45, 
        l_v = 0.450, # 5
        
        n_channels_C = 70, 
        n_channels_H = 69, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 60*np.pi/180, # !!! arbitrary, 
        fooling = 0, 
        n_plates = 140, 
        plate_cond = 45, 
        w = 0.243, # 23
        
        charge_compute = True,
        
        Flow_Type = 'CounterFlow', n_disc = 50)     

    # HX.set_parameters(flow_type='CounterFlow',          # Counter-flow: no LMTD correction (F = 1)
    #                         n_disc=3,                        # Number of regular cells (phase-change cells are added)
    #                         A_H=17.8, A_C=17.8,             # Heat transfer area [m^2]
    #                         V_C=16.63*1e-3, V_H=16.87*1e-3,       # Total volume of each side [m^3] (à vérifier quel côté est lequel!!)
    #                         n_plates = 140,
    #                         n_channels_H = 69,
    #                         n_channels_C = 70,
    #                         L_plate = 0.45,                  # Longueur H
    #                         W_plate = 0.1635,                # Longueur D
    #                         L_eff = 0.4485,                  # Longueur C
    #                         plate_pitch = 0.004,            
    #                         corrugation_amplitude = corrugation_amplitude,  
    #                         corrugation_pitch = corrugation_pitch,    #0.00573,    # Estimation based on enlargement factor
    #                         chevron_angle = 60*(np.pi/180), # Estimation of the chevron angle (no ideal actually + check nomenclature! could be very different depending of the correlation)
    #                         k_wall = 45,                    # ???
    #                         R_fouling = 0.0001,
    #                         roughness = 0)  
    
    Corr_H = {"1P" : "Martin", "2P" : "Longo"}
    Corr_C = {"1P" : "Martin", "2P" : "Cooper"}

    Corr_H_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    Corr_C_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    
elif case_study == 'Evaporator_R1233ZDE':


    # Models inputs
    HX.set_inputs(
        fluid_H='Water', T_su_H=70+273.15, P_su_H=1e5, m_dot_H=2,
        fluid_C='R1233ZDE', T_su_C=20+273.15,  P_su_C=CP.PropsSI('P', 'T', 55.4+273.15, 'Q', 1, 'R1233ZDE'), m_dot_C=0.5
    )

    HX.set_parameters(
        A_c = 17.8, 
        A_h = 17.8, 
        h = 0.346, 
        l = 0.525, 
        l_v = 0.450, # 5
        
        n_channels_C = 70, 
        n_channels_H = 69, 
        
        t_plates = 0.0005, # !!! arbitrary
        beta = 30*np.pi/180, # !!! arbitrary, 
        fooling = 0, 
        n_plates = 140, 
        plate_cond = 45, 
        w = 0.243, # 23
        
        charge_compute = True,
    
        Flow_Type = 'CounterFlow', n_disc = 50)

    Corr_H = {"1P" : "Martin", "2P" : "Longo"}
    Corr_C = {"1P" : "Martin", "2P" : "Cooper"}

    Corr_H_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    Corr_C_DP = {'1P': 'Martin', '2P':'Lockhart_Martinelli'}
    
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
HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

HX.set_void_fraction('Zivi')


"Solve the component"
HX.solve()

# print(f"Heat transfer rate      : {HX.Q_dot:10.1f} W")
# print(f"Hot exhaust temperature : {HX.ex_H.T - 273.15:10.2f} °C")
# print(f"Cold exhaust temperature: {HX.ex_C.T - 273.15:10.2f} °C")
# print(f"Hot exhaust pressure    : {HX.ex_H.p:10.0f} Pa")
# print(f"Cold exhaust pressure   : {HX.ex_C.p:10.0f} Pa")
# print(f"Hot side charge         : {HX.M_H * 1e3:10.2f} g")
# print(f"Cold side charge        : {HX.M_C * 1e3:10.2f} g")
# print(f" Subcooling             : {HX.ex_H.SC}")
# dp_h = HX.su_H.p-HX.ex_H.p
# print(f"Hot pressure drop       : {dp_h:10.0f} Pa")
# dp_c = HX.su_C.p-HX.ex_C.p
# print(f"Cold pressure drop      : {dp_c:10.0f} Pa")

# HX.plot_cells()