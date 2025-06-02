# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 13:54:18 2023

@author: Basile
"""

import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI
from modules.HP_tube_model import thermosyphon_model
from modules.Airflow import nozzle
from modules.P_max_steel_pipes import P_max_adm
from post_process import res_encoding

"""
Code used to vary parametersf a heat pipe bank and collect results
"""

encode = 0 # Decides if the results shall be encoded for later postt processing and optimization

"1) Thermosyphon design parameters"

# Sizing [m] : Oustide diameter
D_o_vect = np.array([20/1000]) # np.linspace(40/1000, 20/1000, 6) # np.linspace(20/1000,50/1000,16) # np.array([38/1000]) # [m]
t = 3.5/1000  # [m] : Thickness

L_pt = np.array([2]) # np.linspace(2, 6, 9) #  np.array([5]) #  np.linspace(0.6,2.6,11)

# Thermal
# Thermal conductivity of thermosyphon pipes (steel)
k_pipe = 42.55 # [W/(m*K)]
foul = 0.2 # fouling factor 

T_lim_tube = 350 + 273.15 # [K]

# Other
F_r = 0.6  # filling ratio
beta = np.pi/2  # [rad] : Inclination of the thermosyphon
geo = 'annular'  # thermosyphon geometry

"2) Tube matrix parameters"

arrang = 'Inline'  # 'Staggered'
pitch_ratio = 2 # 2 # 1.5 # 1.25 
pitch_ratio_2 = 2 # 2 # 1.5 # 1.25 

L_HTX = np.array([3.1]) # np.linspace(2, 6, 6) #  np.array([6]) # [m] : Heat exchanger length
W_HTX = np.array([4.5]) # np.linspace(2, 6, 9) # np.array([3.2, 3.8, 5, 5.6, 6.2, 6.8, 7.4, 8]) # [m] : Heat exchanger width

# Done : All but 6

var_array_1 = D_o_vect
var_array_2 = L_pt
var_array_3 = L_HTX
var_array_4 = W_HTX

Q_dot_tot = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
DP_air = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
DP_oil = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
u_gas_fume_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
V_max_air_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
sonic_vel_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
N_tubes_row_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
P_in_max = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
T_in_max = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
P_max_adm_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
R_tot_tube_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
R_e_o_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
R_c_o_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
DP_gas_nozzle_vect = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
N_T = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
coef_chosen = np.zeros([len(var_array_1), len(var_array_2), len(var_array_3), len(var_array_4)])
M_tubes = np.zeros([len(var_array_1), len(var_array_2)])

coef_vect = np.linspace(0.9,0.1,21) # np.array([0.6]) # 

Q_dot_tot_MW = 0
Q_dot_tot_MW_2 = 0

for m in range(len(var_array_4)):
    for l in range(len(var_array_3)):
        for k in range(len(var_array_2)):
            print('Progress (1) :',  m+1, '/', len(var_array_4),' || Progress (2) :', l+1, '/', len(var_array_3), ' || Progress (3) :', k+1, '/', len(var_array_2) )
            
            # Lines that ignore certain iterations as it is not exxpected that the ones that are following will yield good results
            
            if k > 0:
                if np.isnan(Q_dot_tot_MW_2) == 0 and Q_dot_tot_MW_2 < 5:
                    print("skip SS - L_pt : ",Q_dot_tot_MW_2)
                    Q_dot_tot_MW_2 = np.nan
                    continue
            for j in range(len(var_array_1)):
                if j > 0 :
                    if np.isnan(Q_dot_tot_MW) == 0 and Q_dot_tot_MW < 6.5:
                        print("skip - D_o : ",Q_dot_tot_MW)
                        Q_dot_tot_MW = np.nan
                        continue
                    if np.isnan(Q_dot_tot_MW) == 0 and Q_dot_tot_MW > 8:
                      print("break H - D_o : ",Q_dot_tot_MW)
                      Q_dot_tot_MW = np.nan
                      break
                  
                for c in range(len(coef_vect)):
                    # try:
                        coef = coef_vect[c]
                        coef = 0.62
                        L_a = L_pt[k]/20  # 0.05  # [m] : adiabatic zone length
                        L_e = coef*(L_pt[k]-L_a) # 0.75 # np.linspace(0.1, 2, 20)  # 1.3 # 0.756  # [m] : evap zone length
                        L_c = (1-coef)*(L_pt[k]-L_a) # 0.75 # np.linspace(0.1, 2, 20)  # 0.75 # 0.17  # [m] : cond zone length
                        
                        L_eff = L_a + (L_e + L_c)/2  # [m] : effective thermosyphon length
                        
                        D_o = D_o_vect[j] # !!! Attention au "0"
                    
                        if arrang == 'Inline':
                            S_T = pitch_ratio * D_o  # + 2*h_fins [m] : Tube pitch (écartement)
                            S_L = pitch_ratio * D_o  # + 2*h_fins [m] : Tube pitch
                        else:
                            S_T = pitch_ratio * D_o  # + 2*h_fins [m] : Tube pitch (écartement)
                            S_L = pitch_ratio_2 * D_o  # + 2*h_fins [m] : Tube pitch
                    
                        N_col = int(np.floor(L_HTX[l]/(S_T))-1)
                        
                        L_M = 1.08*(S_T*S_L-0.785*D_o**2)/D_o  # mean beam length
                    
                        D_i = D_o - 2*t  # [m] : Inside diameter
                    
                        N_tubes_row = np.floor((W_HTX[m] - 0.012)/S_T)
                        N_tubes_row_vect[j] = N_tubes_row
                    
                        # Evaporator
                        W_casing = W_HTX[m]  # [m] : cross-sectional width of thermosyphon evaporator
                        # [m] : cross-section height of thermosyphon evaporator
                        h_casing = L_e + 0.016
                        A_casing = W_casing*h_casing  # m^2 : cross-section area of thermosyphon evaporator
                    
                        N_tubes = N_tubes_row*np.ones(N_col)
                        N_T[j] = sum(N_tubes)
                    
                        # Vector for containing
                        T_gas_fume = np.zeros(N_col+1)
                        P_gas_fume = np.zeros(N_col+1)
                        
                        T_cyclo = np.zeros(N_col+1)
                        h_cyclo = np.zeros(N_col+1)
                        
                        Q_dot_row = np.zeros(N_col)
                        Q_dot_tube = np.zeros(N_col)
                        T_v_e = np.zeros(N_col)
                        
                        T_wc_o = np.zeros(N_col)
                        T_we_o = np.zeros(N_col)
                        
                        P_v_c = np.zeros(N_col)
                        P_v_e = np.zeros(N_col)
                        
                        m_dot_v = np.zeros(N_col)
                        
                        R_tot_tube = np.zeros(N_col)  
                        R_c_o = np.zeros(N_col)  
                        R_e_o = np.zeros(N_col)  
                        
                        "2.1) Total tube mass"
                        # Tube volume
                        L_tot = L_e + L_c + L_a # [m]
                        A_mat_tube = np.pi*((D_o/2)**2 - (D_i/2)**2) # [m^2]
                        Vol_cyl = L_tot*A_mat_tube # [m^3]
                        
                        # Ends volume
                        Vol_end = 4/3 * np.pi * ((D_o/2)**3 - (D_i/2)**3) # [m^3]
                        Vol_tube = Vol_cyl + 2*Vol_end # [m^3]
                        
                        Vol_water = F_r*np.pi*(D_i/2)**2*L_tot
                        
                        # For all tubes
                        rho_water = 1000 # [kg/m^3]
                        rho_steel = 8050 # [kg/m^3]
                        
                        Vol_tot_tube = Vol_tube*N_T[j][k] # [m^3]
                        Vol_tot_water = Vol_water*N_T[j][k] # [m^3]

                        M_tubes[j][k] = Vol_tot_tube*rho_steel # [kg]
                        M_water = Vol_tot_water*rho_water # [kg]

                        M_tot_tube = M_tubes[j][k] + M_water # [kg]
                        
                        "3) External fluids data"
                    
                        "3.1) Used fluids inside and outside the thermosyphon"
                    
                        fluid_thermosiphon = 'Water'
                        fluid_ev = 'Air'
                        fluid_cd = 'Cyclopentane'
                    
                        "3.2) Fume gases"
                    
                        D_chimney = 2.2  # [m] : Fume gas chimney diameter
                        A_chimney = (np.pi/4)*D_chimney**2
                    
                        # Thermo state related :
                        # [K] : fume gas supply temperature (Given by Sintef)
                        T_gas_fume_in = 465 + 273.15
                        # Pa : fume gas supply pressure (Given by Sintef)
                        P_gas_fume_in = 101.8*1e3
                        cp_gas_fume_in = PropsSI('C', 'T', T_gas_fume_in,
                                                 'P', P_gas_fume_in, fluid_ev)
                    
                        # T_gas_fume_out = 150 + 273.15  # [K] : Constraint given by Sintef
                    
                        # Components partial pressures
                        p_CO2 = 0.176  # CO2 partial pressure in fume gases : Given by Sintef
                        p_H2O = 0.101  # H2O partial pressure in fume gases : Given by Sintef
                    
                        # Flow rate related
                        m_dot_gas_fume = 28.37  # [kg/s]  : Given by Sintef
                        V_dot_gas_fume = m_dot_gas_fume / \
                            PropsSI('D', 'T', T_gas_fume_in, 'P', P_gas_fume_in, fluid_ev)
                        u_gas_fume_chimney = V_dot_gas_fume/A_chimney
                    
                        "Nozzle computation - HTX Inlet"
                        (u_gas_fume, T_gas_fume[0], P_gas_fume[0]) = nozzle(
                            m_dot_gas_fume, T_gas_fume_in, cp_gas_fume_in, u_gas_fume_chimney, P_gas_fume_in, A_chimney, A_casing)
                            
                        DP_nozzle_in = P_gas_fume_in - P_gas_fume[0]
                    
                        "3.3) Evaporating Cyclopentane"
                    
                        P_sat_cyclo = 3100*1e3 # [Pa]
                        T_sat_cyclo = PropsSI('T','P',P_sat_cyclo,'Q',0,fluid_cd)
                        
                        T_cyclo[0] = 273.15 + 240 # [K]
                        M_dot_cyclo = 13.76 # [kg/s]
                        
                        D_oil_pipe = 0.1  # [m] : oil pipe diameter
                        A_oil_pipe = (np.pi*D_oil_pipe**2)/4  # [m^2] : oil pipe area
                    
                        H_casing = L_c + 0.02  # [m] height of oil divergent
                        A_casing = H_casing*W_casing
                    
                        rho_cyclo_in = PropsSI('D', 'P',P_sat_cyclo,'T',T_cyclo[0],fluid_cd)  # [kg/m^3] : oil density at 110°C
                        V_dot_cyclo = M_dot_cyclo/rho_cyclo_in  # [m^3/s] : Oil volumetric flow rate
                    
                        u_cyclo_in = V_dot_cyclo/A_casing  # [m/s] : Oil velocity in the casing
                        
                        "4) Thermosyphon Simulation"
                    
                        u_gas_fume_vect[j] = u_gas_fume
                        
                        vect_init = (7*101325, 7*101325, 1000, 1000, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 0, 513.15, 635.15, 101325,50) # First initial conditions
                        # P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
                        
                        DP_a = 0
                        
                        for i in range(N_col):
                            
                            R_tot_tube_mean = 0
                            R_e_o_mean = 0
                            R_c_o_mean = 0
                            
                            (Q_dot_tube[i], T_wc_o[i], T_we_o[i], P_v_e[i], T_v_e[i], m_dot_v[i], P_gas_fume[i+1], V_max_air, a_gas_fume, R_tot_tube[i], R_e_o[i], R_c_o[i]) = thermosyphon_model(beta, D_i, D_o, F_r, fluid_cd, fluid_ev, fluid_thermosiphon, geo, k_pipe, L_a, L_c, L_e, L_M, m_dot_gas_fume, M_dot_cyclo, N_col, p_CO2, p_H2O, P_sat_cyclo, P_gas_fume[i], S_L, S_T, T_cyclo[i], T_gas_fume[i], u_gas_fume, u_cyclo_in, A_casing, W_casing, arrang, h_cyclo[i], foul, vect_init)  # call of the model for one tube

                            if i == 0:
                                P_in_max[j][k][l][m] = P_v_e[i]
                                T_in_max[j][k][l][m] = T_v_e[i]
                                P_max_adm_vect[j][k][l][m] = P_max_adm(D_o_vect[j]*1000, T_we_o[i], 0) # !!! attention au "0"
                            
                            if T_in_max[j][k][l][m] > T_lim_tube:
                                break
                                
                            DP_air[j][k][l][m] = DP_air[j][k][l][m] + P_gas_fume[i]-P_gas_fume[i+1]
                    
                            Q_dot_row[i] = Q_dot_tube[i]*N_tubes[i]
                    
                            cp_g = PropsSI('C', 'P', P_gas_fume[i], 'T', T_gas_fume[i], fluid_ev)
                            T_gas_fume[i+1] = T_gas_fume[i] - Q_dot_row[i]/(cp_g*m_dot_gas_fume)
                            
                            if i == 0:
                                h_cyclo[i] = PropsSI('H', 'P', P_sat_cyclo, 'T', T_cyclo[i], fluid_cd)
                                
                            h_cyclo_row = h_cyclo[i] - Q_dot_row[i]/M_dot_cyclo
                            h_cyclo[i+1] = h_cyclo_row
                            
                            T_cyclo[i+1] = PropsSI('T', 'P', P_sat_cyclo, 'H', h_cyclo_row, fluid_cd)
                            
                            R_tot_tube_mean = R_tot_tube_mean + R_tot_tube[i]
                            R_e_o_mean = R_e_o_mean + R_e_o[i]
                            R_c_o_mean = R_c_o_mean + R_c_o[i]
                            
                            vect_init = (P_v_e[i],P_v_e[i],1000,Q_dot_tube[i],T_v_e[i],T_v_e[i],T_wc_o[i],T_wc_o[i],T_wc_o[i],T_we_o[i],T_we_o[i],T_we_o[i], 0, T_cyclo[i+1], T_gas_fume[i+1], P_gas_fume[i+1],R_tot_tube[i])                            
                                        # P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
            
                            # print("T_gas_fume",T_gas_fume[i+1])
                            # print("T_loop_oil",T_loop_oil[i+1])
                            
                            print("Progress : ",i + 1,"/",N_col)
                            
                        R_tot_tube_vect[j][k][l][m] = R_tot_tube_mean/N_col
                        R_e_o_vect[j][k][l][m] = R_e_o_mean/N_col
                        R_c_o_vect[j][k][l][m] = R_c_o_mean/N_col
                        coef_chosen[j][k][l][m] = coef
                        
                        Q_dot_tot[j][k][l][m] = sum(Q_dot_row)
                        Q_dot_tot_MW = Q_dot_tot[j][k][l][m]/1e6
                        Q_dot_tot_MW2 = Q_dot_tot[j][k][l][m]/1e6
                    
                        V_max_air_vect[j][k][l][m] = V_max_air
                        sonic_vel_vect[j][k][l][m] = a_gas_fume
                    
                        N_tube_tot = sum(N_tubes)
                    
                        T_bar_we_o = sum(T_we_o)/len(T_we_o)
                        T_bar_wc_o = sum(T_wc_o)/len(T_wc_o)
                                            
                        rho_fg_out = PropsSI('D', 'T',T_gas_fume[-1],'P',P_gas_fume[-1],'air')
                        
                        V_dot_gas_fume_out = m_dot_gas_fume/rho_fg_out
                        u_gas_fume_out = V_dot_gas_fume_out/A_casing
                        
                        "Nozzle computation - HTX Inlet"
                        (u_gas_fume_out, T_gas_fume_out, P_gas_fume_out) = nozzle(
                            m_dot_gas_fume, T_gas_fume[-1], cp_g, u_gas_fume_out, P_gas_fume[-1], A_casing, A_chimney)
                        
                        # print((P_gas_fume_out - P_gas_fume[-1]))
                        
                        DP_nozzle_out = P_gas_fume_out - P_gas_fume[-1]
                        
                        DP_air[j][k][l][m] = DP_a + DP_nozzle_in + DP_nozzle_out
                
                        if T_in_max[j][k][l][m] > T_lim_tube:
                            if c == coef_vect[-1]:
                                Q_dot_tot[j][k][l][m] = np.nan
                                Q_dot_tot_MW = np.nan
                                Q_dot_tot_MW_2 = np.nan
                                DP_air[j][k][l][m] = np.nan
                                DP_oil[j][k][l][m] = np.nan
                                u_gas_fume_vect[j][k][l][m] = np.nan
                                V_max_air_vect[j][k][l][m] = np.nan
                                sonic_vel_vect[j][k][l][m] = np.nan
                                N_tubes_row_vect[j][k][l][m] = np.nan
                                P_in_max[j][k][l][m] = np.nan
                                T_in_max[j][k][l][m] = np.nan
                                P_max_adm_vect[j][k][l][m] = np.nan
                                R_tot_tube_vect[j][k][l][m] = np.nan
                                R_e_o_vect[j][k][l][m] = np.nan
                                R_c_o_vect[j][k][l][m] = np.nan
                                DP_gas_nozzle_vect[j][k][l][m] = np.nan
                                N_T[j][k][l][m] = np.nan
                                coef_chosen[j][k][l][m] = np.nan
                            else:
                                continue
                        else:
                            break
                            
                    # except:
                    #     Q_dot_tot[j][k][l][m] = np.nan
                    #     Q_dot_tot_MW = np.nan
                    #     Q_dot_tot_MW_2 = np.nan
                    #     DP_air[j][k][l][m] = np.nan
                    #     DP_oil[j][k][l][m] = np.nan
                    #     u_gas_fume_vect[j][k][l][m] = np.nan
                    #     V_max_air_vect[j][k][l][m] = np.nan
                    #     sonic_vel_vect[j][k][l][m] = np.nan
                    #     N_tubes_row_vect[j][k][l][m] = np.nan
                    #     P_in_max[j][k][l][m] = np.nan
                    #     T_in_max[j][k][l][m] = np.nan
                    #     P_max_adm_vect[j][k][l][m] = np.nan
                    #     R_tot_tube_vect[j][k][l][m] = np.nan
                    #     R_e_o_vect[j][k][l][m] = np.nan
                    #     R_c_o_vect[j][k][l][m] = np.nan
                    #     DP_gas_nozzle_vect[j][k][l][m] = np.nan
                    #     N_T[j][k][l][m] = np.nan
                    #     coef_chosen[j][k][l][m] = np.nan
                
                print('Progress 1 :', j+1, '/', len(var_array_1))

    
    if encode == 1:                        
        "5) Encode the results"
        
        W_HTX_string = str(var_array_4[m]).replace('.','_')
        
        res_encoding(Q_dot_tot[:,:,:,m],'W'+ W_HTX_string +'_Q_dot',D_o_vect,L_HTX,L_pt)
        res_encoding(DP_air[:,:,:,m],'W'+ W_HTX_string +'_DP_air',D_o_vect,L_HTX,L_pt)
        res_encoding(DP_oil[:,:,:,m],'W'+ W_HTX_string +'_DP_oil',D_o_vect,L_HTX,L_pt)
        res_encoding(P_in_max[:,:,:,m],'W'+ W_HTX_string +'_P_in_max',D_o_vect,L_HTX,L_pt)
        res_encoding(T_in_max[:,:,:,m],'W'+ W_HTX_string +'_T_in_max',D_o_vect,L_HTX,L_pt)
        res_encoding(P_max_adm_vect[:,:,:,m],'W'+ W_HTX_string +'_P_max_adm',D_o_vect,L_HTX,L_pt)
        res_encoding(R_tot_tube_vect[:,:,:,m],'W'+ W_HTX_string +'_R_tube',D_o_vect,L_HTX,L_pt)
        res_encoding(R_e_o_vect[:,:,:,m],'W'+ W_HTX_string +'_R_e',D_o_vect,L_HTX,L_pt)
        res_encoding(R_c_o_vect[:,:,:,m],'W'+ W_HTX_string +'_R_c',D_o_vect,L_HTX,L_pt)
        res_encoding(N_T[:,:,:,m],'W'+ W_HTX_string +'_N_T',D_o_vect,L_HTX,L_pt)
        res_encoding(coef_chosen[:,:,:,m],'W'+ W_HTX_string +'_coef',D_o_vect,L_HTX,L_pt)
