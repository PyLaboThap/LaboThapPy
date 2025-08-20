# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:04:10 2025

@author: Basile
"""

from machine.examples.CO2_Transcritical_Circuits.CO2_Transcritical_circuit import basic_CO2_TC, REC_CO2_TC
from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

"1) Initiate the heat pump cycle"

case_study = "Plant"

if case_study == "Test_Bench":
    T_cold_source = 0.1+273.15
    T_hot_source = 110+273.15
    
    eta_is_exp = 0.7
    eta_gh = 0.9
    eta_rec = 0.7
    eta_is_pp = 0.9
    
    m_dot = 0.04 # 0.08
    
    PPTD_cd = 10
    SC_cd = 0.1
    
    P_high = 100*1e5
    P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
    P_crit_CO2 = PropsSI('PCRIT','CO2')
    
    P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)  
    
    HSource = MassConnector()
    HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = m_dot*1) # 0.1 # 62.5
    
    CSource = MassConnector()
    CSource.set_properties(fluid = 'Water', T = T_cold_source, p = 5e5, m_dot = m_dot*100) # 1000 # 625000
    
    CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot) # 0.16 # 100
    
    CO2_TC.solve()
    
    eta = (CO2_TC.components['Expander'].model.W_exp.W_dot-CO2_TC.components['Pump'].model.W_pp.W_dot)/CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
    
    T_hw_out = CO2_TC.components['GasHeater'].model.ex_H.T
    
    print(f"eta : {eta}")
    print(f"T_hw_out : {T_hw_out}")

else:
    T_cold_source = 15+273.15
    T_hot_source = 130+273.15
    
    eta_is_exp = 0.9
    eta_gh = 0.95
    eta_rec = 0.9
    eta_is_pp = 0.95
    
    m_dot = 100 # 0.08
    
    Pinch_min_GH = 5
    Pinch_min_REC = 0
    PPTD_cd = 5
    SC_cd = 0.1
    
    P_high = 140*1e5
    P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
    P_crit_CO2 = PropsSI('PCRIT','CO2')
    
    P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)  
    
    HSource = MassConnector()
    HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = m_dot*0.75) # 0.1 # 62.5
    
    CSource = MassConnector()
    CSource.set_properties(fluid = 'Water', T = T_cold_source, p = 5e5, m_dot = m_dot*100) # 1000 # 625000
    
    CO2_TC = REC_CO2_TC(HSource, CSource.T, Pinch_min_GH, Pinch_min_REC, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot, mute_print_flag = 0) # 0.16 # 100
    
    CO2_TC.solve()
    
    eta = (CO2_TC.components['Expander'].model.W_exp.W_dot-CO2_TC.components['Pump'].model.W_pp.W_dot)/CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
    
    T_hw_out = CO2_TC.components['GasHeater'].model.ex_H.T
    
    print(f"eta : {eta}")
    print(f"T_hw_out : {T_hw_out}")

    # CO2_TC.Ts_plot()

# "2) Sensitivity Analysis of high pressure"

# P_high_vec = np.linspace(80,200,25)*1e5
# eta_vec = np.zeros(len(P_high_vec))
# W_exp_vec = np.zeros(len(P_high_vec))
# W_pp_vec = np.zeros(len(P_high_vec))
# T_out_hw = np.zeros(len(P_high_vec))
# eta_gh_vec = np.zeros(len(P_high_vec))
# eta_rec_vec = np.zeros(len(P_high_vec))
# Q_dot_gh_vec = np.zeros(len(P_high_vec))

# for i in range(len(P_high_vec)):
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high_vec[i], m_dot)
#     # CO2_TC = basic_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, PPTD_cd, SC_cd, P_low_guess, P_high_vec[i], m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             # if res <= 1e-2:
#             eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#             W_exp_vec[i] = W_dot_exp
#             W_pp_vec[i] = W_dot_pp
#             Q_dot_gh_vec[i] = Q_dot_gh
#             T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#             eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#             eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']

# plt.figure()
# plt.plot(P_high_vec, eta_vec)
# plt.xlabel("P_high")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, Q_dot_gh_vec)
# plt.xlabel("P_high")
# plt.ylabel("Q_dot_gh")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, W_exp_vec)
# plt.xlabel("P_high")
# plt.ylabel("W_dot_exp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, W_pp_vec)
# plt.xlabel("P_high")
# plt.ylabel("W_dot_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, T_out_hw)
# plt.xlabel("P_high")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, eta_gh_vec, color='r', label='GH')
# # plt.plot(P_high_vec, eta_rec_vec, color='b', label='Rec')
# plt.xlabel("P_high")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "3) Sensitivity Analysis of Heat Source Temperature"

# T_HSink_vec = np.linspace(110,150,26) + 273.15
# eta_vec = np.zeros(len(T_HSink_vec))
# W_exp_vec = np.zeros(len(T_HSink_vec))
# W_pp_vec = np.zeros(len(T_HSink_vec))
# T_out_hw = np.zeros(len(T_HSink_vec))
# eta_gh_vec = np.zeros(len(T_HSink_vec))
# eta_rec_vec = np.zeros(len(T_HSink_vec))

# for i in range(len(T_HSink_vec)):
    
#     HSource.set_T(T_HSink_vec[i])
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']       

# plt.figure()
# plt.plot(T_HSink_vec-273.15, eta_vec)
# plt.xlabel("T_high")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_HSink_vec-273.15, W_exp_vec)
# plt.xlabel("T_high")
# plt.ylabel("W_dot_exp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_HSink_vec-273.15, W_pp_vec)
# plt.xlabel("T_high")
# plt.ylabel("W_dot_pp")
# plt.grid(
# )
# plt.show()
# plt.figure()
# plt.plot(T_HSink_vec-273.15, T_out_hw)
# plt.xlabel("T_high")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_HSink_vec-273.15, eta_gh_vec, color='r', label='GH')
# # plt.plot(P_high_vec, eta_rec_vec, color='b', label='Rec')
# plt.xlabel("T_high")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "4) Sensitivity Analysis of Heat Source Temperature"

# T_Source_vec = np.linspace(0,20,11) + 273.15
# eta_vec = np.zeros(len(T_Source_vec))
# W_exp_vec = np.zeros(len(T_Source_vec))
# W_pp_vec = np.zeros(len(T_Source_vec))
# T_out_hw = np.zeros(len(T_Source_vec))
# eta_gh_vec = np.zeros(len(T_Source_vec))
# eta_rec_vec = np.zeros(len(T_Source_vec))

# for i in range(len(T_Source_vec)):
    
#     P_sat_T_CSource = PropsSI('P', 'T', T_Source_vec[i],'Q',0.5,'CO2')

#     P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)  
    
#     CSource.set_T(T_Source_vec[i])
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']    
                
# plt.figure()
# plt.plot(T_Source_vec-273.15, eta_vec)
# plt.xlabel("T_cold")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, W_exp_vec)
# plt.xlabel("T_cold")
# plt.ylabel("W_dot_exp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, W_pp_vec)
# plt.xlabel("T_cold")
# plt.ylabel("W_dot_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, T_out_hw)
# plt.xlabel("T_cold")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, eta_gh_vec, color='r', label='GH')
# # plt.plot(P_high_vec, eta_rec_vec, color='b', label='Rec')
# plt.xlabel("T_cold")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "5) Sensitivity Analysis of expander efficiency"

# eta_exp_vec = np.linspace(0.6,0.9,31)
# eta_vec = np.zeros(len(eta_exp_vec))
# W_exp_vec = np.zeros(len(eta_exp_vec))
# W_pp_vec = np.zeros(len(eta_exp_vec))
# T_out_hw = np.zeros(len(eta_exp_vec))
# eta_gh_vec = np.zeros(len(eta_exp_vec))
# eta_rec_vec = np.zeros(len(eta_exp_vec))

# for i in range(len(eta_exp_vec)):
    
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_exp_vec[i], eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']   

# plt.figure()
# plt.plot(eta_exp_vec, eta_vec)
# plt.xlabel("eta_exp")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_exp_vec, W_exp_vec)
# plt.xlabel("eta_exp")
# plt.ylabel("W_dot_exp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_exp_vec, W_pp_vec)
# plt.xlabel("eta_exp")
# plt.ylabel("W_dot_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_exp_vec, T_out_hw)
# plt.xlabel("eta_exp")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_exp_vec, eta_gh_vec, color='r', label='GH')
# # plt.plot(P_high_vec, eta_rec_vec, color='b', label='Rec')
# plt.xlabel("eta_exp")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "6) Sensitivity Analysis of PPTD"

# PPTD_vec = np.linspace(1,10,19)
# eta_vec = np.zeros(len(PPTD_vec))
# W_exp_vec = np.zeros(len(PPTD_vec))
# W_pp_vec = np.zeros(len(PPTD_vec))
# T_out_hw = np.zeros(len(PPTD_vec))
# eta_gh_vec = np.zeros(len(PPTD_vec))
# eta_rec_vec = np.zeros(len(PPTD_vec))

# for i in range(len(PPTD_vec)):
    
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_vec[i], SC_cd, P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']  

# plt.figure()
# plt.plot(PPTD_vec, eta_vec)
# plt.xlabel("PPTD_cd")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, W_exp_vec)
# plt.xlabel("PPTD_cd")
# plt.ylabel("W_dot_exp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, W_pp_vec)
# plt.xlabel("PPTD_cd")
# plt.ylabel("W_dot_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, T_out_hw)
# plt.xlabel("PPTD_cd")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, eta_gh_vec, color='r', label='GH')
# # plt.plot(P_high_vec, eta_rec_vec, color='b', label='Rec')
# plt.xlabel("PPTD_cd")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "7) Sensitivity Analysis of SH"

# SC_vec = np.linspace(0.5,10,19)
# eta_vec = np.zeros(len(SC_vec))
# W_exp_vec = np.zeros(len(SC_vec))
# W_pp_vec = np.zeros(len(SC_vec))
# T_out_hw = np.zeros(len(SC_vec))
# eta_gh_vec = np.zeros(len(SC_vec))
# eta_rec_vec = np.zeros(len(SC_vec))

# for i in range(len(SC_vec)):
    
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_vec[i], P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']  
            

# plt.figure()
# plt.plot(SC_vec, eta_vec)
# plt.xlabel("SC")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SC_vec, W_exp_vec)
# plt.xlabel("SC")
# plt.ylabel("W_exp_vec")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SC_vec, W_pp_vec)
# plt.xlabel("SC")
# plt.ylabel("W_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SC_vec, T_out_hw)
# plt.xlabel("SC")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SC_vec, eta_gh_vec, color='r', label='gc')
# plt.plot(SC_vec, eta_rec_vec, color='b', label='rec')
# plt.xlabel("SC")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "8) Sensitivity Analysis of heat sink flowrate"

# m_dot_HS_vec = np.linspace(0.05,1,20)*100
# eta_vec = np.zeros(len(m_dot_HS_vec))
# W_exp_vec = np.zeros(len(m_dot_HS_vec))
# W_pp_vec = np.zeros(len(m_dot_HS_vec))
# T_out_hw = np.zeros(len(m_dot_HS_vec))
# eta_gh_vec = np.zeros(len(m_dot_HS_vec))
# eta_rec_vec = np.zeros(len(m_dot_HS_vec))

# for i in range(len(m_dot_HS_vec)):
    
#     HSource.set_m_dot(m_dot_HS_vec[i])
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec, PPTD_cd, SC_cd, P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']  
            

# plt.figure()
# plt.plot(m_dot_HS_vec, eta_vec)
# plt.xlabel("m_dot_H")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, W_exp_vec)
# plt.xlabel("m_dot_H")
# plt.ylabel("W_exp_vec")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, W_pp_vec)
# plt.xlabel("m_dot_H")
# plt.ylabel("W_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, T_out_hw)
# plt.xlabel("m_dot_H")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, eta_gh_vec, color='r', label='gc')
# plt.plot(m_dot_HS_vec, eta_rec_vec, color='b', label='rec')
# plt.xlabel("m_dot_H")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "9) Sensitivity Analysis of heat source flowrate"

# m_dot_HS_vec = np.linspace(0.5,10,20)
# COP_vec = np.zeros(len(m_dot_HS_vec))
# W_cp_vec = np.zeros(len(m_dot_HS_vec))
# T_out_hw = np.zeros(len(m_dot_HS_vec))
# eta_gc_vec = np.zeros(len(m_dot_HS_vec))
# eta_IHX_vec = np.zeros(len(m_dot_HS_vec))

# for i in range(len(m_dot_HS_vec)):
    
#     CSource.set_m_dot(m_dot_HS_vec[i])
#     CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, PPTD_ev, SH_ev, P_low_guess, P_high)
    
#     CO2_HP.solve()

#     Q_dot_gc = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot
#     Q_dot_ev = CO2_HP.components['Evaporator'].model.Q_dot.Q_dot
#     W_dot_cp = CO2_HP.components['Compressor'].model.W_mec.W_dot

#     if CO2_HP.converged:
#         if Q_dot_gc >= 0 and W_dot_cp >= 0 and Q_dot_ev >= 0:
#             res = Q_dot_gc - W_dot_cp - Q_dot_ev
#             res = res/Q_dot_gc
            
#             if res <= 1e-3:
#                 COP_vec[i] = Q_dot_gc/W_dot_cp
#                 W_cp_vec[i] = W_dot_cp
#                 T_out_hw[i] =  CO2_HP.components['GasCooler'].model.ex_C.T
#                 eta_gc_vec[i] = CO2_HP.components['GasCooler'].model.params['eta']
#                 eta_IHX_vec[i] = CO2_HP.components['IHX'].model.params['eta']
            

# plt.figure()
# plt.plot(m_dot_HS_vec, COP_vec)
# plt.xlabel("m_dot_source")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, W_cp_vec)
# plt.xlabel("m_dot_source")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, T_out_hw)
# plt.xlabel("m_dot_source")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, eta_gc_vec, color='r', label='gc')
# plt.plot(m_dot_HS_vec, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("m_dot_source")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "10) Sensitivity Analysis of Rec efficiency"

# eta_rec_vec = np.linspace(0.3,1,15)
# eta_vec = np.zeros(len(eta_rec_vec))
# W_exp_vec = np.zeros(len(eta_rec_vec))
# W_pp_vec = np.zeros(len(eta_rec_vec))
# T_out_hw = np.zeros(len(eta_rec_vec))
# eta_gh_vec = np.zeros(len(eta_rec_vec))
# eta_rec_vec = np.zeros(len(eta_rec_vec))

# for i in range(len(eta_rec_vec)):
    
#     CO2_TC = REC_CO2_TC(HSource, CSource, eta_is_pp, eta_is_exp, eta_gh, eta_rec_vec[i], PPTD_cd, SC_cd, P_low_guess, P_high, m_dot)
    
#     CO2_TC.solve()

#     Q_dot_gh = CO2_TC.components['GasHeater'].model.Q_dot.Q_dot
#     Q_dot_cd = CO2_TC.components['Condenser'].model.Q_dot.Q_dot
#     W_dot_pp = CO2_TC.components['Pump'].model.W_pp.W_dot
#     W_dot_exp = CO2_TC.components['Expander'].model.W_exp.W_dot

#     if CO2_TC.converged:
#         if Q_dot_gh >= 0 and Q_dot_cd >= 0 and W_dot_exp >= 0 and W_dot_pp >= 0:
#             res = Q_dot_gh - Q_dot_cd - W_dot_exp + W_dot_pp
#             res = res/Q_dot_gh
            
#             if res <= 1e-2:
#                 eta_vec[i] = (W_dot_exp-W_dot_pp)/Q_dot_gh
#                 W_exp_vec[i] = W_dot_exp
#                 W_pp_vec[i] = W_dot_pp
#                 T_out_hw[i] = CO2_TC.components['GasHeater'].model.ex_H.T
#                 eta_gh_vec[i] = CO2_TC.components['GasHeater'].model.params['eta']
#                 eta_rec_vec[i] = CO2_TC.components['Recuperator'].model.params['eta']   

# plt.figure()
# plt.plot(eta_rec_vec, eta_vec)
# plt.xlabel("eta_exp")
# plt.ylabel("eta")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_rec_vec, W_exp_vec)
# plt.xlabel("eta_exp")
# plt.ylabel("W_dot_exp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_rec_vec, W_pp_vec)
# plt.xlabel("eta_exp")
# plt.ylabel("W_dot_pp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_rec_vec, T_out_hw)
# plt.xlabel("eta_exp")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_rec_vec, eta_gh_vec, color='r', label='GH')
# # plt.plot(P_high_vec, eta_rec_vec, color='b', label='Rec')
# plt.xlabel("eta_exp")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()
