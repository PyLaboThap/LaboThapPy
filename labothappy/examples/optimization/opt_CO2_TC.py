# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:04:10 2025

@author: Basile
"""

from examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import basic_CO2_HP, IHX_CO2_HP, Flash_CO2_HP_Series_CP
from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

"1) Initiate the heat pump cycle"

T_cold_source = 0+273.15
T_hot_source = 15+273.15

eta_is_cp = 0.7
eta_gc = 0.9
eta_IHX = 0.3

PPTD_ev = 5
SH_ev = 0.1

P_high = 140*1e5
P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
P_crit_CO2 = PropsSI('PCRIT','CO2')

P_low_guess = 0.8*min(P_sat_T_CSource,P_crit_CO2)

HSource = MassConnector()
HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = 62.5) # 0.1 # 62.5

CSource = MassConnector()
CSource.set_properties(fluid = 'R22', T = T_cold_source, p = 5e5, m_dot = 625000) # 1000 # 625000

CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, PPTD_ev, SH_ev, P_low_guess, P_high, 100) # 0.16 # 100

CO2_HP.solve()

COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/CO2_HP.components['Compressor'].model.W_mec.W_dot

T_hw_out = CO2_HP.components['GasCooler'].model.ex_C.T

T_sat_su_cp = PropsSI('T', 'P', CO2_HP.components['Compressor'].model.su.p, 'Q', 0, 'CO2')
SH_cp = CO2_HP.components['Compressor'].model.su.T - T_sat_su_cp

print(f"COP : {COP}")
print(f"T_hw_out : {T_hw_out}")
print(f"CP_SH : {SH_cp}")

# "2) Sensitivity Analysis of high pressure"

# P_high_vec = np.linspace(80,200,25)*1e5
# COP_vec = np.zeros(len(P_high_vec))
# W_cp_vec = np.zeros(len(P_high_vec))
# T_out_hw = np.zeros(len(P_high_vec))
# eta_gc_vec = np.zeros(len(P_high_vec))
# eta_IHX_vec = np.zeros(len(P_high_vec))

# for i in range(len(P_high_vec)):
#     CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, PPTD_ev, SH_ev, P_low_guess, P_high_vec[i])
    
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
# plt.plot(P_high_vec, COP_vec)
# plt.xlabel("P_cp_ex")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, W_cp_vec)
# plt.xlabel("P_cp_ex")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, T_out_hw)
# plt.xlabel("P_cp_ex")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, eta_gc_vec, color='r', label='gc')
# plt.plot(P_high_vec, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("P_cp_ex")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "3) Sensitivity Analysis of Heat Sink Temperature"

# T_HSink_vec = np.linspace(5,30,26) + 273.15
# COP_vec = np.zeros(len(T_HSink_vec))
# W_cp_vec = np.zeros(len(T_HSink_vec))
# T_out_hw = np.zeros(len(T_HSink_vec))
# eta_gc_vec = np.zeros(len(T_HSink_vec))
# eta_IHX_vec = np.zeros(len(T_HSink_vec))

# for i in range(len(T_HSink_vec)):
    
#     HSource.set_T(T_HSink_vec[i])
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
# plt.plot(T_HSink_vec-273.15, COP_vec)
# plt.xlabel("T_Sink")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_HSink_vec-273.15, W_cp_vec)
# plt.xlabel("T_Sink")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_HSink_vec-273.15, T_out_hw)
# plt.xlabel("T_Sink")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_HSink_vec-273.15, eta_gc_vec, color='r', label='gc')
# plt.plot(T_HSink_vec-273.15, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("T_Sink")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "4) Sensitivity Analysis of Heat Source Temperature"

# T_Source_vec = np.linspace(0,30,31) + 273.15
# COP_vec = np.zeros(len(T_Source_vec))
# W_cp_vec = np.zeros(len(T_Source_vec))
# T_out_hw = np.zeros(len(T_Source_vec))
# eta_gc_vec = np.zeros(len(T_Source_vec))
# eta_IHX_vec = np.zeros(len(T_Source_vec))

# for i in range(len(T_Source_vec)):
    
#     HSource.set_T(T_Source_vec[i])
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
# plt.plot(T_Source_vec-273.15, COP_vec)
# plt.xlabel("T_HSource")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, W_cp_vec)
# plt.xlabel("T_HSource")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, T_out_hw)
# plt.xlabel("T_HSource")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(T_Source_vec-273.15, eta_gc_vec, color='r', label='gc')
# plt.plot(T_Source_vec-273.15, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("T_HSource")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "5) Sensitivity Analysis of Heat Source Temperature"

# eta_cp_vec = np.linspace(0.6,0.9,31)
# COP_vec = np.zeros(len(eta_cp_vec))
# W_cp_vec = np.zeros(len(eta_cp_vec))
# T_out_hw = np.zeros(len(eta_cp_vec))
# eta_gc_vec = np.zeros(len(eta_cp_vec))
# eta_IHX_vec = np.zeros(len(eta_cp_vec))

# eta_carnot =  np.zeros(len(eta_cp_vec)) 


# for i in range(len(eta_cp_vec)):
    
#     CO2_HP = IHX_CO2_HP(HSource, CSource, eta_cp_vec[i], eta_gc, eta_IHX, PPTD_ev, SH_ev, P_low_guess, P_high)
    
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
#                 eta_carnot[i] = 0.5*(1 - T_cold_source/T_out_hw[i])

# plt.figure()
# plt.plot(eta_cp_vec, COP_vec)
# plt.xlabel("eta_cp")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_cp_vec, W_cp_vec)
# plt.xlabel("eta_cp")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_cp_vec, T_out_hw)
# plt.xlabel("eta_cp")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_cp_vec, eta_gc_vec, color='r', label='gc')
# plt.plot(eta_cp_vec, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("eta_cp")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()


# plt.figure()
# plt.plot(eta_cp_vec, COP_vec*eta_carnot, color='r', label='gc')
# plt.xlabel("eta_cp")
# plt.ylabel("eta_carnot")
# plt.grid()
# plt.show()

# "6) Sensitivity Analysis of PPTD"

# PPTD_vec = np.linspace(1,10,19)
# COP_vec = np.zeros(len(PPTD_vec))
# W_cp_vec = np.zeros(len(PPTD_vec))
# T_out_hw = np.zeros(len(PPTD_vec))
# eta_gc_vec = np.zeros(len(PPTD_vec))
# eta_IHX_vec = np.zeros(len(PPTD_vec))

# for i in range(len(PPTD_vec)):
    
#     CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, PPTD_vec[i], SH_ev, P_low_guess, P_high)
    
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
# plt.plot(PPTD_vec, COP_vec)
# plt.xlabel("PPTD")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, W_cp_vec)
# plt.xlabel("PPTD")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, T_out_hw)
# plt.xlabel("PPTD")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(PPTD_vec, eta_gc_vec, color='r', label='gc')
# plt.plot(PPTD_vec, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("PPTD")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "7) Sensitivity Analysis of SH"

# SH_vec = np.linspace(0.5,10,19)
# COP_vec = np.zeros(len(SH_vec))
# W_cp_vec = np.zeros(len(SH_vec))
# T_out_hw = np.zeros(len(SH_vec))
# eta_gc_vec = np.zeros(len(SH_vec))
# eta_IHX_vec = np.zeros(len(SH_vec))

# for i in range(len(SH_vec)):
    
#     CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, PPTD_ev, SH_vec[i], P_low_guess, P_high)
    
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
# plt.plot(SH_vec, COP_vec)
# plt.xlabel("SH")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SH_vec, W_cp_vec)
# plt.xlabel("SH")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SH_vec, T_out_hw)
# plt.xlabel("SH")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(SH_vec, eta_gc_vec, color='r', label='gc')
# plt.plot(SH_vec, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("SH")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()

# "8) Sensitivity Analysis of heat sink flowrate"

# m_dot_HS_vec = np.linspace(0.05,1,20)
# COP_vec = np.zeros(len(m_dot_HS_vec))
# W_cp_vec = np.zeros(len(m_dot_HS_vec))
# T_out_hw = np.zeros(len(m_dot_HS_vec))
# eta_gc_vec = np.zeros(len(m_dot_HS_vec))
# eta_IHX_vec = np.zeros(len(m_dot_HS_vec))

# for i in range(len(m_dot_HS_vec)):
    
#     HSource.set_m_dot(m_dot_HS_vec[i])
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
# plt.xlabel("m_dot_sink")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, W_cp_vec)
# plt.xlabel("m_dot_sink")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, T_out_hw)
# plt.xlabel("m_dot_sink")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(m_dot_HS_vec, eta_gc_vec, color='r', label='gc')
# plt.plot(m_dot_HS_vec, eta_IHX_vec, color='b', label='IHX')
# plt.xlabel("m_dot_sink")
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

# "10) Sensitivity Analysis of IHX efficiency"

# eta_IHX_vec = np.linspace(0.3,1,15)
# COP_vec = np.zeros(len(eta_IHX_vec))
# W_cp_vec = np.zeros(len(eta_IHX_vec))
# T_out_hw = np.zeros(len(eta_IHX_vec))
# eta_gc_vec = np.zeros(len(eta_IHX_vec))

# for i in range(len(eta_IHX_vec)):
    
#     CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX_vec[i], PPTD_ev, SH_ev, P_low_guess, P_high)
    
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
# plt.plot(eta_IHX_vec, COP_vec)
# plt.xlabel("eta_IHX_vec")
# plt.ylabel("COP")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_IHX_vec, W_cp_vec)
# plt.xlabel("eta_IHX_vec")
# plt.ylabel("W_dot_cp")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_IHX_vec, T_out_hw)
# plt.xlabel("eta_IHX_vec")
# plt.ylabel("T_out_hw")
# plt.grid()
# plt.show()

# plt.figure()
# plt.plot(eta_IHX_vec, eta_gc_vec, color='r', label='gc')
# plt.xlabel("eta_IHX_vec")
# plt.ylabel("eta_HX")
# plt.grid()
# plt.show()
