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
T_hot_source = 50+273.15

P_high = 120*1e5
P_mid = 50*1e5
P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
P_crit_CO2 = PropsSI('PCRIT','CO2')

P_low_guess = 0.8*min(P_sat_T_CSource,P_crit_CO2)

HSource = MassConnector()
HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 3e5, m_dot = 0.1)

CSource = MassConnector()
CSource.set_properties(fluid = 'R22', T = T_cold_source, p = 5e5, m_dot = 1)

CO2_HP = Flash_CO2_HP_Series_CP(HSource, CSource, 0.7, 0.9, 3, 3, P_low_guess, P_mid, P_high)

CO2_HP.solve()

# "2) Sensitivity Analysis of high pressure"

# P_high_vec = np.linspace(100,200,21)*1e5
# COP_vec = np.zeros(len(P_high_vec))
# W_cp_vec = np.zeros(len(P_high_vec))

# for i in range(len(P_high_vec)):
#     CO2_HP = IHX_CO2_HP(HSource, CSource, 0.7, 0.95, 0.8, 5, 3, P_low_guess, P_high_vec[i])
    
#     CO2_HP.solve()

#     Q_dot_gc = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot
#     Q_dot_ev = CO2_HP.components['Evaporator'].model.Q_dot.Q_dot
#     W_dot_cp = CO2_HP.components['Compressor'].model.W_dot.W_dot

#     if CO2_HP.converged:
#         if Q_dot_gc >= 0 and W_dot_cp >= 0 and Q_dot_ev >= 0:
#             res = Q_dot_gc - W_dot_cp - Q_dot_ev
#             res = res/CO2_HP.components['GasCooler'].model.Q_dot.Q_dot
            
#             if res <= 1e-3:
#                 COP_vec[i] = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/CO2_HP.components['Compressor'].model.W_dot.W_dot
#                 W_cp_vec[i] = CO2_HP.components['Compressor'].model.W_dot.W_dot

# plt.figure()
# plt.plot(P_high_vec, COP_vec)
# plt.xlabel("P_cp_ex")
# plt.ylabel("COP")
# plt.show()

# plt.figure()
# plt.plot(P_high_vec, W_cp_vec)
# plt.xlabel("P_cp_ex")
# plt.ylabel("W_dot_cp")
# plt.show()

# "3) Sensitivity Analysis of Heat Sink Temperature"

# T_HSource_vec = np.linspace(40,200,21) + 273.15
# COP_vec = np.zeros(len(P_high_vec))
# W_cp_vec = np.zeros(len(P_high_vec))

# Q_dot_gc = np.zeros(len(P_high_vec))
# Q_dot_ev = np.zeros(len(P_high_vec))
# W_dot_cp = np.zeros(len(P_high_vec))
# res = np.zeros(len(P_high_vec))

# for i in range(len(P_high_vec)):
    
#     HSource.set_T(T_HSource_vec[i])
#     CO2_HP = IHX_CO2_HP(HSource, CSource, 0.7, 0.95, 0.8, 5, 3, P_low_guess, P_high)
    
#     CO2_HP.solve()

#     Q_dot_gc[i] = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot
#     Q_dot_ev[i] = CO2_HP.components['Evaporator'].model.Q_dot.Q_dot
#     W_dot_cp[i] = CO2_HP.components['Compressor'].model.W_dot.W_dot

#     if CO2_HP.converged:
#         if Q_dot_gc[i] >= 0 and W_dot_cp[i] >= 0 and Q_dot_ev[i] >= 0:
#             res[i] = Q_dot_gc[i] - W_dot_cp[i] - Q_dot_ev[i]
#             res[i] = res[i]/Q_dot_gc[i]
            
#             if res[i] <= 1e-3:
#                 COP_vec[i] = Q_dot_gc[i]/W_dot_cp[i]
#                 W_cp_vec[i] = W_dot_cp[i]
            

# plt.figure()
# plt.plot(T_HSource_vec - 273.15, COP_vec)
# plt.xlabel("T_h")
# plt.ylabel("COP")
# plt.show()

# plt.figure()
# plt.plot(T_HSource_vec - 273.15, W_cp_vec)
# plt.xlabel("T_h")
# plt.ylabel("W_dot_cp")
# plt.show()

