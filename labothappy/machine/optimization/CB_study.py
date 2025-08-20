# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 16:18:35 2025

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

W_dot_cp = 1e6 # W
W_dot_exp = 1e6 # W

V_st = 50 # m3
T_st_vec = 273.15 + np.array([100,110,120,130,140,150])
P_st = 10*1e5
loss_per_cycle_st = 0.05

RTE = []
t_ch = []
t_dis = []
t_ratio = []
Q_loss_st = []
st_dens = []

for i in range(len(T_st_vec)):
    T_st = T_st_vec[i]
    
    rho_st = PropsSI('D', 'T', T_st, 'P', P_st, 'Water')
    m_st = V_st*rho_st
    
    h_st = PropsSI('H', 'T', T_st, 'P', P_st, 'Water')
    h_st_empty = PropsSI('H', 'T', 288.15, 'P', P_st, 'Water')
    
    Dh_st = h_st - h_st_empty
    
    Q_st_tot = Dh_st*m_st/3600 # Wh
    st_dens.append(Q_st_tot/V_st)
        
    COP_vec = np.array([4.697249244448596, 4.3895555255356085, 4.050117334016632, 4.026683891241342, 3.7880968812034888, 3.5983675455471418])
    eta_vec = np.array([0.09086802588180352, 0.10128587354579623, 0.11146548353387467, 0.1214949488715974, 0.13104562249758014, 0.140919652539438])
    
    Q_dot_ch = COP_vec[i]*W_dot_cp
    t_ch.append(Q_st_tot/Q_dot_ch) # h
    
    # Charge process
    E_cp = W_dot_cp * t_ch[i]
    Q_ch = COP_vec[0]*E_cp
    
    # Storage losses
    Q_loss_st.append((loss_per_cycle_st)*Q_ch)
    
    Q_dis = Q_ch*(1-loss_per_cycle_st) # 5%/cycle
    
    # Discharge Process
    Q_dot_dis = W_dot_exp/eta_vec[i]
    t_dis.append(Q_dis/Q_dot_dis) # h
    
    E_exp = W_dot_exp*t_dis[i]
    
    # RTE
    RTE.append(E_exp/E_cp)
    t_ratio.append(t_dis[i]/t_ch[i])

# Convert temperatures to Celsius for better x-axis readability
T_st_C = T_st_vec - 273.15

# Plot 1: Round-Trip Efficiency (RTE)
plt.figure()
plt.plot(T_st_C, np.array(RTE)*100, marker='o')
plt.xlabel('Storage Temperature [°C]')
plt.ylabel('Round Trip Efficiency [%]')
plt.title('Round Trip Efficiency vs Storage Temperature')
plt.grid(True)

# Plot 2: Charge Time
plt.figure()
plt.plot(T_st_C, t_ch, marker='o')
plt.xlabel('Storage Temperature [°C]')
plt.ylabel('Charge Time [h]')
plt.title('Charge Time vs Storage Temperature')
plt.grid(True)

# Plot 3: Discharge Time
plt.figure()
plt.plot(T_st_C, t_dis, marker='o')
plt.xlabel('Storage Temperature [°C]')
plt.ylabel('Discharge Time [h]')
plt.title('Discharge Time vs Storage Temperature')
plt.grid(True)

# Plot 4: Charge/Discharge Time Ratio
plt.figure()
plt.plot(T_st_C, t_ratio, marker='o')
plt.xlabel('Storage Temperature [°C]')
plt.ylabel('Charge/Discharge Time Ratio')
plt.title('Discharge/Charge Time Ratio vs Storage Temperature')
plt.grid(True)

# Plot 5: Storage Heat Loss per Cycle
plt.figure()
plt.plot(T_st_C, np.array(Q_loss_st)/1e6, marker='o')
plt.xlabel('Storage Temperature [°C]')
plt.ylabel('Heat Loss per Cycle [MWh]')
plt.title('Storage Heat Loss vs Storage Temperature')
plt.grid(True)

# Plot 6: Storage Energy Density
plt.figure()
plt.plot(T_st_C, np.array(st_dens)/1e6, marker='o')
plt.xlabel('Storage Temperature [°C]')
plt.ylabel('Energy Density [MWh/m³]')
plt.title('Storage Energy Density vs Temperature')
plt.grid(True)

plt.show()

