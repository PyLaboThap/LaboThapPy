# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:00:30 2023

@author: Basile
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from matplotlib.lines import Line2D

"""
Max allowable stress for A-106C steel with respect to temperature
A-106C : Seamless Carbon Steel Pipe for High-Temperature Service

https://www.engineeringtoolbox.com/temperature-allowable-stresses-pipes-d_1338.html
"""

def P_max_adm(D_o_in,T_tube_in,plot_boolean):
    """
    Inputs
    ----------
        - D_o_in : Input outer diameter
        - T_tube_in : Input tube temperature
        - plot_boolean : Chooses if a plot shall be done
    
    Outputs
    -------
        - P_max_calc : Maximum allowable pressure
        
    Reference
    ---------
    2007 ASME BPV Code 
    
    """
    T_S_interp = np.array([0, 93.33, 204.444, 315.556, 371.111,
                          398.889, 426.667]) + 273.15  # [K] : Temperature vector
    # [MPa] : Max stress with respect to temperature vector
    S_interp = np.array([158.57942, 158.57942, 158.57942,
                        134.44777, 131.00039, 103.42136, 82.737088])
    
    S_fun = interp1d(T_S_interp, S_interp, kind='linear')
    
    """
    
    Max allowable pressure depending on pipe outside diameter and thickness
    If under critical pressure, associated stauration temperature
    
    """
    
    e = 0  # Because welded head
    
    "Compute P_max for inputs"
    
    S_tube_calc = S_fun(T_tube_in)  # [MPa]
    
    t = 3.5 # [mm]
    D_i_in = D_o_in - 2*t  # [mm]
    
    P_max_1 = S_tube_calc*((2*t - 0.01*D_o_in - 2*e)/(D_o_in - (t-0.005*D_o_in-e)))
    P_max_bouch_1 = (4.8*t*S_tube_calc)/(5*D_i_in)
    
    if P_max_1 > P_max_bouch_1:
        P_max_calc = P_max_bouch_1
    else:
        P_max_calc = P_max_1
    
    "Plot total curve"
    
    if plot_boolean == 1:
        T_tube = 373.15 + 273.15  # [K]
        S_tube = S_fun(T_tube)  # [MPa]
        D_o_vect = np.linspace(18, 58, 101)  # [mm]
        P_max = np.zeros(len(D_o_vect))
        P_max_bouch = np.zeros(len(D_o_vect))
        T_max_sat_cel = np.zeros(len(D_o_vect))
    
        P_crit_w = PropsSI('PCRIT', '',0,'',0,'water')*1e-5
        P_sat_310 = PropsSI('P', 'Q',1,'T',310+273.15,'water')*1e-5
        P_sat_240 = PropsSI('P', 'Q',1,'T',240+273.15,'water')*1e-5
    
        for i in range(len(D_o_vect)):
            D_o = D_o_vect[i]
            D_i = D_o - 2*t  # [mm]
        
            P_max[i] = S_tube*((2*t - 0.01*D_o - 2*e)/(D_o - (t-0.005*D_o-e)))
            P_max_bouch[i] = (4.8*t*S_tube)/(5*D_i)
        
            if P_max[i] > P_max_bouch[i]:
                P_max_fin = P_max_bouch[i]
            else:
                P_max_fin = P_max[i]
                
            if P_max_fin*10 < P_crit_w:
                T_max_sat = PropsSI('T', 'P', P_max_fin*1e6, 'Q', 1, 'water')
                T_max_sat_cel[i] = T_max_sat - 273.15  # [°C]
            else:
                T_max_sat_cel[i] = np.nan
                
        "Plot Results"
        
        # Create Line2D objects and labels for each line
        line1 = Line2D([0], [0], color='blue', label='P$_{max,tube}$')
        line2 = Line2D([0], [0], color='black', label='P$_{max,head}$')
        line3 = Line2D([0], [0], color='orange', label='P$_{c,w}$')
        line4 = Line2D([0], [0], color='red', label='P$_{sat,310}$')
        
        plt.figure()
        
        plt.plot(D_o_vect,P_max*10, color = 'blue')
        plt.plot(D_o_vect,P_max_bouch*10, color = 'black')
        plt.plot([D_o_vect[0],D_o_vect[-1]],[P_crit_w, P_crit_w],color = 'orange')
        plt.plot([D_o_vect[0],D_o_vect[-1]],[P_sat_310, P_sat_310],color = 'red')
        
        plt.grid()
        plt.xlabel("D$_o$ [mm]", fontsize=18)
        plt.ylabel("P$_{op,max} [bar]$", fontsize=18)
        plt.xlim([D_o_vect[1],D_o_vect[-1]])
        plt.ylim([0,700])
        
        x_ticks = np.linspace(D_o_vect[0],D_o_vect[-1],11)
        y_ticks = np.linspace(0,700,8)
        
        plt.xticks(x_ticks, rotation = 45, fontsize=16)
        plt.yticks(y_ticks, fontsize = 18)
        
        # Create a legend with a custom bounding box anchor
        # The bbox_to_anchor parameter takes a tuple (x, y) for the anchor position.
        # In this example, it's positioned in the top-right corner of the plot.
        plt.legend(handles=[line1, line2, line3, line4], loc='upper center', bbox_to_anchor=(0.5, 1.3),ncol = 4, fontsize = 14)
    
    return P_max_calc

