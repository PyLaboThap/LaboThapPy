# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:07:43 2024

@author: Basile
"""

import __init__ 

from connector.mass_connector import MassConnector
from component.steady_state.turbomachinery.turbine.polyn_isentropic_eff.simulation_model import Turb_polyn_eff
from modules.c_turb_polyn_geom import Geometry_polyn_turb
from CoolProp.CoolProp import PropsSI   
import numpy as np 
import matplotlib.pyplot as plt

p_in = np.array([10 , 15 , 20, 25, 30, 35]) # 0.8
T_in = np.array([173  , 195.6, 213.1, 227.6, 240, 250.8])
m_dot_tune = np.array([4.74 , 7.02 , 9.29, 11.57, 13.855, 16.165])

m_dot_model = np.zeros(len(T_in))

#%%
"--------- 1) DECAGONE Turbine ------------------------------------------------------------------------------------------"

"Create turb object"
Turb = Turb_polyn_eff()

"Geometry"

Turb_geom = Geometry_polyn_turb()
Turb_geom.set_parameters("DECAGONE_TURB") 

"Set params"
Turb.set_parameters(
                    D_inlet = Turb_geom.D_inlet, N_turb_rated = Turb_geom.N_turb_rated, turb_voltage = Turb_geom.turb_voltage, turb_phases = Turb_geom.turb_phases,
                    eta_max_motor = Turb_geom.eta_max_motor, W_dot_el_rated = Turb_geom.W_dot_el_rated, eta_m = Turb_geom.eta_m, eta_is_coefs = Turb_geom.eta_is_coefs,
                    eta_is_coefs_red = Turb_geom.eta_is_coefs_red, A_th = Turb_geom.A_th 
                    )

for i in range(len(p_in)): # Loop for model flowrate verification (verify the input vs the ouput flowrate)

    "Set inputs"
    Turb.set_inputs(
                    su_fluid = 'Cyclopentane',
                    su_T = 273.15 + T_in[i],
                    su_p = p_in[i]*1e5,
                    ex_p = 0.8*1e5,
                    N_rot = 50
    )
    
    Turb.solve()
    m_dot_model[i] = Turb.m_dot

plt.figure()
plt.plot(m_dot_tune,m_dot_model)
plt.plot([0,20], [0,20])
plt.grid()
plt.show()

"""

Outputs shall be : 

W_dot = 2062 kW
T_out = 133.8 K
P_out = 0.8 bar                    

"""