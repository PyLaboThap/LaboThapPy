# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:43:39 2023

@author: Basile
"""

import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.component.heat_exchanger.steady_state.heat_pipe_based.simulation_model import HP_HTX

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

HP_HTX = HP_HTX()

# Inputs
HP_HTX.set_inputs(
                  su_C_fluid = 'Cyclopentane',
                  su_C_T = 434,
                  su_C_m_dot = 13.76,
                  su_C_p = 31*1e5,
                   
                  su_H_fluid = 'Air',
                  su_H_T = 465+273.15,
                  su_H_m_dot = 28.37,
                  su_H_p = 101.8*1e3,                  
                  )

# Params
HP_HTX.set_parameters(
                      p_CO2 = 0.176, p_H2O = 0.101, beta = np.pi/2, D_o = 20/1000, t = 3.5/1000, 
                      F_r = 0.6, k_pipe = 42.55, geo = "annular", H_core = 3, L_core = 5, W_core = 3, Bank_side = 'H',
                      coef_evap = 0.64, foul = 0.2, arrang = "inline", pitch_T = 2, pitch_L = 2, D_chimney = 2.2
                      )

# Solve
HP_HTX.solve(20,300,230 + 273.15, 245 + 273.15)
