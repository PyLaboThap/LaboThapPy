# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:43:39 2023

@author: Basile
"""

import numpy as np
from component.heat_exchanger.hex_thermosyphon import HP_HTX

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

HP_HTX = HP_HTX()

# Inputs
HP_HTX.set_inputs(
                  fluid_C = 'Cyclopentane',
                  T_su_C = 434,
                  m_dot_C = 13.76,
                  P_su_C = 31*1e5,
                   
                  fluid_H = 'Air',
                  T_su_H = 465+273.15,
                  m_dot_H = 28.37,
                  P_su_H = 101.8*1e3,                  
                  )

# Params
HP_HTX.set_parameters(
                      p_CO2 = 0.176, p_H2O = 0.101, beta = np.pi/2, D_o = 20/1000, t = 3.5/1000, 
                      F_r = 0.6, k_pipe = 42.55, geo = "annular", H_core = 3, L_core = 2.5, W_core = 3, Bank_side = 'H',
                      coef_evap = 0.64, foul = 0.2, arrang = "inline", pitch_T = 2, pitch_L = 2, D_chimney = 2.2, HP_fluid = 'Water'
                      )

# Solve
HP_HTX.solve(230 + 273.15, 245 + 273.15, n_it_max = 20, print_flag = 1)
HP_HTX.plot_disc()

