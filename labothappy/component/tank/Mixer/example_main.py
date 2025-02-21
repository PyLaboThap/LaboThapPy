# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from simulation_model import Mixer

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Instanciation"

Mixer = Mixer(n_inlets = 2)
 
Mixer.set_inputs(
    su_1_T = 50 + 273.15, # K
    su_1_p = 2*1e5, # Pa
    su_1_m_dot = 1, # kg/s
    su_1_fluid = 'Water',
    
    su_2_T = 100 + 273.15, # K
    su_2_p = 2*1e5, # Pa
    su_2_m_dot = 1, # kg/s
    su_2_fluid = 'Water'
    )

Mixer.solve()
