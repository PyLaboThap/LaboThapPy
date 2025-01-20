# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from simulation_model import Spliter

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Geometry"

Spliter = Spliter(outlet_repartition = [0.3,0.4,0.3])

Spliter.set_inputs(
    su_T=10 + 273.15,  # Temperature in Kelvin
    su_m_dot=13.8,     # Mass flow rate in kg/s
    su_p=0.8 * 1e5,    # Pressure in Pa
    su_fluid="Cyclopentane"
)

Spliter.solve()