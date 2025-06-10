# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

import os
import sys

# Add labohappy root to the path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from component.tank.tank_spliter import Spliter

# 1) Data ------------------------------------------------------------------------------------------

# Create splitter with specified outlet repartition
spliter = Spliter(outlet_repartition=[0.3, 0.4, 0.3])

# Set inputs
spliter.set_inputs(
    su_T=10 + 273.15,        # Temperature in Kelvin
    su_m_dot=13.8,           # Mass flow rate in kg/s
    su_p=0.8 * 1e5,          # Pressure in Pa
    su_fluid="Cyclopentane"  # Working fluid
)

# Solve
spliter.solve()

# You can also print results to verify
for i in range(len(spliter.outlet_repartition)):
    outlet = getattr(spliter, f"ex_{i+1}")
    print(f"Outlet {i+1}: m_dot = {outlet.m_dot} kg/s, p = {outlet.p} Pa, h = {outlet.h} J/kg")
