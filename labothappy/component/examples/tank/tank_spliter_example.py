# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from labothappy.component.tank.tank_spliter import TankSpliter

# 1) Data ------------------------------------------------------------------------------------------

# Create splitter with specified outlet repartition
spliter = TankSpliter(outlet_repartition=[0.3, 0.4, 0.3])

# Set inputs
spliter.set_inputs(
    T_su=10 + 273.15,        # Temperature in Kelvin
    m_dot=13.8,           # Mass flow rate in kg/s
    P_su=0.8 * 1e5,          # Pressure in Pa
    fluid="Cyclopentane"  # Working fluid
)

# Solve
spliter.solve()

# You can also print results to verify
for i in range(len(spliter.outlet_repartition)):
    outlet = getattr(spliter, f"ex_{i+1}")
    print(f"Outlet {i+1}: m_dot = {outlet.m_dot} kg/s, p = {outlet.p} Pa, h = {outlet.h} J/kg")

# Plot States
fig = spliter.plot_thermo_states()
fig.show()

