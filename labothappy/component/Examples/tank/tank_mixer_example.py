# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from component.tank.tank_mixer import Mixer 
from CoolProp.CoolProp import PropsSI


"--------- 1) Data ------------------------------------------------------------------------------------------"

"Instanciation"

Mixer = Mixer(n_inlets = 2)
 
Mixer.set_inputs(
    T_su_1 = 50 + 273.15, # K
    p_su_1 = 2*1e5, # Pa
    m_dot_su_1 = 1, # kg/s
    fluid_su_1 = 'Water',
    
    T_su_2 = 100 + 273.15, # K
    p_su_2 = 2*1e5, # Pa
    m_dot_su_2 = 1, # kg/s
    fluid_su_2 = 'Water'
    )

"--------- 2) Solve ------------------------------------------------------------------------------------------"
Mixer.solve()


"--------- 3) Results ------------------------------------------------------------------------------------------"
print("\n=== Mixer Output Results ===")
print(f"Outlet fluid: {Mixer.ex.fluid}")
print(f"Outlet pressure: {Mixer.ex.p:.2f} Pa")
print(f"Outlet mass flow rate: {Mixer.ex.m_dot:.2f} kg/s")
print(f"Outlet enthalpy: {Mixer.ex.h:.2f} J/kg")

try:
    T_out = PropsSI('T', 'P', Mixer.ex.p, 'H', Mixer.ex.h, Mixer.ex.fluid)
    print(f"Outlet temperature: {T_out - 273.15:.2f} °C")
except:
    print("Could not calculate outlet temperature.")