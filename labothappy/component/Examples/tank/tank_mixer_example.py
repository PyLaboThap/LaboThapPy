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
    su_1_T = 50 + 273.15, # K
    su_1_p = 2*1e5, # Pa
    su_1_m_dot = 1, # kg/s
    su_1_fluid = 'Water',
    
    su_2_T = 100 + 273.15, # K
    su_2_p = 2*1e5, # Pa
    su_2_m_dot = 1, # kg/s
    su_2_fluid = 'Water'
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