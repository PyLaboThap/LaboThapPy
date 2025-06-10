# -*- coding: utf-8 -*-
"""
Created on Wed June 04 14:42:00 2024

@author: Elise
"""

from labothappy.component.valve.isenthalpic_valve_P_ex import IsenthalpicValve_P_ex 
from CoolProp.CoolProp import PropsSI


"--------- 1) Data ------------------------------------------------------------------------------------------"

"Instanciation"

Valve = IsenthalpicValve_P_ex()

Valve.set_inputs(
    T_su = 50 + 273.15, # K
    P_su = 3*1e5, # Pa
    P_ex = 1e5, # kg/s
    fluid = 'R1233zd(E)',

    )

"--------- 2) Solve ------------------------------------------------------------------------------------------"
Valve.solve()


"--------- 3) Results ------------------------------------------------------------------------------------------"
Valve.print_results()