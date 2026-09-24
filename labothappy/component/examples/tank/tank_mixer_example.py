# -*- coding: utf-8 -*-
# LaboThApPy: An open-source tool for advanced thermodynamic cycle simulation, designed for researchers and engineers.
#
# Copyright (C) <2025> <Université catholique de Louvain (UCLouvain), Belgique
#                       Université de Liège (ULiège), Belgique
#                       Université de Mons (UMONS), Belgique>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# List of the contributors to the development of LaboThApPy: see AUTHORS file.
# Description and complete License: see NOTICE & LICENSE files.

"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from labothappy.component.tank.tank_mixer import TankMixer 
from CoolProp.CoolProp import PropsSI


"--------- 1) Data ------------------------------------------------------------------------------------------"

"Instanciation"

Mixer = TankMixer(n_inlets = 2)
 
Mixer.set_inputs(
    T_su_1 = 50 + 273.15, # K
    P_su_1 = 2*1e5, # Pa
    m_dot_su_1 = 1, # kg/s
    fluid_su_1 = 'Water',
    
    T_su_2 = 100 + 273.15, # K
    P_su_2 = 2*1e5, # Pa
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
    
fig = Mixer.plot_thermo_states()
fig.show()

