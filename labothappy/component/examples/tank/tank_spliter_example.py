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

