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
Created on Tue Dec 19 14:43:39 2023

@author: Samuel Gendebien
"""

# import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.component.tank.tank_LV_separator import TankLVSeparator

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

LV_Separator = TankLVSeparator()

# Inputs
LV_Separator.set_inputs(
                  fluid = 'R22',
                  x_su = 0.5,
                  P_su = 100000,
                  m_dot = 14,
                  )

# Params
LV_Separator.set_parameters()

# Solve
LV_Separator.solve()
LV_Separator.print_results()
LV_Separator.print_states_connectors()

fig = LV_Separator.plot_thermo_states()
fig.show()

