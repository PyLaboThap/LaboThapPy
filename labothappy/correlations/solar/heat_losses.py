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
Created on Tue Jul 30 14:32:39 2024

@author: Basile

source : 
Semi-empirical correlation to model heat losses 
along solar parabolic trough collectors 
Rémi Dickes, Vincent Lemort and Sylvain Quoilin 
"""

import numpy as np


def heat_losses(self, k):
    "Calibrated for soponova_microcsp collector"        

    T_amb = self.Q_amb.T_cold - 273.15 # °C
    T_htf = self.T[k] - 273.15 # °C

    Terms = []
    Terms.append(self.params['a'][0]) 
    Terms.append(self.params['a'][1]*(T_htf - T_amb))
    Terms.append(self.params['a'][2]*(T_htf - T_amb)**2) 
    Terms.append(self.DNI*np.cos(self.Theta)*(self.params['a'][3]*T_htf**2))
    Terms.append(self.DNI*np.cos(self.Theta)*(self.params['a'][4]*np.sqrt(self.v_wind)))
    Terms.append(self.params['a'][5]*T_htf**3)
    Terms.append(self.v_wind*self.params['a'][6])
    Terms.append(self.v_wind*self.params['a'][7]*(T_htf - T_amb))
    Terms.append(np.sqrt(self.v_wind)*self.params['a'][8])
    Terms.append(np.sqrt(self.v_wind)*self.params['a'][9]*(T_htf - T_amb))

    Terms_np = np.array(Terms) # W/m (of line collector)

    return sum(Terms_np)