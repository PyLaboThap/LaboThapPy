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

import heat_losses
import numpy as np


def Q_dot_abs(self, k):

        # Sun absorbed power
        eta_opt = self.params['eta_opt']*self.params['alpha_r']*self.params['refl_m']*self.params['envel_tau']
        Q_dot_sun_raw = self.DNI*np.cos(self.Theta)*self.params['W']# *eta_opt # self.L_disc
        Q_dot_sun = Q_dot_sun_raw*eta_opt
        Opt_losses = Q_dot_sun_raw*(1-eta_opt)

        # Heat loss power
        HL = self.heat_losses(k)
        Q_dot_loss = HL #*self.L_disc
        Q_dot_abs = Q_dot_sun - Q_dot_loss

        return Q_dot_abs
