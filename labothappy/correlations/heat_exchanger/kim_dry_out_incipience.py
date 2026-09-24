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
Created on Mon Sep 13 11:21:49 2021

@author: jvega
"""


def kim_dry_out_incipience(G, q, Dh, P_star, rho_l, rho_v, mu_l, sigma, i_fg):
    """
    Inputs
    ------
    ?
    
    Outputs
    -------
    ?
    
    Reference
    ---------
    ?
    
    """
    
    Bo = q/G/i_fg
    We_lo = (Dh*G**2)/(rho_l*sigma)
    Ca = (mu_l*G)/(rho_l*sigma)
    x_di = 1.4*(We_lo**0.03)*(P_star**0.08) - 15*(Bo**0.15)*(Ca**0.35)*(rho_v/rho_l)**0.06
    return x_di