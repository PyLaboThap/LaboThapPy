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
Created on Thu Jul 10 15:21:37 2025

@author: Basile
"""

import numpy as np

def balje_binsley(Re, t_TE, chord, t_blade, xhi1, xhi2):
    """
    Balje-Binsley : Profile pressure losses for axial turbines
    
    Inputs
    ------
    Re : Reynolds number [-]
    t_TE : Trailing Edge thickness [m]
    chord : Blade chord [m]
    t_blade : Blade thickness [m]
    xhi2 : Outlet stagger angle [rad]
    
    Output(s)
    ---------
    Yp : Profile pressure loss coefficient
    
    DP_loss = Y * 1/2 * v2**2 * rho_out
    
    References
    ----------    
    DESIGN METHODOLOGY OF SUPERCRITICAL CO2 BRAYTON CYCLE TURBOMACHINERIES
    Jekyoung Lee, Jeong Ik Lee, Yoonhan Ahn, Hojoon Yoon
    2012    
    """
    
    xhi2 = abs(xhi2 + xhi1)/2
    
    H_TE = 1.4 + 300/Re**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
    theta = 0.036*chord/Re**0.2 # Boundary layer momentum thickness : Empirical equation for turbulent plate
    theta_star = theta/(t_blade*np.sin(xhi2))
    
    A = 1-(1+H_TE)*theta_star-t_TE/t_blade # 1-(1+H_TE)*theta-t_TE/t_blade
    B = 1-H_TE*theta_star-t_TE/t_blade
    
    num_Yp = (np.cos(xhi2)**2 * A**2) / B**2 + (np.sin(xhi2)**2) * B**2
    den_Yp = 1 + 2 * (np.sin(xhi2)**2) * (B**2 - A)
    Yp = abs(1- num_Yp/den_Yp)
    
    # print(f"H_TE_2 : {H_TE}")
    # print(f"theta_2 : {theta_star}")
    # print(f"t_TE : {t_TE}")
    # print(f"t_blade : {t_blade}")
    # print(f"A_2 : {A}")
    # print(f"B_2 : {B}")
    # print(f"stagger_2 : {xhi2}")
    
    return Yp

def kacker_okaapu(AR, solidity, alpha1, alpha2, beta1):
    """
    Kacker-Okaapu : Secondary pressure losses for axial turbines
    
    Inputs
    ------
    AR : Aspect Ratio [-]
    solidity : Solidity [-]
    alpha1 : Absolute (if stator) or relative (then beta1, if rotor) inlet flow angle
    alpha2 : Absolute (if stator) or relative (then beta1, if rotor) outlet flow angle
    beta1 : Relative (if stator) or absolute (then alpha1, if rotor) inlet flow angle

    Output(s)
    ---------
    Ys : Secondary pressure loss coefficient
    
    DP_loss = Y * 1/2 * v2**2 * rho_out
    
    References
    ----------    
    DESIGN METHODOLOGY OF SUPERCRITICAL CO2 BRAYTON CYCLE TURBOMACHINERIES
    Jekyoung Lee, Jeong Ik Lee, Yoonhan Ahn, Hojoon Yoon
    2012    
    """
    
    # 5.2) Kacker-Okaapu : Secondary pressure losses
    Z = solidity*(alpha1-alpha2)/np.cos(alpha2) # Loading Factor
    Ys = abs(0.0334*1/AR*(np.cos(alpha2)/np.cos(beta1))*Z)
    
    return Ys
