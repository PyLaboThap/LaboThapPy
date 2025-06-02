# -*- coding: utf-8 -*-
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
