# -*- coding: utf-8 -*-
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