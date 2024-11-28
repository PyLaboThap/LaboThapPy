# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from Port.Mass_connector import Mass_connector
from CoolProp.CoolProp import PropsSI
from Modules.LV_sep_Geom import LV_sep_Geom
from LV_separator import LV_separator

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Cyclopentane Su"
C5_su = Mass_connector()

# Set the fluid
C5_su.set_fluid("Cyclopentane")

# Set other properties
C5_su.set_m_dot(13.8)  # Example mass flow rate [kg/s]
C5_su.set_T(20 + 273.15)  # Example Temperature [K]

# C5_su.set_x(0.6) # Example quality [-]
C5_su.set_p(0.7*1e5)  # Example Pressure [Pa]

"Pressure Drop"

DP_H_ON = True
DP_C_ON = True

"Geometry"

ACC_sep_geom = LV_sep_Geom()
ACC_sep_geom.set_parameters("DECAGONE_ACC_5_inputs")

ACC_sep = LV_separator()
ACC_sep.set_parameters(point_su = [C5_su, ], geom = ACC_sep_geom)
 
ACC_sep.solve()
