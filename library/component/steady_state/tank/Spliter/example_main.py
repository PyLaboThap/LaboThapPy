# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from Port.Mass_connector import Mass_connector
from CoolProp.CoolProp import PropsSI
from Modules.Spliter_Geom import Spliter_Geom
from Spliter import Spliter

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Cyclopentane Su"
C5_su = Mass_connector()

# Set the fluid
C5_su.set_fluid("Cyclopentane")

# Set other properties
C5_su.set_m_dot(13.8)  # Example mass flow rate [kg/s]
C5_su.set_T(10 + 273.15)  # Example Temperature [K]

# C5_su.set_x(0.6) # Example quality [-]
C5_su.set_p(0.7*1e5)  # Example Pressure [Pa]

"Geometry"

Spliter_Geom = Spliter_Geom()
Spliter_Geom.set_parameters("DECAGONE_ACC_Spliter")

Spliter = Spliter(geom = Spliter_Geom)
Spliter.set_parameters(point_su = [C5_su])
 
Spliter.solve()

