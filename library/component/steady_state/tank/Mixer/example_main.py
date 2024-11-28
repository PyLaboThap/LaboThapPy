# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from Port.Mass_connector import Mass_connector
from CoolProp.CoolProp import PropsSI
from Modules.Mixer_Geom import Mixer_Geom
from Mixer import Mixer

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

"Cyclopentane Su 2"

C5_su_2 = Mass_connector()

# Set the fluid
C5_su_2.set_fluid("Cyclopentane")

# Set other properties
C5_su_2.set_m_dot(13.8)  # Example mass flow rate [kg/s]
C5_su_2.set_T(20 + 273.15)  # Example Temperature [K]

# C5_su.set_x(0.6) # Example quality [-]
C5_su_2.set_p(0.7*1e5)  # Example Pressure [Pa]

"Geometry"

Mixer_geom = Mixer_Geom()
Mixer_geom.set_parameters("DECAGONE_ACC_Mixer")

Mixer = Mixer(geom = Mixer_geom)
Mixer.set_parameters(point_su = [C5_su, C5_su_2])
 
Mixer.solve()
