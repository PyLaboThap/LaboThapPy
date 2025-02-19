# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:42:00 2024

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from Modules.LV_sep_Geom import LV_sep_Geom
from LV_separator import LV_Separator

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Cyclopentane Su"
Separator = LV_Separator()

# Inputs
Separator.set_inputs(
                  su_fluid = 'R22',
                  h_su=223873,
                  p_su = 100000,
                  m_dot_su = 14,                  
                  )

Separator.solve()
