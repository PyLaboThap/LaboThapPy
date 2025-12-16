# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:43:39 2023

@author: Samuel Gendebien
"""

# import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.component.tank.tank_LV_separator import TankLVSeparator

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

LV_Separator = TankLVSeparator()

# Inputs
LV_Separator.set_inputs(
                  fluid = 'R22',
                  x_su = 0.5,
                  P_su = 100000,
                  m_dot = 14,
                  )

# Params
LV_Separator.set_parameters()

# Solve
LV_Separator.solve()
LV_Separator.print_results()
LV_Separator.print_states_connectors()
