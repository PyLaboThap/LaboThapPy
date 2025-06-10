# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:43:39 2023

@author: Samuel Gendebien
"""

import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI

from component.tank.tank_LV_separator import LV_Separator

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

LV_Separator = LV_Separator()

# Inputs
LV_Separator.set_inputs(
                  su_fluid = 'R22',
                  x_su = 0.5,
                  p_su = 100000,
                  m_dot_su = 14,
                  )

# Params
LV_Separator.set_parameters()

# Solve
LV_Separator.solve()
LV_Separator.print_results()
LV_Separator.print_states_connectors()
