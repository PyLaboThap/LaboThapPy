# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:43:39 2023

@author: Samuel Gendebien
"""

import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI

# from labothappy.connector.mass_connector import MassConnector
from component.tank.LV_Separator import LV_Separator

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

LV_Separator = LV_Separator()

# Inputs
LV_Separator.set_inputs(
                  su_fluid = 'R22',
                  T_su=350,
                  p_su = 100000,
                  m_dot_su = 14,
                  
                  )

# Params
LV_Separator.set_parameters()

# Solve
LV_Separator.solve()
LV_Separator.print_results()
LV_Separator.print_states_connectors()
