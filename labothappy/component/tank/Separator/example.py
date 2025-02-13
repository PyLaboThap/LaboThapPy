# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:43:39 2023

@author: Samuel Gendebien
"""

import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI

# from labothappy.connector.mass_connector import MassConnector
from component.tank.Separator.LV_Separator_SG import LV_Separator_SG

"-----------------------------------------------------------  TEST   ----------------------------------------------------------------"

LV_Separator_SG = LV_Separator_SG()

# Inputs
LV_Separator_SG.set_inputs(
                  su_fluid = 'R22',
                  su_T=350,
                  su_p = 100000,
                  su_m_dot = 14,
                  
                  )

# Params
LV_Separator_SG.set_parameters()

# Solve
LV_Separator_SG.solve()
LV_Separator_SG.print_results()
LV_Separator_SG.print_states_connectors()
