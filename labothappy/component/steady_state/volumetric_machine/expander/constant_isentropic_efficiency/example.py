# -*- coding: utf-8 -*-
"""
Created on Aug 03 21:31:37 2023

@author: Elise
"""

import __init__
from component.steady_state.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff

import numpy as np

# Example usage
EXP = ExpanderCstEff()
EXP.print_setup()

# "If the inputs are not set directly BUT throught the connectors"
EXP.su.set_properties(P=955214.9, T=374.18, fluid='R134a')
EXP.ex.set_properties(P=293940.1)
EXP.set_parameters(eta_is=0.8)
EXP.print_setup()

EXP.solve()
EXP.print_results()
