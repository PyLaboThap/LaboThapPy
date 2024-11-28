"""
Created on Aug 03 21:31:37 2023

@author: Elise
"""
import __init__
from component.steady_state.volumetric_machine.compressor.constant_isentropic_efficiency.simulation_model import CompressorCstEff

import numpy as np

# Example usage
CP = CompressorCstEff()
CP.print_setup()

# "If the inputs are not set directly BUT through the connectors"
CP.su.set_properties(P=319296.5575177148, T=331.033964665788, fluid='R1233ZDE')
CP.ex.set_properties(P=606240.1433176235, fluid='R1233ZDE')
CP.set_parameters(eta_is=0.8)
CP.print_setup()

CP.solve()
CP.print_results()
CP.plot_component_comp_cst_eff()
CP.plot_connectors_comp_cst_eff()
