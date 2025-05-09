"""
Created on Aug 03 21:31:37 2023

@author: Elise
"""
import __init__
from component.volumetric_machine.compressor.steady_state.constant_isentropic_efficiency.simulation_model import CompressorCstEff
from CoolProp.CoolProp import PropsSI
import numpy as np

# Example usage
CP = CompressorCstEff()
# CP.print_setup()

T_in_cp = 273.15 - 2
SH = 1
P_out_cp = 110*1e5

P_in_cp = PropsSI('P', 'T', T_in_cp- SH, 'Q', 0, 'CO2')

# "If the inputs are not set directly BUT through the connectors"
CP.su.set_properties(P=P_in_cp, T=T_in_cp, fluid='CO2', m_dot = 0.08)
CP.ex.set_properties(P=P_out_cp, fluid='CO2')
CP.set_parameters(eta_is=0.8)
CP.print_setup()

CP.solve()
CP.print_results()
# CP.plot_component_comp_cst_eff()
# CP.plot_connectors_comp_cst_eff()
