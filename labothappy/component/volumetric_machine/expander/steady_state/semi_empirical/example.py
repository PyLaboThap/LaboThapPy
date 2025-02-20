# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:39:37 2023

@author: Elise
"""

import __init__
from simulation_model import ExpanderSE

import numpy as np

# Example usage
expander = ExpanderSE()

"If the inputs are not set directly BUT throught the connectors"
expander.su.set_fluid('R1233zd(E)')
expander.ex.set_fluid('R1233zd(E)')

# # Set properties for su connector
expander.su.set_p(4*1e5)
# expander.su.set_m_dot(0.166667)
expander.su.set_T(273.15+70)  # You need to set su.h appropriately

# # Set properties for ex connector
expander.ex.set_p(1.1*1e5)

# Set rotational speed
expander.W_exp.set_N(6000)

# Set ambient temperature
expander.Q_amb.set_T_cold(293)

# Setting inputs
# expander.set_inputs(
#     N_rot=3000,
#     T_amb=298.15,
#     su_p=1e5,
#     su_T=300,
#     ex_p=1e4,
#     su_fluid='R134a'  # Make sure to include fluid information
# )

expander.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
            A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)

expander.set_inputs(N_rot=6000, T_amb=293)

# Solve the expander component
expander.solve()

# expander.print_setup()
# expander.print_results()

# print(expander.defined)  # Should print True if the component was successfully solved


# m-su:


# Fluid: R1233zd(E)
# Mass flow rate: 0.449[kg/s]
# Volume flow rate: 151.86253311317586[m^3/h]
# Temperature: 315.74539914862277[K]
# Pressure: 200000.0[Pa]
# Spec. enthalpy: 465664.22392372205[J/kg]
# Spec. entropy: 1887.6846063307426[J/kg/K]
# Mass density: 10.643836678236987[kg/m^3]
# Quality: -1.0[-]


