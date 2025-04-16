# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:39:37 2023

@author: Elise
"""
from labothappy.component.volumetric_machine.expander.expander_semi_empirical import ExpanderSE

import numpy as np

"Example of a semi-empirical expander component"
# Inputs: N_rot, T_amb, P_su, h_su, P_ex, fluid 
#         OR m_dot, T_amb, P_su, h_su, P_ex, fluid
# Parameters: AU_amb, AU_su_n, AU_ex_n, d_su1, m_dot_n, A_leak, W_dot_loss_0, alpha, C_loss, 
# rv_in, V_s, mode

#-------------------------------------------------------------------------------------------#
"First case: m_dot known and N_rot unknown"
#-------------------------------------------------------------------------------------------#
# # Create an instance of the expander component
# expander = ExpanderSE()
# "1. Inputs set through connectors"
# # Set properties for su connector
# expander.su.set_properties(P=4*1e5, T=273.15+70, fluid='R1233zd(E)')
# expander.ex.set_properties(P=1.1*1e5)

# # Set properties for ex connector
# expander.ex.set_properties(P=1.1*1e5)

# # Set rotational speed
# expander.W_mec.set_N(6000)

# # Set ambient temperature
# expander.Q_amb.set_T_cold(293)

# "2. Inputs set directly"
# # expander.set_inputs(
# #     N_rot=6000,
# #     T_amb=298.15,
# #     P_su=400000,
# #     h_su=485571,
# #     P_ex=110000,
# #     fluid='R1233zd(E)'  # Make sure to include fluid information
# # )

# expander.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
#             A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712,
#             mode = 'N_rot')

# # Check the setup
# expander.print_setup()

# # Solve the expander component
# expander.solve()

# # Print the results
# expander.print_results()

#-------------------------------------------------------------------------------------------#
"Second case: N_rot known and m_dot unknown"
#-------------------------------------------------------------------------------------------#
# Create an instance of the expander component
expander = ExpanderSE()
"1. Inputs set through connectors"
# # Set properties for su connector
# expander.su.set_fluid('R1233zd(E)')
# expander.su.set_p(4*1e5)
# expander.su.set_T(273.15+70)  # Sets h_su -> OK for inputs
# expander.su.set_m_dot(0.1)

# # Set properties for ex connector
# expander.ex.set_p(1.1*1e5)

# # Set ambient temperature
# expander.Q_amb.set_T_cold(293)

"2. Inputs set directly"
expander.set_inputs(
    m_dot=0.1,
    T_amb=298.15,
    P_su=400000,
    h_su=485571,
    P_ex=110000,
    fluid='R1233zd(E)'  # Make sure to include fluid information
)

expander.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
            A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712,
            mode = 'm_dot')

# Check the setup
expander.print_setup()

# Solve the expander component
expander.solve()

# Print the results
expander.print_results()




