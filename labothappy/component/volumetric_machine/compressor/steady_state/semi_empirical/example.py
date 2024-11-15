import __init__
from component.steady_state.volumetric_machine.compressor.semi_empirical.simulation_model import CompressorSE

import numpy as np

# Example usage
CP = CompressorSE()
CP.print_setup()

"If the inputs are not set directly BUT throught the connectors"
CP.su.set_fluid('R1233ZDE')

# Set properties for su connector
CP.su.set_p(319296.5575177148)
CP.su.set_T(331.033964665788)  # You need to set su.h appropriately

# Set properties for ex connector
CP.ex.set_p(606240.1433176235)

# Set rotational speed
CP.W_cp.set_N(6000)

# Set ambient temperature
CP.Q_amb.set_T_cold(293)

# Setting inputs
# CP.set_inputs(
#     N_rot=6000,
#     T_amb=293,
#     su_p=319296.5575177148,
#     su_T=331.033964665788,
#     ex_p=606240.1433176235,
#     su_fluid='R1233ZDE'  # Make sure to include fluid information
# )

# Setting parameters
CP.set_parameters(
    AU_amb=9.96513290e+00, AU_su_n=1.02359773e+01, AU_ex_n=2.24133147e+00, d_ex=1.82304791e-02, m_dot_n=0.1, 
    A_leak=3.66336680e-07, W_dot_loss_0=9.05482168e-01, alpha=3.22395090e-03, C_loss=1.11169710e-061, rv_in=1.7, V_s=1.17889079e-04
)

CP.print_setup()
# Solve the expander component
CP.solve()
CP.print_results()

# print(expander.defined)  # Should print True if the component was successfully solved



