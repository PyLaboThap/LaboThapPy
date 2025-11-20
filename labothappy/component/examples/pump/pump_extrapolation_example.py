import __init__
from component.pump.pump_extrapolation import PumpExtrapolationModel

# Instantiate Pump
Pump = PumpExtrapolationModel()

# Set Input
Pump.set_inputs(
    fluid="Cyclopentane",
    T_su=32 + 273.15,
    P_su=1 * 1e5,
    P_ex=31.5 * 1e5,
    N_rot=3000,
)

import numpy as np

"Parameters"
Pump.set_parameters(
    N_rot_rated=0.3,
    min_flowrate=1,
    rated_flowrate=5,
    max_flowrate=10,
    PI_rated=5,
    D_p=4,
    V_dot_curve=np.array([4, 4]),
    Delta_H_curve=np.array([4, 4]),
    eta_is_curve=np.array([0.85, 0.85]),
    NPSH_r_curve=np.array([5, 5]),
    eta_m=0.99,
    eta_max_motor=0.99,
    W_dot_el_rated=5,
)

Pump.solve()
print(Pump.W_dot_wf)
