import __init__
from component.turbomachinery.pump_extrapolation import PumpExtrapolationModel

# Instantiate Pump
Pump = PumpExtrapolationModel()

# Set Input
Pump.set_inputs(
    su_fluid="Cyclopentane",
    su_T=32 + 273.15,
    su_p=1 * 1e5,
    ex_p=31.5 * 1e5,
    Omega_pp=3000,
)

import numpy as np

"Parameters"
Pump.set_parameters(
    Omega_rated=0.3,
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
