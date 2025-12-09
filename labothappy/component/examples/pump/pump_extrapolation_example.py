from labothappy.component.pump.pump_extrapolation import PumpExtrapolation
import numpy as np
# Instantiate Pump
Pump = PumpExtrapolation()

V_dot_curve = np.array([0.0, 1.0, 2.5, 4.0, 6.0, 8.0, 10.0])   # m3/h
Delta_H_curve = np.array([40.0, 38.0, 34.0, 28.0, 20.0, 12.0, 5.0])  # m (head falls with flow)
eta_is_curve = np.array([0.40, 0.68, 0.82, 0.88, 0.84, 0.74, 0.60])  # isentropic eff peaks near mid-flow
NPSH_r_curve = np.array([3.5, 2.8, 2.0, 1.6, 1.8, 2.8, 4.5])  # m, increases again near max flow

# Example rated numbers (consistent units)
N_rot_rated = 3000     # rpm
N_rot_actual = 3000    # rpm (set in inputs)
W_dot_el_rated = 5000  # W (example)
eta_m = 0.98           # mechanical efficiency of pump (shaft)
eta_max_motor = 0.95   # max electrical motor efficiency (fraction)
# Set Input
Pump.set_inputs(
    fluid="Cyclopentane",
    T_su=32 + 273.15, # K
    P_su=1 * 1e5, # Pa
    P_ex=31.5 * 1e5, # Pa
    N_rot=3000, # rpm
)

# import numpy as np

"Parameters"
Pump.set_parameters(
    N_rot_rated=N_rot_rated,
    min_flowrate=1,
    rated_flowrate=5,
    max_flowrate=10,
    PI_rated=5,
    D_p=4,
    V_dot_curve=V_dot_curve,
    Delta_H_curve=Delta_H_curve,
    eta_is_curve=eta_is_curve,
    NPSH_r_curve=NPSH_r_curve,
    eta_m=0.98,
    eta_max_motor=0.95,
    W_dot_el_rated=W_dot_el_rated,
)

Pump.solve()
print(Pump.W_dot_wf)
