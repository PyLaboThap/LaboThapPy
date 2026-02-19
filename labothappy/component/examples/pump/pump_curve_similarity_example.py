from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity
import CoolProp.CoolProp as CP
import numpy as np

# Example characteristic curves and parameters
V_dot_curve = np.array([20, 30, 40, 50, 60, 70, 80])   # m3/h
Delta_H_curve = np.array([57, 55, 52, 49, 45, 42, 36])  # m (head falls with flow)
eta_is_curve = np.array([0.45, 0.59, 0.69, 0.75, 0.79, 0.79, 0.75])  # eff peaks near mid-flow
NPSH_r_curve = np.array([1.1, 1.1, 1.4, 1.8, 2, 3, 4.7])  # m, increases again near max flow
N_rated = 2900 # RPM
# Reference point: water at 20°C and 1 atm for a rated speed of 2900 RPM

# If you have eta_curve -> you can find back W_dot_curve or vice versa

PUMP = PumpCurveSimilarity()

# Set Parameters
PUMP.set_parameters(
    V_dot_curve = V_dot_curve,
    Delta_H_curve = Delta_H_curve,
    eta_is_curve = eta_is_curve,
    NPSH_r_curve = NPSH_r_curve,
    N_rot_rated = N_rated,
    mode = "P_M",  # Mode can be "M_N", "P_M", or "P_N"
)

# Set Inputs
PUMP.set_inputs(
    P_su=1.3e5,  # Suction pressure in Pascals
    T_su=275.15+15,  # Suction temperature in Kelvin
    P_ex=2.4e5,  # Exhaust pressure in Pascals
    # N_rot=1589,  # Rotational speed in RPM
    m_dot = 0.36,  # Mass flow rate in kg/s
    fluid="R1233zd(E)",  # Actual fluid type
)

PUMP.solve()
PUMP.print_results()
rho_curve = CP.PropsSI("D", "T", 293.15, "P", 101325, "Water")  # Density of water at 20°C and 1 atm
PUMP.plot_characteristic_curves(
    speeds_to_plot=[1450, 1750, 2900, 3500], rho_curve=rho_curve
)
