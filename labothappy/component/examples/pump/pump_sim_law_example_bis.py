from labothappy.component.pump.pump_similarity_laws_bis import PumpSimilarityLaws
import CoolProp.CoolProp as CP
import numpy as np

# Example characteristic curves and parameters
# V_dot_curve = np.array([20, 30, 40, 50, 60, 70, 80])   # m3/h
# Delta_H_curve = np.array([57, 55, 52, 49, 45, 42, 36])  # m (head falls with flow)
# eta_is_curve = np.array([0.45, 0.59, 0.69, 0.75, 0.79, 0.79, 0.75])  # eff peaks near mid-flow
# NPSH_curve = np.array([1.1, 1.1, 1.4, 1.8, 2, 3, 4.7])  # m, increases again near max flow
# W_dot_curve = np.array([7100, 8000, 8800, 9100, 9800, 10100, 10900])  # W

V_dot_curve = np.array([0, 0.5, 1, 1.5, 2, 3, 3.3])   # m3/h
Delta_H_curve = np.array([66, 56, 46, 36, 26, 6, 0])  # m (head falls with flow)

W_dot_hyd = (CP.PropsSI("D", "T", 293.15, "P", 101325, "Water")) * 9.81 * Delta_H_curve * (V_dot_curve / 3600)  # W
eta_is_curve = W_dot_hyd / np.array([1030, 900, 830, 760, 620, 550, 500])  # efficiency based on hydraulic and input power
# eta_is_curve = np.array([0.45, 0.59, 0.69, 0.75, 0.79, 0.79, 0.75])  # eff peaks near mid-flow
NPSH_curve = np.array([0, 3, 3, 3.5, 5, 6.8, 0])  # m, increases again near max flow
N_rated = 2900 # RPM
fluid_curve = "Water"
density_curve = CP.PropsSI("D", "T", 293.15, "P", 101325, fluid_curve)  # Density of water at 20°C and 1 atm


# If you have eta_curve -> you can find back W_dot_curve or vice versa

# Reference point: water at 20°C and 1 atm for a rated speed of 

PUMP = PumpSimilarityLaws()

# Set Inputs
# PUMP.set_inputs(
#     fluid = 'Water',
#     T_su = 20 + 273.15,  # K
#     P_su = 1e5,          # Pa
#     P_ex = 3e5,          # Pa
#     N_rot = 2000,        # rpm
# )

PUMP.set_inputs(
    P_su=1.3e5,  # Suction pressure in Pascals
    T_su=275.15+15,  # Suction temperature in Kelvin
    # P_ex=2.4e5,  # Exhaust pressure in Pascals
    N_rot=1589,  # Rotational speed in RPM
    m_dot = 0.36,  # Mass flow rate in kg/s
    fluid="R1233zd(E)",  # Actual fluid type
)


# Set Parameters
PUMP.set_parameters(
    V_dot_curve = V_dot_curve,
    Delta_H_curve = Delta_H_curve,
    eta_is_curve = eta_is_curve,
    NPSH_curve = NPSH_curve,
    # W_dot_curve = W_dot_curve,
    N_rot_rated = N_rated,
    # fluid_curve = fluid_curve,
    # density_curve = density_curve,
    mode = "M_N",  # Mode can be "M_N", "P_M", or "P_N"
)

PUMP.solve()
PUMP.print_results()

# curves_head = [(0, 66), (0.5, 56), (1, 46), (1.5, 36), (2, 26), (3, 6), (3.3, 0)] # ([m^3/h], [m])
# curves_power = [(0, 1030), (1, 900), (1.5, 830), (2, 760), (3, 620), (3.5, 550)]  # ([m^3/h], [W])
# curves_fluid = "Water"
# speed_ref = 2900 # RPM
# # curves_rho = 997  # Density of the curve fluid at reference conditions (kg/m^3)
# curves_rho = CP.PropsSI("D", "T", 293.15, "P", 101325, "Water")  # Density of water at 20°C and 1 atm

# # Create the model
# model = PumpSimilarityLaws()

# # Set the parameters for the model
# model.set_parameters(
#     curves_head=curves_head,
#     curves_power=curves_power,
#     curves_fluid=curves_fluid,
#     speed_ref=speed_ref,
#     curves_rho=curves_rho,
#     mode="P_M",  # Mode can be "M_N", "P_M", or "D_N"
# )

# # Set the inputs for the simulation
# model.set_inputs(
#     P_su=1.3e5,  # Suction pressure in Pascals
#     T_su=275.15+15,  # Suction temperature in Kelvin
#     P_ex=2.4e5,  # Exhaust pressure in Pascals
#     # N_rot=1100,  # Rotational speed in RPM
#     m_dot = 0.36,  # Mass flow rate in kg/s
#     fluid="R1233zd(E)",  # Actual fluid type
# )

# # Solve the system
# model.solve()
# model.print_results()

# # Plot the characteristic curves for the given speeds and flow range
# model.plot_characteristic_curves(
#     speeds_to_plot=[1450, 1750, 2900, 3500], flow_range=(0, 3.3), n_points=50
# )
