from labothappy.component.pump.pump_similarity_laws import PumpSimilarityLaws
import CoolProp.CoolProp as CP

# Example characteristic curves and parameters
curves_head = [(0, 66), (0.5, 56), (1, 46), (1.5, 36), (2, 26), (3, 6), (3.3, 0)] # ([m^3/h], [m])
curves_power = [(0, 1030), (1, 900), (1.5, 830), (2, 760), (3, 620), (3.5, 550)]  # ([m^3/h], [W])
curves_fluid = "Water"
speed_ref = 2900 # RPM
# curves_rho = 997  # Density of the curve fluid at reference conditions (kg/m^3)
curves_rho = CP.PropsSI("D", "T", 293.15, "P", 101325, "Water")  # Density of water at 20°C and 1 atm

# Create the model
model = PumpSimilarityLaws()

# Set the parameters for the model
model.set_parameters(
    curves_head=curves_head,
    curves_power=curves_power,
    curves_fluid=curves_fluid,
    speed_ref=speed_ref,
    curves_rho=curves_rho,
    mode="P_M",  # Mode can be "M_N", "P_M", or "D_N"
)

# Set the inputs for the simulation
model.set_inputs(
    P_su=1.3e5,  # Suction pressure in Pascals
    T_su=275.15+15,  # Suction temperature in Kelvin
    P_ex=2.4e5,  # Exhaust pressure in Pascals
    # N_rot=1100,  # Rotational speed in RPM
    m_dot = 0.36,  # Mass flow rate in kg/s
    fluid="R1233zd(E)",  # Actual fluid type
)

# Solve the system
model.solve()
model.print_results()

# Plot the characteristic curves for the given speeds and flow range
model.plot_characteristic_curves(
    speeds_to_plot=[1450, 1750, 2900, 3500], flow_range=(0, 3.3), n_points=50
)
