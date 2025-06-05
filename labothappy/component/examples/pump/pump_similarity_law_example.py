import __init__ 

from labothappy.component.pump.pump_similarity_laws import PumpSimilarityLaws    

# Example characteristic curves and parameters
curves_head = [(0, 66), (0.5, 56), (1, 46), (1.5, 36), (2, 26), (3, 6), (3.3, 0)]
curves_power = [(0, 1.03), (1, 0.9), (1.5, 0.83), (2, 0.76), (3, 0.62), (3.5, 0.55)]
curves_fluid = 'Water'
speed_ref = 2900

# Create the model
model = PumpSimilarityLaws()

# Set the parameters for the model
model.set_parameters(
    curves_head=curves_head,
    curves_power=curves_power,
    curves_fluid=curves_fluid,
    speed_ref=speed_ref
)

# Set the inputs for the simulation
model.set_inputs(
    P_su=2e5,  # Suction pressure in Pascals
    T_su=300,  # Suction temperature in Kelvin
    P_ex=6e5,  # Exhaust pressure in Pascals
    N_rot=2500,  # Rotational speed in RPM
    fluid='R1233zd(E)'  # Actual fluid type
)

# # Solve the system
model.solve()

# # Plot the characteristic curves for the given speeds and flow range
model.plot_characteristic_curves(speeds_to_plot=[1450, 1750, 2900, 3500], flow_range=(0, 3.3), n_points=50)
