import __init__ 

import matplotlib.pyplot as plt
from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from modules.geometry_cross_flow_fins import GeometryCrossFlowFins
from simulation_model import CrossFlowTubeAndFinsHTX

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Cyclopentane Su"
ACC = CrossFlowTubeAndFinsHTX()

T_in_air = 12.5 + 273.15
P_in_air = 101325

rho_in_air = PropsSI('D', 'T',T_in_air, 'P', P_in_air, 'Air')

V_dot = 256.4/2 # m^3/s
m_dot_air = V_dot*rho_in_air

# Evaporator Case
ACC.set_inputs(
    # First fluid
    Hsu_fluid = "Cyclopentane",
    Hsu_T = 53.6 + 273.15, # K
    Hsu_p = 0.7*1e5, # Pa
    Hsu_m_dot = 13.8/2, # kg/s

    # Second fluid
    Csu_fluid = "Air",
    Csu_T = T_in_air, # K
    Csu_p = P_in_air, # Pa
    Csu_m_dot = m_dot_air, # kg/s  # Make sure to include fluid information
)

"Geometry"

HX_geom = GeometryCrossFlowFins()
HX_geom.set_parameters("DECAGONE_ACC_5Bundle") # DECAGONE_ACC_1Bundle # DECAGONE_ACC

ACC.set_parameters(
    A_finned = HX_geom.A_finned, A_flow = HX_geom.A_flow, A_in_tot = HX_geom.A_in_tot, A_out_tot = HX_geom.A_out_tot, A_unfinned = HX_geom.A_unfinned, 
    B_V_tot = HX_geom.B_V_tot, Fin_OD = HX_geom.Fin_OD, Fin_per_m = HX_geom.Fin_per_m, Fin_t = HX_geom.Fin_t, Fin_type = HX_geom.Fin_type,
    Finned_tube_flag = HX_geom.Finned_tube_flag, fouling = HX_geom.fouling, h = HX_geom.h, k_fin = HX_geom.k_fin,
    L = HX_geom.L, n_passes = HX_geom.n_passes, n_rows = HX_geom.n_rows, n_tubes = HX_geom.n_tubes, pitch = HX_geom.pitch, pitch_ratio = HX_geom.pitch_ratio,
    T_V_tot = HX_geom.T_V_tot, tube_arrang = HX_geom.tube_arrang, Tube_cond = HX_geom.Tube_cond, Tube_L = HX_geom.Tube_L, Tube_OD = HX_geom.Tube_OD, 
    Tube_t = HX_geom.Tube_t, w = HX_geom.w, 
 
    Fin_Side = 'C', H_DP_ON = True, C_DP_ON = True, n_disc = 100)
 
ACC.solve()

# Plot the matrix using imshow
plt.pcolor(ACC.T_matrix)

# Add a colorbar to indicate color values
plt.colorbar()

# Show the plot
plt.show()