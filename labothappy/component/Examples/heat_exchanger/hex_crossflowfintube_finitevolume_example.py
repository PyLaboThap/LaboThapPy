import __init__ 

import matplotlib.pyplot as plt
from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from toolbox.geometries.heat_exchanger.geometry_tube_and_fins_hx import CFTubeAndFinsGeom
from component.heat_exchanger.hex_crossflowfintube_finitevolume import CrossFlowTubeAndFinsHTX

"--------- 1) Data ------------------------------------------------------------------------------------------"

"Cyclopentane Su"
ACC = CrossFlowTubeAndFinsHTX()

T_in_air = 100 + 273.15
P_in_air = 101325

rho_in_air = PropsSI('D', 'T',T_in_air, 'P', P_in_air, 'Air')

V_dot = 256.4/2 # m^3/s
m_dot_air = V_dot*rho_in_air

# Evaporator Case
ACC.set_inputs(
    # First fluid
    Hsu_fluid = "Cyclopentane",
    Hsu_T = 53 + 273.15, # K
    Hsu_p = 116645.96, # Pa
    # Hsu_x = 0.5, # Quality 
    Hsu_m_dot = 13.8/2, # kg/s

    # Second fluid
    Csu_fluid = "Air",
    Csu_T = T_in_air, # K
    Csu_p = P_in_air, # Pa
    Csu_m_dot = m_dot_air, # kg/s  # Make sure to include fluid information
)

"Geometry"
HX_geom = CFTubeAndFinsGeom()
HX_geom.set_parameters("DECAGONE_ACC") # DECAGONE_ACC_1Bundle # DECAGONE_ACC

ACC.set_parameters(
    A_finned = HX_geom.A_finned, A_flow = HX_geom.A_flow, A_in_tot = HX_geom.A_in_tot, A_out_tot = HX_geom.A_out_tot, A_unfinned = HX_geom.A_unfinned, 
    B_V_tot = HX_geom.B_V_tot, Fin_OD = HX_geom.Fin_OD, Fin_per_m = HX_geom.Fin_per_m, Fin_t = HX_geom.Fin_t, Fin_type = HX_geom.Fin_type,
    Finned_tube_flag = HX_geom.Finned_tube_flag, fouling = HX_geom.fouling, h = HX_geom.h, k_fin = HX_geom.k_fin,
    L = HX_geom.L, Tube_pass = HX_geom.n_passes, n_rows = HX_geom.n_rows, n_tubes = HX_geom.n_tubes, pitch = HX_geom.pitch, pitch_ratio = HX_geom.pitch_ratio,
    T_V_tot = HX_geom.T_V_tot, tube_arrang = HX_geom.tube_arrang, Tube_cond = HX_geom.Tube_cond, Tube_L = HX_geom.Tube_L, Tube_OD = HX_geom.Tube_OD, 
    Tube_t = HX_geom.Tube_t, w = HX_geom.w, 
 
    Fin_Side = 'C', H_DP_ON = True, C_DP_ON = True, n_disc = 100)
 
ACC.solve()

# Plot the matrix using imshow for the tube side
plt.pcolor(ACC.T_matrix)
plt.colorbar()
plt.show()

#Plot the pressure drop evolution
import matplotlib.pyplot as plt
import numpy as np

# Slice the matrix
tube_side_mask = np.zeros_like(ACC.P_matrix, dtype=bool)
bank_side_mask = np.zeros_like(ACC.P_matrix, dtype=bool)

tube_side_mask[1::2, :] = True
bank_side_mask[0::2, :] = True

tube_side_values = np.ma.masked_where(~tube_side_mask, ACC.P_matrix)
bank_side_values = np.ma.masked_where(~bank_side_mask, ACC.P_matrix)

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Plot both matrices
tube_cmap = ax.imshow(tube_side_values, cmap='viridis', aspect='auto', origin='lower')
bank_cmap = ax.imshow(bank_side_values, cmap='plasma', aspect='auto', origin='lower')

# Add colorbars
cbar_tube = plt.colorbar(tube_cmap, ax=ax, orientation='vertical', pad=0.02)
cbar_tube.set_label('Tube Side Pressure [Pa]')

cbar_bank = plt.colorbar(bank_cmap, ax=ax, orientation='vertical', pad=0.12)
cbar_bank.set_label('Bank Side Pressure [Pa]')

# Labels
ax.set_title('Pressure Distribution: Tube Side vs Bank Side')
ax.set_xlabel('Disc Index')
ax.set_ylabel('Row Index')

plt.tight_layout()
plt.show()

ACC.print_setup()
print("=== Heat Exchanger outlet conditions ===")
print(f"ex bundle: fluid_bundle = {ACC.B_su.fluid}, T_ex_bundle ={(int(round(ACC.B_ex.T)))}[K], h_ex_bundle={(int(round(ACC.B_ex.h)))}[J/kg-K], p_ex_bundle={(int(round(ACC.B_ex.p)))}[Pa], m_dot_ex_bundle={(int(round(ACC.B_su.m_dot)))}[kg/s]")
print("======")
print(f"ex tube: fluid={ACC.T_su.fluid}, T_ex_tube ={(int(round(ACC.T_ex.T)))}[K], h_ex_tube={(int(round(ACC.T_ex.h)))}[J/kg-K], p_ex_tube={(int(round(ACC.T_ex.p)))}[Pa], m_dot_tube={(int(round(ACC.T_su.m_dot)))}[kg/s]")
