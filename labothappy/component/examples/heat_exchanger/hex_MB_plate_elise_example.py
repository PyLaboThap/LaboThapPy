"""
Example of use of the moving-boundary plate heat exchanger model (HexMBPlate),
"""

from labothappy.toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP
from labothappy.component.heat_exchanger.hex_MB_plate_elise import HexMBPlate


"HX instantiation"

HX = HexMBPlate('Plate')


case_study = "C5_COND"

if case_study == "Simple":
    # Simple case where pressure drops nad heat transfer coefficients are user-defined, 
    # and the geometry is a simple plate pack
    # The void fraction correlation
    # Condenser: cyclopentane condensed by water
    HX.set_inputs(
        fluid_H='Cyclopentane', T_su_H=139 + 273.15, P_su_H=0.77e5, m_dot_H=0.014,   # K, Pa, kg/s
        fluid_C='Water',        T_su_C=12 + 273.15,  P_su_C=5e5,    m_dot_C=0.2,
    )

    HX.set_parameters(flow_type='CounterFlow',          # Counter-flow: no LMTD correction (F = 1)
                      n_disc=50,                        # Number of regular cells (phase-change cells are added)
                      A_H=0.752, A_C=0.752,             # Heat transfer area [m^2]
                      V_C=0.9*1e-3, V_H=0.9*1e-3)       # Total volume of each side [m^3]

    HX.set_htc(htc_type='user_defined',
             htc_user_h={'liquid': 500, 'vapor': 200, 'two-phase': 2000, 'vapor-wet': 1000},
             htc_user_c={'liquid': 5000})

    HX.set_dp(dp_type='user_defined', dp_user_h=20e3, dp_user_c=50e3)

if case_study == "C5_COND":
    # Same as above, but with correlations for the heat transfer coefficients

    HX.set_inputs(
        fluid_H='Cyclopentane', T_su_H=139 + 273.15, P_su_H=0.77e5, m_dot_H=0.014,   # K, Pa, kg/s
        fluid_C='Water',        T_su_C=12 + 273.15,  P_su_C=5e5,    m_dot_C=0.2,
    )

    HX.set_parameters(flow_type='CounterFlow',          # Counter-flow: no LMTD correction (F = 1)
                      n_disc=50,                        # Number of regular cells (phase-change cells are added)
                      A_H=0.752, A_C=0.752,             # Heat transfer area [m^2]
                      V_C=0.9*1e-3, V_H=0.9*1e-3,       # Total volume of each side [m^3]
                      n_plates = 24,
                      n_channels_H = 12,
                      n_channels_C = 11,
                      L_plate = 0.324,
                      W_plate = 0.094,
                      L_eff = 0.2686,
                      plate_pitch = 0.004,
                      corrugation_amplitude = 0.002,
                      corrugation_pitch = 0.008,
                      chevron_angle = 0.4363,
                      k_wall = 45,
                      R_fouling = 0.0001)

    HX.set_htc(htc_type='correlation',
             htc_corr_h={'single-phase': 'Martin', 'two-phase': 'Longo_cond'},
             htc_corr_c={'single-phase': 'Martin', 'two-phase': 'Martin'})


    HX.set_dp(dp_type='user_defined', dp_user_h=20e3, dp_user_c=50e3)

if case_study == "C5_EVAP":
    # Evaporator: cyclopentane evaporated by a thermal oil
    HX.set_inputs(
        fluid_H='INCOMP::T66',  T_su_H=243 + 273.15, P_su_H=5e5,    m_dot_H=0.4,
        fluid_C='Cyclopentane', T_su_C=41 + 273.15,  P_su_C=31.5e5, m_dot_C=0.014,
    )
    
    HX.set_parameters(flow_type='CounterFlow',          # Counter-flow: no LMTD correction (F = 1)
                          n_disc=50,                        # Number of regular cells (phase-change cells are added)
                          A_H=0.752, A_C=0.752,             # Heat transfer area [m^2]
                          V_C=0.9*1e-3, V_H=0.9*1e-3,       # Total volume of each side [m^3]
                          n_plates = 24,
                          n_channels_H = 12,
                          n_channels_C = 11,
                          L_plate = 0.324,
                          W_plate = 0.094,
                          L_eff = 0.2686,
                          plate_pitch = 0.004,
                          corrugation_amplitude = 0.002,
                          corrugation_pitch = 0.008,
                          chevron_angle = 0.4363,
                          k_wall = 45,
                          R_fouling = 0.0001)
    
    HX.set_htc(htc_type='correlation',
                htc_corr_h={'single-phase': 'Martin', 'two-phase': 'Martin'},
                htc_corr_c={'single-phase': 'Martin', 'two-phase': 'Cooper'})

    HX.set_dp(dp_type='user_defined', dp_user_h=20e3, dp_user_c=50e3)

HX.solve()

print(f"Heat transfer rate      : {HX.Q_dot:10.1f} W")
print(f"Hot exhaust temperature : {HX.ex_H.T - 273.15:10.2f} °C")
print(f"Cold exhaust temperature: {HX.ex_C.T - 273.15:10.2f} °C")
print(f"Hot exhaust pressure    : {HX.ex_H.p:10.0f} Pa")
print(f"Cold exhaust pressure   : {HX.ex_C.p:10.0f} Pa")
print(f"Hot side charge         : {HX.charge_H * 1e3:10.2f} g")
print(f"Cold side charge        : {HX.charge_C * 1e3:10.2f} g")

HX.plot_cells()
HX.print_states_connectors()

# def plate_parameters(geom):
#     """Translate a PlateGeomSWEP geometry into the parameter names of HexMBPlate."""
#     params = dict(
#         A_c=geom.A_c,                           # Heat transfer area, cold side [m^2]
#         A_h=geom.A_h,                           # Heat transfer area, hot side [m^2]
#         pack_height=geom.h,                     # Height of the plate pack [m]
#         plate_width=geom.w,                     # Plate width [m]
#         plate_length=geom.l,                    # Plate length [m]
#         port_distance=geom.l_v,                 # Distance between the ports [m]
#         t_casing=geom.casing_t,                 # Casing thickness [m]
#         A_cs_c=geom.C_CS, A_cs_h=geom.H_CS,     # Cross-section of one channel [m^2]
#         D_h_c=geom.C_Dh, D_h_h=geom.H_Dh,       # Hydraulic diameter of one channel [m]
#         n_channels_c=geom.C_n_canals, n_channels_h=geom.H_n_canals,   # Number of channels [-]
#         t_channel_c=geom.C_canal_t, t_channel_h=geom.H_canal_t,       # Channel thickness [m]
#         V_c=geom.C_V_tot, V_h=geom.H_V_tot,     # Internal volume of each side [m^3]
#         n_plates=geom.n_plates,                 # Number of plates [-]
#         t_plate=geom.t_plates,                  # Plate thickness [m]
#         k_plate=geom.plate_cond,                # Plate thermal conductivity [W/(m K)]
#         corrugation_pitch=geom.plate_pitch_co,  # Corrugation pitch [m]
#         chevron_angle=geom.chevron_angle,       # Chevron angle [rad]
#         fouling=geom.fooling,                   # Fouling factor, hot side [-]
#     )
#     # Only some geometries define these (needed by the Martin-Holger, Amalfi and Shah correlations)
#     for name_geom, name_model in (('amplitude', 'amplitude'), ('w_v', 'port_width'), ('phi', 'enlargement_factor')):
#         if hasattr(geom, name_geom):
#             params[name_model] = getattr(geom, name_geom)
#     return params


# def fix_casing(geom, t_casing):
#     """Set the casing thickness and recompute the channel dimensions, with the same
#     formulas as geometry_plate_hx_swep.py."""
#     geom.casing_t = t_casing
#     for side, n_channels in (('C', geom.C_n_canals), ('H', geom.H_n_canals)):
#         t = ((geom.h - 2 * t_casing) - geom.n_plates * geom.t_plates) / (2 * n_channels)   # Channel thickness [m]
#         setattr(geom, f'{side}_canal_t', t)
#         setattr(geom, f'{side}_CS', t * (geom.w - 2 * t_casing))                         # Channel cross-section [m^2]
#         setattr(geom, f'{side}_Dh', 4 * t * geom.w / (2 * t + 2 * geom.w))               # Hydraulic diameter [m]
#     return geom


# "HX instantiation"

# HX = HexMBPlate('Plate')

# case_study = "C5_COND"

# if case_study == "C5_COND":
#     # Condenser: cyclopentane condensed by water
#     HX.set_inputs(
#         fluid_H='Cyclopentane', T_su_H=139 + 273.15, P_su_H=0.77e5, m_dot_H=0.014,   # K, Pa, kg/s
#         fluid_C='Water',        T_su_C=12 + 273.15,  P_su_C=5e5,    m_dot_C=0.2,
#     )
#     geom = PlateGeomSWEP()
#     geom.set_parameters("B20Hx24/1P")

#     htc_corr_h = {'single_phase': 'Gnielinski', 'two-phase': 'Han_cond_BPHEX'}
#     htc_corr_c = {'single_phase': 'Gnielinski', 'two-phase': 'Han_Boiling_BPHEX_HTC'}
#     dp_user_h, dp_user_c = 20e3, 50e3     # Pa

# elif case_study == "C5_EVAP":
#     # Evaporator: cyclopentane evaporated by a thermal oil
#     HX.set_inputs(
#         fluid_H='INCOMP::T66',  T_su_H=243 + 273.15, P_su_H=5e5,    m_dot_H=0.4,
#         fluid_C='Cyclopentane', T_su_C=41 + 273.15,  P_su_C=31.5e5, m_dot_C=0.014,
#     )
#     geom = PlateGeomSWEP()
#     geom.set_parameters("B20Hx24/1P")

#     # The old example used 'Boiling_curve' (pool boiling on a tube, needs 'Tube_OD'),
#     # which does not apply to plates: a plate boiling correlation is used instead.
#     htc_corr_h = {'single_phase': 'Gnielinski', 'two-phase': 'Han_cond_BPHEX'}
#     htc_corr_c = {'single_phase': 'Gnielinski', 'two-phase': 'Han_Boiling_BPHEX_HTC'}
#     dp_user_h, dp_user_c = 20e3, 50e3

# elif case_study == "C5_DSH":
#     # Desuperheater: superheated cyclopentane vapour cooled by water (see note above)
#     HX.set_inputs(
#         fluid_H='Cyclopentane', T_su_H=205 + 273.15, P_su_H=1e5, m_dot_H=0.014,
#         fluid_C='Water',        T_su_C=12 + 273.15,  P_su_C=4e5, m_dot_C=0.2,
#     )
#     geom = PlateGeomSWEP()
#     geom.set_parameters("B35TM0x10/1P")
#     fix_casing(geom, 0.005)   # The stored 0.05 m casing is thicker than the plate pack (see note above)

#     htc_corr_h = {'single_phase': 'Gnielinski', 'two-phase': 'Han_cond_BPHEX'}
#     htc_corr_c = {'single_phase': 'Gnielinski', 'two-phase': 'Han_Boiling_BPHEX_HTC'}
#     dp_user_h, dp_user_c = 0.0, 0.0

# elif case_study == "R245fa_COND":
#     # Condenser: R245fa condensed by water
#     HX.set_inputs(
#         fluid_H='R245fa',     T_su_H=334.38, P_su_H=2.1e5, m_dot_H=0.4197,
#         fluid_C='Water',      T_su_C=15 + 273.15, P_su_C=2e5, m_dot_C=2.6,
#     )
#     geom = PlateGeomSWEP()
#     geom.set_parameters("P200THx140/1P_Condenser")

#     htc_corr_h = {'single_phase': 'martin_holger_plate_HTC', 'two-phase': 'Han_cond_BPHEX'}
#     htc_corr_c = {'single_phase': 'water_plate_HTC', 'two-phase': 'Han_Boiling_BPHEX_HTC'}
#     dp_user_h, dp_user_c = 20e3, 50e3

# "Parameters"

# HX.set_parameters(flow_type='CounterFlow',          # Counter-flow: no LMTD correction (F = 1)
#                   n_disc=50,                        # Number of regular cells (phase-change cells are added)
#                   void_fraction_model='Zivi',
#                   A_h=geom.A_h, A_c=geom.A_c,
#                   V_c=geom.C_V_tot, V_h=geom.H_V_tot)       # Void fraction model for the charge

# "Heat transfer and pressure drops"

# # HX.set_htc(htc_type='correlation', htc_corr_h=htc_corr_h, htc_corr_c=htc_corr_c)

# # Other possibilities:
# HX.set_htc(htc_type='user_defined',
#              htc_user_h={'liquid': 500, 'vapor': 200, 'two-phase': 2000, 'vapor-wet': 1000},
#              htc_user_c={'liquid': 5000})
# #   HX.set_dp()                                        # No pressure drops
# #   HX.set_dp(dp_type='correlation_disc', dp_corr_h=..., dp_corr_c=...)

# HX.set_dp(dp_type='user_defined', dp_user_h=dp_user_h, dp_user_c=dp_user_c)

# "Solve"

# HX.solve()

# print(f"Heat transfer rate      : {HX.Q_dot:10.1f} W")
# print(f"Hot exhaust temperature : {HX.ex_H.T - 273.15:10.2f} °C")
# print(f"Cold exhaust temperature: {HX.ex_C.T - 273.15:10.2f} °C")
# print(f"Hot exhaust pressure    : {HX.ex_H.p:10.0f} Pa")
# print(f"Cold exhaust pressure   : {HX.ex_C.p:10.0f} Pa")
# print(f"Hot side charge         : {HX.charge_h * 1e3:10.2f} g")
# print(f"Cold side charge        : {HX.charge_c * 1e3:10.2f} g")

# HX.plot_cells()



# # /!\ PRESSURE DROP CORRELATION NOT WORKING!! R1233ZDE NOT WORKING BC OF VISCOSITY PROPERTY!!!