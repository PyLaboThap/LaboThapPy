"""
Example of use of the moving-boundary plate heat exchanger model (HexMBPlate),
"""

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