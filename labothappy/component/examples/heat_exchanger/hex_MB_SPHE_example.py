# -*- coding: utf-8 -*-
"""
Spiral plate heat exchanger (SPHE) with HexMBChargeSensitive, geometry resolved by hex_geometry.py.
Cells are numbered from the core on both streams: the cold stream enters at the core, the hot one at the periphery.
"""
from labothappy.component.heat_exchanger.hex_MB import HexMB
from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.plot_SPHE import plot_spiral, plot_incidence, plot_contacts

HX = HexMB('SPHE')

case_study = "Condenser"

# -------------------------------------------------------------------------------------------------------------
# Correlation names (SPHE_HTC_dispatcher):
#   1P   : 'Picon_Nunez' (was "SPHE"), 'Rajavel' (water), 'Curved_rect_duct', 'Rect_duct'
#   2P   : hot  (condensation) 'Zakeralhoseini'
#          cold (evaporation)  'Zhao' (was "Zhao_SPHE"), 'Nariai'          (None / 'auto': Zhao)
# Pressure drop: not available yet for SPHE -> set_DP() (no pressure drop) until the DP bank exists.
# -------------------------------------------------------------------------------------------------------------

if case_study == "Condenser":
    HX.set_inputs(
        fluid_H='Propane', T_su_H=80 + 273.15, P_su_H=22e5, m_dot_H=0.02,
        fluid_C='Water', T_su_C=50 + 273.15, P_su_C=5e5, m_dot_C=0.3)

    Corr_H = {"1P": "Picon_Nunez", "2P": "Zakeralhoseini"}
    Corr_C = {"1P": "Rajavel", "2P": None}

    Corr_H_DP = {"1P" : "Picon_Nunez", "2P" : "Choi"}
    Corr_C_DP = {"1P" : "Picon_Nunez", "2P" : "Choi"}
    
    params = dict(b_C=0.005, b_H=0.005, D_ext=0.39, H=0.39, r0_in=1e-3, t=0.002,
                  inner_ch='cold', charge_compute = False, k_wall = 50)

    params = {'b_C' : 0.005, 
              'b_H' : 0.005, 
              'D_ext' : 0.4, 
              'H' : 0.4, 
              'r0_in' : 1e-3, 
              't' : 0.002,
              'inner_ch' : 'cold', 
              'charge_compute' : False}

elif case_study == "Evaporator":
    HX.set_inputs(
        fluid_H='Water', T_su_H=10 + 273.15, P_su_H=5e5, m_dot_H=2,
        fluid_C='Propane', h_su_C=353418.2920873372, P_su_C=5.4e5, m_dot_C=0.02)

    Corr_H = {"1P": "Rajavel", "2P": None}
    Corr_C = {"1P": "Picon_Nunez", "2P": "Zhao"}

    Corr_H_DP = {"1P" : "Picon_Nunez", "2P" : "Choi"}
    Corr_C_DP = {"1P" : "Picon_Nunez", "2P" : "Choi"}
    
    params = dict(b_C=0.00143, b_H=0.00174, D_ext=0.18, H=0.1, r0_in=2.5e-2, t=0.002,
                  inner_ch='hot', charge_compute = False, k_wall = 50)

# Geometry: give ONE of D_ext (outer diameter), n_turns (number of turns) or A_sp (exchange area of one spiral).
# Optional: k_wall (wall conduction), R_f_H / R_f_C [m2 K/W], n_parallel. (n_series > 1 is not available.)
HX.set_parameters(**params, Flow_Type="CounterFlow", n_disc=100)


# UD_H_HTC = {'Liquid': 5000,
#             'Vapor' : 2000,
#             'Two-Phase' : 10000,
#             'Vapor-wet' : 10000,
#             'Dryout' : 10000,
#             'Transcritical' : 5000}

# UD_C_HTC = {'Liquid': 5000,
#             'Vapor' : 2000,
#             'Two-Phase' : 10000,
#             'Vapor-wet' : 10000,
#             'Dryout' : 10000,
#             'Transcritical' : 5000}

HX.set_htc(htc_type = 'Correlation_Disc', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 28
# HX.set_htc(htc_type = 'User-Defined', UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC) # 'User-Defined' or 'Correlation'

# HX.set_DP()                       # no pressure drop (SPHE pressure drop correlations: to come)
HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

HX.solve()
HX.plot_cells()

# # %%
# plot_spiral(HX)
# # W, B = plot_incidence(HX)
# # ax, W2, segs = plot_contacts(HX)

# print(f"Q = {HX.Q_dot/1e3:.2f} kW | N = {HX.params['n_turns']:.2f} turns | A = {HX.A_h:.3f} m2 | "
#       f"L_H = {HX.params['L_H']:.2f} m, L_C = {HX.params['L_C']:.2f} m")

