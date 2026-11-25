# -*- coding: utf-8 -*-
"""
Spiral channel heat exchanger (SCHE) with HexMBChargeSensitive, geometry resolved by hex_geometry.py.
Cells are numbered from the core on both streams: the cold stream enters at the core, the hot one at the periphery.
"""
from labothappy.component.heat_exchanger.hex_MB import HexMB
from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.plot_SPHE import plot_spiral, plot_incidence, plot_contacts

HX = HexMB('SCHE')

case_study = "Evaporator"

# Correlation names (curved pipe dispatcher, evaluated on the elliptical channel with D = D_h):
#   1P : 'Sigalotti' (was "curved_pipe_htc")         2P hot (condensation) : 'Zakeralhoseini'
#   2P cold (evaporation) : 'Zhao' (was "Zhao_SPHE") or 'Nariai'; None / 'auto' = Zhao (laminar) / Nariai (turbulent)
# Pressure drop: not available yet for SCHE -> set_DP().

if case_study == "Condenser":
    HX.set_inputs(
        fluid_H='Propane', T_su_H=80 + 273.15, P_su_H=22e5, m_dot_H=0.02,
        fluid_C='Water', T_su_C=50 + 273.15, P_su_C=5e5, m_dot_C=0.3)
    
    Corr_H = {"1P": "Sigalotti", "2P": "Zakeralhoseini"}
    Corr_C = {"1P": "Sigalotti", "2P": None}
    params = dict(w_ch_C=0.00463, h_ch_C=0.00463, w_ch_H=0.00178, h_ch_H=0.00178,
                  D_ext=0.13, H=0.13, r0_in=3e-2, P_H=0.00496, P_V=0.00496, inner_ch='hot',
                  charge_compute = False)
    n_disc = 50

    Corr_H_DP = {"1P" : "Churchill", "2P" : "Choi"}
    Corr_C_DP = {"1P" : "Churchill", "2P" : "Choi"}

elif case_study == "Evaporator":
    HX.set_inputs(
        fluid_H='Water', T_su_H=10 + 273.15, P_su_H=5e5, m_dot_H=0.3,
        fluid_C='Propane', h_su_C=353418.2920873372, P_su_C=5.4e5, m_dot_C=0.02)
    
    Corr_H = {"1P": "Sigalotti", "2P": None}
    Corr_C = {"1P": "Sigalotti", "2P": "Zhao"}
    params = dict(w_ch_C=0.00444, h_ch_C=0.00444, w_ch_H=0.00495, h_ch_H=0.00495,
                  D_ext=0.13, H=0.071, r0_in=3e-2, P_H=0.00706, P_V=0.00706, inner_ch='hot',
                  charge_compute = False)
    
    n_disc = 50

    Corr_H_DP = {"1P" : "Churchill", "2P" : "Choi"}
    Corr_C_DP = {"1P" : "Churchill", "2P" : "Choi"}

# Geometry: give ONE of D_ext or n_turns. Optional: exchange_frac, k_wall, t_wall, n_parallel (n_series > 1 not available).
HX.set_parameters(**params, Flow_Type="CounterFlow", n_disc=n_disc)
HX.set_htc(htc_type='Correlation', Corr_H=Corr_H, Corr_C=Corr_C)
# HX.set_DP()
HX.set_DP(DP_type="Correlation_Disc", Corr_C=Corr_C_DP, Corr_H=Corr_H_DP)

#%%

HX.solve()
# HX.plot_cells()

# plot_spiral(HX)

# # print(f"Q = {HX.Q_dot/1e3:.2f} kW | N = {HX.params['n_turns']:.2f} turns | A_h = {HX.A_h:.3f} m2, A_c = {HX.A_c:.3f} m2")

