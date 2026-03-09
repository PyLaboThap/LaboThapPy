# -*- coding: utf-8 -*-
"""
ORC test case — circuit_mix
============================
Identical to the original circuit_rec example, adapted to use MixedCircuit.
The solve() method is now pluggable:

    orc.solve(method='successive_substitution')  # baseline
    orc.solve(method='wegstein')                 # Wegstein acceleration
    orc.solve(method='fsolve')                   # Newton-Powell (scipy fsolve)
    orc.solve(method='hybr')                     # Newton-Powell (scipy root)
    orc.solve(method='lm')                       # Levenberg-Marquardt
    orc.solve(method='broyden1')                 # quasi-Newton
    orc.solve(method='anderson')                 # Anderson mixing
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit_mix import MixedCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.expander.expander_semi_empirical import ExpanderSE
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity
from labothappy.toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP

# -------- 1) Instantiate Circuit --------
fluid = 'R1233zd(E)'
orc = MixedCircuit(fluid)

orc.mute_print()

# -------- 2) Create components --------
Expander   = ExpanderSE()
Condenser  = HexMBChargeSensitive('Plate')
Pump       = PumpCurveSimilarity()
Evaporator = HexMBChargeSensitive('Plate')

# -------- 3) Set component parameters --------

# Expander
Expander.set_parameters(
    AU_amb=9.3, AU_su_n=4.75, AU_ex_n=17.7, d_su1=6.48e-3, m_dot_n=0.1,
    A_leak=9.99e-06, W_dot_loss_0=2.37e+1, alpha=1.16e-1, C_loss=1.13,
    rv_in=1.7, V_s=2*0.0000712, mode='P_M'
)

# Evaporator
evaporator_geom = PlateGeomSWEP()
evaporator_geom.set_parameters("P200THx140/1P_Evaporator")

Evaporator.set_parameters(
    A_c=evaporator_geom.A_c, A_h=evaporator_geom.A_h,
    h=evaporator_geom.h, l=evaporator_geom.l,
    l_v=evaporator_geom.l_v, w_v=evaporator_geom.w_v,
    C_CS=evaporator_geom.C_CS, C_Dh=evaporator_geom.C_Dh,
    C_V_tot=evaporator_geom.C_V_tot, C_canal_t=evaporator_geom.C_canal_t,
    C_n_canals=evaporator_geom.C_n_canals,
    H_CS=evaporator_geom.H_CS, H_Dh=evaporator_geom.H_Dh,
    H_V_tot=evaporator_geom.H_V_tot, H_canal_t=evaporator_geom.H_canal_t,
    H_n_canals=evaporator_geom.H_n_canals,
    casing_t=evaporator_geom.casing_t, chevron_angle=evaporator_geom.chevron_angle,
    fooling=evaporator_geom.fooling, n_plates=evaporator_geom.n_plates,
    plate_cond=evaporator_geom.plate_cond, plate_pitch_co=evaporator_geom.plate_pitch_co,
    t_plates=evaporator_geom.t_plates, w=evaporator_geom.w,
    amplitude=evaporator_geom.amplitude, phi=evaporator_geom.phi,
    Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0,
)
Evaporator.set_htc(
    htc_type="Correlation",
    Corr_H={"1P": "water_plate_HTC"},
    Corr_C={"1P": "martin_holger_plate_HTC", "2P": "amalfi_plate_HTC"}
)
Evaporator.set_DP()

# Condenser
condenser_geom = PlateGeomSWEP()
condenser_geom.set_parameters("P200THx140/1P_Condenser")

Condenser.set_parameters(
    A_c=condenser_geom.A_c, A_h=condenser_geom.A_h,
    h=condenser_geom.h, l=condenser_geom.l,
    l_v=condenser_geom.l_v, w_v=condenser_geom.w_v,
    C_CS=condenser_geom.C_CS, C_Dh=condenser_geom.C_Dh,
    C_V_tot=condenser_geom.C_V_tot, C_canal_t=condenser_geom.C_canal_t,
    C_n_canals=condenser_geom.C_n_canals,
    H_CS=condenser_geom.H_CS, H_Dh=condenser_geom.H_Dh,
    H_V_tot=condenser_geom.H_V_tot, H_canal_t=condenser_geom.H_canal_t,
    H_n_canals=condenser_geom.H_n_canals,
    casing_t=condenser_geom.casing_t, chevron_angle=condenser_geom.chevron_angle,
    fooling=condenser_geom.fooling, n_plates=condenser_geom.n_plates,
    plate_cond=condenser_geom.plate_cond, plate_pitch_co=condenser_geom.plate_pitch_co,
    t_plates=condenser_geom.t_plates, w=condenser_geom.w,
    amplitude=condenser_geom.amplitude, phi=condenser_geom.phi,
    Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=10,
)
Condenser.set_htc(
    htc_type="Correlation",
    Corr_H={"1P": "martin_holger_plate_HTC", "2P": "shah_condensation_plate_HTC"},
    Corr_C={"1P": "water_plate_HTC"}
)
Condenser.set_DP()

# Pump
V_dot_curve    = np.array([20, 30, 40, 50, 60, 70, 80])
Delta_H_curve  = np.array([57, 55, 52, 49, 45, 42, 36])
eta_is_curve   = np.array([0.45, 0.59, 0.69, 0.75, 0.79, 0.79, 0.75])
NPSH_r_curve   = np.array([1.1, 1.1, 1.4, 1.8, 2, 3, 4.7])
N_rated        = 2900  # RPM

Pump.set_parameters(
    V_dot_curve=V_dot_curve,
    Delta_H_curve=Delta_H_curve,
    eta_is_curve=eta_is_curve,
    NPSH_r_curve=NPSH_r_curve,
    N_rot_rated=N_rated,
    mode="P_M",
)

# -------- 4) Add components to circuit --------
orc.add_component(Expander,   "Expander")
orc.add_component(Condenser,  "Condenser")
orc.add_component(Pump,       "Pump")
orc.add_component(Evaporator, "Evaporator")

# -------- 5) Link components --------
orc.link_components("Expander",   "m-ex",   "Condenser",  "m-su_H")
orc.link_components("Condenser",  "m-ex_H", "Pump",       "m-su")
orc.link_components("Pump",       "m-ex",   "Evaporator", "m-su_C")
orc.link_components("Evaporator", "m-ex_C", "Expander",   "m-su")

# -------- 6) Add fluid sources --------
CD_source    = MassConnector('Water')
T_su_w_cd    = 5 + 273.15
P_su_w_cd    = 2e5
m_dot_w_cd   = 3.0   # kg/s

EV_source    = MassConnector('Water')
T_su_w_ev    = 70 + 273.15
P_su_w_ev    = 2e5
m_dot_w_ev   = 2.5   # kg/s

orc.add_source("CD_Water", CD_source, orc.components["Condenser"],  "m-su_C")
orc.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd,
                           m_dot=m_dot_w_cd, target="CD_Water")

orc.add_source("EV_Water", EV_source, orc.components["Evaporator"], "m-su_H")
orc.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev,
                           m_dot=m_dot_w_ev, target="EV_Water")

# -------- 7) Set cycle guesses --------
m_dot_ref  = 0.4   # kg/s
SC_cd      = 5     # K
N_exp      = 6000  # RPM
T_amb      = 293   # K

P_LP_guess = PropsSI("P", "T", T_su_w_cd + 10, "Q", 0, fluid)
P_HP_guess = PropsSI("P", "T", T_su_w_ev - 10, "Q", 1, fluid)

orc.set_cycle_guess(target="Pump:su",      m_dot=m_dot_ref, SC=SC_cd, p=P_LP_guess)
orc.set_cycle_guess(target="Pump:ex",      p=P_HP_guess)
orc.set_cycle_guess(target="Expander:W",   N_rot=N_exp)
orc.set_cycle_guess(target="Expander:Q_amb", T_amb=T_amb)
orc.set_cycle_guess(target="Expander:ex",  p=P_LP_guess)

# -------- 8) Solve — swap method here for comparison --------
METHOD = 'wegstein'   # <-- change to compare: 'successive_substitution',
                      #     'wegstein', 'fsolve', 'lm', 'broyden1', 'anderson'

orc.solve(method=METHOD)

print(f"\n[{METHOD}] Converged: {orc.converged}")
print(f"  P_HP = {Expander.su.p:.2f} Pa")
print(f"  P_LP = {Expander.ex.p:.2f} Pa")
print(f"  T_su_expander = {Expander.su.T - 273.15:.2f} °C")
print(f"  m_dot = {Pump.su.m_dot:.4f} kg/s")
print(f"  W_exp = {Expander.W.W_dot:.2f} W")
print(f"  W_pump = {Pump.W.W_dot:.2f} W")
print(f"  Q_ev = {Evaporator.Q.Q_dot:.2f} W")
print(f"  Q_cd = {Condenser.Q.Q_dot:.2f} W")

