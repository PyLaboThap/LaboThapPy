# -*- coding: utf-8 -*-
"""
Example of use of the TeeJunction component (diverging tee, EES / Crane TP-410 losses).

Four cases:
    1. Reproduction of the EES help example (water, v = 10 m/s, smooth tee)
    2. Subcooled R134a liquid line splitting to two circuits, with refrigerant charge
    3. Two-phase R134a after an expansion valve, reducing tee (12 mm run, 8 mm branch)
    4. Sweep of the flow split: effect of velocity_basis ('combined' vs 'own')
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.component.pipe.tee_junction import TeeJunction   # adapt to where tee_junction.py is placed


# ============================================================================
# CASE 1 - Reproduction of the EES example
# ============================================================================
# EES help pages (K_tee_straight / K_tee_branched): RelRough = 0, V = 10 m/s,
# rho = 1000 kg/m3  ->  DELTAP = 6028 Pa (run) and 18083 Pa (branch)

print("\n#### CASE 1 - EES example (water, 10 m/s) ####")

D = 20e-3                                      # [m]
P_su = 3e5                                     # [Pa]
T_su = 20 + 273.15                             # [K]
h_su = PropsSI('H', 'P', P_su, 'T', T_su, 'Water')
rho = PropsSI('D', 'P', P_su, 'T', T_su, 'Water')
m_dot = rho * 10 * np.pi * D**2 / 4            # mass flow giving v = 10 m/s

tee_ees = TeeJunction(outlet_repartition=[0.5, 0.5])   # [run, branch]
tee_ees.set_parameters(D=D)                    # defaults: roughness = 0, velocity_basis = 'combined'

tee_ees.su.set_fluid('Water')
tee_ees.su.set_p(P_su)
tee_ees.su.set_h(h_su)
tee_ees.su.set_m_dot(m_dot)

tee_ees.solve()
tee_ees.print_results()

# Scale the result to rho = 1000 kg/m3 to compare with the EES numbers
print(f"dP run    (rho = 1000) : {tee_ees.DP[0] * 1000 / rho:.0f} Pa   (EES: 6028 Pa)")
print(f"dP branch (rho = 1000) : {tee_ees.DP[1] * 1000 / rho:.0f} Pa   (EES: 18083 Pa)")


# ============================================================================
# CASE 2 - Subcooled liquid line (R134a), split to two circuits, with charge
# ============================================================================

print("\n#### CASE 2 - R134a liquid line, 60/40 split, charge ####")

fluid = 'R134a'
P_su = 10e5                                    # [Pa]
T_su = 30 + 273.15                             # [K]  (subcooled at 10 bar)
h_su = PropsSI('H', 'P', P_su, 'T', T_su, fluid)

tee_liq = TeeJunction(outlet_repartition=[0.6, 0.4])
tee_liq.set_parameters(
    D=10e-3,               # run inner diameter [m]
    roughness=1.5e-6,      # drawn copper tube [m]
    L_run=35e-3,           # internal length inlet -> run outlet [m]
    L_branch=20e-3,        # internal length of the branch, from the run wall [m]
)

tee_liq.su.set_fluid(fluid)
tee_liq.su.set_p(P_su)
tee_liq.su.set_h(h_su)
tee_liq.su.set_m_dot(0.05)                     # [kg/s]

tee_liq.solve()
tee_liq.print_results()

# Accessing the results
print(f"P at run outlet    : {tee_liq.ex_run.p / 1e5:.5f} bar")
print(f"P at branch outlet : {tee_liq.ex_branch.p / 1e5:.5f} bar")
print(f"Charge breakdown   : run {tee_liq.charge_detail['run'] * 1e3:.3f} g, "
      f"branch {tee_liq.charge_detail['branch'] * 1e3:.3f} g")


# ============================================================================
# CASE 3 - Two-phase flow after an expansion valve, reducing tee
# ============================================================================

print("\n#### CASE 3 - Two-phase R134a, reducing tee 12 mm / 8 mm ####")

P_su = 3e5                                     # [Pa]
x_su = 0.25                                    # inlet vapour quality [-]
h_su = PropsSI('H', 'P', P_su, 'Q', x_su, fluid)

tee_2ph = TeeJunction(outlet_repartition=[0.7, 0.3])
tee_2ph.set_parameters(
    D=12e-3,
    D_branch=8e-3,                 # reducing tee: sudden contraction added on the branch
    roughness=1.5e-6,
    velocity_basis='own',          # each K uses the velocity of its own outlet stream
    L_run=40e-3,
    L_branch=30e-3,
    void_fraction_model='Zivi',    # for the charge (the pressure drop uses the homogeneous density)
)

tee_2ph.su.set_fluid(fluid)
tee_2ph.su.set_p(P_su)
tee_2ph.su.set_h(h_su)
tee_2ph.su.set_m_dot(0.03)

tee_2ph.solve()
tee_2ph.print_results()

d = tee_2ph.DP_detail
print(f"Branch: tee {d['tee'][1]:.1f} Pa + reducer {d['reducer'][1]:.1f} Pa "
      f"(K_reducer = {d['K_reducer']:.3f})")


# ============================================================================
# CASE 4 - Sweep of the flow split: 'combined' vs 'own' velocity basis
# ============================================================================

print("\n#### CASE 4 - Flow-split sweep (R134a liquid, D = 10 mm, 0.05 kg/s) ####")

h_su = PropsSI('H', 'P', 10e5, 'T', 30 + 273.15, fluid)
fractions_branch = np.linspace(0.1, 0.9, 9)

print(f"{'x_branch':>9} | {'dP_run comb.':>12} {'dP_br comb.':>12} | {'dP_run own':>11} {'dP_br own':>10}  [Pa]")
for x_b in fractions_branch:
    results = {}
    for basis in ('combined', 'own'):
        tee = TeeJunction(outlet_repartition=[1 - x_b, x_b])
        tee.set_parameters(D=10e-3, roughness=1.5e-6, velocity_basis=basis)
        tee.print_flag = False
        tee.su.set_fluid(fluid)
        tee.su.set_p(10e5)
        tee.su.set_h(h_su)
        tee.su.set_m_dot(0.05)
        tee.solve()
        results[basis] = tee.DP
    print(f"{x_b:9.2f} | {results['combined'][0]:12.1f} {results['combined'][1]:12.1f} | "
          f"{results['own'][0]:11.1f} {results['own'][1]:10.1f}")

# With 'combined' the losses do not depend on the split (EES/Crane behaviour).
# With 'own' they scale with the square of each stream's mass flow.
