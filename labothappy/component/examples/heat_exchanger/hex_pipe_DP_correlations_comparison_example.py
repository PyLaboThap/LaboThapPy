"""
Comparison of the pipe pressure drop correlations of labothappy.correlations.pressure_drop,
used as a plausibility / regression check.

Interface (pipe_DP.py):
    pipe_DP(AS, geom, cond, flow_state, choices)
        flow_state : 'SP' (single phase), 'EVAP' or 'COND'
        choices    : {'SP': name} / {'EVAP': name} / {'COND': name}, or None for the automatic choice
    geom : D_t_i [m] (tube inner diameter), L [m], K [m] (roughness), n_t [-] (tubes in parallel)
    cond : m_dot [kg/s] (or G [kg/m2/s] of ONE tube) and the fluid state:
           single phase rho [kg/m3], mu [Pa s] (+ p [Pa], h [J/kg] for Cheng_CO2); two-phase p [Pa], x [-]
    Returns the pressure drop over the length L [Pa] (friction + gravity, no acceleration).

Four checks:
    A) Single-phase Darcy friction factor vs Reynolds number (smooth pipe),
       for all single-phase correlations except 'Cheng_CO2'.
    B) Single-phase pressure drop vs mass flow rate for water in a pipe,
       using the pipe_DP dispatcher (flow_state 'SP').
    C) Two-phase frictional pressure drop vs vapor quality for R134a
       ('Friedel' vs 'MSH', flow_state 'EVAP').
    D) Supercritical CO2 pressure drop vs temperature across the
       pseudo-boiling region ('Cheng_CO2').
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from labothappy.correlations.pressure_drop.pipe_DP_correlation_bank import (
    friction_factor_churchill,
    friction_factor_swamee_jain,
    friction_factor_haaland,
    friction_factor_konakov,
    friction_factor_petukhov)

from labothappy.correlations.pressure_drop.pipe_DP import (
    pipe_DP
)

fig, axs = plt.subplots(2, 2, figsize=(12, 9))

"--------- A) Single-phase friction factor vs Reynolds number (smooth pipe) ---------------------------"

Re_vec = np.logspace(2.5, 7, 200)
K_smooth = 0.0
d_hyd_ref = 0.02  # m (only sets K/d_hyd, irrelevant here since K=0)

f_correlations = {
    "Churchill":     [friction_factor_churchill(K_smooth, d_hyd_ref, Re) for Re in Re_vec],
    "Swamee-Jain":   [friction_factor_swamee_jain(K_smooth, d_hyd_ref, Re) for Re in Re_vec],
    "Haaland":       [friction_factor_haaland(K_smooth, d_hyd_ref, Re) for Re in Re_vec],
    "Konakov":       [friction_factor_konakov(Re) for Re in Re_vec],
    "Petukhov":      [friction_factor_petukhov(Re) for Re in Re_vec],
}

ax = axs[0, 0]
for name, f_vec in f_correlations.items():
    ax.plot(Re_vec, f_vec, label=name)
ax.axvline(4000, color='k', linestyle=':', linewidth=1, label='Re = 4000 (turbulent limit)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Reynolds number [-]')
ax.set_ylabel('Darcy friction factor f [-]')
ax.set_title('A) Friction factor vs Re (smooth pipe)\nSwamee-Jain/Konakov/Petukhov are turbulent-only')
ax.legend(fontsize=8)
ax.grid(True, which='both', alpha=0.3)

print("A) Friction factor at Re=1e5 (smooth pipe):")
idx_1e5 = np.argmin(np.abs(Re_vec - 1e5))
for name, f_vec in f_correlations.items():
    print(f"   {name:<14s}: f = {f_vec[idx_1e5]:.5f}")

"--------- B) Single-phase pressure drop vs mass flow rate (water) --------------------------------------"

fluid_1p = 'Water'
AS_1p = CP.AbstractState("HEOS", fluid_1p)
AS_1p.update(CP.PT_INPUTS, 2e5, 20 + 273.15)  # 2 bar, 20 degC (liquid)

pipe_geom_1p = {'D_t_i': 0.02, 'L': 10.0, 'K': 0.0, 'n_t': 1}  # smooth pipe, 20 mm ID, 10 m long
m_dot_vec = np.linspace(0.02, 2.0, 80)  # kg/s

correlations_1P = ['Churchill', 'Swamee_Jain', 'Haaland', 'Konakov', 'Petukhov']
dP_1P = {corr: [] for corr in correlations_1P}

for m_dot in m_dot_vec:
    # fluid state of the single-phase flow (constant here: 2 bar, 20 degC)
    cond_1p = {'m_dot': m_dot, 'rho': AS_1p.rhomass(), 'mu': AS_1p.viscosity()}
    for corr in correlations_1P:
        dP_1P[corr].append(
            pipe_DP(AS_1p, pipe_geom_1p, cond_1p, flow_state="SP", choices={"SP": corr})
        )

ax = axs[0, 1]
for corr in correlations_1P:
    ax.plot(m_dot_vec, np.array(dP_1P[corr]) * 1e-3, label=corr)
ax.set_xlabel('Mass flow rate [kg/s]')
ax.set_ylabel('Pressure drop [kPa]')
ax.set_title('B) Water pipe pressure drop vs mass flow rate\n(D=20 mm, L=10 m, smooth)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

idx_mid = len(m_dot_vec) // 2
print(f"\nB) Pressure drop at m_dot={m_dot_vec[idx_mid]:.3f} kg/s:")
for corr in correlations_1P:
    print(f"   {corr:<14s}: dP = {dP_1P[corr][idx_mid]:9.1f} Pa")

"--------- C) Two-phase frictional pressure drop vs vapor quality (R134a) -------------------------------"

fluid_2p = 'R134a'
T_sat = 10 + 273.15  # K
P_sat = CP.PropsSI('P', 'T', T_sat, 'Q', 0, fluid_2p)
AS_2p = CP.AbstractState("HEOS", fluid_2p)

pipe_geom_2p = {'D_t_i': 0.01, 'L': 3.0, 'K': 0.0, 'n_t': 1}  # 10 mm ID tube, 3 m long
m_dot_2p = 0.05  # kg/s

x_vec = np.linspace(0.02, 0.98, 60)
dP_friedel, dP_msh = [], []

for x in x_vec:
    # two-phase state: saturation pressure and vapor quality (properties computed from AS)
    cond_2p = {'m_dot': m_dot_2p, 'p': P_sat, 'x': x}
    dP_friedel.append(
        pipe_DP(AS_2p, pipe_geom_2p, cond_2p, flow_state="EVAP", choices={"EVAP": 'Friedel'})
    )
    dP_msh.append(
        pipe_DP(AS_2p, pipe_geom_2p, cond_2p, flow_state="EVAP", choices={"EVAP": 'MSH'})
    )

ax = axs[1, 0]
ax.plot(x_vec, np.array(dP_friedel) * 1e-3, label='Friedel')
ax.plot(x_vec, np.array(dP_msh) * 1e-3, label='MSH (Muller-Steinhagen & Heck)')
ax.set_xlabel('Vapor quality x [-]')
ax.set_ylabel('Frictional pressure drop [kPa]')
ax.set_title(f'C) R134a two-phase dP vs quality\n(T_sat={T_sat - 273.15:.0f} degC, D=10 mm, L=3 m)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

idx_mid_x = len(x_vec) // 2
print(f"\nC) Two-phase dP at x={x_vec[idx_mid_x]:.2f} (R134a, T_sat={T_sat - 273.15:.0f} degC):")
print(f"   Friedel        : dP = {dP_friedel[idx_mid_x]:9.1f} Pa")
print(f"   MSH            : dP = {dP_msh[idx_mid_x]:9.1f} Pa")

"--------- D) Supercritical CO2 pressure drop vs temperature (pseudo-boiling region) ---------------------"

fluid_co2 = 'CO2'
AS_co2 = CP.AbstractState("HEOS", fluid_co2)
P_co2 = 90e5  # Pa (supercritical, P_crit(CO2) = 73.8 bar)
T_co2_vec = np.linspace(0 + 273.15, 60 + 273.15, 100)  # K

pipe_geom_co2 = {'D_t_i': 0.008, 'L': 2.0, 'K': 0.0, 'n_t': 1}
m_dot_co2 = 0.05  # kg/s

dP_co2 = []
for T in T_co2_vec:
    AS_co2.update(CP.PT_INPUTS, P_co2, T)
    # supercritical state: Cheng_CO2 also needs the pressure and the enthalpy
    cond_co2 = {'m_dot': m_dot_co2, 'p': P_co2, 'h': AS_co2.hmass(),
                'rho': AS_co2.rhomass(), 'mu': AS_co2.viscosity()}
    dP_co2.append(
        pipe_DP(AS_co2, pipe_geom_co2, cond_co2, flow_state="SP", choices={"SP": 'Cheng_CO2'})
    )

ax = axs[1, 1]
ax.plot(T_co2_vec - 273.15, np.array(dP_co2) * 1e-3)
ax.set_xlabel('Temperature [degC]')
ax.set_ylabel('Pressure drop [kPa]')
ax.set_title(f'D) Supercritical CO2 dP vs T\n(P={P_co2 * 1e-5:.0f} bar, D=8 mm, L=2 m)\nCheng_CO2 correlation, crosses pseudo-boiling region')
ax.grid(True, alpha=0.3)

print(f"\nD) Supercritical CO2 dP range over T in [{T_co2_vec[0] - 273.15:.0f}, {T_co2_vec[-1] - 273.15:.0f}] degC:")
print(f"   min = {min(dP_co2):9.1f} Pa, max = {max(dP_co2):9.1f} Pa")

plt.tight_layout()
plt.show()
