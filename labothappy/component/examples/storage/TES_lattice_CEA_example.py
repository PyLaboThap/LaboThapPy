# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2026

@author: titouanjanod

Example of the PCM-TES lattice model (TES_lattice_CEA.py): standalone simulation (or sizing) of the TES alone,
then all the figures and the saving of the results.
Same cases as the __main__ of titouanjanod_perso/WP3/dynamic_modelling/lattice_TES_duCEA.py:
uncomment the case to simulate.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Absolute path to the LaboThapPy folder (4 levels above this file), if labothappy is not installed
LABOTHAPPY_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
if LABOTHAPPY_PATH not in sys.path:
    sys.path.append(LABOTHAPPY_PATH)

from labothappy.component.storage.TES_lattice_CEA import TES_simulation, create_result_folder, use_svg_inline


# Initilising variables if not needed because simulation and not sizing
Q_TES_predesign = None
T_target = None
T_HTF_min = None
T_HTF_max = None


SIM = TES_simulation()


""" TES Yvernault """
# fluid = "water"
# SIM.set_parameters(
#         T_PCM = 74 + 273.15,
#         cp_liq = 2.2e3,
#         Dh_latent = 238e3,
#         cp_sol = 2.4e3,
#         rho_PCM_sol = 939,
#         rho_PCM_liq = 858,
#         lambda_eq = 8
#         )
# Q_TES_predesign = 120*3600*1000 # J
# T_HTF_min = 20 + 273.15 # °C
# T_HTF_max = 85 + 273.15 # °C
# """ charging"""
# T_su_HTF = 85 +273.15# °C
# m_dot_tot = 525 / 3600 # kg/s
# p_HTF_init = 3e5 # Pa
# T_HTF_init = 85 +273.15 # °C
# T_PCM_eq_init = 40 +273.15 # °C
# T_target = 80 + 273.15 # °C

# """ discharging """
# T_su_HTF = 65 + 273.15# °C
# m_dot_tot = 266 / 3600 # kg/s
# p_HTF_init = 3e5 # Pa
# T_HTF_init = 65 +273.15 # °C
# T_PCM_eq_init = 84 +273.15 # °C
# T_target = 65 + 273.15 # °C


# n_box_parallel = 1
# n_box_serie = 1
# n_pass_box = 16
# n_segment_box = 32 # number of segment per length of the box
# Dt_ch_dis = 8*3600 # Time of charge or discharge, in seconds



""" TES Torrecid """
# fluid = "cyclopentane"
# SIM.set_parameters(
#         T_PCM = 168 + 273.15,
#         Dh_latent = 270e3,
#         cp_sol = 1.4e3,
#         cp_liq = 2.7e3,
#         rho_PCM_sol = 1520,
#         rho_PCM_liq = 1520,
#         lambda_eq = 16.5
#         )

# Q_TES_predesign = 10919.2*3600*1000 # J
# T_HTF_min = 31.8 + 273.15 # °C
# T_HTF_max = 190 + 273.15 # °C

# T_su_HTF = 190 +273.15# °C
# m_dot_tot = 14.34 # kg/s

# p_HTF_init = 35.81e5 # Pa
# T_HTF_init = 190 +273.15 # °C
# T_PCM_eq_init = 31.8 +273.15 # °C

# T_target = 172 + 273.15 # °C

# # n_box_parallel = 8
# # n_box_serie = 37
# n_pass_box = 1
# n_segment_box = 2



""" T_PCM = 92°C """
# fluid = "cyclopentane"
# SIM.set_parameters(
#         T_PCM = 97.2 + 273.15,
#         cp_sol = 1.5e3,
#         cp_liq = 2.8e3,
#         Dh_latent = 185e3,
#         rho_PCM_sol = 1414,
#         rho_PCM_liq = 1414,
#         lambda_eq = 8
#         )

# Q_TES_predesign = 10919.2*3600*1000 # J
# T_HTF_min = 140 + 273.15 # °C
# T_HTF_max = 190 + 273.15 # °C

# T_su_HTF = 190 +273.15# °C
# m_dot_tot = 14.34 # kg/s

# p_HTF_init = 35.81e5 # Pa
# T_HTF_init = 190 +273.15 # °C
# T_PCM_eq_init = 160 +273.15 # °C

# T_target = 172 + 273.15 # °C

# n_box_parallel = 8
# n_box_serie = 37
# n_pass_box = 1
# n_segment_box = 16



""" T_PCM = 150°C """
# fluid = "cyclopentane"
# SIM.set_parameters(
#         T_PCM = 150 + 273.15,
#         cp_sol = 1.7e3,
#         cp_liq = 2.5e3,
#         Dh_latent = 173e3,
#         rho_PCM_sol = 1620,
#         rho_PCM_liq = 1620,
#         lambda_eq = 8
#         )

# Q_TES_predesign = 10919.2*3600*1000 # J
# T_HTF_min = 140 + 273.15 # °C
# T_HTF_max = 190 + 273.15 # °C

# T_su_HTF = 190 +273.15# °C
# m_dot_tot = 14.34 # kg/s

# p_HTF_init = 35.81e5 # Pa
# T_HTF_init = 190 +273.15 # °C
# T_PCM_eq_init = 160 +273.15 # °C

# T_target = 172 + 273.15 # °C



""" T_PCM = 186°C """
# fluid = "cyclopentane"
# SIM.set_parameters(
#         T_PCM = 186 + 273.15,
#         cp_sol = 1.7e3,
#         cp_liq = 2.7e3,
#         Dh_latent = 358e3,
#         rho_PCM_sol = 1470,
#         rho_PCM_liq = 1470,
#         lambda_eq = 8
#         )
# """ charging - REC """
# Q_TES_predesign = 8097*3600*1000 # J
# T_HTF_min = 85.78 + 273.15 # °C
# T_HTF_max = 226 + 273.15 # °C
# T_su_HTF = 226 +273.15# °C
# m_dot_tot = 5.112 # kg/s
# p_HTF_init = 40.57e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 172.5 +273.15 # °C
# T_target = 196 + 273.15 # °C
# """ charging - no REC """
# Q_TES_predesign = 10069*3600*1000 # J
# T_HTF_min = 22.12 + 273.15 # °C
# T_HTF_max = 223 + 273.15 # °C
# T_su_HTF = 223 +273.15# °C
# m_dot_tot = 7.307 # kg/s
# p_HTF_init = 40e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 22.12 +273.15 # °C
# T_target = 196 + 273.15 # °C
# """ discharging - REC """
# Q_TES_predesign = 8097*3600*1000 # J
# T_su_HTF = 85.78 +273.15# °C
# m_dot_tot = 4.11 # kg/s
# p_HTF_init = 35.57e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 196 +273.15 # °C
# T_target = 172.5 + 273.15 # °C
# T_HTF_min = T_su_HTF
# T_HTF_max = T_PCM_eq_init
# """ discharging - no REC """
# Q_TES_predesign = 10069*3600*1000 # J
# T_su_HTF = 22.12 +273.15# °C
# m_dot_tot = 3.776 # kg/s
# p_HTF_init = 36.57e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 196 +273.15 # °C
# T_target = 153.3 + 273.15 # °C
# T_HTF_min = T_su_HTF
# T_HTF_max = T_PCM_eq_init

# n_pass_box = 1
# n_segment_box = 2
# n_box_parallel = 9
# n_box_serie = 23



""" T_PCM = 224°C """
# fluid = "cyclopentane"
# SIM.set_parameters(
#         T_PCM = 135 + 273.15,
#         cp_sol = 1.8e3,
#         cp_liq = 2.8e3,
#         Dh_latent = 260e3,
#         rho_PCM_sol = 1750,
#         rho_PCM_liq = 1750,
#         lambda_eq = 8
#         )
# """Charging with rec """
# Q_TES_predesign = 7971*3600*1000 # J
# T_HTF_min = 85.78 + 273.15 # °C
# T_HTF_max = 236 + 273.15 # °C
# T_su_HTF = 236 +273.15# °C
# m_dot_tot = 18.17 # kg/s
# p_HTF_init = 50e5 # Pa
# T_HTF_init = 236 +273.15 # °C
# T_PCM_eq_init = 85.78 +273.15 # °C
# T_target = 229 + 273.15 # °C
# """Charging no rec """
# Q_TES_predesign = 8801*3600*1000 # J
# T_HTF_min = 85.78 + 273.15 # °C
# T_HTF_max = 236 + 273.15 # °C
# T_su_HTF = 236 +273.15# °C
# m_dot_tot = 20.06 # kg/s
# p_HTF_init = 50e5 # Pa
# T_HTF_init = 236 +273.15 # °C
# T_PCM_eq_init = 85.78 +273.15 # °C
# T_target = 229 + 273.15 # °C
# """ discharging - REC """
# Q_TES_predesign = 7971*3600*1000 # J
# T_su_HTF = 85.78 +273.15 # °C
# m_dot_tot = 4.093 # kg/s
# p_HTF_init = 42e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 229 +273.15 # °C
# T_target = 171.6 + 273.15 # °C
# T_HTF_min = T_su_HTF
# T_HTF_max = T_PCM_eq_init
# """ discharging - no REC """
# Q_TES_predesign = 8801*3600*1000 # J
# T_su_HTF = 22.12 +273.15# °C
# m_dot_tot = 3.644 # kg/s
# p_HTF_init = 42e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 229 +273.15 # °C
# T_target = 142.7 + 273.15 # °C
# T_HTF_min = T_su_HTF
# T_HTF_max = T_PCM_eq_init

# n_pass_box = 1
# n_segment_box = 4
# n_box_parallel = 10
# n_box_serie = 19


""" SEHRENE experimental campain """
fluid = "water"
SIM.set_parameters(
        T_PCM = 118.9 + 273.15,
        cp_sol = 1.637e3,
        cp_liq = 2.846e3,
        Dh_latent = 335.5e3,
        rho_PCM_sol = 1448.4,
        rho_PCM_liq = 1281.5,
        lambda_eq = 8
        )
T_su_HTF = 135 +273.15 # °C
m_dot_tot = 700 / 3600 # kg/s
p_HTF_init = 6e5 # Pa
T_HTF_init = T_su_HTF
T_PCM_eq_init = 80 +273.15 # °C


n_subsegment_cell = 5 # number of HTF sub-segments per cell
n_pass_box = 4
n_segment_box = 16
n_box_parallel = 1
n_box_serie = 1



# """ Zorlu use case """
# fluid = "cyclopentane"
# SIM.set_parameters(
#         T_PCM = 135 + 273.15,
#         cp_sol = 1043,
#         cp_liq = 763.3,
#         Dh_latent = 221e3,
#         rho_PCM_sol = 1200,
#         rho_PCM_liq = 1200,
#         lambda_eq = 12
#         )
# """ discharging - REC """
# # Q_TES_predesign = 7971*3600*1000 # J
# T_su_HTF = 96.1 +273.15 # °C
# m_dot_tot = 62.74 # kg/s
# p_HTF_init = 22.258e5 # Pa
# T_HTF_init = T_su_HTF
# T_PCM_eq_init = 170.4 +273.15 # °C
# T_target = 124.8 + 273.15 # °C
# T_HTF_min = T_su_HTF
# T_HTF_max = T_PCM_eq_init

# n_pass_box = 1
# n_segment_box = 2
# n_box_parallel = 31
# n_box_serie = 200



SIM.set_parameters(n_segment_box = n_segment_box, n_pass_box = n_pass_box,
                   n_box_parallel = n_box_parallel, n_box_serie = n_box_serie,
                   n_subsegment_cell = n_subsegment_cell)
# SIM.set_parameters(htc_correlation = 'Hausen') # default: 'Gnielinski'


name_simu = "Erythritol comparaison Bénédicte 4"
test_solo = True # i.e. simu or design?
bypass_on = False
Dt_ch_dis = 9*3600 # Time of charge or discharge, in seconds
tau = 20 # time constant / time step, [s]

save_results = False
result_root = r"C:\Users\titou\OneDrive - Universite de Liege\!SEHRENE\WP3\TES modelling\Results"

# Experimental data of the CEA (B. Champel) on the figures: loaded with plots_tests_CEA.py of titouanjanod_perso
plot_exp_CEA = False
name_CEA_test = "2026-02-10 CH (80) 135 700kgh"
CEA_SCRIPT_DIR = r"C:\Users\titou\Documents\GitHub\titouanjanod_perso\TES_related"


if Q_TES_predesign is None:
    Q_TES_predesign = 7971 * 3600 * 1000  # J
    T_target = 171.6 + 273.15             # °C
    T_HTF_min = 25 + 273.15 # °C
    T_HTF_max = 100 + 273.15 # °C


use_svg_inline()

if save_results:
    result_folder = create_result_folder(result_root, name_simu)
else:
    result_folder = None


if test_solo==False:
    best_config = SIM.TES_sizing(T_target, Q_TES_predesign, T_HTF_min, T_HTF_max, Dt_ch_dis, tau, m_dot_tot, T_HTF_init,
                                 p_HTF_init, T_PCM_eq_init, fluid, bypass_on)

    print("best_config= ", best_config)

    n_box_serie = best_config["N_series"]
    n_box_parallel = best_config["N_parallel"]
    SIM.set_parameters(n_box_serie = n_box_serie, n_box_parallel = n_box_parallel)

    if save_results:
        with open(result_folder / f"best_config-{name_simu}.txt", "w") as f:

            f.write(f"N_parallel = {best_config['N_parallel']}\n")
            f.write(f"N_series = {best_config['N_series']}\n")
            f.write(f"N_box = {best_config['N_box']}\n")
            f.write(f"margin = {best_config['margin']}\n")
            f.write(f"index = {best_config['index']}\n")
            f.write(f"target_reached = {best_config['target_reached']}\n")


history = SIM.simulate(Dt_ch_dis, tau, m_dot_tot, T_HTF_init, p_HTF_init, T_PCM_eq_init, fluid, T_target, bypass_on)


print("Model type: 1D (no radial consideration)")
print("-----------")
print(f"tau= {tau} sec")
print(f"Nt= {len(history.time)} time discretisations")
print(f"Nx= {SIM.TES.Nx} axial discretisations")
print("-----------")
print(f"Q_tot = {history.Q_tot/3.6e6:.2f} kWh")
print(f"Solving time= {history.solving_time:.2f}")


exp = None
if plot_exp_CEA:
    if CEA_SCRIPT_DIR not in sys.path:
        sys.path.append(CEA_SCRIPT_DIR)
    from plots_tests_CEA import load_CEA_test, compute_derived
    exp = compute_derived(load_CEA_test(name_CEA_test))


figs = history.plot_all(save_dir=result_folder, name_simu=name_simu, exp=exp)

if save_results:
    history.save_results(result_folder, name_simu)

plt.show()
