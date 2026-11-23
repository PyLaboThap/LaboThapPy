# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2026

@author: titouanjanod

Example of the PCM-TES lattice model (TES_lattice_CEA.py): simulation of the TES alone, with a constant HTF inlet,
or sizing of the number of modules. Then all the figures, with the measurements if available, and the saving of
the results.

Choose the campaign and the settings in "User choices", the data of each campaign are in "Campaigns".
"""

import os
import sys
import matplotlib.pyplot as plt

# Absolute path to the LaboThapPy folder (4 levels above this file), if labothappy is not installed
LABOTHAPPY_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
if LABOTHAPPY_PATH not in sys.path:
    sys.path.append(LABOTHAPPY_PATH)

from labothappy.component.storage.TES_lattice_CEA import TES_simulation, create_result_folder 

# =============================================================================
# User choices
# =============================================================================
# "Yvernault_charge", "Yvernault_discharge", "CEA_SoC_erythritol", "Zorlu_UC" or "Torrecid_charge_sizing"
campaign = "CEA_SoC_erythritol"

do_sizing = False # True: sizing of the number of modules (only for the campaigns with sizing data), then simulation
bypass_on = False # True: part of the HTF bypasses the TES, to keep the mixed outlet at T_target

duration_simu_forced = None # [s], if not None: replaces the duration of the campaign
dt = 20 # time step [s]
n_segment_module = 32 # number of cells along one module
n_subsegment_cell = 5 # number of HTF sub-segments per cell (the PCM stays at one temperature per cell)

save_results = False
run_name = "Erythritol comparaison Bénédicte 5" # name of the results folder, and suffix of the saved files
result_root = r"C:\Users\titou\OneDrive - Universite de Liege\!SEHRENE\WP3\TES modelling\Results"
show_plots = True

SIM = TES_simulation()
# optional, default: 'Gnielinski'
SIM.set_parameters(htc_correlation = 'Hausen')


# =============================================================================
# Experimental data on the figures (optional)
# =============================================================================
# plot_all(exp=exp_data) adds the measurements (crosses) on the figures of the model.
# exp_data is a pandas DataFrame or a dict of arrays, with these columns:
#   t_h          time from the start of the test [h]
#   Q_dot_kW     heat rate exchanged between the HTF and the TES, > 0 in charge and in discharge [kW]
#                (Q_dot_smooth_kW is used instead if present)
#   Q_TES_kWh    energy exchanged since the start, > 0 in charge and in discharge [kWh]
#   T_HTF_ex_C   HTF outlet temperature [°C]
#   T_PCM_min_C, T_PCM_avg_C, T_PCM_max_C   (optional) min / mean / max of the PCM thermocouples [°C]
# With your own data, e.g. a csv file whose columns are renamed as above:
#   import pandas as pd
#   exp_data = pd.read_csv("my_test.csv")
exp_data = None
plot_exp = True # plot the measurements of the campaign, if it has some


# =============================================================================
# Campaigns
# =============================================================================
if campaign in ("Yvernault_charge", "Yvernault_discharge"):
    # Tests of Yvernault, L., Vesin, S., & Pouvreau, J. (2026) https://doi.org/10.1016/j.est.2026.120714
    fluid = "water"
    SIM.set_parameters(
            T_PCM = 74 + 273.15,
            cp_liq = 2.2e3,
            Dh_latent = 238e3,
            cp_sol = 2.4e3,
            rho_PCM_sol = 939,
            rho_PCM_liq = 858,
            lambda_eq = 8
            )
    p_HTF = 3e5 # Pa
    duration_simu = 8*3600 # s

    if campaign == "Yvernault_charge":
        T_HTF_su = 85 + 273.15 # K
        m_dot_tot = 525 / 3600 # kg/s
        T_PCM_eq_init = 40 + 273.15 # K
        T_target = 80 + 273.15 # K
    else:
        T_HTF_su = 65 + 273.15 # K
        m_dot_tot = 266 / 3600 # kg/s
        T_PCM_eq_init = 84 + 273.15 # K
        T_target = 65 + 273.15 # K

    n_pass_module = 16
    n_module_parallel = 1
    n_module_serie = 1


elif campaign == "CEA_SoC_erythritol":
    # CEA Grenoble campaign for the State of Charge (SoC) measurement (B. Champel), SEHRENE deliverable D4.7
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
    T_HTF_su = 135 + 273.15 # K
    m_dot_tot = 700 / 3600 # kg/s
    p_HTF = 6e5 # Pa
    T_PCM_eq_init = 80 + 273.15 # K
    T_target = None # no bypass in these tests
    duration_simu = 9*3600 # s

    n_pass_module = 4
    n_module_parallel = 1
    n_module_serie = 1

    # measurements, loaded with titouanjanod_perso/TES_related/plots_tests_CEA.py
    if plot_exp:
        CEA_SCRIPT_DIR = r"C:\Users\titou\Documents\GitHub\titouanjanod_perso\TES_related"
        if CEA_SCRIPT_DIR not in sys.path:
            sys.path.append(CEA_SCRIPT_DIR)
        from plots_tests_CEA import load_CEA_test, compute_derived
        exp_data = compute_derived(load_CEA_test("2026-02-10 CH (80) 135 700kgh"))


elif campaign == "Zorlu_UC":
    # Discharge of the TES, initial data of SEHRENE deliverable D3.3
    fluid = "cyclopentane"
    SIM.set_parameters(
            T_PCM = 135 + 273.15,
            cp_sol = 1043,
            cp_liq = 763.3,
            Dh_latent = 221e3,
            rho_PCM_sol = 1200,
            rho_PCM_liq = 1200,
            lambda_eq = 12
            )
    T_HTF_su = 96.1 + 273.15 # K
    m_dot_tot = 62.74 # kg/s
    p_HTF = 22.258e5 # Pa
    T_PCM_eq_init = 170.4 + 273.15 # K
    T_target = 124.8 + 273.15 # K
    duration_simu = 9*3600 # s

    # sizing data
    Q_TES_predesign = 7971*3600*1000 # J
    T_HTF_min = T_HTF_su
    T_HTF_max = T_PCM_eq_init

    # replaced by the result of the sizing if do_sizing
    n_pass_module = 1
    n_module_parallel = 31
    n_module_serie = 200


elif campaign == "Torrecid_charge_sizing":
    # Charge, initial data of SEHRENE deliverable D3.3
    fluid = "cyclopentane"
    SIM.set_parameters(
            T_PCM = 168 + 273.15,
            Dh_latent = 270e3,
            cp_sol = 1.4e3,
            cp_liq = 2.7e3,
            rho_PCM_sol = 1520,
            rho_PCM_liq = 1520,
            lambda_eq = 16.5
            )
    T_HTF_su = 190 + 273.15 # K
    m_dot_tot = 14.34 # kg/s
    p_HTF = 35.81e5 # Pa
    T_PCM_eq_init = 31.8 + 273.15 # K
    T_target = 172 + 273.15 # K
    duration_simu = 9*3600 # s

    # sizing data
    Q_TES_predesign = 10919.2*3600*1000 # J
    T_HTF_min = 31.8 + 273.15 # K
    T_HTF_max = 190 + 273.15 # K

    # replaced by the result of the sizing if do_sizing
    n_pass_module = 1
    n_module_parallel = 8
    n_module_serie = 37

else:
    raise ValueError(f"Unknown campaign: '{campaign}'")


if duration_simu_forced is not None:
    duration_simu = duration_simu_forced

# "module" in this example = "box" in the model
SIM.set_parameters(n_segment_box = n_segment_module, n_pass_box = n_pass_module,
                   n_box_parallel = n_module_parallel, n_box_serie = n_module_serie,
                   n_subsegment_cell = n_subsegment_cell)


# =============================================================================
# Sizing (optional) and simulation
# =============================================================================

if save_results:
    result_folder = create_result_folder(result_root, run_name)
else:
    result_folder = None

if do_sizing:
    if campaign not in ("Zorlu_UC", "Torrecid_charge_sizing"):
        raise ValueError(f"No sizing data for the campaign '{campaign}'")

    best_config = SIM.TES_sizing(T_target, Q_TES_predesign, T_HTF_min, T_HTF_max, duration_simu, dt, m_dot_tot,
                                 T_HTF_su, p_HTF, T_PCM_eq_init, fluid, bypass_on)
    print("best_config = ", best_config)

    SIM.set_parameters(n_box_serie = best_config["N_series"], n_box_parallel = best_config["N_parallel"])

    if save_results:
        with open(result_folder / f"best_config-{run_name}.txt", "w") as f:
            for key in ("N_parallel", "N_series", "N_box", "margin", "index", "target_reached"):
                f.write(f"{key} = {best_config[key]}\n")

history = SIM.simulate(duration_simu, dt, m_dot_tot, T_HTF_su, p_HTF, T_PCM_eq_init, fluid, T_target, bypass_on)


# =============================================================================
# Results
# =============================================================================
print("Model type: 1D (no radial consideration)")
print("-----------")
print(f"dt = {dt} s")
print(f"Nt = {len(history.time)} time steps")
print(f"N_cells = {SIM.TES.Nx} cells along the HTF flow")
print(f"n_subsegment_cell = {n_subsegment_cell} HTF sub-segments per cell")
print("-----------")
print(f"Q_tot = {history.Q_tot/3.6e6:.2f} kWh")
print(f"Solving time = {history.solving_time:.2f} s")

if not plot_exp:
    exp_data = None

figs = history.plot_all(save_dir=result_folder, name_simu=run_name, exp=exp_data)

if save_results:
    history.save_results(result_folder, run_name)

if show_plots:
    plt.show()
