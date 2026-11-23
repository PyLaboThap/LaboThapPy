# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2026

@author: titouanjanod

Example of the PCM-TES lattice model (TES_lattice_CEA.py): simulation of the TES alone, with a constant HTF supply
or with the measured HTF supply of a CEA test (m_dot(t), T_HTF,su(t)), or sizing of the number of modules.
Then all the figures, with the measurements if available, and the saving of the results.

Choose the campaign and the settings in "User choices", the data of each campaign are in "Campaigns".
"""

import os
import re
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Absolute path to the LaboThapPy folder (4 levels above this file), if labothappy is not installed
LABOTHAPPY_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
if LABOTHAPPY_PATH not in sys.path:
    sys.path.append(LABOTHAPPY_PATH)

from labothappy.component.storage.TES_lattice_CEA import TES_simulation, create_result_folder
from labothappy.component.storage.TES_lattice_CEA import skip_exp as skip_exp_default

# =============================================================================
# User choices
# =============================================================================
# "Yvernault_charge", "Yvernault_discharge", "CEA_SoC_erythritol", "Zorlu_UC" or "Torrecid_charge_sizing"
campaign = "CEA_SoC_erythritol"

do_sizing = False # True: sizing of the number of modules (only for the campaigns with sizing data), then simulation
bypass_on = False # True: part of the HTF bypasses the TES, to keep the mixed outlet at T_target

# -----------------------------------------------------------------------------
# HTF supply (= what enters the TES): constant or measured
# -----------------------------------------------------------------------------
# "constant": T_HTF_su and m_dot_tot stay the same during the whole simulation. Their values are in "Campaigns"
#             below. This is the classical run, it works for every campaign.
#             For "CEA_SoC_erythritol", the crosses on the figures are the test "2026-02-10 CH (80) 135 700kgh"
#             (same 135 °C and 700 kg/h as the constant values).
# "measured": the model is fed, at each time step, with the m_dot(t) and T_HTF,su(t) MEASURED during the CEA test
#             CEA_test (data every 5 s, linearly interpolated at each dt). The PCM starts uniform at the measured
#             T_PCM,avg. Only for campaign = "CEA_SoC_erythritol", without sizing and without bypass.
HTF_supply = "measured"

# ---- the settings below are only used if HTF_supply = "measured" ----

# name of the xlsx file of the test (without .xlsx), in the folder DATA_DIR of plots_tests_CEA.py
CEA_test = "2026-02-11 CH (80) 135 1500kgh"

# Moment where the simulation starts:
#   "raw"  : at the FIRST ROW of the file. The CEA first brought the PCM to 80 °C (charges) before the real test:
#            this pre-heating is simulated too, from a PCM that is ~uniform (e.g. ~52 °C for the 700 kg/h test).
#            The start of the real test is drawn as a dotted vertical line on the figures vs time.
#   "test" : at the START OF THE REAL TEST (pre-heating removed, t = 0 as in plots_tests_CEA.py). The PCM is then
#            NOT uniform anymore in reality (warmer near the tubes), but it is taken uniform at T_PCM,avg in the model.
measured_start = "raw"

# Zoom on what happens BEFORE the start of the real test (only with measured_start = "raw"):
# two extra figures, Q_dot(t) and T(t), from t = 0 to (start of the real test + zoom_after_start_h).
zoom_preheating = False
zoom_after_start_h = 0 # [h] how long after the start of the real test the zoom goes on
skip_exp_zoom = 20 # on the zoomed figures, one measured point out of skip_exp_zoom (1 = every point, i.e. every 5 s)

# Warning if the PCM is not uniform at the start of the simulation (the uniform initial T_PCM is then rough):
dT_PCM_uniform_max = 5 # [K] max allowed T_PCM,max - T_PCM,min at the first row

# -----------------------------------------------------------------------------
# Time and space discretisation
# -----------------------------------------------------------------------------
duration_simu_forced = 8 *3600 # [s] (!) in SECONDS, e.g. 8*3600 for 8 h. If None: duration of the campaign / of the test
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
# Checks
# =============================================================================
if HTF_supply == "measured":
    if campaign != "CEA_SoC_erythritol":
        raise ValueError(f"HTF_supply = 'measured' only with the campaign 'CEA_SoC_erythritol', not '{campaign}'")
    if do_sizing or bypass_on:
        raise ValueError("HTF_supply = 'measured': no sizing and no bypass (do_sizing = bypass_on = False)")
    if measured_start not in ("raw", "test"):
        raise ValueError(f"Unknown measured_start: '{measured_start}' ('raw' or 'test')")
elif HTF_supply != "constant":
    raise ValueError(f"Unknown HTF_supply: '{HTF_supply}' ('constant' or 'measured')")


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
    # constant supply: nominal values of the test "2026-02-10 CH (80) 135 700kgh" (not used if HTF_supply = "measured")
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
    #   constant supply: test with the same nominal values as above, for the crosses on the figures
    #   measured supply: test CEA_test, it gives the inputs of the model AND the crosses on the figures
    if plot_exp or HTF_supply == "measured":
        CEA_SCRIPT_DIR = r"C:\Users\titou\Documents\GitHub\titouanjanod_perso\TES_related"
        if CEA_SCRIPT_DIR not in sys.path:
            sys.path.append(CEA_SCRIPT_DIR)
        from plots_tests_CEA import load_CEA_test, compute_derived

        if HTF_supply == "measured":
            test_name = CEA_test
            # the rows before the start of the real test (pre-heating) are kept only with measured_start = "raw"
            cut_preheating = (measured_start == "test")
        else:
            test_name = "2026-02-10 CH (80) 135 700kgh"
            cut_preheating = True # the constant run starts at the start of the real test
        # exp_data: DataFrame with t_s, t_h, m_dot_kg_s, T_HTF_su_C, T_HTF_ex_C, T_PCM_min/avg/max_C, Q_dot_kW,
        # Q_dot_smooth_kW (moving average over 60 s), Q_TES_kWh (integral of the raw Q_dot from t = 0)...
        exp_data = compute_derived(load_CEA_test(test_name, cut_preheating=cut_preheating))

    if HTF_supply == "measured":
        # initial state of the model: PCM uniform at the mean of the 48 thermocouples at t = 0
        # (below the melting zone, h is linear in T -> same stored energy as the real non-uniform PCM)
        T_PCM_eq_init = exp_data["T_PCM_avg_C"].iloc[0] + 273.15 # K
        dT_PCM_init = exp_data["T_PCM_max_C"].iloc[0] - exp_data["T_PCM_min_C"].iloc[0] # K
        if dT_PCM_init > dT_PCM_uniform_max:
            warnings.warn(f"{CEA_test}: T_PCM,max - T_PCM,min = {dT_PCM_init:.1f} K at the start of the simulation "
                          f"(> {dT_PCM_uniform_max} K), the uniform initial PCM temperature is a rough assumption")
        duration_simu = exp_data["t_s"].iloc[-1] # s

        # title of the figures: nominal values of the test (file name) and date (one test)
        flow_nom = re.search(r"(\d+)\s*kgh", CEA_test)
        title_parts = [exp_data.attrs["mode"].capitalize(),
                       rf"$T_{{HTF,in}}$ = {exp_data.attrs['T_HTF_in_nom_C']:g} °C"]
        if flow_nom:
            title_parts.append(f"{flow_nom.group(1)} kg/h")
        title_parts.append(exp_data.attrs["date"])
        title_measured = " – ".join(title_parts)
        # start of the real charge / discharge in the time of the simulation [h]
        mode_test = exp_data.attrs["mode"]
        t_start_test_h = exp_data.attrs["t_start_raw_s"]/3600 if measured_start == "raw" else 0.0


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

if duration_simu < dt:
    raise ValueError(f"duration_simu = {duration_simu} s < dt = {dt} s: no time step to compute. "
                     f"(duration_simu_forced is in seconds, e.g. 8*3600 for 8 h)")

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

if HTF_supply == "measured":
    # m_dot(t) and T_HTF,su(t) of the test, interpolated at each time step
    history = SIM.simulate_measured(exp_data["t_s"], exp_data["m_dot_kg_s"], exp_data["T_HTF_su_C"] + 273.15,
                                    dt, p_HTF, T_PCM_eq_init, fluid, duration_simu)
    history.title = title_measured
else:
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

if HTF_supply == "measured":
    # energy exchanged with the HTF since t = 0, model vs measurements [kWh] (> 0 in charge and in discharge)
    sign = 1 if history.get_mode() == "charge" else -1
    # Q_dot_tot[n] is the step n -> n+1: energy at time[n] = sum of the steps before
    Q_TES_model_kWh = sign*np.concatenate(([0], history.Q_cum[:-1]))/3.6e6
    t_end_h = history.time[-1]/3600

    def Q_TES_at(t_h):
        return (np.interp(t_h, history.time/3600, Q_TES_model_kWh),
                np.interp(t_h, exp_data["t_h"], exp_data["Q_TES_kWh"]))

    print("-----------")
    print(f"Measured supply: {CEA_test}, start = '{measured_start}'")
    print(f"T_PCM,init = {T_PCM_eq_init-273.15:.2f} °C (uniform, measured T_PCM,min/max = "
          f"{exp_data['T_PCM_min_C'].iloc[0]:.2f} / {exp_data['T_PCM_max_C'].iloc[0]:.2f} °C)")
    Q_start = Q_TES_at(t_start_test_h)
    Q_end = Q_TES_at(t_end_h)
    if measured_start == "raw":
        print(f"Q_TES at the start of the {mode_test} (t = {t_start_test_h:.2f} h): "
              f"model = {Q_start[0]:.2f} kWh, exp. = {Q_start[1]:.2f} kWh")
    print(f"Q_TES at the end (t = {t_end_h:.2f} h): model = {Q_end[0]:.2f} kWh, exp. = {Q_end[1]:.2f} kWh")
    print(f"Q_TES during the {mode_test} only: model = {Q_end[0]-Q_start[0]:.2f} kWh, "
          f"exp. = {Q_end[1]-Q_start[1]:.2f} kWh")

if not plot_exp:
    exp_data = None

figs = history.plot_all(exp=exp_data)

def autoscale_y_window(ax, x_max, margin=0.05):
    """y limits of ax from the curves and points with x <= x_max (call before adding the vertical lines)."""
    y = []
    for line in ax.get_lines():
        x_d, y_d = np.asarray(line.get_xdata(), float), np.asarray(line.get_ydata(), float)
        y.append(y_d[x_d <= x_max])
    for coll in ax.collections:
        xy = np.asarray(coll.get_offsets(), float)
        if xy.size:
            y.append(xy[xy[:, 0] <= x_max, 1])
    y = np.concatenate(y) if y else np.array([])
    y = y[np.isfinite(y)]
    if y.size:
        dy = margin * max(y.max() - y.min(), 1e-9)
        ax.set_ylim(y.min() - dy, y.max() + dy)


# measured supply from the first row of the file: zoom on the pre-heating, start of the real charge / discharge
if HTF_supply == "measured" and measured_start == "raw":
    keys_time = ["Q_vs_t", "T_vs_t"]

    if zoom_preheating:
        t_zoom_h = t_start_test_h + zoom_after_start_h
        history.skip_exp = skip_exp_zoom
        figs["Q_vs_t_zoom"] = history.plot_Q_vs_t(exp=exp_data)
        figs["T_vs_t_zoom"] = history.plot_T_vs_t(exp=exp_data)
        history.skip_exp = skip_exp_default
        for key in ("Q_vs_t_zoom", "T_vs_t_zoom"):
            for ax in figs[key].axes:
                ax.set_xlim(0, t_zoom_h)
                autoscale_y_window(ax, t_zoom_h)
        keys_time += ["Q_vs_t_zoom", "T_vs_t_zoom"]

    color_start = "#3A3A3A"
    for key in keys_time:
        ax = figs[key].axes[0]
        ax.axvline(t_start_test_h, color=color_start, lw=1.5, ls=":")
        ax.text(t_start_test_h, 1.0, f" start of the {mode_test}",
                transform=ax.get_xaxis_transform(), rotation=90, va="top", ha="right",
                color=color_start, fontsize=14)
    for ax in figs["heat_map"].axes[:2]:
        ax.axhline(t_start_test_h, color="white", lw=1.5, ls=":")

if save_results:
    for key, fig in figs.items():
        fig.savefig(Path(result_folder) / f"{key}-{run_name}.svg", bbox_inches="tight")
    history.save_results(result_folder, run_name)

if show_plots:
    plt.show()
