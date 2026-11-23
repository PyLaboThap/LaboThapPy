# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2026

@author: titouanjanod

Example of the PCM-TES lattice model (TES_lattice.py): simulation of the TES alone with a constant HTF supply
(T_HTF,su, m_dot), optionally the sizing of the number of modules, then all the figures and the saving of the results.

To run it: choose a campaign and the settings in "User choices", then run the file.
The data of each campaign (PCM, HTF supply, number of modules...) are in "Campaigns" further down.
To simulate your own case: copy the template at the top of "Campaigns".
"Advanced" (measured HTF supply of the CEA tests) needs the CEA data files: an ordinary user can skip it.
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

from labothappy.component.storage.TES_lattice import TES_simulation, create_result_folder
from labothappy.component.storage.TES_lattice import skip_exp as skip_exp_default

# =============================================================================
# User choices
# =============================================================================
# Campaign: "Yvernault_charge", "Yvernault_discharge", "CEA_SoC_erythritol", "Zorlu_UC" or "Torrecid_charge_sizing"
# or your own case (data in "Campaigns" below, with a template to copy)
campaign = "CEA_SoC_erythritol"

# Options
do_sizing = False   # True: sizing of the number of modules first (campaigns with sizing data: "Zorlu_UC",
                    # "Torrecid_charge_sizing", or yours)
bypass_on = False   # True: part of the HTF bypasses the TES, to keep the mixed outlet at T_target
heat_losses = True  # losses to the ambient: Q_dot_loss = UA_amb * n_module * (T_PCM_mean - T_amb),
                    # UA_amb of ONE module given in "Campaigns" (0 if unknown). False: no loss
T_amb = 20 + 273.15 # [K] ambient temperature (only used if heat_losses)

# Discretisation
duration_simu_forced = 9*3600 # [s] (!) in SECONDS. None: duration of the campaign (or of the measured test)
dt = 20                       # [s] time step
n_segment_module = 16         # number of cells along one module
n_subsegment_cell = 5         # HTF sub-segments per cell (the PCM stays at one temperature per cell)

SIM = TES_simulation()
SIM.set_parameters(htc_correlation = 'Hausen') # HTF-side heat transfer: 'Gnielinski' (default) or 'Hausen'

# Outputs
show_plots = True
save_results = False # True: figures (svg) + results (npz, txt) in a new folder in result_root
result_root = r"your_path" # your folder
run_name = "Erythritol_700kgh_135C_from80C" # name of the results folder, and suffix of the saved files

# Measurements on the figures (crosses), optional
# exp_data: pandas DataFrame or dict of arrays with the columns
#   t_h [h], Q_dot_kW (HTF side, > 0 in charge and discharge), Q_TES_kWh, T_HTF_ex_C,
#   and optionally T_PCM_min_C, T_PCM_avg_C, T_PCM_max_C (thermocouples). Q_dot_smooth_kW is used if present.
# e.g. exp_data = pandas.read_csv("my_test.csv") with the columns renamed as above.
exp_data = None
plot_exp = False # plot the measurements of the campaign, if it has some ("CEA_SoC_erythritol": needs the CEA files)


# =============================================================================
# Advanced: measured HTF supply of the CEA tests (needs the CEA data files and plots_tests_CEA.py)
# =============================================================================
# HTF_supply = "constant": T_HTF_su and m_dot_tot of "Campaigns" during the whole simulation (works for every campaign).
#              For "CEA_SoC_erythritol", the crosses are the test "2026-02-10 CH (80) 135 700kgh" (same 135 °C, 700 kg/h).
# HTF_supply = "measured": the model is fed at each time step with m_dot(t) and T_HTF,su(t) measured during CEA_test
#              (every 5 s, linearly interpolated), PCM uniform at the measured T_PCM,avg at the start.
#              Only for "CEA_SoC_erythritol", without sizing and without bypass.
HTF_supply = "constant"

# ---- only used if HTF_supply = "measured" ----
CEA_test = "2026-02-10 CH (80) 135 700kgh" # xlsx file of the test (without .xlsx), in DATA_DIR of plots_tests_CEA.py

# Start of the simulation:
#   "raw"  : first row of the file, i.e. the pre-heating of the PCM to 80 °C by the CEA is simulated too
#            (start of the real test = dotted vertical line on the figures vs time)
#   "test" : start of the real test (pre-heating removed, t = 0 as in plots_tests_CEA.py). The real PCM is then
#            not uniform anymore (warmer near the tubes), but it is taken uniform at T_PCM,avg in the model.
measured_start = "test"

# Zoom on the pre-heating (only with "raw"): Q_dot(t) and T(t) from t = 0 to (start of the test + zoom_after_start_h)
zoom_preheating = False
zoom_after_start_h = 0 # [h]
skip_exp_zoom = 20     # one measured point out of skip_exp_zoom on the zoom (1 = every 5 s)

dT_PCM_uniform_max = 5 # [K] warning if T_PCM,max - T_PCM,min > this at the start (uniform initial PCM too rough)


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
# To add YOUR case: add an "elif" block below (e.g. copy "Yvernault_charge"), with a new name in campaign = "...".
# The module is the CEA one (0.5 m x 0.5 m x 1 m, 16 tubes DN10, aluminium Kelvin-cell lattice); another geometry
# can be given with SIM.set_parameters(De_tube=..., thk_tube=..., L_tube=..., V_TES=..., ...), see TES_step.__init__.
#
#   elif campaign == "My_case":
#       fluid = "water"                  # HTF, CoolProp name
#       SIM.set_parameters(              # PCM + lattice (equivalent material)
#               T_PCM = ... + 273.15,    # melting temperature [K]
#               cp_sol = ..., cp_liq = ...,            # [J/(kg.K)]
#               Dh_latent = ...,                       # [J/kg]
#               rho_PCM_sol = ..., rho_PCM_liq = ...,  # [kg/m3]
#               lambda_eq = ...                        # conductivity of PCM + lattice [W/(m.K)], 8 for the CEA module
#               )
#       T_HTF_su = ... + 273.15     # HTF supply temperature [K] (> T_PCM_eq_init: charge, < : discharge)
#       m_dot_tot = ... / 3600      # HTF flow rate entering the whole TES [kg/s]
#       p_HTF = ...e5               # HTF pressure [Pa]
#       T_PCM_eq_init = ... + 273.15 # initial PCM temperature, uniform [K]
#       duration_simu = ...*3600    # [s] (replaced by duration_simu_forced if not None)
#       n_pass_module = 4           # passes in one module (16 tubes: 4 branches of 4 tubes in series)
#       n_module_parallel = 1       # modules in parallel / in series
#       n_module_serie = 1
#       UA_amb_module = 0           # [W/K] heat losses of ONE module (fitted on a standby test), 0 if unknown
#       T_target = None             # [K] only for bypass_on or do_sizing: HTF outlet temperature to keep / to reach
#
#       # only for do_sizing (with T_target above):
#       Q_TES_predesign = ...*3600*1000 # energy to store [J]
#       T_HTF_min = ...             # lowest / highest HTF temperature of the cycle [K]: the stored energy of one
#       T_HTF_max = ...             # module is computed between them, to estimate the number of modules
#
# The measured HTF supply ("Advanced") only works with "CEA_SoC_erythritol".

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
    p_HTF = 5e5 # Pa
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

    # Yvernault et al.: U = 1.06 W/(m2.K) on A = 2.5 m2, fitted on a standby period
    UA_amb_module = 1.06 * 2.5 # W/K


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
    p_HTF = 5e5 # Pa
    T_PCM_eq_init = 80 + 273.15 # K
    T_target = None # no bypass in these tests
    duration_simu = 9*3600 # s

    n_pass_module = 4
    n_module_parallel = 1
    n_module_serie = 1

    # CEA analysis of a 2.5-day standby after a charge (3 periods, "A" = 4.75 to 6 W/K): mean value.
    # Fitted with T_shell - T_ext, here used with T_PCM_mean - T_amb (T_shell not available to us): to be refitted.
    UA_amb_module = 5.3 # W/K

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

    UA_amb_module = 0 # W/K, unknown: no heat loss


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

    UA_amb_module = 0 # W/K, unknown: no heat loss

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

# heat losses: UA_amb of ONE module, the model multiplies it by the number of modules
SIM.set_parameters(UA_amb = UA_amb_module if heat_losses else 0)
T_amb_simu = T_amb if heat_losses else None


# =============================================================================
# Sizing (optional) and simulation
# =============================================================================

if save_results:
    result_folder = create_result_folder(result_root, run_name)
else:
    result_folder = None

if do_sizing:
    # sizing data of the campaign, see the template in "Campaigns"
    missing = [name for name in ("Q_TES_predesign", "T_HTF_min", "T_HTF_max") if name not in globals()]
    if T_target is None:
        missing.append("T_target")
    if missing:
        raise ValueError(f"do_sizing: missing sizing data for the campaign '{campaign}': {', '.join(missing)}")

    best_config = SIM.TES_sizing(T_target, Q_TES_predesign, T_HTF_min, T_HTF_max, duration_simu, dt, m_dot_tot,
                                 T_HTF_su, p_HTF, T_PCM_eq_init, fluid, bypass_on, T_amb=T_amb_simu)
    print("best_config = ", best_config)

    SIM.set_parameters(n_box_serie = best_config["N_series"], n_box_parallel = best_config["N_parallel"])

    if save_results:
        with open(result_folder / f"best_config-{run_name}.txt", "w") as f:
            for key in ("N_parallel", "N_series", "N_box", "margin", "index", "target_reached"):
                f.write(f"{key} = {best_config[key]}\n")

if HTF_supply == "measured":
    # m_dot(t) and T_HTF,su(t) of the test, interpolated at each time step
    history = SIM.simulate_measured(exp_data["t_s"], exp_data["m_dot_kg_s"], exp_data["T_HTF_su_C"] + 273.15,
                                    dt, p_HTF, T_PCM_eq_init, fluid, duration_simu, T_amb=T_amb_simu)
    history.title = title_measured
else:
    history = SIM.simulate(duration_simu, dt, m_dot_tot, T_HTF_su, p_HTF, T_PCM_eq_init, fluid, T_target, bypass_on,
                           T_amb=T_amb_simu)


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
print(f"Q_tot = {history.Q_tot/3.6e6:.2f} kWh (given by the HTF)")
if SIM.TES.UA_amb > 0:
    print(f"Q_loss = {history.Q_loss_tot/3.6e6:.2f} kWh (lost to the ambient, UA_amb = {SIM.TES.UA_amb} W/K per module, "
          f"T_amb = {T_amb-273.15:.1f} °C)")
    print(f"Q_stored = Q_tot - Q_loss = {(history.Q_tot-history.Q_loss_tot)/3.6e6:.2f} kWh")
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
