# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2026

@author: titouanjanod

Example of the PCM-TES lattice model (TES_lattice.py): simulation of the TES alone with a constant HTF supply
(T_HTF,su, m_dot), optionally the sizing of the number of modules, then all the figures and the saving of the results.

To run it: choose a campaign and the settings in "User choices", then run the file.
The data of each campaign (PCM, HTF supply, number of modules...) are in "Campaigns" further down.
To simulate your own case: copy the template at the top of "Campaigns".
"""

import os
import sys
import matplotlib.pyplot as plt
from pathlib import Path

# Absolute path to the LaboThapPy folder (4 levels above this file), if labothappy is not installed
LABOTHAPPY_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
if LABOTHAPPY_PATH not in sys.path:
    sys.path.append(LABOTHAPPY_PATH)

from labothappy.component.storage.TES_lattice import TES_simulation, create_result_folder

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
duration_simu_forced = None # [s] (!) in SECONDS. None: duration of the campaign
dt = 20                       # [s] time step
n_segment_module = 16         # number of cells along one module
n_subsegment_cell = 5         # HTF sub-segments per cell (the PCM stays at one temperature per cell)

SIM = TES_simulation()
SIM.set_parameters(htc_correlation = 'Hausen') # HTF-side heat transfer: 'Gnielinski' (default) or 'Hausen'

# Outputs
show_plots = True
save_results = False # True: figures (svg) + results (npz, txt) in a new folder in result_root
result_root = r"C:\Users\titou\OneDrive - Universite de Liege\!SEHRENE\TES modelling\Results" # your folder
run_name = "Yvernault_charge_newMASSES" # name of the results folder, and suffix of the saved files

# Measurements on the figures (crosses), optional
# exp_data: pandas DataFrame or dict of arrays with the columns
#   t_h [h], Q_dot_kW (HTF side, > 0 in charge and discharge), Q_TES_kWh, T_HTF_ex_C,
#   and optionally T_PCM_min_C, T_PCM_avg_C, T_PCM_max_C (thermocouples). Q_dot_smooth_kW is used if present.
# e.g. exp_data = pandas.read_csv("my_test.csv") with the columns renamed as above. None: model only.
# Campaigns "Yvernault_charge" / "Yvernault_discharge": if None, the digitized measurements of the paper are loaded
# (titouanjanod_perso/TES_related/Yvernault_to_dataframe.py).
exp_data = None
skip_exp = 1           # only one measured point out of skip_exp is plotted (1: all; e.g. 60 for data every 5 s)
dt_grid_Yvernault = 300 # [s] time step of the measurements of Yvernault after interpolation (None: digitized times)


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
            lambda_eq = 8,
            lambda_tube = 45
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

    # measurements of the paper (digitized with WebPlotDigitizer), interpolated on a common time grid
    if exp_data is None:
        TES_RELATED_DIR = r"C:\Users\titou\Documents\GitHub\titouanjanod_perso\TES_related"
        if TES_RELATED_DIR not in sys.path:
            sys.path.append(TES_RELATED_DIR)
        from Yvernault_to_dataframe import load_Yvernault
        exp_data = load_Yvernault("charge" if campaign == "Yvernault_charge" else "discharge",
                                  dt_grid_s=dt_grid_Yvernault)


elif campaign == "CEA_SoC_erythritol":
    # CEA Grenoble campaign for the State of Charge (SoC) measurement (B. Champel), SEHRENE deliverable D4.7 https://doi.org/10.5281/zenodo.24496277
    # (nominal values of the test "2026-02-10 CH (80) 135 700kgh")
    fluid = "water"
    SIM.set_parameters(
            T_PCM = 118.9 + 273.15,
            cp_sol = 1.637e3,
            cp_liq = 2.846e3,
            Dh_latent = 335.5e3,
            rho_PCM_sol = 1448.4,
            rho_PCM_liq = 1281.5,
            lambda_eq = 8.7
            )
    T_HTF_su = 135 + 273.15 # K
    m_dot_tot = 700 / 3600 # kg/s
    p_HTF = 5e5 # Pa
    T_PCM_eq_init = 80 + 273.15 # K
    T_target = None # no bypass in these tests
    duration_simu = 9*3600 # s

    n_pass_module = 4
    n_module_parallel = 1
    n_module_serie = 1

    # fitted on a standby of the module (solid PCM, quasi-steady), with T_PCM_mean - T_amb
    UA_amb_module = 5.4 # W/K


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

history.skip_exp = skip_exp
figs = history.plot_all(exp=exp_data)

if save_results:
    for key, fig in figs.items():
        fig.savefig(Path(result_folder) / f"{key}-{run_name}.svg", bbox_inches="tight")
    history.save_results(result_folder, run_name, exp=exp_data)

if show_plots:
    plt.show()
