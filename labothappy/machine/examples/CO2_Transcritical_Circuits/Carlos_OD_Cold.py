# -*- coding: utf-8 -*-
"""
CO2 Transcritical Cycle — Off-Design Study  (6-map version - COLD SOURCE SWEPT)
===============================================================================
Generates SIX separate figures:

  Series A  — T_su_w_cd swept, m_dot_w_cd = 4542 kg/s (nominal)
      A1 : T_su_w_cd = 0.1 °C   (nominal)
      A2 : T_su_w_cd = 5.0 °C   
      A3 : T_su_w_cd = 10.0 °C

  Series B  — m_dot_w_cd swept, T_su_w_cd = 0.1 °C (nominal)
      B1 : m_dot_w_cd = 4542 × 0.85 = 3860.7 kg/s
      B2 : m_dot_w_cd = 4542        = 4542.0 kg/s  (nominal)
      B3 : m_dot_w_cd = 4542 × 1.15 = 5223.3 kg/s

Hot source is kept constant at nominal: T_su_w_gh = 150 °C, m_dot_w_gh = 192.24 kg/s.

Each figure contains three contour maps:
    1. Net thermal efficiency   η  [ - ]
    2. Net electrical power     W_net  [kW]
    3. High-side pressure       P_high [bar]
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from joblib import Parallel, delayed
from tqdm import tqdm
import os
from contextlib import redirect_stdout, redirect_stderr

from labothappy.machine.circuit_fpi import CircuitFPI
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.expander.turbine_mean_line_Aungier import AxialTurbineMeanLine
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity

os.environ["NUMEXPR_MAX_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"]     = "1"
os.environ["MKL_NUM_THREADS"]     = "1"
os.environ["OPENBLAS_NUM_THREADS"]= "1"

# =============================================================================
# GLOBAL EFFICIENCY CORRECTION FACTORS
# =============================================================================
ETA_MOTOR    = 0.95
ETA_SEC_PUMP = 0.80
SC_TARGET    = 0.5   

MAP_SAVE_PATH = r"C:\Users\Basile\Desktop\Travail\Thèse\Travail\WP1\Turbomachines\Save Maps\turb_map.parquet"

# =============================================================================
# CONSTANT HOT SOURCE NOMINAL VALUES
# =============================================================================
T_GH_NOM_K  = 150.0 + 273.15
MDOT_GH_NOM = 192.24


# =============================================================================
# HELPER — net power & efficiency
# =============================================================================
def compute_net_performance(W_turb_shaft, W_pump_shaft, Q_gh):
    W_turb_elec = W_turb_shaft * ETA_MOTOR
    W_pump_elec = W_pump_shaft / (ETA_SEC_PUMP * ETA_MOTOR)
    W_net       = W_turb_elec - W_pump_elec
    eta         = W_net / Q_gh if Q_gh > 0 else np.nan
    return W_net, eta


# =============================================================================
# CIRCUIT BUILDER  (T_su_w_cd and m_dot_w_cd are now parameters)
# =============================================================================
def CO2_OD_TC(mdot_CO2_guess, P_high_guess,
              N_pp, N_exp,
              T_su_w_cd,
              m_dot_w_cd,
              SC_target=SC_TARGET,
              damping_it=0.3,
              mute_print=1):
    """
    Build and configure the CO2 transcritical cycle circuit.
    """
    fluid  = 'CO2'
    CO2_TC = CircuitFPI(fluid)
    if mute_print:
        CO2_TC.mute_print()

    # ------------------------------------------------------------------
    # Components
    # ------------------------------------------------------------------
    Turbine     = AxialTurbineMeanLine(fluid)
    GasHeater   = HexMBChargeSensitive('Shell&Tube')
    Recuperator = HexMBChargeSensitive('PCHE')
    Pump        = PumpCurveSimilarity()
    Condenser   = HexMBChargeSensitive('Shell&Tube')

    # ------------------------------------------------------------------
    # Turbine parameters (Truncated for brevity, kept exactly as yours)
    # ------------------------------------------------------------------
    Turb_params = {
        'mdot_rated'                : 318.437021666738,
        'Wdot_rated'                : 15244151.612636594,
        'N_rot_rated'               : 2864.775003723627,
        'total_to_static_efficiency': 0.8851473369105862,
        'DP_rated'                  : 2.93,
        'n_stages'                  : 7,
        'p0_su'                     : 15309670.5,
        'T0_su'                     : 406.4,
        'p_ex'                      : 5220928,
        'r_m'                       : 0.20956432565203412,
        'delta_tip'                 : 0.0004,
        'N_lw'                      : 0,
        'D_lw'                      : 0,
        'e_blade'                   : 2e-06,
        'stator': {
            'h_blade_S': [0.024581131784973807, 0.02731115991390685, 0.030535241526643685, 0.0343428949996286, 0.03884110173212556,  0.04415656370209263, 0.05044235012970171,  0.05396446226956401],
            'chord_S'  : [0.017387241421509592, 0.018189646227485024, 0.019157515021682635, 0.020309513815586707, 0.021666371275761282, 0.023250883041521515, 0.025089330758218745, 0.026110623182803244],
            'xhi_S1'   : [-0.5493774236565352]*8,
            'xhi_S2'   : [1.1268112059950612]*8,
            'pitch_S'  : [0.014118294382457852, 0.014769840363216227, 0.015555741716329778, 0.01649115508675976, 0.01759291197811346,  0.01887952225855718, 0.020372326403956538, 0.021201607297417914],
            'o_S'      : [0.00407007938170583,  0.004257909709554279, 0.00448447255112079,  0.004754137325719167, 0.005071756287732223, 0.0054426655373106475, 0.005873017193716411, 0.006112085665777279],
            't_TE_S'   : [0.0005]*8,
            't_blade_S': [0.0017387241421509593, 0.0018189646227485025, 0.0019157515021682636, 0.002030951381558671, 0.0021666371275761284, 0.0023250883041521517, 0.0025089330758218745, 0.0026110623182803248],
            'n_blade_S': [93, 89, 85, 80, 75, 70, 65, 62],
            'R_c_S'    : [0.02593658923532629,  0.02713353838605762, 0.02857731056016104,  0.03029574989120514, 0.03231977738022937,  0.034683397336424615, 0.03742581415244761,  0.03894927848261039],
        },
        'rotor': {
            'h_blade_R': [0.02599623192810525,  0.028975662485097392, 0.03249353115560146,  0.03664835013019067, 0.041556335411602884, 0.04735747516384924, 0.054218932168004906, None],
            'chord_R'  : [0.017804442744401783, 0.018691926981905842, 0.019754111574023726, 0.0210108023910845, 0.02248355410347207,  0.02419698360942641, 0.02617885482072475,  None],
            'xhi_R1'   : [0.5033968740474793]*7 + [None],
            'xhi_R2'   : [-1.1544956444218015]*7 + [None],
            'pitch_R'  : [0.015148494029818542, 0.015903589253319174, 0.01680732419089942,  0.017876550204477332, 0.01912960658162921,  0.02058743803494024, 0.022273666839919944, None],
            'o_R'      : [0.004776675516045738, 0.005014774752793548, 0.005299743578132708, 0.005636895621769747, 0.006032013691274512, 0.006491702145839364, 0.007023409643056108, None],
            't_TE_R'   : [0.0005]*7 + [None],
            't_blade_R': [0.0019275069515591894, 0.0020395300073626607, 0.002174650407702253,  0.0022904090712424882, 0.0024358345376652743, 0.002590844225248941, 0.002758858045024202,  None],
            'n_blade_R': [87, 83, 78, 74, 69, 64, 59, None],
            'R_c_R'    : [0.02655892943743007,  0.02788279177780649, 0.029467255051191132, 0.03134186372129275, 0.03353877094100647,  0.036094697795777626, 0.03905105977860335,  None],
        },
    }

    Turbine.set_parameters(
        r_m        = Turb_params['r_m'],
        nStages    = Turb_params['n_stages'],
        mdot_rated = Turb_params['mdot_rated'],
        DP_rated   = Turb_params['DP_rated'],
        damping    = 0.2,
        delta_tip  = Turb_params['delta_tip'],
        N_lw       = Turb_params['N_lw'],
        D_lw       = Turb_params['D_lw'],
        e_blade    = Turb_params['e_blade'],
        solve_type = 'map',
        map_mode   = "P_N",
    )
    
    Turbine.set_stage_parameters(**Turb_params['stator'], **Turb_params['rotor'])
    Turbine.load_map_df(MAP_SAVE_PATH)

    # ------------------------------------------------------------------
    # Gas Heater
    # ------------------------------------------------------------------
    params_GH = {
        'n_series': 1, 'n_parallel': 3,
        'foul_t': 0.000176, 'foul_s': 0.000176,
        'tube_cond': 20, 'Shell_Side': 'H', 'Flow_Type': 'Shell&Tube',
        'n_disc': 20, 'Shell_ID': 1.8288, 'Tube_L': 5.832,
        'Tube_OD': 0.009525, 'Tube_t': 0.0005588,
        'central_spacing': 0.37, 'Tube_pass': 1,
        'n_tubes': 16164, 'pitch_ratio': 1.33,
        'tube_layout': 0, 'Baffle_cut': 34.973,
    }
    GasHeater.set_parameters(**{k: v for k, v in params_GH.items()})
    GasHeater.set_htc(htc_type="Correlation_Disc",
                      Corr_H={"SC": "Shell_Kern_HTC", "1P": "Shell_Kern_HTC", "2P": "Shell_Kern_HTC"},
                      Corr_C={"SC": "Gnielinski",     "1P": "Gnielinski",     "2P": "Flow_boiling"})
    GasHeater.set_DP(DP_type="Correlation_Disc",
                     Corr_H={"SC": "Shell_Kern_DP",  "1P": "Shell_Kern_DP",  "2P": "Shell_Kern_DP"},
                     Corr_C={"SC": "Gnielinski_DP",  "1P": "Gnielinski_DP",  "2P": "Choi_DP"})

    # ------------------------------------------------------------------
    # Condenser
    # ------------------------------------------------------------------
    params_CD = {
        'n_series': 1, 'n_parallel': 10.0,
        'foul_t': 0.000176, 'foul_s': 0.000176,
        'tube_cond': 20, 'Shell_Side': 'C', 'Flow_Type': 'Shell&Tube',
        'n_disc': 20, 'Shell_ID': 1.6764, 'Tube_L': 6.8,
        'Tube_OD': 0.0254, 'Tube_t': 0.0012446,
        'central_spacing': 1.7, 'Tube_pass': 1.0,
        'cross_passes': 3, 'n_tubes': 2106,
        'pitch_ratio': 1.25, 'tube_layout': 0.0, 'Baffle_cut': 18.0,
    }
    Condenser.set_parameters(**{k: v for k, v in params_CD.items()})
    Condenser.set_htc(htc_type="Correlation_Disc",
                      Corr_H={"SC": "Gnielinski",    "1P": "Gnielinski",    "2P": "Thome_Condensation"},
                      Corr_C={"SC": "Shell_Kern_HTC","1P": "Shell_Kern_HTC","2P": "Shell_Kern_HTC"})
    Condenser.set_DP(DP_type="Correlation_Disc",
                     Corr_H={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
                     Corr_C={"SC": "Shell_Kern_DP", "1P": "Shell_Kern_DP", "2P": "Shell_Kern_DP"})

    # ------------------------------------------------------------------
    # Pump
    # ------------------------------------------------------------------
    D_H_rated   = 1211
    V_dot_rated = 1298
    eta_rated   = 0.867
    N_pp_rated  = 2900

    V_dot_curve   = np.array([76.2, 91.0, 113.7, 136.3, 158.9, 181.8, 204.7, 227.6, 250.2, 272.9, 295.5]) * V_dot_rated / 250.2
    eta_is_curve  = np.array([40.4, 47.2,  56.0,  62.6,  68.4,  72.8, 76.1,  78.5,  79.2,  78.7,  76.8]) * eta_rated / 79.2
    Delta_H_curve = np.array([1406.5, 1394.7, 1370.9, 1342.4, 1309.0, 1270.8, 1232.6, 1177.6, 1105.7, 1014.6,  918.7]) * D_H_rated / 1105.7
    NPSH_r_curve  = np.array([1.1, 1.1, 1.25, 1.4, 1.6, 1.8, 1.9, 2.0, 3.0, 3.85, 4.7])

    Pump.set_parameters(
        V_dot_curve   = V_dot_curve,
        Delta_H_curve = Delta_H_curve,
        eta_is_curve  = eta_is_curve,
        NPSH_r_curve  = NPSH_r_curve,
        N_rot_rated   = N_pp_rated,
        mode          = "M_N",
    )

    # ------------------------------------------------------------------
    # Recuperator
    # ------------------------------------------------------------------
    params_REC = {
        'alpha': 32.62, 'D_c': 2.42e-3,
        'C_V_tot': 1, 'H_V_tot': 1, 'k_cond': 60,
        'L_c': 0.7432303013776589, 'N_c': 736, 'N_p': 563,
        'R_p': 1, 't_2': 0.0012282802564224898,
        't_3': 0.0009428803890487963,
    }
    Recuperator.set_parameters(**params_REC, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=20)
    Recuperator.set_htc(htc_type='Correlation_Disc',
                        Corr_H={"1P": "Gnielinski", "SC": "Gnielinski", "2P": "Thome_Condensation"},
                        Corr_C={"1P": "Gnielinski", "SC": "Gnielinski", "2P": "Flow_boiling"})
    Recuperator.set_DP(DP_type="Correlation_Disc",
                       Corr_H={"SC": "Darcy_Weisbach", "1P": "Darcy_Weisbach", "2P": "Choi_DP"},
                       Corr_C={"SC": "Darcy_Weisbach", "1P": "Darcy_Weisbach", "2P": "Choi_DP"})

    # ------------------------------------------------------------------
    # Assemble circuit
    # ------------------------------------------------------------------
    CO2_TC.add_component(Turbine,     "Turbine")
    CO2_TC.add_component(Condenser,   "Condenser")
    CO2_TC.add_component(Pump,        "Pump")
    CO2_TC.add_component(GasHeater,   "GasHeater")
    CO2_TC.add_component(Recuperator, "Recuperator")

    CO2_TC.link_components("Turbine",    "m-ex",   "Recuperator", "m-su_H")
    CO2_TC.link_components("Recuperator","m-ex_H", "Condenser",   "m-su_H")
    CO2_TC.link_components("Condenser",  "m-ex_H", "Pump",        "m-su")
    CO2_TC.link_components("Pump",       "m-ex",   "Recuperator", "m-su_C")
    CO2_TC.link_components("Recuperator","m-ex_C", "GasHeater",   "m-su_C")
    CO2_TC.link_components("GasHeater",  "m-ex_C", "Turbine",     "m-su")

    # ------------------------------------------------------------------
    # Boundary conditions  ← T_su_w_cd and m_dot_w_cd are injected here
    # ------------------------------------------------------------------
    CD_source  = MassConnector('Water')
    GH_source  = MassConnector('Water')

    CO2_TC.add_source("CD_Water", CD_source,
                      CO2_TC.components["Condenser"], "m-su_C")
    CO2_TC.set_source_properties(T=T_su_w_cd, fluid='Water',      # ← parametric
                                 P=5e5, m_dot=m_dot_w_cd,         # ← parametric
                                 target="CD_Water")

    CO2_TC.add_source("GH_Water", GH_source,
                      CO2_TC.components["GasHeater"], "m-su_H")
    CO2_TC.set_source_properties(T=T_GH_NOM_K, fluid='Water',     # ← CONSTANT NOMINAL
                                 P=5e5, m_dot=MDOT_GH_NOM,        # ← CONSTANT NOMINAL
                                 target="GH_Water")

    # ------------------------------------------------------------------
    # Component inputs
    # ------------------------------------------------------------------
    Pump.set_inputs(fluid="CO2", N_rot=N_pp)
    Turbine.set_inputs(N_rot=N_exp)

    # ------------------------------------------------------------------
    # Initial guesses (Updated for parametric cold side)
    # ------------------------------------------------------------------
    P_LP_guess     = PropsSI("P", "T", T_su_w_cd + 10, "Q", 0, 'CO2')
    T_sat_LP       = PropsSI("T", "P", P_LP_guess, "Q", 0.5, 'CO2')
    h_SC_guess     = PropsSI("H", "P", P_LP_guess, "T", T_sat_LP - 5, 'CO2')
    h_su_exp_guess = PropsSI("H", "P", P_high_guess, "T", T_GH_NOM_K - 10, 'CO2')

    P_LP_min = PropsSI("P", "T", T_su_w_cd + 1, "Q", 0, 'CO2')

    CO2_TC.set_cycle_guess(target="Pump:su",
                            m_dot=mdot_CO2_guess, h=h_SC_guess, p=P_LP_guess)
    CO2_TC.set_cycle_guess(target="Pump:ex",    p=P_high_guess)
    CO2_TC.set_cycle_guess(target="Turbine:su",
                            m_dot=mdot_CO2_guess, p=P_high_guess, h=h_su_exp_guess)
    CO2_TC.set_cycle_guess(target="Turbine:ex", p=P_LP_guess)

    CO2_TC.set_iteration_variable(
        it_var        = 'Turbine:ex-p',
        objective     = 'Pump:su-SC',
        target_value  = SC_target,
        obj_type      = "Target_val",
        damping_factor= damping_it,
    )

    return CO2_TC


# =============================================================================
# PLOTTING — paper-quality contour maps
# =============================================================================
def plot_and_save(eta_2D, Wnet_2D, Phigh_2D,
                  N_pp_sweep, N_exp_sweep,
                  N_PP_NOMINAL, N_EXP_NOMINAL,
                  title_top, filename):
    """
    Draw three filled-contour maps and save to <filename>.
    Typography and layout are tuned for journal / paper inclusion.
    """
    N_exp_grid, N_pp_grid = np.meshgrid(N_exp_sweep, N_pp_sweep)

    FS_sup    = 20
    FS_title  = 17
    FS_label  = 14
    FS_tick   = 12
    FS_cbar   = 13
    FS_clabel = 12
    FS_legend = 12

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    fig.suptitle(title_top, fontsize=FS_sup, fontweight='bold', y=1.01)

    maps_cfg = [
        {'data': eta_2D,   'title': 'Net Thermal Efficiency',
         'label': r'$\eta_\mathrm{net}$  [ - ]',
         'cmap': 'plasma',   'fmt': '%.3f', 'n_levels': 20, 'n_lines': 8},
        {'data': Wnet_2D,  'title': 'Net Electrical Power',
         'label': r'$\dot{W}_\mathrm{net}$  [kW]',
         'cmap': 'RdYlGn_r',   'fmt': '%.0f', 'n_levels': 20, 'n_lines': 8},
        {'data': Phigh_2D, 'title': 'High-Side Pressure',
         'label': r'$P_\mathrm{high}$  [bar]',
         'cmap': 'viridis',  'fmt': '%.1f', 'n_levels': 20, 'n_lines': 8},
    ]

    for ax, cfg in zip(axes, maps_cfg):
        data = cfg['data']

        cf = ax.contourf(N_exp_grid, N_pp_grid, data,
                         levels=cfg['n_levels'], cmap=cfg['cmap'])
        cbar = fig.colorbar(cf, ax=ax, shrink=0.88, pad=0.03)
        cbar.set_label(cfg['label'], fontsize=FS_cbar, labelpad=6)
        cbar.ax.tick_params(labelsize=FS_tick)
        cbar.locator   = ticker.MaxNLocator(nbins=6)
        cbar.update_ticks()

        cs = ax.contour(N_exp_grid, N_pp_grid, data,
                        levels=cfg['n_lines'], colors='white',
                        linewidths=0.7, alpha=0.65)
        ax.clabel(cs, inline=True, fontsize=FS_clabel,
                  fmt=cfg['fmt'], use_clabeltext=True)

        ax.plot(N_EXP_NOMINAL, N_PP_NOMINAL,
                'r*', markersize=11, zorder=5,
                label=f'Nominal\n({N_EXP_NOMINAL:.0f}, {N_PP_NOMINAL:.0f}) RPM')

        ax.set_xlabel(r'$N_\mathrm{exp}$  [RPM]', fontsize=FS_label, labelpad=5)
        ax.set_ylabel(r'$N_\mathrm{pp}$  [RPM]',  fontsize=FS_label, labelpad=5)
        ax.set_title(cfg['title'], fontsize=FS_title, fontweight='bold', pad=8)

        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5, integer=True))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, integer=True))
        ax.tick_params(axis='both', labelsize=FS_tick, which='both',
                       direction='in', length=4)

        ax.legend(fontsize=FS_legend, loc='lower right',
                  framealpha=0.7, edgecolor='white')

    plt.tight_layout(rect=[0, 0, 1, 1])
    # plt.savefig(filename, dpi=180, bbox_inches='tight')
    plt.show()
    
    plt.close(fig)
    print(f"  Saved → {filename}")


# =============================================================================
# OFF-DESIGN SWEEP
# =============================================================================
if __name__ == "__main__":

    # ---- Reference / nominal values -----------------------------------
    MDOT_GUESS       = 280
    P_HIGH_GUESS_NOM = 120.0e5   # Pa — nominal high-side pressure guess

    N_PP_NOMINAL  = 2900.0
    N_EXP_NOMINAL = 2864.775

    T_CD_NOM    = 0.1 + 273.15   # K
    MDOT_CD_NOM = 4542.0         # kg/s

    # ---- Grid definition  (±15 %, 20 × 20) ----------------------------
    N_POINTS = 5
    N_pp_sweep  = np.linspace(0.85 * N_PP_NOMINAL,  1 * N_PP_NOMINAL,  N_POINTS)
    N_exp_sweep = np.linspace(0.85 * N_EXP_NOMINAL, 1.15 * N_EXP_NOMINAL, N_POINTS)

    # ---- 6 study cases (COLD SOURCE SWEPT) ----------------------------
    #  Series A: vary T_su_w_cd  (m_dot_w_cd = MDOT_CD_NOM)
    #  Series B: vary m_dot_w_cd (T_su_w_cd  = T_CD_NOM)
    CASES = [
        # tag,  T_cd [K],               mdot_cd [kg/s],       filename
        ("A4", 15 + 273.15,            MDOT_CD_NOM,          "CO2_OD_A4_Tcd_15.png"),
        # ("A5", 20.0 + 273.15,            MDOT_CD_NOM,          "CO2_OD_A5_Tcd_20.png"),
        #("A3", 10.0 + 273.15,           MDOT_CD_NOM,          "CO2_OD_A3_Tcd_10_0.png"),
        #("B1", T_CD_NOM,                MDOT_CD_NOM * 0.85,   "CO2_OD_B1_mdot_cd_-15pct.png"),
        #("B2", T_CD_NOM,                MDOT_CD_NOM,          "CO2_OD_B2_mdot_cd_nom.png"),
        #("B3", T_CD_NOM,                MDOT_CD_NOM * 1.15,   "CO2_OD_B3_mdot_cd_+15pct.png"),
    ]

    # ---- Parallel single-point solver ---------------------------------
    def _try_solve(N_pp, N_exp, T_cd, mdot_cd, P_high_guess,
                   damping, max_iter, verbose=False):
        """
        Build a fresh circuit with the given damping and solve it.
        Returns the circuit object if converged, None otherwise.
        """
        cyc = CO2_OD_TC(
            mdot_CO2_guess = MDOT_GUESS,
            P_high_guess   = P_high_guess,
            N_pp           = N_pp,
            N_exp          = N_exp,
            T_su_w_cd      = T_cd,
            m_dot_w_cd     = mdot_cd,
            SC_target      = SC_TARGET,
            damping_it     = damping,
            mute_print     = 1, # 0 if verbose else 1,
        )
        # cyc.solve(method='wegstein', max_iter=max_iter, tol=1e-3)
        cyc.solve(method='wegstein', max_iter=max_iter, tol=1e-3)
        return cyc if cyc.converged else None

    # Damping ladder: start gentle, escalate if needed
    DAMPING_LADDER = [0.1, 0.3, 0.5]
    MAX_ITER       = 80

    def evaluate_point(i, j, N_pp, N_exp, T_cd, mdot_cd, P_high_guess):
        nan_row = (i, j, np.nan, np.nan, np.nan)
        
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull), redirect_stderr(fnull):
                cyc = None
                for damp in DAMPING_LADDER:
                    # try:
                        cyc = _try_solve(N_pp, N_exp, T_cd, mdot_cd,
                                         P_high_guess, damp, MAX_ITER)
                        
                        if cyc is not None:
                            break   # converged — stop retrying
                            
                    # except Exception:
                    #     cyc = None
                    #     continue   # try next damping level

        if cyc is None:
            print(f"  [NC all dampings] N_pp={N_pp:.0f} | N_exp={N_exp:.0f}")
            return nan_row

        W_turb = cyc.components['Turbine'].model.W.W_dot
        W_pump = cyc.components['Pump'].model.W.W_dot
        Q_gh   = cyc.components['GasHeater'].model.Q.Q_dot
        P_high = cyc.components['Turbine'].model.su.p

        W_net, eta = compute_net_performance(W_turb, W_pump, Q_gh)
        return i, j, eta, W_net / 1e3, P_high / 1e5

    def diagnose_point(T_cd, mdot_cd, P_high_guess, N_pp, N_exp):
        """
        Solve ONE point in series trying each damping level in turn.
        Full verbosity on the first attempt that succeeds (or the last that fails).
        """
        T_cd_C = T_cd - 273.15
        print(f"\n{'#'*65}")
        print(f"  DIAGNOSTIC POINT")
        print(f"  T_cd      = {T_cd_C:.1f} °C  ({T_cd:.2f} K)")
        print(f"  m_dot_cd  = {mdot_cd:.2f} kg/s")
        print(f"  P_high_0  = {P_high_guess/1e5:.1f} bar  (initial guess)")
        print(f"  N_pp      = {N_pp:.0f} RPM  |  N_exp = {N_exp:.0f} RPM")
        print(f"  SC_TARGET = {SC_TARGET} K")
        print(f"  Damping ladder: {DAMPING_LADDER}")
        print(f"{'#'*65}\n")

        cyc = None
        for damp in DAMPING_LADDER:
            print(f"  --- Trying damping = {damp} ---")
            try:
                cyc = _try_solve(N_pp, N_exp, T_cd, mdot_cd,
                             P_high_guess, damp, MAX_ITER,
                             verbose=True)
            except Exception as exc:
                import traceback
                print(f"  *** EXCEPTION (damping={damp}): "
                      f"{type(exc).__name__}: {exc}")
                traceback.print_exc()
                cyc = None

            if cyc is not None:
                print(f"\n  ✓ Converged at damping = {damp}")
                break
            else:
                print(f"  ✗ Did not converge at damping = {damp}\n")

        if cyc is None:
            print("\n  ✗ All damping levels failed.")
        else:
            W_turb = cyc.components['Turbine'].model.W.W_dot
            W_pump = cyc.components['Pump'].model.W.W_dot
            Q_gh   = cyc.components['GasHeater'].model.Q.Q_dot
            P_high = cyc.components['Turbine'].model.su.p
            W_net, eta = compute_net_performance(W_turb, W_pump, Q_gh)
            print(f"  η_net  = {eta:.4f}")
            print(f"  W_net  = {W_net/1e3:.1f} kW")
            print(f"  P_high = {P_high/1e5:.2f} bar")
        print(f"{'#'*65}\n")

    tasks = [(i, j, Npp, Nex)
             for i, Npp  in enumerate(N_pp_sweep)
             for j, Nex  in enumerate(N_exp_sweep)]

    # ---- Main loop over the 6 cases -----------------------------------
    for tag, T_cd, mdot_cd, fname in CASES:

        T_cd_C = T_cd - 273.15

        # P_HIGH_GUESS remains mostly dictated by the constant hot source in a Transcritical cycle.
        # We can safely use the nominal P_HIGH_GUESS as the starting point here.
        P_HIGH_GUESS = P_HIGH_GUESS_NOM
        print(f"\n{'='*65}")
        print(f"  Case {tag}: T_cd = {T_cd_C:.1f} °C | ṁ_cd = {mdot_cd:.2f} kg/s")
        print(f"  P_HIGH_GUESS = {P_HIGH_GUESS/1e5:.1f} bar")
        print(f"  Grid: {N_POINTS}×{N_POINTS} = {len(tasks)} points")
        print(f"{'='*65}")

        print(f"\n  >>> Running diagnostic at nominal (N_pp={N_PP_NOMINAL:.0f}, "
              f"N_exp={N_EXP_NOMINAL:.0f}) RPM ...")
        diagnose_point(T_cd, mdot_cd, P_HIGH_GUESS,
                       N_PP_NOMINAL, N_EXP_NOMINAL)
        print("  >>> Diagnostic done. Launching full sweep ...\n")

        results = Parallel(n_jobs=-2, backend='loky', verbose=0)(
            delayed(evaluate_point)(i, j, Npp, Nex, T_cd, mdot_cd, P_HIGH_GUESS)
            for i, j, Npp, Nex in tqdm(tasks, desc=f"Case {tag}", unit="pt")
        )
        
        print("A")
        
        eta_2D   = np.full((N_POINTS, N_POINTS), np.nan)
        Wnet_2D  = np.full((N_POINTS, N_POINTS), np.nan)
        Phigh_2D = np.full((N_POINTS, N_POINTS), np.nan)

        for i, j, eta, W_net_kW, P_high_bar in results:
            eta_2D[i, j]   = eta
            Wnet_2D[i, j]  = W_net_kW
            Phigh_2D[i, j] = P_high_bar

        print("B")

        # Determine which quantity varies for the subtitle
        if tag.startswith("A"):
            subtitle = (
                f"T$_{{cold}}$ = {T_cd_C:.1f} °C  |  "
                r"$\dot{m}_{cold}$" + f" = {mdot_cd:.2f} kg/s  |  "
                f"$N_{{pp}}$ & $N_{{exp}}$ ±15 %"
            )
        else:
            subtitle = (
                f"T$_{{cold}}$ = {T_cd_C:.1f} °C  |  "
                r"$\dot{m}_{cold}$" + f" = {mdot_cd:.2f} kg/s  ({mdot_cd/MDOT_CD_NOM*100:.0f} %)  |  "
                f"$N_{{pp}}$ & $N_{{exp}}$ ±15 %"
            )

        print("C")

        title_top = (
            f"CO$_2$ Transcritical Cycle — Off-Design Map  (Case {tag})\n"
            + subtitle
        )

        plot_and_save(
            eta_2D, Wnet_2D, Phigh_2D,
            N_pp_sweep, N_exp_sweep,
            N_PP_NOMINAL, N_EXP_NOMINAL,
            title_top, fname,
        )
        
        print("D")

        # Quick summary at nearest nominal node
        i_nom = np.argmin(np.abs(N_pp_sweep  - N_PP_NOMINAL))
        j_nom = np.argmin(np.abs(N_exp_sweep - N_EXP_NOMINAL))
        print(f"  Nominal node → η={eta_2D[i_nom,j_nom]:.4f} | "
              f"W_net={Wnet_2D[i_nom,j_nom]:.1f} kW | "
              f"P_HP={Phigh_2D[i_nom,j_nom]:.2f} bar")

    print("\n✓  All 6 maps saved.")

