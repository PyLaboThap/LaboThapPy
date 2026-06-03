# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 12:50:00 2025

@author: Basile
"""

import matplotlib.pyplot as plt
import numpy as np

from labothappy.toolbox.turbomachinery.mean_line_axial_turbine_mapping import map_plot, plot_power_eta_vs_mdot, filter_sparse_by_proximity
from labothappy.component.expander.turbine_mean_line_Aungier import AxialTurbineMeanLine, generate_map_processes

if __name__ == "__main__":

    solve_type    = "mean_line"   # "map" or "mean_line"
    map_mode      = "M_N"   # "M_N", "P_N", or "P_M"
    case_study    = "TCO2_ORC"
    map_generation = True

    # -------------------------------------------------------------------------
    # Operating point — set all three so any mode can pick what it needs
    # -------------------------------------------------------------------------
    m_dot_op = 318.437   # kg/s
    N_rot_op = 2864.775  # RPM
    P_ex_op  = 5220928   # Pa

    #%%
    if case_study == "Salah_Case":
        Turb_OD = AxialTurbineMeanLine('CO2')

        Turb_OD.set_inputs(
            m_dot = 655.18,
            P_su  = 25000000.0,
            T_su  = 923,
            N_rot = 1506.9946780513785,
            fluid = 'CO2',
        )
        
        Turb_OD.set_parameters(
            r_m        = 0.261423771889,
            nStages    = 12,
            mdot_rated = 655.18,
            DP_rated   = 2.5,
            damping    = 0.5,
            delta_tip  = 0.0004,
            N_lw       = 0,
            D_lw       = 0,
            e_blade    = 2e-06,
        )
        
        Turb_OD.set_stage_parameters(
            h_blade_S = [0.05893535333, 0.06254061127, 0.06644686421, 0.07068270152, 0.07527979324, 0.08027326565, 0.08570212832, 0.09160976073, 0.09804446776, 0.105060115, 0.1127168568, 0.1210819728, 0.1254195395],
            chord_S   = [0.008645525688, 0.009066896935, 0.009520522543, 0.01000898242, 0.01053510451, 0.01110199118, 0.0117130487, 0.0123720203, 0.0130830232, 0.01385059029, 0.01467971707, 0.01557591465, 0.01603830059],
            xhi_S1    = [-0.6455466816]*13,
            xhi_S2    = [1.146678265]*13,
            pitch_S   = [0.006966193979, 0.007305716867, 0.007671228935, 0.008064808966, 0.008488735591, 0.008945508564, 0.009437872521, 0.009968843587, 0.01054173924, 0.0111602119, 0.01182828672, 0.01255040431, 0.01292297508],
            o_S       = [0.001744270233, 0.001829283609, 0.001920804434, 0.002019353216, 0.002125500504, 0.002239872211, 0.00236315557, 0.002496105791, 0.002639553539, 0.002794413345, 0.00296169307, 0.003142504605, 0.003235792862],
            t_TE_S    = [0.0005]*13,
            t_blade_S = [0.003025933991, 0.003173413927, 0.00333218289, 0.003503143847, 0.003687286579, 0.003885696912, 0.004099567045, 0.004330207105, 0.00457905812, 0.004847706601, 0.005137900974, 0.005451570129, 0.005613405206],
            n_blade_S = [236, 225, 214, 204, 194, 184, 174, 165, 156, 147, 139, 131, 127],
            R_c_S     = [0.01352981312, 0.01418923794, 0.01489913922, 0.01566355437, 0.01648690899, 0.01737405813, 0.01833033244, 0.01936158987, 0.02047427367, 0.02167547758, 0.02297301932, 0.02437552349, 0.02509913423],
            h_blade_R = [0.06085456338, 0.06461845758, 0.06869830659, 0.07312430505, 0.07792992342, 0.08315231023, 0.08883274986, 0.09501718485, 0.1017568126, 0.1091087686, 0.1171369103, 0.1259127187, None],
            chord_R   = [0.008870626618, 0.009309055087, 0.009781089003, 0.01028943789, 0.01083707309, 0.01142725573, 0.01206356808, 0.01274994866, 0.01349073167, 0.01429069144, 0.01515509243, 0.01608974592, None],
            xhi_R1    = [0.5693595771]*12 + [None],
            xhi_R2    = [-1.17766032]*12  + [None],
            pitch_R   = [0.007763880051, 0.008147607851, 0.008560748305, 0.009005672883, 0.009484982197, 0.01000153051, 0.01055845315, 0.01115919723, 0.01180755622, 0.01250770875, 0.01326426248, 0.01408230362, None],
            o_R       = [0.00229560867, 0.002409068546, 0.002531225097, 0.002662779514, 0.002804500743, 0.002957232725, 0.003121902507, 0.003299529328, 0.003491234828, 0.00369825454, 0.003921950849, 0.004163827632, None],
            t_TE_R    = [0.0005]*12 + [None],
            t_blade_R = [0.003104719316, 0.00325816928, 0.003423381151, 0.003601303262, 0.00379297558, 0.003999539504, 0.004222248827, 0.00446248203, 0.004721756085, 0.005001742002, 0.005304282352, 0.00563141107, None],
            n_blade_R = [212, 202, 192, 182, 173, 164, 156, 147, 139, 131, 124, 117, None],
            R_c_R     = [0.01388208476, 0.01456820328, 0.01530691263, 0.01610245309, 0.01695947464, 0.0178830808, 0.01887887765, 0.01995302875, 0.02111231693, 0.02236421375, 0.02371695786, 0.02517964358, None],
        )

    #%%
    elif case_study == "TCO2_ORC":
        Turb_OD = AxialTurbineMeanLine('CO2')

        Turb_OD.set_inputs(
            m_dot = m_dot_op,
            P_su  = 15309670,
            T_su  = 406.4,
            fluid = 'CO2',
            N_rot = N_rot_op,
            P_ex  = P_ex_op,
        )

        turb_params = {'type': 'Axial Turbine',
            'mdot_rated': 318.437021666738,
            'Wdot_rated': 15244151.612636594,
            'N_rot_rated': 2864.775003723627,
            'total_to_static_efficiency': 0.8851473369105862,
            'DP_rated': 2.93,
            'n_stages': 7,
            'p0_su': 15309670.5,
            'T0_su': 406.4,
            'p_ex': 5220928,
            'r_m': 0.20956432565203412,
            'delta_tip': 0.0004,
            'N_lw': 0,
            'D_lw': 0,
            'e_blade': 2e-06,
            'stator': {
                'h_blade_S': [0.024581131784973807, 0.02731115991390685, 0.030535241526643685, 0.0343428949996286, 0.03884110173212556, 0.04415656370209263, 0.05044235012970171, 0.05396446226956401],
                'chord_S':   [0.017387241421509592, 0.018189646227485024, 0.019157515021682635, 0.020309513815586707, 0.021666371275761282, 0.023250883041521515, 0.025089330758218745, 0.026110623182803244],
                'xhi_S1':    [-0.5493774236565352]*8,
                'xhi_S2':    [1.1268112059950612]*8,
                'pitch_S':   [0.014118294382457852, 0.014769840363216227, 0.015555741716329778, 0.01649115508675976, 0.01759291197811346, 0.01887952225855718, 0.020372326403956538, 0.021201607297417914],
                'o_S':       [0.00407007938170583, 0.004257909709554279, 0.00448447255112079, 0.004754137325719167, 0.005071756287732223, 0.0054426655373106475, 0.005873017193716411, 0.006112085665777279],
                't_TE_S':    [0.0005]*8,
                't_blade_S': [0.0017387241421509593, 0.0018189646227485025, 0.0019157515021682636, 0.002030951381558671, 0.0021666371275761284, 0.0023250883041521517, 0.0025089330758218745, 0.0026110623182803248],
                'n_blade_S': [93, 89, 85, 80, 75, 70, 65, 62],
                'R_c_S':     [0.02593658923532629, 0.02713353838605762, 0.02857731056016104, 0.03029574989120514, 0.03231977738022937, 0.034683397336424615, 0.03742581415244761, 0.03894927848261039],
            },
            'rotor': {
                'h_blade_R': [0.02599623192810525, 0.028975662485097392, 0.03249353115560146, 0.03664835013019067, 0.041556335411602884, 0.04735747516384924, 0.054218932168004906, None],
                'chord_R':   [0.017804442744401783, 0.018691926981905842, 0.019754111574023726, 0.0210108023910845, 0.02248355410347207, 0.02419698360942641, 0.02617885482072475, None],
                'xhi_R1':    [0.5033968740474793]*7 + [None],
                'xhi_R2':    [-1.1544956444218015]*7 + [None],
                'pitch_R':   [0.015148494029818542, 0.015903589253319174, 0.01680732419089942, 0.017876550204477332, 0.01912960658162921, 0.02058743803494024, 0.022273666839919944, None],
                'o_R':       [0.004776675516045738, 0.005014774752793548, 0.005299743578132708, 0.005636895621769747, 0.006032013691274512, 0.006491702145839364, 0.007023409643056108, None],
                't_TE_R':    [0.0005]*7 + [None],
                't_blade_R': [0.0019275069515591894, 0.0020395300073626607, 0.002174650407702253, 0.0022904090712424882, 0.0024358345376652743, 0.002590844225248941, 0.002758858045024202, None],
                'n_blade_R': [87, 83, 78, 74, 69, 64, 59, None],
                'R_c_R':     [0.02655892943743007, 0.02788279177780649, 0.029467255051191132, 0.03134186372129275, 0.03353877094100647, 0.036094697795777626, 0.03905105977860335, None],
            }
        }

        Turb_OD.set_parameters(
            r_m         = turb_params['r_m'],
            nStages     = turb_params['n_stages'],
            mdot_rated  = turb_params['mdot_rated'],
            DP_rated    = turb_params['DP_rated'],
            N_rot_rated = turb_params['N_rot_rated'],
            damping     = 0.3,
            delta_tip   = 0.0004,
            N_lw        = 0,
            D_lw        = 0,
            e_blade     = 2e-06,
            solve_type  = solve_type,
        )
        Turb_OD.set_stage_parameters(
            h_blade_S = turb_params['stator']['h_blade_S'],
            chord_S   = turb_params['stator']['chord_S'],
            xhi_S1    = turb_params['stator']['xhi_S1'],
            xhi_S2    = turb_params['stator']['xhi_S2'],
            pitch_S   = turb_params['stator']['pitch_S'],
            o_S       = turb_params['stator']['o_S'],
            t_TE_S    = turb_params['stator']['t_TE_S'],
            t_blade_S = turb_params['stator']['t_blade_S'],
            n_blade_S = turb_params['stator']['n_blade_S'],
            R_c_S     = turb_params['stator']['R_c_S'],
            h_blade_R = turb_params['rotor']['h_blade_R'],
            chord_R   = turb_params['rotor']['chord_R'],
            xhi_R1    = turb_params['rotor']['xhi_R1'],
            xhi_R2    = turb_params['rotor']['xhi_R2'],
            pitch_R   = turb_params['rotor']['pitch_R'],
            o_R       = turb_params['rotor']['o_R'],
            t_TE_R    = turb_params['rotor']['t_TE_R'],
            t_blade_R = turb_params['rotor']['t_blade_R'],
            n_blade_R = turb_params['rotor']['n_blade_R'],
            R_c_R     = turb_params['rotor']['R_c_R'],
        )

    #%%
    # MAP_SAVE_PATH = r"Your_Path.parquet"   # <-- set your path here
    MAP_SAVE_PATH = r"C:\Users\Basile\Desktop\Travail\Thèse\Travail\WP1\Turbomachines\Save Maps\turb_map.parquet"
    
    # -------------------------------------------------------------------------
    # Map generation (optional)
    # -------------------------------------------------------------------------
    if map_generation:
        df_map = generate_map_processes(
            Turb_OD,
            m_grid      = np.linspace(0.6*Turb_OD.params['mdot_rated'], 1.4*Turb_OD.params['mdot_rated'], 30),
            N_grid      = np.linspace(0.3*Turb_OD.params['N_rot_rated'], 1.5*Turb_OD.params['N_rot_rated'], 30),
            max_workers = -2,
        )
    
        Turb_OD.save_map_df(MAP_SAVE_PATH, df_map)
        df_reloaded = Turb_OD.load_map_df(MAP_SAVE_PATH)
    
        assert len(df_reloaded) == len(df_map), "Row count mismatch after reload!"
        print(f"[check] DataFrame round-trip OK  ({len(df_reloaded)} rows)")
    
        df_ok = df_reloaded[df_reloaded["converged"] == True].dropna(
            subset=["W_dot", "eta_is", "P_ex_calc"]
        )
    
        N_GRID  = 120
        m_fine  = np.linspace(df_ok["m_dot"].min(), df_ok["m_dot"].max(), N_GRID)
        N_fine  = np.linspace(df_ok["N_rot"].min(), df_ok["N_rot"].max(), N_GRID)
        MM, NN  = np.meshgrid(m_fine, N_fine)
        m_flat  = MM.ravel()
        N_flat  = NN.ravel()
    
        interp   = Turb_OD.map_interpolator
        res      = interp.query_batch(m_flat, N_flat)
    
        W_grid   = res["W_dot"].reshape(N_GRID, N_GRID)    / 1e6
        eta_grid = res["eta_is"].reshape(N_GRID, N_GRID)
        Pex_grid = res["P_ex_calc"].reshape(N_GRID, N_GRID) / 1e5
    
        hull_mask = interp.hull_mask(m_flat, N_flat).reshape(N_GRID, N_GRID)
        W_grid    = np.where(hull_mask, W_grid,   np.nan)
        eta_grid  = np.where(hull_mask, eta_grid, np.nan)
        Pex_grid  = np.where(hull_mask, Pex_grid, np.nan)
    
        TARGETS  = [
            ("W_dot [MW]", df_ok["W_dot"]    / 1e6, W_grid,   "viridis"),
            ("eta_is [–]", df_ok["eta_is"],          eta_grid, "plasma"),
            ("P_ex [bar]", df_ok["P_ex_calc"] / 1e5, Pex_grid, "cividis"),
        ]
        N_LEVELS = 20
    
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(13, 13), constrained_layout=True)
        fig.suptitle("Performance map  –  raw data (left)  vs  interpolated reconstruction (right)",
                     fontsize=13, fontweight="bold")
    
        for row, (label, raw_vals, grid_vals, cmap) in enumerate(TARGETS):
            vmin   = float(np.nanmin(grid_vals))
            vmax   = float(np.nanmax(grid_vals))
            levels = np.linspace(vmin, vmax, N_LEVELS + 1)
    
            ax_raw = axes[row, 0]
            sc = ax_raw.scatter(df_ok["m_dot"], df_ok["N_rot"],
                                c=raw_vals, cmap=cmap, vmin=vmin, vmax=vmax,
                                s=40, edgecolors="k", linewidths=0.4, zorder=3)
            ax_raw.set_title(f"{label}  –  raw map points", fontsize=10)
            ax_raw.set_xlabel("ṁ  [kg/s]")
            ax_raw.set_ylabel("N  [rpm]")
            fig.colorbar(sc, ax=ax_raw, label=label, pad=0.02)
    
            ax_int = axes[row, 1]
            cf = ax_int.contourf(MM, NN, grid_vals, levels=levels, cmap=cmap, extend="neither")
            ax_int.contour(MM, NN, grid_vals, levels=levels, colors="white", linewidths=0.35, alpha=0.5)
            ax_int.scatter(df_ok["m_dot"], df_ok["N_rot"],
                           c=raw_vals, cmap=cmap, vmin=vmin, vmax=vmax,
                           s=25, edgecolors="k", linewidths=0.5, zorder=4)
            ax_int.set_title(f"{label}  –  interpolated (method={interp.method!r}, n={interp.n_points} pts)",
                             fontsize=10)
            ax_int.set_xlabel("ṁ  [kg/s]")
            ax_int.set_ylabel("N  [rpm]")
            fig.colorbar(cf, ax=ax_int, label=label, pad=0.02)
    
        plt.show()
    
        res_nodes = interp.query_batch(df_ok["m_dot"].to_numpy(), df_ok["N_rot"].to_numpy())
        err_W     = (res_nodes["W_dot"]     - df_ok["W_dot"].to_numpy())    / 1e6
        err_eta   = (res_nodes["eta_is"]    - df_ok["eta_is"].to_numpy())
        err_Pex   = (res_nodes["P_ex_calc"] - df_ok["P_ex_calc"].to_numpy()) / 1e5
    
        print("\n── Interpolation residuals at map nodes ──────────────────────────────")
        print(f"  W_dot    :  max |err| = {np.abs(err_W).max():.4f} MW     RMS = {np.sqrt((err_W**2).mean()):.4f} MW")
        print(f"  eta_is   :  max |err| = {np.abs(err_eta).max():.6f}      RMS = {np.sqrt((err_eta**2).mean()):.6f}")
        print(f"  P_ex_calc:  max |err| = {np.abs(err_Pex).max():.4f} bar  RMS = {np.sqrt((err_Pex**2).mean()):.4f} bar")
        print("──────────────────────────────────────────────────────────────────────")
    
        df_clean = filter_sparse_by_proximity(df_map, rp_col=None, group_by='N_rot',
                                              rp_tol_rel=0.6, m_tol_rel=0.6, min_neighbors=2)
    
        fig, ax = map_plot(
            df_clean, rp_col='RP_calc',
            use_grid=True, nx=600, ny=600,
            min_circle_ratio=0.01, max_area_factor=50.0, long_edge_q=1,
            fill_holes=True, hole_method='nearest', hole_smooth_sigma=0.8,
            smooth_sigma=0.6, show_points=True,
            levels=24, focus_high=True, max_iso_speeds=4,
            figsize=(9, 6), dpi=220,
        )
        plt.show()
    
        _ = plot_power_eta_vs_mdot(df_map, speeds=None, max_lines=5)
        plt.show()     
        
        # -------------------------------------------------------------------------
        # Diagnose P_N inversion quality — does P_ex_calc actually vary with m_dot?
        # -------------------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle("P_N inversion diagnosis", fontsize=13, fontweight="bold")
    
        N_lines   = np.linspace(interp.N_range[0], interp.N_range[1], 7)
        cmap      = plt.cm.viridis
        colors    = [cmap(i / (len(N_lines) - 1)) for i in range(len(N_lines))]
        N_rated  = turb_params['N_rot_rated']   # needed for the diagnosis labels
    
        for color, N_target in zip(colors, N_lines):
            # Grab raw map points near this N_rot (±5% band)
            mask = np.abs(df_ok["N_rot"] - N_target) / N_target < 0.05
            sub  = df_ok[mask].sort_values("m_dot")
            if len(sub) < 2:
                continue
            label = f"N = {N_target/N_rated*100:.0f}% N$_{{r}}$"
    
            # Left: raw (m_dot, P_ex_calc) — shows if P_ex varies with m_dot
            axes[0].plot(sub["m_dot"], sub["P_ex_calc"] / 1e5,
                         'o-', color=color, lw=1.5, ms=4, label=label)
    
            # Right: inversion result — query_PN over the P_ex range of this N slice
            P_sweep = np.linspace(sub["P_ex_calc"].min(), sub["P_ex_calc"].max(), 40)
            m_inv   = [interp.query_PN(N_target, P)["m_dot"] for P in P_sweep]
            axes[1].plot(P_sweep / 1e5, m_inv,
                         color=color, lw=1.8, label=label)
    
        axes[0].set_xlabel("$\\dot{m}$  [kg/s]",      fontsize=11)
        axes[0].set_ylabel("$P_{ex,calc}$  [bar]",     fontsize=11)
        axes[0].set_title("Raw map: $P_{ex}$ vs $\\dot{m}$ at fixed N\n"
                          "(slope = sensitivity; flat = degenerate inverse)", fontsize=10)
        axes[0].grid(True, alpha=0.3)
    
        axes[1].set_xlabel("$P_{ex}$  [bar]",          fontsize=11)
        axes[1].set_ylabel("$\\dot{m}$ recovered  [kg/s]", fontsize=11)
        axes[1].set_title("Interpolated inverse: $\\dot{m}$(N, $P_{ex}$)\n"
                          "(should match left panel slope)", fontsize=10)
        axes[1].grid(True, alpha=0.3)
    
        handles, labels_leg = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels_leg, loc='lower center', ncol=4,
                   bbox_to_anchor=(0.5, -0.13), fontsize=9, framealpha=0.8)
        plt.tight_layout()
        plt.show()
    
        # Print the actual P_ex range per N slice to quantify degeneracy
        print("\n── P_ex spread per iso-speed slice ───────────────────────────────────")
        for N_target in N_lines:
            mask  = np.abs(df_ok["N_rot"] - N_target) / N_target < 0.05
            sub   = df_ok[mask]
            if len(sub) < 2:
                continue
            dP    = sub["P_ex_calc"].max() - sub["P_ex_calc"].min()
            dm    = sub["m_dot"].max()     - sub["m_dot"].min()
            print(f"  N={N_target:7.1f} rpm  |  ΔP_ex={dP/1e5:6.2f} bar  |  "
                  f"Δm_dot={dm:6.2f} kg/s  |  "
                  f"sensitivity={dP/dm/1e5:.3f} bar/(kg/s)")
        print("──────────────────────────────────────────────────────────────────────")
            
    #%%
    # -------------------------------------------------------------------------
    # Solve
    # -------------------------------------------------------------------------

    if solve_type == "map":
        Turb_OD.load_map_df(MAP_SAVE_PATH)

        if map_mode == "M_N":
            
            Turb_OD.set_inputs(m_dot=m_dot_op, N_rot=N_rot_op)
            Turb_OD.set_parameters(map_mode='M_N')
            Turb_OD.solve()

            print(f"[M_N] m_dot={m_dot_op:.3f} kg/s | N_rot={N_rot_op:.2f} RPM")
            print(f"  → P_ex  = {Turb_OD.P_ex_calc/1e5:.3f} bar")

        elif map_mode == "P_N":
            Turb_OD.set_inputs(P_ex=P_ex_op, N_rot=N_rot_op)
            Turb_OD.set_parameters(map_mode='P_N')
            Turb_OD.solve()

            print(f"[P_N] N_rot={N_rot_op:.2f} RPM | P_ex={P_ex_op/1e5:.3f} bar")
            print(f"  → m_dot = {Turb_OD.inputs['m_dot']:.3f} kg/s")

        elif map_mode == "P_M":
            Turb_OD.set_inputs(m_dot=m_dot_op, P_ex=P_ex_op)
            Turb_OD.set_parameters(map_mode='P_M')
            Turb_OD.solve()

            print(f"[P_M] m_dot={m_dot_op:.3f} kg/s | P_ex={P_ex_op/1e5:.3f} bar")
            print(f"  → N_rot = {Turb_OD.inputs['N_rot']:.2f} RPM")

        print(f"  W_dot   = {Turb_OD.W_dot/1e6:.3f} MW")
        print(f"  eta_is  = {Turb_OD.eta_is:.4f}")

    else:
        Turb_OD.solve()
        print(f"[mean_line] W_dot={Turb_OD.W_dot/1e6:.3f} MW | eta_is={Turb_OD.eta_is:.4f}")

   
