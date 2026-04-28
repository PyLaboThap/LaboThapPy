

from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector

from labothappy.correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model
from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.toolbox.turbomachinery.mean_line_axial_turbine_mapping import map_plot, map_plot_clean, plot_power_eta_vs_mdot, filter_sparse_by_proximity

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, RBFInterpolator
from scipy.spatial import Delaunay
import pyswarms as ps

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import warnings

warnings.filterwarnings("ignore")


class MapInterpolationError(RuntimeError):
    """Raised when map-based interpolation cannot produce a valid result."""


class AxialTurbineMeanLine(BaseComponent):

    def __init__(self, fluid):
        super().__init__()
        self.inputs  = {}
        self.params  = {}
        self.fluid   = fluid
        self.AS      = CP.AbstractState('HEOS', fluid)
        self.stages  = []
        self.Vel_Tri_Last_Stage = {}
        
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W  = WorkConnector()
        
        self.Dh0_stage_guess    = 0
        self.map_interpolator   = None

    # =========================================================================
    #  Inner class: stage  —  all state stored in plain dicts
    # =========================================================================

    class stage:

        # State keys stored per position index (1, 2, 3)
        _STATE_KEYS = ('H', 'S', 'P', 'D', 'A', 'V')

        def __init__(self, fluid):
            self.AS = CP.AbstractState('HEOS', fluid)

            # Replace DataFrames with nested dicts: total_states[pos][key]
            self.total_states  = {
                1: dict.fromkeys(self._STATE_KEYS, np.nan),
                2: dict.fromkeys(self._STATE_KEYS, np.nan),
                3: dict.fromkeys(self._STATE_KEYS, np.nan),
            }
            self.static_states = {
                1: dict.fromkeys(self._STATE_KEYS, np.nan),
                2: dict.fromkeys(self._STATE_KEYS, np.nan),
                3: dict.fromkeys(self._STATE_KEYS, np.nan),
            }

            self.eta_is_R = None
            self.eta_is_S = None
            self.A_flow_S = None
            self.A_flow_R = None
            self.h_blade_S = None
            self.h_blade_R = None
            self.chord_S   = None
            self.chord_R   = None
            self.stage     = None
            self.AR        = None
            self.xhi_S1    = None
            self.xhi_S2    = None
            self.xhi_R1    = None
            self.xhi_R2    = None
            self.Vel_Tri_R = {}
            self.Vel_Tri_S = {}

        # ------------------------------------------------------------------
        #  State updaters — write directly into dicts, no .loc[] overhead
        # ------------------------------------------------------------------

        def _write_state(self, target_dict, position):
            """Write current AS state into target_dict[position]."""
            d = target_dict[position]
            d['H'] = self.AS.hmass()
            d['S'] = self.AS.smass()
            d['P'] = self.AS.p()
            d['D'] = self.AS.rhomass()
            try:
                d['A'] = self.AS.speed_sound()
            except Exception:
                d['A'] = -1.0
            d['V'] = self.AS.viscosity()

        def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            self._write_state(self.total_states, position)

        def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            self._write_state(self.static_states, position)

    # =========================================================================
    #  Inner class: MapInterpolator  (unchanged)
    # =========================================================================

    class MapInterpolator:
        _DEFAULT_TARGETS = ["W_dot", "eta_is", "P_ex_calc"]

        def __init__(
            self,
            df,
            *,
            method        = "linear",
            m_dot_col     = "m_dot",
            N_rot_col     = "N_rot",
            targets       = None,
            normalize     = True,
            fallback_nn   = True,
            min_converged_frac = 0.2,
        ):
            if method not in ("linear", "rbf", "nearest"):
                raise ValueError(f"method must be 'linear', 'rbf', or 'nearest'; got {method!r}")

            self.method       = method
            self.m_dot_col    = m_dot_col
            self.N_rot_col    = N_rot_col
            self.targets      = targets if targets is not None else list(self._DEFAULT_TARGETS)
            self.normalize    = normalize
            self.fallback_nn  = fallback_nn

            conv_col = df.get("converged", pd.Series(True, index=df.index))
            df_ok    = df[conv_col == True].copy()
            frac     = len(df_ok) / max(len(df), 1)

            if frac < min_converged_frac:
                raise MapInterpolationError(
                    f"Only {frac:.1%} of map points converged "
                    f"(threshold {min_converged_frac:.1%})."
                )

            df_ok = df_ok.dropna(subset=[m_dot_col, N_rot_col] + self.targets)
            if len(df_ok) < 4:
                raise MapInterpolationError(
                    f"Only {len(df_ok)} usable points after filtering – need at least 4."
                )

            m_raw = df_ok[m_dot_col].to_numpy(dtype=float)
            N_raw = df_ok[N_rot_col].to_numpy(dtype=float)

            self._m_scale = (float(m_raw.min()), float(m_raw.max()))
            self._N_scale = (float(N_raw.min()), float(N_raw.max()))

            pts = self._normalize(m_raw, N_raw)

            self._primary = {}
            self._nn      = {}

            for col in self.targets:
                if col not in df_ok.columns:
                    warnings.warn(f"MapInterpolator: target column {col!r} not found – skipping.")
                    continue
                vals = df_ok[col].to_numpy(dtype=float)
                if method == "linear":
                    self._primary[col] = LinearNDInterpolator(pts, vals)
                elif method == "rbf":
                    self._primary[col] = RBFInterpolator(pts, vals, kernel="thin_plate_spline", smoothing=0.0)
                else:
                    self._primary[col] = NearestNDInterpolator(pts, vals)
                self._nn[col] = NearestNDInterpolator(pts, vals)

            try:
                self._hull = Delaunay(pts)
            except Exception:
                self._hull = None

            self._pts_raw = np.column_stack([m_raw, N_raw])
            self.n_points = len(df_ok)
            self.m_range  = self._m_scale
            self.N_range  = self._N_scale

        def _normalize(self, m, N):
            m, N = np.asarray(m, float), np.asarray(N, float)
            if not self.normalize:
                return np.column_stack([m, N])
            m_lo, m_hi = self._m_scale
            N_lo, N_hi = self._N_scale
            m_n = (m - m_lo) / max(m_hi - m_lo, 1e-12)
            N_n = (N - N_lo) / max(N_hi - N_lo, 1e-12)
            return np.column_stack([m_n, N_n])

        def hull_mask(self, m_arr, N_arr):
            pts = self._normalize(np.asarray(m_arr, float), np.asarray(N_arr, float))
            if self._hull is None:
                return np.ones(len(pts), dtype=bool)
            return self._hull.find_simplex(pts) >= 0

        def query(self, m_dot, N_rot):
            res = self.query_batch(np.array([m_dot], dtype=float), np.array([N_rot], dtype=float))
            return {k: float(v[0]) for k, v in res.items()}

        def query_batch(self, m_arr, N_arr):
            m_arr = np.asarray(m_arr, dtype=float)
            N_arr = np.asarray(N_arr, dtype=float)
            pts   = self._normalize(m_arr, N_arr)
            out   = {}
            for col, interp in self._primary.items():
                vals     = np.asarray(interp(pts), dtype=float).ravel()
                nan_mask = ~np.isfinite(vals)
                if nan_mask.any():
                    if self.fallback_nn:
                        vals[nan_mask] = self._nn[col](pts[nan_mask]).ravel()
                    else:
                        raise MapInterpolationError(
                            f"{nan_mask.sum()} query point(s) outside hull for {col!r}."
                        )
                out[col] = vals
            return out

        def save(self, path):
            with open(path, "wb") as fh:
                pickle.dump(self, fh, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"[MapInterpolator] Saved to {path}  ({self.n_points} pts)")

        @staticmethod
        def load(path):
            with open(path, "rb") as fh:
                obj = pickle.load(fh)
            if not hasattr(obj, "_primary"):
                raise TypeError(f"Expected MapInterpolator, got {type(obj)}")
            return obj

        def plot_coverage(self, ax=None, color="steelblue", alpha=0.55):
            if ax is None:
                _, ax = plt.subplots(figsize=(6, 4))
            ax.scatter(self._pts_raw[:, 0], self._pts_raw[:, 1],
                       s=14, c=color, alpha=alpha, label=f"{self.n_points} pts")
            ax.set_xlabel(f"{self.m_dot_col} [kg/s]")
            ax.set_ylabel(f"{self.N_rot_col} [rpm]")
            ax.set_title(
                f"MapInterpolator  –  method={self.method!r}  |  "
                f"m ∈ [{self.m_range[0]:.1f}, {self.m_range[1]:.1f}]  "
                f"N ∈ [{self.N_range[0]:.0f}, {self.N_range[1]:.0f}]"
            )
            ax.legend(fontsize=8)
            return ax

        def __repr__(self):
            return (
                f"MapInterpolator(method={self.method!r}, "
                f"n_points={self.n_points}, "
                f"m_range={self.m_range}, "
                f"N_range={self.N_range}, "
                f"targets={self.targets})"
            )

    # =========================================================================
    #  Data Handling
    # =========================================================================

    def set_stage_parameters(self, **parameters):
        n_stages = max(len(arr) for arr in parameters.values())
        if not hasattr(self, "stages"):
            self.stages = []
        while len(self.stages) < n_stages:
            self.stages.append(self.stage(self.fluid))
        for key, array in parameters.items():
            for i in range(len(array)):
                setattr(self.stages[i], key, array[i])

    def _snapshot_from_machine(self):
        base_inputs = dict(
            P_su  = float(self.inputs['P_su']),
            T_su  = float(self.inputs['T_su']),
            P_ex  = float(self.inputs['P_ex']),
            m_dot = float(self.inputs['m_dot']),
            N_rot = float(self.inputs['N_rot']),
            fluid = self.fluid,
        )
        base_params = dict(self.params)
        stage_params = {}
        import numbers
        def _ok(v): return (v is None) or isinstance(v, (numbers.Number, np.floating, np.integer))
        if self.stages:
            keys      = set().union(*(vars(st).keys() for st in self.stages))
            blacklist = {
                'AS', 'total_states', 'static_states',
                'Vel_Tri_S', 'Vel_Tri_R',
                'eta_is_R', 'eta_is_S',
                'M1_S', 'M2_S', 'M2_R', 'M3_R',
                'Y_vec_S', 'Y_vec_R',
                'delta_S', 'delta_R',
                'beta_g_S', 'beta_g_R',
            }
            for k in sorted(keys - blacklist):
                vals = [getattr(st, k, None) for st in self.stages]
                if any(_ok(v) for v in vals) and all(_ok(v) for v in vals):
                    stage_params[k] = vals
        return base_inputs, base_params, stage_params

    #%% Stage solving
    
    def _wegstein_solve(self, f, x0, tol=1e-6, max_iter=100, q_min=-5.0, q_max=0.0):
        """
        Wegstein Fixed-Point Solving approach used to solve the stages systems
        
        First iterations use regular FPI with damping factors and then switches to Wegstein when necessary history has been computed.
        """
        
        x_in       = np.array(x0, dtype=float)
        x_out_prev = None
        x_in_prev  = None

        for c in range(max_iter):
            x_out = np.array(f(x_in), dtype=float)
            res   = np.sum(np.abs((x_in - x_out) / np.where(np.abs(x_out) > 1e-30, x_out, 1.0)))

            if res < tol:
                f(x_out)
                return x_out

            if x_in_prev is None:
                x_new = (1.0 - self.params['damping']) * x_in + self.params['damping'] * x_out
            else:
                dx   = x_in  - x_in_prev
                df_  = x_out - x_out_prev
                mask = np.abs(dx) > 1e-30
                q    = np.where(mask, np.divide(df_, dx, where=mask, out=np.zeros_like(dx)), 0.0)
                q    = np.clip(q, q_min, q_max)
                w    = np.where(mask, q / (q - 1.0), self.params['damping'])
                x_new = w * x_in + (1.0 - w) * x_out

            x_in_prev  = x_in.copy()
            x_out_prev = x_out.copy()
            x_in       = x_new

        raise RuntimeError(
            f"Wegstein failed to converge after {max_iter} iterations "
            f"(final residual: {res:.3e})"
        )

    #%%  Map interface  (unchanged public API)

    def load_map(self, df, *, method="linear", targets=None, normalize=True,
                 fallback_nn=True, min_converged_frac=0.2, verbose=True):
        self.map_interpolator = self.MapInterpolator(
            df, method=method, targets=targets, normalize=normalize,
            fallback_nn=fallback_nn, min_converged_frac=min_converged_frac,
        )
        if verbose:
            print(
                f"[load_map] {self.map_interpolator.n_points} pts loaded  "
                f"| method={method!r}  "
                f"| m ∈ [{self.map_interpolator.m_range[0]:.1f}, "
                f"{self.map_interpolator.m_range[1]:.1f}] kg/s  "
                f"| N ∈ [{self.map_interpolator.N_range[0]:.0f}, "
                f"{self.map_interpolator.N_range[1]:.0f}] rpm"
            )
        return self.map_interpolator

    def load_map_from_file(self, path, verbose=True):
        self.map_interpolator = self.MapInterpolator.load(path)
        if verbose:
            print(f"[load_map_from_file] Loaded {self.map_interpolator}")
        return self.map_interpolator

    def save_map_df(self, path, df, fmt="parquet"):
        from pathlib import Path as _Path
        p = _Path(path)
        if not p.suffix:
            p = p.with_suffix("." + fmt)
        fmt = fmt.lower()
        if fmt == "parquet":
            df.to_parquet(p, index=False)
        elif fmt == "csv":
            df.to_csv(p, index=False)
        else:
            raise ValueError(f"fmt must be 'parquet' or 'csv'; got {fmt!r}")
        print(f"[save_map_df] {len(df)} rows saved to '{p}'  (fmt={fmt!r})")

    def load_map_df(self, path, fmt="auto"):
        from pathlib import Path as _Path
        p = _Path(path)
        if fmt == "auto":
            fmt = p.suffix.lstrip(".").lower()
            if fmt not in ("parquet", "csv"):
                raise ValueError(f"Cannot infer format from extension {p.suffix!r}.")
        if fmt == "parquet":
            df = pd.read_parquet(p)
        elif fmt == "csv":
            df = pd.read_csv(p)
            if "converged" in df.columns:
                df["converged"] = df["converged"].astype(bool)
        else:
            raise ValueError(f"fmt must be 'parquet' or 'csv'; got {fmt!r}")
        print(f"[load_map_df] {len(df)} rows loaded from '{p}'")
        self.load_map(df, verbose=True)
        return df

    def solve_from_map(self, m_dot=None, N_rot=None, *, strict_hull=False):
        if self.map_interpolator is None:
            raise MapInterpolationError("No performance map loaded. Call turbine.load_map(df_map) first.")

        _m = float(m_dot if m_dot is not None else self.inputs["m_dot"])
        _N = float(N_rot if N_rot is not None else self.inputs["N_rot"])

        if strict_hull:
            inside = self.map_interpolator.hull_mask(np.array([_m]), np.array([_N]))[0]
            if not inside:
                raise MapInterpolationError(
                    f"Operating point (m_dot={_m:.3f} kg/s, N_rot={_N:.1f} rpm) "
                    "is outside the map convex hull."
                )

        result = self.map_interpolator.query(_m, _N)

        self.W_dot     = result.get("W_dot",     float("nan"))
        self.eta_is    = result.get("eta_is",    float("nan"))
        self.P_ex_calc = result.get("P_ex_calc", float("nan"))

        if self.ex is not None:
            try:
                self.ex.p = self.P_ex_calc
            except AttributeError:
                pass

        self._from_map = True

        if m_dot is not None: self.inputs["m_dot"] = _m
        if N_rot is not None: self.inputs["N_rot"] = _N

        self.AS.update(CP.PSmass_INPUTS, self.P_ex_calc, self.su.s)
        h_ex_is  = self.AS.hmass()
        self.h_ex = self.su.h - self.eta_is * (self.su.h - h_ex_is)
        self.solved = True
        self.update_connectors(self.P_ex_calc, self.h_ex)

    # =========================================================================
    #  Loss Models  — only dict access changed, all physics identical
    # =========================================================================

    def stator_blade_row_system(self, x):
        stage = self.stages[self.curr_stage_index]
        ss    = stage.static_states
        ts    = stage.total_states

        h_static_out = x[0] * 1e5
        p_static_out = x[1] * 1e5

        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)

        stage.Vel_Tri_S['u']  = u  = self.u
        A_flow                     = stage.h_blade_S * (2 * np.pi * self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot'] / (ss[2]['D'] * A_flow)

        if self.curr_stage_index == 0:
            stage.Vel_Tri_S['alpha1'] = alpha1 = stage.xhi_S1
            stage.Vel_Tri_S['vu1']    = vu1    = vm * np.tan(alpha1)
        else:
            wu1 = np.tan(stage.Vel_Tri_S['beta1']) * vm
            stage.Vel_Tri_S['wu1']    = wu1
            stage.Vel_Tri_S['vu1']    = vu1    = wu1 + u
            stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1 / vm)

        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(vm**2 + vu1**2)
        stage.M1_S = v1 / ss[1]['A']

        hin  = ss[1]['H']
        h0in = hin + (vu1**2 + vm**2) / 2
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, ss[1]['S'], 1)

        h02 = h0in
        arg_v2 = max(0.0, 2 * (h02 - h_static_out))
        stage.Vel_Tri_R['v2'] = v2 = np.sqrt(arg_v2)
        stage.M2_S = v2 / ss[2]['A']
        stage.Vel_Tri_S['alpha2'] = np.arctan2(np.sqrt(v2**2 - vm**2), vm)

        stage.Re_s  = stage.chord_S * (ss[2]['D'] * vm) / ss[2]['V']
        stage.AR_S  = stage.h_blade_S / stage.chord_S
        stage.beta_g_S = np.arcsin(stage.o_S / stage.pitch_S)

        stage.Y_vec_S = aungier_loss_model(
            stage.Vel_Tri_S['alpha1'], stage.Vel_Tri_S['alpha2'],
            stage.beta_g_S * 180 / np.pi, stage.xhi_S1, stage.chord_S,
            0, self.params['D_lw'], self.params['e_blade'],
            stage.h_blade_S, ss[2]['V'],
            stage.M1_S, stage.M2_S, self.params['N_lw'],
            stage.R_c_S, ss[2]['D'], stage.pitch_S,
            stage.t_blade_S, stage.t_TE_S, vm, v2, 1,
        )

        self.compute_deviation_stator(stage)
        alpha2_calc = stage.xhi_S2 + stage.delta_S
        v2_new      = vm / np.cos(alpha2_calc)

        Y     = stage.Y_vec_S['Y_tot']
        p0_out = (ts[1]['P'] + Y * p_static_out) / (1 + Y)
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)

        sout = ts[2]['S']
        hout = ts[2]['H'] - v2_new**2 / 2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)

        pout_calc = ss[2]['P']
        self.AS.update(CP.PSmass_INPUTS, pout_calc, ss[1]['S'])
        hout_s = self.AS.hmass()
        stage.eta_is_S = (ss[1]['H'] - ss[2]['H']) / (ss[1]['H'] - hout_s)

        return np.array([hout, pout_calc]) * 1e-5

    def rotor_blade_row_system(self, x):
        stage = self.stages[self.curr_stage_index]
        ss    = stage.static_states
        ts    = stage.total_states

        h_static_out = x[0] * 1e5
        p_static_out = x[1] * 1e5

        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)

        stage.Vel_Tri_R['u']   = u  = self.u
        A_flow                      = stage.h_blade_R * (2 * np.pi * self.params['r_m'])
        stage.Vel_Tri_R['vm']  = vm = self.inputs['m_dot'] / (ss[3]['D'] * A_flow)
        stage.Vel_Tri_R['vu2'] = vu2 = vm * np.tan(stage.Vel_Tri_R['alpha2'])

        wu2 = vu2 - u
        stage.Vel_Tri_R['wu2'] = wu2
        stage.Vel_Tri_R['w2']  = w2 = np.sqrt(wu2**2 + vm**2)
        stage.M2_R              = w2 / ss[2]['A']
        stage.Vel_Tri_R['beta2'] = np.arctan(wu2 / vm)

        hin  = ss[2]['H']
        h0in = hin + w2**2 / 2
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, ss[2]['S'], 2)

        h03 = ts[2]['H']
        arg_w3 = max(0.0, 2 * (h03 - h_static_out))
        stage.Vel_Tri_R['w3'] = w3 = np.sqrt(arg_w3)
        stage.M3_R = w3 / ss[3]['A']
        stage.Vel_Tri_R['beta3'] = -np.arccos(vm / w3)

        stage.Re_r   = stage.chord_R * (ss[3]['D'] * vm) / ss[3]['V']
        stage.AR_R   = stage.h_blade_R / stage.chord_R
        stage.beta_g_R = np.arcsin(stage.o_R / stage.pitch_R)

        stage.Y_vec_R = aungier_loss_model(
            -stage.Vel_Tri_R['beta2'], -stage.Vel_Tri_R['beta3'],
            stage.beta_g_R * 180 / np.pi, -stage.xhi_R1, stage.chord_R,
            self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'],
            stage.h_blade_R, ss[3]['V'],
            stage.M2_R, stage.M3_R, self.params['N_lw'],
            stage.R_c_R, ss[3]['D'], stage.pitch_R,
            stage.t_blade_R, stage.t_TE_R, vm, w3, 1,
        )

        self.compute_deviation_rotor(stage)
        beta3_calc = stage.xhi_R2 + stage.delta_R
        w3_new     = vm / np.cos(beta3_calc)

        Y      = stage.Y_vec_R['Y_tot']
        p0_out = (ts[2]['P'] + Y * p_static_out) / (1 + Y)
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)

        sout = ts[3]['S']
        hout = h0in - w3_new**2 / 2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 3)

        pout_calc = ss[3]['P']
        self.AS.update(CP.PSmass_INPUTS, pout_calc, ss[2]['S'])
        hout_s    = self.AS.hmass()
        stage.eta_is_R = (ss[2]['H'] - ss[3]['H']) / (ss[2]['H'] - hout_s)

        return np.array([hout, pout_calc]) * 1e-5

    def last_blade_row_system(self, x):
        stage = self.stages[-1]
        ss    = stage.static_states
        ts    = stage.total_states

        h_static_out = x[0] * 1e5
        p_static_out = x[1] * 1e5

        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)

        stage.Vel_Tri_S['u']  = u  = self.u
        A_flow                     = stage.h_blade_S * (2 * np.pi * self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot'] / (ss[2]['D'] * A_flow)

        wu1 = np.tan(stage.Vel_Tri_S['beta1']) * vm
        stage.Vel_Tri_S['wu1']    = wu1
        stage.Vel_Tri_S['vu1']    = vu1    = wu1 + u
        stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1 / vm)

        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(vm**2 + vu1**2)
        stage.M1_S = v1 / ss[1]['A']

        hin  = ss[1]['H']
        h0in = hin + (vu1**2 + vm**2) / 2
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, ss[1]['S'], 1)

        h02 = h0in
        arg_v2 = max(0.0, 2 * (h02 - h_static_out))
        stage.Vel_Tri_R['v2'] = v2 = np.sqrt(arg_v2)
        stage.M2_S = v2 / ss[2]['A']
        stage.Vel_Tri_S['alpha2'] = np.arctan2(np.sqrt(v2**2 - vm**2), vm)

        alpha2 = stage.Vel_Tri_S['alpha2']
        solidity = stage.chord_S / stage.pitch_S
        a        = 0.0117

        D_e = (np.cos(alpha2) / np.cos(alpha1)) * (
            1.12
            + a * (alpha1 - stage.xhi_S1)
            + 0.61 * np.cos(alpha1)**2 / solidity * (np.tan(alpha1) - np.tan(alpha2))
        )

        P_cst  = np.cos(alpha2) / 2 * solidity * (v1 / v2)**2
        Yp     = 0.004 * (1 + 3.1 * (D_e - 1)**2 + 0.4 * (D_e - 1)**8) / P_cst
        EW_Cst = np.cos((alpha1 + alpha2) / 2)**3 / np.cos(alpha1)**2
        AR_S   = stage.h_blade_S / stage.chord_S
        Yew    = 0.02 * (solidity / AR_S) / EW_Cst

        DP_loss = (Yp + Yew) * (v1**2) * ss[1]['D'] / 2
        p0_out  = ts[1]['P'] - DP_loss
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)

        sout = ts[2]['S']
        hout = h0in - v2**2 / 2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)

        pout_calc = ss[2]['P']
        self.AS.update(CP.PSmass_INPUTS, pout_calc, ss[1]['S'])
        hout_s    = self.AS.hmass()
        stage.eta_is_S = (ss[1]['H'] - ss[2]['H']) / (ss[1]['H'] - hout_s)

        return np.array([hout, pout_calc]) * 1e-5

    # =========================================================================
    #  Deviation models  (unchanged)
    # =========================================================================

    def compute_deviation_stator(self, stage):
        delta_0S = (
            np.arcsin(
                (stage.o_S / stage.pitch_S)
                * (1 + (1 - stage.o_S / stage.pitch_S) * (2 * stage.beta_g_S / np.pi)**2)
            )
            - stage.beta_g_S
        )
        if stage.M2_S <= 0.5:
            stage.delta_S = delta_0S
        else:
            X = 2 * stage.M2_S - 1
            stage.delta_S = delta_0S * (1 - 10*X**3 + 15*X**4 - 6*X**5)

    def compute_deviation_rotor(self, stage):
        delta_0R = (
            np.arcsin(
                (stage.o_R / stage.pitch_R)
                * (1 + (1 - stage.o_R / stage.pitch_R) * (2 * stage.beta_g_R / np.pi)**2)
            )
            - abs(stage.beta_g_R)
        )
        if stage.M3_R <= 0.5:
            stage.delta_R = delta_0R
        else:
            X = 2 * stage.M3_R - 1
            stage.delta_R = delta_0R * (1 - 10*X**3 + 15*X**4 - 6*X**5)

    # =========================================================================
    #  Flow computations  (unchanged)
    # =========================================================================

    def computeVelTriangle(self):
        self.Vel_Tri['vu2OverU'] = (2*(1-self.R) + self.psi) / 2
        self.Vel_Tri['vu3OverU'] = (2*(1-self.R) - self.psi) / 2
        self.Vel_Tri['vmOverU']  = self.phi

        self.Vel_Tri['wu2OverU'] = self.Vel_Tri['vu2OverU'] - 1
        self.Vel_Tri['wu3OverU'] = self.Vel_Tri['vu3OverU'] - 1

        self.Vel_Tri['v2OverU'] = np.sqrt(self.Vel_Tri['vu2OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['w2OverU'] = np.sqrt(self.Vel_Tri['wu2OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['v3OverU'] = np.sqrt(self.Vel_Tri['vu3OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['w3OverU'] = np.sqrt(self.Vel_Tri['wu3OverU']**2 + self.Vel_Tri['vmOverU']**2)

        self.Vel_Tri['alpha1'] = self.Vel_Tri['alpha3'] = np.arctan(self.Vel_Tri['vu3OverU'] / self.Vel_Tri['vmOverU'])
        self.Vel_Tri['alpha2'] = np.arctan(self.Vel_Tri['vu2OverU'] / self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta1']  = self.Vel_Tri['beta3']  = np.arctan(self.Vel_Tri['wu3OverU'] / self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta2']  = np.arctan(self.Vel_Tri['wu2OverU'] / self.Vel_Tri['vmOverU'])

    def computeVelTriangleLastStage(self):
        self.Vel_Tri_Last_Stage['u']      = self.Vel_Tri['u']
        self.Vel_Tri_Last_Stage['vu2']    = 0
        self.Vel_Tri_Last_Stage['vu1']    = self.Vel_Tri['vu3']
        self.Vel_Tri_Last_Stage['vm']     = self.Vel_Tri['vm']
        self.Vel_Tri_Last_Stage['wu2']    = self.Vel_Tri_Last_Stage['vu2'] - self.Vel_Tri_Last_Stage['u']
        self.Vel_Tri['v2']                = np.sqrt(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w2']                = np.sqrt(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w3']                = np.sqrt(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri_Last_Stage['alpha1'] = self.Vel_Tri['alpha3']
        self.Vel_Tri_Last_Stage['alpha2'] = 0
        self.Vel_Tri_Last_Stage['beta1']  = self.Vel_Tri['beta3']
        self.Vel_Tri_Last_Stage['beta2']  = np.arctan(self.Vel_Tri['u'] / self.Vel_Tri['vm'])

    def computeBladeRow(self, stage_index, row_type):
        stage = self.stages[stage_index]
        self.curr_stage_index = stage_index
        ss = stage.static_states

        if 'P_ex' not in self.inputs:
            RP_1_row = 5 ** (1 / (2 * self.nStages))
        else:
            RP_1_row = (self.inputs['P_su'] / self.inputs['P_ex']) ** (1 / (2 * self.nStages))

        if row_type == 'S':
            h_out_guess = (ss[1]['H'] - self.Dh0_stage_guess / 2) if self.Dh0_stage_guess != 0 else ss[1]['H'] * 0.99
            pout_guess  = ss[1]['P'] / RP_1_row
            x0          = np.array([h_out_guess, pout_guess]) * 1e-5
            x_out       = self._wegstein_solve(self.stator_blade_row_system, x0)
            self.stator_blade_row_system(x_out)
        else:
            h_out_guess = (ss[2]['H'] - self.Dh0_stage_guess / 2) if self.Dh0_stage_guess != 0 else ss[2]['H'] * 0.99
            pout_guess  = ss[2]['P'] / RP_1_row
            x0          = np.array([h_out_guess, pout_guess]) * 1e-5
            x_out       = self._wegstein_solve(self.rotor_blade_row_system, x0)
            self.rotor_blade_row_system(x_out)
            self.compute_deviation_rotor(stage)

    def computeRepeatingStages(self):
        self.nStages = self.params['nStages']

        for i in range(int(self.nStages)):
            if i == 0:
                self.computeBladeRow(i, 'S')
                self.compute_deviation_stator(self.stages[i])
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']
                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']
                self.Dh0_stage_guess = (
                    self.stages[i].total_states[1]['H'] - self.stages[i].total_states[3]['H']
                )
            else:
                # Copy state dict by value — no DataFrame .loc needed
                self.stages[i].static_states[1] = dict(self.stages[i-1].static_states[3])
                self.computeBladeRow(i, 'S')
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']
                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

    def computeLastStage(self):
        stage = self.stages[-1]
        stage.static_states[1] = dict(self.stages[-2].static_states[3])

        if 'P_ex' not in self.inputs:
            RP_1_row = 5 ** (1 / (2 * self.nStages))
        else:
            RP_1_row = (self.inputs['P_su'] / self.inputs['P_ex']) ** (1 / (2 * self.nStages))

        h_out_guess = stage.static_states[1]['H'] - self.Dh0_stage_guess / 2
        pout_guess  = stage.static_states[1]['P'] / RP_1_row
        x0          = np.array([h_out_guess, pout_guess]) * 1e-5

        self._wegstein_solve(self.last_blade_row_system, x0)

    # =========================================================================
    #  Main solver
    # =========================================================================

    def solve(self, map_case=False):
        
        if self.params.get('solve_type') == "map":
            self.solve_from_map()
            return

        self.omega_rads = 2 * np.pi * self.inputs['N_rot'] / 60
        self.u          = self.omega_rads * self.params['r_m']

        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)

        self.computeRepeatingStages()
        self.computeLastStage()

        hin  = self.stages[0].total_states[1]['H']
        hout = self.stages[-1].static_states[2]['H']

        self.AS.update(
            CP.PSmass_INPUTS,
            self.stages[-1].static_states[2]['P'],
            self.stages[0].static_states[1]['S'],
        )
        hout_s = self.AS.hmass()

        self.W_dot    = self.inputs['m_dot'] * (hin - hout)
        self.eta_is   = (hin - hout) / (hin - hout_s)
        self._from_map = False
        self.solved    = True

        if not map_case:
            self.update_connectors(
                self.stages[-1].static_states[2]['P'],
                self.stages[-1].static_states[2]['H'],
            )

    def update_connectors(self, P_ex, h_ex):
        self.ex.reset()
        self.ex.set_p(P_ex)
        self.ex.set_h(h_ex)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_fluid(self.su.fluid)

        self.W.set_W_dot(self.W_dot)
        self.W.set_N_rot(self.inputs['N_rot'])

# =============================================================================
#  Parallel map generation helpers  (unchanged public API)
# =============================================================================

def _eval_point_from_snapshot(m, N, base_inputs, base_params, stage_params):
    try:
        turb = AxialTurbineMeanLine(base_inputs['fluid'])
        turb.set_inputs(
            m_dot = float(m),
            P_su  = base_inputs['P_su'],
            T_su  = base_inputs['T_su'],
            N_rot = float(N),
            fluid = base_inputs['fluid'],
            P_ex  = base_inputs['P_ex'],
        )
        turb.set_parameters(**base_params)
        if stage_params:
            turb.set_stage_parameters(**stage_params)

        turb.solve(map_case=True)

        P_ex_calc = float(turb.stages[-1].static_states[2]['P'])
        RP_calc   = turb.inputs['P_su'] / P_ex_calc if P_ex_calc else np.nan
        RP_target = turb.inputs['P_su'] / turb.inputs['P_ex'] if turb.inputs.get('P_ex') else np.nan

        if turb.eta_is < 0.3:
            return dict(
                m_dot=float(m), N_rot=float(N),
                P_su=float(base_inputs.get('P_su', np.nan)),
                T_su=float(base_inputs.get('T_su', np.nan)),
                P_ex_target=float(base_inputs.get('P_ex', np.nan)),
                P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
                W_dot=np.nan, eta_is=np.nan, converged=False,
            )

        return dict(
            m_dot       = float(m),
            N_rot       = float(N),
            P_su        = float(turb.inputs['P_su']),
            T_su        = float(turb.inputs['T_su']),
            P_ex_target = float(turb.inputs['P_ex']),
            P_ex_calc   = P_ex_calc,
            RP_target   = RP_target,
            RP_calc     = RP_calc,
            W_dot       = float(getattr(turb, 'W_dot', np.nan)),
            eta_is      = float(getattr(turb, 'eta_is', np.nan)),
            converged   = True,
            note        = "",
        )

    except Exception as e:
        return dict(
            m_dot=float(m), N_rot=float(N),
            P_su=float(base_inputs.get('P_su', np.nan)),
            T_su=float(base_inputs.get('T_su', np.nan)),
            P_ex_target=float(base_inputs.get('P_ex', np.nan)),
            P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
            W_dot=np.nan, eta_is=np.nan, converged=False, note=str(e),
        )


import os, sys
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
import joblib
from contextlib import contextmanager
from joblib.parallel import BatchCompletionCallBack


@contextmanager
def tqdm_joblib(tqdm_object):
    class TqdmBatchCompletionCallback(BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)
    old_cb = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_cb
        tqdm_object.close()


def generate_map_processes(machine, m_grid, N_grid, max_workers=-1, desc="Operation map"):
    os.environ.setdefault("OMP_NUM_THREADS",     "1")
    os.environ.setdefault("MKL_NUM_THREADS",     "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")

    base_inputs, base_params, stage_params = machine._snapshot_from_machine()
    tasks = [(m, N) for N in N_grid for m in m_grid]
    total = len(tasks)

    with tqdm(total=total, desc=desc, unit="pt",
              dynamic_ncols=True, miniters=1, mininterval=0,
              ascii=True, file=sys.stdout) as bar, tqdm_joblib(bar):
        results = Parallel(n_jobs=max_workers, backend="loky", prefer="processes")(
            delayed(_eval_point_from_snapshot)(m, N, base_inputs, base_params, stage_params)
            for (m, N) in tasks
        )

    return pd.DataFrame(results).sort_values(['N_rot', 'm_dot'], ignore_index=True)

if __name__ == "__main__":
    
    case_study = "TCO2_ORC"
    
    if case_study == "Salah_Case":
        Turb_OD = AxialTurbineMeanLine('CO2')
        
        Turb_OD.set_inputs(
              m_dot = 655.18,
              P_su = 25000000.0,
              T_su = 923,
              N_rot = 1506.9946780513785,
              fluid = 'CO2', 
              P_ex = 100*1e5
            )
        
        Turb_OD.set_parameters(
            r_m = 0.261423771889,
            nStages = 12,
            mdot_rated = 655.18,
            DP_rated = 2.5,
            damping = 0.5,
            delta_tip = 0.0004,
            N_lw = 0,
            D_lw = 0,
            e_blade = 2e-06
            )
        
        Turb_OD.set_stage_parameters(
            h_blade_S = [0.05893535333, 0.06254061127, 0.06644686421, 0.07068270152, 0.07527979324, 0.08027326565, 0.08570212832, 0.09160976073, 0.09804446776, 0.105060115, 0.1127168568, 0.1210819728, 0.1254195395],
            chord_S = [0.008645525688, 0.009066896935, 0.009520522543, 0.01000898242, 0.01053510451, 0.01110199118, 0.0117130487, 0.0123720203, 0.0130830232, 0.01385059029, 0.01467971707, 0.01557591465, 0.01603830059],
            xhi_S1 = [-0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816],
            xhi_S2 = [1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265],
            pitch_S = [0.006966193979, 0.007305716867, 0.007671228935, 0.008064808966, 0.008488735591, 0.008945508564, 0.009437872521, 0.009968843587, 0.01054173924, 0.0111602119, 0.01182828672, 0.01255040431, 0.01292297508],
            o_S = [0.001744270233, 0.001829283609, 0.001920804434, 0.002019353216, 0.002125500504, 0.002239872211, 0.00236315557, 0.002496105791, 0.002639553539, 0.002794413345, 0.00296169307, 0.003142504605, 0.003235792862],
            t_TE_S = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
            t_blade_S = [0.003025933991, 0.003173413927, 0.00333218289, 0.003503143847, 0.003687286579, 0.003885696912, 0.004099567045, 0.004330207105, 0.00457905812, 0.004847706601, 0.005137900974, 0.005451570129, 0.005613405206],
            n_blade_S = [236, 225, 214, 204, 194, 184, 174, 165, 156, 147, 139, 131, 127],
            R_c_S = [0.01352981312, 0.01418923794, 0.01489913922, 0.01566355437, 0.01648690899, 0.01737405813, 0.01833033244, 0.01936158987, 0.02047427367, 0.02167547758, 0.02297301932, 0.02437552349, 0.02509913423],
        
            h_blade_R = [0.06085456338, 0.06461845758, 0.06869830659, 0.07312430505, 0.07792992342, 0.08315231023, 0.08883274986, 0.09501718485, 0.1017568126, 0.1091087686, 0.1171369103, 0.1259127187, None],
            chord_R = [0.008870626618, 0.009309055087, 0.009781089003, 0.01028943789, 0.01083707309, 0.01142725573, 0.01206356808, 0.01274994866, 0.01349073167, 0.01429069144, 0.01515509243, 0.01608974592, None],
            xhi_R1 = [0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, None],
            xhi_R2 = [-1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, None],
            pitch_R = [0.007763880051, 0.008147607851, 0.008560748305, 0.009005672883, 0.009484982197, 0.01000153051, 0.01055845315, 0.01115919723, 0.01180755622, 0.01250770875, 0.01326426248, 0.01408230362, None],
            o_R = [0.00229560867, 0.002409068546, 0.002531225097, 0.002662779514, 0.002804500743, 0.002957232725, 0.003121902507, 0.003299529328, 0.003491234828, 0.00369825454, 0.003921950849, 0.004163827632, None],
            t_TE_R = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
            t_blade_R = [0.003104719316, 0.00325816928, 0.003423381151, 0.003301303262, 0.00379297558, 0.003999539504, 0.004222248827, 0.00446248203, 0.004721756085, 0.005001742002, 0.005304282352, 0.00563141107, None],
            n_blade_R = [212, 202, 192, 182, 173, 164, 156, 147, 139, 131, 124, 117, None],
            R_c_R = [0.01388208476, 0.01456820328, 0.01530691263, 0.01610245309, 0.01695947464, 0.0178830808, 0.01887887765, 0.01995302875, 0.02111231693, 0.02236421375, 0.02371695786, 0.02517964358, None],
            )
        
    elif case_study == "TCO2_ORC":
        Turb_OD = AxialTurbineMeanLine('CO2')
        
        Turb_OD.set_inputs(
            m_dot = 318.437,
            P_su = 15309670,
            T_su = 406.4,
            fluid = 'CO2', 
            N_rot = 2864.775003723627,
            P_ex = 5220928,
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
            'stator': {'h_blade_S': [0.024581131784973807,
              0.02731115991390685,
              0.030535241526643685,
              0.0343428949996286,
              0.03884110173212556,
              0.04415656370209263,
              0.05044235012970171,
              0.05396446226956401],
              'chord_S': [0.017387241421509592,
              0.018189646227485024,
              0.019157515021682635,
              0.020309513815586707,
              0.021666371275761282,
              0.023250883041521515,
              0.025089330758218745,
              0.026110623182803244],
              'xhi_S1': [-0.5493774236565352,
              -0.5493774236565352,
              -0.5493774236565352,
              -0.5493774236565352,
              -0.5493774236565352,
              -0.5493774236565352,
              -0.5493774236565352,
              -0.5493774236565352],
              'xhi_S2': [1.1268112059950612,
              1.1268112059950612,
              1.1268112059950612,
              1.1268112059950612,
              1.1268112059950612,
              1.1268112059950612,
              1.1268112059950612,
              1.1268112059950612],
              'pitch_S': [0.014118294382457852,
              0.014769840363216227,
              0.015555741716329778,
              0.01649115508675976,
              0.01759291197811346,
              0.01887952225855718,
              0.020372326403956538,
              0.021201607297417914],
              'o_S': [0.00407007938170583,
              0.004257909709554279,
              0.00448447255112079,
              0.004754137325719167,
              0.005071756287732223,
              0.0054426655373106475,
              0.005873017193716411,
              0.006112085665777279],
              't_TE_S': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
              't_blade_S': [0.0017387241421509593,
              0.0018189646227485025,
              0.0019157515021682636,
              0.002030951381558671,
              0.0021666371275761284,
              0.0023250883041521517,
              0.0025089330758218745,
              0.0026110623182803248],
              'n_blade_S': [93, 89, 85, 80, 75, 70, 65, 62],
              'R_c_S': [0.02593658923532629,
              0.02713353838605762,
              0.02857731056016104,
              0.03029574989120514,
              0.03231977738022937,
              0.034683397336424615,
              0.03742581415244761,
              0.03894927848261039]},
            'rotor': {'h_blade_R': [0.02599623192810525,
              0.028975662485097392,
              0.03249353115560146,
              0.03664835013019067,
              0.041556335411602884,
              0.04735747516384924,
              0.054218932168004906,
              None],
              'chord_R': [0.017804442744401783,
              0.018691926981905842,
              0.019754111574023726,
              0.0210108023910845,
              0.02248355410347207,
              0.02419698360942641,
              0.02617885482072475,
              None],
              'xhi_R1': [0.5033968740474793,
              0.5033968740474793,
              0.5033968740474793,
              0.5033968740474793,
              0.5033968740474793,
              0.5033968740474793,
              0.5033968740474793,
              None],
              'xhi_R2': [-1.1544956444218015,
              -1.1544956444218015,
              -1.1544956444218015,
              -1.1544956444218015,
              -1.1544956444218015,
              -1.1544956444218015,
              -1.1544956444218015,
              None],
              'pitch_R': [0.015148494029818542,
              0.015903589253319174,
              0.01680732419089942,
              0.017876550204477332,
              0.01912960658162921,
              0.02058743803494024,
              0.022273666839919944,
              None],
              'o_R': [0.004776675516045738,
              0.005014774752793548,
              0.005299743578132708,
              0.005636895621769747,
              0.006032013691274512,
              0.006491702145839364,
              0.007023409643056108,
              None],
              't_TE_R': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
              't_blade_R': [0.0019275069515591894,
              0.0020395300073626607,
              0.002174650407702253,
              0.0022904090712424882,
              0.0024358345376652743,
              0.002590844225248941,
              0.002758858045024202,
              None],
              'n_blade_R': [87, 83, 78, 74, 69, 64, 59, None],
              'R_c_R': [0.02655892943743007,
              0.02788279177780649,
              0.029467255051191132,
              0.03134186372129275,
              0.03353877094100647,
              0.036094697795777626,
              0.03905105977860335,
              None]}}
        
        Turb_OD.set_parameters(
            r_m = turb_params['r_m'],
            nStages = turb_params['n_stages'],
            mdot_rated = turb_params['mdot_rated'],
            DP_rated = turb_params['DP_rated'],
            N_rot_rated = turb_params['N_rot_rated'], # RPM 
            damping = 0.3,
            delta_tip = 0.0004,
            N_lw = 0,
            D_lw = 0,
            e_blade = 2e-06
            )
        
        Turb_OD.set_stage_parameters(
            # --- Stator ---
            h_blade_S  = turb_params['stator']['h_blade_S'],
            chord_S    = turb_params['stator']['chord_S'],
            xhi_S1     = turb_params['stator']['xhi_S1'],
            xhi_S2     = turb_params['stator']['xhi_S2'],
            pitch_S    = turb_params['stator']['pitch_S'],
            o_S        = turb_params['stator']['o_S'],
            t_TE_S     = turb_params['stator']['t_TE_S'],
            t_blade_S  = turb_params['stator']['t_blade_S'],
            n_blade_S  = turb_params['stator']['n_blade_S'],
            R_c_S      = turb_params['stator']['R_c_S'],
            # --- Rotor ---
            h_blade_R  = turb_params['rotor']['h_blade_R'],
            chord_R    = turb_params['rotor']['chord_R'],
            xhi_R1     = turb_params['rotor']['xhi_R1'],
            xhi_R2     = turb_params['rotor']['xhi_R2'],
            pitch_R    = turb_params['rotor']['pitch_R'],
            o_R        = turb_params['rotor']['o_R'],
            t_TE_R     = turb_params['rotor']['t_TE_R'],
            t_blade_R  = turb_params['rotor']['t_blade_R'],
            n_blade_R  = turb_params['rotor']['n_blade_R'],
            R_c_R      = turb_params['rotor']['R_c_R'],
        )
        
        Turb_OD.solve()
        
    # -------------------------------------------------------------------------
    # 1) Generate the performance map
    # -------------------------------------------------------------------------
    df_map = generate_map_processes(
        Turb_OD,
        m_grid=np.linspace(0.6*Turb_OD.params['mdot_rated'], 1.4*Turb_OD.params['mdot_rated'], 30),
        N_grid=np.linspace(0.3*Turb_OD.params['N_rot_rated'], 1.5*Turb_OD.params['N_rot_rated'], 30),
        max_workers=-2
    )

    # -------------------------------------------------------------------------
    # 2) Save the raw DataFrame, then reload it from disk
    #    → set MAP_SAVE_PATH to whatever directory you want
    # -------------------------------------------------------------------------
    MAP_SAVE_PATH = r"C:\Users\basil\OneDrive\Desktop\Travail\Thèse\WP1\Turbine Maps\turb_map.parquet"   # <-- set your path here

    Turb_OD.save_map_df(MAP_SAVE_PATH, df_map)
    df_reloaded = Turb_OD.load_map_df(MAP_SAVE_PATH)   # also calls load_map()

    assert len(df_reloaded) == len(df_map), "Row count mismatch after reload!"
    print(f"[check] DataFrame round-trip OK  ({len(df_reloaded)} rows)")

    # -------------------------------------------------------------------------
    # 3) Interpolation validation: raw map vs interpolated reconstruction
    #
    #    Strategy: build a fine regular grid that covers the convex hull of
    #    the converged points, query the interpolator on every node, then plot
    #    raw data (scatter) and interpolated surface (filled contour) side by
    #    side for W_dot, eta_is and P_ex_calc.
    # -------------------------------------------------------------------------
    df_ok = df_reloaded[df_reloaded["converged"] == True].dropna(
        subset=["W_dot", "eta_is", "P_ex_calc"]
    )

    # Fine query grid (denser than the original map)
    N_GRID  = 120
    m_fine  = np.linspace(df_ok["m_dot"].min(),  df_ok["m_dot"].max(),  N_GRID)
    N_fine  = np.linspace(df_ok["N_rot"].min(),  df_ok["N_rot"].max(),  N_GRID)
    MM, NN  = np.meshgrid(m_fine, N_fine)           # shape (N_GRID, N_GRID)
    m_flat  = MM.ravel()
    N_flat  = NN.ravel()

    # Query the interpolator on the grid (NaN outside hull → shown as white)
    interp  = Turb_OD.map_interpolator
    res     = interp.query_batch(m_flat, N_flat)

    W_grid   = res["W_dot"].reshape(N_GRID, N_GRID)   / 1e6   # → MW
    eta_grid = res["eta_is"].reshape(N_GRID, N_GRID)
    Pex_grid = res["P_ex_calc"].reshape(N_GRID, N_GRID) / 1e5  # → bar

    # Mask out-of-hull cells so they don't pollute colour scaling
    hull_mask = interp.hull_mask(m_flat, N_flat).reshape(N_GRID, N_GRID)
    W_grid   = np.where(hull_mask, W_grid,   np.nan)
    eta_grid = np.where(hull_mask, eta_grid, np.nan)
    Pex_grid = np.where(hull_mask, Pex_grid, np.nan)

    # ---- figure layout: 3 rows × 2 columns ---------------------------------
    #   col 0 = raw scattered data   |   col 1 = interpolated contour
    TARGETS = [
        ("W_dot [MW]",    df_ok["W_dot"]    / 1e6,  W_grid,   "viridis"),
        ("eta_is [–]",    df_ok["eta_is"],           eta_grid, "plasma"),
        ("P_ex [bar]",    df_ok["P_ex_calc"] / 1e5,  Pex_grid, "cividis"),
    ]
    N_LEVELS = 20

    fig, axes = plt.subplots(
        nrows=3, ncols=2,
        figsize=(13, 13),
        constrained_layout=True,
    )
    fig.suptitle(
        "Performance map  –  raw data (left)  vs  interpolated reconstruction (right)",
        fontsize=13, fontweight="bold"
    )

    for row, (label, raw_vals, grid_vals, cmap) in enumerate(TARGETS):

        vmin = float(np.nanmin(grid_vals))
        vmax = float(np.nanmax(grid_vals))
        levels = np.linspace(vmin, vmax, N_LEVELS + 1)

        # ---- left: scatter of raw converged points -------------------------
        ax_raw = axes[row, 0]
        sc = ax_raw.scatter(
            df_ok["m_dot"], df_ok["N_rot"],
            c=raw_vals, cmap=cmap,
            vmin=vmin, vmax=vmax,
            s=40, edgecolors="k", linewidths=0.4, zorder=3
        )
        ax_raw.set_title(f"{label}  –  raw map points", fontsize=10)
        ax_raw.set_xlabel("ṁ  [kg/s]")
        ax_raw.set_ylabel("N  [rpm]")
        fig.colorbar(sc, ax=ax_raw, label=label, pad=0.02)

        # ---- right: filled contour of interpolated grid --------------------
        ax_int = axes[row, 1]
        cf = ax_int.contourf(
            MM, NN, grid_vals,
            levels=levels, cmap=cmap, extend="neither"
        )
        # overlay iso-lines
        ax_int.contour(
            MM, NN, grid_vals,
            levels=levels, colors="white", linewidths=0.35, alpha=0.5
        )
        # overlay original scatter so position accuracy is visible
        ax_int.scatter(
            df_ok["m_dot"], df_ok["N_rot"],
            c=raw_vals, cmap=cmap,
            vmin=vmin, vmax=vmax,
            s=25, edgecolors="k", linewidths=0.5, zorder=4
        )
        ax_int.set_title(
            f"{label}  –  interpolated  "
            f"(method={interp.method!r}, n={interp.n_points} pts)",
            fontsize=10
        )
        ax_int.set_xlabel("ṁ  [kg/s]")
        ax_int.set_ylabel("N  [rpm]")
        fig.colorbar(cf, ax=ax_int, label=label, pad=0.02)

    plt.show()

    # ---- point-by-point residual table  ------------------------------------
    #   Re-query the interpolator exactly at the raw map nodes and compare.
    res_nodes = interp.query_batch(df_ok["m_dot"].to_numpy(),
                                    df_ok["N_rot"].to_numpy())

    err_W   = (res_nodes["W_dot"]     - df_ok["W_dot"].to_numpy())   / 1e6
    err_eta = (res_nodes["eta_is"]    - df_ok["eta_is"].to_numpy())
    err_Pex = (res_nodes["P_ex_calc"] - df_ok["P_ex_calc"].to_numpy()) / 1e5

    print("\n── Interpolation residuals at map nodes ──────────────────────────────")
    print(f"  W_dot    :  max |err| = {np.abs(err_W).max():.4f} MW   "
          f"  RMS = {np.sqrt((err_W**2).mean()):.4f} MW")
    print(f"  eta_is   :  max |err| = {np.abs(err_eta).max():.6f}    "
          f"  RMS = {np.sqrt((err_eta**2).mean()):.6f}")
    print(f"  P_ex_calc:  max |err| = {np.abs(err_Pex).max():.4f} bar "
          f"  RMS = {np.sqrt((err_Pex**2).mean()):.4f} bar")
    print("──────────────────────────────────────────────────────────────────────")

    # -------------------------------------------------------------------------
    # 4) Standard map plot (unchanged workflow)
    # -------------------------------------------------------------------------
    df_clean = filter_sparse_by_proximity(df_map, rp_col=None, group_by='N_rot',
                                          rp_tol_rel=0.6, m_tol_rel=0.6, min_neighbors=2)

    fig, ax = map_plot(
        df_clean, rp_col='RP_calc',
        use_grid=True, nx=600, ny=600,
        min_circle_ratio=0.01,
        max_area_factor=50.0,
        long_edge_q=1,
        fill_holes=True, hole_method='nearest', hole_smooth_sigma=0.8,
        smooth_sigma=0.6,
        show_points=True,
        levels=24, focus_high=True, max_iso_speeds=4,
        figsize=(9,6), dpi=220
    )
    plt.show()

    _ = plot_power_eta_vs_mdot(df_map, speeds=None, max_lines=5)
    plt.show()
    
    

