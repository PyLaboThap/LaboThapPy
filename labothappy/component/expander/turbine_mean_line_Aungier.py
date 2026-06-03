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

        self.Dh0_stage_guess  = 0
        self.map_interpolator = None

    # =========================================================================
    #  Inner class: stage
    # =========================================================================

    class stage:

        _STATE_KEYS = ('H', 'S', 'P', 'D', 'A', 'V')

        def __init__(self, fluid):
            self.AS = CP.AbstractState('HEOS', fluid)

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

        def _write_state(self, target_dict, position):
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
    #  Inner class: MapInterpolator
    # =========================================================================

    class MapInterpolator:
        _DEFAULT_TARGETS = ["W_dot", "eta_is", "P_ex_calc"]

        def __init__(
            self,
            df,
            *,
            method             = "linear",
            m_dot_col          = "m_dot",
            N_rot_col          = "N_rot",
            targets            = None,
            normalize          = True,
            fallback_nn        = True,
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

            # ---- forward interpolants: (m_dot, N_rot) -> targets ------------
            pts_fwd = self._normalize_2d(m_raw, N_raw, self._m_scale, self._N_scale)

            self._primary = {}
            self._nn      = {}

            for col in self.targets:
                if col not in df_ok.columns:
                    warnings.warn(f"MapInterpolator: target column {col!r} not found – skipping.")
                    continue
                vals = df_ok[col].to_numpy(dtype=float)
                self._primary[col] = self._build_interp(method, pts_fwd, vals)
                self._nn[col]      = NearestNDInterpolator(pts_fwd, vals)

            try:
                self._hull = Delaunay(pts_fwd)
            except Exception:
                self._hull = None

            # ---- inverse interpolants for P_N and P_M modes -----------------
            self._inv_primary = {}
            self._inv_nn      = {}
            self._inv_hull    = {}

            if "P_ex_calc" in df_ok.columns:
                P_raw         = df_ok["P_ex_calc"].to_numpy(dtype=float)
                self._P_scale = (float(P_raw.min()), float(P_raw.max()))

                # P_N mode: (N_rot, DP=P_su-P_ex) -> m_dot  [DP replaces P_ex as input]
                if "DP_calc" in df_ok.columns:
                    DP_raw         = df_ok["DP_calc"].to_numpy(dtype=float)
                    self._DP_scale = (float(DP_raw.min()), float(DP_raw.max()))
                else:
                    DP_raw         = df_ok["P_su"].to_numpy(dtype=float) - P_raw
                    self._DP_scale = (float(DP_raw.min()), float(DP_raw.max()))

                pts_PN = self._normalize_2d(N_raw, DP_raw, self._N_scale, self._DP_scale)
                self._inv_primary["PN_to_m"] = self._build_interp(method, pts_PN, m_raw)
                self._inv_nn["PN_to_m"]      = NearestNDInterpolator(pts_PN, m_raw)
                try:    self._inv_hull["PN"] = Delaunay(pts_PN)
                except: self._inv_hull["PN"] = None

                # P_M mode unchanged: (m_dot, P_ex) -> N_rot
                pts_PM = self._normalize_2d(m_raw, P_raw, self._m_scale, self._P_scale)
                self._inv_primary["PM_to_N"] = self._build_interp(method, pts_PM, N_raw)
                self._inv_nn["PM_to_N"]      = NearestNDInterpolator(pts_PM, N_raw)
                try:    self._inv_hull["PM"] = Delaunay(pts_PM)
                except: self._inv_hull["PM"] = None

            else:
                self._P_scale  = (0.0, 1.0)
                self._DP_scale = (0.0, 1.0)

            self._pts_raw = np.column_stack([m_raw, N_raw])
            self.n_points = len(df_ok)
            self.m_range  = self._m_scale
            self.N_range  = self._N_scale
            self.P_range  = self._P_scale
            self.DP_range = self._DP_scale

        # ------------------------------------------------------------------
        #  Private helpers
        # ------------------------------------------------------------------

        @staticmethod
        def _build_interp(method, pts, vals):
            if method == "linear":
                return LinearNDInterpolator(pts, vals)
            elif method == "rbf":
                return RBFInterpolator(pts, vals, kernel="thin_plate_spline", smoothing=0.0)
            else:
                return NearestNDInterpolator(pts, vals)

        @staticmethod
        def _normalize_2d(a, b, scale_a, scale_b):
            a = np.asarray(a, float)
            b = np.asarray(b, float)
            a_n = (a - scale_a[0]) / max(scale_a[1] - scale_a[0], 1e-12)
            b_n = (b - scale_b[0]) / max(scale_b[1] - scale_b[0], 1e-12)
            return np.column_stack([a_n, b_n])

        def _normalize(self, m, N):
            return self._normalize_2d(m, N, self._m_scale, self._N_scale)

        def _query_interp(self, primary, nn, pts):
            vals     = np.asarray(primary(pts), dtype=float).ravel()
            nan_mask = ~np.isfinite(vals)
            if nan_mask.any():
                if self.fallback_nn:
                    vals[nan_mask] = nn(pts[nan_mask]).ravel()
                else:
                    raise MapInterpolationError("Query point(s) outside hull.")
            return vals

        # ------------------------------------------------------------------
        #  Public query API
        # ------------------------------------------------------------------

        def hull_mask(self, m_arr, N_arr):
            pts = self._normalize(np.asarray(m_arr, float), np.asarray(N_arr, float))
            if self._hull is None:
                return np.ones(len(pts), dtype=bool)
            return self._hull.find_simplex(pts) >= 0

        def query(self, m_dot, N_rot):
            """M_N mode: (m_dot, N_rot) -> all targets."""
            res = self.query_batch(np.array([m_dot], float), np.array([N_rot], float))
            return {k: float(v[0]) for k, v in res.items()}

        def query_batch(self, m_arr, N_arr):
            """M_N mode vectorised."""
            pts = self._normalize(np.asarray(m_arr, float), np.asarray(N_arr, float))
            out = {}
            for col, interp in self._primary.items():
                out[col] = self._query_interp(interp, self._nn[col], pts)
            return out

        def query_PN(self, N_rot, P_ex, P_su):
            """P_N mode: (N_rot, P_ex) -> m_dot, W_dot, eta_is.
            Internally uses DP = P_su - P_ex as the interpolation coordinate.
            """
            DP    = P_su - P_ex
            pts   = self._normalize_2d(
                np.array([N_rot], float), np.array([DP], float),
                self._N_scale, self._DP_scale,
            )
            m_dot = float(self._query_interp(
                self._inv_primary["PN_to_m"], self._inv_nn["PN_to_m"], pts
            )[0])
            fwd = self.query(m_dot, N_rot)
            return {"m_dot": m_dot, "W_dot": fwd["W_dot"],
                    "eta_is": fwd["eta_is"], "P_ex_calc": P_ex}

        def query_PM(self, m_dot, P_ex):
            """P_M mode: (m_dot, P_ex) -> N_rot, W_dot, eta_is."""
            pts   = self._normalize_2d(
                np.array([m_dot], float), np.array([P_ex], float),
                self._m_scale, self._P_scale,
            )
            N_rot = float(self._query_interp(
                self._inv_primary["PM_to_N"], self._inv_nn["PM_to_N"], pts
            )[0])
            fwd = self.query(m_dot, N_rot)
            return {"N_rot": N_rot, "W_dot": fwd["W_dot"],
                    "eta_is": fwd["eta_is"], "P_ex_calc": P_ex}

        # ------------------------------------------------------------------
        #  Persistence
        # ------------------------------------------------------------------

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
                f"P_range={self.P_range}, "
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
        base_params  = dict(self.params)
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

    # =========================================================================
    #  Stage solving
    # =========================================================================

    def _wegstein_solve(self, f, x0, tol=1e-6, max_iter=100, q_min=-5.0, q_max=0.0):
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

    # =========================================================================
    #  Map interface
    # =========================================================================

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
                f"| m in [{self.map_interpolator.m_range[0]:.1f}, "
                f"{self.map_interpolator.m_range[1]:.1f}] kg/s  "
                f"| N in [{self.map_interpolator.N_range[0]:.0f}, "
                f"{self.map_interpolator.N_range[1]:.0f}] rpm  "
                f"| P in [{self.map_interpolator.P_range[0]/1e5:.2f}, "
                f"{self.map_interpolator.P_range[1]/1e5:.2f}] bar"
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

    def solve_from_map(self, m_dot=None, N_rot=None, P_ex=None, *,
                       mode="M_N", strict_hull=False):
        """
        Resolve turbine from interpolated performance map.

        Modes
        -----
        M_N : (m_dot, N_rot) → P_ex, W_dot, eta_is   [default]
        P_N : (N_rot, P_ex)  → m_dot, W_dot, eta_is
        P_M : (m_dot, P_ex)  → N_rot, W_dot, eta_is
        """
        if self.map_interpolator is None:
            raise MapInterpolationError(
                "No performance map loaded. Call turbine.load_map(df_map) first."
            )

        if mode == "M_N":
            _m = float(m_dot if m_dot is not None else self.inputs["m_dot"])
            _N = float(N_rot if N_rot is not None else self.inputs["N_rot"])
            if strict_hull:
                inside = self.map_interpolator.hull_mask(np.array([_m]), np.array([_N]))[0]
                if not inside:
                    raise MapInterpolationError(
                        f"Operating point (m_dot={_m:.3f}, N_rot={_N:.1f}) outside hull."
                    )
            result = self.map_interpolator.query(_m, _N)
            
            self.inputs["m_dot"] = _m
            self.inputs["N_rot"] = _N
            self.m_dot_ex = _m
            
        elif mode == "P_N":
            _N   = float(N_rot if N_rot is not None else self.inputs["N_rot"])
            _P   = float(P_ex  if P_ex  is not None else self.inputs.get("P_ex", float("nan")))
            _Psu = float(self.inputs.get("P_su", float("nan")))
            result = self.map_interpolator.query_PN(_N, _P, _Psu)   # <-- add _Psu
            
            self.inputs["N_rot"] = _N
            self.m_dot_ex = result["m_dot"]
            
        elif mode == "P_M":
            _m = float(m_dot if m_dot is not None else self.inputs["m_dot"])
            _P = float(P_ex  if P_ex  is not None else self.inputs.get("P_ex", float("nan")))
            result = self.map_interpolator.query_PM(_m, _P)

            self.inputs["m_dot"] = _m
            self.W.set_N_rot(result["N_rot"])
            self.m_dot_ex = _m
            
        else:
            raise ValueError(f"mode must be 'M_N', 'P_N', or 'P_M'; got {mode!r}")

        self.eta_is    = result.get("eta_is",    float("nan"))
        self.P_ex_calc = result.get("P_ex_calc", float("nan"))

        if self.ex is not None:
            try:
                self.ex.p = self.P_ex_calc
            except AttributeError:
                pass

        self._from_map = True
        self.AS.update(CP.PSmass_INPUTS, self.P_ex_calc, self.su.s)
        h_ex_is   = self.AS.hmass()
        self.h_ex = self.su.h - self.eta_is * (self.su.h - h_ex_is)
        self.W_dot     = self.m_dot_ex * (self.su.h - self.h_ex)
        
        self.solved = True
        self.update_connectors(self.P_ex_calc, self.h_ex, self.m_dot_ex)

    # =========================================================================
    #  Loss Models
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
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(arg_v2)
        stage.M2_S = v2 / ss[2]['A']
        stage.Vel_Tri_S['alpha2'] = np.arctan2(np.sqrt(max(0.0, v2**2 - vm**2)), vm)

        stage.Re_s     = stage.chord_S * (ss[2]['D'] * vm) / ss[2]['V']
        stage.AR_S     = stage.h_blade_S / stage.chord_S
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

        Y      = stage.Y_vec_S['Y_tot']
        p0_out = (ts[1]['P'] + Y * p_static_out) / (1 + Y)
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)

        sout = ts[2]['S']
        hout = ts[2]['H'] - v2_new**2 / 2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)

        pout_calc = ss[2]['P']
        self.AS.update(CP.PSmass_INPUTS, pout_calc, ss[1]['S'])
        hout_s    = self.AS.hmass()
        stage.eta_is_S = (ss[1]['H'] - ss[2]['H']) / (ss[1]['H'] - hout_s)

        return np.array([hout, pout_calc]) * 1e-5

    def rotor_blade_row_system(self, x):
        stage = self.stages[self.curr_stage_index]
        ss    = stage.static_states
        ts    = stage.total_states

        h_static_out = x[0] * 1e5
        p_static_out = x[1] * 1e5

        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)

        stage.Vel_Tri_R['u']   = u   = self.u
        A_flow                       = stage.h_blade_R * (2 * np.pi * self.params['r_m'])
        stage.Vel_Tri_R['vm']  = vm  = self.inputs['m_dot'] / (ss[3]['D'] * A_flow)
        stage.Vel_Tri_R['vu2'] = vu2 = vm * np.tan(stage.Vel_Tri_R['alpha2'])

        wu2 = vu2 - u
        stage.Vel_Tri_R['wu2']  = wu2
        stage.Vel_Tri_R['w2']   = w2 = np.sqrt(wu2**2 + vm**2)
        stage.M2_R               = w2 / ss[2]['A']
        stage.Vel_Tri_R['beta2'] = np.arctan(wu2 / vm)

        hin  = ss[2]['H']
        h0in = hin + w2**2 / 2
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, ss[2]['S'], 2)

        h03    = ts[2]['H']
        arg_w3 = max(0.0, 2 * (h03 - h_static_out))
        stage.Vel_Tri_R['w3'] = w3 = np.sqrt(arg_w3)
        stage.M3_R = w3 / ss[3]['A']
        stage.Vel_Tri_R['beta3'] = -np.arccos(np.clip(vm / (w3 + 1e-30), -1.0, 1.0))

        stage.Re_r     = stage.chord_R * (ss[3]['D'] * vm) / ss[3]['V']
        stage.AR_R     = stage.h_blade_R / stage.chord_R
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

        h02    = h0in
        arg_v2 = max(0.0, 2 * (h02 - h_static_out))
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(arg_v2)
        stage.M2_S = v2 / ss[2]['A']
        stage.Vel_Tri_S['alpha2'] = np.arctan2(np.sqrt(max(0.0, v2**2 - vm**2)), vm)

        alpha2   = stage.Vel_Tri_S['alpha2']
        solidity = stage.chord_S / stage.pitch_S
        a        = 0.0117

        D_e = (np.cos(alpha2) / np.cos(alpha1)) * (
            1.12
            + a * (alpha1 - stage.xhi_S1)
            + 0.61 * np.cos(alpha1)**2 / solidity * (np.tan(alpha1) - np.tan(alpha2))
        )

        P_cst  = np.cos(alpha2) / 2 * solidity * (v1 / (v2 + 1e-30))**2
        Yp     = 0.004 * (1 + 3.1 * (D_e - 1)**2 + 0.4 * (D_e - 1)**8) / (P_cst + 1e-30)
        EW_Cst = np.cos((alpha1 + alpha2) / 2)**3 / np.cos(alpha1)**2
        AR_S   = stage.h_blade_S / stage.chord_S
        Yew    = 0.02 * (solidity / AR_S) / (EW_Cst + 1e-30)

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
    #  Deviation models
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
    #  Flow computations
    # =========================================================================

    def computeVelTriangle(self):
        self.Vel_Tri['vu2OverU'] = (2*(1-self.R) + self.psi) / 2
        self.Vel_Tri['vu3OverU'] = (2*(1-self.R) - self.psi) / 2
        self.Vel_Tri['vmOverU']  = self.phi
        self.Vel_Tri['wu2OverU'] = self.Vel_Tri['vu2OverU'] - 1
        self.Vel_Tri['wu3OverU'] = self.Vel_Tri['vu3OverU'] - 1
        self.Vel_Tri['v2OverU']  = np.sqrt(self.Vel_Tri['vu2OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['w2OverU']  = np.sqrt(self.Vel_Tri['wu2OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['v3OverU']  = np.sqrt(self.Vel_Tri['vu3OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['w3OverU']  = np.sqrt(self.Vel_Tri['wu3OverU']**2 + self.Vel_Tri['vmOverU']**2)
        self.Vel_Tri['alpha1']   = self.Vel_Tri['alpha3'] = np.arctan(self.Vel_Tri['vu3OverU'] / self.Vel_Tri['vmOverU'])
        self.Vel_Tri['alpha2']   = np.arctan(self.Vel_Tri['vu2OverU'] / self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta1']    = self.Vel_Tri['beta3']  = np.arctan(self.Vel_Tri['wu3OverU'] / self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta2']    = np.arctan(self.Vel_Tri['wu2OverU'] / self.Vel_Tri['vmOverU'])

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
            self.solve_from_map(mode = self.params['map_mode'])
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

        self.W_dot     = self.inputs['m_dot'] * (hin - hout)
        self.eta_is    = (hin - hout) / (hin - hout_s)
        self._from_map = False
        self.solved    = True

        if not map_case:
            self.update_connectors(
                self.stages[-1].static_states[2]['P'],
                self.stages[-1].static_states[2]['H'],
                self.su.m_dot,
            )

    def update_connectors(self, P_ex, h_ex, m_dot_ex):
        
        self.ex.reset()
        
        self.ex.set_p(P_ex)
        self.ex.set_h(h_ex)
        self.ex.set_m_dot(m_dot_ex)
        self.ex.set_fluid(self.su.fluid)
        self.W.set_W_dot(self.W_dot)
        self.W.set_N_rot(self.inputs['N_rot'])


# =============================================================================
#  Parallel map generation helpers
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
                DP_calc = np.nan,
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
            DP_calc     = float(turb.inputs['P_su']) - P_ex_calc,
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
            DP_calc = np.nan,
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
    os.environ.setdefault("OMP_NUM_THREADS",      "1")
    os.environ.setdefault("MKL_NUM_THREADS",      "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS",  "1")

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
