from labothappy.connector.mass_connector import MassConnector
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

# =============================================================================
#  MapInterpolationError  –  module-level so it can be caught without
#  instantiating the turbine class first.
# =============================================================================

class MapInterpolationError(RuntimeError):
    """Raised when map-based interpolation cannot produce a valid result."""


# =============================================================================
#  AxialTurbineMeanLine
# =============================================================================

class AxialTurbineMeanLine(BaseComponent):

    def __init__(self, fluid):
        super().__init__()
        
        # Inputs
        self.inputs = {}
        
        # Params
        self.params = {}  

        # Abstract State 
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Blade Dictionnary
        self.stages = []

        # Velocity Triangle Data
        self.Vel_Tri_Last_Stage = {}

        self.su = MassConnector()
        self.ex = MassConnector()
        
        self.Dh0_stage_guess = 0

        # Map-interpolation slot (populated by load_map / set_map)
        self.map_interpolator = None

    # =========================================================================
    #  Inner class: stage
    # =========================================================================
    
    class stage(object):
        
        def __init__(self, fluid):
            self.total_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3])
            self.static_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3])
            self.AS = CP.AbstractState('HEOS', fluid)
            
            self.eta_is_R = None
            self.eta_is_S = None
            
            self.A_flow_S = None
            self.A_flow_R = None
            
            self.h_blade_S = None
            self.h_blade_R = None
            
            self.chord_S = None
            self.chord_R = None
            
            self.stage = None
            self.AR = None
            
            self.xhi_S1 = None
            self.xhi_S2 = None
            
            self.xhi_R1 = None
            self.xhi_R2 = None
            
            self.Vel_Tri_R = {}
            self.Vel_Tri_S = {}
            
        def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.total_states.loc[position, 'H'] = self.AS.hmass()
            self.total_states.loc[position, 'S'] = self.AS.smass()
            self.total_states.loc[position, 'P'] = self.AS.p()
            self.total_states.loc[position, 'D'] = self.AS.rhomass()
            try:
                self.total_states.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                self.total_states.loc[position, 'A'] = -1
            self.total_states.loc[position, 'V'] = self.AS.viscosity()        
            
            return
        
        def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.static_states.loc[position, 'H'] = self.AS.hmass()
            self.static_states.loc[position, 'S'] = self.AS.smass()
            self.static_states.loc[position, 'P'] = self.AS.p()
            self.static_states.loc[position, 'D'] = self.AS.rhomass()
            try:
                self.static_states.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                self.static_states.loc[position, 'A'] = -1
            self.static_states.loc[position, 'V'] = self.AS.viscosity()         

            return

    # =========================================================================
    #  Inner class: MapInterpolator
    # =========================================================================

    class MapInterpolator:
        """
        2-D scattered interpolant over a turbine performance map.

        Build one from any DataFrame produced by ``generate_map_processes``
        (or ``generate_map_m_dot_N_rot``), then attach it to the turbine with
        ``turbine.load_map(df)`` or assign it directly via
        ``turbine.map_interpolator = turbine.MapInterpolator(df)``.

        Parameters
        ----------
        df : pd.DataFrame
            Must contain ``m_dot_col``, ``N_rot_col``, and all columns
            listed in ``targets``.
        method : {'linear', 'rbf', 'nearest'}
            Interpolation strategy:

            * ``'linear'``  – ``LinearNDInterpolator`` (Delaunay triangulation).
              Exact at data nodes, returns NaN outside the convex hull
              (handled automatically when ``fallback_nn=True``).
            * ``'rbf'``     – ``RBFInterpolator`` with thin-plate-spline kernel.
              Smooth, extrapolates beyond the hull, slightly slower to build.
            * ``'nearest'`` – ``NearestNDInterpolator``. Fastest, always
              returns a value.

        targets : list[str], optional
            Columns to interpolate.
            Default: ``['W_dot', 'eta_is', 'P_ex_calc']``.
        normalize : bool
            Scale inputs to [0, 1] before fitting (strongly recommended for
            ``'rbf'`` where coordinate magnitudes matter).
        fallback_nn : bool
            Fall back silently to nearest-neighbour when the primary
            interpolant returns NaN (i.e. out-of-hull for ``'linear'``).
        min_converged_frac : float
            Fraction of converged rows required; raises
            ``MapInterpolationError`` if too few points converged.
        """

        _DEFAULT_TARGETS = ["W_dot", "eta_is", "P_ex_calc"]

        def __init__(
            self,
            df: pd.DataFrame,
            *,
            method: str = "linear",
            m_dot_col: str = "m_dot",
            N_rot_col: str = "N_rot",
            targets=None,
            normalize: bool = True,
            fallback_nn: bool = True,
            min_converged_frac: float = 0.2,
        ):
            if method not in ("linear", "rbf", "nearest"):
                raise ValueError(
                    f"method must be 'linear', 'rbf', or 'nearest'; got {method!r}"
                )

            self.method          = method
            self.m_dot_col       = m_dot_col
            self.N_rot_col       = N_rot_col
            self.targets         = targets if targets is not None else list(self._DEFAULT_TARGETS)
            self.normalize       = normalize
            self.fallback_nn     = fallback_nn

            # ---- keep only converged rows -----------------------------------
            conv_col = df.get("converged", pd.Series(True, index=df.index))
            df_ok = df[conv_col == True].copy()

            frac = len(df_ok) / max(len(df), 1)
            if frac < min_converged_frac:
                raise MapInterpolationError(
                    f"Only {frac:.1%} of map points converged "
                    f"(threshold {min_converged_frac:.1%}). "
                    "Generate a denser map before building an interpolator."
                )

            df_ok = df_ok.dropna(subset=[m_dot_col, N_rot_col] + self.targets)
            if len(df_ok) < 4:
                raise MapInterpolationError(
                    f"Only {len(df_ok)} usable points after filtering – need at least 4."
                )

            # ---- raw coordinates -------------------------------------------
            m_raw = df_ok[m_dot_col].to_numpy(dtype=float)
            N_raw = df_ok[N_rot_col].to_numpy(dtype=float)

            self._m_scale = (float(m_raw.min()), float(m_raw.max()))
            self._N_scale = (float(N_raw.min()), float(N_raw.max()))

            pts = self._normalize(m_raw, N_raw)

            # ---- build interpolants ----------------------------------------
            self._primary = {}
            self._nn      = {}

            for col in self.targets:
                if col not in df_ok.columns:
                    warnings.warn(
                        f"MapInterpolator: target column {col!r} not found – skipping."
                    )
                    continue
                vals = df_ok[col].to_numpy(dtype=float)

                if method == "linear":
                    self._primary[col] = LinearNDInterpolator(pts, vals)
                elif method == "rbf":
                    self._primary[col] = RBFInterpolator(
                        pts, vals, kernel="thin_plate_spline", smoothing=0.0
                    )
                else:
                    self._primary[col] = NearestNDInterpolator(pts, vals)

                self._nn[col] = NearestNDInterpolator(pts, vals)

            # ---- convex-hull tester (for hull_mask / strict_hull) ----------
            try:
                self._hull = Delaunay(pts)
            except Exception:
                self._hull = None

            # ---- public metadata -------------------------------------------
            self._pts_raw = np.column_stack([m_raw, N_raw])
            self.n_points = len(df_ok)
            self.m_range  = self._m_scale
            self.N_range  = self._N_scale

        # ------------------------------------------------------------------
        # Private helpers
        # ------------------------------------------------------------------

        def _normalize(self, m, N):
            m, N = np.asarray(m, float), np.asarray(N, float)
            if not self.normalize:
                return np.column_stack([m, N])
            m_lo, m_hi = self._m_scale
            N_lo, N_hi = self._N_scale
            m_n = (m - m_lo) / max(m_hi - m_lo, 1e-12)
            N_n = (N - N_lo) / max(N_hi - N_lo, 1e-12)
            return np.column_stack([m_n, N_n])

        # ------------------------------------------------------------------
        # Public query API
        # ------------------------------------------------------------------

        def hull_mask(self, m_arr, N_arr):
            """Boolean array: True where ``(m, N)`` lies inside the convex hull."""
            pts = self._normalize(np.asarray(m_arr, float), np.asarray(N_arr, float))
            if self._hull is None:
                return np.ones(len(pts), dtype=bool)
            return self._hull.find_simplex(pts) >= 0

        def query(self, m_dot: float, N_rot: float) -> dict:
            """
            Interpolate all targets at a single ``(m_dot, N_rot)`` point.

            Returns
            -------
            dict[str, float]
                Keys match ``self.targets``.

            Raises
            ------
            MapInterpolationError
                If point is outside the hull and ``fallback_nn=False``.
            """
            res = self.query_batch(
                np.array([m_dot], dtype=float),
                np.array([N_rot],  dtype=float),
            )
            return {k: float(v[0]) for k, v in res.items()}

        def query_batch(self, m_arr, N_arr) -> dict:
            """
            Vectorised version of :meth:`query`.

            Parameters
            ----------
            m_arr, N_arr : array_like, shape (n,)

            Returns
            -------
            dict[str, np.ndarray]  –  shape (n,) for each target.
            """
            m_arr = np.asarray(m_arr, dtype=float)
            N_arr = np.asarray(N_arr, dtype=float)
            pts   = self._normalize(m_arr, N_arr)

            out = {}
            for col, interp in self._primary.items():
                vals = np.asarray(interp(pts), dtype=float).ravel()

                nan_mask = ~np.isfinite(vals)
                if nan_mask.any():
                    if self.fallback_nn:
                        vals[nan_mask] = self._nn[col](pts[nan_mask]).ravel()
                    else:
                        raise MapInterpolationError(
                            f"{nan_mask.sum()} query point(s) are outside the map "
                            f"convex hull for target {col!r}. "
                            "Set fallback_nn=True or restrict the operating range."
                        )
                out[col] = vals

            return out

        # ------------------------------------------------------------------
        # Persistence
        # ------------------------------------------------------------------

        def save(self, path) -> None:
            """Pickle the interpolator to *path*."""
            with open(path, "wb") as fh:
                pickle.dump(self, fh, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"[MapInterpolator] Saved to {path}  ({self.n_points} pts)")

        @staticmethod
        def load(path):
            """Load a pickled ``MapInterpolator`` from *path*."""
            with open(path, "rb") as fh:
                obj = pickle.load(fh)
            if not hasattr(obj, "_primary"):
                raise TypeError(f"Expected MapInterpolator, got {type(obj)}")
            return obj

        # ------------------------------------------------------------------
        # Diagnostics
        # ------------------------------------------------------------------

        def plot_coverage(self, ax=None, color="steelblue", alpha=0.55):
            """Scatter-plot the support points used to build the interpolant."""
            if ax is None:
                _, ax = plt.subplots(figsize=(6, 4))
            ax.scatter(
                self._pts_raw[:, 0], self._pts_raw[:, 1],
                s=14, c=color, alpha=alpha, label=f"{self.n_points} pts"
            )
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
        """
        Assign stage parameters from arrays.
        If a stage doesn't exist, instantiate it.
        """
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
            P_su = float(self.inputs['P_su']),
            T_su = float(self.inputs['T_su']),
            P_ex = float(self.inputs['P_ex']),
            m_dot = float(self.inputs['m_dot']),
            N_rot = float(self.inputs['N_rot']),
            fluid = self.fluid
        )
        base_params = dict(self.params)

        stage_params = {}
        import numbers
        def _ok(v): return (v is None) or isinstance(v, (numbers.Number, np.floating, np.integer))
        if self.stages:
            keys = set().union(*(vars(st).keys() for st in self.stages))
            blacklist = {'AS','total_states','static_states','Vel_Tri_S','Vel_Tri_R',
                          'eta_is_R','eta_is_S','M1_S','M2_S','M2_R','M3_R',
                          'Y_vec_S','Y_vec_R','delta_S','delta_R','beta_g_S','beta_g_R'}
            for k in sorted(keys - blacklist):
                vals = [getattr(st, k, None) for st in self.stages]
                if any((_ok(v) for v in vals)) and all((_ok(v) for v in vals)):
                    stage_params[k] = vals
        return base_inputs, base_params, stage_params

    # =========================================================================
    #  Map-interpolation interface
    # =========================================================================

    def load_map(
        self,
        df: pd.DataFrame,
        *,
        method: str = "linear",
        targets=None,
        normalize: bool = True,
        fallback_nn: bool = True,
        min_converged_frac: float = 0.2,
        verbose: bool = True,
    ):
        """
        Build a :class:`MapInterpolator` from *df* and attach it to this
        turbine instance.  After calling this method, :meth:`solve_from_map`
        becomes available.

        Parameters
        ----------
        df : pd.DataFrame
            Output of :func:`generate_map_processes` or
            :meth:`generate_map_m_dot_N_rot`.
        method : {'linear', 'rbf', 'nearest'}
            Interpolation method – see :class:`MapInterpolator` for details.
        targets : list[str], optional
            Which columns to interpolate.
            Default: ``['W_dot', 'eta_is', 'P_ex_calc']``.
        normalize : bool
            Normalise ``(m_dot, N_rot)`` to [0,1] before fitting.
        fallback_nn : bool
            Use nearest-neighbour for out-of-hull queries.
        min_converged_frac : float
            Minimum acceptable fraction of converged rows in *df*.
        verbose : bool
            Print a summary line.

        Returns
        -------
        MapInterpolator
            The constructed interpolator (also stored as
            ``self.map_interpolator``).

        Examples
        --------
        >>> df_map = generate_map_processes(Turb_OD, m_grid=..., N_grid=...)
        >>> Turb_OD.load_map(df_map, method='rbf')
        >>> Turb_OD.set_inputs(m_dot=480, ...)
        >>> Turb_OD.solve_from_map()
        >>> print(Turb_OD.W_dot, Turb_OD.eta_is)
        """
        self.map_interpolator = self.MapInterpolator(
            df,
            method=method,
            targets=targets,
            normalize=normalize,
            fallback_nn=fallback_nn,
            min_converged_frac=min_converged_frac,
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

    def load_map_from_file(self, path, verbose: bool = True):
        """
        Restore a previously pickled :class:`MapInterpolator` and attach it.

        Parameters
        ----------
        path : str or Path
            File written by ``turbine.map_interpolator.save(path)``.
        verbose : bool

        Returns
        -------
        MapInterpolator
        """
        self.map_interpolator = self.MapInterpolator.load(path)
        if verbose:
            print(f"[load_map_from_file] Loaded {self.map_interpolator}")
        return self.map_interpolator

    def save_map_df(self, path, df: pd.DataFrame, fmt: str = "parquet") -> None:
        """
        Persist a performance-map DataFrame to disk.

        Parameters
        ----------
        path : str or Path
            Destination file.  The extension is added automatically when
            ``fmt`` is ``'parquet'`` or ``'csv'`` and *path* has none.
        df : pd.DataFrame
            Any DataFrame returned by :func:`generate_map_processes` or
            :meth:`generate_map_m_dot_N_rot`.
        fmt : {'parquet', 'csv'}
            Storage format.

            * ``'parquet'`` *(default)* – compact binary, preserves dtypes
              exactly (bool ``converged`` column stays bool, floats stay
              float64).  Requires ``pyarrow`` or ``fastparquet`` (both are
              common scipy-stack packages).
            * ``'csv'``   – plain text, universally readable, slightly
              larger on disk and reloads all columns as strings until cast.

        Examples
        --------
        >>> df_map = generate_map_processes(Turb_OD, m_grid=..., N_grid=...)
        >>> Turb_OD.save_map_df("turb_map.parquet", df_map)
        >>> # or
        >>> Turb_OD.save_map_df("turb_map.csv", df_map, fmt="csv")
        """
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

    def load_map_df(self, path, fmt: str = "auto") -> pd.DataFrame:
        """
        Reload a performance-map DataFrame from disk and immediately build
        (or rebuild) the :class:`MapInterpolator` attached to this instance.

        This is the counterpart of :meth:`save_map_df`.  After calling it,
        :meth:`solve_from_map` is ready to use without running the expensive
        mean-line sweep again.

        Parameters
        ----------
        path : str or Path
            File previously written by :meth:`save_map_df`.
        fmt : {'auto', 'parquet', 'csv'}
            ``'auto'`` *(default)* infers the format from the file extension.

        Returns
        -------
        pd.DataFrame
            The reloaded map (also passed to :meth:`load_map` internally).

        Examples
        --------
        >>> # First session – generate and save
        >>> df_map = generate_map_processes(Turb_OD, m_grid=..., N_grid=...)
        >>> Turb_OD.save_map_df("turb_map.parquet", df_map)

        >>> # Later session – reload and solve instantly
        >>> df_map = Turb_OD.load_map_df("turb_map.parquet")
        >>> Turb_OD.set_inputs(m_dot=480, P_su=140e5, T_su=394.15,
        ...                    N_rot=2100, fluid='CO2', P_ex=39.8e5)
        >>> Turb_OD.solve_from_map()
        """
        from pathlib import Path as _Path
        p = _Path(path)

        if fmt == "auto":
            fmt = p.suffix.lstrip(".").lower()
            if fmt not in ("parquet", "csv"):
                raise ValueError(
                    f"Cannot infer format from extension {p.suffix!r}. "
                    "Pass fmt='parquet' or fmt='csv' explicitly."
                )

        if fmt == "parquet":
            df = pd.read_parquet(p)
        elif fmt == "csv":
            df = pd.read_csv(p)
            # Restore bool dtype for the 'converged' column if present
            if "converged" in df.columns:
                df["converged"] = df["converged"].astype(bool)
        else:
            raise ValueError(f"fmt must be 'parquet' or 'csv'; got {fmt!r}")

        print(f"[load_map_df] {len(df)} rows loaded from '{p}'")
        self.load_map(df, verbose=True)
        return df

    def solve_from_map(
        self,
        m_dot: float = None,
        N_rot: float = None,
        *,
        strict_hull: bool = False,
    ) -> None:
        """
        Resolve the turbine at the current operating point using the
        pre-computed performance map instead of running the full mean-line
        solver.

        Produces the same output attributes as :meth:`solve`:

        * ``self.W_dot``      – shaft power [W]
        * ``self.eta_is``     – isentropic efficiency [–]
        * ``self.P_ex_calc``  – interpolated outlet static pressure [Pa]

        Also sets ``self._from_map = True`` so downstream code can detect
        which solver was used.

        Parameters
        ----------
        m_dot : float, optional
            Override the mass-flow rate [kg/s].  If omitted,
            ``self.inputs['m_dot']`` is used.
        N_rot : float, optional
            Override the rotational speed [rpm].  If omitted,
            ``self.inputs['N_rot']`` is used.
        strict_hull : bool
            If True, raise :class:`MapInterpolationError` when the query
            point lies outside the map's convex hull.  If False (default),
            the nearest hull point is used transparently.

        Raises
        ------
        MapInterpolationError
            When no map has been loaded yet, or interpolation cannot
            produce a finite result.

        Notes
        -----
        Call :meth:`load_map` (or assign ``self.map_interpolator`` directly)
        before using this method.

        Examples
        --------
        >>> df_map = generate_map_processes(Turb_OD, m_grid=..., N_grid=...)
        >>> Turb_OD.load_map(df_map)          # once
        >>>
        >>> # Fast repeated evaluations (e.g. in a system-level optimiser)
        >>> for m in np.linspace(300, 600, 50):
        ...     Turb_OD.set_inputs(m_dot=m, N_rot=2000, ...)
        ...     Turb_OD.solve_from_map()
        ...     print(m, Turb_OD.W_dot, Turb_OD.eta_is)
        """
        
        if self.map_interpolator is None:
            raise MapInterpolationError(
                "No performance map loaded. "
                "Call  turbine.load_map(df_map)  first."
            )

        _m = float(m_dot if m_dot is not None else self.inputs["m_dot"])
        _N = float(N_rot if N_rot is not None else self.inputs["N_rot"])

        # Optional strict convex-hull check
        if strict_hull:
            inside = self.map_interpolator.hull_mask(
                np.array([_m]), np.array([_N])
            )[0]
            if not inside:
                raise MapInterpolationError(
                    f"Operating point (m_dot={_m:.3f} kg/s, N_rot={_N:.1f} rpm) "
                    "is outside the map convex hull. "
                    "Extend the map range or set strict_hull=False."
                )

        result = self.map_interpolator.query(_m, _N)

        self.W_dot     = result.get("W_dot",     float("nan"))
        self.eta_is    = result.get("eta_is",    float("nan"))
        self.P_ex_calc = result.get("P_ex_calc", float("nan"))

        # Propagate outlet pressure to the exhaust connector if present
        if self.ex is not None:
            try:
                self.ex.p = self.P_ex_calc
            except AttributeError:
                pass

        self._from_map = True

        # Sync overrides back to inputs dict
        if m_dot is not None:
            self.inputs["m_dot"] = _m
        if N_rot is not None:
            self.inputs["N_rot"] = _N

        self.AS.update(CP.PSmass_INPUTS, self.P_ex_calc, self.su.s)
        h_ex_is = self.AS.hmass()
        
        self.h_ex = self.su.h - self.eta_is*(self.su.h - h_ex_is)
        
        self.solved = True
        self.update_connectors(self.P_ex_calc, self.h_ex)
        
        return

    # =========================================================================
    #  Loss Models
    # =========================================================================

    def stator_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]

        h_static_out = x[0]*1e5
        p_static_out = x[1]*1e5
                
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        stage.Vel_Tri_S['u'] = u = self.u
        
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][2]*A_flow)
        
        if self.curr_stage_index == 0:
            stage.Vel_Tri_S['alpha1'] = alpha1 = stage.xhi_S1
            stage.Vel_Tri_S['vu1'] = vu1 = vm*np.tan(alpha1)
        else:
            stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
            stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u 
            stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
        
        
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][1]

        hin = stage.static_states['H'][1]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        h02 = h0in
        
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][2]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)

        stage.Re_s = stage.chord_S*(stage.static_states['D'][2]*vm)/(stage.static_states['V'][2])
        stage.AR_S = stage.h_blade_S/stage.chord_S
                
        stage.beta_g_S = np.arcsin(stage.o_S/stage.pitch_S)
        
        stage.Y_vec_S = aungier_loss_model(stage.Vel_Tri_S['alpha1'], stage.Vel_Tri_S['alpha2'], stage.beta_g_S*180/np.pi, stage.xhi_S1, stage.chord_S, 
                               0, self.params['D_lw'], self.params['e_blade'], stage.h_blade_S, stage.static_states['V'][2], 
                               stage.M1_S, stage.M2_S, self.params['N_lw'], stage.R_c_S, stage.static_states['D'][2], stage.pitch_S, stage.t_blade_S, stage.t_TE_S,
                               vm, v2,1)
        
        self.compute_deviation_stator(stage)
        alpha2_calc = stage.xhi_S2 + stage.delta_S
        
        v2_new = vm/np.cos(alpha2_calc)
        
        Y = stage.Y_vec_S['Y_tot']
        p0_out = (stage.total_states['P'][1] + Y * p_static_out)/(1+Y)

        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = stage.total_states['H'][2]-v2_new**2/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
                
        pout_calc = stage.static_states['P'][2]

        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        return np.array([hout, pout_calc])*1e-5

    def rotor_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]
        
        [h_static_out, p_static_out] = x*1e5
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)
        
        stage.Vel_Tri_R['u'] = u = self.u
        
        A_flow = stage.h_blade_R*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_R['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][3]*A_flow)
        stage.Vel_Tri_R['vu2'] = vu2 = vm*np.tan(stage.Vel_Tri_R['alpha2'])    
        
        stage.Vel_Tri_R['wu2'] = wu2 = vu2 - u
        stage.Vel_Tri_R['w2'] = w2 = np.sqrt(wu2**2 + vm**2)
        stage.M2_R = w2/stage.static_states['A'][2]
        stage.Vel_Tri_R['beta2'] = np.arctan(wu2/vm)
        
        hin = stage.static_states['H'][2]
        h0in = hin + (w2**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        h03 = stage.total_states['H'][2]
        stage.Vel_Tri_R['w3'] = w3 = np.sqrt(2*(h03 - h_static_out))
        stage.M3_R = w3/stage.static_states['A'][3]
        stage.Vel_Tri_R['beta3'] = -np.arccos(vm/w3)
                
        stage.Re_r = stage.chord_R*(stage.static_states['D'][3]*vm)/(stage.static_states['V'][3])
        stage.AR_R = stage.h_blade_R/stage.chord_R
                
        stage.beta_g_R =  np.arcsin(stage.o_R/stage.pitch_R)
                
        stage.Y_vec_R = aungier_loss_model(-stage.Vel_Tri_R['beta2'], -stage.Vel_Tri_R['beta3'], stage.beta_g_R*180/np.pi, -stage.xhi_R1, stage.chord_R, 
                               self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.static_states['V'][3], 
                               stage.M2_R, stage.M3_R, self.params['N_lw'], stage.R_c_R, stage.static_states['D'][3], stage.pitch_R, stage.t_blade_R, stage.t_TE_R,
                               vm, w3,1)

        self.compute_deviation_rotor(stage)
        beta3_calc = stage.xhi_R2 + stage.delta_R

        w3_new = vm/np.cos(beta3_calc)

        Y = stage.Y_vec_R['Y_tot']
                
        p0_out = (stage.total_states['P'][2] + Y * p_static_out)/(1+Y)

        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.total_states['S'][3]
        
        hout = h0in-(w3_new**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 3)
        
        pout_calc = stage.static_states['P'][3]

        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][2])
        hout_s = self.AS.hmass()

        stage.eta_is_R = (stage.static_states['H'][2]-stage.static_states['H'][3])/(stage.static_states['H'][2]-hout_s)
        
        return np.array([hout, pout_calc])*1e-5

    def last_blade_row_system(self, x):
        [h_static_out, p_static_out] = x*1e5
        
        stage = self.stages[-1]
        
        stage.Vel_Tri_S['u'] = u = self.u
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][2]*A_flow)
        
        stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
        stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u 
        stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
        
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][1]

        hin = stage.static_states['H'][1]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)    
                
        h02 = h0in
        
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][2]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)
        
        AR_S = stage.h_blade_S/stage.chord_S
        solidity = (stage.chord_S/stage.pitch_S)
        
        a = 0.0117

        D_e = (np.cos(alpha2)/np.cos(alpha1))*(1.12+a*(alpha1 - stage.xhi_S1)+0.61*np.cos(alpha1)**2 / solidity * (np.tan(alpha1)-np.tan(alpha2)))
        
        P_cst = np.cos(alpha2)/2 * solidity * (v1/v2)**2
        
        Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
        EW_Cst = np.cos((alpha1+alpha2)/2)**3 / np.cos(alpha1)**2

        Yew = 0.02*(solidity/AR_S)/EW_Cst

        DP_loss = (Yp+Yew)*(v1**2)*stage.static_states['D'][1]/2
        p0_out = stage.total_states['P'][1]-DP_loss
                
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(v2**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
        
        pout_calc = stage.static_states['P'][2]

        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        return np.array([hout, pout_calc])*1e-5

    def compute_deviation_stator(self, stage):
        
        delta_0S = np.arcsin((stage.o_S/stage.pitch_S)*(1+(1-stage.o_S/stage.pitch_S)*(2*stage.beta_g_S/np.pi)**2)) - stage.beta_g_S
        
        if stage.M2_S <= 0.5:
            stage.delta_S = delta_0S
        else:
            X = 2*stage.M2_S-1
            stage.delta_S = delta_0S*(1-10*X**3 + 15*X**4 - 6*X**5)
        
        return 

    def compute_deviation_rotor(self, stage):
                
        delta_0R = np.arcsin((stage.o_R/stage.pitch_R)*(1+(1-stage.o_R/stage.pitch_R)*(2*stage.beta_g_R/np.pi)**2)) - abs(stage.beta_g_R)
        
        if stage.M3_R <= 0.5:
            stage.delta_R = delta_0R
        else:
            X = 2*stage.M3_R-1
            stage.delta_R = delta_0R*(1-10*X**3 + 15*X**4 - 6*X**5)

        return 

    # =========================================================================
    #  Flow Computations
    # =========================================================================

    def computeVelTriangle(self):

        self.Vel_Tri['vu2OverU'] = (2*(1-self.R) + self.psi)/2
        self.Vel_Tri['vu3OverU'] = (2*(1-self.R) - self.psi)/2
        self.Vel_Tri['vmOverU']  = self.phi
        
        self.Vel_Tri['wu2OverU']  = self.Vel_Tri['vu2OverU'] - 1
        self.Vel_Tri['wu3OverU']  = self.Vel_Tri['vu3OverU'] - 1

        self.Vel_Tri['v2OverU']  = np.sqrt(self.Vel_Tri['vu2OverU']*self.Vel_Tri['vu2OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['w2OverU']  = np.sqrt(self.Vel_Tri['wu2OverU']*self.Vel_Tri['wu2OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['v3OverU']  = np.sqrt(self.Vel_Tri['vu3OverU']*self.Vel_Tri['vu3OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['w3OverU']  = np.sqrt(self.Vel_Tri['wu3OverU']*self.Vel_Tri['wu3OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])

        self.Vel_Tri['alpha1'] = self.Vel_Tri['alpha3'] = np.arctan(self.Vel_Tri['vu3OverU']/self.Vel_Tri['vmOverU'])
        self.Vel_Tri['alpha2'] = np.arctan(self.Vel_Tri['vu2OverU']/self.Vel_Tri['vmOverU'])

        self.Vel_Tri['beta1'] = self.Vel_Tri['beta3'] = np.arctan(self.Vel_Tri['wu3OverU']/self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta2'] = np.arctan(self.Vel_Tri['wu2OverU']/self.Vel_Tri['vmOverU'])
        
        return 
    
    def computeVelTriangleLastStage(self):

        self.Vel_Tri_Last_Stage['u'] = self.Vel_Tri['u']
        self.Vel_Tri_Last_Stage['vu2'] = 0
        self.Vel_Tri_Last_Stage['vu1'] = self.Vel_Tri['vu3']
        self.Vel_Tri_Last_Stage['vm']  = self.Vel_Tri['vm']
        
        self.Vel_Tri_Last_Stage['wu2'] = self.Vel_Tri_Last_Stage['vu2'] - self.Vel_Tri_Last_Stage['u']
        self.Vel_Tri['v2'] = np.sqrt(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w2'] = np.sqrt(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w3'] = np.sqrt(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)

        self.Vel_Tri_Last_Stage['alpha1'] = self.Vel_Tri['alpha3'] 
        self.Vel_Tri_Last_Stage['alpha2'] = 0

        self.Vel_Tri_Last_Stage['beta1'] = self.Vel_Tri['beta3']
        self.Vel_Tri_Last_Stage['beta2'] = np.arctan(self.Vel_Tri['u']/self.Vel_Tri['vm'])
        
        return 
    
    def computeBladeRow(self, stage_index, row_type):
        stage = self.stages[stage_index]
        
        self.curr_stage_index = stage_index
               
        if row_type == 'S':
            
            if 'P_ex' not in self.inputs:
                RP_1_row = 5**(1/(2*self.nStages)) 
            else:
                RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
            
            if self.Dh0_stage_guess != 0:
                h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2    
            else:
                h_out_guess = stage.static_states['H'][1]*0.99
                
            pout_guess = stage.static_states['P'][1]/RP_1_row

            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
            
            res = 1
            x_in = x0_disc
            c = 0
            
            while res > 1e-8:
                if c > 1000:
                    raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")
                x_out = self.stator_blade_row_system(x_in)
                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                c += 1
                
            self.stator_blade_row_system(x_out)

        else:

            if 'P_ex' not in self.inputs:
                RP_1_row = 5**(1/(2*self.nStages)) 
            else:
                RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
            
            if self.Dh0_stage_guess != 0:
                h_out_guess = stage.static_states['H'][2] - self.Dh0_stage_guess/2    
            else:
                h_out_guess = stage.static_states['H'][2]*0.99
            
            pout_guess = stage.static_states['P'][2]/RP_1_row

            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

            res = 1
            x_in = x0_disc
            c = 0
            
            while res > 1e-8:
                if c > 1000:
                    raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")
                x_out = self.rotor_blade_row_system(x_in) 
                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                c += 1
            
            self.rotor_blade_row_system(x_out)
            self.compute_deviation_rotor(stage)

        return
            
    def computeRepeatingStages(self):
                
        self.nStages = self.params['nStages']
        
        for i in range(int(self.nStages)):
                    
            if i == 0:
                self.computeBladeRow(i, 'S')
                
                self.compute_deviation_stator(self.stages[i])
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']
                
                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

                self.Dh0_stage_guess = self.stages[i].total_states['H'][1] - self.stages[i].total_states['H'][3]

            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
                self.computeBladeRow(i, 'S')
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']

                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

        return
    
    def computeLastStage(self):
        stage = self.stages[-1]
        
        stage.static_states.loc[1] = self.stages[-2].static_states.loc[3]
    
        if 'P_ex' not in self.inputs:
            RP_1_row = 5**(1/(2*self.nStages)) 
        else:
            RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
                
        h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2  
        pout_guess = stage.static_states['P'][1]/RP_1_row

        x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

        res = 1
        x_in = x0_disc
        c = 0
            
        while res > 1e-8:
            if c > 1000:
                raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")
            x_out = self.last_blade_row_system(x_in) 
            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
            c += 1
        
        self.last_blade_row_system(x_out)
        
        return

    # =========================================================================
    #  Map generation
    # =========================================================================
    
    def generate_map_m_dot_N_rot(
        self,
        m_dot_grid=None,
        N_rot_grid=None,
        *,
        m_dot_range=None,
        N_rot_range=None,
        fixed_P_su=None,
        fixed_T_su=None,
        fixed_P_ex=None,
        per_point_hook=None,
        max_retries=2,
        mach_limit=1.2,
        pressure_tol=0.02,
        verbose=False
        ):
        """
        Build an operation map by sweeping mass flow and speed.
    
        Parameters
        ----------
        m_dot_grid : iterable of float
            Mass-flow values [kg/s]. If None, use m_dot_range.
        N_rot_grid : iterable of float
            Rotational speeds [rpm]. If None, use N_rot_range.
        m_dot_range : (mmin, mmax, n)
            Range spec if m_dot_grid is None.
        N_rot_range : (nmin, nmax, n)
            Range spec if N_rot_grid is None.
        fixed_P_su, fixed_T_su, fixed_P_ex : float or None
            If given, overrides the current inputs for the whole map (Pa, K, Pa).
        per_point_hook : callable(self) -> None
            Called right before solve() for each (ṁ, N) point.
        max_retries : int
            Retries on convergence failure.
        mach_limit : float
            Warn if any stage exit Mach exceeds this value.
        pressure_tol : float
            Relative tolerance on outlet pressure.
        verbose : bool
            Print progress.
    
        Returns
        -------
        pandas.DataFrame
        """
        import numpy as _np
        import pandas as _pd
        
        if m_dot_grid is None:
            if m_dot_range is None:
                m0 = self.inputs.get('m_dot', 1.0)
                m_dot_grid = _np.linspace(0.7*m0, 1.3*m0, 9)
            else:
                m_dot_grid = _np.linspace(*m_dot_range)
        else:
            m_dot_grid = _np.array(list(m_dot_grid), dtype=float)
    
        if N_rot_grid is None:
            if N_rot_range is None:
                N0 = self.inputs.get('N_rot', 1000.0)
                N_rot_grid = _np.linspace(0.6*N0, 1.2*N0, 9)
            else:
                N_rot_grid = _np.linspace(*N_rot_range)
        else:
            N_rot_grid = _np.array(list(N_rot_grid), dtype=float)
    
        _P_su0 = self.inputs.get('P_su', None)
        _T_su0 = self.inputs.get('T_su', None)
        _P_ex0 = self.inputs.get('P_ex', None)
    
        rows = []
        total_pts = len(m_dot_grid) * len(N_rot_grid)
        idx = 0
    
        for N in N_rot_grid:
            for m in m_dot_grid:
                idx += 1
                if verbose:
                    print(f"[{idx}/{total_pts}] N={N:.2f} rpm, ṁ={m:.3f} kg/s")
    
                self.set_inputs(
                    m_dot = float(m),
                    P_su  = float(fixed_P_su if fixed_P_su is not None else _P_su0),
                    T_su  = float(fixed_T_su if fixed_T_su is not None else _T_su0),
                    N_rot = float(N),
                    fluid = self.fluid,
                    P_ex  = float(fixed_P_ex if fixed_P_ex is not None else _P_ex0),
                )
    
                hook_note = ""
                if per_point_hook is not None:
                    try:
                        per_point_hook(self)
                    except Exception as e:
                        hook_note = f"per_point_hook failed: {e}"
    
                converged = False
                notes = []
                for attempt in range(max_retries + 1):
                    try:
                        if attempt > 0 and 'damping' in self.params:
                            self.params['damping'] = min(0.8, self.params['damping'] * 1.5)
                        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
                        self.solve()
                        converged = True
                        break
                    except Exception as e:
                        notes.append(f"attempt {attempt}: {e}")
    
                if converged:
                    try:
                        P_su = self.inputs['P_su']
                        T_su = self.inputs['T_su']
                        P_ex_target = self.inputs['P_ex']
                        P_ex_calc = float(self.stages[-1].static_states['P'][2])
                        W_dot = getattr(self, 'W_dot', _np.nan)
                        eta   = getattr(self, 'eta_is', _np.nan)
                        RP_target = P_su / P_ex_target if P_ex_target else _np.nan
                        RP_calc   = P_su / P_ex_calc   if _np.isfinite(P_ex_calc) and P_ex_calc else _np.nan
                        PR = RP_target
                        machs = []
                        for st in self.stages:
                            for label in ('M2_S','M3_R'):
                                if hasattr(st, label):
                                    val = getattr(st, label)
                                    if val is not None:
                                        machs.append(val)
                        mach_warn = bool(len([x for x in machs if _np.isfinite(x) and x > mach_limit]) > 0)
                        pressure_warn = False
                        if P_ex_target and _np.isfinite(P_ex_calc):
                            rel_err = abs(P_ex_calc - P_ex_target)/P_ex_target
                            pressure_warn = rel_err > pressure_tol
                            if pressure_warn:
                                notes.append(f"P_ex mismatch {rel_err:.2%}")
                    except Exception as e:
                        P_su = self.inputs.get('P_su', _np.nan)
                        T_su = self.inputs.get('T_su', _np.nan)
                        P_ex_target = self.inputs.get('P_ex', _np.nan)
                        P_ex_calc = _np.nan
                        PR = _np.nan
                        W_dot = _np.nan
                        eta = _np.nan
                        mach_warn = False
                        pressure_warn = True
                        notes.append(f"post-process error: {e}")
                else:
                    P_su = self.inputs.get('P_su', _np.nan)
                    T_su = self.inputs.get('T_su', _np.nan)
                    P_ex_target = self.inputs.get('P_ex', _np.nan)
                    P_ex_calc = _np.nan
                    PR = _np.nan
                    W_dot = _np.nan
                    eta = _np.nan
                    mach_warn = False
                    pressure_warn = True
    
                rows.append(dict(
                    m_dot=float(m),
                    N_rot=float(N),
                    P_su=float(P_su) if P_su is not None else _np.nan,
                    T_su=float(T_su) if T_su is not None else _np.nan,
                    P_ex_target=float(P_ex_target) if P_ex_target is not None else _np.nan,
                    P_ex_calc=float(P_ex_calc) if _np.isfinite(P_ex_calc) else _np.nan,
                    PR=float(PR) if _np.isfinite(PR) else _np.nan,
                    W_dot=float(W_dot) if _np.isfinite(W_dot) else _np.nan,
                    eta_is=float(eta) if _np.isfinite(eta) else _np.nan,
                    converged=bool(converged),
                    mach_warn=bool(mach_warn),
                    pressure_warn=bool(pressure_warn),
                    notes="; ".join([hook_note] + notes)
                ))
    
        self.set_inputs(
            m_dot=self.inputs['m_dot'],
            P_su=_P_su0,
            T_su=_T_su0,
            N_rot=self.inputs['N_rot'],
            fluid=self.fluid,
            P_ex=_P_ex0
        )
    
        df = _pd.DataFrame(rows)
        df.sort_values(by=['N_rot','m_dot'], inplace=True, ignore_index=True)
        return df

    # =========================================================================
    #  Mean-line solver
    # =========================================================================
    
    def solve(self, map_case = False):
        
        if 'solve_type' in self.params:
            if self.params['solve_type']:
                self.solve_from_map()
                return 
        
        self.omega_rads = 2*np.pi*self.inputs['N_rot']/60
        self.u = self.omega_rads*self.params['r_m']
        
        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
        
        self.computeRepeatingStages()
        
        self.computeLastStage()
        
        hin = self.stages[0].total_states['H'][1]
        hout = self.stages[-1].static_states['H'][2]
        
        self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])

        hout_s = self.AS.hmass()
        
        self.W_dot = self.inputs['m_dot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - hout_s)

        self._from_map = False
        
        self.solved = True
        
        if not map_case:
            self.update_connectors(self.stages[-1].static_states['P'][2], self.stages[-1].static_states['H'][2])
        
        return 
    
    def update_connectors(self, P_ex, h_ex):
        self.ex.reset()
    
        self.ex.set_p(P_ex)
        self.ex.set_h(h_ex)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_fluid(self.su.fluid)
        
        return
    
# =============================================================================
#  Parallel map generation helpers  (module-level, unchanged)
# =============================================================================

def _eval_point_from_snapshot(m, N, base_inputs, base_params, stage_params):
    try:
        turb = AxialTurbineMeanLine(base_inputs['fluid'])
        turb.set_inputs(
            m_dot=float(m),
            P_su=base_inputs['P_su'],
            T_su=base_inputs['T_su'],
            N_rot=float(N),
            fluid=base_inputs['fluid'],
            P_ex=base_inputs['P_ex'],
        )
        turb.set_parameters(**base_params)
        
        if stage_params:
            turb.set_stage_parameters(**stage_params)

        turb.solve(map_case=True)

        P_ex_calc = float(turb.stages[-1].static_states['P'][2])
        RP_calc = turb.inputs['P_su'] / P_ex_calc if P_ex_calc else np.nan
        RP_target = turb.inputs['P_su'] / turb.inputs['P_ex'] if turb.inputs.get('P_ex') else np.nan

        if turb.eta_is < 0.3:
            return dict(
                m_dot=float(m), N_rot=float(N),
                P_su=float(base_inputs.get('P_su', np.nan)),
                T_su=float(base_inputs.get('T_su', np.nan)),
                P_ex_target=float(base_inputs.get('P_ex', np.nan)),
                P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
                W_dot=np.nan, eta_is=np.nan, converged=False
            )
            
        else:
            return dict(
                m_dot=float(m),
                N_rot=float(N),
                P_su=float(turb.inputs['P_su']),
                T_su=float(turb.inputs['T_su']),
                P_ex_target=float(turb.inputs['P_ex']),
                P_ex_calc=P_ex_calc,
                RP_target=RP_target,
                RP_calc=RP_calc,
                W_dot=float(getattr(turb, 'W_dot', np.nan)),
                eta_is=float(getattr(turb, 'eta_is', np.nan)),
                converged=True,
                note=""
            )

    except Exception as e:
        return dict(
            m_dot=float(m), N_rot=float(N),
            P_su=float(base_inputs.get('P_su', np.nan)),
            T_su=float(base_inputs.get('T_su', np.nan)),
            P_ex_target=float(base_inputs.get('P_ex', np.nan)),
            P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
            W_dot=np.nan, eta_is=np.nan, converged=False, note=str(e)
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
    """Context manager to patch joblib to report into tqdm progress bar."""
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
    """
    Multiprocessing map sweep using joblib (loky backend).

    Parameters
    ----------
    machine : AxialTurbineMeanLine
    m_grid : array-like  –  mass-flow values [kg/s]
    N_grid : array-like  –  rotational speeds [rpm]
    max_workers : int    –  -1 = all cores
    desc : str           –  tqdm label

    Returns
    -------
    pd.DataFrame  –  one row per (ṁ, N), sorted by [N_rot, m_dot]
    """
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")

    base_inputs, base_params, stage_params = machine._snapshot_from_machine()

    tasks = [(m, N) for N in N_grid for m in m_grid]
    total = len(tasks)

    with tqdm(total=total, desc=desc, unit="pt",
              dynamic_ncols=True, miniters=1, mininterval=0,
              ascii=True, file=sys.stdout) as bar, tqdm_joblib(bar):
        results = Parallel(
            n_jobs=max_workers,
            backend="loky",
            prefer="processes",
        )(
            delayed(_eval_point_from_snapshot)(m, N, base_inputs, base_params, stage_params)
            for (m, N) in tasks
        )

    return pd.DataFrame(results).sort_values(['N_rot', 'm_dot'], ignore_index=True)


# =============================================================================
#  __main__  –  demo / regression test
# =============================================================================

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
        
        # Turb.set_inputs(
        #     mdot = 318.437021666738, # kg/s
        #     W_dot = 15*1e6, # W : 
        #     p0_su = 15309670.5, # Pa
        #     T0_su = 406.4, # K
        #     p_ex = 5220928, # 5742510, # Pa
        #     )
            
        # turb_params = {'type': 'Axial Turbine',
        #     'mdot_rated': 318.437021666738,
        #     'Wdot_rated': 15244151.612636594,
        #     'N_rot_rated': 2864.775003723627,
        #     'total_to_static_efficiency': 0.8851473369105862,
        #     'DP_rated': 2.93,
        #     'n_stages': 7,
        #     'p0_su': 15309670.5,
        #     'T0_su': 406.4,
        #     'p_ex': 5220928,
        #     'r_m': 0.20956432565203412,
        #     'delta_tip': 0.0004,
        #     'N_lw': 0,
        #     'D_lw': 0,
        #     'e_blade': 2e-06,
        #     'stator': {'h_blade_S': [0.024581131784973807,
        #       0.02731115991390685,
        #       0.030535241526643685,
        #       0.0343428949996286,
        #       0.03884110173212556,
        #       0.04415656370209263,
        #       0.05044235012970171,
        #       0.05396446226956401],
        #      'chord_S': [0.017387241421509592,
        #       0.018189646227485024,
        #       0.019157515021682635,
        #       0.020309513815586707,
        #       0.021666371275761282,
        #       0.023250883041521515,
        #       0.025089330758218745,
        #       0.026110623182803244],
        #      'xhi_S1': [-0.5493774236565352,
        #       -0.5493774236565352,
        #       -0.5493774236565352,
        #       -0.5493774236565352,
        #       -0.5493774236565352,
        #       -0.5493774236565352,
        #       -0.5493774236565352,
        #       -0.5493774236565352],
        #      'xhi_S2': [1.1268112059950612,
        #       1.1268112059950612,
        #       1.1268112059950612,
        #       1.1268112059950612,
        #       1.1268112059950612,
        #       1.1268112059950612,
        #       1.1268112059950612,
        #       1.1268112059950612],
        #      'pitch_S': [0.014118294382457852,
        #       0.014769840363216227,
        #       0.015555741716329778,
        #       0.01649115508675976,
        #       0.01759291197811346,
        #       0.01887952225855718,
        #       0.020372326403956538,
        #       0.021201607297417914],
        #      'o_S': [0.00407007938170583,
        #       0.004257909709554279,
        #       0.00448447255112079,
        #       0.004754137325719167,
        #       0.005071756287732223,
        #       0.0054426655373106475,
        #       0.005873017193716411,
        #       0.006112085665777279],
        #      't_TE_S': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
        #      't_blade_S': [0.0017387241421509593,
        #       0.0018189646227485025,
        #       0.0019157515021682636,
        #       0.002030951381558671,
        #       0.0021666371275761284,
        #       0.0023250883041521517,
        #       0.0025089330758218745,
        #       0.0026110623182803248],
        #      'n_blade_S': [93, 89, 85, 80, 75, 70, 65, 62],
        #      'R_c_S': [0.02593658923532629,
        #       0.02713353838605762,
        #       0.02857731056016104,
        #       0.03029574989120514,
        #       0.03231977738022937,
        #       0.034683397336424615,
        #       0.03742581415244761,
        #       0.03894927848261039]},
        #     'rotor': {'h_blade_R': [0.02599623192810525,
        #       0.028975662485097392,
        #       0.03249353115560146,
        #       0.03664835013019067,
        #       0.041556335411602884,
        #       0.04735747516384924,
        #       0.054218932168004906,
        #       None],
        #      'chord_R': [0.017804442744401783,
        #       0.018691926981905842,
        #       0.019754111574023726,
        #       0.0210108023910845,
        #       0.02248355410347207,
        #       0.02419698360942641,
        #       0.02617885482072475,
        #       None],
        #      'xhi_R1': [0.5033968740474793,
        #       0.5033968740474793,
        #       0.5033968740474793,
        #       0.5033968740474793,
        #       0.5033968740474793,
        #       0.5033968740474793,
        #       0.5033968740474793,
        #       None],
        #      'xhi_R2': [-1.1544956444218015,
        #       -1.1544956444218015,
        #       -1.1544956444218015,
        #       -1.1544956444218015,
        #       -1.1544956444218015,
        #       -1.1544956444218015,
        #       -1.1544956444218015,
        #       None],
        #      'pitch_R': [0.015148494029818542,
        #       0.015903589253319174,
        #       0.01680732419089942,
        #       0.017876550204477332,
        #       0.01912960658162921,
        #       0.02058743803494024,
        #       0.022273666839919944,
        #       None],
        #      'o_R': [0.004776675516045738,
        #       0.005014774752793548,
        #       0.005299743578132708,
        #       0.005636895621769747,
        #       0.006032013691274512,
        #       0.006491702145839364,
        #       0.007023409643056108,
        #       None],
        #      't_TE_R': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
        #      't_blade_R': [0.0019275069515591894,
        #       0.0020395300073626607,
        #       0.002174650407702253,
        #       0.0022904090712424882,
        #       0.0024358345376652743,
        #       0.002590844225248941,
        #       0.002758858045024202,
        #       None],
        #      'n_blade_R': [87, 83, 78, 74, 69, 64, 59, None],
        #      'R_c_R': [0.02655892943743007,
        #       0.02788279177780649,
        #       0.029467255051191132,
        #       0.03134186372129275,
        #       0.03353877094100647,
        #       0.036094697795777626,
        #       0.03905105977860335,
        #       None]}}
        
        turb_params = {'type': 'Axial Turbine',
         'mdot_rated': 318.437021666738,
         'Wdot_rated': 15254770.535339229,
         'N_rot_rated': 2531.510432924823,
         'total_to_static_efficiency': 0.8881461131749432,
         'DP_rated': 2.93,
         'n_stages': 4,
         'p0_su': 15309670.5,
         'T0_su': 406.4,
         'p_ex': 5220928,
         'r_m': 0.321831367853914,
         'delta_tip': 0.0004,
         'N_lw': 0,
         'D_lw': 0,
         'e_blade': 2e-06,
         'stator': {'h_blade_S': [0.01639689254849712,
           0.01991071092020615,
           0.02457275090429655,
           0.03076015473971693,
           0.034272877736767765],
          'chord_S': [0.0041660122936226195,
           0.0045541416421869966,
           0.005073048254277887,
           0.0057468794397966495,
           0.006122580068148005],
          'xhi_S1': [-0.5558768237656709,
           -0.5558768237656709,
           -0.5558768237656709,
           -0.5558768237656709,
           -0.5558768237656709],
          'xhi_S2': [1.2232747686058605,
           1.2232747686058605,
           1.2232747686058605,
           1.2232747686058605,
           1.223275174515224],
          'pitch_S': [0.004391986404224528,
           0.004801168783399159,
           0.005348222086351282,
           0.006058603428737446,
           0.006454682925263108],
          'o_S': [0.0014613537156267446,
           0.0015975017213675359,
           0.001779523773199484,
           0.0020158902640450382,
           0.0021476785235384197],
          't_TE_S': [0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
          't_blade_S': [0.0014581043027679166,
           0.0015939495747654487,
           0.0017755668889972603,
           0.0020114078039288274,
           0.002142903023851802],
          'n_blade_S': [460, 421, 378, 334, 313],
          'R_c_S': [0.00650111641451165,
           0.007106797315350143,
           0.007916557838291198,
           0.008968080174769136,
           0.009554365895926587]},
         'rotor': {'h_blade_R': [0.01799734346049472,
           0.02202340945274219,
           0.027361923783276284,
           0.03445184918650644,
           None],
          'chord_R': [0.004343392521587047,
           0.004791997871406995,
           0.005381998077472141,
           0.006140443385180993,
           None],
          'xhi_R1': [0.6393865201089375,
           0.6393865201089375,
           0.6393865201089375,
           0.6393865201089375,
           None],
          'xhi_R2': [-1.2277087317706656,
           -1.2277087317706656,
           -1.2277087317706656,
           -1.2277087317706656,
           None],
          'pitch_R': [0.004237171385136019,
           0.004674805732487493,
           0.0052503770118362296,
           0.005990273933204133,
           None],
          'o_R': [0.0012136938206303944,
           0.0013390496428988508,
           0.0015039160694796739,
           0.0017158518728125028,
           None],
          't_TE_R': [0.0005, 0.0005, 0.0005, 0.0005, None],
          't_blade_R': [0.0015201873825554664,
           0.0016771992549924483,
           0.0018836993271152493,
           0.002149155184813347,
           None],
          'n_blade_R': [477, 433, 385, 338, None],
          'R_c_R': [0.00677792056926526,
           0.007477975057298509,
           0.008398678059084426,
           0.009582241834689071,
           None]},
         'CAPEX': {'Turbine': 838803.4741547115,
          'Alternator': 478951.6442829673,
          'Installation': 461214.2914531875,
          'Total': 1778969.4098908664}}
        
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
    MAP_SAVE_PATH = r"C:\Users\Basile\Desktop\Travail\Thèse\Travail\WP1\Turbomachines\Save Maps\turb_map.parquet"   # <-- set your path here

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
    
    