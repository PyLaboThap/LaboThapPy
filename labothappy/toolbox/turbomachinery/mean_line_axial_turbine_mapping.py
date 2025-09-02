# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 17:35:53 2025

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.ticker import MaxNLocator
from matplotlib.path import Path
from scipy.spatial import ConvexHull

#%%

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.tri import TriAnalyzer, LinearTriInterpolator
from scipy.ndimage import gaussian_filter

def map_plot_clean(
    df,
    *,
    rp_col=None,                 # auto: PR > RP > RP_calc > RP_target
    m_col='m_dot',
    n_col='N_rot',
    z_col='eta_is',
    only_converged=True,
    drop_na=True,
    rp_clip=(1.2, 6.0),          # absolute RP bounds to keep (None to disable)
    qclip_m=(0.01, 0.99),        # quantile clip on m_dot
    qclip_rp=None,               # optional quantile clip on RP within kept data
    levels=24,
    nx=500, ny=500,
    min_circle_ratio=0.01,       # mask skinny triangles
    long_edge_q=0.92,            # mask triangles with an edge longer than this quantile
    smooth_sigma=0.0,            # Gaussian blur (0 => no smoothing)
    cmap='viridis',
    iso_speeds_max=4,            # draw up to N distinct speeds for clarity
    fig=None, ax=None, dpi=150
):
    # --- pick RP column ---
    if rp_col is None:
        for cand in ('PR', 'RP', 'RP_calc', 'RP_target'):
            if cand in df.columns:
                rp_col = cand
                break
        if rp_col is None:
            raise ValueError("No RP/PR column found. Provide rp_col or add one of: PR, RP, RP_calc, RP_target.")

    d = df.copy()

    # --- basic filtering ---
    if only_converged and 'converged' in d.columns:
        d = d[d['converged']]

    if drop_na:
        d = d[np.isfinite(d[m_col]) & np.isfinite(d[rp_col]) & np.isfinite(d[z_col])]

    if rp_clip is not None:
        lo, hi = rp_clip
        d = d[(d[rp_col] >= lo) & (d[rp_col] <= hi)]

    if qclip_m is not None:
        qlo, qhi = d[m_col].quantile(qclip_m[0]), d[m_col].quantile(qclip_m[1])
        d = d[(d[m_col] >= qlo) & (d[m_col] <= qhi)]

    if qclip_rp is not None:
        qlo, qhi = d[rp_col].quantile(qclip_rp[0]), d[rp_col].quantile(qclip_rp[1])
        d = d[(d[rp_col] >= qlo) & (d[rp_col] <= qhi)]

    if len(d) < 10:
        raise ValueError("Not enough points to plot after filtering.")

    x = d[m_col].to_numpy(float)
    y = d[rp_col].to_numpy(float)
    z = d[z_col].to_numpy(float)

    # --- triangulation + cleaning masks ---
    tri = mtri.Triangulation(x, y)

    # 1) remove very flat/skinny triangles
    mask_flat = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio=min_circle_ratio)

    # 2) remove long-edge triangles (they span gaps)
    # compute max edge length per triangle
    p = np.column_stack([x, y])
    tris = tri.triangles
    e01 = np.linalg.norm(p[tris[:, 0]] - p[tris[:, 1]], axis=1)
    e12 = np.linalg.norm(p[tris[:, 1]] - p[tris[:, 2]], axis=1)
    e20 = np.linalg.norm(p[tris[:, 2]] - p[tris[:, 0]], axis=1)
    emax = np.maximum(np.maximum(e01, e12), e20)
    thr = np.quantile(emax, long_edge_q)
    mask_long = emax > thr

    tri.set_mask(mask_flat | mask_long)

    # --- grid only inside the triangulation ---
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    Xi, Yi = np.meshgrid(xi, yi)

    interp = LinearTriInterpolator(tri, z)
    Zi = interp(Xi, Yi)

    # mask grid points outside any triangle
    finder = tri.get_trifinder()
    mask_out = finder(Xi, Yi) == -1
    Zi_masked = np.ma.masked_array(Zi, mask=mask_out | ~np.isfinite(Zi))

    if smooth_sigma and smooth_sigma > 0:
        zi_filled = Zi_masked.filled(np.nan)
        zi_filled = gaussian_filter(np.nan_to_num(zi_filled, nan=0.0), sigma=smooth_sigma)
        Zi_masked = np.ma.masked_array(zi_filled, mask=Zi_masked.mask)

    # --- plot ---
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi)

    cont = ax.contourf(Xi, Yi, Zi_masked, levels=levels, cmap=cmap)
    cbar = fig.colorbar(cont, ax=ax)
    cbar.set_label("Isentropic efficiency, η")

    # scatter the kept points
    ax.scatter(x, y, s=12, c='k', alpha=0.35, lw=0)

    # iso-speed polylines (up to iso_speeds_max distinct speeds)
    if n_col in d.columns and iso_speeds_max:
        speeds = np.unique(d[n_col].round(0))
        # pick 4 evenly spaced speeds over the available set
        pick_idx = np.linspace(0, len(speeds) - 1, num=min(iso_speeds_max, len(speeds)), dtype=int)
        picked = speeds[pick_idx]

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, N in enumerate(picked):
            g = d[np.isfinite(d[n_col]) & (d[n_col].round(0) == N)]
            if len(g) < 5:
                continue
            g = g.sort_values(m_col)
            ax.plot(g[m_col], g[rp_col], lw=2.2, color=colors[i % len(colors)], label=f"{int(N)} rpm")

        if len(picked):
            ax.legend(title="Iso-speeds", ncol=2, frameon=False)

    ax.set_xlabel(f"Mass flow, $\dot m$ [kg/s]")
    ax.set_ylabel(f"Pressure ratio, $R_P$ ({rp_col})")
    ax.set_title("Efficiency map: $\eta(\dot m, RP)$ — cleaned triangulation (no extrapolation)")
    ax.grid(True, alpha=0.25)

    return fig, ax

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.ticker import MaxNLocator

def map_plot(
    df, m_col='m_dot', rp_col='RP_calc', eta_col='eta_is', n_col='N_rot',
    levels=20, focus_high=True,
    figsize=(8,6), dpi=200, max_iso_speeds=4,
    # artifact controls for triangulation
    min_circle_ratio=0.01, max_area_factor=8.0, refine_subdiv=0,
    long_edge_q=None,              # e.g. 0.95 to drop triangles with very long edge
    show_points=True,
    # render on dense grid using LinearTriInterpolator
    use_grid=False, nx=400, ny=400,
    # smoothing and hole filling
    smooth_sigma=None,             # Gaussian sigma on gridded field
    fill_holes=True,               # fill masked interior holes inside hull
    hole_method='nearest',         # 'nearest' only (simple & robust)
    hole_smooth_sigma=0.5          # light blur after filling
):
    d = df.copy()

    # ensure RP_calc exists
    if rp_col not in d.columns and {'P_su','P_ex_calc'}.issubset(d.columns):
        d[rp_col] = d['P_su'] / d['P_ex_calc']

    # keep converged & finite
    if 'converged' in d.columns:
        d = d[d['converged'] == True]
    d = d[[m_col, rp_col, eta_col, n_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if d.empty:
        raise ValueError("No valid data to plot.")

    x = d[m_col].to_numpy(float)
    y = d[rp_col].to_numpy(float)
    z = d[eta_col].to_numpy(float)

    # efficiency levels for the filled contour (color shading)
    vmin, vmax = np.nanmin(z), np.nanmax(z)
    if isinstance(levels, int):
        if focus_high and vmax > vmin:
            u = np.linspace(0, 1, levels)
            gamma = 3.0
            eta_levels = vmin + (1 - (1 - u)**gamma) * (vmax - vmin)
        else:
            eta_levels = np.linspace(vmin, vmax, levels)
    else:
        eta_levels = np.asarray(levels, float)

    # Contour lines: coarse + dense near the top
    line_levels = np.unique(
        np.concatenate([
            MaxNLocator(nbins=10).tick_values(vmin, vmax),
            np.arange(0.91, 0.951, 0.01)
        ])
    )

    # --- triangulation (two copies: clean + full hull) ---
    tri_full = mtri.Triangulation(x, y)  # ONLY for in-hull mask
    tri = mtri.Triangulation(x, y)

    # mask skinny triangles
    analyzer = mtri.TriAnalyzer(tri)
    try:
        mask_flat = analyzer.circle_ratio_mask(min_circle_ratio=min_circle_ratio)
    except AttributeError:
        mask_flat = analyzer.get_flat_tri_mask(min_circle_ratio=min_circle_ratio)

    # mask very large-area triangles (span gaps)
    cells = tri.triangles
    xa, ya = x[cells[:,0]], y[cells[:,0]]
    xb, yb = x[cells[:,1]], y[cells[:,1]]
    xc, yc = x[cells[:,2]], y[cells[:,2]]
    area = 0.5 * np.abs(xa*(yb - yc) + xb*(yc - ya) + xc*(ya - yb))
    med_area = np.median(area)
    mask_area = area > (max_area_factor * med_area)

    # optional: mask triangles with an unusually long edge
    if long_edge_q is not None:
        p = np.column_stack([x, y])
        e01 = np.linalg.norm(p[cells[:,0]] - p[cells[:,1]], axis=1)
        e12 = np.linalg.norm(p[cells[:,1]] - p[cells[:,2]], axis=1)
        e20 = np.linalg.norm(p[cells[:,2]] - p[cells[:,0]], axis=1)
        emax = np.maximum(np.maximum(e01, e12), e20)
        edge_thr = np.quantile(emax, long_edge_q)
        mask_long = emax > edge_thr
    else:
        mask_long = np.zeros(len(cells), dtype=bool)

    tri.set_mask(mask_flat | mask_area | mask_long)

    # optional refinement for tri rendering
    if refine_subdiv and refine_subdiv > 0:
        ref = mtri.UniformTriRefiner(tri)
        tri_f, z_f = ref.refine_field(z, subdiv=int(refine_subdiv))
    else:
        tri_f, z_f = tri, z

    # figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if use_grid:
        # Interpolate on a regular grid
        lin = mtri.LinearTriInterpolator(tri, z)
        xi = np.linspace(x.min(), x.max(), nx)
        yi = np.linspace(y.min(), y.max(), ny)
        XI, YI = np.meshgrid(xi, yi)
        ZI = lin(XI, YI)

        # mask OUTSIDE the full triangulation hull (concave via TriFinder)
        finder_full = tri_full.get_trifinder()
        outside = finder_full(XI, YI) == -1

        # mark interior holes (inside hull but NaN)
        holes = (~outside) & ~np.isfinite(ZI)

        # fill interior holes if requested
        if fill_holes and np.any(holes):
            if hole_method == 'nearest':
                # nearest-neighbor fill from scattered points
                try:
                    from scipy.interpolate import NearestNDInterpolator
                    nn = NearestNDInterpolator(np.c_[x, y], z)
                    ZI[holes] = nn(XI[holes], YI[holes])
                except Exception:
                    # fallback: griddata
                    from scipy.interpolate import griddata
                    ZI[holes] = griddata(np.c_[x, y], z, (XI[holes], YI[holes]), method='nearest')
            # optional light blur to remove patchy look
            if hole_smooth_sigma and hole_smooth_sigma > 0:
                from scipy.ndimage import gaussian_filter
                keep_mask = (~outside) & np.isfinite(ZI)
                Zfill = np.where(keep_mask, ZI, 0.0)
                W = keep_mask.astype(float)
                Znum = gaussian_filter(Zfill, sigma=hole_smooth_sigma)
                Zden = gaussian_filter(W, sigma=hole_smooth_sigma)
                good = Zden > 1e-12
                ZI = np.where(good, Znum / Zden, np.nan)

        # final mask: outside hull remains masked
        ZI = np.ma.array(ZI, mask=outside | ~np.isfinite(ZI))

        # optional overall smoothing (mask-aware)
        if smooth_sigma is not None and smooth_sigma > 0:
            from scipy.ndimage import gaussian_filter
            mask = ZI.mask.copy()
            Zfill = ZI.filled(np.nan)
            valid = np.isfinite(Zfill).astype(float)
            Znum = gaussian_filter(np.nan_to_num(Zfill, nan=0.0), sigma=smooth_sigma)
            Zden = gaussian_filter(valid, sigma=smooth_sigma)
            Zsmooth = np.where(Zden > 0, Znum / Zden, np.nan)
            ZI = np.ma.array(Zsmooth, mask=mask | ~np.isfinite(Zsmooth))

        cf = ax.contourf(XI, YI, ZI, levels=eta_levels, cmap="viridis")
        c  = ax.contour (XI, YI, ZI, levels=line_levels, colors="k", linewidths=0.5, alpha=0.45)

    else:
        # direct plotting on (cleaned/refined) triangulation
        cf = ax.tricontourf(tri_f, z_f, levels=eta_levels, cmap="viridis")
        c  = ax.tricontour (tri_f, z_f,
                            levels=MaxNLocator(nbins=10).tick_values(vmin, vmax),
                            colors="k", linewidths=0.5, alpha=0.45)

    ax.clabel(c, inline=True, fmt="%.3f", fontsize=7)
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Isentropic efficiency, η")

    if show_points:
        ax.scatter(x, y, s=10, c="k", alpha=0.35, zorder=3)

    # iso-speed lines (up to max_iso_speeds)
    speeds = np.sort(d[n_col].unique())
    if len(speeds) > max_iso_speeds:
        targets = np.linspace(speeds.min(), speeds.max(), max_iso_speeds)
        chosen = sorted(set(min(speeds, key=lambda s: abs(s - t)) for t in targets))
    else:
        chosen = list(speeds)
    for Nval, grp in d.groupby(n_col):
        if Nval not in chosen: continue
        g = grp.sort_values(m_col)
        ax.plot(g[m_col].values, g[rp_col].values, '-', lw=1.4, alpha=0.95, label=f"{Nval:.0f} rpm")

    ax.set_xlabel("Mass flow, $\dot m$ [kg/s]")
    ax.set_ylabel("Pressure ratio, $RP_{calc}$")
    ax.set_title("Efficiency map: η($\dot m$, RP) — cleaned triangulation" + (" (gridded)" if use_grid else ""))
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(title="Iso-speeds", loc="lower center", bbox_to_anchor=(0.5, 1.08),
              ncol=min(4, len(chosen)), frameon=False)

    fig.tight_layout()
    return fig, ax



#%%
def pick_speeds(d, n=5):
    """Pick up to n representative N_rot values (evenly across the unique range)."""
    speeds = np.sort(d["N_rot"].unique())
    if len(speeds) <= n:
        return speeds.tolist()
    idx = np.linspace(0, len(speeds)-1, n).round().astype(int)
    return speeds[idx].tolist()

def plot_power_eta_vs_mdot(df_map, speeds=None, max_lines=5, figsize=(7.5, 5), dpi=200):
    # Keep only converged, finite values
    d = df_map.copy()
    if "converged" in d.columns:
        d = d[d["converged"] == True]
    cols_needed = {"m_dot", "N_rot", "W_dot", "eta_is"}
    d = d[list(cols_needed & set(d.columns))].replace([np.inf, -np.inf], np.nan).dropna()
    if d.empty:
        raise ValueError("No valid rows in df_map to plot (check converged/W_dot/eta_is).")

    # Choose which speeds to plot
    chosen = speeds if speeds is not None else pick_speeds(d, n=max_lines)

    # --- Plot W_dot = f(m_dot)
    fig1, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    for N in chosen:
        g = d[d["N_rot"] == N].sort_values("m_dot")
        if g.empty: 
            continue
        ax1.plot(g["m_dot"], g["W_dot"], lw=1.6, label=f"{N:.0f} rpm")
    ax1.set_xlabel(r"Mass flow, $\dot m$ [kg/s]")
    ax1.set_ylabel(r"Shaft power, $\dot W$ [W]")
    ax1.set_title(r"$\dot W$ vs $\dot m$ at selected iso-speeds")
    ax1.grid(True, ls=":", alpha=0.5)
    ax1.legend(title="Iso-speed", ncol=min(4, len(chosen)), frameon=False)
    fig1.tight_layout()

    # --- Plot eta = f(m_dot)
    fig2, ax2 = plt.subplots(figsize=figsize, dpi=dpi)
    for N in chosen:
        g = d[d["N_rot"] == N].sort_values("m_dot")
        if g.empty:
            continue
        ax2.plot(g["m_dot"], g["eta_is"], lw=1.6, label=f"{N:.0f} rpm")
    ax2.set_xlabel(r"Mass flow, $\dot m$ [kg/s]")
    ax2.set_ylabel(r"Isentropic efficiency, $\eta$ [-]")
    ax2.set_title(r"$\eta$ vs $\dot m$ at selected iso-speeds")
    ax2.grid(True, ls=":", alpha=0.5)
    ax2.legend(title="Iso-speed", ncol=min(4, len(chosen)), frameon=False)
    fig2.tight_layout()

    return (fig1, ax1), (fig2, ax2)

# --- Filter sparse points by proximity in (RP, m_dot) -------------------------
def filter_sparse_by_proximity(
    df,
    rp_col=None,           # auto-detect if None: PR > RP > RP_calc > RP_target
    m_col='m_dot',
    group_by='N_rot',      # set to None to consider all speeds together
    rp_tol_rel=0.02,       # neighborhood half-width as fraction of RP range
    m_tol_rel=0.05,        # neighborhood half-width as fraction of m_dot range
    min_neighbors=3,       # require at least this many neighbors (excluding self)
    return_flag=False      # if True, returns (df_filtered, mask_kept)
):
    import numpy as _np
    import pandas as _pd

    # pick the RP/PR column
    if rp_col is None:
        for cand in ('PR', 'RP', 'RP_calc', 'RP_target'):
            if cand in df.columns:
                rp_col = cand
                break
        if rp_col is None:
            raise ValueError("No RP/PR column found. Provide rp_col or add one of: PR, RP, RP_calc, RP_target.")

    # work only on rows with finite values
    finite_mask = _np.isfinite(df[rp_col].to_numpy()) & _np.isfinite(df[m_col].to_numpy())
    work = df.loc[finite_mask].copy()

    # helper to mark non-isolated within one group
    def _mark_group(g):
        vals = g[[rp_col, m_col]].to_numpy(float)
        if len(vals) == 0:
            g['_keep'] = False
            return g
        # scales from group ranges (avoid zero range)
        rp_min, rp_max = _np.nanmin(vals[:,0]), _np.nanmax(vals[:,0])
        m_min,  m_max  = _np.nanmin(vals[:,1]), _np.nanmax(vals[:,1])
        rp_rng = max(rp_max - rp_min, _np.finfo(float).eps)
        m_rng  = max(m_max  - m_min,  _np.finfo(float).eps)
        rp_tol = rp_tol_rel * rp_rng
        m_tol  = m_tol_rel  * m_rng

        # pairwise “rectangular” neighborhood: |ΔRP|<=rp_tol AND |Δm|<=m_tol
        drp = _np.abs(vals[:,None,0] - vals[None,:,0])
        dm  = _np.abs(vals[:,None,1] - vals[None,:,1])
        neigh = (drp <= rp_tol) & (dm <= m_tol)

        counts = neigh.sum(axis=1) - 1  # exclude self
        g['_keep'] = counts >= int(min_neighbors)
        return g

    if group_by is None:
        work = _mark_group(work)
    else:
        work = work.groupby(group_by, as_index=False, group_keys=False).apply(_mark_group)

    keep_idx = work.index[work['_keep'].to_numpy()]
    kept_mask_global = _np.zeros(len(df), dtype=bool)
    kept_mask_global[finite_mask] = df.index[finite_mask].isin(keep_idx)

    df_filtered = df.loc[kept_mask_global].reset_index(drop=True)
    if return_flag:
        return df_filtered, kept_mask_global
    return df_filtered

