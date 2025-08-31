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

def map_plot(
    df, m_col='m_dot', rp_col='RP_calc', eta_col='eta_is', n_col='N_rot',
    levels=20, focus_high=True,
    figsize=(8,6), dpi=200, max_iso_speeds=4,
    # artifact controls for triangulation
    min_circle_ratio=0.02, max_area_factor=6.0, refine_subdiv=0,
    show_points=True,
    # new: render on dense grid using LinearTriInterpolator
    use_grid=False, nx=400, ny=400, smooth_sigma=None):
    
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
            gamma = 3.0   # tweak to densify near max
            base_levels = vmin + (1 - (1 - u)**gamma) * (vmax - vmin)
        else:
            base_levels = np.linspace(vmin, vmax, levels)
    else:
        base_levels = np.asarray(levels, float)
    
    eta_levels = base_levels  # used for filled contours
    
    # ---------------------------------------
    # Contour (line) levels for isolines:
    # keep coarse round ticks AND force extra ones in 0.91–0.95 band
    line_levels = np.unique(
        np.concatenate([
            MaxNLocator(nbins=10).tick_values(vmin, vmax),  # coarse auto levels
            np.arange(0.91, 0.951, 0.01)                    # fine band
        ])
    )

    # --- triangulation + masks (always build; also used for grid interpolation) ---
    tri = mtri.Triangulation(x, y)
    analyzer = mtri.TriAnalyzer(tri)
    try:
        mask_bad = analyzer.circle_ratio_mask(min_circle_ratio=min_circle_ratio)
    except AttributeError:
        # older Matplotlib fallback
        mask_bad = analyzer.get_flat_tri_mask(min_circle_ratio=min_circle_ratio)

    # area-based mask to avoid bridging sparse gaps
    cells = tri.triangles
    xa, ya = x[cells[:,0]], y[cells[:,0]]
    xb, yb = x[cells[:,1]], y[cells[:,1]]
    xc, yc = x[cells[:,2]], y[cells[:,2]]
    area = 0.5 * np.abs(xa*(yb - yc) + xb*(yc - ya) + xc*(ya - yb))
    med_area = np.median(area)
    mask_area = area > (max_area_factor * med_area)

    tri.set_mask(mask_bad | mask_area)

    # optional refinement for tri rendering
    if refine_subdiv and refine_subdiv > 0:
        ref = mtri.UniformTriRefiner(tri)
        tri_f, z_f = ref.refine_field(z, subdiv=int(refine_subdiv))
    else:
        tri_f, z_f = tri, z

    # figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if use_grid:
        # Piecewise-linear interpolation on dense regular grid
        lin = mtri.LinearTriInterpolator(tri, z)

        xi = np.linspace(x.min(), x.max(), nx)
        yi = np.linspace(y.min(), y.max(), ny)
        XI, YI = np.meshgrid(xi, yi)
        ZI = lin(XI, YI)

        # mask outside convex hull (SciPy)
        pts = np.c_[x, y]
        hull = ConvexHull(pts)
        hull_path = Path(pts[hull.vertices])
        mask_out = ~hull_path.contains_points(np.c_[XI.ravel(), YI.ravel()]).reshape(ZI.shape)
        ZI = np.ma.array(ZI, mask=mask_out | ~np.isfinite(ZI))

        # optional gentle smoothing (mask-aware)
        if smooth_sigma is not None:
            from scipy.ndimage import gaussian_filter
            mask = ZI.mask.copy()
            Zfill = ZI.filled(np.nan)
            valid = np.isfinite(Zfill).astype(float)
            Znum = gaussian_filter(np.nan_to_num(Zfill, nan=0.0), sigma=smooth_sigma)
            Zden = gaussian_filter(valid, sigma=smooth_sigma)
            Zsmooth = np.where(Zden > 0, Znum / Zden, np.nan)
            ZI = np.ma.array(Zsmooth, mask=mask | ~np.isfinite(Zsmooth))

        cf = ax.contourf(XI, YI, ZI, levels=eta_levels, cmap="viridis")
        c  = ax.contour (XI, YI, ZI,
                         levels=line_levels,
                         colors="k", linewidths=0.5, alpha=0.45)
    else:
        # plot directly on (cleaned/refined) triangulation
        cf = ax.tricontourf(tri_f, z_f, levels=eta_levels, cmap="viridis")
        c  = ax.tricontour (tri_f, z_f,
                            levels=MaxNLocator(nbins=10).tick_values(vmin, vmax),
                            colors="k", linewidths=0.5, alpha=0.45)

    ax.clabel(c, inline=True, fmt="%.3f", fontsize=7)
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Isentropic efficiency, η")

    # show computed points
    if show_points:
        ax.scatter(x, y, s=10, c="k", alpha=0.35, zorder=3)

    # iso-speed lines (max 4)
    speeds = np.sort(d[n_col].unique())
    if len(speeds) > max_iso_speeds:
        targets = np.linspace(speeds.min(), speeds.max(), max_iso_speeds)
        chosen = sorted(set(min(speeds, key=lambda s: abs(s - t)) for t in targets))
    else:
        chosen = list(speeds)

    for Nval, grp in d.groupby(n_col):
        if Nval not in chosen:
            continue
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

