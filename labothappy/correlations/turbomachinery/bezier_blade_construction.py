"""
meridional_channel.py
---------------------
Construction et visualisation du canal méridional (hub & tip)
d'un impeller centrifuge dans le plan (x, r).

Méthode Aungier (Centrifugal Compressors, 2000)
───────────────────────────────────────────────
La courbe est une Bézier cubique dont les points de contrôle intérieurs
sont calculés analytiquement pour satisfaire simultanément :
  - tangente imposée en P1  (angle alpha1 → Y1' = tan(alpha1))
  - tangente imposée en P3  (angle alpha3 → Y3' = tan(alpha3))
  - passage exact par un point intermédiaire (X2, Y2)

Les coefficients M1, M2, M3 (dérivées secondes aux nœuds) sont calculés
par les formules d'Aungier, puis convertis en bras de contrôle Bézier.

Segments droits optionnels de longueur L1 (inlet) et L3 (outlet) peuvent
précéder / suivre la Bézier, dans la direction des tangentes imposées.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from math import comb


# ─────────────────────────────────────────────
#  Bézier d'ordre n  (Bernstein général)
# ─────────────────────────────────────────────

def _bernstein(n: int, k: int, t: np.ndarray) -> np.ndarray:
    """B_k^n(t) = C(n,k) * t^k * (1-t)^(n-k)."""
    return comb(n, k) * (t ** k) * ((1 - t) ** (n - k))


def _bezier_order_n(P: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Évalue une Bézier d'ordre n (n+1 points de contrôle P).

    P : (n+1, 2)  t : (m,)  →  retourne (m, 2)
    """
    n   = len(P) - 1
    pts = np.zeros((len(t), 2))
    for k in range(n + 1):
        pts += _bernstein(n, k, t)[:, None] * P[k]
    return pts


# ─────────────────────────────────────────────
#  Formules d'Aungier : M1, M2, M3
# ─────────────────────────────────────────────

def _aungier_bezier_ctrl(
    X1: float, Y1: float, dY1: float,
    X2: float, Y2: float,
    X3: float, Y3: float, dY3: float,
) -> np.ndarray:
    """
    Calcule les 4 points de contrôle d'une Bézier cubique selon Aungier
    (Centrifugal Compressors, 2000).

    La cubique passe par (X1,Y1) et (X3,Y3) avec tangentes imposées
    dY1 = Y1' = tan(alpha1) et dY3 = Y3' = tan(alpha3),
    et passe exactement par le point intermédiaire (X2, Y2).

    Formules Aungier
    ────────────────
    M2 =   6*(Y3-Y2) / [(X3-X2)*(X3-X1)]
         - 6*(Y2-Y1) / [(X2-X1)*(X3-X1)]
         + 2*(dY1-dY3)  / (X3-X1)

    M1 =   3*(Y2-Y1)/(X2-X1)²  -  3*dY1/(X2-X1)  -  M2/2

    M3 = - 3*(Y3-Y2)/(X3-X2)²  +  3*dY3/(X3-X2)  -  M2/2

    Les bras de contrôle intérieurs de la Bézier cubique sont ensuite :
        P1 = P0 + (X2-X1)/3 * [1, dY1 + M1*(X2-X1)/3]  (approx. linéaire)
        P2 = P3 - (X3-X2)/3 * [1, dY3 + M3*(X3-X2)/3]

    Note : M1, M2, M3 sont les dérivées secondes d(²Y)/dX² aux trois
    nœuds — ils servent à positionner P1 et P2 de façon à honorer
    simultanément les deux tangentes ET le point de passage (X2,Y2).

    Paramètres
    ----------
    X1, Y1, dY1 : point de départ et tangente (dY/dX)
    X2, Y2      : point intermédiaire de passage
    X3, Y3, dY3 : point d'arrivée et tangente (dY/dX)

    Retourne
    --------
    ctrl : (4, 2) ndarray  [P0, P1, P2, P3]
    """
    # ── Dérivées secondes aux nœuds (formules Aungier) ────────────────
    M2 = (  6.0 * (Y3 - Y2) / ((X3 - X2) * (X3 - X1))
           - 6.0 * (Y2 - Y1) / ((X2 - X1) * (X3 - X1))
           + 2.0 * (dY1 - dY3) / (X3 - X1) )

    M1 = (  3.0 * (Y2 - Y1) / (X2 - X1) ** 2
           - 3.0 * dY1       / (X2 - X1)
           - M2 / 2.0 )

    M3 = ( -3.0 * (Y3 - Y2) / (X3 - X2) ** 2
           + 3.0 * dY3       / (X3 - X2)
           - M2 / 2.0 )

    # ── Conversion en points de contrôle Bézier cubique ───────────────
    # Pour une Bézier cubique B(t) avec t ∈ [0,1] paramétré par X :
    #   B'(0)  = 3*(P1-P0) / (X3-X1)  →  P1 = P0 + h13/3 * [1, dY1']
    #   B'(1)  = 3*(P3-P2) / (X3-X1)  →  P2 = P3 - h13/3 * [1, dY3']
    # où dY1' et dY3' tiennent compte des corrections de courbure M1/M3.
    h12 = X2 - X1
    h23 = X3 - X2
    h13 = X3 - X1

    # Pentes corrigées aux extrémités (dérivée première incluant M)
    dY1_corr = dY1 + M1 * h12 / 3.0
    dY3_corr = dY3 - M3 * h23 / 3.0

    P0 = np.array([X1, Y1])
    P3 = np.array([X3, Y3])
    P1 = P0 + (h13 / 3.0) * np.array([1.0, dY1_corr])
    P2 = P3 - (h13 / 3.0) * np.array([1.0, dY3_corr])

    return np.array([P0, P1, P2, P3]), (M1, M2, M3)


# ─────────────────────────────────────────────
#  Segments droits optionnels
# ─────────────────────────────────────────────

def _straight_segment(
    origin: np.ndarray,
    angle_deg: float,
    length: float,
    n: int,
    direction: float = 1.0,
) -> np.ndarray:
    """
    Génère n points sur un segment droit depuis `origin`.

    direction = +1 : depuis origin vers l'avant
    direction = -1 : depuis origin vers l'arrière (segment d'outlet)
    """
    dx = direction * length * np.cos(np.radians(angle_deg))
    dr = direction * length * np.sin(np.radians(angle_deg))
    s  = np.linspace(0.0, 1.0, n)
    xs = origin[0] + s * dx
    rs = origin[1] + s * dr
    return np.column_stack([xs, rs])


# ─────────────────────────────────────────────
#  Clipping à r = r_target
# ─────────────────────────────────────────────

def _clip_contour_at_r(pts: np.ndarray, r_target: float) -> tuple[np.ndarray, float]:
    """
    Clip un contour (m, 2) en (x, r) à r = r_target.
    Suppose r monotone croissant.
    Retourne (contour clippé, x à l'intersection).
    """
    r   = pts[:, 1]
    x   = pts[:, 0]
    idx = np.clip(np.searchsorted(r, r_target), 1, len(r) - 1)
    r0, r1 = r[idx - 1], r[idx]
    x0, x1 = x[idx - 1], x[idx]
    frac    = (r_target - r0) / (r1 - r0) if (r1 - r0) > 1e-14 else 0.0
    x_cross = float(x0 + frac * (x1 - x0))
    return np.vstack([pts[:idx], [[x_cross, r_target]]]), x_cross


# ─────────────────────────────────────────────
#  Construction du canal méridional
# ─────────────────────────────────────────────

def build_meridional_channel(
    # Points de passage
    X1: float, Y1: float,      # point de départ  (inlet de la Bézier)
    X2: float, Y2: float,      # point intermédiaire imposé
    X3: float, Y3: float,      # point d'arrivée  (outlet de la Bézier)
    # Tangentes aux extrémités [degrés]
    alpha1_deg: float,         # angle à l'inlet  → dY1 = tan(alpha1)
    alpha3_deg: float,         # angle à l'outlet → dY3 = tan(alpha3)
    # Segments droits optionnels
    L1: float = 0.0,           # longueur du segment droit à l'inlet  (0 = absent)
    L3: float = 0.0,           # longueur du segment droit à l'outlet (0 = absent)
    # Clipping outlet
    r_out: float | None = None,  # si fourni, clip le contour à r = r_out
    # Discrétisation
    n_points: int = 600,
) -> dict:
    """
    Construit un contour méridional par la méthode d'Aungier :
    segment droit L1 + Bézier cubique (passant par X2,Y2) + segment droit L3.

    Les tangentes sont imposées aux deux extrémités de la Bézier.
    Les coefficients M1, M2, M3 sont calculés analytiquement pour garantir
    le passage exact par (X2, Y2).

    Paramètres
    ----------
    X1, Y1       : point de départ de la Bézier  (= fin du segment L1)
    X2, Y2       : point intermédiaire de passage
    X3, Y3       : point d'arrivée de la Bézier  (= début du segment L3)
    alpha1_deg   : angle de tangente à l'inlet [°]   (dY/dX = tan(alpha1))
    alpha3_deg   : angle de tangente à l'outlet [°]
    L1           : longueur du segment droit inlet  [même unité que X, Y]
    L3           : longueur du segment droit outlet
    r_out        : rayon de clipping (optionnel)
    n_points     : nombre de points de discrétisation

    Retourne
    --------
    dict :
        "contour"    – (m, 2)   contour complet (L1 + Bézier + L3, clippé si r_out)
        "bezier"     – (m, 2)   partie Bézier seule
        "ctrl"       – (4, 2)   points de contrôle Bézier
        "M1","M2","M3"          dérivées secondes aux nœuds (scalaires)
        "seg_inlet"  – (m, 2) | None   segment droit inlet
        "seg_outlet" – (m, 2) | None   segment droit outlet
        "x_out"      – float | None    x à l'intersection r_out
    """
    dY1 = np.tan(np.radians(alpha1_deg))
    dY3 = np.tan(np.radians(alpha3_deg))

    # ── Points de contrôle Bézier ──────────────────────────────────────
    ctrl, (M1, M2, M3) = _aungier_bezier_ctrl(
        X1, Y1, dY1,
        X2, Y2,
        X3, Y3, dY3,
    )

    # ── Courbe Bézier ──────────────────────────────────────────────────
    t_bez  = np.linspace(0.0, 1.0, n_points)
    bezier = _bezier_order_n(ctrl, t_bez)

    # ── Segments droits optionnels ─────────────────────────────────────
    n_seg  = max(10, n_points // 20)

    seg_inlet = None
    if L1 > 0.0:
        end_inlet = np.array([X1, Y1])
        seg_inlet = _straight_segment(end_inlet, alpha1_deg, -L1, n_seg, direction=-1.0)
        seg_inlet = seg_inlet[::-1]   # ordre : loin → P1

    seg_outlet = None
    if L3 > 0.0:
        start_outlet = np.array([X3, Y3])
        seg_outlet   = _straight_segment(start_outlet, alpha3_deg, L3, n_seg)

    # ── Assemblage ────────────────────────────────────────────────────
    parts = []
    if seg_inlet  is not None: parts.append(seg_inlet)
    parts.append(bezier)
    if seg_outlet is not None: parts.append(seg_outlet)
    contour = np.vstack(parts)

    # ── Clipping ──────────────────────────────────────────────────────
    x_out = None
    if r_out is not None:
        contour, x_out = _clip_contour_at_r(contour, r_out)

    return {
        "contour":    contour,
        "bezier":     bezier,
        "ctrl":       ctrl,
        "M1": M1, "M2": M2, "M3": M3,
        "seg_inlet":  seg_inlet,
        "seg_outlet": seg_outlet,
        "x_out":      x_out,
        "X1": X1, "Y1": Y1,
        "X2": X2, "Y2": Y2,
        "X3": X3, "Y3": Y3,
        "alpha1_deg": alpha1_deg,
        "alpha3_deg": alpha3_deg,
    }


# ─────────────────────────────────────────────
#  Visualisation
# ─────────────────────────────────────────────

def plot_meridional(
    hub:        dict,
    tip:        dict,
    show_ctrl:  bool = True,
    show_annot: bool = True,
    title:      str  = "Vue méridionale — Méthode Aungier",
    ax:         plt.Axes | None = None,
) -> plt.Figure:
    """
    Trace les contours hub et tip dans le plan méridional (r horizontal,
    x axial, axe y inversé — inlet en haut).

    hub, tip : dicts retournés par build_meridional_channel()
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    else:
        fig = ax.get_figure()

    def _plot_contour(geo, color, label):
        cnt  = geo["contour"]
        ctrl = geo["ctrl"]
        ax.plot(cnt[:, 1], cnt[:, 0], color=color, lw=2.5,
                label=label, zorder=4)
        # point intermédiaire imposé
        ax.plot(geo["Y2"], geo["X2"], "D", color=color,
                ms=6, zorder=6, label=f"{label} — point X2,Y2")
        if show_ctrl:
            ax.plot(ctrl[:, 1], ctrl[:, 0], "o--", color=color,
                    lw=0.9, ms=5, alpha=0.5, zorder=5)

    _plot_contour(hub, "#8B4513", "Hub")
    _plot_contour(tip, "#2c7bb6", "Tip (shroud)")

    # Passage fluide
    h_cnt = hub["contour"]
    t_cnt = tip["contour"]
    pass_r = np.concatenate([h_cnt[:, 1], t_cnt[::-1, 1]])
    pass_x = np.concatenate([h_cnt[:, 0], t_cnt[::-1, 0]])
    ax.fill(pass_r, pass_x, color="#dceefb", alpha=0.5, zorder=1)
    ax.fill(pass_r, pass_x, hatch="///", facecolor="none",
            edgecolor="#aac8e0", lw=0.4, alpha=0.4, zorder=2)

    # Inlet
    ax.plot([hub["Y1"], tip["Y1"]], [hub["X1"], tip["X1"]],
            color="#e67e22", lw=2.0, ls="--", label="Inlet", zorder=3)
    # Outlet
    ax.plot([hub["contour"][-1, 1], tip["contour"][-1, 1]],
            [hub["contour"][-1, 0], tip["contour"][-1, 0]],
            color="#27ae60", lw=2.0, ls="--", label="Outlet", zorder=3)

    # Axe de rotation
    ax.axvline(0, color="#ccc", lw=0.8, ls=":", zorder=0)

    if show_annot:
        for geo, color in [(hub, "#8B4513"), (tip, "#2c7bb6")]:
            ax.annotate(
                f"M1={geo['M1']:.2f}\nM2={geo['M2']:.2f}\nM3={geo['M3']:.2f}",
                xy=(geo["Y2"], geo["X2"]),
                xytext=(geo["Y2"] + 0.002, geo["X2"]),
                fontsize=7, color=color,
                arrowprops=dict(arrowstyle="->", color=color, lw=0.8),
            )

    ax.set_xlabel("r  — radial  [m]")
    ax.set_ylabel("x  — axial  [m]")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(left=0)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.8)
    ax.grid(True, ls=":", alpha=0.35)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────
#  Exemple d'utilisation
# ─────────────────────────────────────────────

if __name__ == "__main__":

    hub = build_meridional_channel(
        X1=0.000, Y1=0.003,
        X2=0.005, Y2=0.012,
        X3=0.011, Y3=0.0293,
        alpha1_deg=5.0,
        alpha3_deg=88.0,
        L1=0.001, L3=0.01,
        r_out=0.0293,
        n_points=600,
    )

    tip = build_meridional_channel(
        X1=0.000, Y1=0.012,
        X2=0.004, Y2=0.020,
        X3=0.009, Y3=0.0293,
        alpha1_deg=3.0,
        alpha3_deg=85.0,
        L1=0.0, L3=0.0,
        r_out=0.0293,
        n_points=600,
    )

    print(f"Hub  — M1={hub['M1']:.3f}  M2={hub['M2']:.3f}  M3={hub['M3']:.3f}")
    print(f"Tip  — M1={tip['M1']:.3f}  M2={tip['M2']:.3f}  M3={tip['M3']:.3f}")

    fig = plot_meridional(hub, tip, show_ctrl=True, show_annot=True)

    plt.show()