"""
casey_robinson_meridional.py
============================

Méthode de Casey & Robinson (ASME GT2008-50561) — section "Meridional flow channel".
Construit la géométrie méridionale d'un impeller centrifuge (plan z–r) en trois
patches Bézier raccordés C¹ :

       Inlet duct  ─►  Impeller  ─►  Vaneless diffuser
                   ◀ LE ▶          ◀ TE ▶

Méthodologie
────────────
1. Chaque patch = une Bézier hub + une Bézier shroud.
   Les coins (LE/TE) sont des dicts partagés entre patches → C⁰ garanti.

2. Hub de l'impeller   : degré 4 (3 points internes)
   Shroud de l'impeller: degré 3 (2 points internes), élevé à 4 pour matcher le hub.
   Inlet & diffuser    : degré 3.

3. Points internes placés le long des tangentes imposées aux coins,
   à une fraction "shape factor" ∈ [0,1] de la longueur LE→TE.

4. Continuité C¹ aux raccords : même tangente au TE du patch i et au LE du patch i+1.

5. Élévation de degré exacte :
       P'_k = (k/(n+1)) * P_{k-1} + (1 - k/(n+1)) * P_k
"""

from __future__ import annotations
from math import comb
import numpy as np
import matplotlib.pyplot as plt


# ════════════════════════════════════════════════════════════════════════════
#  BÉZIER : primitives
# ════════════════════════════════════════════════════════════════════════════

def bezier_eval(P: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Input(s)
    ----------
    - P : coordinates of the points to interpolate
    - t is the curve parcour position vector (between 0 and 1 and translates the 
    evolution along the curve) - its number of elements represents the discretization level

    Objective
    ------------
    bezier interpolation - Interpolates data using:
    
    B(t) = Σ C(n,k) * t^k * (1-t)^(n-k) * P_k

    -> n is the interpolation degree, equal to the number of provided points to interp - 1
    -> k is an iterator representing each of the provided points

    Output(s)
    ------------
    returns out : Contains the coordinates of interpolated values (= number of disc in t)
    Interpolation being done based on the P values
    """
    n   = len(P) - 1
    out = np.zeros((len(t), 2))
    for k in range(n + 1):
        w    = comb(n, k) * (t ** k) * ((1 - t) ** (n - k))
        out += w[:, None] * P[k]
    return out


def elevate_degree(P: np.ndarray) -> np.ndarray:
    """
    Elevates the degree of the bezier polynomial without modifying the curve.
    Identity : P'_k = (k/(n+1)) * P_{k-1} + (1 - k/(n+1)) * P_k
    """
    n      = len(P) - 1
    new    = np.zeros((n + 2, 2))
    new[0] = P[0]
    new[-1]= P[-1]
    for k in range(1, n + 1):
        a      = k / (n + 1)
        new[k] = a * P[k - 1] + (1 - a) * P[k]
    return new


# ════════════════════════════════════════════════════════════════════════════
#  CONSTRUCTEURS DE DICTS
# ════════════════════════════════════════════════════════════════════════════

def make_corner(z: float, r: float) -> dict:
    """
    Crée un coin de patch.

    Retourne
    --------
    {"z": float, "r": float}
    """
    return {"z": z, "r": r}


def corner_xy(c: dict) -> np.ndarray:
    """Coordonnées (z, r) d'un coin sous forme (2,) ndarray."""
    return np.array([c["z"], c["r"]])


def make_patch(
    name: str,
    LE_hub:     dict, LE_shroud:  dict,
    TE_hub:     dict, TE_shroud:  dict,
    tan_LE_hub:    float = 0.0,
    tan_LE_shroud: float = 0.0,
    tan_TE_hub:    float = np.pi / 2,
    tan_TE_shroud: float = np.pi / 2,
    sf_hub:    tuple = (0.5, 0.5, 0.5),
    sf_shroud: tuple = (0.5, 0.5),
) -> dict:
    """
    Crée un patch méridional.

    Paramètres
    ----------
    name          : identifiant ("inlet", "impeller", "diffuser")
    LE_hub/shroud : coins d'entrée (dicts make_corner)
    TE_hub/shroud : coins de sortie
    tan_LE/TE_*   : angles de tangente [rad depuis +z]
                    0 = axial pur,  π/2 = radial pur
    sf_hub        : 3 shape factors hub   (points P1, P2, P3 internes)
    sf_shroud     : 2 shape factors shroud (points P1, P2 internes)

    Retourne
    --------
    dict avec toutes les clés ci-dessus + "ctrl_hub" et "ctrl_shroud"
    initialisés à None (remplis par build_channel).
    """
    return {
        "name":           name,
        "LE_hub":         LE_hub,
        "LE_shroud":      LE_shroud,
        "TE_hub":         TE_hub,
        "TE_shroud":      TE_shroud,
        "tan_LE_hub":     tan_LE_hub,
        "tan_LE_shroud":  tan_LE_shroud,
        "tan_TE_hub":     tan_TE_hub,
        "tan_TE_shroud":  tan_TE_shroud,
        "sf_hub":         sf_hub,
        "sf_shroud":      sf_shroud,
        "ctrl_hub":       None,
        "ctrl_shroud":    None,
    }


# ════════════════════════════════════════════════════════════════════════════
#  CONSTRUCTION DES POINTS DE CONTRÔLE
# ════════════════════════════════════════════════════════════════════════════

def _tangent_vec(angle: float) -> np.ndarray:
    return np.array([np.cos(angle), np.sin(angle)])


def _patch_length(LE: dict, TE: dict) -> float:
    return float(np.hypot(TE["z"] - LE["z"], TE["r"] - LE["r"]))


def build_hub_ctrl(p: dict) -> np.ndarray:
    """
    Hub : Bézier degré 4 — 5 points de contrôle, 3 internes.

    P0 = LE_hub
    P1 = P0 + sf[0] * L * tan_LE          (bras LE)
    P2 = P1 + sf[1] * (P3 - P1)           (point milieu)
    P3 = P4 - sf[2] * L * tan_TE          (bras TE)
    P4 = TE_hub
    """
    L  = _patch_length(p["LE_hub"], p["TE_hub"])
    t0 = _tangent_vec(p["tan_LE_hub"])
    t1 = _tangent_vec(p["tan_TE_hub"])
    P0 = corner_xy(p["LE_hub"])
    P4 = corner_xy(p["TE_hub"])
    P1 = P0 + p["sf_hub"][0] * L * t0
    P3 = P4 - p["sf_hub"][2] * L * t1
    P2 = P1 + p["sf_hub"][1] * (P3 - P1)
    return np.array([P0, P1, P2, P3, P4])


def build_shroud_ctrl(p: dict) -> np.ndarray:
    """
    Shroud : Bézier degré 3 élevée à degré 4.

    P0 = LE_shroud
    P1 = P0 + sf[0] * L * tan_LE
    P2 = P3 - sf[1] * L * tan_TE
    P3 = TE_shroud
    → élévation au degré 4 (identité exacte)
    """
    L  = _patch_length(p["LE_shroud"], p["TE_shroud"])
    t0 = _tangent_vec(p["tan_LE_shroud"])
    t1 = _tangent_vec(p["tan_TE_shroud"])
    P0 = corner_xy(p["LE_shroud"])
    P3 = corner_xy(p["TE_shroud"])
    P1 = P0 + p["sf_shroud"][0] * L * t0
    P2 = P3 - p["sf_shroud"][1] * L * t1
    return elevate_degree(np.array([P0, P1, P2, P3]))


def discretise_patch(p: dict, n_points: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """Construit et discrétise hub & shroud. Retourne (hub_pts, shroud_pts) (n,2)."""
    p["ctrl_hub"]    = build_hub_ctrl(p)
    p["ctrl_shroud"] = build_shroud_ctrl(p)
    t                = np.linspace(0.0, 1.0, n_points)
    return bezier_eval(p["ctrl_hub"], t), bezier_eval(p["ctrl_shroud"], t)


# ════════════════════════════════════════════════════════════════════════════
#  ASSEMBLAGE MULTI-PATCH
# ════════════════════════════════════════════════════════════════════════════

def enforce_C1_between(p_prev: dict, p_next: dict) -> None:
    """Impose C¹ au raccord : tangentes TE de p_prev → LE de p_next."""
    p_next["tan_LE_hub"]    = p_prev["tan_TE_hub"]
    p_next["tan_LE_shroud"] = p_prev["tan_TE_shroud"]


def build_channel(patches: list[dict], n_per_patch: int = 200) -> dict:
    """
    Assemble les patches en un canal méridional complet.

    Retourne
    --------
    {
      "hub_full"    : (N, 2)  hub complet,
      "shroud_full" : (N, 2)  shroud complet,
      "patches"     : liste de dicts par patch (hub, shroud, ctrl_hub, ctrl_shroud, name)
    }
    """
    for i in range(len(patches) - 1):
        enforce_C1_between(patches[i], patches[i + 1])

    hub_chunks, shroud_chunks, patch_data = [], [], []

    for k, p in enumerate(patches):
        hub_pts, shroud_pts = discretise_patch(p, n_per_patch)
        if k > 0:                      # évite de doubler le point au raccord
            hub_pts    = hub_pts[1:]
            shroud_pts = shroud_pts[1:]
        hub_chunks.append(hub_pts)
        shroud_chunks.append(shroud_pts)
        patch_data.append({
            "name":        p["name"],
            "hub":         hub_pts,
            "shroud":      shroud_pts,
            "ctrl_hub":    p["ctrl_hub"],
            "ctrl_shroud": p["ctrl_shroud"],
        })

    return {
        "hub_full":    np.vstack(hub_chunks),
        "shroud_full": np.vstack(shroud_chunks),
        "patches":     patch_data,
    }


# ════════════════════════════════════════════════════════════════════════════
#  VISUALISATION
# ════════════════════════════════════════════════════════════════════════════

PATCH_COLORS = {
    "inlet":    ("#d3d1c7", "Inlet duct"),
    "impeller": ("#85b7eb", "Impeller"),
    "diffuser": ("#9fe1cb", "Vaneless diffuser"),
}


def annotate_stations(ax: plt.Axes, channel: dict) -> None:
    patches = channel["patches"]
    stations = []
    if len(patches) >= 1:
        stations.append(("0", patches[0]["hub"][0],  patches[0]["shroud"][0]))
    if len(patches) >= 2:
        stations.append(("1", patches[1]["hub"][0],  patches[1]["shroud"][0]))
        stations.append(("2", patches[1]["hub"][-1], patches[1]["shroud"][-1]))
    if len(patches) >= 3:
        stations.append(("3", patches[2]["hub"][-1], patches[2]["shroud"][-1]))
    for label, hub_pt, sh_pt in stations:
        z_mid = 0.5 * (hub_pt[0] + sh_pt[0])
        r_top = max(hub_pt[1], sh_pt[1]) + 0.012
        ax.text(z_mid, r_top, label, ha="center", va="bottom",
                fontsize=10, fontweight="bold", color="#c0392b",
                bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec="#c0392b", lw=1.0))


def plot_channel(channel: dict, title: str = "Meridional channel — Casey/Robinson",
                 show_ctrl: bool = True, ax: plt.Axes | None = None) -> plt.Figure:
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 6))
    else:
        fig = ax.get_figure()

    for pd in channel["patches"]:
        color, _ = PATCH_COLORS.get(pd["name"].lower(), ("#eeeeee", pd["name"]))
        poly_z = np.concatenate([pd["hub"][:, 0], pd["shroud"][::-1, 0]])
        poly_r = np.concatenate([pd["hub"][:, 1], pd["shroud"][::-1, 1]])
        ax.fill(poly_z, poly_r, color=color, alpha=0.45, zorder=1, edgecolor="none")

    ax.plot(channel["hub_full"][:, 0],    channel["hub_full"][:, 1],
            "-", color="#444441", lw=2.5, label="Hub",    zorder=4)
    ax.plot(channel["shroud_full"][:, 0], channel["shroud_full"][:, 1],
            "-", color="#185fa5", lw=2.5, label="Shroud", zorder=4)

    if show_ctrl:
        for pd in channel["patches"]:
            for ctrl, color, marker in [
                (pd["ctrl_hub"],    "#444441", "o"),
                (pd["ctrl_shroud"], "#185fa5", "s"),
            ]:
                ax.plot(ctrl[:, 0], ctrl[:, 1], "--", color=color,
                        lw=0.6, alpha=0.4, zorder=5)
                ax.plot(ctrl[1:-1, 0], ctrl[1:-1, 1], marker, color=color,
                        ms=5, alpha=0.7, zorder=5,
                        markeredgecolor="white", markeredgewidth=0.8)

    seen = set()
    for pd in channel["patches"]:
        for corner in [tuple(pd["hub"][0]),    tuple(pd["hub"][-1]),
                       tuple(pd["shroud"][0]), tuple(pd["shroud"][-1])]:
            if corner not in seen:
                seen.add(corner)
                ax.plot(corner[0], corner[1], "k+", ms=10, mew=1.5, zorder=6)

    ax.axhline(0, color="#aaa", lw=0.6, ls=":", zorder=0)
    annotate_stations(ax, channel)

    from matplotlib.patches import Patch as MplPatch
    legend_patches = [
        MplPatch(facecolor=color, alpha=0.45, edgecolor="none", label=label)
        for _, (color, label) in PATCH_COLORS.items()
    ]
    leg1 = ax.legend(handles=legend_patches, loc="upper left",
                     fontsize=8, framealpha=0.9, title="Patches")
    ax.add_artist(leg1)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    ax.set_xlabel("z  — axial  [m]")
    ax.set_ylabel("r  — radial [m]")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, ls=":", alpha=0.35)
    fig.tight_layout()
    return fig


# ════════════════════════════════════════════════════════════════════════════
#  CAS TEST : ECKARDT ROTOR O (1980)
# ════════════════════════════════════════════════════════════════════════════

def eckardt_rotor_O() -> list[dict]:
    """
    3 patches du Eckardt Rotor O  (D2=400 mm, r1h=45 mm, r1s=140 mm, b2=26 mm).
    z = 0 au LE impeller.
    
    Storage of the coordinates of the different geometry elements
    """
    LE_inlet_hub    = make_corner(z=-0.030, r=0.045)
    LE_inlet_shroud = make_corner(z=-0.030, r=0.140)
    LE_imp_hub      = make_corner(z=0.000,  r=0.045)
    LE_imp_shroud   = make_corner(z=0.000,  r=0.140)
    TE_imp_hub      = make_corner(z=0.130,          r=0.200)
    TE_imp_shroud   = make_corner(z=0.130 - 0.026,  r=0.200)
    r_diff_out      = 0.200 * 1.69
    TE_diff_hub     = make_corner(z=0.130,          r=r_diff_out)
    TE_diff_shroud  = make_corner(z=0.130 - 0.026,  r=r_diff_out)

    inlet = make_patch(
        name="inlet",
        LE_hub=LE_inlet_hub,    LE_shroud=LE_inlet_shroud,
        TE_hub=LE_imp_hub,      TE_shroud=LE_imp_shroud,
        tan_LE_hub=0.0,  tan_LE_shroud=0.0,
        tan_TE_hub=0.0,  tan_TE_shroud=0.0,
        sf_hub=(0.4, 0.5, 0.4), sf_shroud=(0.4, 0.4),
    )
    impeller = make_patch(
        name="impeller",
        LE_hub=LE_imp_hub,      LE_shroud=LE_imp_shroud,
        TE_hub=TE_imp_hub,      TE_shroud=TE_imp_shroud,
        tan_LE_hub=0.0,         tan_LE_shroud=0.0,
        tan_TE_hub=np.pi / 2,   tan_TE_shroud=np.pi / 2,
        sf_hub=(0.55, 0.50, 0.55), sf_shroud=(0.45, 0.45),
    )
    diffuser = make_patch(
        name="diffuser",
        LE_hub=TE_imp_hub,      LE_shroud=TE_imp_shroud,
        TE_hub=TE_diff_hub,     TE_shroud=TE_diff_shroud,
        tan_LE_hub=np.pi / 2,   tan_LE_shroud=np.pi / 2,
        tan_TE_hub=np.pi / 2,   tan_TE_shroud=np.pi / 2,
        sf_hub=(0.3, 0.5, 0.3), sf_shroud=(0.3, 0.3),
    )
    return [inlet, impeller, diffuser]

# ════════════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    patches = eckardt_rotor_O()
    channel = build_channel(patches, n_per_patch=200)

    fig = plot_channel(channel,
                        title="Eckardt Rotor O — méthode Casey/Robinson (3 patches Bézier)",
                        show_ctrl=True)
    
    fig2, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    for ax, sf_le in zip(axes, [0.30, 0.55, 0.80]):
        ps = eckardt_rotor_O()
        ps[1]["sf_hub"] = (sf_le, 0.50, 0.55)
        ch = build_channel(ps, n_per_patch=200)
        plot_channel(ch, title=f"sf_hub[0] = {sf_le:.2f}", show_ctrl=True, ax=ax)
    fig2.suptitle("Sensibilité au shape factor du bras LE hub",
                  fontsize=13, fontweight="bold", y=1.02)

    plt.show()