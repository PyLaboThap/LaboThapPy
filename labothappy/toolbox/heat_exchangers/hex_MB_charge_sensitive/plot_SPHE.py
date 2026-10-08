# import numpy as np
# import matplotlib.pyplot as plt


# def plot_spiral(hx, ax=None, n=2000, tol=1e-6):
#     """Vue en plan de l'échangeur à spirale (HTX_Type == 'SPHE').

#     Prend l'objet échangeur lui-même et lit hx.params, tel que rempli par setup_geom :
#         t, C_canal_t, H_canal_t, r0_in, D_ext, N, radial_limits, theta_end.

#     Utilisable en fonction  : plot_spiral(hx)
#     ou collée comme méthode : def plot_spiral(self, ax=None, ...)  puis  hx.plot_spiral()
#     """
#     if getattr(hx, 'HTX_Type', None) != 'SPHE':
#         raise ValueError(f"plot_spiral : HTX_Type = {getattr(hx, 'HTX_Type', None)!r}, attendu 'SPHE'.")

#     P = hx.params
#     t = P['t']
#     p = P['C_canal_t'] + P['H_canal_t'] + 2 * t        # pas radial (non stocké dans params : recalculé)
#     a = p / (2 * np.pi)
#     theta_max = 2 * np.pi * P['N']

#     def r(element, th, pos=0.5):
#         """Rayon de l'élément à l'angle th ; pos = 0 face intérieure, 0.5 médiane, 1 face extérieure."""
#         lo, hi = P['radial_limits'][element]
#         return P['r0_in'] + a * th + lo + pos * (hi - lo)

#     if ax is None:
#         _, ax = plt.subplots(figsize=(7, 7))

#     # canaux puis tôles (ordre de tracé)
#     styles = {
#         "hot":    ("tab:red",  0.35, "Canal chaud"),
#         "cold":   ("tab:blue", 0.35, "Canal froid"),
#         "sheet1": ("0.25",     1.0,  "Tôles"),
#         "sheet2": ("0.25",     1.0,  None),
#     }
#     for el, (color, alpha, label) in styles.items():
#         th = np.linspace(0.0, P['theta_end'][el], n)
#         r_in, r_out = r(el, th, 0.0), r(el, th, 1.0)
#         x = np.r_[r_out * np.cos(th), (r_in * np.cos(th))[::-1]]
#         y = np.r_[r_out * np.sin(th), (r_in * np.sin(th))[::-1]]
#         ax.fill(x, y, color=color, alpha=alpha, lw=0, label=label)

#     # cercle de contrôle D_ext
#     phi = np.linspace(0, 2 * np.pi, 400)
#     R_ext = P['D_ext'] / 2
#     ax.plot(R_ext * np.cos(phi), R_ext * np.sin(phi), "k--", lw=0.6, label="D_ext")

#     # contrôle : la face externe de la tôle 2 doit tomber sur D_ext
#     R_out = float(r("sheet2", theta_max, 1.0))
#     if abs(R_out - R_ext) > tol:
#         msg = (f"⚠ face externe de la tôle 2 : R = {R_out*1e3:.2f} mm ≠ D_ext/2 = {R_ext*1e3:.2f} mm "
#                f"(écart {1e3*(R_out - R_ext):+.2f} mm)")
#         print(msg)
#         ax.text(0.5, 0.01, msg, transform=ax.transAxes, ha="center", va="bottom", fontsize=8, color="crimson")

#     ax.set_aspect("equal")
#     ax.set_xlabel("x [m]")
#     ax.set_ylabel("y [m]")
#     ax.set_title(f"Spirale – p = {p*1e3:.1f} mm, N = {P['N']:.2f} tours, canal intérieur : {P['inner_channel']}")
#     ax.legend(loc="upper right", fontsize=8)
#     return ax

#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors


def plot_spiral(hx, ax=None, n=2000, w=None, values_c=None, values_h=None,
                flow_from="inner", flow_from_hot=None, label_cells=True, tol=1e-6):
    """Vue en plan de l'échangeur à spirale (HTX_Type == 'SPHE'), avec la discrétisation des DEUX fluides.

    Lit hx.params (rempli par setup_geom) : t, C_canal_t, H_canal_t, r0_in, D_ext, N,
    radial_limits, theta_end, inner_channel (+ L_cold, L_hot si présents, pour contrôle).

    Discrétisation : w (par défaut hx.w) = parts de longueur thermique, somme = 1, identiques côté chaud
    et côté froid. Pour chaque canal, la frontière de la cellule k est placée à l'abscisse curviligne
    cumsum(w)[k] · L_canal le long de sa ligne médiane, puis ramenée à l'angle θ_k (inversion de s(θ), Newton).

    w            : vecteur des parts ; None -> hx.w ; absent -> tracé sans discrétisation.
    values_c/h   : une valeur par cellule (ex. T_froid, T_chaud) pour colorer les cellules ; None -> teintes alternées.
    flow_from    : 'inner' = cellule 0 côté noyau (θ = 0) ; 'outer' = cellule 0 côté extérieur (côté froid,
                   et côté chaud si flow_from_hot n'est pas précisé).
    flow_from_hot: idem pour le chaud si son sens diffère.
    label_cells  : numérote les cellules (0 = premier élément de w) si <= 40 cellules.

    Utilisable en fonction : plot_spiral(hx)   ou collée comme méthode : def plot_spiral(self, ...).
    """
    if getattr(hx, 'HTX_Type', None) != 'SPHE':
        raise ValueError(f"plot_spiral : HTX_Type = {getattr(hx, 'HTX_Type', None)!r}, attendu 'SPHE'.")

    P = hx.params
    t = P['t']
    p = P['C_canal_t'] + P['H_canal_t'] + 2 * t        # pas radial (non stocké dans params : recalculé)
    a = p / (2 * np.pi)
    theta_max = 2 * np.pi * P['N']

    def r(element, th, pos=0.5):
        """Rayon de l'élément à l'angle th ; pos = 0 face intérieure, 0.5 médiane, 1 face extérieure."""
        lo, hi = P['radial_limits'][element]
        return P['r0_in'] + a * np.asarray(th) + lo + pos * (hi - lo)

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    def band(element, th0, th1, **kw):
        th = np.linspace(th0, th1, max(int(n * (th1 - th0) / theta_max), 2))
        r_in, r_out = r(element, th, 0.0), r(element, th, 1.0)
        x = np.r_[r_out * np.cos(th), (r_in * np.cos(th))[::-1]]
        y = np.r_[r_out * np.sin(th), (r_in * np.sin(th))[::-1]]
        ax.fill(x, y, lw=0, **kw)

    # --- discrétisation : frontières θ_k de chaque canal ---
    F = lambda rr: 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a   # ∫√(r²+a²)dr / a

    def boundaries(el, flow):
        th_end = P['theta_end'][el]
        if th_end <= 0:
            return None
        s_of = lambda th: F(r(el, th)) - F(r(el, 0.0))                    # abscisse curviligne
        L = float(s_of(th_end))
        key = 'L_cold' if el == 'cold' else 'L_hot'
        if key in P and abs(P[key] - L) > 1e-6 * L:
            print(f"⚠ {key}(params) = {P[key]:.4f} m ≠ longueur recalculée {L:.4f} m")
        frac = np.r_[0.0, np.cumsum(w)]
        frac[-1] = 1.0
        if flow == "outer":
            frac = 1.0 - frac
        s_k = frac * L
        th = s_k / L * th_end                                              # point de départ
        for _ in range(6):                                                 # Newton : ds/dθ = √(r²+a²)
            th = th - (s_of(th) - s_k) / np.sqrt(r(el, th)**2 + a**2)
        return np.clip(th, 0.0, th_end)

    if w is None:
        w = getattr(hx, 'w', None)
    th_k = {}
    if w is not None:
        w = np.asarray(w, float).ravel()
        if abs(w.sum() - 1.0) > 1e-6:
            print(f"⚠ somme(w) = {w.sum():.6f} ≠ 1 : w normalisé pour le tracé")
            w = w / w.sum()
        th_k = {"cold": boundaries("cold", flow_from),
                "hot": boundaries("hot", flow_from_hot or flow_from)}
    elif values_c is not None or values_h is not None:
        raise ValueError("values_c / values_h fournis mais aucun vecteur w (hx.w absent).")

    # --- tracé : canaux (avec cellules) puis tôles ---
    cfg = {
        "hot":  dict(color="tab:red",  cmap=cm.Reds,  values=values_h, txt="darkred", name="chaud", alpha=0.30),
        "cold": dict(color="tab:blue", cmap=cm.Blues, values=values_c, txt="navy",    name="froid", alpha=0.35),
    }
    for el, c in cfg.items():
        tk = th_k.get(el)
        if tk is None:
            band(el, 0.0, P['theta_end'][el], color=c['color'], alpha=c['alpha'], label=f"Canal {c['name']}")
            continue

        vals = c['values']
        if vals is not None:
            vals = np.asarray(vals, float).ravel()
            if len(vals) != len(w):
                raise ValueError(f"values pour le {c['name']} : {len(vals)} valeurs pour {len(w)} cellules.")
            norm = colors.Normalize(vals.min(), vals.max())
        for k in range(len(w)):
            lo_, hi_ = sorted((tk[k], tk[k + 1]))
            if vals is not None:
                band(el, lo_, hi_, color=c['cmap'](0.25 + 0.75 * norm(vals[k])), alpha=0.9)
            else:
                band(el, lo_, hi_, color=c['color'], alpha=0.55 if k % 2 == 0 else 0.25)
        ax.fill([], [], color=c['color'], alpha=0.4, label=f"Canal {c['name']} ({len(w)} cellules)")
        if vals is not None:
            sm = cm.ScalarMappable(norm=colors.Normalize(vals.min(), vals.max()),
                                   cmap=colors.LinearSegmentedColormap.from_list(
                                       "_", c['cmap'](np.linspace(0.25, 1.0, 64))))
            ax.figure.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, label=f"valeur par cellule ({c['name']})")

        # frontières de cellules (traits radiaux à travers le canal)
        ri, ro = r(el, tk, 0.0), r(el, tk, 1.0)
        for k in range(len(tk)):
            ax.plot([ri[k] * np.cos(tk[k]), ro[k] * np.cos(tk[k])],
                    [ri[k] * np.sin(tk[k]), ro[k] * np.sin(tk[k])],
                    color="k", lw=0.6, label=f"_bnd_{el}")

        if label_cells and len(w) <= 40:
            for k in range(len(w)):
                th_m = 0.5 * (tk[k] + tk[k + 1])
                rm = r(el, th_m)
                ax.text(rm * np.cos(th_m), rm * np.sin(th_m), str(k), fontsize=6,
                        ha="center", va="center", color=c['txt'])

    band("sheet1", 0.0, P['theta_end']['sheet1'], color="0.25", label="Tôles")
    band("sheet2", 0.0, P['theta_end']['sheet2'], color="0.25")

    # cercle de contrôle D_ext
    phi = np.linspace(0, 2 * np.pi, 400)
    R_ext = P['D_ext'] / 2
    ax.plot(R_ext * np.cos(phi), R_ext * np.sin(phi), "k--", lw=0.6, label="D_ext")

    # contrôle : la face externe de la tôle 2 doit tomber sur D_ext
    R_out = float(r("sheet2", theta_max, 1.0))
    if abs(R_out - R_ext) > tol:
        msg = (f"⚠ face externe de la tôle 2 : R = {R_out*1e3:.2f} mm ≠ D_ext/2 = {R_ext*1e3:.2f} mm "
               f"(écart {1e3*(R_out - R_ext):+.2f} mm)")
        print(msg)
        ax.text(0.5, 0.01, msg, transform=ax.transAxes, ha="center", va="bottom", fontsize=8, color="crimson")

    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"Spirale – p = {p*1e3:.1f} mm, N = {P['N']:.2f} tours, canal intérieur : {P['inner_channel']}")
    ax.legend(loc="upper right", fontsize=8)
    return ax