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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_incidence(HX, w=None, tol_incid=1e-9, ax=None):
    """
    Matrice d'incidence géométrique W[i,j] = fraction de l'aire totale d'échange
    entre la cellule chaude i et la cellule froide j.
    Convention : cellule 0 = la plus proche du centre, pour le chaud ET le froid
    (numérotation par θ croissant des deux côtés).
    """
    P = HX.params
    H, n = P['H'], P['n_disc']
    rl, te = P['radial_limits'], P['theta_end']

    # --- géométrie ---------------------------------------------------------
    p = rl['hot'][1] - rl['sheet1'][0]            # pas (m/tour)
    a = p / (2*np.pi)
    rc1 = 0.5*sum(rl['sheet1'])
    rc2 = 0.5*sum(rl['sheet2'])
    th_h, th_c = te['hot'], te['cold']

    if w is None:
        w = getattr(HX, 'w', None)
        
    w = np.full(n, 1.0/n) if w is None else np.asarray(w, float)
    w = w / w.sum()

    # aire de tôle cumulée (arc exact)
    def G(theta, rc):
        rho = rc + a*theta
        return H*0.5*(rho*np.sqrt(rho**2 + a**2) + a**2*np.arcsinh(rho/a))/a

    G2 = lambda t: G(t, rc2) - G(0.0, rc2)          # sheet2 : chaud(θ) | froid(θ)
    G1 = lambda t: G(t, rc1) - G(2*np.pi, rc1)      # sheet1 : froid(θ) | chaud(θ-2π)

    # abscisse d'aire cumulée de chaque côté (croissante avec θ)
    tg_h = np.linspace(0, th_h, 40001)
    tg_c = np.linspace(0, th_c, 40001)
    S_H = G2(tg_h) + G1(tg_h + 2*np.pi)
    S_C = G2(np.minimum(tg_c, th_h)) + G1(np.maximum(tg_c, 2*np.pi))
    A_tot = S_H[-1]

    # frontières des cellules : mêmes w, cumul depuis le centre des deux côtés
    u = np.concatenate([[0.0], np.cumsum(w)])*A_tot
    thH = np.interp(u, S_H, tg_h)
    thC = np.interp(u, S_C, tg_c)

    # --- recouvrements -----------------------------------------------------
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            # (1) sheet2 : chaud(θ) en face de froid(θ)
            lo, hi = max(thH[i], thC[j]), min(thH[i+1], thC[j+1])
            if hi > lo:
                W[i, j] += G2(hi) - G2(lo)
            # (2) sheet1 : chaud(θ) en face de froid(θ+2π)
            lo, hi = max(thH[i], thC[j] - 2*np.pi), min(thH[i+1], thC[j+1] - 2*np.pi)
            if hi > lo:
                W[i, j] += G1(hi + 2*np.pi) - G1(lo + 2*np.pi)
    W /= A_tot
    B = W > tol_incid

    # --- tracé -------------------------------------------------------------
    if ax is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))
    else:
        ax1, ax2 = ax
        fig = ax1.figure

    Wm = np.ma.masked_where(~B, W)
    im = ax1.imshow(Wm, origin="upper", cmap="viridis",
                    norm=LogNorm(vmin=Wm.min(), vmax=Wm.max()))
    ax1.set_title("$w_{ij}$ (fraction de l'aire totale)")
    fig.colorbar(im, ax=ax1, fraction=0.046)
    ax2.imshow(B, origin="upper", cmap="Greys")
    ax2.set_title(f"Incidence : {B.sum()} contacts / {n*n}")

    for axx in (ax1, ax2):
        axx.set_xlabel("cellule froide j (0 = intérieur)")
        axx.set_ylabel("cellule chaude i (0 = intérieur)")
        axx.set_xticks(range(n)); axx.set_yticks(range(n))

    fig.tight_layout()
    plt.show()
    return W, B


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors


def plot_contacts(hx, ax=None, w=None, flow_from="inner", flow_from_hot=None,
                  label_pairs=True, fontsize=5, min_label_len=6e-3,
                  n_pts=2000, cmap="PRGn"):
    """Contacts entre cellules chaudes (i) et froides (j) dans l'échangeur à spirale.

    Chaque canal est découpé en sous-segments là où le partenaire change : un segment chaud est borné
    par les frontières des cellules chaudes ET par celles des cellules froides en regard
    (même tour via sheet2, tour suivant via sheet1), et inversement pour le froid.

    Géométrie (cf. radial_limits) : à l'angle θ on rencontre dans l'ordre
        froid(θ) | sheet2(θ) | chaud(θ) | sheet1(θ+2π) | froid(θ+2π)
    donc le chaud en θ touche le froid en θ (via sheet2) et le froid en θ+2π (via sheet1).

    Affichage : tôles colorées par j - i (0 = même numéro), étiquette "ij" = chaud i / froid j,
    posée dans le canal chaud : près de la face intérieure -> contact via sheet2,
    près de la face extérieure -> contact via sheet1.

    Retourne (ax, W, segs) : W[i, j] = fraction de l'aire d'échange totale, segs = liste de
    (θ_a, θ_b, i, j_sheet2, j_sheet1) en angle du canal chaud.
    """
    if getattr(hx, 'HTX_Type', None) != 'SPHE':
        raise ValueError(f"plot_contacts : HTX_Type = {getattr(hx, 'HTX_Type', None)!r}, attendu 'SPHE'.")

    P = hx.params
    H = P['H']
    p = P['C_canal_t'] + P['H_canal_t'] + 2 * P['t']
    a = p / (2 * np.pi)
    theta_max = 2 * np.pi * P['N']
    TWO_PI = 2 * np.pi
    th_h, th_c = P['theta_end']['hot'], P['theta_end']['cold']

    def r(element, th, pos=0.5):
        lo, hi = P['radial_limits'][element]
        return P['r0_in'] + a * np.asarray(th) + lo + pos * (hi - lo)

    F = lambda rr: 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a

    # --- vecteur w -------------------------------------------------------------
    if w is None:
        w = getattr(hx, 'w', None)
    if w is None:
        raise ValueError("plot_contacts : w absent (argument ou hx.w).")
    w = np.asarray(w, float).ravel()
    if abs(w.sum() - 1.0) > 1e-6:
        print(f"⚠ somme(w) = {w.sum():.6f} ≠ 1 : w normalisé")
        w = w / w.sum()
    nc = len(w)

    # --- frontières des cellules (identique à plot_spiral) -----------------------
    def boundaries(el, flow):
        th_end = P['theta_end'][el]
        s_of = lambda th: F(r(el, th)) - F(r(el, 0.0))
        L = float(s_of(th_end))
        frac = np.r_[0.0, np.cumsum(w)]
        frac[-1] = 1.0
        if flow == "outer":
            frac = 1.0 - frac
        s_k = frac * L
        th = s_k / L * th_end
        for _ in range(6):
            th = th - (s_of(th) - s_k) / np.sqrt(r(el, th)**2 + a**2)
        return np.clip(th, 0.0, th_end)

    tkC = boundaries("cold", flow_from)
    tkH = boundaries("hot", flow_from_hot or flow_from)
    ivC = (np.minimum(tkC[:-1], tkC[1:]), np.maximum(tkC[:-1], tkC[1:]))   # (lo, hi) par cellule
    ivH = (np.minimum(tkH[:-1], tkH[1:]), np.maximum(tkH[:-1], tkH[1:]))

    def find(iv, x):
        k = np.nonzero((iv[0] <= x) & (x < iv[1]))[0]
        return int(k[0]) if len(k) else -1

    def cuts(*arrays, upper):
        x = np.concatenate(arrays)
        x = np.unique(np.r_[0.0, x[(x > 0) & (x < upper)], upper])
        keep = np.r_[True, np.diff(x) > 1e-9]                               # supprime les segments dégénérés
        return x[keep]

    cuts_h = cuts(tkH, tkC, tkC - TWO_PI, upper=th_h)    # sous-découpage du canal chaud
    cuts_c = cuts(tkC, tkH, tkH + TWO_PI, upper=th_c)    # sous-découpage du canal froid

    # --- segments de contact + matrice W ----------------------------------------
    segs, W = [], np.zeros((nc, nc))
    for t0, t1 in zip(cuts_h[:-1], cuts_h[1:]):
        m = 0.5 * (t0 + t1)
        i, j2, j1 = find(ivH, m), find(ivC, m), find(ivC, m + TWO_PI)
        segs.append((t0, t1, i, j2, j1))
        if i >= 0 and j2 >= 0:      # sheet2 : chaud(θ) | froid(θ)
            W[i, j2] += H * (F(r("sheet2", t1)) - F(r("sheet2", t0)))
        if i >= 0 and j1 >= 0:      # sheet1 : chaud(θ) | froid(θ+2π)
            W[i, j1] += H * (F(r("sheet1", t1 + TWO_PI)) - F(r("sheet1", t0 + TWO_PI)))
    A_tot = W.sum()
    W = W / A_tot

    # --- tracé -------------------------------------------------------------------
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    def band(element, th0, th1, **kw):
        th = np.linspace(th0, th1, max(int(n_pts * (th1 - th0) / theta_max), 2))
        r_in, r_out = r(element, th, 0.0), r(element, th, 1.0)
        x = np.r_[r_out * np.cos(th), (r_in * np.cos(th))[::-1]]
        y = np.r_[r_out * np.sin(th), (r_in * np.sin(th))[::-1]]
        ax.fill(x, y, lw=0, **kw)

    def radial(el, th, **kw):
        th = np.asarray(th)
        ri, ro = r(el, th, 0.0), r(el, th, 1.0)
        for k in range(len(th)):
            ax.plot([ri[k] * np.cos(th[k]), ro[k] * np.cos(th[k])],
                    [ri[k] * np.sin(th[k]), ro[k] * np.sin(th[k])], **kw)

    # canaux : une teinte alternée par cellule
    for el, col, iv in (("hot", "tab:red", ivH), ("cold", "tab:blue", ivC)):
        for k in range(nc):
            band(el, iv[0][k], iv[1][k], color=col, alpha=0.40 if k % 2 == 0 else 0.18)
        radial(el, cuts_h if el == "hot" else cuts_c, color="0.35", lw=0.5, ls=":")   # sous-découpage
        radial(el, tkH if el == "hot" else tkC, color="k", lw=0.9)                    # frontières de cellules

    # tôles : base grise (zones sans contact chaud/froid), puis morceaux colorés par j - i
    band("sheet1", 0.0, P['theta_end']['sheet1'], color="0.80")
    band("sheet2", 0.0, P['theta_end']['sheet2'], color="0.80")

    pairs = [(i, j) for (_, _, i, j2, j1) in segs for j in (j2, j1) if j >= 0]
    dmax = max(1, max(abs(j - i) for i, j in pairs))
    norm = colors.Normalize(-dmax, dmax)
    cmap_ = plt.get_cmap(cmap)

    for t0, t1, i, j2, j1 in segs:
        if j2 >= 0:
            band("sheet2", t0, t1, color=cmap_(norm(j2 - i)))
        if j1 >= 0:
            band("sheet1", t0 + TWO_PI, t1 + TWO_PI, color=cmap_(norm(j1 - i)))

    # étiquettes "ij" dans le canal chaud
    if label_pairs:
        sep = "" if nc <= 10 else ","
        for t0, t1, i, j2, j1 in segs:
            tm = 0.5 * (t0 + t1)
            if r("hot", tm) * (t1 - t0) < min_label_len:
                continue
            for j, pos in ((j2, 0.22), (j1, 0.78)):      # 0.22 : côté sheet2 ; 0.78 : côté sheet1
                if j < 0:
                    continue
                rm = r("hot", tm, pos)
                ax.text(rm * np.cos(tm), rm * np.sin(tm), f"{i}{sep}{j}", fontsize=fontsize,
                        ha="center", va="center", color="k")

    sm = cm.ScalarMappable(norm=norm, cmap=cmap_)
    cb = ax.figure.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, label="j − i  (froid − chaud ; 0 = même numéro)")
    if dmax <= 10:
        cb.set_ticks(np.arange(-dmax, dmax + 1))

    phi = np.linspace(0, TWO_PI, 400)
    R_ext = P['D_ext'] / 2
    ax.plot(R_ext * np.cos(phi), R_ext * np.sin(phi), "k--", lw=0.6, label="D_ext")
    ax.fill([], [], color="tab:red", alpha=0.4, label="Canal chaud")
    ax.fill([], [], color="tab:blue", alpha=0.4, label="Canal froid")
    ax.fill([], [], color="0.8", label="Tôle sans contact")

    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"Contacts chaud i / froid j – {nc} cellules, A_tot = {A_tot:.4f} m², "
                 f"part de la diagonale = {np.trace(W):.1%}")
    ax.legend(loc="upper right", fontsize=8)

    print(f"Part de l'aire sur la diagonale (même numéro) : {np.trace(W):.1%}")
    print(f"Écart marginales vs w : lignes {np.abs(W.sum(1) - w).max():.1e}, "
          f"colonnes {np.abs(W.sum(0) - w).max():.1e}")
    
    return ax, W, segs
