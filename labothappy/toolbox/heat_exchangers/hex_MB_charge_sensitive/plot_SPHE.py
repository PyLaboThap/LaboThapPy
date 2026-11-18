import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.colors import LogNorm

# determine_cell_overlap_SPHE doit être importable ici (même module ou import)
# from cell_overlap_sphe import determine_cell_overlap_SPHE

SPIRAL_TYPES = ('SPHE', 'SCHE')
TWO_PI = 2 * np.pi


def _check_spiral(hx, fname):
    if getattr(hx, 'HTX_Type', None) not in SPIRAL_TYPES:
        raise ValueError(f"{fname} : HTX_Type = {getattr(hx, 'HTX_Type', None)!r}, attendu 'SPHE' ou 'SCHE'.")


def _pitch(P):
    """Pas radial : pitch_H pour le SCHE, C_canal_t + H_canal_t + 2t pour le SPHE (identiques via les alias)."""
    return P['pitch_H'] if 'pitch_H' in P else P['C_canal_t'] + P['H_canal_t'] + 2 * P['t']


def _wall_label(hx):
    """SPHE : tôles spiralées ; SCHE : nervures entre canaux (vue en plan d'une couche)."""
    return "Nervures" if hx.HTX_Type == 'SCHE' else "Tôles"


def _sphe_like_limits(P):
    """
    Limites radiales fictives pour tracer un SCHE en vue en plan comme un SPHE :
    deux bandes de canal séparées par deux parois minces, sur un pas pitch_H.
    Épaisseurs des canaux proportionnelles à leurs largeurs. Purement visuel (non exact).
    """
    p = _pitch(P)
    t = min(P.get('t_wall_min', 0.1 * p), 0.2 * p)            # épaisseur de paroi tracée
    wc, wh = P['C_canal_w'], P['H_canal_w']
    tc = (p - 2 * t) * wc / (wc + wh)
    th = (p - 2 * t) - tc
    base = min(P['radial_limits']['cold'][0], P['radial_limits']['hot'][0])
    if P['inner_channel'] == "cold":
        # froid(θ) | sheet2(θ) | chaud(θ) | sheet1(θ+2π) | froid(θ+2π)
        lim = {'cold': (0.0, tc), 'sheet2': (tc, tc + t), 'hot': (tc + t, tc + t + th), 'sheet1': (-t, 0.0)}
    else:
        # froid(θ-2π) | sheet1(θ) | chaud(θ) | sheet2(θ) | froid(θ)
        lim = {'hot': (0.0, th), 'sheet2': (th, th + t), 'cold': (th + t, th + t + tc), 'sheet1': (-t, 0.0)}
    return {k: (base + lo, base + hi) for k, (lo, hi) in lim.items()}


def _sheet1_orientation(P):
    """
    Retourne (shift_c, shift_s1, pos_s2, pos_s1) :
      - le chaud en θ touche le froid en θ via sheet2, et le froid en θ + shift_c via sheet1
        (sheet1 évaluée à θ + shift_s1) ;
      - pos_s2 / pos_s1 : position relative dans le canal chaud de la face portée par sheet2 / sheet1.
    """
    if P['inner_channel'] == "cold":
        # froid(θ) | sheet2(θ) | chaud(θ) | sheet1(θ+2π) | froid(θ+2π)
        return TWO_PI, TWO_PI, 0.22, 0.78
    # chaud intérieur : froid(θ-2π) | sheet1(θ) | chaud(θ) | sheet2(θ) | froid(θ)
    return -TWO_PI, 0.0, 0.78, 0.22


# =============================================================================
def plot_SCHE_section(hx, ax=None, n_turns=3, n_layers=3, first_turn=None, annotate=True, n_pts=60):
    """Coupe radiale (plan r-z, à θ = 0) des canaux d'un SCHE, comme la cellule unitaire de l'article.

    Convention de la cellule (hauteur pitch_V, largeur pitch_H) :
      - canaux elliptiques complets : largeur totale *_canal_w, hauteur totale *_canal_h
      - canal froid centré à z = (l + 1/4)·pitch_V, canal chaud à z = (l + 3/4)·pitch_V
        (rangées décalées d'une demi-cellule, comme radialement)
      - les deux canaux sont décalés radialement selon radial_limits (même pas pitch_H)
    pitch_H / pitch_V = écarts centre à centre entre deux canaux du même fluide ;
    chaud et froid sont décalés de (pitch_H/2, pitch_V/2).

    n_turns, n_layers : fenêtre tracée (None -> tous les tours / toutes les couches)
    first_turn        : premier tour tracé (None -> tours du milieu de la spirale)
    annotate          : cotes (largeurs / hauteurs des canaux, pitch_H, pitch_V, épaisseur de plaque)

    """
    _check_spiral(hx, "plot_SCHE_section")
    if hx.HTX_Type != 'SCHE':
        raise ValueError("plot_SCHE_section : réservé au SCHE.")

    P = hx.params
    p, pv = _pitch(P), P['pitch_V']
    a = p / TWO_PI
    mm = 1e3
    n_ch = P.get('n_ch', int(P['H'] // pv))           # canaux de chaque fluide empilés
    n_plates = P.get('n_plates', 2 * n_ch + 1)       # plaques : centre chaud -> centre froid, + 2 plaques de bord

    # tours présents dans la coupe θ = 0 pour chaque canal : θ = 2πk <= theta_end
    k_max = {el: int(np.floor(P['theta_end'][el] / TWO_PI + 1e-9)) for el in ('hot', 'cold')}
    n_tot = max(k_max.values()) + 1
    if n_turns is None or n_turns >= n_tot:
        k0, k1 = 0, n_tot
    else:
        k0 = (n_tot - n_turns) // 2 if first_turn is None else int(first_turn)
        k1 = min(k0 + n_turns, n_tot)
    l1 = n_ch if n_layers is None else min(n_layers, n_ch)

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))

    def r_center(el, k):
        return P['r0_in'] + a * TWO_PI * k + np.mean(P['radial_limits'][el])

    def ellipse(xc, zc, w, h):
        phi = np.linspace(0, TWO_PI, 2 * n_pts)
        return xc + 0.5 * w * np.cos(phi), zc + 0.5 * h * np.sin(phi)

    # tours effectivement tracés pour chaque fluide
    ks = {el: [k for k in range(k0, k1) if k <= k_max[el]] for el in ('cold', 'hot')}

    # métal : fond de la fenêtre (de la nervure avant le premier canal à celle après le dernier)
    t_rib = 0.25 * p                                   # marge de métal aux bords de la fenêtre
    edges = []
    for el, wch in (('cold', P['C_canal_w']), ('hot', P['H_canal_w'])):
        for k in ks[el]:
            xc = r_center(el, k)
            edges += [xc - 0.5 * wch, xc + 0.5 * wch]
    x_lo, x_hi = min(edges) - 0.5 * t_rib, max(edges) + 0.5 * t_rib
    # bas du métal : sous le bas des canaux de la première couche (+ marge)
    z_lo = min(0.0, min(0.25 * pv - 0.5 * P['C_canal_h'], 0.75 * pv - 0.5 * P['H_canal_h']) - 0.1 * pv)
    ax.add_patch(plt.Rectangle((x_lo * mm, z_lo * mm), (x_hi - x_lo) * mm, (l1 * pv - z_lo) * mm,
                               color="0.72", lw=0, zorder=0))

    cfg = {"cold": dict(w=P['C_canal_w'], h=P['C_canal_h'], z0=0.25, color="tab:blue", name="froid"),
           "hot":  dict(w=P['H_canal_w'], h=P['H_canal_h'], z0=0.75, color="tab:red",  name="chaud")}
    for el, c in cfg.items():
        for k in ks[el]:
            xc = r_center(el, k)
            for l in range(l1):
                x, z = ellipse(xc, (l + c['z0']) * pv, c['w'], c['h'])
                ax.fill(x * mm, z * mm, color=c['color'], alpha=0.55, lw=0.8, ec="k", zorder=2)

    # plans de liaison entre plaques : ils passent par le centre des canaux
    # (une plaque va du centre d'un canal chaud au centre du canal froid suivant)
    for l in range(l1):
        for z0 in (cfg['cold']['z0'], cfg['hot']['z0']):
            ax.plot([x_lo * mm, x_hi * mm], [(l + z0) * pv * mm] * 2, color="0.25", lw=0.7, ls="--", zorder=3)

    # --- cotes ---------------------------------------------------------------------
    ext = 0.04 * max(p, pv)                            # dépassement des lignes d'attache
    arrow = dict(arrowstyle="<->", lw=0.8, color="k", shrinkA=0, shrinkB=0)
    kw_ext = dict(color="0.2", lw=0.5, zorder=1.5)     # lignes d'attache (sous les canaux)
    kw_txt = dict(fontsize=7.5, zorder=6, linespacing=1.0, multialignment="center", textcoords="offset points",
                  bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85))

    def cote_h(x1, x2, z1, z2, zl, txt, below=False):
        """Cote horizontale de x1 à x2 sur la ligne z = zl ; attaches depuis (x1, z1) et (x2, z2). Unités : m."""
        s = -1 if below else 1
        for x, z in ((x1, z1), (x2, z2)):
            ax.plot([x * mm] * 2, [z * mm, (zl + s * ext) * mm], **kw_ext)
        ax.annotate("", xy=(x1 * mm, zl * mm), xytext=(x2 * mm, zl * mm), arrowprops=arrow, zorder=5)
        ax.annotate(txt, xy=(0.5 * (x1 + x2) * mm, zl * mm), xytext=(0, -3 if below else 3),
                    ha="center", va="top" if below else "bottom", **kw_txt)

    def cote_v(z1, z2, x1, x2, xl, txt, left=True):
        """Cote verticale de z1 à z2 sur la ligne x = xl ; attaches depuis (x1, z1) et (x2, z2). Unités : m."""
        s = -1 if left else 1
        for x, z in ((x1, z1), (x2, z2)):
            ax.plot([x * mm, (xl + s * ext) * mm], [z * mm] * 2, **kw_ext)
        ax.annotate("", xy=(xl * mm, z1 * mm), xytext=(xl * mm, z2 * mm), arrowprops=arrow, zorder=5)
        ax.annotate(txt, xy=(xl * mm, 0.5 * (z1 + z2) * mm), xytext=(4 * s, 0),
                    ha="right" if left else "left", va="center", **kw_txt)

    if annotate:
        wc, hc = cfg['cold']['w'], cfg['cold']['h']
        wh, hh = cfg['hot']['w'], cfg['hot']['h']
        z_cold = lambda l: (l + cfg['cold']['z0']) * pv
        z_hot = lambda l: (l + cfg['hot']['z0']) * pv
        xl = x_lo - 0.08 * p                             # ligne de cotes à gauche
        xr = x_hi + 0.08 * p                             # ligne de cotes à droite

        # largeur du canal froid : sous la première couche
        if ks['cold']:
            xc = r_center('cold', ks['cold'][0])
            cote_h(xc - 0.5 * wc, xc + 0.5 * wc, z_cold(0), z_cold(0), z_lo - 0.2 * pv,
                   f"w froid = {wc*mm:.2f} mm", below=True)

        # largeur du canal chaud : au-dessus de la dernière couche
        z_top = l1 * pv + 0.15 * pv
        if ks['hot']:
            xc = r_center('hot', ks['hot'][0])
            cote_h(xc - 0.5 * wh, xc + 0.5 * wh, z_hot(l1 - 1), z_hot(l1 - 1), z_top,
                   f"w chaud = {wh*mm:.2f} mm")
            z_top += 0.35 * pv

        # pitch_H : entre deux centres de canaux froids successifs (dernière couche)
        kc = k0 if k0 <= k_max['cold'] else k0 - 1
        xa, xb = r_center('cold', kc), r_center('cold', kc + 1)
        cote_h(xa, xb, z_cold(l1 - 1), z_cold(l1 - 1), z_top, f"pitch_H = {p*mm:.2f} mm")

        # pitch_V : entre deux centres de canaux froids superposés (à gauche)
        x_att = r_center('cold', ks['cold'][0]) if ks['cold'] else x_lo
        cote_v(z_cold(0), z_cold(1), x_att, x_lo if l1 < 2 else x_att, xl, f"pitch_V\n{pv*mm:.2f} mm")

        # épaisseur d'une plaque : entre deux plans de liaison (à gauche, dernière couche)
        if l1 >= 2:
            cote_v(z_cold(l1 - 1), z_hot(l1 - 1), x_lo, x_lo, xl, f"e plaque\n{0.5*pv*mm:.2f} mm")

        # hauteurs des canaux (à droite) : froid sur la première couche, chaud sur la dernière
        if ks['cold']:
            xc = r_center('cold', ks['cold'][-1])
            cote_v(z_cold(0) - 0.5 * hc, z_cold(0) + 0.5 * hc, xc, xc, xr,
                   f"h froid\n{hc*mm:.2f} mm", left=False)
        if ks['hot']:
            xc = r_center('hot', ks['hot'][-1])
            l_h = l1 - 1
            x_line = xr if l1 >= 2 else xr + 0.35 * max(p, pv)
            cote_v(z_hot(l_h) - 0.5 * hh, z_hot(l_h) + 0.5 * hh, xc, xc, x_line,
                   f"h chaud\n{hh*mm:.2f} mm", left=False)


    marg = max(0.5 * p, 1.2 * pv)
    ax.set_aspect("equal")
    ax.set_xlim(x_lo * mm - marg * mm, x_hi * mm + marg * mm)
    ax.set_ylim((z_lo - 0.5 * pv) * mm, (l1 + 0.9) * pv * mm)
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("z [mm]")
    ax.set_title("Coupe radiale des canaux (θ = 0)")
    return ax

# =============================================================================
def plot_spiral(hx, ax=None, n=2000, w=None, values_c=None, values_h=None,
                flow_from="inner", flow_from_hot=None, label_cells=True, tol=1e-6,
                section=True, ax_section=None, **section_kw):
    """Vue en plan de l'échangeur à spirale (SPHE ou SCHE), avec la discrétisation des DEUX fluides.

    SPHE : tôles spiralées de hauteur H.
    SCHE : vue en plan d'UNE couche ; les n_plates couches empilées sont identiques.
           Les bandes grises sont les nervures entre canaux (largeur e_H / 2 de part et d'autre).

    Lit hx.params (rempli par setup_geom) : t, C_canal_t, H_canal_t (ou pitch_H), r0_in, D_ext, N,
    radial_limits, theta_end, inner_channel (+ L_cold, L_hot si présents, pour contrôle).

    w            : vecteur des parts ; None -> hx.w ; absent -> tracé sans discrétisation.
    values_c/h   : une valeur par cellule pour colorer les cellules ; None -> teintes alternées.
    flow_from    : 'inner' = cellule 0 côté noyau (θ = 0) ; 'outer' = cellule 0 côté extérieur.
    flow_from_hot: idem pour le chaud si son sens diffère.
    label_cells  : numérote les cellules si <= 40 cellules.
    section      : SCHE uniquement, ajoute la coupe radiale des canaux (plot_SCHE_section)
    ax_section   : axe pour la coupe ; None -> créé à côté de la vue en plan
    section_kw   : transmis à plot_SCHE_section (n_turns, n_layers, first_turn, annotate)
    """
    _check_spiral(hx, "plot_spiral")

    P = hx.params
    p = _pitch(P)
    a = p / TWO_PI
    theta_max = TWO_PI * P['N']
    with_section = section and hx.HTX_Type == 'SCHE'
    approx = hx.HTX_Type == 'SCHE'
    if approx:                                        # SCHE tracé comme un SPHE (approximation visuelle)
        P = dict(P, radial_limits=_sphe_like_limits(P))

    def r(element, th, pos=0.5):
        """Rayon de l'élément à l'angle th ; pos = 0 face intérieure, 0.5 médiane, 1 face extérieure."""
        lo, hi = P['radial_limits'][element]
        return P['r0_in'] + a * np.asarray(th) + lo + pos * (hi - lo)

    if ax is None:
        if with_section and ax_section is None:
            _, (ax, ax_section) = plt.subplots(1, 2, figsize=(14, 6.5),
                                               gridspec_kw=dict(width_ratios=[1.1, 1]))
        else:
            _, ax = plt.subplots(figsize=(7, 7))
    elif with_section and ax_section is None:
        _, ax_section = plt.subplots(figsize=(7, 5))

    def band(element, th0, th1, **kw):
        th = np.linspace(th0, th1, max(int(n * (th1 - th0) / theta_max), 2))
        r_in, r_out = r(element, th, 0.0), r(element, th, 1.0)
        x = np.r_[r_out * np.cos(th), (r_in * np.cos(th))[::-1]]
        y = np.r_[r_out * np.sin(th), (r_in * np.sin(th))[::-1]]
        ax.fill(x, y, lw=0, **kw)

    F = lambda rr: 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a

    def boundaries(el, flow):
        th_end = P['theta_end'][el]
        if th_end <= 0:
            return None
        s_of = lambda th: F(r(el, th)) - F(r(el, 0.0))
        L = float(s_of(th_end))
        key = 'L_cold' if el == 'cold' else 'L_hot'
        if not approx and key in P and abs(P[key] - L) > 1e-6 * L:
            print(f"⚠ {key}(params) = {P[key]:.4f} m ≠ longueur recalculée {L:.4f} m")
        frac = np.r_[0.0, np.cumsum(w)]
        frac[-1] = 1.0
        if flow == "outer":
            frac = 1.0 - frac
        s_k = frac * L
        th = s_k / L * th_end
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

    band("sheet1", 0.0, P['theta_end']['sheet1'], color="0.25", label=_wall_label(hx))
    band("sheet2", 0.0, P['theta_end']['sheet2'], color="0.25")

    phi = np.linspace(0, TWO_PI, 400)
    R_ext = P['D_ext'] / 2
    ax.plot(R_ext * np.cos(phi), R_ext * np.sin(phi), "k--", lw=0.6, label="D_ext")

    R_out = float(r("sheet2", theta_max, 1.0))
    if not approx and abs(R_out - R_ext) > tol:
        msg = (f"⚠ face externe de sheet2 : R = {R_out*1e3:.2f} mm ≠ D_ext/2 = {R_ext*1e3:.2f} mm "
               f"(écart {1e3*(R_out - R_ext):+.2f} mm)")
        print(msg)
        ax.text(0.5, 0.01, msg, transform=ax.transAxes, ha="center", va="bottom", fontsize=8, color="crimson")

    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    title = f"Spirale – p = {p*1e3:.2f} mm, N = {P['N']:.2f} tours, canal intérieur : {P['inner_channel']}"
    if hx.HTX_Type == 'SCHE':
        title += f"\nVue d'une couche – {P.get('n_plates', '?')} couches de pas {P['pitch_V']*1e3:.2f} mm"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)

    if with_section:
        plot_SCHE_section(hx, ax=ax_section, **section_kw)
        ax.plot([P['r0_in'], R_ext], [0, 0], color="k", lw=1.2, ls="-.", label="_coupe")   # trace de la coupe
        ax.figure.tight_layout()
    return ax



# =============================================================================
def plot_incidence(HX, w=None, tol_incid=1e-9, ax=None):
    """
    Matrice d'incidence géométrique W[i,j] = fraction de l'aire totale d'échange
    entre la cellule chaude i et la cellule froide j (SPHE ou SCHE).
    Convention : cellule 0 = la plus proche du centre, pour le chaud ET le froid.

    W est calculée par determine_cell_overlap_SPHE, c'est-à-dire exactement la matrice
    utilisée par le solveur (frontières de cellules en longueur de canal, les deux
    orientations de inner_channel gérées).
    """
    _check_spiral(HX, "plot_incidence")

    if w is not None:                                  # w imposé : calcul temporaire, hx.w restauré
        w_saved = getattr(HX, 'w', None)
        HX.w = np.asarray(w, float)
        try:
            W = determine_cell_overlap_SPHE(HX)
        finally:
            HX.w = w_saved
    else:
        if getattr(HX, 'w', None) is None:
            HX.w = np.full(HX.params['n_disc'], 1.0 / HX.params['n_disc'])
        W = determine_cell_overlap_SPHE(HX)

    n = W.shape[0]
    B = W > tol_incid

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
        if n <= 40:
            axx.set_xticks(range(n)); axx.set_yticks(range(n))

    fig.tight_layout()
    plt.show()
    return W, B


# =============================================================================
def plot_contacts(hx, ax=None, w=None, flow_from="inner", flow_from_hot=None,
                  label_pairs=True, fontsize=5, min_label_len=6e-3,
                  n_pts=2000, cmap="PRGn"):
    """Contacts entre cellules chaudes (i) et froides (j) dans l'échangeur à spirale (SPHE ou SCHE).

    Chaque canal est découpé en sous-segments là où le partenaire change.

    Géométrie (cf. radial_limits), selon inner_channel :
        'cold' : froid(θ) | sheet2(θ) | chaud(θ) | sheet1(θ+2π) | froid(θ+2π)
        'hot'  : froid(θ-2π) | sheet1(θ) | chaud(θ) | sheet2(θ) | froid(θ)
    Le chaud en θ touche donc toujours le froid en θ via sheet2, et le froid en θ ± 2π via sheet1.

    SCHE : sheet1/sheet2 sont les deux parois d'échange du canal chaud ; l'aire affichée
    utilise H_contact (= n_plates · P_ex / 2) au lieu de H.

    Retourne (ax, W, segs) : W[i, j] = fraction de l'aire d'échange totale, segs = liste de
    (θ_a, θ_b, i, j_sheet2, j_sheet1) en angle du canal chaud.
    """
    _check_spiral(hx, "plot_contacts")

    P = hx.params
    H = P.get('H_contact', P['H'])                  # hauteur d'échange par paroi
    p = _pitch(P)
    a = p / TWO_PI
    theta_max = TWO_PI * P['N']
    th_h, th_c = P['theta_end']['hot'], P['theta_end']['cold']
    shift_c, shift_s1, pos_s2, pos_s1 = _sheet1_orientation(P)

    def r(element, th, pos=0.5):
        lo, hi = P['radial_limits'][element]
        return P['r0_in'] + a * np.asarray(th) + lo + pos * (hi - lo)

    F = lambda rr: 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a

    if w is None:
        w = getattr(hx, 'w', None)
    if w is None:
        raise ValueError("plot_contacts : w absent (argument ou hx.w).")
    w = np.asarray(w, float).ravel()
    if abs(w.sum() - 1.0) > 1e-6:
        print(f"⚠ somme(w) = {w.sum():.6f} ≠ 1 : w normalisé")
        w = w / w.sum()
    nc = len(w)

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
    ivC = (np.minimum(tkC[:-1], tkC[1:]), np.maximum(tkC[:-1], tkC[1:]))
    ivH = (np.minimum(tkH[:-1], tkH[1:]), np.maximum(tkH[:-1], tkH[1:]))

    def find(iv, x):
        k = np.nonzero((iv[0] <= x) & (x < iv[1]))[0]
        return int(k[0]) if len(k) else -1

    def cuts(*arrays, upper):
        x = np.concatenate(arrays)
        x = np.unique(np.r_[0.0, x[(x > 0) & (x < upper)], upper])
        keep = np.r_[True, np.diff(x) > 1e-9]
        return x[keep]

    cuts_h = cuts(tkH, tkC, tkC - shift_c, upper=th_h)    # sous-découpage du canal chaud
    cuts_c = cuts(tkC, tkH, tkH + shift_c, upper=th_c)    # sous-découpage du canal froid

    # --- segments de contact + matrice W ----------------------------------------
    segs, W = [], np.zeros((nc, nc))
    for t0, t1 in zip(cuts_h[:-1], cuts_h[1:]):
        m = 0.5 * (t0 + t1)
        i, j2, j1 = find(ivH, m), find(ivC, m), find(ivC, m + shift_c)
        segs.append((t0, t1, i, j2, j1))
        if i >= 0 and j2 >= 0:      # sheet2 : chaud(θ) | froid(θ)
            W[i, j2] += H * (F(r("sheet2", t1)) - F(r("sheet2", t0)))
        if i >= 0 and j1 >= 0:      # sheet1 : chaud(θ) | froid(θ ± 2π)
            W[i, j1] += H * (F(r("sheet1", t1 + shift_s1)) - F(r("sheet1", t0 + shift_s1)))
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

    for el, col, iv in (("hot", "tab:red", ivH), ("cold", "tab:blue", ivC)):
        for k in range(nc):
            band(el, iv[0][k], iv[1][k], color=col, alpha=0.40 if k % 2 == 0 else 0.18)
        radial(el, cuts_h if el == "hot" else cuts_c, color="0.35", lw=0.5, ls=":")
        radial(el, tkH if el == "hot" else tkC, color="k", lw=0.9)

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
            band("sheet1", t0 + shift_s1, t1 + shift_s1, color=cmap_(norm(j1 - i)))

    # étiquettes "ij" dans le canal chaud, du côté de la paroi concernée
    if label_pairs:
        sep = "" if nc <= 10 else ","
        for t0, t1, i, j2, j1 in segs:
            tm = 0.5 * (t0 + t1)
            if r("hot", tm) * (t1 - t0) < min_label_len:
                continue
            for j, pos in ((j2, pos_s2), (j1, pos_s1)):
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
    ax.fill([], [], color="0.8", label=f"{_wall_label(hx)} sans contact")

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