import numpy as np

def determine_cell_overlap_STHE(w_cells, n_passes, shell_side):
    w_cells = np.asarray(w_cells, dtype=float)
    n_cells = w_cells.size
    total_w = w_cells.sum()

    # Shell-side relative positions
    shell_edges = np.empty(n_cells + 1)
    shell_edges[0] = 0.0
    shell_edges[1:] = np.cumsum(w_cells) / total_w
    shell_edges[-1] = 1.0

    shell_start = shell_edges[:-1]
    shell_end   = shell_edges[1:]

    # Tube-side (unfolded over passes)
    tube_edges = np.empty(n_cells + 1)
    tube_edges[0] = 0.0
    tube_edges[1:] = np.cumsum(w_cells) * n_passes / total_w

    overlap_matrix = np.zeros((n_cells, n_cells), dtype=float)

    for i in range(n_cells):
        start = tube_edges[i]
        end   = tube_edges[i + 1]

        p_start = int(start)
        p_end   = int(np.ceil(end))

        for p in range(p_start, p_end):
            # relative to current pass
            rel_start = start - p
            rel_end   = min(end - p, 1.0)

            if p & 1:  # reverse pass (faster than % 2)
                sub_start = 1.0 - rel_end
                sub_end   = 1.0 - rel_start
            else:
                sub_start = rel_start
                sub_end   = rel_end

            # vectorized overlap
            overlap = np.minimum(sub_end, shell_end) - np.maximum(sub_start, shell_start)
            overlap_matrix[i] += np.maximum(overlap, 0.0)

    # If the hot side is currently on the tube side, transpose it
    matrix = overlap_matrix
    # if shell_side != 'H':
    #     matrix = overlap_matrix.T  # ensure j corresponds to hot side
    
    return total_w * matrix / n_passes

#%%

TWO_PI = 2 * np.pi

def _find_cells(lo, hi, x):
    """
    Version vectorisée de find() : pour chaque x, indice k de la cellule telle que
    lo[k] <= x < hi[k], ou -1 si aucune.
    Les cellules forment une partition de [0, θ_end] : après tri par (lo, hi),
    la cellule candidate est la dernière dont lo <= x (searchsorted). En cas de
    lo égaux (cellule de largeur nulle), le tri secondaire sur hi place la
    cellule non dégénérée en dernier, donc c'est elle qui est retenue.
    """
    order = np.lexsort((hi, lo))
    lo_s, hi_s = lo[order], hi[order]
    idx = np.searchsorted(lo_s, x, side='right') - 1
    idx_c = np.clip(idx, 0, len(lo) - 1)
    ok = (idx >= 0) & (x < hi_s[idx_c])
    return np.where(ok, order[idx_c], -1)


def determine_cell_overlap_SPHE(hx, flow_from="inner"):
    """
    Matrice W (nc x nc) des surfaces de contact normalisées entre cellules
    chaudes (lignes) et froides (colonnes), somme(W) = 1.

    Version vectorisée : même résultat que la version à boucle, sans boucle
    Python sur les segments.
    """
    if getattr(hx, 'HTX_Type', None) not in ('SPHE', 'SCHE'):
        raise ValueError(f"determine_cell_overlap_SPHE : HTX_Type = {getattr(hx, 'HTX_Type', None)!r}, attendu 'SPHE' ou 'SCHE'.")

    P = hx.params
    H = P['H']
    a = (P['b_C'] + P['b_H'] + 2 * P['t']) / TWO_PI
    r0 = P['r0_in']
    RL = P['radial_limits']
    th_h = P['theta_end']['hot']

    def r(element, th, pos=0.5):
        lo, hi = RL[element]
        return r0 + a * th + lo + pos * (hi - lo)

    def F(rr):
        return 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a

    # --- vecteur w -------------------------------------------------------------
    w = getattr(hx, 'w', None)
    if w is None:
        raise ValueError("plot_contacts : w absent (argument ou hx.w).")
    w = np.asarray(w, float).ravel()
    if abs(w.sum() - 1.0) > 1e-6:
        w = w / w.sum()
    nc = len(w)

    # --- frontières des cellules (identique à plot_spiral) -----------------------
    frac = np.r_[0.0, np.cumsum(w)]
    frac[-1] = 1.0
    if flow_from == "outer":
        frac = 1.0 - frac

    def boundaries(el):
        th_end = P['theta_end'][el]
        F0 = F(r(el, 0.0))
        L = float(F(r(el, th_end)) - F0)
        s_k = frac * L
        th = s_k / L * th_end
        for _ in range(6):                                   # Newton sur s(θ)
            th = th - (F(r(el, th)) - F0 - s_k) / np.sqrt(r(el, th)**2 + a**2)
        return np.clip(th, 0.0, th_end)

    tkC = boundaries("cold")
    tkH = boundaries("hot")
    loC, hiC = np.minimum(tkC[:-1], tkC[1:]), np.maximum(tkC[:-1], tkC[1:])
    loH, hiH = np.minimum(tkH[:-1], tkH[1:]), np.maximum(tkH[:-1], tkH[1:])

    # --- sheet1 orientation depends on which channel is inside ------------------
    if P['inner_ch'] == "cold":
        # hot is outer: hot(θ) | sheet1(θ+2π) | cold(θ+2π)
        shift_c, shift_s1 = TWO_PI, TWO_PI
    else:
        # hot is inner: cold(θ-2π) | sheet1(θ) | hot(θ)
        shift_c, shift_s1 = -TWO_PI, 0.0

    # --- sub-division of the hot channel -----------------------------------------
    x = np.concatenate((tkH, tkC, tkC - shift_c))
    x = np.unique(np.r_[0.0, x[(x > 0) & (x < th_h)], th_h])
    x = x[np.r_[True, np.diff(x) > 1e-9]]

    t0, t1 = x[:-1], x[1:]
    m = 0.5 * (t0 + t1)

    i  = _find_cells(loH, hiH, m)
    j2 = _find_cells(loC, hiC, m)                # sheet2 : hot(θ) | cold(θ)
    j1 = _find_cells(loC, hiC, m + shift_c)      # sheet1 : hot(θ) | cold(θ ± 2π)

    W = np.zeros((nc, nc))

    ok = (i >= 0) & (j2 >= 0)
    dA = H * (F(r("sheet2", t1[ok])) - F(r("sheet2", t0[ok])))
    np.add.at(W, (i[ok], j2[ok]), dA)

    ok = (i >= 0) & (j1 >= 0)
    dA = H * (F(r("sheet1", t1[ok] + shift_s1)) - F(r("sheet1", t0[ok] + shift_s1)))
    np.add.at(W, (i[ok], j1[ok]), dA)

    return W / W.sum()

#%%

import numpy as np

def determine_cell_mean_curvature_radius_SPHE(hx, channel="hot", flow_from="inner", n_sub=10):
    """
    Calcule le rayon de courbure moyen et le rapport de courbure delta (Dh / (2 * rho))
    de chaque cellule de discrétisation pour un échangeur spiralé (SPHE).

    Paramètres:
    - hx : objet échangeur contenant .params et .w
    - channel : 'hot' (canal chaud) ou 'cold' (canal froid)
    - flow_from : 'inner' ou 'outer'
    - n_sub : nombre de points de la méthode des trapèzes par cellule (10 par défaut)

    Retourne:
    - delta_cells : tableau numpy des rapports de courbure par cellule (sans dimension)
    - rho_mean : tableau numpy des rayons de courbure moyens par cellule [m]

    Version vectorisée : toutes les cellules sont intégrées en une fois.
    """
    if getattr(hx, 'HTX_Type', None) not in ('SPHE', 'SCHE'):
        raise ValueError(f"determine_cell_mean_curvature_radius_SPHE : HTX_Type = {getattr(hx, 'HTX_Type', None)!r}, attendu 'SPHE'.")

    P = hx.params
    H = P['H']
    a = (P['b_C'] + P['b_H'] + 2 * P['t']) / (2 * np.pi)
    r0 = P['r0_in']

    # --- sélection du canal -----------------------------------------------------
    ch = channel.lower()
    if ch in ('hot', 'chaud'):
        el = "hot"
        canal_t = P['b_H']
        Dh = P.get('Dh_H', 2 * H * canal_t / (H + canal_t))
    elif ch in ('cold', 'froid'):
        el = "cold"
        canal_t = P['b_C']
        Dh = P.get('Dh_C', 2 * H * canal_t / (H + canal_t))
    else:
        raise ValueError("Le paramètre 'channel' doit être 'hot' ou 'cold'.")

    lo, hi = P['radial_limits'][el]
    r_off = r0 + 0.5 * (lo + hi)                 # rayon au centre du canal à θ = 0

    def r(th):
        return r_off + a * th

    def F(rr):
        return 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a

    # --- vecteur w -------------------------------------------------------------
    w = getattr(hx, 'w', None)
    if w is None:
        raise ValueError("w absent (argument ou hx.w).")
    w = np.asarray(w, float).ravel()
    if abs(w.sum() - 1.0) > 1e-6:
        w = w / w.sum()

    # --- frontières des cellules du canal choisi uniquement ---------------------
    th_end = P['theta_end'][el]
    F0 = F(r(0.0))
    L = float(F(r(th_end)) - F0)
    frac = np.r_[0.0, np.cumsum(w)]
    frac[-1] = 1.0
    if flow_from == "outer":
        frac = 1.0 - frac
    s_k = frac * L
    tk = s_k / L * th_end
    for _ in range(6):                           # Newton sur s(θ)
        tk = tk - (F(r(tk)) - F0 - s_k) / np.sqrt(r(tk)**2 + a**2)
    tk = np.clip(tk, 0.0, th_end)

    # --- rayon de courbure de la spirale d'Archimède ---------------------------
    def curvature_radius(th):
        rr = r(th)
        return (rr**2 + a**2)**1.5 / (rr**2 + 2 * a**2)

    t0, t1 = tk[:-1], tk[1:]
    dt = t1 - t0

    # Grille (cellules x n_sub) et méthode des trapèzes à pas constant :
    # moyenne = (somme des points - (premier + dernier) / 2) / (n_sub - 1)
    theta_sub = t0[:, None] + dt[:, None] * np.linspace(0.0, 1.0, n_sub)
    rho_sub = curvature_radius(theta_sub)
    rho_mean = (rho_sub.sum(axis=1) - 0.5 * (rho_sub[:, 0] + rho_sub[:, -1])) / (n_sub - 1)

    # Cellules de largeur nulle : rayon au point milieu
    degen = np.abs(dt) < 1e-12
    if degen.any():
        rho_mean[degen] = curvature_radius(0.5 * (t0[degen] + t1[degen]))

    delta_cells = Dh / (2.0 * rho_mean)
    return delta_cells, rho_mean

#%%

if __name__ == "__main__":
    # ===============================
    # Example usage
    # ===============================
    
    case_ex = 'Wron'

    shell_side = 'C'

    if case_ex == "Easy":
        
        w_cells = np.array(
            [0.51267984, 0.2861556 , 0.20116456])
        
        n_disc = 3
        n_passes = 4
    
    elif case_ex == "Wrong":
        
        w_cells= np.array([6.65073121, 0.65695199, 0.36510396, 0.20022861, 0.13874025, 0.10644548,
         0.0865282,  0.07301715, 0.06325175, 0.05586614])
    
        n_disc = 10
        n_passes = 2

    elif case_ex == "0.5":
        
        w_cells= np.array([0.5])
    
        n_disc = 2
        n_passes = 1

    else:
        
        w_cells = np.array([
            0.19295366, 0.15062147, 0.12378974, 0.105168,
            0.09147944, 0.08099766, 0.07271404, 0.06600607,
            0.06046302, 0.0558069
        ])
        
        n_disc = 10
        n_passes = 2
    
    over_matrix = determine_cell_overlap_STHE(w_cells, n_passes, shell_side)
    
    import matplotlib.pyplot as plt
    
    # Plot heat map
    plt.figure(figsize=(8,6))
    im = plt.imshow(over_matrix, origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(im, label='Overlap fraction')
    plt.xlabel('Shell cell index')
    plt.ylabel('Tube cell index')
    plt.title(f'Heat Map of Tube-Shell Cell Overlaps ({n_passes}-Pass HX)')
    
    # Add text labels
    n_tube, n_shell = over_matrix.shape
    for i in range(n_tube):
        for j in range(n_shell):
            val = over_matrix[i, j]
            plt.text(j, i, f'{val:.3f}', ha='center', va='center',
                      color='white' if val > over_matrix.max()/2 else 'black')
    
    plt.tight_layout()
    plt.show()
    
    print(np.sum(over_matrix))
    