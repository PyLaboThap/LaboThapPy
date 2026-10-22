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

def determine_cell_overlap_SPHE(hx, flow_from="inner"):
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
    w = getattr(hx, 'w', None)
    if w is None:
        raise ValueError("plot_contacts : w absent (argument ou hx.w).")
    w = np.asarray(w, float).ravel()
    if abs(w.sum() - 1.0) > 1e-6:
        # print(f"⚠ somme(w) = {w.sum():.6f} ≠ 1 : w normalisé")
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
    tkH = boundaries("hot", flow_from)
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
    
    return W

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
    
    over_matrix = determine_cell_overlap(w_cells, n_passes, shell_side)
    
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
    