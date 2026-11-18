# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 10:03:56 2026

@author: Basile
"""

import numpy as np
from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import determine_cell_overlap_STHE
from scipy.optimize import brentq

def compute_cell_x_intersection(cell1, cell2):
    """
    Compute the intersection of two x-intervals.

    Parameters
    ----------
    cell1 : array-like of length 2
        First x-interval [x1, x2]
    cell2 : array-like of length 2
        Second x-interval [x3, x4]

    Returns
    -------
    np.ndarray or None
        [x_start, x_end] if intersection exists, otherwise None
    """
    a1, a2 = np.sort(cell1)
    b1, b2 = np.sort(cell2)

    x_start = max(a1, b1)
    x_end = min(a2, b2)

    if x_start <= x_end:
        return np.array([x_start, x_end])
    else:
        return None

def compute_cell_x_position(w_start, w_end, n_passes, sum_w):
    """
    Compute tube cell positions mapped to 0-1 coordinates per pass.
    Automatically decomposes the cell if it crosses pass turns.

    Parameters
    ----------
    w_start, w_end : float
        Global relative coordinates (0 → 1), w_start < w_end
    n_passes : int
        Number of tube passes

    Returns
    -------
    cells : np.ndarray, shape (N, 2)
        x_start, x_end of sub-cells in pass-local coordinates (0→1 forward, 1→0 backward)
    """
    if w_start >= w_end:
        raise ValueError("w_start must be < w_end")

    # Pass turn locations in global w
    turns = np.linspace(0.0, sum_w, int(n_passes) + 1)

    # Split the cell at pass boundaries
    cuts = [w_start]
    cuts += [t for t in turns[1:-1] if w_start < t < w_end]
    cuts.append(w_end)

    cells = []

    for ws, we in zip(cuts[:-1], cuts[1:]):
        # Determine pass index
        p = np.searchsorted(turns, ws, side='right') - 1
        p = min(p, n_passes - 1)

        # Local coordinate in this pass
        local_ws = (ws - turns[p]) / (turns[p + 1] - turns[p])
        local_we = (we - turns[p]) / (turns[p + 1] - turns[p])

        # Map to pass-local 0-1 coordinates
        if p % 2 == 0:
            # forward pass
            x_start, x_end = local_ws, local_we
        else:
            # backward pass
            x_start, x_end = 1 - local_ws, 1 - local_we

        cells.append([x_start, x_end])

    return np.array(cells)

def interpolate_temperature_xrel(T_start, T_end, x_start, x_end, x_point):
    """
    Linearly interpolate temperature along x_rel (works for forward or backward).
    """
    if x_end == x_start:
        return T_start
    
    if x_end < x_start:
        save = x_end
        x_end = x_start
        x_start = save
    
    fraction = (x_point - x_start) / (x_end - x_start)
    return T_start + fraction * (T_end - T_start)

def compute_LMTD_STHE(HX, debug = 0):
    
    Tj_int_couples = []
    Tk_int_couples = []
        
    w_LMTD = HX.w
    
    # if sum(HX.w) <= 1:
    #     w_LMTD = HX.w
    # else:
    #     w_LMTD = HX.w/sum(HX.w)

    # overlap_matrix = HX.overlap_matrix
    
    w_sum = sum(w_LMTD)
    # w_sum = sum(HX.w)
    
    # overlap_matrix = determine_cell_overlap(w_LMTD, HX.params['Tube_pass'], HX.params['Shell_Side'])
    overlap_matrix = HX.overlap_matrix
    
    LMTD_matrix = np.zeros([len(HX.hvec_c)-1, len(HX.hvec_h)-1])
    
    w = w_LMTD
    w_cum = np.insert(np.cumsum(w),0,0)
    
    if debug:
        print(f"overlap_LMTD : {overlap_matrix}")
        print(f"w_LMTD : {w_LMTD}")
        print(f"w_cum : {w_cum}")

    for k in range(len(HX.hvec_c)-1):
        LMTD_equiv = 0
        UA_tot = 0
        
        contact_j = np.where(overlap_matrix[k, :] > 0.0)[0]
        
        if debug:
            print("-------------------------")
                
            print(f"Tube cell : {k}")
            print(f"contacts shell cells : {contact_j}")
        
        if HX.params['Shell_Side'] == 'H':
            Tki = HX.Tvec_c[k]
            Tko = HX.Tvec_c[k+1]
        else:
            Tki = HX.Tvec_h[k]
            Tko = HX.Tvec_h[k+1]
        
        
        for j in contact_j:
            
            if debug:
                print("----")
                
                print(f"Contact with shell cell {j}")
            
            if HX.params['Shell_Side'] == 'H':
                Tji = HX.Tvec_h[j+1]
                Tjo = HX.Tvec_h[j]
            else:
                Tji = HX.Tvec_c[j+1]
                Tjo = HX.Tvec_c[j]
                
            "1) Determine spatial intersections"
            intersections = []
            x_rel_k = []
            
            w_start_k = w_cum[k]
            w_end_k = w_cum[k+1]
            
            w_start_j = w_cum[j]
            w_end_j = w_cum[j+1]
            
            x_rel_j = np.array([w_start_j, w_end_j])/w_sum
            x_rel_k = compute_cell_x_position(w_start_k, w_end_k, HX.params['Tube_pass'], w_sum)
            
            if debug:
                print(f"x_rel_k : {x_rel_k}")
                print(f"x_rel_j : {x_rel_j}")
            
            valid_intersections = []
            valid_x_rel_k = []
            
            for x_relk in x_rel_k:
                inter = compute_cell_x_intersection(x_relk, x_rel_j)
                if inter is not None:
                    valid_intersections.append(inter)
                    valid_x_rel_k.append(x_relk)
            
            # Then use these in the loop
            intersections = valid_intersections
            x_rel_k = valid_x_rel_k

            if debug:
                print(f"intersections : {intersections}")
                
            "2) Get the temperature values of the intersections" 
            
            if debug:
                print(f"Tki : {Tki}, Tko : {Tko}")
                print(f"Tji : {Tji}, Tjo : {Tjo}")
            
            # Linearly interpolate each temperature inside its cell at the intersection location
            
            for i in range(len(intersections)):
                intersection = intersections[i]
                x_relk = x_rel_k[i]
                
                x_int_start, x_int_end = intersection
                
                if debug:
                    print(f"----")
                    print(f"x_relk : {x_relk}")
                    print(f"x_rel_j : {x_rel_j}")
    
                    print(f"x_int_start : {x_int_start}, x_int_end : {x_int_end}")
                
                # Tube temperatures at intersection endpoints
                Tki_int = interpolate_temperature_xrel(Tki, Tko, x_relk[0], x_relk[1], x_int_start)
                Tko_int = interpolate_temperature_xrel(Tki, Tko, x_relk[0], x_relk[1], x_int_end)
            
                # Shell temperatures at intersection endpoints
                Tji_int = interpolate_temperature_xrel(Tjo, Tji, x_rel_j[0], x_rel_j[1], x_int_start)
                Tjo_int = interpolate_temperature_xrel(Tjo, Tji, x_rel_j[0], x_rel_j[1], x_int_end)
                
                if debug:
                    print(f"Tki : {Tki_int}, Tko : {Tko_int}")
                    print(f"Tji : {Tji_int}, Tjo : {Tjo_int}")
                
                if HX.params['Shell_Side'] == 'H':
                    Thi = Tjo_int
                    Tho = Tji_int
                    
                    Tci = Tki_int
                    Tco = Tko_int
                else:
                    Tci = Tjo_int
                    Tco = Tji_int
                    
                    Thi = Tki_int
                    Tho = Tko_int
                   
                deltaT1 = max(Thi - Tco, 1e-6)
                deltaT2 = max(Tho - Tci, 1e-6)
                
                # print(f"deltaT1 : {deltaT1}")
                # print(f"deltaT2 : {deltaT2}")
                
                # Avoid division by zero
                if deltaT1 == deltaT2:
                    LMTD_inter = deltaT1  # ΔT constant → LMTD = ΔT
                else:
                    LMTD_inter =  (deltaT1 - deltaT2) / np.log(deltaT1 / deltaT2)
                
                LMTD_equiv += LMTD_inter*HX.UA_matrix[k][j]
                UA_tot += HX.UA_matrix[k][j]
                
                if debug:    
                    print(f"LMTD_inter : {LMTD_inter}")
                    print(f"UA_weight : {HX.UA_matrix[j][k]}")
                    print(f"UA_tot : {UA_tot}")
                
                # if LMTD_equiv > 0:
                Tj_int_couples.append([Tji_int, Tjo_int])
                Tk_int_couples.append([Tki_int, Tko_int])
                
            if UA_tot != 0:
                LMTD_matrix[k][j] = LMTD_equiv/UA_tot
            else:
                LMTD_matrix[k][j] = 0
                # print(f"LMTD : {LMTD_equiv/UA_tot}")

    # Compute LMTD_vector
    LMTD_vector = np.zeros(len(HX.hvec_c)-1)

    for k in range(len(HX.hvec_c)-1):
        UA_tot = 0
        
        for j in range(len(HX.hvec_h)-1):            
            LMTD_vector[k] += LMTD_matrix[k][j]*HX.UA_matrix[k][j]
            UA_tot += HX.UA_matrix[k][j]
            
        if UA_tot != 0:
            LMTD_vector[k] = LMTD_vector[k]/UA_tot
        else:
            LMTD_vector[k] = 0
            
    if debug:
        print(f"--------------------")
        print(f"--------------------")
    
        import matplotlib.pyplot as plt
        Tj_unique = np.unique(Tj_int_couples, axis = 0)
        Tk_unique = np.unique(Tk_int_couples, axis = 0)
    
        print("----------------------")
    
        print("Hot")
        
        print(Tj_unique)
        
        print("Cold")
    
        print(Tk_unique)
        
        x = [0]
        Th = [Tj_unique[0, 0]]
        Tc = [Tk_unique[0, 0]]
        
        for (Thi, Tho), (Tci, Tco) in zip(Tj_unique, Tk_unique):
            x.append(x[-1] + 1)
            Th.append(Tho)
            Tc.append(Tco)
        
        x = np.array(x)
        Th = np.array(Th)
        Tc = np.array(Tc)
        
        plt.figure()
        plt.plot(x, Th, '-r', label='Hot side')
        plt.plot(x, Tc, '-b', label='Cold side')
        
        plt.xlabel("Segment index (cumulative)")
        plt.ylabel("Temperature")
        plt.title("Sequential temperature profiles")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        print(f"--------------------")
        print(f"--------------------")
    
    return LMTD_matrix, LMTD_vector

#%%

def sphe_segments(thH, thC, th_h, th_c):
    """Segments (i, j, shift, t0, t1) en angle du canal chaud ; thH, thC croissants (cellule 0 = noyau)."""
    TWO_PI = 2*np.pi
    x = np.unique(np.r_[0.0, thH, thC, thC - TWO_PI, th_h])
    x = x[(x >= 0) & (x <= th_h)]
    x = x[np.r_[True, np.diff(x) > 1e-9]]
    segs = []
    for t0, t1 in zip(x[:-1], x[1:]):
        m = 0.5*(t0 + t1)
        i = np.searchsorted(thH, m) - 1
        for shift in (0.0, TWO_PI):            # sheet2 / sheet1
            if m + shift < th_c:
                j = np.searchsorted(thC, m + shift) - 1
                segs.append((i, j, shift, t0, t1))
    return segs

TWO_PI = 2 * np.pi

def compute_LMTD_SPHE(HX, U=None, debug=False):
    """
    LMTD par paire (i chaud, j froid) pour un échangeur spirale, à partir de HX seul.

    Lit : HX.params, HX.w, HX.Tvec_h, HX.Tvec_c (n+1 valeurs chacun)
    Convention supposée : cellule 0 côté noyau pour les deux fluides ; Tvec_x[k] = T à la borne
    intérieure (θ le plus petit) de la cellule k, Tvec_x[k+1] = T à sa borne extérieure.
    U : None -> on répartit HX.UA_matrix[i, j] (i chaud, j froid) au prorata de l'aire des segments
        scalaire ou callable (i, j) -> U : UA = U · aire du segment
    Retourne LMTD_matrix[i, j], LMTD_hot[i], LMTD_cold[j], UA[i, j]

    Version vectorisée : mêmes résultats et mêmes messages que la version à boucles.
    """
    P = HX.params
    H = P['H']
    a = (P['C_canal_t'] + P['H_canal_t'] + 2 * P['t']) / TWO_PI
    r0 = P['r0_in']
    RL = P['radial_limits']
    th_h, th_c = P['theta_end']['hot'], P['theta_end']['cold']

    w = np.asarray(HX.w, float).ravel()
    w = w / w.sum()
    n = len(w)

    Th = np.asarray(HX.Tvec_h, float).ravel()
    Tc = np.asarray(HX.Tvec_c, float).ravel()
    if len(Th) != n + 1 or len(Tc) != n + 1:
        raise ValueError(f"Tvec_h/Tvec_c : {len(Th)}/{len(Tc)} valeurs, attendu {n+1}")

    # --- géométrie ---------------------------------------------------------
    def r(el, th, pos=0.5):
        lo, hi = RL[el]
        return r0 + a * th + lo + pos * (hi - lo)

    def F(rr):
        return 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a

    frac = np.r_[0.0, np.cumsum(w)]
    frac[-1] = 1.0

    def boundaries(el):
        th_end = P['theta_end'][el]
        F0 = F(r(el, 0.0))
        L = float(F(r(el, th_end)) - F0)
        s_k = frac * L
        th = s_k / L * th_end
        for _ in range(6):                                   # Newton sur s(θ)
            th = th - (F(r(el, th)) - F0 - s_k) / np.sqrt(r(el, th)**2 + a**2)
        return np.clip(th, 0.0, th_end)

    thH, thC = boundaries('hot'), boundaries('cold')

    # sheet1 orientation depends on which channel is inside
    if P['inner_channel'] == "cold":     # hot outer: hot(θ) | sheet1(θ+2π) | cold(θ+2π)
        shift_c, shift_s1 = TWO_PI, TWO_PI
    else:                                # hot inner: cold(θ-2π) | sheet1(θ) | hot(θ)
        shift_c, shift_s1 = -TWO_PI, 0.0

    x = np.unique(np.r_[0.0, thH, thC, thC - shift_c, th_h])
    x = x[(x >= 0) & (x <= th_h)]
    x = x[np.r_[True, np.diff(x) > 1e-9]]

    T0s, T1s = x[:-1], x[1:]
    m = 0.5 * (T0s + T1s)
    i_seg = np.clip(np.searchsorted(thH, m) - 1, 0, n - 1)
    k_seg = np.arange(len(m))

    I, J, SH, T0, T1, K = [], [], [], [], [], []
    for shift in (0.0, shift_c):
        mc = m + shift
        ok = (mc > 0) & (mc < th_c)          # lower bound needed for the -2π case
        I.append(i_seg[ok])
        J.append(np.clip(np.searchsorted(thC, mc[ok]) - 1, 0, n - 1))
        SH.append(np.full(ok.sum(), shift))
        T0.append(T0s[ok]); T1.append(T1s[ok]); K.append(k_seg[ok])
        
    I, J, SH, T0, T1, K = map(np.concatenate, (I, J, SH, T0, T1, K))
    order = np.lexsort((SH, K))                 # même ordre que la boucle d'origine
    i, j, sh, t0, t1 = I[order], J[order], SH[order], T0[order], T1[order]

    # --- aires des segments ------------------------------------------------
    s2 = sh == 0.0
    A_seg = np.empty(len(i))
    A_seg[s2]  = H * (F(r('sheet2', t1[s2]))  - F(r('sheet2', t0[s2])))
    A_seg[~s2] = H * (F(r('sheet1', t1[~s2] + shift_s1)) - F(r('sheet1', t0[~s2] + shift_s1)))
    
    A_pair = np.zeros((n, n))
    np.add.at(A_pair, (i, j), A_seg)

    # --- ΔT aux extrémités des segments -------------------------------------
    def Tint(nodes, edges, k, th):              # interpolation sur les bornes COMPLÈTES de la cellule
        f = (th - edges[k]) / (edges[k + 1] - edges[k])
        return nodes[k] + f * (nodes[k + 1] - nodes[k])

    dT0 = Tint(Th, thH, i, t0) - Tint(Tc, thC, j, t0 + sh)
    dT1 = Tint(Th, thH, i, t1) - Tint(Tc, thC, j, t1 + sh)

    bad = (dT0 <= 0) | (dT1 <= 0)
    warn = {k: f"⚠ croisement de températures, segment ({i[k]},{j[k]}), ΔT = {dT0[k]:.3g}/{dT1[k]:.3g}"
            for k in np.nonzero(bad)[0]}
    
    # if not debug:
    #     for msg in warn.values():
    #         print(msg)
            
    dT0 = np.where(bad, np.maximum(dT0, 1e-6), dT0)
    dT1 = np.where(bad, np.maximum(dT1, 1e-6), dT1)

    with np.errstate(divide='ignore', invalid='ignore'):
        same = np.abs(dT0 - dT1) < 1e-9 * np.maximum(dT0, dT1)
        lm = np.where(same, dT0, (dT0 - dT1) / np.log(dT0 / dT1))

    # --- UA par segment -----------------------------------------------------
    if U is None:
        UA_user = getattr(HX, 'UA_matrix', None)
        if UA_user is None:
            raise ValueError("U non fourni et HX.UA_matrix absent.")
        UA_user = np.asarray(UA_user, float)      # supposé [i chaud, j froid]
        Ap = A_pair[i, j]
        ua = np.divide(UA_user[i, j] * A_seg, Ap, out=np.zeros_like(A_seg), where=Ap > 0)
    elif callable(U):
        # un appel par paire (i, j) distincte au lieu d'un appel par segment
        pairs, inv = np.unique(np.stack((i, j), axis=1), axis=0, return_inverse=True)
        U_pair = np.array([U(int(p), int(q)) for p, q in pairs], float)
        ua = U_pair[inv.ravel()] * A_seg
    else:
        ua = U * A_seg

    UA = np.zeros((n, n))
    QL = np.zeros((n, n))
    np.add.at(UA, (i, j), ua)
    np.add.at(QL, (i, j), ua * lm)

    if debug:                                   # même ordre que l'original : avertissement puis ligne debug
        for k in range(len(i)):
            if k in warn:
                print(warn[k])
            print(f"({i[k]},{j[k]}) shift={sh[k]:.2f} θ=[{t0[k]:.3f},{t1[k]:.3f}] "
                  f"ΔT={dT0[k]:.2f}/{dT1[k]:.2f} LMTD={lm[k]:.2f}")

    LMTD_matrix = np.divide(QL, UA, out=np.zeros_like(QL), where=UA > 0)
    s1, s0 = UA.sum(1), UA.sum(0)
    LMTD_hot = np.divide(QL.sum(1), s1, out=np.zeros(n), where=s1 > 0)
    LMTD_cold = np.divide(QL.sum(0), s0, out=np.zeros(n), where=s0 > 0)
    return LMTD_matrix, LMTD_hot, LMTD_cold, UA
 
#%% ---------------------------------------------------------------------------
# Géométrie
 
def _sphe_geometry_arrays(HX, w):
    """
    Cœur vectorisé de sphe_geometry.
    Retourne thH, thC et les segments sous forme de tableaux
    (i, j, shift, t0, t1), dans le même ordre que la version liste :
    segment par segment, shift = 0 puis shift = 2π.
    """
    P = HX.params
    a = (P['C_canal_t'] + P['H_canal_t'] + 2 * P['t']) / TWO_PI
    n = len(w)
 
    def r(el, th):
        lo, hi = P['radial_limits'][el]
        return P['r0_in'] + a * np.asarray(th) + 0.5 * (lo + hi)
 
    def F(rr):
        return 0.5 * (rr * np.sqrt(rr**2 + a**2) + a**2 * np.arcsinh(rr / a)) / a
 
    def boundaries(el):
        th_end = P['theta_end'][el]
        F0 = F(r(el, 0.0))
        s_of = lambda th: F(r(el, th)) - F0
        L = float(s_of(th_end))
        frac = np.r_[0.0, np.cumsum(w)]
        frac[-1] = 1.0
        s_k = frac * L
        th = s_k / L * th_end
        for _ in range(6):                                   # Newton sur s(θ)
            th = th - (s_of(th) - s_k) / np.sqrt(r(el, th)**2 + a**2)
        return np.clip(th, 0.0, th_end)
 
    thH, thC = boundaries('hot'), boundaries('cold')
    th_h, th_c = P['theta_end']['hot'], P['theta_end']['cold']
 
    x = np.unique(np.r_[0.0, thH, thC, thC - TWO_PI, th_h])
    x = x[(x >= 0) & (x <= th_h)]
    x = x[np.r_[True, np.diff(x) > 1e-9]]
 
    t0, t1 = x[:-1], x[1:]
    m = 0.5 * (t0 + t1)
    i = np.clip(np.searchsorted(thH, m) - 1, 0, n - 1)
    seg = np.arange(len(m))
 
    I, J, SH, T0, T1, K = [], [], [], [], [], []
    for shift in (0.0, TWO_PI):
        ok = (m + shift) < th_c
        j = np.clip(np.searchsorted(thC, m[ok] + shift) - 1, 0, n - 1)
        I.append(i[ok]); J.append(j); SH.append(np.full(ok.sum(), shift))
        T0.append(t0[ok]); T1.append(t1[ok]); K.append(seg[ok])
 
    I, J, SH, T0, T1, K = map(np.concatenate, (I, J, SH, T0, T1, K))
    order = np.lexsort((SH, K))          # même ordre que la boucle d'origine
 
    return thH, thC, I[order], J[order], SH[order], T0[order], T1[order]
 
 
def sphe_geometry(HX, w):
    """
    Frontières angulaires des cellules (chaud, froid) et segments de contact.
 
    w : parts de longueur des cellules (somme = 1), identiques chaud/froid.
    Retourne thH, thC (n+1 valeurs croissantes) et segs = liste de
    (i chaud, j froid, shift, t0, t1), angles exprimés dans le repère du canal chaud
    (shift = 0 : contact via sheet2 ; shift = 2π : contact via sheet1, froid un tour plus loin).
 
    Interface inchangée (liste de tuples) pour les fonctions de tracé.
    """
    thH, thC, i, j, sh, t0, t1 = _sphe_geometry_arrays(HX, np.asarray(w, float))
    segs = [(int(a), int(b), float(c), float(d), float(e))
            for a, b, c, d, e in zip(i, j, sh, t0, t1)]
    return thH, thC, segs
 
 
#%% ---------------------------------------------------------------------------
# Cache de géométrie
 
def _freeze(v):
    """Rend hashable les dicts/listes/arrays de paramètres."""
    if isinstance(v, dict):
        return tuple(sorted((k, _freeze(x)) for k, x in v.items()))
    if isinstance(v, (list, tuple, np.ndarray)):
        return tuple(_freeze(x) for x in np.asarray(v, dtype=object).ravel())
    return v
 
 
def _sphe_geometry_cached(HX, w):
    """
    Cache de la géométrie par (paramètres géométriques, w), stocké sur l'objet HX.
    La clé couvre tout ce que lit _sphe_geometry_arrays : changer un de ces
    paramètres sur le même HX invalide automatiquement l'entrée.
    """
    P = HX.params
    key = (
        P['C_canal_t'], P['H_canal_t'], P['t'], P['r0_in'],
        _freeze(P['radial_limits']), _freeze(P['theta_end']),
        len(w), w.tobytes(),
    )
    cache = HX.__dict__.setdefault('_sphe_geom_cache', {})
    if key not in cache:
        cache[key] = _sphe_geometry_arrays(HX, w)
    return cache[key]
 
 
#%% ---------------------------------------------------------------------------
# ΔT minimal pour une puissance Q
 
def sphe_min_dT(HX, Q, pvec_c=None, pvec_h=None, w_guess=None):
    """
    ΔT minimal (T_chaud - T_froid) sur tous les segments de contact pour une puissance Q.
    Met à jour les vecteurs de HX via calculate_cell_boundaries(Q).
    pvec_c / pvec_h : conservés pour compatibilité d'appel (non utilisés).
    Retourne (dT_min, (i, j, shift, θ)).
    """
    HX.calculate_cell_boundaries(Q)
 
    Th = np.asarray(HX.Tvec_h, float)
    Tc = np.asarray(HX.Tvec_c, float)
    n = len(Th) - 1
 
    if w_guess is not None and len(w_guess) == n:
        w = np.asarray(w_guess, float)
        w = w / w.sum()
    else:
        w = np.full(n, 1.0 / n)
 
    thH, thC, i, j, shift, t0, t1 = _sphe_geometry_cached(HX, w)
 
    # ΔT linéaire par morceaux : min aux extrémités. Ordre entrelacé
    # (seg0 t0, seg0 t1, seg1 t0, ...) comme la boucle d'origine.
    th  = np.column_stack((t0, t1)).ravel()
    ii  = np.repeat(i, 2)
    jj  = np.repeat(j, 2)
    sh  = np.repeat(shift, 2)
    thc = th + sh
 
    # Interpolation sur les bornes COMPLÈTES de chaque cellule
    fH = (th  - thH[ii]) / (thH[ii + 1] - thH[ii])
    fC = (thc - thC[jj]) / (thC[jj + 1] - thC[jj])
 
    dT = (Th[ii] + fH * (Th[ii + 1] - Th[ii])) - (Tc[jj] + fC * (Tc[jj + 1] - Tc[jj]))
 
    k = int(np.argmin(dT))
    return float(dT[k]), (int(ii[k]), int(jj[k]), float(sh[k]), float(th[k]))
 
 
#%% ---------------------------------------------------------------------------
# Pincement interne
 
def internal_pinching_SPHE(HX, stream='hot', pvec_c=None, pvec_h=None,
                           Q_hi=None, dT_min=0.0, w_guess=None,
                           rtol=1e-5, max_iter=50):
    """
    Qmax géométrique : plus grand Q tel que ΔT >= dT_min sur tous les segments de contact.
 
    stream  : 'hot' fait le calcul ; 'cold' renvoie None (le pincement ne dépend pas du fluide,
              ce qui évite le doublon dans solve_hx).
    Q_hi    : borne haute (défaut HX.Qmax_ext).
    w_guess : positions des cellules (défaut uniformes ; HX.w de l'itération précédente
              donne un pincement plus fidèle).
    rtol    : précision relative sur Qmax (rapportée à Q_hi).
    Retourne Qmax (= Q_hi si le pincement géométrique n'est pas actif avant Q_hi),
    ou None si stream == 'cold' ou si les deux fluides sont supercritiques.
    Effets : HX.Qmax_geo, HX.pinch_dT, HX.pinch_where renseignés ; HX positionné sur Qmax.
    """
    if stream not in ('hot', 'cold'):
        raise ValueError(stream)
    if stream == 'cold':
        return None
    if getattr(HX, 'SC_h', False) and getattr(HX, 'SC_c', False):
        return None
 
    if Q_hi is None:
        Q_hi = getattr(HX, 'Qmax_ext', None)
        if Q_hi is None:
            raise ValueError("Q_hi absent et HX.Qmax_ext non défini (lancer external_pinching avant).")
 
    # ---- évaluations mises en cache ----
    cache = {}
    last = [None]                         # dernier Q évalué = état courant de HX
 
    def eval_dT(Q):
        if Q not in cache:
            cache[Q] = sphe_min_dT(HX, Q, pvec_c, pvec_h, w_guess)
            last[0] = Q
        return cache[Q]
 
    def f(Q):
        return eval_dT(Q)[0] - dT_min
 
    # ---- bornes ----
    Q_lo = 1.0
    if f(Q_lo) < 0:
        raise ValueError("ΔT < dT_min même à Q = 0 : vérifier T_hi > T_ci.")
 
    if f(Q_hi) >= 0:
        Q_star = Q_hi                     # pas de pincement géométrique avant Q_hi
    else:
        brentq(f, Q_lo, Q_hi, xtol=rtol * Q_hi, maxiter=max_iter)
        # Comme l'ancienne bissection (qui renvoyait lo) : plus grand Q évalué
        # du côté faisable (ΔT >= dT_min).
        Q_star = max(Q for Q, v in cache.items() if v[0] - dT_min >= 0)
 
    # ---- repositionner HX sur Q_star seulement si nécessaire ----
    if last[0] != Q_star:
        cache.pop(Q_star)
        eval_dT(Q_star)
 
    HX.pinch_dT, HX.pinch_where = cache[Q_star]
    HX.Qmax_geo = Q_star
    return Q_star
 

