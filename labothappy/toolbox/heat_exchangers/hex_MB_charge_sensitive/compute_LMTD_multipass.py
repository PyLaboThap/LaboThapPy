# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 10:03:56 2026

@author: Basile
"""

import numpy as np
from toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import determine_cell_overlap

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

def determine_LMTD_multipass(HX):
    
    Tj_int_couples = []
    Tk_int_couples = []
    
    debug = 0
    
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
        
        # Tcell_k = HX.Tvec_c
        
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

# -*- coding: utf-8 -*-
# """
# Created on Sat Jan 10 10:03:56 2026

# @author: Basile
# """

# import numpy as np
# from toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import determine_cell_overlap

# def compute_cell_x_intersection(cell1, cell2):
#     a1, a2 = cell1 if cell1[0] <= cell1[1] else (cell1[1], cell1[0])
#     b1, b2 = cell2 if cell2[0] <= cell2[1] else (cell2[1], cell2[0])

#     x_start = max(a1, b1)
#     x_end   = min(a2, b2)
#     return None if x_start > x_end else np.array([x_start, x_end])

# def compute_cell_x_intersection(cell1, cell2):
#     x_start = max(cell1[0], cell2[0])
#     x_end   = min(cell1[1], cell2[1])
#     return None if x_start > x_end else np.array([x_start, x_end])

# def compute_cell_x_position(w_start, w_end, n_passes, sum_w):
#     """
#     Compute tube cell positions mapped to 0-1 coordinates per pass.
#     Automatically decomposes the cell if it crosses pass turns.

#     Parameters
#     ----------
#     w_start, w_end : float
#         Global relative coordinates (0 → 1), w_start < w_end
#     n_passes : int
#         Number of tube passes

#     Returns
#     -------
#     cells : np.ndarray, shape (N, 2)
#         x_start, x_end of sub-cells in pass-local coordinates (0→1 forward, 1→0 backward)
#     """
#     if w_start >= w_end:
#         raise ValueError("w_start must be < w_end")

#     # Pass turn locations in global w
#     turns = np.linspace(0.0, sum_w, n_passes + 1)

#     # Split the cell at pass boundaries
#     cuts = [w_start]
#     cuts += [t for t in turns[1:-1] if w_start < t < w_end]
#     cuts.append(w_end)

#     cells = []

#     for ws, we in zip(cuts[:-1], cuts[1:]):
#         # Determine pass index
#         p = np.searchsorted(turns, ws, side='right') - 1
#         p = min(p, n_passes - 1)

#         # Local coordinate in this pass
#         local_ws = (ws - turns[p]) / (turns[p + 1] - turns[p])
#         local_we = (we - turns[p]) / (turns[p + 1] - turns[p])

#         # Map to pass-local 0-1 coordinates
#         if p % 2 == 0:
#             # forward pass
#             x_start, x_end = local_ws, local_we
#         else:
#             # backward pass
#             x_start, x_end = 1 - local_ws, 1 - local_we

#         cells.append([x_start, x_end])

#     return np.array(cells)

# def compute_cell_x_position(w_start, w_end, n_passes, turns, sum_w):
#     if w_start >= w_end:
#         raise ValueError("w_start must be < w_end")

#     cuts = [w_start]
#     cuts += [t for t in turns[1:-1] if w_start < t < w_end]
#     cuts.append(w_end)

#     cells = []

#     for ws, we in zip(cuts[:-1], cuts[1:]):
#         p = np.searchsorted(turns, ws, side="right") - 1
#         p = min(p, n_passes - 1)

#         denom = turns[p + 1] - turns[p]
#         local_ws = (ws - turns[p]) / denom
#         local_we = (we - turns[p]) / denom

#         if p % 2 == 0:
#             x0, x1 = local_ws, local_we
#         else:
#             x0, x1 = 1 - local_ws, 1 - local_we

#         cells.append((x0, x1) if x0 < x1 else (x1, x0))

#     return np.array(cells, dtype=np.float64)

# def compute_cell_x_position(w_start, w_end, n_passes, turns, sum_w):
#     if w_start >= w_end:
#         raise ValueError("w_start must be < w_end")

#     # Use boolean mask instead of list comprehension
#     turns_inside = turns[1:-1]
#     mask = (turns_inside > w_start) & (turns_inside < w_end)
#     cuts = np.concatenate(([w_start], turns_inside[mask], [w_end])).astype(np.float32)

#     n_cuts = len(cuts) - 1
#     cells = np.empty((n_cuts, 2), dtype=np.float32)

#     for i in range(n_cuts):
#         ws, we = cuts[i], cuts[i + 1]
#         p = np.searchsorted(turns, ws, side="right") - 1
#         p = min(p, n_passes - 1)

#         denom = turns[p + 1] - turns[p]
#         local_ws = (ws - turns[p]) / denom
#         local_we = (we - turns[p]) / denom

#         if p % 2 == 0:
#             x0, x1 = local_ws, local_we
#         else:
#             x0, x1 = 1 - local_ws, 1 - local_we

#         # Assign directly to array
#         cells[i, 0] = min(x0, x1)
#         cells[i, 1] = max(x0, x1)

#     return cells

# def compute_cell_x_position(w_start, w_end, n_passes, turns, sum_w):
#     """
#     Compute relative positions of tube cell within multipass heat exchanger.
#     Fully vectorized using NumPy.
#     Returns array of shape (n_segments, 2) in float32.
#     """
#     if w_start >= w_end:
#         raise ValueError("w_start must be < w_end")

#     # 1) Select turns that fall inside this tube cell
#     turns_inside = turns[1:-1]  # ignore first and last
#     mask = (turns_inside > w_start) & (turns_inside < w_end)
#     cuts = np.concatenate(([w_start], turns_inside[mask], [w_end])).astype(np.float32)

#     # 2) Vectorized differences
#     ws = cuts[:-1]
#     we = cuts[1:]
#     n_cuts = len(ws)

#     # 3) Determine which pass each segment belongs to
#     pass_idx = np.searchsorted(turns, ws, side='right') - 1
#     pass_idx = np.minimum(pass_idx, n_passes - 1)

#     # 4) Compute local positions inside each pass
#     denom = turns[pass_idx + 1] - turns[pass_idx]
#     local_ws = (ws - turns[pass_idx]) / denom
#     local_we = (we - turns[pass_idx]) / denom

#     # 5) Handle reversed passes
#     even_pass = (pass_idx % 2 == 0)
#     x0 = np.where(even_pass, local_ws, 1 - local_ws)
#     x1 = np.where(even_pass, local_we, 1 - local_we)

#     # 6) Ensure x0 < x1
#     cells = np.stack([np.minimum(x0, x1), np.maximum(x0, x1)], axis=1).astype(np.float32)

#     return cells

# def interpolate_temperature_xrel(T_start, T_end, x_start, x_end, x_point):
#     """
#     Linearly interpolate temperature along x_rel (works for forward or backward).
#     """
#     if x_end == x_start:
#         return T_start
    
#     if x_end < x_start:
#         save = x_end
#         x_end = x_start
#         x_start = save
    
#     fraction = (x_point - x_start) / (x_end - x_start)
#     return T_start + fraction * (T_end - T_start)

# def determine_LMTD_multipass(HX, debug=0):
#     """
#     Optimized multipass LMTD computation.
#     Uses float32, precomputed slopes, and minimal loops.
#     """
#     # Preallocate lists for temperature couples (optional)
#     Tj_int_couples = []
#     Tk_int_couples = []

#     w_cum = np.insert(HX.w_cumsum.astype(np.float32), 0, 0.0)
#     w_sum = w_cum[-1]
#     n_passes = int(HX.params['Tube_pass'])
#     turns = np.linspace(0.0, w_sum, n_passes + 1, dtype=np.float32)

#     overlap_matrix = HX.overlap_matrix
#     len_hvec = len(HX.hvec_c) - 1

#     LMTD_matrix = np.zeros((len_hvec, len_hvec), dtype=np.float32)

#     for k in range(len_hvec):
#         LMTD_equiv = np.float32(0.0)
#         UA_tot = np.float32(0.0)

#         contact_j = np.where(overlap_matrix[k, :] > 0.0)[0]

#         # Tube temperatures
#         if HX.params['Shell_Side'] == 'H':
#             Tki, Tko = HX.Tvec_c[k], HX.Tvec_c[k+1]
#         else:
#             Tki, Tko = HX.Tvec_h[k], HX.Tvec_h[k+1]

#         w_start_k, w_end_k = w_cum[k], w_cum[k+1]
#         x_rel_k_all = compute_cell_x_position(w_start_k, w_end_k, n_passes, turns, w_sum).astype(np.float32)

#         # Precompute slope for tube linear interpolation
#         slopes_tube = np.array([(x[1]-x[0], Tko-Tki) for x in x_rel_k_all], dtype=np.float32)

#         for j in contact_j:
#             # Shell temperatures
#             if HX.params['Shell_Side'] == 'H':
#                 Tji, Tjo = HX.Tvec_h[j+1], HX.Tvec_h[j]
#             else:
#                 Tji, Tjo = HX.Tvec_c[j+1], HX.Tvec_c[j]

#             w_start_j, w_end_j = w_cum[j], w_cum[j+1]
#             x_rel_j = np.array([w_start_j, w_end_j], dtype=np.float32) / w_sum

#             # Determine intersections
#             valid_intersections = []
#             valid_x_rel_k = []

#             for x_relk in x_rel_k_all:
#                 inter = compute_cell_x_intersection(x_relk, x_rel_j)
#                 if inter is not None:
#                     valid_intersections.append(inter)
#                     valid_x_rel_k.append(x_relk)

#             intersections = valid_intersections
#             x_rel_k = valid_x_rel_k

#             # Precompute tube slopes for intersections
#             for idx, (intersection, x_relk) in enumerate(zip(intersections, x_rel_k)):
#                 x_int_start, x_int_end = intersection

#                 denom_tube = x_relk[1] - x_relk[0]
#                 slope_tube = (Tko - Tki) / denom_tube if denom_tube != 0 else 0.0
#                 Tki_int = Tki + slope_tube * (x_int_start - x_relk[0])
#                 Tko_int = Tki + slope_tube * (x_int_end - x_relk[0])

#                 denom_shell = x_rel_j[1] - x_rel_j[0]
#                 slope_shell = (Tji - Tjo) / denom_shell if denom_shell != 0 else 0.0
#                 Tjo_int = Tjo + slope_shell * (x_int_start - x_rel_j[0])
#                 Tji_int = Tjo + slope_shell * (x_int_end - x_rel_j[0])

#                 if HX.params['Shell_Side'] == 'H':
#                     Thi, Tho = Tjo_int, Tji_int
#                     Tci, Tco = Tki_int, Tko_int
#                 else:
#                     Tci, Tco = Tjo_int, Tji_int
#                     Thi, Tho = Tki_int, Tko_int

#                 deltaT1 = max(Thi - Tco, 1e-6)
#                 deltaT2 = max(Tho - Tci, 1e-6)

#                 if deltaT1 == deltaT2:
#                     LMTD_inter = deltaT1
#                 else:
#                     LMTD_inter = (deltaT1 - deltaT2) / np.log(deltaT1 / deltaT2)

#                 LMTD_equiv += LMTD_inter * HX.UA_matrix[k, j]
#                 UA_tot += HX.UA_matrix[k, j]

#                 # Optional: save couples
#                 Tj_int_couples.append([Tji_int, Tjo_int])
#                 Tk_int_couples.append([Tki_int, Tko_int])

#             LMTD_matrix[k, j] = LMTD_equiv / UA_tot if UA_tot != 0 else 0.0

#     UA_sum = HX.UA_matrix.sum(axis=1)
#     LMTD_vector = np.divide(
#         (LMTD_matrix * HX.UA_matrix).sum(axis=1),
#         UA_sum,
#         out=np.zeros_like(UA_sum, dtype=np.float32),
#         where=UA_sum != 0
#     )

#     return LMTD_matrix, LMTD_vector