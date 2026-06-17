import numpy as np

def determine_cell_overlap(w_cells, n_passes, shell_side):
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
    