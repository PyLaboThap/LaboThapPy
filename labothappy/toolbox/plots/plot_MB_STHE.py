import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

def plot_MB_STHE(HX, H_total=0.3, L_total=1.0, props=None, legend_range=None, save_path = None):
    """
    Trace le rectangle de la shell et les rectangles représentant les tubes,
    avec les lignes rouges pour les cellules et lignes bleues pour les passes.
    
    Les lignes bleues ne sont tracées qu'aux positions réelles de fin de cellule.
    """

    w = np.array(HX.w)
    w = w/np.sum(w)

    n_passes = HX.params["Tube_pass"]
    
    # Marges et dimensions des tubes
    spacing_frac = 0.3/n_passes
    tube_spacing = H_total * spacing_frac
    tube_height = (H_total - (n_passes + 1) * tube_spacing) / n_passes
    if tube_height <= 0:
        raise ValueError("Pas assez de hauteur pour le nombre de passes demandé")
    
    # Figure
    fig, ax = plt.subplots(figsize=(8, H_total*15))
    
    "1) Rectangle de la shell"
    shell_rect = plt.Rectangle((0,0), L_total, H_total, edgecolor='black', facecolor='lightgray')
    ax.add_patch(shell_rect)
    
    if HX.params['Shell_Side'] == 'H':
        tube_color = 'blue'
        shell_color = 'red'
    else:
        tube_color = 'red'
        shell_color = 'blue'
    
    # Choisir la couleur des lignes selon si on trace une heatmap ou pas
    line_color = 'black' if props == 'htc' else None  # noir si heatmap
    tube_line_color = 'black' if props == 'htc' else tube_color
    cell_line_color = 'black' if props == 'htc' else shell_color
    
    # Lignes de séparation verticale des cellules
    cum_w = np.cumsum(w)
    for x in cum_w * L_total:
        ax.plot([x, x], [0, H_total], color=cell_line_color, lw=1, zorder=12)
    
    "2) Tracer les tubes"
    y_positions = tube_spacing + np.arange(n_passes)*(tube_height + tube_spacing)
    for y_start in y_positions:
        tube_rect = plt.Rectangle((0, y_start), L_total, tube_height,
                                  edgecolor='black', facecolor='white', zorder=10)
        ax.add_patch(tube_rect)
    
    # Lignes bleues pour les passages de cellules
    pass_limits = np.linspace(0, 1, n_passes + 1)
    cell_start = 0.0
    for w_i in w:
        cell_end = cell_start + w_i
        
        # Intersection avec chaque passe
        for p in range(n_passes):
            pass_start, pass_end = pass_limits[p], pass_limits[p+1]
            overlap_start = max(cell_start, pass_start)
            overlap_end = min(cell_end, pass_end)
            
            if overlap_start >= overlap_end:
                continue  # pas d'intersection
            
            # On ne trace que la ligne bleue si c'est **la fin réelle de la cellule**
            draw_end = overlap_end == cell_end
            
            if not draw_end:
                continue
            
            # Position horizontale relative sur le tube
            x_pos = (overlap_end - pass_start) / (pass_end - pass_start) * L_total
            
            # Inversion pour passes paires
            if p % 2 == 1:
                x_pos = L_total - x_pos
            
            y_start_tube = y_positions[p]
            ax.plot([x_pos, x_pos], [y_start_tube, y_start_tube + tube_height],
                    color=tube_line_color, lw=1, zorder=20)
        
        cell_start = cell_end
    
    
    "3) Dynamic tube inlet/outlet pipes"

    pipe_length = 0.25 * L_total
    pipe_height = tube_height * 0.8
    tube_bottom = y_positions[0]
    tube_top = y_positions[-1] + tube_height
    
    # Inlet always at bottom-left
    pipe_in = plt.Rectangle(
        (-pipe_length, tube_bottom + tube_height*0.1),
        pipe_length,
        pipe_height,
        edgecolor='black',
        facecolor='white',
    )
    ax.add_patch(pipe_in)
    
    # Outlet
    if n_passes % 2 == 0:
        # Even number of passes → outlet top-left
        pipe_out = plt.Rectangle(
            (-pipe_length, tube_top - pipe_height - tube_height*0.1),
            pipe_length,
            pipe_height,
            edgecolor='black',
            facecolor='white',
        )
    else:
        # Odd number of passes → outlet top-right
        pipe_out = plt.Rectangle(
            (L_total, tube_top - pipe_height - tube_height*0.1),
            pipe_length,
            pipe_height,
            edgecolor='black',
            facecolor='white',
        )
    ax.add_patch(pipe_out)
    
    # Arrows showing flow
    arrow_length = 0.2 * pipe_length
    arrow_width = 0.02 * H_total
    arrow_head_width = 0.06 * H_total
    arrow_head_length = 0.08 * pipe_length
    
    # Inlet arrow → right
    ax.arrow(
        -pipe_length + 0.05*pipe_length,
        tube_bottom + tube_height*0.1 + pipe_height/2,
        arrow_length,
        0,
        width=arrow_width,
        head_width=arrow_head_width,
        head_length=arrow_head_length,
        fc='black',
        ec='black'
    )
    
    # Outlet arrow → left or right depending on location
    if n_passes % 2 == 0:
        # Even → outlet top-left → arrow left
        ax.arrow(
            -pipe_length + 0.35*pipe_length,
            tube_top - pipe_height - tube_height*0.1 + pipe_height/2,
            -arrow_length,
            0,
            width=arrow_width,
            head_width=arrow_head_width,
            head_length=arrow_head_length,
            fc='black',
            ec='black'
        )
    else:
        # Odd → outlet top-right → arrow right
        ax.arrow(
            L_total + pipe_length - 0.35*pipe_length,
            tube_top - pipe_height - tube_height*0.1 + pipe_height/2,
            arrow_length,
            0,
            width=arrow_width,
            head_width=arrow_head_width,
            head_length=arrow_head_length,
            fc='black',
            ec='black'
        )
    
    "4) Demi-cercles de chaque coté de shell"
    R = H_total / 2
    
    # Demi-cercle gauche
    left_semi = Wedge(
        center=(0, H_total / 2),
        r=R,
        theta1=90,
        theta2=270,
        edgecolor='black',
        facecolor='white',
        zorder=1
    )
    ax.add_patch(left_semi)
    
    # Demi-cercle droit
    right_semi = Wedge(
        center=(L_total, H_total / 2),
        r=R,
        theta1=-90,
        theta2=90,
        edgecolor='black',
        facecolor='white',
        zorder=1
    )
    ax.add_patch(right_semi)
    
    # Ligne horizontale de séparation des demi-cercles
    y_mid = H_total / 2 - 0.002
    
    "5) Horizontal line for collector separation"
    
    def draw_collectors(ax, side, center_x, center_y, R, n_passes, lw=2):
        if n_passes <= 1:
            return
    
        eps = 0.01 * R
    
        # Vertical positions between tube passes
        y_positions = np.linspace(
            center_y - R,
            center_y + R,
            n_passes + 1
        )[1:-1]
    
        for i, y in enumerate(y_positions):
            pass_number = i + 1  # separator is after this pass
    
            # Decide whether to draw separator
            if side == 'left' and pass_number % 2 == 1:
                draw = True
            elif side == 'right' and pass_number % 2 == 0:
                draw = True
            else:
                draw = False
    
            if not draw:
                continue
    
            dy = y - center_y
            dx = np.sqrt(max(R**2 - dy**2, 0.0))
    
            if side == 'left':
                x_start = center_x - dx + eps
                x_end   = center_x - eps
            else:  # right
                x_start = center_x + eps
                x_end   = center_x + dx - eps
    
            ax.plot(
                [x_start, x_end],
                [y, y],
                color='black',
                lw=lw,
                zorder=2
            )
        return


    draw_collectors(
        ax,
        side='left',
        center_x=0,
        center_y=y_mid,
        R=R,
        n_passes=n_passes
    )

    draw_collectors(
        ax,
        side='right',
        center_x=L_total,
        center_y=y_mid,
        R=R,
        n_passes=n_passes
        )
    
    "6) Tuyaux Entrées coté shell"
    # Taille des carrés
    square_size = 0.1 * L_total  # tu peux ajuster la taille
    
    # Carré en haut à droite
    square_top_right = plt.Rectangle(
        (L_total - 1.2*square_size, H_total),  # position coin inférieur gauche
        square_size,
        square_size*0.6,
        edgecolor='black',
        facecolor='lightgrey',
    )
    ax.add_patch(square_top_right)
    
    # Carré en bas à gauche
    square_bottom_left = plt.Rectangle(
        (0.2*square_size, -square_size*0.6),  # position coin inférieur gauche
        square_size,
        square_size*0.6,
        edgecolor='black',
        facecolor='lightgrey',
    )
    ax.add_patch(square_bottom_left)

    # Flèche carré haut-droite (vers le bas)
    ax.arrow(
        L_total - square_size/2 - 0.02,   # x : centre du carré
        H_total + square_size/2 + 0*square_size,  # y : légèrement au-dessus pour que la flèche reste dans le carré
        0,                                # dx = 0 (flèche verticale)
        -0.3*square_size,                 # dy = vers le bas
        width=0.08*square_size,
        head_width=0.2*square_size,
        head_length=0.1*square_size,
        fc='black',
        ec='black',
    )
    
    # Flèche carré bas-gauche (vers le bas)
    ax.arrow(
        square_size/2 + 0.02,             # x : centre du carré
        -0.1*square_size,   # y : légèrement au-dessus du coin supérieur pour que la flèche reste dans le carré
        0,                                # dx = 0
        -0.3*square_size,                 # dy = vers le bas
        width=0.08*square_size,
        head_width=0.2*square_size,
        head_length=0.1*square_size,
        fc='black',
        ec='black',
    )

    # Ajustements figure
    ax.set_xlim(-R - 0.2*L_total, L_total + R + 0.2*L_total)
    # ax.set_ylim(-0.3*H_total, 1.3*H_total)
    ax.axis('off')
    
    "7) Propriété thermo"
    
    if props == 'htc':
        # Valeurs des tubes
        prop_values_tubes = np.array(HX.alpha_c)
        prop_values_shell = np.array(HX.alpha_h)
    
        if len(prop_values_tubes) != len(w) or len(prop_values_shell) != len(w):
            raise ValueError("HX.alpha_c et HX.alpha_h doivent avoir la même longueur que HX.w")
    
        # --- Normalisation commune et colormap ---
        all_values = np.concatenate([prop_values_tubes, prop_values_shell])
        norm = mcolors.Normalize(
            vmin=legend_range[0] if legend_range is not None else np.min([HX.alpha_c, HX.alpha_h]),
            vmax=legend_range[1] if legend_range is not None else np.max([HX.alpha_c, HX.alpha_h]),
            clip=True
        )
        viridis = plt.cm.viridis
        cmap = LinearSegmentedColormap.from_list(
            'viridis_light', viridis(np.linspace(0.3, 1.0, 256))
        )    
        # --- Heatmap des tubes ---
        cell_start = 0.0
                
        for idx, w_i in enumerate(w):
            cell_end = cell_start + w_i
            
            for p in range(n_passes):
                pass_start, pass_end = pass_limits[p], pass_limits[p+1]
                overlap_start = max(cell_start, pass_start)
                overlap_end = min(cell_end, pass_end)
                if overlap_start >= overlap_end:
                    continue
                
                frac_start = (overlap_start - pass_start) / (pass_end - pass_start)
                frac_end   = (overlap_end   - pass_start) / (pass_end - pass_start)
                
                if p % 2 == 1:
                    frac_start, frac_end = 1 - frac_end, 1 - frac_start
                
                x0 = frac_start * L_total
                width = (frac_end - frac_start) * L_total
                y0 = y_positions[p]
                
                rect = plt.Rectangle(
                    (x0, y0),
                    width,
                    tube_height-0.002,
                    facecolor=cmap(norm(prop_values_tubes[idx])),
                    edgecolor='none',
                    zorder=15
                )
                ax.add_patch(rect)
            
            cell_start = cell_end
    
        # --- Heatmap de la shell ---
        cell_start = 0.0
        for idx, w_i in enumerate(w):
            cell_end = cell_start + w_i
            x0 = cell_start * L_total
            width = w_i * L_total
    
            rect_shell = plt.Rectangle(
                (x0, 0),
                width,
                H_total-0.005,
                facecolor=cmap(norm(prop_values_shell[idx])),
                edgecolor='none',
                zorder=5
            )
            ax.add_patch(rect_shell)
            cell_start = cell_end
    
        # --- Ajouter la légende ---
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # nécessaire pour colorbar
        
        cbar = plt.colorbar(
            sm,
            ax=ax,
            fraction=0.03,
            pad=0.02,
        )
        
        cbar.set_label(
            r"Heat Transfer Coefficient [$\mathrm{W\,m^{-2}\,K^{-1}}$]",
            fontsize=13
        )
        
        if legend_range is not None:
            cbar.set_ticks(np.linspace(legend_range[0], legend_range[1], 5))
     
    elif props is None:
        """
        Default coloring when no property is shown:
        - Shell side: red gradient (hot)
        - Tube side: blue gradient (cold)
        """

        # Normalized axial coordinate based on cells
        cell_centers = np.cumsum(w) - 0.5 * w
        norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

        if HX.params['Shell_Side'] == 'H':
            cmap_shell = plt.cm.Reds
            cmap_tube = plt.cm.Blues_r
        else:
            cmap_shell = plt.cm.Blues
            cmap_tube = plt.cm.Reds 
            
        def lightened_norm(x, vmin=0.2, vmax=0.6):
            """Compress normalized values into a lighter colormap range."""
            return vmin + (vmax - vmin) * x

        # ---- Tube side (cold – blue gradient) ----
        cell_start = 0.0
        for idx, w_i in enumerate(w):
            cell_end = cell_start + w_i

            for p in range(n_passes):
                pass_start, pass_end = pass_limits[p], pass_limits[p+1]
                overlap_start = max(cell_start, pass_start)
                overlap_end   = min(cell_end, pass_end)

                if overlap_start >= overlap_end:
                    continue

                frac_start = (overlap_start - pass_start) / (pass_end - pass_start)
                frac_end   = (overlap_end   - pass_start) / (pass_end - pass_start)

                if p % 2 == 1:
                    frac_start, frac_end = 1 - frac_end, 1 - frac_start

                x0 = frac_start * L_total
                width = (frac_end - frac_start) * L_total
                y0 = y_positions[p]

                tube_margin = 0.08 * tube_height  # NEW

                rect = plt.Rectangle(
                    (x0 + 0.002, y0),
                    width - 0.002,
                    tube_height-0.002,
                    facecolor=cmap_tube(lightened_norm(norm(cell_centers[idx]))),
                    edgecolor='none',
                    zorder=15
                )

                ax.add_patch(rect)

            cell_start = cell_end

        # ---- Shell side (hot – red gradient) ----
        cell_start = 0.0
        for idx, w_i in enumerate(w):
            x0 = cell_start * L_total
            width = w_i * L_total

            rect_shell = plt.Rectangle(
                (x0+0.002, 0),
                width-0.005,
                H_total-0.002,
                facecolor=cmap_shell(lightened_norm(norm(cell_centers[idx]))),
                edgecolor='none',
                zorder=5
            )
            ax.add_patch(rect_shell)

            cell_start += w_i

        
    "8) Save and return"
        
    if save_path is not None:
        # Save as SVG
        plt.savefig(save_path, format="svg", bbox_inches="tight")
        
    plt.show()
    return plt

