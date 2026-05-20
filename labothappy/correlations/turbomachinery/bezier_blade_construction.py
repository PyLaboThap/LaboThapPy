# """
# blade_bezier.py
# ---------------
# Reconstruction d'un profil de pale 2D (cartésien x, y) via des courbes de
# Bézier cubiques indépendantes pour le pressure side et le suction side.

# Paramètres géométriques d'entrée
# ─────────────────────────────────
#   x_le, y_le          : coordonnées du leading edge
#   x_te, y_te          : coordonnées du trailing edge
#   beta1_deg           : angle de pale à l'entrée  [degrés, mesuré depuis l'axe x]
#   beta2_deg           : angle de pale à la sortie [degrés, mesuré depuis l'axe x]
#   d1, d2              : longueurs des bras de contrôle aux extrémités (camber)
#   t_max_ss            : épaisseur max du suction side  (normalisée par la corde)
#   t_max_ps            : épaisseur max du pressure side (normalisée par la corde)
#   s_tmax_ss           : position relative (0–1) de t_max sur le SS
#   s_tmax_ps           : position relative (0–1) de t_max sur le PS
#   t_te_ss, t_te_ps    : demi-épaisseur de bord de fuite SS / PS (normalisée)
#   d1_t, d2_t          : bras de contrôle de la Bézier d'épaisseur (entrée/sortie)
#   n_points            : nombre de points de discrétisation (défaut 200)

# Sorties  (dict retourné par build_blade)
# ────────────────────────────────────────
#   "camber"      – (n, 2) ndarray  points de la camber line
#   "ss"          – (n, 2) ndarray  points du suction side
#   "ps"          – (n, 2) ndarray  points du pressure side
#   "ctrl_camber" – (4, 2) ndarray  points de contrôle de la camber line
#   "ctrl_ss"     – (4, 2) ndarray  points de contrôle du suction side
#   "ctrl_ps"     – (4, 2) ndarray  points de contrôle du pressure side
#   "chord"       – float           longueur de corde LE→TE
# """

# from __future__ import annotations

# import numpy as np
# import matplotlib.pyplot as plt


# # ─────────────────────────────────────────────
# #  Évaluation d'une Bézier cubique
# # ─────────────────────────────────────────────

# def _bezier_cubic(P: np.ndarray, t: np.ndarray) -> np.ndarray:
#     """
#     Évalue une Bézier cubique définie par 4 points de contrôle P (4×2)
#     pour un vecteur de paramètres t ∈ [0, 1].

#     Retourne un tableau (n, 2).
#     """
#     t = t[:, None]          # (n, 1) pour broadcast
#     B = (  (1 - t)**3       * P[0]
#          + 3*(1-t)**2 * t   * P[1]
#          + 3*(1-t)   * t**2 * P[2]
#          +             t**3 * P[3] )
#     return B                # (n, 2)


# def _bezier_derivative(P: np.ndarray, t: np.ndarray) -> np.ndarray:
#     """Dérivée première dB/dt d'une Bézier cubique."""
#     t = t[:, None]
#     dB = (  3*(1-t)**2       * (P[1] - P[0])
#            + 6*(1-t)*t       * (P[2] - P[1])
#            + 3*t**2           * (P[3] - P[2]) )
#     return dB               # (n, 2)


# # ─────────────────────────────────────────────
# #  Construction de la camber line
# # ─────────────────────────────────────────────

# def _build_camber(
#     x_le: float, y_le: float,
#     x_te: float, y_te: float,
#     beta1_deg: float, beta2_deg: float,
#     d1: float, d2: float,
#     n: int,
# ) -> tuple[np.ndarray, np.ndarray]:
#     """
#     Construit la camber line par une Bézier cubique.

#     Les tangentes en P0 et P3 sont imposées par les angles de pale β1, β2.
#     d1 et d2 contrôlent la longueur des bras de contrôle (degrés de liberté).

#     Retourne (points (n,2), ctrl_pts (4,2)).
#     """
#     beta1 = np.radians(beta1_deg)
#     beta2 = np.radians(beta2_deg)

#     P0 = np.array([x_le, y_le])
#     P3 = np.array([x_te, y_te])

#     # Directions tangentes unitaires
#     tan1 = np.array([np.cos(beta1), np.sin(beta1)])
#     tan2 = np.array([np.cos(beta2), np.sin(beta2)])

#     # Points de contrôle intérieurs
#     P1 = P0 + d1 * tan1
#     P2 = P3 - d2 * tan2          # on remonte depuis TE

#     ctrl = np.array([P0, P1, P2, P3])
#     t = np.linspace(0.0, 1.0, n)
#     pts = _bezier_cubic(ctrl, t)
#     return pts, ctrl


# # ─────────────────────────────────────────────
# #  Distribution d'épaisseur via Bézier cubique
# # ─────────────────────────────────────────────

# def _thickness_bezier(
#     t_le: float,
#     t_max: float,
#     s_tmax: float,
#     t_te: float,
#     d1_t: float,
#     d2_t: float,
#     n: int,
# ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
#     """
#     Construit une distribution d'épaisseur scalaire e(s), s ∈ [0,1],
#     par une Bézier cubique dans l'espace (s, e).

#     Points de contrôle :
#       Q0 = (0,       t_le)           bord d'attaque
#       Q1 = (d1_t,    t_max)          contrôle la pente à l'entrée
#       Q2 = (s_tmax,  t_max)          pic d'épaisseur (tangente horizontale)
#       Q3 = (1,       t_te)           bord de fuite

#     Retourne (s_vals (n,), e_vals (n,), ctrl_pts (4,2)).

#     Note : la sortie est en coordonnées (s, e) — à mapper ensuite
#            sur la normale à la camber line.
#     """
#     Q0 = np.array([0.0,     t_le  ])
#     Q1 = np.array([d1_t,    t_max ])
#     Q2 = np.array([s_tmax,  t_max ])
#     Q3 = np.array([1.0,     t_te  ])

#     ctrl = np.array([Q0, Q1, Q2, Q3])
#     t_param = np.linspace(0.0, 1.0, n)
#     pts = _bezier_cubic(ctrl, t_param)   # (n, 2) : (s, e)

#     # On ré-échantillonne uniformément en s ∈ [0,1] par interpolation
#     s_raw = pts[:, 0]
#     e_raw = pts[:, 1]
#     s_uniform = np.linspace(0.0, 1.0, n)
#     e_uniform = np.interp(s_uniform, s_raw, e_raw)

#     return s_uniform, e_uniform, ctrl


# # ─────────────────────────────────────────────
# #  Projection de l'épaisseur sur la normale
# # ─────────────────────────────────────────────

# def _offset_curve(
#     camber: np.ndarray,
#     ctrl_camber: np.ndarray,
#     thickness: np.ndarray,
#     sign: float,
#     chord: float,
# ) -> tuple[np.ndarray, np.ndarray]:
#     """
#     Déplace la camber line de ±thickness (normalisée par la corde)
#     dans la direction normale locale.

#     sign = +1  → suction side (toward low-pressure)
#     sign = -1  → pressure side

#     Retourne (offset_pts (n,2), ctrl_approx (4,2)).
#     """
#     n = camber.shape[0]
#     t_param = np.linspace(0.0, 1.0, n)

#     # Dérivée tangentielle
#     dB = _bezier_derivative(ctrl_camber, t_param)       # (n, 2)
#     norms = np.linalg.norm(dB, axis=1, keepdims=True)
#     norms = np.where(norms < 1e-12, 1e-12, norms)       # évite /0
#     tangents = dB / norms                               # (n, 2) unitaires

#     # Normale = rotation 90° de la tangente
#     # [cos, sin] → [-sin, cos] (vers suction) ou [sin, -cos] (vers pressure)
#     normals = sign * np.column_stack([-tangents[:, 1], tangents[:, 0]])

#     # Décalage
#     offset = camber + normals * (thickness * chord)[:, None]

#     # Points de contrôle approximatifs (LE et TE exacts, intermédiaires estimés)
#     ctrl_approx = np.array([
#         offset[0],
#         offset[n // 3],
#         offset[2 * n // 3],
#         offset[-1],
#     ])

#     return offset, ctrl_approx


# # ─────────────────────────────────────────────
# #  Longueurs caractéristiques
# # ─────────────────────────────────────────────

# def meridional_length(geo: dict) -> float:
#     """
#     Longueur méridionale de la pale = longueur d'arc de la camber line.

#     Calculée par sommation des segments discrets sur geo["camber"].
#     Converge à < 1e-5 de précision relative pour n_points >= 200.

#     Returns
#     -------
#     float : longueur d'arc  [même unité que x, y]
#     """
#     pts  = geo["camber"]                             # (n, 2)
#     segs = np.diff(pts, axis=0)                      # (n-1, 2)
#     return float(np.sum(np.linalg.norm(segs, axis=1)))


# def axial_length(geo: dict) -> float:
#     """
#     Longueur axiale de la pale = projection de la camber line sur l'axe x.

#     Pour une géométrie cartésienne 2D, c'est simplement x_te - x_le.
#     Les points de contrôle intermédiaires n'affectent pas cette valeur
#     car la courbe de Bézier est garantie de passer par P0 et P3.

#     Returns
#     -------
#     float : longueur axiale  [même unité que x]
#     """
#     return float(geo["camber"][-1, 0] - geo["camber"][0, 0])


# # ─────────────────────────────────────────────
# #  Fonction principale
# # ─────────────────────────────────────────────

# def build_blade(
#     # Géométrie globale
#     x_le: float = 0.0,
#     y_le: float = 0.0,
#     x_te: float = 1.0,
#     y_te: float = 0.0,
#     # Angles de pale [degrés]
#     beta1_deg: float = 50.0,
#     beta2_deg: float = 20.0,
#     # Bras de contrôle camber line
#     d1: float = 0.35,
#     d2: float = 0.35,
#     # Épaisseur suction side
#     t_max_ss: float = 0.08,
#     s_tmax_ss: float = 0.35,
#     t_te_ss: float = 0.005,
#     d1_t_ss: float = 0.15,
#     d2_t_ss: float = 0.65,
#     # Épaisseur pressure side
#     t_max_ps: float = 0.06,
#     s_tmax_ps: float = 0.40,
#     t_te_ps: float = 0.005,
#     d1_t_ps: float = 0.15,
#     d2_t_ps: float = 0.65,
#     # Discrétisation
#     n_points: int = 300,
# ) -> dict:
#     """
#     Reconstruit un profil de pale 2D cartésien via des Bézier cubiques.

#     Paramètres
#     ----------
#     x_le, y_le      : coordonnées du leading edge
#     x_te, y_te      : coordonnées du trailing edge
#     beta1_deg       : angle de pale à l'entrée [°] depuis l'axe x
#     beta2_deg       : angle de pale à la sortie [°] depuis l'axe x
#     d1, d2          : longueurs des bras de contrôle de la camber line
#                       (normalisées par la longueur LE→TE, typiquement 0.2–0.5)
#     t_max_ss/ps     : épaisseur maximale SS/PS normalisée par la corde
#     s_tmax_ss/ps    : position curviligne (0–1) du max d'épaisseur SS/PS
#     t_te_ss/ps      : demi-épaisseur au bord de fuite SS/PS (normalisée)
#     d1_t_ss/ps      : bras de contrôle Bézier épaisseur côté LE (0–1)
#     d2_t_ss/ps      : bras de contrôle Bézier épaisseur côté TE (0–1)
#     n_points        : nombre de points de discrétisation

#     Retourne
#     --------
#     dict avec les clés :
#         "camber"      – (n, 2) ndarray
#         "ss"          – (n, 2) ndarray  suction side
#         "ps"          – (n, 2) ndarray  pressure side
#         "ctrl_camber" – (4, 2) ndarray
#         "ctrl_ss"     – (4, 2) ndarray
#         "ctrl_ps"     – (4, 2) ndarray
#         "chord"       – float
#     """
#     # 1. Camber line
#     camber, ctrl_camber = _build_camber(
#         x_le, y_le, x_te, y_te,
#         beta1_deg, beta2_deg,
#         d1, d2, n_points,
#     )

#     # Corde réelle (distance LE–TE)
#     chord = float(np.linalg.norm(np.array([x_te - x_le, y_te - y_le])))

#     # 2. Distributions d'épaisseur
#     _, e_ss, ctrl_ss = _thickness_bezier(
#         t_le=0.0, t_max=t_max_ss, s_tmax=s_tmax_ss,
#         t_te=t_te_ss, d1_t=d1_t_ss, d2_t=d2_t_ss, n=n_points,
#     )
#     _, e_ps, ctrl_ps = _thickness_bezier(
#         t_le=0.0, t_max=t_max_ps, s_tmax=s_tmax_ps,
#         t_te=t_te_ps, d1_t=d1_t_ps, d2_t=d2_t_ps, n=n_points,
#     )

#     # 3. Projection normale sur la camber line
#     ss, ctrl_ss = _offset_curve(camber, ctrl_camber, e_ss, +1.0, chord)
#     ps, ctrl_ps = _offset_curve(camber, ctrl_camber, e_ps, -1.0, chord)

#     return {
#         "camber":      camber,
#         "ss":          ss,
#         "ps":          ps,
#         "ctrl_camber": ctrl_camber,
#         "ctrl_ss":     ctrl_ss,
#         "ctrl_ps":     ctrl_ps,
#         "chord":       chord,
#     }


# # ─────────────────────────────────────────────
# #  Visualisation
# # ─────────────────────────────────────────────

# def plot_blade(
#     geo: dict,
#     show_camber: bool = True,
#     show_ctrl: bool = True,
#     title: str = "Profil de pale — Bézier cubique",
#     ax: plt.Axes | None = None,
# ) -> plt.Figure:
#     """
#     Trace le profil complet (SS, PS, camber line, points de contrôle).

#     Paramètres
#     ----------
#     geo         : dict retourné par build_blade()
#     show_camber : afficher la camber line
#     show_ctrl   : afficher les polygones de contrôle Bézier
#     title       : titre du graphique
#     ax          : axes matplotlib existants (optionnel)

#     Retourne
#     --------
#     Figure matplotlib.
#     """
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(11, 5))
#     else:
#         fig = ax.get_figure()

#     camber      = geo["camber"]
#     ss          = geo["ss"]
#     ps          = geo["ps"]
#     ctrl_camber = geo["ctrl_camber"]
#     ctrl_ss     = geo["ctrl_ss"]
#     ctrl_ps     = geo["ctrl_ps"]

#     # ── Remplissage du profil ──────────────────────────────────────────
#     contour_x = np.concatenate([ss[:, 0], ps[::-1, 0]])
#     contour_y = np.concatenate([ss[:, 1], ps[::-1, 1]])
#     ax.fill(contour_x, contour_y, color="#d0e8f5", alpha=0.6, zorder=1)

#     # ── Suction side ──────────────────────────────────────────────────
#     ax.plot(ss[:, 0], ss[:, 1],
#             color="#1a6fa8", lw=2.0, label="Suction side", zorder=3)

#     # ── Pressure side ─────────────────────────────────────────────────
#     ax.plot(ps[:, 0], ps[:, 1],
#             color="#c0392b", lw=2.0, label="Pressure side", zorder=3)

#     # ── Camber line ───────────────────────────────────────────────────
#     if show_camber:
#         ax.plot(camber[:, 0], camber[:, 1],
#                 color="#2ecc71", lw=1.4, ls="--", label="Camber line", zorder=4)

#     # ── Polygones de contrôle ─────────────────────────────────────────
#     if show_ctrl:
#         for ctrl, color, label in [
#             (ctrl_camber, "#2ecc71", "Ctrl camber"),
#             (ctrl_ss,     "#1a6fa8", "Ctrl SS"),
#             (ctrl_ps,     "#c0392b", "Ctrl PS"),
#         ]:
#             ax.plot(ctrl[:, 0], ctrl[:, 1],
#                     "o--", color=color, lw=0.8, ms=5, alpha=0.6,
#                     label=label, zorder=5)

#     # ── Marqueurs LE / TE ─────────────────────────────────────────────
#     le = camber[0]
#     te = camber[-1]
#     ax.plot(*le, "ko", ms=6, zorder=6)
#     ax.plot(*te, "ks", ms=6, zorder=6)
#     ax.annotate("LE", le, textcoords="offset points", xytext=(-14, 6), fontsize=9)
#     ax.annotate("TE", te, textcoords="offset points", xytext=(5, 6),   fontsize=9)

#     # ── Axes & légende ────────────────────────────────────────────────
#     ax.set_aspect("equal")
#     ax.set_xlabel("x  [m]")
#     ax.set_ylabel("y  [m]")
#     ax.set_title(title, fontsize=12, fontweight="bold")
#     ax.legend(loc="upper right", fontsize=8, framealpha=0.8)
#     ax.grid(True, ls=":", alpha=0.4)
#     fig.tight_layout()
#     return fig


# # ─────────────────────────────────────────────
# #  Visualisation du passage
# # ─────────────────────────────────────────────

# def plot_passage(
#     geo: dict,
#     pitch: float,
#     show_camber: bool = True,
#     show_ctrl: bool = False,
#     show_lengths: bool = True,
#     title: str = "Passage de pale — Bézier cubique",
#     ax: plt.Axes | None = None,
# ) -> plt.Figure:
#     """
#     Trace le passage complet entre deux pales consécutives.

#     La pale de référence est celle de geo. La pale voisine est obtenue
#     par translation de +pitch selon y (direction tangentielle).
#     Le passage fluide entre le SS de la pale de référence et le PS
#     de la pale voisine est colorié.

#     Paramètres
#     ----------
#     geo          : dict retourné par build_blade()
#     pitch        : espacement entre pales [même unité que x, y]
#     show_camber  : afficher les camber lines
#     show_ctrl    : afficher les polygones de contrôle Bézier
#     show_lengths : annoter les longueurs méridionale et axiale
#     title        : titre du graphique
#     ax           : axes matplotlib existants (optionnel)

#     Retourne
#     --------
#     Figure matplotlib.
#     """
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(12, 7))
#     else:
#         fig = ax.get_figure()

#     camber = geo["camber"]
#     ss     = geo["ss"]
#     ps     = geo["ps"]

#     # Pale voisine : translation de +pitch en y
#     ss_next     = ss     + np.array([0.0, pitch])
#     ps_next     = ps     + np.array([0.0, pitch])
#     camber_next = camber + np.array([0.0, pitch])

#     # ── Passage fluide ────────────────────────────────────────────────
#     # Délimité par : SS de la pale ref (bas) + PS de la pale voisine (haut)
#     # Contour : SS LE→TE, puis PS_next TE→LE
#     pass_x = np.concatenate([ss[:, 0],      ps_next[::-1, 0]])
#     pass_y = np.concatenate([ss[:, 1],      ps_next[::-1, 1]])
#     ax.fill(pass_x, pass_y, color="#e8f4e8", alpha=0.7, zorder=1,
#             label="Passage fluide")

#     # ── Pale de référence ─────────────────────────────────────────────
#     blade_x = np.concatenate([ss[:, 0], ps[::-1, 0]])
#     blade_y = np.concatenate([ss[:, 1], ps[::-1, 1]])
#     ax.fill(blade_x, blade_y, color="#d0e8f5", alpha=0.8, zorder=2)
#     ax.plot(ss[:, 0], ss[:, 1], color="#1a6fa8", lw=2.0,
#             label="Suction side", zorder=4)
#     ax.plot(ps[:, 0], ps[:, 1], color="#c0392b", lw=2.0,
#             label="Pressure side", zorder=4)

#     # ── Pale voisine ──────────────────────────────────────────────────
#     blade_next_x = np.concatenate([ss_next[:, 0], ps_next[::-1, 0]])
#     blade_next_y = np.concatenate([ss_next[:, 1], ps_next[::-1, 1]])
#     ax.fill(blade_next_x, blade_next_y, color="#d0e8f5", alpha=0.8, zorder=2)
#     ax.plot(ss_next[:, 0], ss_next[:, 1], color="#1a6fa8", lw=2.0,
#             zorder=4)
#     ax.plot(ps_next[:, 0], ps_next[:, 1], color="#c0392b", lw=2.0,
#             zorder=4)

#     # ── Camber lines ──────────────────────────────────────────────────
#     if show_camber:
#         ax.plot(camber[:, 0],      camber[:, 1],      color="#2ecc71",
#                 lw=1.4, ls="--", label="Camber line", zorder=5)
#         ax.plot(camber_next[:, 0], camber_next[:, 1], color="#2ecc71",
#                 lw=1.4, ls="--", zorder=5)

#     # ── Polygones de contrôle ─────────────────────────────────────────
#     if show_ctrl:
#         for ctrl, color in [
#             (geo["ctrl_camber"],              "#2ecc71"),
#             (geo["ctrl_ss"],                  "#1a6fa8"),
#             (geo["ctrl_ps"],                  "#c0392b"),
#             (geo["ctrl_camber"] + [0, pitch], "#2ecc71"),
#             (geo["ctrl_ss"]     + [0, pitch], "#1a6fa8"),
#             (geo["ctrl_ps"]     + [0, pitch], "#c0392b"),
#         ]:
#             ax.plot(ctrl[:, 0], ctrl[:, 1], "o--", color=color,
#                     lw=0.8, ms=4, alpha=0.5, zorder=6)

#     # ── Marqueurs LE / TE ─────────────────────────────────────────────
#     for cam, offset_label in [(camber, ""), (camber_next, "")]:
#         le, te = cam[0], cam[-1]
#         ax.plot(*le, "ko", ms=5, zorder=7)
#         ax.plot(*te, "ks", ms=5, zorder=7)
#     ax.annotate("LE", camber[0], textcoords="offset points",
#                 xytext=(-18, 4), fontsize=8)
#     ax.annotate("TE", camber[-1], textcoords="offset points",
#                 xytext=(5, 4), fontsize=8)

#     # ── Annotations de longueur ───────────────────────────────────────
#     if show_lengths:
#         L_mer = meridional_length(geo)
#         L_ax  = axial_length(geo)

#         # Flèche longueur axiale (sous la pale)
#         y_arrow = ps[:, 1].min() - 0.04 * pitch
#         x0, x1  = camber[0, 0], camber[-1, 0]
#         ax.annotate("", xy=(x1, y_arrow), xytext=(x0, y_arrow),
#                     arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2))
#         ax.text((x0 + x1) / 2, y_arrow - 0.03 * pitch,
#                 f"$L_{{ax}}$ = {L_ax:.3f} m",
#                 ha="center", va="top", fontsize=8, color="#555")

#         # Flèche longueur méridionale (le long de la camber)
#         # Positionnée à mi-passage en y, décalée à gauche
#         x_mid  = camber[:, 0].mean()
#         y_mid  = camber[:, 1].mean()
#         ax.annotate(f"$L_{{mer}}$ = {L_mer:.3f} m",
#                     xy=(x_mid, y_mid),
#                     xytext=(x_mid - 0.12 * L_ax, y_mid + 0.08 * pitch),
#                     fontsize=8, color="#2ecc71",
#                     arrowprops=dict(arrowstyle="->", color="#2ecc71", lw=1.0))

#         # Flèche pitch (à droite)
#         x_pitch = camber[:, 0].max() + 0.05 * L_ax
#         y_bot   = ps[0, 1]
#         y_top   = ps_next[0, 1]
#         ax.annotate("", xy=(x_pitch, y_top), xytext=(x_pitch, y_bot),
#                     arrowprops=dict(arrowstyle="<->", color="#888", lw=1.2))
#         ax.text(x_pitch + 0.02 * L_ax, (y_bot + y_top) / 2,
#                 f"pitch = {pitch:.3f} m",
#                 ha="left", va="center", fontsize=8, color="#888")

#     # ── Axes & légende ────────────────────────────────────────────────
#     ax.set_aspect("equal")
#     ax.set_xlabel("x  [m]")
#     ax.set_ylabel("y  [m]")
#     ax.set_title(title, fontsize=12, fontweight="bold")
#     ax.legend(loc="upper right", fontsize=8, framealpha=0.8)
#     ax.grid(True, ls=":", alpha=0.35)
#     fig.tight_layout()
#     return fig


# # ─────────────────────────────────────────────
# #  Canal méridional (hub & tip)
# # ─────────────────────────────────────────────

# def build_meridional_channel(
#     # Inlet — plan axial : hub et tip partagent x_in, diffèrent en r
#     x_in:        float = 0.0,
#     r_hub_in:    float = 0.03,
#     r_tip_in:    float = 0.10,
#     # Outlet — plan radial : hub et tip partagent r_out, diffèrent en x
#     r_out:       float = 0.25,
#     x_hub_out:   float = 0.08,
#     x_tip_out:   float = 0.04,
#     # Angles de tangente hub [degrés depuis l'axe x] : entrée, sortie
#     alpha_hub_in:  float = 10.0,
#     alpha_hub_out: float = 82.0,
#     # Angles de tangente tip [degrés depuis l'axe x] : entrée, sortie
#     alpha_tip_in:  float =  5.0,
#     alpha_tip_out: float = 88.0,
#     # Longueurs des bras de contrôle (normalisées par la distance axiale de chaque contour)
#     d_hub_in:  float = 0.50,
#     d_hub_out: float = 0.50,
#     d_tip_in:  float = 0.40,
#     d_tip_out: float = 0.40,
#     # Discrétisation
#     n_points: int = 300,
# ) -> dict:
#     """
#     Construit les contours hub et tip du canal méridional d'un impeller
#     centrifuge par deux courbes de Bézier cubiques dans le plan (x, r).

#     Convention géométrique
#     ──────────────────────
#     Inlet  — plan axial    : x = x_in   (constant), r ∈ [r_hub_in, r_tip_in]
#     Outlet — plan radial   : r = r_out  (constant), x ∈ [x_tip_out, x_hub_out]

#     Les deux contours ont donc :
#       Hub  : (x_in, r_hub_in) → (x_hub_out, r_out)
#       Tip  : (x_in, r_tip_in) → (x_tip_out, r_out)

#     Paramètres
#     ----------
#     x_in                : abscisse axiale de l'inlet
#     r_hub_in, r_tip_in  : rayons hub / tip à l'inlet
#     r_out               : rayon de l'outlet (commun hub et tip)
#     x_hub_out           : abscisse axiale du hub à l'outlet  (> x_tip_out)
#     x_tip_out           : abscisse axiale du tip à l'outlet
#     alpha_hub_in/out    : angles de tangente du hub [°] depuis l'axe x
#     alpha_tip_in/out    : angles de tangente du tip [°] depuis l'axe x
#     d_hub_in/out        : longueurs bras de contrôle hub (× Δx du contour)
#     d_tip_in/out        : longueurs bras de contrôle tip (× Δx du contour)
#     n_points            : nombre de points de discrétisation

#     Retourne
#     --------
#     dict :
#         "hub"        – (n, 2) ndarray  contour hub  dans (x, r)
#         "tip"        – (n, 2) ndarray  contour tip  dans (x, r)
#         "ctrl_hub"   – (4, 2) ndarray  points de contrôle hub
#         "ctrl_tip"   – (4, 2) ndarray  points de contrôle tip
#         "x_in"       – float
#         "r_hub_in"   – float
#         "r_tip_in"   – float
#         "r_out"      – float
#         "x_hub_out"  – float
#         "x_tip_out"  – float
#     """
#     def _contour(x0, r0, x1, r1, a_in, a_out, d_in, d_out):
#         P0 = np.array([x0, r0])
#         P3 = np.array([x1, r1])
#         dx = abs(x1 - x0) if abs(x1 - x0) > 1e-12 else abs(r1 - r0)
#         a0 = np.radians(a_in)
#         a1 = np.radians(a_out)
#         tan_in  = np.array([np.cos(a0), np.sin(a0)])
#         tan_out = np.array([np.cos(a1), np.sin(a1)])
#         P1 = P0 + d_in  * dx * tan_in
#         P2 = P3 - d_out * dx * tan_out
#         ctrl = np.array([P0, P1, P2, P3])
#         t    = np.linspace(0.0, 1.0, n_points)
#         pts  = _bezier_cubic(ctrl, t)
#         return pts, ctrl

#     hub, ctrl_hub = _contour(
#         x_in, r_hub_in, x_hub_out, r_out,
#         alpha_hub_in, alpha_hub_out, d_hub_in, d_hub_out,
#     )
#     tip, ctrl_tip = _contour(
#         x_in, r_tip_in, x_tip_out, r_out,
#         alpha_tip_in, alpha_tip_out, d_tip_in, d_tip_out,
#     )

#     return {
#         "hub":       hub,
#         "tip":       tip,
#         "ctrl_hub":  ctrl_hub,
#         "ctrl_tip":  ctrl_tip,
#         "x_in":      x_in,
#         "r_hub_in":  r_hub_in,
#         "r_tip_in":  r_tip_in,
#         "r_out":     r_out,
#         "x_hub_out": x_hub_out,
#         "x_tip_out": x_tip_out,
#     }


# def plot_meridional(
#     channel: dict,
#     show_ctrl:  bool = True,
#     show_annot: bool = True,
#     title: str = "Vue méridionale — Hub & Tip",
#     ax: plt.Axes | None = None,
# ) -> plt.Figure:
#     """
#     Vue méridionale (r horizontal, x axial inversé) d'un impeller centrifuge.

#     Inlet  : segment vertical   r ∈ [r_hub_in, r_tip_in] à x = x_in
#     Outlet : segment horizontal x ∈ [x_tip_out, x_hub_out] à r = r_out
#     """
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(8, 7))
#     else:
#         fig = ax.get_figure()

#     hub      = channel["hub"]        # (n,2) en (x, r)
#     tip      = channel["tip"]
#     ctrl_hub = channel["ctrl_hub"]
#     ctrl_tip = channel["ctrl_tip"]

#     x_in      = channel["x_in"]
#     r_hub_in  = channel["r_hub_in"]
#     r_tip_in  = channel["r_tip_in"]
#     r_out     = channel["r_out"]
#     x_hub_out = channel["x_hub_out"]
#     x_tip_out = channel["x_tip_out"]

#     def _r(pts): return pts[:, 1]
#     def _x(pts): return pts[:, 0]

#     # ── Passage fluide ────────────────────────────────────────────────
#     # Contour fermé dans le plan (r, x) :
#     #   hub LE→TE  + outlet (x_hub_out→x_tip_out à r=r_out) + tip TE→LE
#     outlet_r = np.array([r_out, r_out])
#     outlet_x = np.array([x_hub_out, x_tip_out])

#     pass_r = np.concatenate([_r(hub), outlet_r, _r(tip)[::-1]])
#     pass_x = np.concatenate([_x(hub), outlet_x, _x(tip)[::-1]])
#     ax.fill(pass_r, pass_x, color="#dceefb", alpha=0.7, zorder=1,
#             label="Passage fluide")
#     ax.fill(pass_r, pass_x, hatch="///", facecolor="none",
#             edgecolor="#aac8e0", lw=0.4, alpha=0.4, zorder=2)

#     # ── Contours hub / tip ────────────────────────────────────────────
#     ax.plot(_r(hub), _x(hub), color="#8B4513", lw=2.5, label="Hub",          zorder=4)
#     ax.plot(_r(tip), _x(tip), color="#2c7bb6", lw=2.5, label="Tip (shroud)", zorder=4)

#     # ── Inlet : vertical (r varie, x constant) ────────────────────────
#     ax.plot([r_hub_in, r_tip_in], [x_in, x_in],
#             color="#e67e22", lw=2.0, ls="--", label="Inlet", zorder=3)

#     # ── Outlet : horizontal (x varie, r constant) ─────────────────────
#     ax.plot([r_out, r_out], [x_hub_out, x_tip_out],
#             color="#27ae60", lw=2.0, ls="--", label="Outlet", zorder=3)

#     # ── Axe de rotation ───────────────────────────────────────────────
#     ax.axvline(0, color="#ccc", lw=0.8, ls=":", zorder=0)
#     ax.text(0.002, x_in, "axe", ha="left", va="center", fontsize=7, color="#aaa")

#     # ── Polygones de contrôle ─────────────────────────────────────────
#     if show_ctrl:
#         for ctrl, color in [(ctrl_hub, "#8B4513"), (ctrl_tip, "#2c7bb6")]:
#             ax.plot(ctrl[:, 1], ctrl[:, 0], "o--", color=color,
#                     lw=0.9, ms=5, alpha=0.55, zorder=5)

#     # ── Annotations ───────────────────────────────────────────────────
#     if show_annot:
#         dr     = r_out
#         margin = 0.03 * dr

#         # r_hub_in
#         ax.annotate("", xy=(r_hub_in, x_in), xytext=(0, x_in),
#                     arrowprops=dict(arrowstyle="<->", color="#8B4513", lw=1.0))
#         ax.text(r_hub_in / 2, x_in - margin,
#                 f"$r_{{hub,in}}$={r_hub_in:.3f} m",
#                 ha="center", va="bottom", fontsize=7.5, color="#8B4513")

#         # r_tip_in
#         ax.annotate("", xy=(r_tip_in, x_in), xytext=(0, x_in),
#                     arrowprops=dict(arrowstyle="<->", color="#2c7bb6", lw=1.0))
#         ax.text(r_tip_in / 2, x_in - margin * 2.5,
#                 f"$r_{{tip,in}}$={r_tip_in:.3f} m",
#                 ha="center", va="bottom", fontsize=7.5, color="#2c7bb6")

#         # r_out (outlet radius, both hub and tip)
#         ax.annotate("", xy=(r_out, x_tip_out), xytext=(0, x_tip_out),
#                     arrowprops=dict(arrowstyle="<->", color="#555", lw=1.0))
#         ax.text(r_out / 2, x_tip_out - margin,
#                 f"$r_{{out}}$={r_out:.3f} m",
#                 ha="center", va="bottom", fontsize=7.5, color="#555")

#         # Δx axial length — arrow to the right of tip contour
#         r_ann  = max(_r(tip)) + 0.06 * dr
#         ax.annotate("", xy=(r_ann, x_hub_out), xytext=(r_ann, x_in),
#                     arrowprops=dict(arrowstyle="<->", color="#444", lw=1.2))
#         ax.text(r_ann + margin * 0.6, (x_in + x_hub_out) / 2,
#                 f"$\\Delta x$={abs(x_hub_out - x_in):.3f} m",
#                 ha="left", va="center", fontsize=7.5, color="#444")

#     # ── Labels inlet / outlet ─────────────────────────────────────────
#     _m = 0.02 * r_out
#     ax.text((r_hub_in + r_tip_in) / 2, x_in - _m,
#             "Inlet", ha="center", va="bottom",
#             fontsize=8, fontstyle="italic", color="#e67e22")
#     ax.text(r_out + _m, (x_hub_out + x_tip_out) / 2,
#             "Outlet", ha="left", va="center",
#             fontsize=8, fontstyle="italic", color="#27ae60")

#     # ── Axes & légende ────────────────────────────────────────────────
#     ax.set_xlabel("r  — radial  [m]")
#     ax.set_ylabel("x  — axial  [m]")
#     ax.set_title(title, fontsize=12, fontweight="bold")
#     ax.set_xlim(left=0)
#     ax.invert_yaxis()    # inlet en haut, outlet en bas
#     ax.set_aspect("equal")
#     ax.legend(loc="lower right", fontsize=8, framealpha=0.8)
#     ax.grid(True, ls=":", alpha=0.35)
#     fig.tight_layout()
#     return fig


# # ─────────────────────────────────────────────
# #  Exemple d'utilisation
# # ─────────────────────────────────────────────

# if __name__ == "__main__":

#     geo = build_blade(
#         x_le=0.0,  y_le=0.0,
#         x_te=1.0,  y_te=0.0,
#         beta1_deg=55.0,
#         beta2_deg=25.0,
#         d1=0.38,   d2=0.32,
#         t_max_ss=0.09, s_tmax_ss=0.32, t_te_ss=0.004,
#         d1_t_ss=0.12,  d2_t_ss=0.60,
#         t_max_ps=0.07, s_tmax_ps=0.38, t_te_ps=0.004,
#         d1_t_ps=0.12,  d2_t_ps=0.60,
#         n_points=400,
#     )

#     print(f"Corde calculée       : {geo['chord']:.4f} m")
#     print(f"Longueur méridionale : {meridional_length(geo):.4f} m")
#     print(f"Longueur axiale      : {axial_length(geo):.4f} m")

#     PITCH = 0.55

#     # Canal méridional centrifuge : entrée axiale, sortie radiale
#     # Hub et tip partagent r_out (même rayon de sortie), mais des x_out différents
#     channel = build_meridional_channel(
#         x_in=0.0,    r_hub_in=0.003,  r_tip_in=0.011,
#         r_out=0.0218,  x_hub_out=0.010, x_tip_out=0.003,
        
#         alpha_hub_in=44.36,   alpha_hub_out=41.9,
#         alpha_tip_in=44.36,   alpha_tip_out=41.9,
        
#         d_hub_in=0.55, d_hub_out=0.55,
#         d_tip_in=0.45, d_tip_out=0.45,
#         n_points=400,
#     )

#     fig1 = plot_blade(geo, show_camber=True, show_ctrl=True)
#     fig1.savefig("/mnt/user-data/outputs/blade_profile.png",
#                  dpi=150, bbox_inches="tight")

#     fig2 = plot_passage(geo, pitch=PITCH, show_camber=True,
#                         show_ctrl=False, show_lengths=True)
#     fig2.savefig("/mnt/user-data/outputs/blade_passage.png",
#                  dpi=150, bbox_inches="tight")

#     fig3 = plot_meridional(channel, show_ctrl=True, show_annot=True)
#     fig3.savefig("/mnt/user-data/outputs/blade_meridional.png",
#                  dpi=150, bbox_inches="tight")

#     plt.show()
"""
blade_bezier.py
---------------
Reconstruction d'un profil de pale 2D (cartésien x, y) via des courbes de
Bézier cubiques indépendantes pour le pressure side et le suction side.

Paramètres géométriques d'entrée
─────────────────────────────────
  x_le, y_le          : coordonnées du leading edge
  x_te, y_te          : coordonnées du trailing edge
  beta1_deg           : angle de pale à l'entrée  [degrés, mesuré depuis l'axe x]
  beta2_deg           : angle de pale à la sortie [degrés, mesuré depuis l'axe x]
  d1, d2              : longueurs des bras de contrôle aux extrémités (camber)
  t_max_ss            : épaisseur max du suction side  (normalisée par la corde)
  t_max_ps            : épaisseur max du pressure side (normalisée par la corde)
  s_tmax_ss           : position relative (0–1) de t_max sur le SS
  s_tmax_ps           : position relative (0–1) de t_max sur le PS
  t_te_ss, t_te_ps    : demi-épaisseur de bord de fuite SS / PS (normalisée)
  d1_t, d2_t          : bras de contrôle de la Bézier d'épaisseur (entrée/sortie)
  n_points            : nombre de points de discrétisation (défaut 200)

Sorties  (dict retourné par build_blade)
────────────────────────────────────────
  "camber"      – (n, 2) ndarray  points de la camber line
  "ss"          – (n, 2) ndarray  points du suction side
  "ps"          – (n, 2) ndarray  points du pressure side
  "ctrl_camber" – (4, 2) ndarray  points de contrôle de la camber line
  "ctrl_ss"     – (4, 2) ndarray  points de contrôle du suction side
  "ctrl_ps"     – (4, 2) ndarray  points de contrôle du pressure side
  "chord"       – float           longueur de corde LE→TE
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────
#  Évaluation d'une Bézier cubique
# ─────────────────────────────────────────────

def _bezier_cubic(P: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Évalue une Bézier cubique définie par 4 points de contrôle P (4×2)
    pour un vecteur de paramètres t ∈ [0, 1].

    Retourne un tableau (n, 2).
    """
    t = t[:, None]          # (n, 1) pour broadcast
    B = (  (1 - t)**3       * P[0]
         + 3*(1-t)**2 * t   * P[1]
         + 3*(1-t)   * t**2 * P[2]
         +             t**3 * P[3] )
    return B                # (n, 2)


def _bezier_derivative(P: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Dérivée première dB/dt d'une Bézier cubique."""
    t = t[:, None]
    dB = (  3*(1-t)**2       * (P[1] - P[0])
           + 6*(1-t)*t       * (P[2] - P[1])
           + 3*t**2           * (P[3] - P[2]) )
    return dB               # (n, 2)


# ─────────────────────────────────────────────
#  Construction de la camber line
# ─────────────────────────────────────────────

def _build_camber(
    x_le: float, y_le: float,
    x_te: float, y_te: float,
    beta1_deg: float, beta2_deg: float,
    d1: float, d2: float,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Construit la camber line par une Bézier cubique.

    Les tangentes en P0 et P3 sont imposées par les angles de pale β1, β2.
    d1 et d2 contrôlent la longueur des bras de contrôle (degrés de liberté).

    Retourne (points (n,2), ctrl_pts (4,2)).
    """
    beta1 = np.radians(beta1_deg)
    beta2 = np.radians(beta2_deg)

    P0 = np.array([x_le, y_le])
    P3 = np.array([x_te, y_te])

    # Directions tangentes unitaires
    tan1 = np.array([np.cos(beta1), np.sin(beta1)])
    tan2 = np.array([np.cos(beta2), np.sin(beta2)])

    # Points de contrôle intérieurs
    P1 = P0 + d1 * tan1
    P2 = P3 - d2 * tan2          # on remonte depuis TE

    ctrl = np.array([P0, P1, P2, P3])
    t = np.linspace(0.0, 1.0, n)
    pts = _bezier_cubic(ctrl, t)
    return pts, ctrl


# ─────────────────────────────────────────────
#  Distribution d'épaisseur via Bézier cubique
# ─────────────────────────────────────────────

def _thickness_bezier(
    t_le: float,
    t_max: float,
    s_tmax: float,
    t_te: float,
    d1_t: float,
    d2_t: float,
    n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construit une distribution d'épaisseur scalaire e(s), s ∈ [0,1],
    par une Bézier cubique dans l'espace (s, e).

    Points de contrôle :
      Q0 = (0,       t_le)           bord d'attaque
      Q1 = (d1_t,    t_max)          contrôle la pente à l'entrée
      Q2 = (s_tmax,  t_max)          pic d'épaisseur (tangente horizontale)
      Q3 = (1,       t_te)           bord de fuite

    Retourne (s_vals (n,), e_vals (n,), ctrl_pts (4,2)).

    Note : la sortie est en coordonnées (s, e) — à mapper ensuite
           sur la normale à la camber line.
    """
    Q0 = np.array([0.0,     t_le  ])
    Q1 = np.array([d1_t,    t_max ])
    Q2 = np.array([s_tmax,  t_max ])
    Q3 = np.array([1.0,     t_te  ])

    ctrl = np.array([Q0, Q1, Q2, Q3])
    t_param = np.linspace(0.0, 1.0, n)
    pts = _bezier_cubic(ctrl, t_param)   # (n, 2) : (s, e)

    # On ré-échantillonne uniformément en s ∈ [0,1] par interpolation
    s_raw = pts[:, 0]
    e_raw = pts[:, 1]
    s_uniform = np.linspace(0.0, 1.0, n)
    e_uniform = np.interp(s_uniform, s_raw, e_raw)

    return s_uniform, e_uniform, ctrl


# ─────────────────────────────────────────────
#  Projection de l'épaisseur sur la normale
# ─────────────────────────────────────────────

def _offset_curve(
    camber: np.ndarray,
    ctrl_camber: np.ndarray,
    thickness: np.ndarray,
    sign: float,
    chord: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Déplace la camber line de ±thickness (normalisée par la corde)
    dans la direction normale locale.

    sign = +1  → suction side (toward low-pressure)
    sign = -1  → pressure side

    Retourne (offset_pts (n,2), ctrl_approx (4,2)).
    """
    n = camber.shape[0]
    t_param = np.linspace(0.0, 1.0, n)

    # Dérivée tangentielle
    dB = _bezier_derivative(ctrl_camber, t_param)       # (n, 2)
    norms = np.linalg.norm(dB, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1e-12, norms)       # évite /0
    tangents = dB / norms                               # (n, 2) unitaires

    # Normale = rotation 90° de la tangente
    # [cos, sin] → [-sin, cos] (vers suction) ou [sin, -cos] (vers pressure)
    normals = sign * np.column_stack([-tangents[:, 1], tangents[:, 0]])

    # Décalage
    offset = camber + normals * (thickness * chord)[:, None]

    # Points de contrôle approximatifs (LE et TE exacts, intermédiaires estimés)
    ctrl_approx = np.array([
        offset[0],
        offset[n // 3],
        offset[2 * n // 3],
        offset[-1],
    ])

    return offset, ctrl_approx


# ─────────────────────────────────────────────
#  Longueurs caractéristiques
# ─────────────────────────────────────────────

def meridional_length(geo: dict) -> float:
    """
    Longueur méridionale de la pale = longueur d'arc de la camber line.

    Calculée par sommation des segments discrets sur geo["camber"].
    Converge à < 1e-5 de précision relative pour n_points >= 200.

    Returns
    -------
    float : longueur d'arc  [même unité que x, y]
    """
    pts  = geo["camber"]                             # (n, 2)
    segs = np.diff(pts, axis=0)                      # (n-1, 2)
    return float(np.sum(np.linalg.norm(segs, axis=1)))


def axial_length(geo: dict) -> float:
    """
    Longueur axiale de la pale = projection de la camber line sur l'axe x.

    Pour une géométrie cartésienne 2D, c'est simplement x_te - x_le.
    Les points de contrôle intermédiaires n'affectent pas cette valeur
    car la courbe de Bézier est garantie de passer par P0 et P3.

    Returns
    -------
    float : longueur axiale  [même unité que x]
    """
    return float(geo["camber"][-1, 0] - geo["camber"][0, 0])


# ─────────────────────────────────────────────
#  Fonction principale
# ─────────────────────────────────────────────

def build_blade(
    # Géométrie globale
    x_le: float = 0.0,
    y_le: float = 0.0,
    x_te: float = 1.0,
    y_te: float = 0.0,
    # Angles de pale [degrés]
    beta1_deg: float = 50.0,
    beta2_deg: float = 20.0,
    # Bras de contrôle camber line
    d1: float = 0.35,
    d2: float = 0.35,
    # Épaisseur suction side
    t_max_ss: float = 0.08,
    s_tmax_ss: float = 0.35,
    t_te_ss: float = 0.005,
    d1_t_ss: float = 0.15,
    d2_t_ss: float = 0.65,
    # Épaisseur pressure side
    t_max_ps: float = 0.06,
    s_tmax_ps: float = 0.40,
    t_te_ps: float = 0.005,
    d1_t_ps: float = 0.15,
    d2_t_ps: float = 0.65,
    # Discrétisation
    n_points: int = 300,
) -> dict:
    """
    Reconstruit un profil de pale 2D cartésien via des Bézier cubiques.

    Paramètres
    ----------
    x_le, y_le      : coordonnées du leading edge
    x_te, y_te      : coordonnées du trailing edge
    beta1_deg       : angle de pale à l'entrée [°] depuis l'axe x
    beta2_deg       : angle de pale à la sortie [°] depuis l'axe x
    d1, d2          : longueurs des bras de contrôle de la camber line
                      (normalisées par la longueur LE→TE, typiquement 0.2–0.5)
    t_max_ss/ps     : épaisseur maximale SS/PS normalisée par la corde
    s_tmax_ss/ps    : position curviligne (0–1) du max d'épaisseur SS/PS
    t_te_ss/ps      : demi-épaisseur au bord de fuite SS/PS (normalisée)
    d1_t_ss/ps      : bras de contrôle Bézier épaisseur côté LE (0–1)
    d2_t_ss/ps      : bras de contrôle Bézier épaisseur côté TE (0–1)
    n_points        : nombre de points de discrétisation

    Retourne
    --------
    dict avec les clés :
        "camber"      – (n, 2) ndarray
        "ss"          – (n, 2) ndarray  suction side
        "ps"          – (n, 2) ndarray  pressure side
        "ctrl_camber" – (4, 2) ndarray
        "ctrl_ss"     – (4, 2) ndarray
        "ctrl_ps"     – (4, 2) ndarray
        "chord"       – float
    """
    # 1. Camber line
    camber, ctrl_camber = _build_camber(
        x_le, y_le, x_te, y_te,
        beta1_deg, beta2_deg,
        d1, d2, n_points,
    )

    # Corde réelle (distance LE–TE)
    chord = float(np.linalg.norm(np.array([x_te - x_le, y_te - y_le])))

    # 2. Distributions d'épaisseur
    _, e_ss, ctrl_ss = _thickness_bezier(
        t_le=0.0, t_max=t_max_ss, s_tmax=s_tmax_ss,
        t_te=t_te_ss, d1_t=d1_t_ss, d2_t=d2_t_ss, n=n_points,
    )
    _, e_ps, ctrl_ps = _thickness_bezier(
        t_le=0.0, t_max=t_max_ps, s_tmax=s_tmax_ps,
        t_te=t_te_ps, d1_t=d1_t_ps, d2_t=d2_t_ps, n=n_points,
    )

    # 3. Projection normale sur la camber line
    ss, ctrl_ss = _offset_curve(camber, ctrl_camber, e_ss, +1.0, chord)
    ps, ctrl_ps = _offset_curve(camber, ctrl_camber, e_ps, -1.0, chord)

    return {
        "camber":      camber,
        "ss":          ss,
        "ps":          ps,
        "ctrl_camber": ctrl_camber,
        "ctrl_ss":     ctrl_ss,
        "ctrl_ps":     ctrl_ps,
        "chord":       chord,
    }


# ─────────────────────────────────────────────
#  Visualisation
# ─────────────────────────────────────────────

def plot_blade(
    geo: dict,
    show_camber: bool = True,
    show_ctrl: bool = True,
    title: str = "Profil de pale — Bézier cubique",
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """
    Trace le profil complet (SS, PS, camber line, points de contrôle).

    Paramètres
    ----------
    geo         : dict retourné par build_blade()
    show_camber : afficher la camber line
    show_ctrl   : afficher les polygones de contrôle Bézier
    title       : titre du graphique
    ax          : axes matplotlib existants (optionnel)

    Retourne
    --------
    Figure matplotlib.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 5))
    else:
        fig = ax.get_figure()

    camber      = geo["camber"]
    ss          = geo["ss"]
    ps          = geo["ps"]
    ctrl_camber = geo["ctrl_camber"]
    ctrl_ss     = geo["ctrl_ss"]
    ctrl_ps     = geo["ctrl_ps"]

    # ── Remplissage du profil ──────────────────────────────────────────
    contour_x = np.concatenate([ss[:, 0], ps[::-1, 0]])
    contour_y = np.concatenate([ss[:, 1], ps[::-1, 1]])
    ax.fill(contour_x, contour_y, color="#d0e8f5", alpha=0.6, zorder=1)

    # ── Suction side ──────────────────────────────────────────────────
    ax.plot(ss[:, 0], ss[:, 1],
            color="#1a6fa8", lw=2.0, label="Suction side", zorder=3)

    # ── Pressure side ─────────────────────────────────────────────────
    ax.plot(ps[:, 0], ps[:, 1],
            color="#c0392b", lw=2.0, label="Pressure side", zorder=3)

    # ── Camber line ───────────────────────────────────────────────────
    if show_camber:
        ax.plot(camber[:, 0], camber[:, 1],
                color="#2ecc71", lw=1.4, ls="--", label="Camber line", zorder=4)

    # ── Polygones de contrôle ─────────────────────────────────────────
    if show_ctrl:
        for ctrl, color, label in [
            (ctrl_camber, "#2ecc71", "Ctrl camber"),
            (ctrl_ss,     "#1a6fa8", "Ctrl SS"),
            (ctrl_ps,     "#c0392b", "Ctrl PS"),
        ]:
            ax.plot(ctrl[:, 0], ctrl[:, 1],
                    "o--", color=color, lw=0.8, ms=5, alpha=0.6,
                    label=label, zorder=5)

    # ── Marqueurs LE / TE ─────────────────────────────────────────────
    le = camber[0]
    te = camber[-1]
    ax.plot(*le, "ko", ms=6, zorder=6)
    ax.plot(*te, "ks", ms=6, zorder=6)
    ax.annotate("LE", le, textcoords="offset points", xytext=(-14, 6), fontsize=9)
    ax.annotate("TE", te, textcoords="offset points", xytext=(5, 6),   fontsize=9)

    # ── Axes & légende ────────────────────────────────────────────────
    ax.set_aspect("equal")
    ax.set_xlabel("x  [m]")
    ax.set_ylabel("y  [m]")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)
    ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────
#  Visualisation du passage
# ─────────────────────────────────────────────

def plot_passage(
    geo: dict,
    pitch: float,
    show_camber: bool = True,
    show_ctrl: bool = False,
    show_lengths: bool = True,
    title: str = "Passage de pale — Bézier cubique",
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """
    Trace le passage complet entre deux pales consécutives.

    La pale de référence est celle de geo. La pale voisine est obtenue
    par translation de +pitch selon y (direction tangentielle).
    Le passage fluide entre le SS de la pale de référence et le PS
    de la pale voisine est colorié.

    Paramètres
    ----------
    geo          : dict retourné par build_blade()
    pitch        : espacement entre pales [même unité que x, y]
    show_camber  : afficher les camber lines
    show_ctrl    : afficher les polygones de contrôle Bézier
    show_lengths : annoter les longueurs méridionale et axiale
    title        : titre du graphique
    ax           : axes matplotlib existants (optionnel)

    Retourne
    --------
    Figure matplotlib.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 7))
    else:
        fig = ax.get_figure()

    camber = geo["camber"]
    ss     = geo["ss"]
    ps     = geo["ps"]

    # Pale voisine : translation de +pitch en y
    ss_next     = ss     + np.array([0.0, pitch])
    ps_next     = ps     + np.array([0.0, pitch])
    camber_next = camber + np.array([0.0, pitch])

    # ── Passage fluide ────────────────────────────────────────────────
    # Délimité par : SS de la pale ref (bas) + PS de la pale voisine (haut)
    # Contour : SS LE→TE, puis PS_next TE→LE
    pass_x = np.concatenate([ss[:, 0],      ps_next[::-1, 0]])
    pass_y = np.concatenate([ss[:, 1],      ps_next[::-1, 1]])
    ax.fill(pass_x, pass_y, color="#e8f4e8", alpha=0.7, zorder=1,
            label="Passage fluide")

    # ── Pale de référence ─────────────────────────────────────────────
    blade_x = np.concatenate([ss[:, 0], ps[::-1, 0]])
    blade_y = np.concatenate([ss[:, 1], ps[::-1, 1]])
    ax.fill(blade_x, blade_y, color="#d0e8f5", alpha=0.8, zorder=2)
    ax.plot(ss[:, 0], ss[:, 1], color="#1a6fa8", lw=2.0,
            label="Suction side", zorder=4)
    ax.plot(ps[:, 0], ps[:, 1], color="#c0392b", lw=2.0,
            label="Pressure side", zorder=4)

    # ── Pale voisine ──────────────────────────────────────────────────
    blade_next_x = np.concatenate([ss_next[:, 0], ps_next[::-1, 0]])
    blade_next_y = np.concatenate([ss_next[:, 1], ps_next[::-1, 1]])
    ax.fill(blade_next_x, blade_next_y, color="#d0e8f5", alpha=0.8, zorder=2)
    ax.plot(ss_next[:, 0], ss_next[:, 1], color="#1a6fa8", lw=2.0,
            zorder=4)
    ax.plot(ps_next[:, 0], ps_next[:, 1], color="#c0392b", lw=2.0,
            zorder=4)

    # ── Camber lines ──────────────────────────────────────────────────
    if show_camber:
        ax.plot(camber[:, 0],      camber[:, 1],      color="#2ecc71",
                lw=1.4, ls="--", label="Camber line", zorder=5)
        ax.plot(camber_next[:, 0], camber_next[:, 1], color="#2ecc71",
                lw=1.4, ls="--", zorder=5)

    # ── Polygones de contrôle ─────────────────────────────────────────
    if show_ctrl:
        for ctrl, color in [
            (geo["ctrl_camber"],              "#2ecc71"),
            (geo["ctrl_ss"],                  "#1a6fa8"),
            (geo["ctrl_ps"],                  "#c0392b"),
            (geo["ctrl_camber"] + [0, pitch], "#2ecc71"),
            (geo["ctrl_ss"]     + [0, pitch], "#1a6fa8"),
            (geo["ctrl_ps"]     + [0, pitch], "#c0392b"),
        ]:
            ax.plot(ctrl[:, 0], ctrl[:, 1], "o--", color=color,
                    lw=0.8, ms=4, alpha=0.5, zorder=6)

    # ── Marqueurs LE / TE ─────────────────────────────────────────────
    for cam, offset_label in [(camber, ""), (camber_next, "")]:
        le, te = cam[0], cam[-1]
        ax.plot(*le, "ko", ms=5, zorder=7)
        ax.plot(*te, "ks", ms=5, zorder=7)
    ax.annotate("LE", camber[0], textcoords="offset points",
                xytext=(-18, 4), fontsize=8)
    ax.annotate("TE", camber[-1], textcoords="offset points",
                xytext=(5, 4), fontsize=8)

    # ── Annotations de longueur ───────────────────────────────────────
    if show_lengths:
        L_mer = meridional_length(geo)
        L_ax  = axial_length(geo)

        # Flèche longueur axiale (sous la pale)
        y_arrow = ps[:, 1].min() - 0.04 * pitch
        x0, x1  = camber[0, 0], camber[-1, 0]
        ax.annotate("", xy=(x1, y_arrow), xytext=(x0, y_arrow),
                    arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2))
        ax.text((x0 + x1) / 2, y_arrow - 0.03 * pitch,
                f"$L_{{ax}}$ = {L_ax:.3f} m",
                ha="center", va="top", fontsize=8, color="#555")

        # Flèche longueur méridionale (le long de la camber)
        # Positionnée à mi-passage en y, décalée à gauche
        x_mid  = camber[:, 0].mean()
        y_mid  = camber[:, 1].mean()
        ax.annotate(f"$L_{{mer}}$ = {L_mer:.3f} m",
                    xy=(x_mid, y_mid),
                    xytext=(x_mid - 0.12 * L_ax, y_mid + 0.08 * pitch),
                    fontsize=8, color="#2ecc71",
                    arrowprops=dict(arrowstyle="->", color="#2ecc71", lw=1.0))

        # Flèche pitch (à droite)
        x_pitch = camber[:, 0].max() + 0.05 * L_ax
        y_bot   = ps[0, 1]
        y_top   = ps_next[0, 1]
        ax.annotate("", xy=(x_pitch, y_top), xytext=(x_pitch, y_bot),
                    arrowprops=dict(arrowstyle="<->", color="#888", lw=1.2))
        ax.text(x_pitch + 0.02 * L_ax, (y_bot + y_top) / 2,
                f"pitch = {pitch:.3f} m",
                ha="left", va="center", fontsize=8, color="#888")

    # ── Axes & légende ────────────────────────────────────────────────
    ax.set_aspect("equal")
    ax.set_xlabel("x  [m]")
    ax.set_ylabel("y  [m]")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.8)
    ax.grid(True, ls=":", alpha=0.35)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────
#  Canal méridional (hub & tip)
# ─────────────────────────────────────────────

def _clip_contour_at_r(pts: np.ndarray, r_target: float) -> tuple[np.ndarray, float]:
    """
    Clip a contour (n, 2) in (x, r) at r = r_target.
    Returns the clipped points and the x value at the crossing.
    Assumes r is monotonically increasing along the contour.
    """
    r = pts[:, 1]
    x = pts[:, 0]
    # Find crossing index
    idx = np.searchsorted(r, r_target)
    idx = np.clip(idx, 1, len(r) - 1)
    # Linear interpolation between idx-1 and idx
    r0, r1 = r[idx - 1], r[idx]
    x0, x1 = x[idx - 1], x[idx]
    frac = (r_target - r0) / (r1 - r0) if (r1 - r0) > 1e-14 else 0.0
    x_cross = float(x0 + frac * (x1 - x0))
    clipped = np.vstack([pts[:idx], [[x_cross, r_target]]])
    return clipped, x_cross


def build_meridional_channel(
    # Inlet — axial plane : x = x_in (constant), r ∈ [r_hub_in, r_tip_in]
    x_in:       float = 0.0,
    r_hub_in:   float = 0.03,
    r_tip_in:   float = 0.10,
    # Outlet — radial plane : r = r_out (constant), x found by intersection
    r_out:      float = 0.25,
    # Axial extent used to size each Bézier before clipping at r_out.
    # Must be large enough that both contours actually reach r_out.
    x_max:      float = 0.15,
    # Tangent angles [degrees from x-axis] at inlet and outlet
    alpha_hub_in:  float = 5.0,
    alpha_hub_out: float = 88.0,
    alpha_tip_in:  float = 3.0,
    alpha_tip_out: float = 85.0,
    # Control arm lengths (× x_max, applied uniformly to size the Bézier)
    d_hub_in:  float = 0.55,
    d_hub_out: float = 0.55,
    d_tip_in:  float = 0.45,
    d_tip_out: float = 0.45,
    # Discretisation
    n_points: int = 600,
) -> dict:
    """
    Construit les contours hub et tip d'un impeller centrifuge dans le plan
    méridional (x, r), puis les clip à r = r_out.

    Workflow
    ────────
    1. Chaque contour est une Bézier cubique partant de (x_in, r_hub/tip_in)
       et ciblant (x_max, r_out) — x_max sert uniquement à dimensionner les
       bras de contrôle et n'apparaît pas dans le résultat final.
    2. La courbe est discrétisée puis clippée à r = r_out par interpolation
       linéaire → donne x_hub_out et x_tip_out automatiquement.

    Paramètres
    ----------
    x_in                : abscisse axiale de l'inlet [m]
    r_hub_in, r_tip_in  : rayons hub / tip à l'inlet [m]
    r_out               : rayon de l'outlet commun [m]
    x_max               : étendue axiale max pour dimensionner la Bézier [m]
                          (doit être > x_hub_out réel ; typiquement 1.5–2× Δr)
    alpha_hub/tip_in    : angle de tangente à l'inlet [°]  (0° = axial pur)
    alpha_hub/tip_out   : angle de tangente à l'outlet [°] (90° = radial pur)
    d_hub/tip_in/out    : longueurs des bras de contrôle (× x_max)
    n_points            : discrétisation de la Bézier avant clipping

    Retourne
    --------
    dict :
        "hub"        – (n, 2)  contour hub clippé dans (x, r)
        "tip"        – (n, 2)  contour tip clippé dans (x, r)
        "ctrl_hub"   – (4, 2)  points de contrôle Bézier hub (avant clipping)
        "ctrl_tip"   – (4, 2)  points de contrôle Bézier tip (avant clipping)
        "x_in"       – float
        "r_hub_in"   – float
        "r_tip_in"   – float
        "r_out"      – float
        "x_hub_out"  – float   ← calculé par intersection
        "x_tip_out"  – float   ← calculé par intersection
    """
    def _make_contour(r_start, a_in_deg, a_out_deg, d_in, d_out):
        # The Bézier runs from (x_in, r_start) toward a virtual endpoint
        # placed well beyond r_out so the curve naturally crosses r_out
        # at an interior point, making the clipping meaningful.
        # The virtual endpoint is placed at (x_max, r_out * 1.3) so that
        # the curve always overshoots r_out before x_max.
        r_virtual = r_out * 1.3
        P0 = np.array([x_in,  r_start  ])
        P3 = np.array([x_max, r_virtual])
        scale = x_max - x_in
        tan_in  = np.array([np.cos(np.radians(a_in_deg)),
                             np.sin(np.radians(a_in_deg))])
        tan_out = np.array([np.cos(np.radians(a_out_deg)),
                             np.sin(np.radians(a_out_deg))])
        P1 = P0 + d_in  * scale * tan_in
        P2 = P3 - d_out * scale * tan_out
        ctrl = np.array([P0, P1, P2, P3])
        t    = np.linspace(0.0, 1.0, n_points)
        pts  = _bezier_cubic(ctrl, t)
        return pts, ctrl

    hub_full, ctrl_hub = _make_contour(
        r_hub_in, alpha_hub_in, alpha_hub_out, d_hub_in, d_hub_out)
    tip_full, ctrl_tip = _make_contour(
        r_tip_in, alpha_tip_in, alpha_tip_out, d_tip_in, d_tip_out)

    # Clip both contours at r = r_out
    hub, x_hub_out = _clip_contour_at_r(hub_full, r_out)
    tip, x_tip_out = _clip_contour_at_r(tip_full, r_out)

    return {
        "hub":       hub,
        "tip":       tip,
        "ctrl_hub":  ctrl_hub,
        "ctrl_tip":  ctrl_tip,
        "x_in":      x_in,
        "r_hub_in":  r_hub_in,
        "r_tip_in":  r_tip_in,
        "r_out":     r_out,
        "x_hub_out": x_hub_out,
        "x_tip_out": x_tip_out,
    }


def plot_meridional(
    channel: dict,
    show_ctrl:  bool = True,
    show_annot: bool = True,
    title: str = "Vue méridionale — Hub & Tip",
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """
    Vue méridionale (r horizontal, x axial inversé) d'un impeller centrifuge.

    Inlet  : segment vertical   r ∈ [r_hub_in, r_tip_in] à x = x_in
    Outlet : segment horizontal x ∈ [x_tip_out, x_hub_out] à r = r_out
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    else:
        fig = ax.get_figure()

    hub      = channel["hub"]        # (n,2) en (x, r)
    tip      = channel["tip"]
    ctrl_hub = channel["ctrl_hub"]
    ctrl_tip = channel["ctrl_tip"]

    x_in      = channel["x_in"]
    r_hub_in  = channel["r_hub_in"]
    r_tip_in  = channel["r_tip_in"]
    r_out     = channel["r_out"]
    x_hub_out = channel["x_hub_out"]
    x_tip_out = channel["x_tip_out"]

    def _r(pts): return pts[:, 1]
    def _x(pts): return pts[:, 0]

    # ── Passage fluide ────────────────────────────────────────────────
    # Contour fermé dans le plan (r, x) :
    #   hub LE→TE  + outlet (x_hub_out→x_tip_out à r=r_out) + tip TE→LE
    outlet_r = np.array([r_out, r_out])
    outlet_x = np.array([x_hub_out, x_tip_out])

    pass_r = np.concatenate([_r(hub), outlet_r, _r(tip)[::-1]])
    pass_x = np.concatenate([_x(hub), outlet_x, _x(tip)[::-1]])
    ax.fill(pass_r, pass_x, color="#dceefb", alpha=0.7, zorder=1,
            label="Passage fluide")
    ax.fill(pass_r, pass_x, hatch="///", facecolor="none",
            edgecolor="#aac8e0", lw=0.4, alpha=0.4, zorder=2)

    # ── Contours hub / tip ────────────────────────────────────────────
    ax.plot(_r(hub), _x(hub), color="#8B4513", lw=2.5, label="Hub",          zorder=4)
    ax.plot(_r(tip), _x(tip), color="#2c7bb6", lw=2.5, label="Tip (shroud)", zorder=4)

    # ── Inlet : vertical (r varie, x constant) ────────────────────────
    ax.plot([r_hub_in, r_tip_in], [x_in, x_in],
            color="#e67e22", lw=2.0, ls="--", label="Inlet", zorder=3)

    # ── Outlet : horizontal (x varie, r constant) ─────────────────────
    ax.plot([r_out, r_out], [x_hub_out, x_tip_out],
            color="#27ae60", lw=2.0, ls="--", label="Outlet", zorder=3)

    # ── Axe de rotation ───────────────────────────────────────────────
    ax.axvline(0, color="#ccc", lw=0.8, ls=":", zorder=0)
    ax.text(0.002, x_in, "axe", ha="left", va="center", fontsize=7, color="#aaa")

    # ── Polygones de contrôle ─────────────────────────────────────────
    if show_ctrl:
        for ctrl, color in [(ctrl_hub, "#8B4513"), (ctrl_tip, "#2c7bb6")]:
            ax.plot(ctrl[:, 1], ctrl[:, 0], "o--", color=color,
                    lw=0.9, ms=5, alpha=0.55, zorder=5)

    # ── Annotations ───────────────────────────────────────────────────
    if show_annot:
        dr     = r_out
        margin = 0.03 * dr

        # r_hub_in
        ax.annotate("", xy=(r_hub_in, x_in), xytext=(0, x_in),
                    arrowprops=dict(arrowstyle="<->", color="#8B4513", lw=1.0))
        ax.text(r_hub_in / 2, x_in - margin,
                f"$r_{{hub,in}}$={r_hub_in:.3f} m",
                ha="center", va="bottom", fontsize=7.5, color="#8B4513")

        # r_tip_in
        ax.annotate("", xy=(r_tip_in, x_in), xytext=(0, x_in),
                    arrowprops=dict(arrowstyle="<->", color="#2c7bb6", lw=1.0))
        ax.text(r_tip_in / 2, x_in - margin * 2.5,
                f"$r_{{tip,in}}$={r_tip_in:.3f} m",
                ha="center", va="bottom", fontsize=7.5, color="#2c7bb6")

        # r_out (outlet radius, both hub and tip)
        ax.annotate("", xy=(r_out, x_tip_out), xytext=(0, x_tip_out),
                    arrowprops=dict(arrowstyle="<->", color="#555", lw=1.0))
        ax.text(r_out / 2, x_tip_out - margin,
                f"$r_{{out}}$={r_out:.3f} m",
                ha="center", va="bottom", fontsize=7.5, color="#555")

        # Δx axial length — arrow to the right of tip contour
        r_ann  = max(_r(tip)) + 0.06 * dr
        ax.annotate("", xy=(r_ann, x_hub_out), xytext=(r_ann, x_in),
                    arrowprops=dict(arrowstyle="<->", color="#444", lw=1.2))
        ax.text(r_ann + margin * 0.6, (x_in + x_hub_out) / 2,
                f"$\\Delta x$={abs(x_hub_out - x_in):.3f} m",
                ha="left", va="center", fontsize=7.5, color="#444")

    # ── Labels inlet / outlet ─────────────────────────────────────────
    _m = 0.02 * r_out
    ax.text((r_hub_in + r_tip_in) / 2, x_in - _m,
            "Inlet", ha="center", va="bottom",
            fontsize=8, fontstyle="italic", color="#e67e22")
    ax.text(r_out + _m, (x_hub_out + x_tip_out) / 2,
            "Outlet", ha="left", va="center",
            fontsize=8, fontstyle="italic", color="#27ae60")

    # ── Axes & légende ────────────────────────────────────────────────
    ax.set_xlabel("r  — radial  [m]")
    ax.set_ylabel("x  — axial  [m]")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(left=0)
    ax.invert_yaxis()    # inlet en haut, outlet en bas
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.8)
    ax.grid(True, ls=":", alpha=0.35)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────
#  Exemple d'utilisation
# ─────────────────────────────────────────────

if __name__ == "__main__":

    geo = build_blade(
        x_le=0.0,  y_le=0.0,
        x_te=1.0,  y_te=0.0,
        beta1_deg=55.0,
        beta2_deg=25.0,
        d1=0.38,   d2=0.32,
        t_max_ss=0.09, s_tmax_ss=0.32, t_te_ss=0.004,
        d1_t_ss=0.12,  d2_t_ss=0.60,
        t_max_ps=0.07, s_tmax_ps=0.38, t_te_ps=0.004,
        d1_t_ps=0.12,  d2_t_ps=0.60,
        n_points=400,
    )

    print(f"Corde calculée       : {geo['chord']:.4f} m")
    print(f"Longueur méridionale : {meridional_length(geo):.4f} m")
    print(f"Longueur axiale      : {axial_length(geo):.4f} m")

    PITCH = 0.55

    channel = build_meridional_channel(
        x_in=0.0,   r_hub_in=0.003, r_tip_in=0.012,
        r_out=0.0293, x_max=0.011,
        alpha_hub_in=5,  alpha_hub_out=88,
        alpha_tip_in=3,  alpha_tip_out=85,
        d_hub_in=0.55, d_hub_out=0.55,
        d_tip_in=0.45, d_tip_out=0.45,
    )
    
    #     # Hub et tip partagent r_out (même rayon de sortie), mais des x_out différents
    #     channel = build_meridional_channel(
    #         x_in=0.0,    r_hub_in=0.003,  r_tip_in=0.011,
    #         r_out=0.0218,  x_hub_out=0.010, x_tip_out=0.003,
            
    #         alpha_hub_in=44.36,   alpha_hub_out=41.9,
    #         alpha_tip_in=44.36,   alpha_tip_out=41.9,
            
    #         d_hub_in=0.55, d_hub_out=0.55,
    #         d_tip_in=0.45, d_tip_out=0.45,
    #         n_points=400,
    #     )
    
    print(f"x_hub_out (computed) : {channel['x_hub_out']:.4f} m")
    print(f"x_tip_out (computed) : {channel['x_tip_out']:.4f} m")

    fig1 = plot_blade(geo, show_camber=True, show_ctrl=True)
    fig1.savefig("/mnt/user-data/outputs/blade_profile.png",
                 dpi=150, bbox_inches="tight")

    fig2 = plot_passage(geo, pitch=PITCH, show_camber=True,
                        show_ctrl=False, show_lengths=True)
    fig2.savefig("/mnt/user-data/outputs/blade_passage.png",
                 dpi=150, bbox_inches="tight")

    fig3 = plot_meridional(channel, show_ctrl=True, show_annot=True)
    fig3.savefig("/mnt/user-data/outputs/blade_meridional.png",
                 dpi=150, bbox_inches="tight")

    plt.show()
