"""
Author : Basile Chaudoir

Sizing par optimisation (PSO) d'un échangeur à canaux spiralés empilés (SCHE),
calqué sur SPHESizingOpt. Le modèle thermo-hydraulique est
HexMB('SCHE') (moving boundary, discrétisé), avec la résolution géométrique de
hex_geometry.py (SCHEGeometry) et les corrélations curved_pipe_HTC / curved_pipe_DP.

Noms de corrélations (clés '1P' / '2P') :
    HTC 1P : 'Sigalotti' (tube courbe), ou un nom de corrélation de tube droit ;
    HTC 2P : chaud (condensation) 'Zakeralhoseini' ; froid (évaporation) 'Zhao', 'Nariai' (None : auto par la turbulence)
    DP  1P : 'Sigalotti', 'Ito' ; DP 2P : corrélations de tubes ('Friedel', 'MSH', 'Choi', ...; None : auto)
             ou 'Cioncolini', 'Ruffell', 'Fsadni'

Géométrie (cf. setup_geom, branche SCHE) :
    - canaux elliptiques complets, largeur w_ch_*, hauteur h_ch_*
    - P_H / P_V : écarts centre à centre entre deux canaux du MÊME fluide
    - chaud et froid en quinconce, décalés de (P_H/2, P_V/2)
    - n_ch = floor(H / P_V) canaux de chaque fluide en parallèle

Deux modes d'objectif (paramètre obj de sizing(), ou params['obj']) :
    obj = 'V' : V_tot = n_parallel * n_series * pi * D_ext^2 / 4 * H      [m^3]
    obj = 'A' : A_tot = n_parallel * n_series * A_unit                    [m^2]
                A_unit = HX.A_h / n_parallel (surface côté chaud d'UNE unité) si disponible, sinon estimation
sous contraintes de pertes de charge (DP_h, DP_c) et, au choix :
    - SH_min : surchauffe minimale en sortie du côté froid (fluide qui s'évapore) [K] ;
    - Q_dot  : puissance minimale [W].
Les deux peuvent être actives en même temps ; None désactive une contrainte.

Surchauffe en sortie : SH = T_out - T_dew(p_out) si h_out >= h_dew.
Si la sortie est encore diphasique, SH = (h_out - h_dew) / cp_v,dew (< 0) :
c'est une « surchauffe équivalente » qui varie continûment avec h_out et
guide l'essaim, alors que T_out - T_dew resterait à 0 pour un corps pur.

Variables de design (ORDER) :
    D_ext       diamètre extérieur                         [m]
    H           hauteur de l'empilement                    [m]
    w_ch_C   largeur du canal froid                     [m]
    h_ch_C   hauteur du canal froid                     [m]
    w_ch_H   largeur du canal chaud                     [m]
    h_ch_H   hauteur du canal chaud                     [m]
    P_H     pas radial (même fluide, centre à centre)  [m]
    P_V     pas vertical (même fluide, centre à centre)[m]
    r0_in       rayon intérieur de départ de la spirale    [m]
    n_parallel  nombre d'unités en parallèle               [-] (entier)
    n_series    nombre d'unités en série                   [-] (entier)

Une variable dont les bornes sont égales (lb == ub) est FIGÉE.
Rapport d'aspect imposé : params C_aspect / H_aspect = h / w (ex. 1.0 = canal
circulaire). La hauteur correspondante n'est alors plus une variable : h = aspect * w.

Pas imposés par la tenue en pression (pitch_mode = 'pressure', défaut) :
    P_H et P_V ne sont plus des variables : ils sont calculés pour chaque
    candidat comme les pas MINIMAUX qui tiennent la pression.
    - chaque canal est assimilé à un tube de diamètre intérieur max(w, h) ;
      son épaisseur requise t_x est la plus petite épaisseur telle que
          Internal_Max_P_carbon_steel(D_o, t, T) >= P_x       (pression du fluide)
          External_Max_P_carbon_steel(D_o, t, T) >= P_autre   (canal vide, autre fluide en charge)
      (labothappy.toolbox.piping.pipe_thickness, D_o = max(w, h) + 2t)
    - chaque canal garde sa propre paroi (pas de paroi partagée) :
          jeu même fluide  >= max(2 t_x,      min_wall)
          paroi chaud-froid >= max(t_C + t_H, min_wall)   (diagonale)
    - pitch_H0 = max(C_w + jeu_C, H_w + jeu_H), pitch_V0 = idem avec les hauteurs,
      puis les deux pas sont multipliés par le même facteur s >= 1 (bisection)
      jusqu'à respecter la paroi diagonale.
    P de calcul = P_design_factor * P_su ; T de calcul = T_design ou max des T d'entrée.
    pitch_mode = 'free' : ancien comportement (P_H, P_V optimisés).

Faisabilité géométrique vérifiée AVANT le calcul HX (pénalité graduée) :
    - jeu entre deux canaux du même fluide >= min_wall (horizontal et vertical)
    - paroi minimale entre canaux chaud et froid voisins (diagonale) >= min_wall
    - au moins min_turns tours, au moins une couche de canaux
"""

import numpy as np
import warnings
from contextlib import contextmanager

import pyswarms as ps
import joblib
from joblib import Parallel, delayed
from joblib.parallel import BatchCompletionCallBack
from tqdm import tqdm

from labothappy.component.base_component import BaseComponent
from labothappy.component.heat_exchanger.hex_MB import HexMB
from labothappy.toolbox.piping.pipe_thickness import Internal_Max_P_carbon_steel, External_Max_P_carbon_steel

warnings.filterwarnings('ignore')


@contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar."""
    class TqdmBatchCompletionCallback(BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_cb = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_cb
        tqdm_object.close()


#%%

def _ellipse(w, h):
    """Section, périmètre (Ramanujan) et diamètre hydraulique d'une ellipse complète w x h."""
    a_, b_ = w / 2.0, h / 2.0
    A = np.pi * a_ * b_
    Per = np.pi * (3*(a_ + b_) - np.sqrt((3*a_ + b_) * (a_ + 3*b_)))
    return A, Per, 4.0 * A / Per


def _diag_wall(Cw, Ch, Hw, Hh, p, pv, n_phi=181):
    """
    Paroi minimale entre un canal froid centré en (0, 0) et un canal chaud
    centré en (p/2, pv/2). Retourne (t_wall_min, chevauchement).
    """
    phi = np.linspace(0.0, 2*np.pi, n_phi)[:-1]
    xc, zc = 0.5*Cw*np.cos(phi), 0.5*Ch*np.sin(phi)
    xh, zh = 0.5*p + 0.5*Hw*np.cos(phi), 0.5*pv + 0.5*Hh*np.sin(phi)
    overlap = (((xc - 0.5*p)/(0.5*Hw))**2 + ((zc - 0.5*pv)/(0.5*Hh))**2 < 1).any() \
              or ((xh/(0.5*Cw))**2 + (zh/(0.5*Ch))**2 < 1).any()
    d = np.sqrt((xc[:, None] - xh[None, :])**2 + (zc[:, None] - zh[None, :])**2).min()
    return (0.0 if overlap else float(d)), bool(overlap)


def outlet_superheat(HX, fluid):
    """
    Surchauffe (équivalente) en sortie du côté froid (fluide qui s'évapore) [K].
    Convention du modèle : le froid entre en k = 0 et sort en k = -1.
    """
    import CoolProp.CoolProp as CP
    h_out = float(np.asarray(HX.hvec_c)[-1])
    p_out = float(np.asarray(HX.pvec_c)[-1])
    T_out = float(np.asarray(HX.Tvec_c)[-1])

    T_dew = CP.PropsSI('T', 'P', p_out, 'Q', 1, fluid)
    h_dew = CP.PropsSI('H', 'P', p_out, 'Q', 1, fluid)
    if h_out >= h_dew:
        return T_out - T_dew
    cp_v = CP.PropsSI('C', 'P', p_out, 'Q', 1, fluid)
    return (h_out - h_dew) / cp_v


_T_ASME = (273.15, 426.667 + 273.15)     # plage du tableau de contraintes admissibles


def required_thickness(D_in, p_int, p_ext, T, n_iter=40):
    """
    Épaisseur minimale d'un tube de diamètre intérieur D_in tenant la pression
    interne p_int et la pression externe p_ext à la température T [K]
    (fonctions de labothappy.toolbox.piping.pipe_thickness, D_o = D_in + 2t).
    Retourne np.inf si aucune épaisseur <= 10 D_in ne convient.
    """
    T = float(np.clip(T, *_T_ASME))      # hors plage : contrainte du bord le plus proche

    def ok(t):
        D_o = D_in + 2.0*t
        return (Internal_Max_P_carbon_steel(D_o, t, T) >= p_int and
                External_Max_P_carbon_steel(D_o, t, T) >= p_ext)

    hi = 0.05 * D_in
    while not ok(hi):
        hi *= 2.0
        if hi > 10.0 * D_in:
            return np.inf
    lo = 0.0
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi


def pressure_pitches(w_ch_C, h_ch_C, w_ch_H, h_ch_H, P_h, P_c, T, min_wall=0.0):
    """
    Pas minimaux (P_H, P_V) tenant la pression, et épaisseurs requises (t_C, t_H).
    Voir l'en-tête du module pour la méthode. Retourne des np.inf si la pression
    ne peut pas être tenue.
    """
    Cw, Ch, Hw, Hh = w_ch_C, h_ch_C, w_ch_H, h_ch_H
    t_C = required_thickness(max(Cw, Ch), P_c, P_h, T)
    t_H = required_thickness(max(Hw, Hh), P_h, P_c, T)
    if not (np.isfinite(t_C) and np.isfinite(t_H)):
        return np.inf, np.inf, t_C, t_H

    g_C, g_H = max(2*t_C, min_wall), max(2*t_H, min_wall)
    g_diag = max(t_C + t_H, min_wall)
    pH0 = max(Cw + g_C, Hw + g_H)
    pV0 = max(Ch + g_C, Hh + g_H)

    def diag_ok(s):
        t_w, ov = _diag_wall(Cw, Ch, Hw, Hh, s*pH0, s*pV0, n_phi=121)
        return (not ov) and t_w >= g_diag

    if not diag_ok(1.0):
        lo, hi = 1.0, 1.5
        while not diag_ok(hi):
            lo, hi = hi, 2.0*hi
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if diag_ok(mid):
                hi = mid
            else:
                lo = mid
        s = hi
    else:
        s = 1.0
    return s*pH0, s*pV0, t_C, t_H


def sche_geometry(D_ext, H, w_ch_C, h_ch_C, w_ch_H, h_ch_H,
                  P_H, P_V, r0_in, inner_ch='hot', min_wall=0.0, exchange_frac=1.0):
    """
    Géométrie simplifiée du SCHE, cohérente avec setup_geom (branche SCHE).

    Retourne un dict avec V_unit (objectif), A_ht (estimation de A_h), N (tours),
    n_ch, t_wall_min, les jeux entre canaux du même fluide et 'viol' :
    somme des violations géométriques RELATIVES (0 = faisable), utilisée comme
    pénalité graduée pour guider l'essaim vers le domaine faisable.
    """
    p, pv = P_H, P_V
    Cw, Ch, Hw, Hh = w_ch_C, h_ch_C, w_ch_H, h_ch_H

    gaps = {'C_gap_H': p - Cw, 'C_gap_V': pv - Ch, 'H_gap_H': p - Hw, 'H_gap_V': pv - Hh}
    viol = 0.0
    for k, g in gaps.items():
        ref = p if k.endswith('H') else pv
        viol += max(min_wall - g, 0.0) / ref

    t_wall_min, overlap = _diag_wall(Cw, Ch, Hw, Hh, p, pv)
    if overlap:
        viol += 1.0
    viol += max(min_wall - t_wall_min, 0.0) / min(p, pv)

    # Nombre de tours (même construction que setup_geom)
    b_inner, b_outer = (Hw, Cw) if inner_ch == 'hot' else (Cw, Hw)
    c_in = 0.5 * b_inner
    g = max(0.5*p - 0.5*(Cw + Hw), 0.0)
    outer_offset = max(b_inner, c_in + 0.25*p + 0.5*g)
    N = (D_ext/2 - r0_in - outer_offset) / p

    n_ch = int(H // pv)

    A_c, Per_c, Dh_c = _ellipse(Cw, Ch)
    A_h, Per_h, Dh_h = _ellipse(Hw, Hh)
    L_ch = np.pi * ((D_ext/2)**2 - r0_in**2) / p              # longueur d'une spirale (approx.)
    A_ht = L_ch * max(n_ch, 0) * exchange_frac * Per_h         # surface côté chaud (approx.)

    return {
        **gaps, 't_wall_min': t_wall_min, 'overlap': overlap, 'viol': viol,
        'n_turns': N, 'n_p': n_ch, 'L_ch': L_ch, 'A_ht': A_ht,
        'Dh_C': Dh_c, 'Dh_H': Dh_h, 'A_flow_C': A_c * n_ch, 'A_flow_H': A_h * n_ch,
        'V_unit': np.pi * D_ext**2 / 4 * H,
    }


def plot_channel_section(HX, n_cols=4, max_layers=12, fName='', dpi=300):
    """
    Coupe (radiale x verticale) du SCHE : canaux elliptiques froids (bleu) et chauds (rouge) en quinconce.
    Froid en (i p, j pv), chaud décalé de (p/2, pv/2) (même construction que _diag_wall).
    HX      : modèle résolu (paramètres w_ch_C/H, h_ch_C/H, P_H, P_V, H, n_p lus dans HX.params)
    n_cols  : nombre de colonnes radiales dessinées (la coupe se répète tous les P_H)
    max_layers : nombre maximal de couches dessinées (n_ch peut être grand ; la hauteur H complète est tracée)
    Les cotes sont en mm ; la paroi minimale chaud-froid est indiquée.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse, Rectangle

    hp = HX.params
    Cw, Ch, Hw, Hh = hp['w_ch_C'], hp['h_ch_C'], hp['w_ch_H'], hp['h_ch_H']
    p, pv, H = hp['P_H'], hp['P_V'], hp['H']
    n_ch = int(hp.get('n_p', H // pv))
    n_l = min(n_ch, max_layers)
    mm = 1e3
    z0 = 0.5 * (H - n_ch * pv)                       # empilement centré en hauteur
    t_w, _ = _diag_wall(Cw, Ch, Hw, Hh, p, pv)

    fig, ax = plt.subplots(figsize=(5, 4.2))
    ax.add_patch(Rectangle((0, 0), n_cols * p * mm, H * mm, fill=False, ec='0.4', lw=1.0))
    for i in range(n_cols):
        for j in range(n_l):
            xc, zc = (i * p + 0.25 * p) * mm, (z0 + j * pv + 0.25 * pv) * mm
            xh, zh = xc + 0.5 * p * mm, zc + 0.5 * pv * mm
            ax.add_patch(Ellipse((xc, zc), Cw * mm, Ch * mm, fc='tab:blue', ec='k', lw=0.6, alpha=0.8))
            ax.add_patch(Ellipse((xh, zh), Hw * mm, Hh * mm, fc='tab:red', ec='k', lw=0.6, alpha=0.8))
    ax.set_xlim(0, n_cols * p * mm)
    ax.set_ylim(0, H * mm if n_l == n_ch else (z0 + n_l * pv + 0.5 * pv) * mm)
    ax.set_aspect('equal')
    ax.set_xlabel('r [mm]')
    ax.set_ylabel('z [mm]')
    ax.set_title(f"froid {Cw*mm:.1f}x{Ch*mm:.1f} mm | chaud {Hw*mm:.1f}x{Hh*mm:.1f} mm\n"
                 f"pas {p*mm:.1f} / {pv*mm:.1f} mm | paroi min {t_w*mm:.2f} mm" + (
                     f" | {n_l}/{n_ch} couches" if n_l < n_ch else f" | {n_ch} couches"), fontsize=8)
    ax.plot([], [], 's', color='tab:blue', label='froid')
    ax.plot([], [], 's', color='tab:red', label='chaud')
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False)
    plt.tight_layout()
    if fName:
        plt.savefig(fName, dpi=dpi)
    plt.show()


#%%

class SCHESizingOpt(BaseComponent):

    ORDER = ['D_ext', 'H', 'w_ch_C', 'h_ch_C', 'w_ch_H', 'h_ch_H',
             'P_H', 'P_V', 'r0_in', 'n_parallel', 'n_series']

    _GEOM_KEYS = ['D_ext', 'H', 'w_ch_C', 'h_ch_C', 'w_ch_H', 'h_ch_H',
                  'P_H', 'P_V', 'r0_in']

    _CORR_KEYS = {
        'H_Corr', 'C_Corr', 'H_DP', 'C_DP',
        'htc_type', 'DP_type',
        'UD_H_HTC', 'UD_C_HTC', 'UD_H_DP', 'UD_C_DP',
    }

    _CONSTRAINT_KEYS = {'Q_dot', 'DP_h', 'DP_c', 'SH_min'}

    # Clés propres à l'optimiseur, jamais transmises au modèle HX
    _OPT_ONLY_KEYS = {'min_turns', 'min_wall', 'PF_Q_dot', 'PF_DP', 'PF_SH', 'obj', 'w_tol',
                      'C_aspect', 'H_aspect', 'pitch_mode', 'P_design_factor', 'T_design',
                      'accept_pinched'}

    _OBJECTIVES = {
        'V': ('V_tot', 'L',  1e3),
        'A': ('A_tot', 'm2', 1.0),
    }

    DEFAULT_BOUNDS = {
        'D_ext':      [0.05, 0.40],
        'H':          [0.05, 0.40],
        'w_ch_C':  [1e-3, 6e-3],
        'h_ch_C':  [1e-3, 6e-3],
        'w_ch_H':  [1e-3, 6e-3],
        'h_ch_H':  [1e-3, 6e-3],
        'P_H':    [2e-3, 1.2e-2],
        'P_V':    [2e-3, 1.2e-2],
        'r0_in':      [5e-3, 5e-2],
        'n_parallel': [1, 1],     # figé par défaut (lb == ub)
        'n_series':   [1, 1],     # figé par défaut (lb == ub)
    }

    DEFAULT_PARAMETERS = {
        'inner_ch': 'hot',
        'n_disc':        50,
        'Flow_Type':     'CounterFlow',
        'k_wall':        16.0,           # conductivité de paroi [W/m.K] (inox) ; None -> paroi négligée
        'exchange_frac': 1.0,            # fraction du périmètre mouillé qui échange
        'H_Corr': {"1P": "Sigalotti", "2P": "Zakeralhoseini"},
        'C_Corr': {"1P": "Sigalotti", "2P": "Zhao"},
        'H_DP':   {"1P": "Sigalotti", "2P": None},
        'C_DP':   {"1P": "Sigalotti", "2P": None},

        # --- propres à l'optimiseur ---
        'obj':       'V',     # 'V' (volume total) ou 'A' (surface d'échange totale)
        'min_turns': 1.0,     # nombre minimal de tours de spirale
        'min_wall':  3e-4,    # paroi minimale entre deux canaux (fabrication / tenue en pression) [m]
        'C_aspect':  None,    # h/w imposé pour le canal froid (1.0 = circulaire), None = libre
        'H_aspect':  None,    # idem canal chaud
        'w_tol':     1e-3,    # écart max |sum(w) - 1| pour accepter un design convergé
        # Un échangeur surdimensionné renvoie le Q du pincement (sum(w) < 1, HX.pinched = True) au lieu d'échouer :
        # accepté si True (surdimensionnement pénalisé par l'objectif, Q insuffisant par la contrainte Q_dot).
        'accept_pinched': True,
        'pitch_mode':      'pressure',  # 'pressure' : pas calculés pour tenir la pression ; 'free' : optimisés
        'P_design_factor': 1.0,         # P de calcul = facteur * P_su (ex. 1.1, ou 1.43 pour l'épreuve PED)
        'T_design':        None,        # T de calcul [K] ; None -> max des températures d'entrée
        'PF_Q_dot':  1000,
        'PF_DP':     100.0,
        'PF_SH':     1000,    # pénalité par unité de violation relative de SH_min
    }

    def __init__(self):
        super().__init__()

        self.HX = HexMB('SCHE')
        self.bounds = {k: list(v) for k, v in self.DEFAULT_BOUNDS.items()}
        self.params = {**self.DEFAULT_PARAMETERS, **getattr(self, 'params', {})}

        self.corr_params = {
            'H_Corr': None, 'C_Corr': None, 'H_DP': None, 'C_DP': None,
            'htc_type': "Correlation", 'DP_type': "Correlation",
            'UD_H_HTC': None, 'UD_C_HTC': None, 'UD_H_DP': None, 'UD_C_DP': None,
        }
        self.Q_dot_constr = self.DP_h_constr = self.DP_c_constr = None
        self.SH_min_constr = None

    #%% ---------------------------------------------------------------------
    # Corrélations / contraintes / bornes

    def set_corr(self, **kwargs):
        """Update incrémental des corrélations (les clés non passées sont conservées)."""
        unknown = set(kwargs) - self._CORR_KEYS
        if unknown:
            raise KeyError(f"Clés de corrélation inconnues : {unknown}")
        self.corr_params.update(kwargs)
        self._apply_corr(self.HX)

    def _apply_corr(self, HX):
        """Ne mute que le HX reçu en argument -- thread-safe."""
        c = self.corr_params

        if c['htc_type'] == "Correlation":
            HX.set_htc(htc_type='Correlation_Disc', Corr_H=c['H_Corr'], Corr_C=c['C_Corr'])
        else:
            HX.set_htc(htc_type='User-Defined', UD_H_HTC=c['UD_H_HTC'], UD_C_HTC=c['UD_C_HTC'])

        if c['DP_type'] == "Correlation":
            HX.set_DP(DP_type='Correlation_Disc', Corr_H=c['H_DP'], Corr_C=c['C_DP'])
        else:
            HX.set_DP(DP_type='User-Defined', UD_H_DP=c['UD_H_DP'], UD_C_DP=c['UD_C_DP'])

    def set_constraints(self, Q_dot=None, DP_h=None, DP_c=None, SH_min=None):
        """
        Q_dot  : puissance minimale [W] (None = pas de contrainte)
        DP_h/c : pertes de charge maximales [Pa]
        SH_min : surchauffe minimale en sortie du côté froid [K]
        """
        self.Q_dot_constr = Q_dot
        self.DP_h_constr = DP_h
        self.DP_c_constr = DP_c
        self.SH_min_constr = SH_min

    def set_bounds(self, **bounds):
        for key, value in bounds.items():
            if key not in self.ORDER:
                raise KeyError(f"Variable de design inconnue : {key}")
            self.bounds[key] = list(value)

    def _apply_deferred_parameters(self):
        """Extrait corrélations et contraintes de self.params et les dispatche."""
        corr_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CORR_KEYS}
        constraint_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CONSTRAINT_KEYS}

        if corr_kwargs:
            self.set_corr(**corr_kwargs)
        if constraint_kwargs:
            self.set_constraints(**constraint_kwargs)

    def _hx_params(self, params):
        """Sous-ensemble de params à transmettre au modèle HX."""
        return {k: v for k, v in params.items() if k not in self._OPT_ONLY_KEYS}

    #%% ---------------------------------------------------------------------
    # Variables libres / figées / liées

    def _tied(self):
        """Hauteurs liées à la largeur par un rapport d'aspect imposé."""
        tied = {}
        if self.params.get('C_aspect') is not None:
            tied['h_ch_C'] = ('w_ch_C', float(self.params['C_aspect']))
        if self.params.get('H_aspect') is not None:
            tied['h_ch_H'] = ('w_ch_H', float(self.params['H_aspect']))
        return tied

    def _pressure_mode(self):
        mode = str(self.params.get('pitch_mode', 'pressure')).lower()
        if mode not in ('pressure', 'free'):
            raise ValueError(f"pitch_mode = {mode!r} : attendu 'pressure' ou 'free'.")
        return mode == 'pressure'

    def _design_conditions(self):
        """(P_h, P_c, T) de calcul mécanique, à partir des entrées."""
        f = float(self.params.get('P_design_factor', 1.0))
        P_h, P_c = f * self.inputs['P_su_H'], f * self.inputs['P_su_C']

        T = self.params.get('T_design')
        if T is None:
            Ts = []
            for s in ('H', 'C'):
                if self.inputs.get(f'T_su_{s}') is not None:
                    Ts.append(self.inputs[f'T_su_{s}'])
                elif self.inputs.get(f'h_su_{s}') is not None:
                    import CoolProp.CoolProp as CP
                    Ts.append(CP.PropsSI('T', 'H', self.inputs[f'h_su_{s}'],
                                         'P', self.inputs[f'P_su_{s}'], self.inputs[f'fluid_{s}']))
            if not Ts:
                raise ValueError("T_design non fourni et aucune température / enthalpie d'entrée.")
            T = max(Ts)
        return P_h, P_c, float(T)

    def _prepare_variables(self):
        tied = self._tied()
        self._dc = self._design_conditions() if self._pressure_mode() else None
        skip = set(tied) | ({'P_H', 'P_V'} if self._dc is not None else set())
        keys = [k for k in self.ORDER if k not in skip]

        lb = np.array([self.bounds[k][0] for k in keys], dtype=float)
        ub = np.array([self.bounds[k][1] for k in keys], dtype=float)
        if np.any(lb > ub):
            bad = [k for k, l, u in zip(keys, lb, ub) if l > u]
            raise ValueError(f"Borne inférieure > borne supérieure pour : {bad}")

        if 'n_series' in keys and ub[keys.index('n_series')] > 1:
            raise ValueError("SCHE : n_series > 1 n'est pas implémenté dans le modèle (utiliser n_parallel) : "
                             "figer n_series = [1, 1].")

        free = lb < ub
        self._free_keys = [k for k, f in zip(keys, free) if f]
        self._fixed_vals = {k: l for k, l, f in zip(keys, lb, free) if not f}
        self._tied_keys = tied

        if not self._free_keys:
            raise ValueError("Toutes les variables sont figées : rien à optimiser.")

        return lb[free], ub[free]

    def _x_to_design(self, x):
        d = dict(self._fixed_vals)
        d.update(zip(self._free_keys, x))
        for k, (ref, aspect) in self._tied_keys.items():
            d[k] = aspect * d[ref]
        if self._dc is not None:                       # pas imposés par la pression
            P_h, P_c, T = self._dc
            d['P_H'], d['P_V'], d['t_C_req'], d['t_H_req'] = pressure_pitches(
                d['w_ch_C'], d['h_ch_C'], d['w_ch_H'], d['h_ch_H'],
                P_h, P_c, T, self.params['min_wall'])
        d['n_parallel'] = max(1, int(np.round(d['n_parallel'])))
        d['n_series'] = max(1, int(np.round(d['n_series'])))
        return d

    def _geometry(self, d):
        return sche_geometry(**{k: d[k] for k in self._GEOM_KEYS},
                             inner_ch=self.params['inner_ch'],
                             min_wall=self.params['min_wall'],
                             exchange_frac=self.params['exchange_frac'])

    #%% ---------------------------------------------------------------------
    # Score

    @staticmethod
    def _heat_transfer_area(HX, geom):
        """Surface côté chaud d'UNE unité : HX.A_h / n_parallel (A_h contient déjà n_parallel) si valide,
        sinon estimation. Retourne (A, source)."""
        A = getattr(HX, 'A_h', None)
        try:
            A = float(A) / max(1, int(getattr(HX, 'params', {}).get('n_parallel', 1)))
            if np.isfinite(A) and A > 0:
                return A, 'model'
        except (TypeError, ValueError):
            pass
        return geom['A_ht'], 'estimate'

    def _w_deviation(self, HX):
        """|sum(w) - 1| ; np.inf si w absent ou non fini (design rejeté).
        Un design surdimensionné (HX.pinched, sum(w) < 1) est accepté (écart 0) si params['accept_pinched']."""
        try:
            w_sum = getattr(HX, 'w_sum', None)
            if w_sum is None:
                w_sum = np.sum(np.asarray(HX.w, float))
            w_sum = float(w_sum)
            if not np.isfinite(w_sum):
                return np.inf
            if self.params.get('accept_pinched', True) and getattr(HX, 'pinched', False) and w_sum < 1.0:
                return 0.0
            return abs(w_sum - 1.0)
        except Exception:
            return np.inf

    def compute_score(self, HX, geom, n_units, Q_dot_constr, DP_h_constr, DP_c_constr):
        """
        PURE : ne lit que self.params (figé pendant l'optim), n'écrit rien.
        score = objectif + PF_Q * viol_Q + PF_DP * (viol_DP_h + viol_DP_c)
        """
        PF_Q_dot = self.params['PF_Q_dot']
        PF_DP = self.params['PF_DP']
        PF_SH = self.params['PF_SH']
        obj = self.params['obj']
        SH_min = self.SH_min_constr

        V_tot = n_units * geom['V_unit']
        A_unit, A_source = self._heat_transfer_area(HX, geom)
        A_tot = n_units * A_unit
        SH = None

        try:
            pen_Q = max(Q_dot_constr - HX.Q.Q_dot, 0) / Q_dot_constr if Q_dot_constr else 0
            pen_DP_h = max(HX.DP_h - DP_h_constr, 0) / DP_h_constr if DP_h_constr else 0
            pen_DP_c = max(HX.DP_c - DP_c_constr, 0) / DP_c_constr if DP_c_constr else 0

            pen_SH = 0
            if SH_min is not None:
                SH = outlet_superheat(HX, self.inputs['fluid_C'])
                pen_SH = max(SH_min - SH, 0) / max(SH_min, 1.0)   # relative (en K si SH_min < 1 K)

            penalty = PF_Q_dot * pen_Q + PF_DP * (pen_DP_h + pen_DP_c) + PF_SH * pen_SH
            objective = V_tot if obj == 'V' else A_tot
            score = objective + penalty

            if not np.isfinite(score):
                raise ValueError("score non fini")

        except Exception:
            return {
                'score': 1e12, 'penalty': 1e12,
                'pen_Q': None, 'pen_DP_h': None, 'pen_DP_c': None, 'pen_SH': None, 'SH': SH,
                'V_tot': V_tot, 'A_tot': A_tot, 'A_source': A_source,
            }

        return {
            'score': score, 'penalty': penalty,
            'pen_Q': pen_Q, 'pen_DP_h': pen_DP_h, 'pen_DP_c': pen_DP_c, 'pen_SH': pen_SH, 'SH': SH,
            'V_tot': V_tot, 'A_tot': A_tot, 'A_source': A_source,
        }

    #%% ---------------------------------------------------------------------
    # Évaluation d'un candidat

    def _build_local_params(self, d):
        return {**self.params, **{k: d[k] for k in self._GEOM_KEYS},
                'n_parallel': d['n_parallel'], 'n_series': d['n_series']}

    def _new_HX(self, params):
        HX = HexMB('SCHE')
        HX.set_inputs(**self.inputs)
        self._apply_corr(HX)
        HX.set_parameters(**self._hx_params(params))
        return HX

    def simulate_HX(self, x):
        """
        THREAD-SAFE : HX local, score local, n'écrit jamais sur self.
        x : vecteur des variables LIBRES (ordre self._free_keys).
        """
        warnings.filterwarnings('ignore')

        d = self._x_to_design(x)
        if not (np.isfinite(d['P_H']) and np.isfinite(d['P_V'])):
            return 1e9                               # pression impossible à tenir pour ces canaux
        geom = self._geometry(d)

        # ---- faisabilité géométrique (avant tout calcul HX) ----
        viol = geom['viol']
        viol += max(self.params['min_turns'] - geom['n_turns'], 0.0)
        viol += 0.0 if geom['n_p'] >= 1 else 1.0
        if viol > 0:
            return 1e6 * (1 + viol)                 # pénalité graduée

        try:
            HX = self._new_HX(self._build_local_params(d))
            HX.solve()
        except Exception:
            return 1e12

        # ---- convergence : sum(w) doit valoir 1 à w_tol près ----
        w_dev = self._w_deviation(HX)
        if not w_dev <= self.params['w_tol']:
            return 1e6 * (1 + min(w_dev, 1e3))

        n_units = d['n_parallel'] * d['n_series']
        result = self.compute_score(HX, geom, n_units,
                                    self.Q_dot_constr, self.DP_h_constr, self.DP_c_constr)
        return result['score']

    #%% ---------------------------------------------------------------------
    # Optimisation

    def sizing(self, n_jobs=-1, backend="threading", chunksize="auto",
               n_particles=30, max_iter=100, patience=None, rtol=1e-3, obj=None):
        """
        obj  : 'V' (volume total) ou 'A' (surface d'échange totale).
               None -> valeur de self.params['obj'] (défaut 'V').
        rtol : amélioration RELATIVE minimale du meilleur score pour remettre
               le compteur de patience à zéro.
        """
        self._apply_deferred_parameters()

        if obj is not None:
            self.params['obj'] = obj
        obj = str(self.params['obj']).upper()
        if obj not in self._OBJECTIVES:
            raise ValueError(f"obj = {obj!r} : attendu 'V' (volume) ou 'A' (surface).")
        self.params['obj'] = obj
        _, unit, scale = self._OBJECTIVES[obj]

        lb, ub = self._prepare_variables()
        D = lb.size

        optimizer = ps.single.GlobalBestPSO(
            n_particles=n_particles, dimensions=D,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=(lb, ub),
        )

        if patience is None:
            patience = max_iter

        pbar = tqdm(total=n_particles * int(max_iter), unit="eval", dynamic_ncols=True,
                    desc=f'SCHE (min {obj})')

        def objective_wrapper(X):
            Xp = [np.clip(xi, lb, ub) for xi in np.asarray(X)]
            costs = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
                delayed(self.simulate_HX)(xi) for xi in Xp
            )
            pbar.update(len(Xp))
            return np.asarray(costs, dtype=float)

        no_improve, best_cost = 0, np.inf

        for i in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            curr = optimizer.swarm.best_cost

            pbar.set_description(f"Iter {i+1} | Best: {curr*scale:.3f} {unit}")
            pbar.refresh()

            if curr < best_cost * (1 - rtol):
                best_cost = curr
                no_improve = 0
            else:
                no_improve += 1

            if no_improve >= patience:
                break

        pbar.close()

        best = self._x_to_design(optimizer.swarm.best_pos)

        # ---- reconstruction finale mono-thread ----
        self.params.update(self._build_local_params(best))
        self.geom = self._geometry(best)

        if self.geom['viol'] > 0 or self.geom['n_turns'] < self.params['min_turns'] or self.geom['n_p'] < 1:
            raise RuntimeError(
                f"SCHE sizing: aucune géométrie faisable trouvée (best_cost={best_cost:.3e}). "
                f"Élargis les bornes de P_H / P_V ou réduis min_wall."
            )

        self.HX = self._new_HX(self.params)
        try:
            self.HX.solve()
        except Exception as e:
            raise RuntimeError(
                f"SCHE sizing did not converge to a feasible geometry "
                f"(best_cost={best_cost:.3e}); final solve failed with: {e}"
            ) from e

        self.w_dev = self._w_deviation(self.HX)
        if not self.w_dev <= self.params['w_tol']:
            raise RuntimeError(
                f"SCHE sizing: no design converged with |sum(w) - 1| <= {self.params['w_tol']:.1e} "
                f"(best design: |sum(w) - 1| = {self.w_dev:.3e}). Widen the bounds, increase "
                f"n_particles / max_iter, or relax w_tol."
            )

        n_units = best['n_parallel'] * best['n_series']
        result = self.compute_score(self.HX, self.geom, n_units,
                                    self.Q_dot_constr, self.DP_h_constr, self.DP_c_constr)

        self.score = result['score']
        self.penalty = result['penalty']
        self.pen_Q = result['pen_Q']
        self.pen_DP_h = result['pen_DP_h']
        self.pen_DP_c = result['pen_DP_c']
        self.pen_SH = result['pen_SH']
        self.SH = result['SH']
        self.V_tot = result['V_tot']
        self.A_tot = result['A_tot']
        self.A_source = result['A_source']

        self.best_design = best
        return best

    #%% ---------------------------------------------------------------------
    # Sorties

    def print_results(self):
        p, g = self.params, self.geom
        hp = getattr(self.HX, 'params', {})          # géométrie exacte calculée par setup_geom
        mm = 1e3

        print(f"\nObjective : min {'V (volume total)' if p['obj'] == 'V' else 'A (surface d échange totale)'}")
        print("\nBest Position")
        print("-------------")
        print(f"D_ext      : {p['D_ext']:.3f} [m]")
        print(f"H          : {p['H']:.3f} [m]")
        print(f"Canal froid: {p['w_ch_C']*mm:.2f} x {p['h_ch_C']*mm:.2f} [mm]")
        print(f"Canal chaud: {p['w_ch_H']*mm:.2f} x {p['h_ch_H']*mm:.2f} [mm]")
        src = " (calculé : pression)" if self._dc is not None else ""
        print(f"P_H        : {p['P_H']*mm:.2f} [mm]{src}")
        print(f"P_V        : {p['P_V']*mm:.2f} [mm]{src}")
        if self._dc is not None:
            b = self.best_design
            P_h, P_c, T = self._dc
            print(f"Pression   : P_h = {P_h/1e5:.1f} bar, P_c = {P_c/1e5:.1f} bar, T = {T-273.15:.1f} °C")
            print(f"t requis   : froid {b['t_C_req']*mm:.3f} [mm], chaud {b['t_H_req']*mm:.3f} [mm]")
        print(f"r0_in      : {p['r0_in']*mm:.1f} [mm]")
        print(f"n_parallel : {p['n_parallel']} [-]")
        print(f"n_series   : {p['n_series']} [-]")

        print("\nGéométrie")
        print("-------------")
        print(f"N (tours)  : {hp.get('n_turns', g['n_turns']):.2f} [-]")
        print(f"n_ch       : {hp.get('n_p', g['n_p'])} canaux / fluide, "
              f"{hp.get('n_plates', 2*g['n_p'] + 1)} plaques de {0.5*p['P_V']*mm:.2f} mm")
        print(f"Dh froid   : {hp.get('Dh_C', g['Dh_C'])*mm:.2f} [mm]")
        print(f"Dh chaud   : {hp.get('Dh_H', g['Dh_H'])*mm:.2f} [mm]")
        if 'L_H' in hp:
            print(f"L chaud    : {hp['L_H']:.2f} [m]   L froid : {hp['L_C']:.2f} [m]")
        else:
            print(f"L canal    : {g['L_ch']:.2f} [m]   (estimation)")
        print(f"paroi min  : {g['t_wall_min']*mm:.2f} [mm] (min {p['min_wall']*mm:.2f} mm)")
        print(f"A_tot      : {self.A_tot:.3f} [m^2] (source : {self.A_source})")
        print(f"V_unit     : {g['V_unit']*1e3:.2f} [L]")
        print(f"V_tot      : {self.V_tot*1e3:.2f} [L]")

        print("\nResults")
        print("-------------")
        print(f"Q_dot      : {self.HX.Q.Q_dot:,.1f} W   (min {self.Q_dot_constr}; pincé : {getattr(self.HX, 'pinched', None)})")
        if self.SH_min_constr is not None:
            print(f"SH sortie  : {self.SH:.2f} K  (min {self.SH_min_constr} K, côté froid)")
        print(f"DP_h       : {self.HX.DP_h:,.1f} Pa  (max {self.DP_h_constr})")
        print(f"DP_c       : {self.HX.DP_c:,.1f} Pa  (max {self.DP_c_constr})")
        print(f"|sum(w)-1| : {self.w_dev:.2e}  (max {self.params['w_tol']:.1e})")
        print(f"penalty    : {self.penalty:.4g}  "
              f"(viol. rel. Q: {self.pen_Q}, DP_h: {self.pen_DP_h}, DP_c: {self.pen_DP_c}, SH: {self.pen_SH})")

    def export_params_dict(self):
        out = {
            "type": "SCHE",
            "Flow_Type": self.params["Flow_Type"],
            "inner_ch": self.params["inner_ch"],

            "Q_dot": self.HX.Q.Q_dot,
            "DP_h": self.HX.DP_h,
            "DP_c": self.HX.DP_c,
            "SH_out": self.SH,
            "SH_min": self.SH_min_constr,

            "fluid_H": self.inputs['fluid_H'],
            "m_dot_H": self.inputs['m_dot_H'],
            "T_su_H": self.inputs.get('T_su_H'),
            "P_su_H": self.inputs['P_su_H'],

            "fluid_C": self.inputs['fluid_C'],
            "m_dot_C": self.inputs['m_dot_C'],
            "T_su_C": self.inputs.get('T_su_C'),
            "P_su_C": self.inputs['P_su_C'],
        }
        out.update({k: self.params[k] for k in self._GEOM_KEYS})
        out.update({
            "n_parallel": self.params["n_parallel"],
            "n_series": self.params["n_series"],
            "k_wall": self.params.get("k_wall"),
            "exchange_frac": self.params.get("exchange_frac"),
            "pitch_mode": self.params.get("pitch_mode"),
            "t_C_req": self.best_design.get('t_C_req'),
            "t_H_req": self.best_design.get('t_H_req'),
            "n_p": self.geom['n_p'],
            "t_wall_min": self.geom['t_wall_min'],
            "obj": self.params["obj"],
            "V_unit": self.geom['V_unit'],
            "V_tot": self.V_tot,
            "A_tot": self.A_tot,
            "A_source": self.A_source,
            "w_dev": self.w_dev,
        })
        return out

#%%

import CoolProp.CoolProp as CP

if __name__ == "__main__":      # indispensable sous Windows avec joblib
  
    HX_opt = SCHESizingOpt()
    case_study = "Evaporator"
 
    if case_study == "Condenser":
 
        # Propane : T_sat(21.17 bar) = 60 °C, refoulement compresseur à 80 °C
        # Q disponible jusqu'à 57 °C (3 K de sous-refroidissement) ≈ 6.4 kW
        HX_opt.set_inputs(
            fluid_H = 'Propane',
            T_su_H = 80 + 273.15,    # K
            P_su_H = 21.17*1e5,      # Pa
            m_dot_H = 0.020,         # kg/s
 
            fluid_C = 'Water',
            T_su_C = 50 + 273.15,    # K
            P_su_C = 5*1e5,          # Pa
            m_dot_C = 0.286,         # kg/s
        )
 
        Q_dot_cstr = 6e3      # W
        DP_h_cstr  = 40e3     # Pa
        DP_c_cstr  = 40e3     # Pa
 
        HX_opt.set_bounds(
            D_ext     = [0.15, 0.25],
            H         = [0.02, 0.2],
            w_ch_C = [0.5e-3, 0.5e-2],
            w_ch_H = [0.5e-3, 0.5e-2],
            P_H   = [2e-3, 1e-2],      # ignoré si pitch_mode = 'pressure'
            P_V   = [2e-3, 1e-2],
            r0_in     = [3e-2, 3e-2],      # figé
        )
 
        HX_opt.set_parameters(
            obj = "A",

            inner_ch = 'hot',
            n_disc = 100,
            C_aspect = 1.0,                # canaux circulaires
            H_aspect = 1.0,
            min_wall = 3e-4,
            pitch_mode = 'pressure',       # pas calculés pour tenir la pression
            P_design_factor = 1.43,        # pression d'épreuve (PED)
            k_wall = 50, # W/(m*K)         # carbon steel
 
            H_Corr = {"1P": "Sigalotti", "2P": "Zakeralhoseini"},   # condensation côté chaud
            C_Corr = {"1P": "Sigalotti", "2P": "Zhao"},
            H_DP   = {"1P": "Sigalotti", "2P": None},
            C_DP   = {"1P": "Sigalotti", "2P": None},
 
            Q_dot = Q_dot_cstr,
            DP_h  = DP_h_cstr,
            DP_c  = DP_c_cstr,
 
            charge_compute = False,
        )
 
    elif case_study == "Evaporator":
 
        # Propane : T_sat(4.745 bar) = 0 °C
        # Entrée après détente isenthalpique du liquide à 57 °C / 21.17 bar (titre ≈ 0.42)
        # Q pour 1 K de surchauffe ≈ 4.5 kW -> l'eau se refroidit d'environ 4.3 K
        h_su_C = CP.PropsSI('H', 'T', 57 + 273.15, 'P', 21.17*1e5, 'Propane')
 
        HX_opt.set_inputs(
            fluid_H = 'Water',
            T_su_H = 10 + 273.15,    # K
            P_su_H = 2*1e5,          # Pa
            m_dot_H = 0.25,          # kg/s
 
            fluid_C = 'Propane',
            h_su_C = h_su_C,         # J/kg
            P_su_C = 4.745*1e5,      # Pa
            m_dot_C = 0.020,         # kg/s
        )
 
        SH_min_cstr = 1.0     # K : surchauffe minimale du propane en sortie (côté froid)
        DP_h_cstr  = 40e3     # Pa
        DP_c_cstr  = 40e3     # Pa
 
        HX_opt.set_bounds(
            D_ext     = [0.15, 0.25],
            H         = [0.02, 0.15],
            w_ch_C = [0.5e-3, 0.5e-2],
            w_ch_H = [0.5e-3, 0.5e-2],
            P_H   = [2e-3, 1e-2],      # ignoré si pitch_mode = 'pressure'
            P_V   = [2e-3, 1e-2],
            r0_in     = [3e-2, 3e-2],      # figé
        )
 
        HX_opt.set_parameters(
            obj = "A",

            inner_ch = 'hot',
            n_disc = 50,
            C_aspect = 1.0,                # canaux circulaires
            H_aspect = 1.0,
            min_wall = 3e-4,
            pitch_mode = 'pressure',       # pas calculés pour tenir la pression
            P_design_factor = 1.43,        # pression d'épreuve (PED)
            k_wall = 50, # W/(m*K)         # carbon steel
 
            H_Corr = {"1P": "Sigalotti", "2P": "Zakeralhoseini"},
            C_Corr = {"1P": "Sigalotti", "2P": "Zhao"},             # évaporation côté froid
            H_DP   = {"1P": "Sigalotti", "2P": None},
            C_DP   = {"1P": "Sigalotti", "2P": None},
 
            SH_min  = SH_min_cstr,         # à la place de Q_dot
            DP_h    = DP_h_cstr,
            DP_c    = DP_c_cstr,
 
            charge_compute = False,
        )
 
    elif case_study == "CO2_GasCooler":
 
        # CO2 supercritique (P > P_crit = 73.8 bar) : monophasique, seules les corrélations '1P' servent
        # Refroidi de 120 °C à 35 °C : Q ≈ 12.2 kW, l'eau s'échauffe d'environ 9.7 K, pincement ≈ 20 K
        HX_opt.set_inputs(
            fluid_H = 'CO2',
            T_su_H = 120 + 273.15,   # K
            P_su_H = 100*1e5,        # Pa
            m_dot_H = 0.05,          # kg/s
 
            fluid_C = 'Water',
            T_su_C = 15 + 273.15,    # K
            P_su_C = 2*1e5,          # Pa
            m_dot_C = 0.3,           # kg/s
        )
 
        Q_dot_cstr = 12e3     # W
        DP_h_cstr  = 100e3    # Pa
        DP_c_cstr  = 40e3     # Pa
 
        HX_opt.set_bounds(
            D_ext     = [0.10, 0.30],
            H         = [0.02, 0.2],
            w_ch_C = [0.5e-3, 0.5e-2],
            w_ch_H = [0.5e-3, 0.4e-2],  # CO2 dense à haute pression : petits canaux
            P_H   = [2e-3, 1e-2],      # ignoré si pitch_mode = 'pressure'
            P_V   = [2e-3, 1e-2],
            r0_in     = [3e-2, 3e-2],      # figé
        )
 
        HX_opt.set_parameters(
            obj = "A",
            
            inner_ch = 'hot',
            n_disc = 100,                  # résoudre le pic de cp pseudo-critique (~45 °C à 100 bar)
            C_aspect = 1.0,                # canaux circulaires
            H_aspect = 1.0,
            min_wall = 3e-4,
            pitch_mode = 'pressure',       # parois épaisses côté CO2 (143 bar de calcul)
            P_design_factor = 1.43,        # pression d'épreuve (PED)
            k_wall = 50, # W/(m*K)         # carbon steel
 
            H_Corr = {"1P": "Sigalotti", "2P": "Zakeralhoseini"},   # 2P jamais utilisé
            C_Corr = {"1P": "Sigalotti", "2P": "Zhao"},             # 2P jamais utilisé
            H_DP   = {"1P": "Sigalotti", "2P": None},
            C_DP   = {"1P": "Sigalotti", "2P": None},
 
            Q_dot = Q_dot_cstr,
            DP_h  = DP_h_cstr,
            DP_c  = DP_c_cstr,
 
            charge_compute = False,
        )
 
    else:
        raise ValueError(f"case_study inconnu : {case_study!r}")
 
    # ------------------------------------------------------------------ sizing
    HX_opt.sizing(n_particles=30, max_iter=60, patience=15, n_jobs=-1, backend="threading")
 
    HX_opt.print_results()
    plot_channel_section(HX_opt.HX)
