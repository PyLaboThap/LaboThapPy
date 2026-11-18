"""
Author : Basile Chaudoir

Sizing par optimisation (PSO) d'un échangeur à plaques spiralées (SPHE),
calqué sur PCHESizingOpt. Le modèle thermo-hydraulique est
HexMBChargeSensitive('SPHE') (moving boundary, discrétisé), avec la résolution
géométrique de hex_geometry.py (SPHEGeometry : D_ext -> N tours, longueurs, A_unit)
et les corrélations SPHE_HTC / SPHE_DP.

Noms de corrélations (clés '1P' / '2P') :
    HTC 1P  : 'Picon_Nunez', 'Rajavel' (eau), 'Curved_rect_duct', 'Rect_duct'
    HTC 2P  : chaud (condensation) 'Zakeralhoseini' ; froid (évaporation) 'Zhao', 'Nariai' (None : Zhao)
    DP  1P  : 'Picon_Nunez' (Minton), 'Curved_rect_duct', 'Rect_duct', 'Sigalotti', 'Ito'
    DP  2P  : corrélations de tubes sur D_h ('Friedel', 'MSH', 'Choi', ...; None : auto), ou
              'Cioncolini', 'Ruffell', 'Fsadni'

Deux modes d'objectif (paramètre obj de sizing(), ou params['obj']) :
    obj = 'V' : minimiser le volume total
                V_tot = n_parallel * n_series * pi * D_ext^2 / 4 * H        [m^3]
    obj = 'A' : minimiser la surface d'échange totale
                A_tot = n_parallel * n_series * A_unit                      [m^2]
                A_unit = surface d'UNE spirale calculée par le modèle (params['A_unit']),
                sinon estimation spirale d'Archimède (2 * L_plate * H).
sous contraintes de puissance (Q_dot) et de pertes de charge (DP_h, DP_c).

Variables de design (ORDER) :
    D_ext       diamètre extérieur de la spirale        [m]
    H           hauteur (largeur des plaques)           [m]
    C_canal_t   épaisseur (gap) du canal froid          [m]
    H_canal_t   épaisseur (gap) du canal chaud          [m]
    r0_in       rayon intérieur de départ de la spirale [m]
    n_parallel  nombre d'unités en parallèle            [-] (entier)
    n_series    nombre d'unités en série                [-] (entier)

Une variable dont les bornes sont égales (lb == ub) est FIGÉE : elle est
retirée du vecteur PSO et injectée telle quelle dans chaque évaluation.
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

def spiral_geometry(D_ext, H, C_canal_t, H_canal_t, r0_in, t):
    """
    Géométrie simplifiée d'une spirale d'Archimède à deux plaques.

    Pas radial par tour : p = C_canal_t + H_canal_t + 2t.
    V_unit = pi D_ext^2 / 4 * H : volume enveloppe d'une unité (objectif).
    n_turns sert à la faisabilité ; L_plate et A_ht sont informatifs.
    """
    R_ext = D_ext / 2
    pitch = C_canal_t + H_canal_t + 2 * t
    n_turns = (R_ext - r0_in) / pitch
    L_plate = np.pi * (R_ext**2 - r0_in**2) / pitch
    A_ht = 2 * L_plate * H
    V_unit = np.pi * D_ext**2 / 4 * H

    return {
        'R_ext': R_ext, 'pitch': pitch, 'n_turns': n_turns,
        'L_plate': L_plate, 'A_ht': A_ht, 'V_unit': V_unit,
    }


#%%

class SPHESizingOpt(BaseComponent):

    ORDER = ['D_ext', 'H', 'C_canal_t', 'H_canal_t', 'r0_in', 'n_parallel', 'n_series']

    _CORR_KEYS = {
        'H_Corr', 'C_Corr', 'H_DP', 'C_DP',
        'htc_type', 'DP_type',
        'UD_H_HTC', 'UD_C_HTC', 'UD_H_DP', 'UD_C_DP',
    }

    _CONSTRAINT_KEYS = {'Q_dot', 'DP_h', 'DP_c'}

    # Clés propres à l'optimiseur, jamais transmises au modèle HX
    _OPT_ONLY_KEYS = {'min_turns', 'PF_Q_dot', 'PF_DP', 'obj', 'w_tol', 'accept_pinched'}

    # Modes d'objectif : clé du résultat, unité affichée, facteur d'affichage
    _OBJECTIVES = {
        'V': ('V_tot', 'L',   1e3),
        'A': ('A_tot', 'm2',  1.0),
    }

    # Bornes par défaut, ordre de grandeur d'un SPHE de petite/moyenne taille.
    # set_bounds() ne fait qu'un update incrémental.
    DEFAULT_BOUNDS = {
        'D_ext':      [0.15, 1.2],
        'H':          [0.10, 1.5],
        'C_canal_t':  [1e-3, 1e-2],
        'H_canal_t':  [1e-3, 1e-2],
        'r0_in':      [1e-2, 1e-1],
        'n_parallel': [1, 1],     # figé par défaut (lb == ub)
        'n_series':   [1, 1],     # figé par défaut (lb == ub)
    }

    DEFAULT_PARAMETERS = {
        't':             0.002,          # épaisseur de plaque [m]
        'inner_channel': 'cold',
        'n_disc':        50,
        'Flow_Type':     'CounterFlow',
        'H_Corr': {"1P": "Picon_Nunez", "2P": "Zakeralhoseini"},
        'C_Corr': {"1P": "Picon_Nunez", "2P": "Zhao"},
        'H_DP':   {"1P": "Picon_Nunez", "2P": None},          # 2P None : Friedel / MSH (auto)
        'C_DP':   {"1P": "Picon_Nunez", "2P": None},

        # --- propres à l'optimiseur ---
        'obj':       'V',    # 'V' (volume total) ou 'A' (surface d'échange totale)
        'min_turns': 1.0,    # nombre minimal de tours de spirale
        'w_tol':     1e-3,   # écart max |sum(w) - 1| pour accepter un design convergé
        # Un échangeur surdimensionné renvoie le Q du pincement (sum(w) < 1, HX.pinched = True) au lieu
        # d'échouer : accepté si True (le surdimensionnement est pénalisé par l'objectif, un Q insuffisant
        # par la contrainte Q_dot).
        'accept_pinched': True,
        # Facteurs de pénalité [unité de l'objectif par unité de violation RELATIVE] :
        # m^3 en mode 'V', m^2 en mode 'A'. Ils doivent rester grands devant la
        # valeur de l'objectif pour que les contraintes soient quasi dures
        # (ex. PF = 1000 : 1 % de Q_dot manquant coûte 10 m^3 ou 10 m^2).
        'PF_Q_dot':  1000,
        'PF_DP':     100.0,
    }

    def __init__(self):
        super().__init__()

        self.HX = HexMB('SPHE')
        self.bounds = {k: list(v) for k, v in self.DEFAULT_BOUNDS.items()}
        self.params = {**self.DEFAULT_PARAMETERS, **getattr(self, 'params', {})}

        self.corr_params = {
            'H_Corr': None, 'C_Corr': None, 'H_DP': None, 'C_DP': None,
            'htc_type': "Correlation", 'DP_type': "Correlation",
            'UD_H_HTC': None, 'UD_C_HTC': None, 'UD_H_DP': None, 'UD_C_DP': None,
        }
        self.Q_dot_constr = self.DP_h_constr = self.DP_c_constr = None

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

    def set_constraints(self, Q_dot=None, DP_h=None, DP_c=None):
        self.Q_dot_constr = Q_dot
        self.DP_h_constr = DP_h
        self.DP_c_constr = DP_c

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
    # Variables libres / figées

    def _prepare_variables(self):
        lb = np.array([self.bounds[k][0] for k in self.ORDER], dtype=float)
        ub = np.array([self.bounds[k][1] for k in self.ORDER], dtype=float)
        if np.any(lb > ub):
            bad = [k for k, l, u in zip(self.ORDER, lb, ub) if l > u]
            raise ValueError(f"Borne inférieure > borne supérieure pour : {bad}")

        i_s = self.ORDER.index('n_series')
        if ub[i_s] > 1:
            raise ValueError("SPHE : n_series > 1 n'est pas implémenté dans le modèle (utiliser n_parallel) : "
                             "figer n_series = [1, 1].")

        free = lb < ub
        self._free_keys = [k for k, f in zip(self.ORDER, free) if f]
        self._fixed_vals = {k: l for k, l, f in zip(self.ORDER, lb, free) if not f}

        if not self._free_keys:
            raise ValueError("Toutes les variables sont figées : rien à optimiser.")

        return lb[free], ub[free]

    def _x_to_design(self, x):
        d = dict(self._fixed_vals)
        d.update(zip(self._free_keys, x))
        d['n_parallel'] = max(1, int(np.round(d['n_parallel'])))
        d['n_series'] = max(1, int(np.round(d['n_series'])))
        return d

    #%% ---------------------------------------------------------------------
    # Score

    @staticmethod
    def _heat_transfer_area(HX, geom):
        """
        Surface d'échange d'UNE spirale : params['A_unit'] écrit par la résolution géométrique
        du modèle (hex_geometry.SPHEGeometry) si elle existe et est valide, sinon l'estimation
        spirale d'Archimède. (HX.A_h contient déjà le facteur n_parallel : il n'est pas utilisé ici.)
        Retourne (A_unit [m^2], source).
        """
        A = getattr(HX, 'params', {}).get('A_unit')
        try:
            A = float(A)
            if np.isfinite(A) and A > 0:
                return A, 'model'
        except (TypeError, ValueError):
            pass
        return geom['A_ht'], 'estimate'

    def _w_deviation(self, HX):
        """
        Écart de fermeture du moving boundary : |sum(w) - 1|.
        Un design bien convergé remplit exactement l'échangeur (sum(w) = 1).
        Retourne np.inf si w est absent ou non fini (design rejeté).
        Un design surdimensionné (HX.pinched, sum(w) < 1) est accepté (écart 0) si params['accept_pinched'].
        """
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
        objectif = V_tot [m^3] si obj == 'V', A_tot [m^2] si obj == 'A'.
        Les violations sont RELATIVES (sans dimension), donc indépendantes des
        unités W / Pa. V_tot et A_tot sont toujours calculés tous les deux.
        """
        PF_Q_dot = self.params['PF_Q_dot']
        PF_DP = self.params['PF_DP']
        obj = self.params['obj']

        V_tot = n_units * geom['V_unit']
        A_unit, A_source = self._heat_transfer_area(HX, geom)
        A_tot = n_units * A_unit

        try:
            pen_Q = max(Q_dot_constr - HX.Q.Q_dot, 0) / Q_dot_constr if Q_dot_constr else 0
            pen_DP_h = max(HX.DP_h - DP_h_constr, 0) / DP_h_constr if DP_h_constr else 0
            pen_DP_c = max(HX.DP_c - DP_c_constr, 0) / DP_c_constr if DP_c_constr else 0

            penalty = PF_Q_dot * pen_Q + PF_DP * (pen_DP_h + pen_DP_c)
            objective = V_tot if obj == 'V' else A_tot
            score = objective + penalty

            if not np.isfinite(score):
                raise ValueError("score non fini")

        except Exception:
            return {
                'score': 1e12, 'penalty': 1e12,
                'pen_Q': None, 'pen_DP_h': None, 'pen_DP_c': None,
                'V_tot': V_tot, 'A_tot': A_tot, 'A_source': A_source,
            }

        return {
            'score': score, 'penalty': penalty,
            'pen_Q': pen_Q, 'pen_DP_h': pen_DP_h, 'pen_DP_c': pen_DP_c,
            'V_tot': V_tot, 'A_tot': A_tot, 'A_source': A_source,
        }

    #%% ---------------------------------------------------------------------
    # Évaluation d'un candidat

    def _build_local_params(self, d):
        return {
            **self.params,
            'D_ext': d['D_ext'], 'H': d['H'],
            'C_canal_t': d['C_canal_t'], 'H_canal_t': d['H_canal_t'],
            'r0_in': d['r0_in'],
            'n_parallel': d['n_parallel'], 'n_series': d['n_series'],
        }

    def simulate_HX(self, x):
        """
        THREAD-SAFE : HX local, score local, n'écrit jamais sur self.
        x : vecteur des variables LIBRES (ordre self._free_keys).
        """
        warnings.filterwarnings('ignore')

        d = self._x_to_design(x)
        t = self.params['t']
        geom = spiral_geometry(d['D_ext'], d['H'], d['C_canal_t'], d['H_canal_t'], d['r0_in'], t)

        # ---- faisabilité : au moins min_turns tours de spirale ----
        if geom['n_turns'] < self.params['min_turns']:
            # pénalité graduée pour guider l'essaim vers le domaine faisable
            return 1e6 * (1 + self.params['min_turns'] - geom['n_turns'])

        local_params = self._build_local_params(d)

        HX = HexMB('SPHE')
        HX.set_inputs(**self.inputs)
        self._apply_corr(HX)
        HX.set_parameters(**self._hx_params(local_params))

        try:
            HX.solve()
        except Exception:
            return 1e12

        # ---- convergence : sum(w) doit valoir 1 à w_tol près ----
        w_dev = self._w_deviation(HX)
        if not w_dev <= self.params['w_tol']:
            # rejeté ; pénalité graduée (finie) pour guider l'essaim
            return 1e6 * (1 + min(w_dev, 1e3))

        n_units = d['n_parallel'] * d['n_series']
        result = self.compute_score(
            HX, geom, n_units,
            self.Q_dot_constr, self.DP_h_constr, self.DP_c_constr,
        )
        return result['score']

    #%% ---------------------------------------------------------------------
    # Optimisation

    def sizing(self, n_jobs=-1, backend="threading", chunksize="auto",
               n_particles=30, max_iter=100, patience=None, rtol=1e-3, obj=None):
        """
        obj  : 'V' (volume total) ou 'A' (surface d'échange totale).
               None -> valeur de self.params['obj'] (défaut 'V').
        rtol : amélioration RELATIVE minimale du meilleur score pour remettre
        le compteur de patience à zéro (indépendante de l'unité de l'objectif).
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

        pbar = tqdm(total=n_particles * int(max_iter), unit="eval", dynamic_ncols=True, desc=f'SPHE (min {obj})')

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

        best_free = optimizer.swarm.best_pos
        best = self._x_to_design(best_free)

        # ---- reconstruction finale mono-thread ----
        self.params.update(self._build_local_params(best))
        self.geom = spiral_geometry(best['D_ext'], best['H'], best['C_canal_t'],
                                    best['H_canal_t'], best['r0_in'], self.params['t'])

        self.HX = HexMB('SPHE')
        self.HX.set_inputs(**self.inputs)
        self._apply_corr(self.HX)
        self.HX.set_parameters(**self._hx_params(self.params))

        try:
            self.HX.solve()
        except Exception as e:
            raise RuntimeError(
                f"SPHE sizing did not converge to a feasible geometry "
                f"(best_cost={best_cost:.3e}); final solve failed with: {e}"
            ) from e

        self.w_dev = self._w_deviation(self.HX)
        if not self.w_dev <= self.params['w_tol']:
            raise RuntimeError(
                f"SPHE sizing: no design converged with |sum(w) - 1| <= {self.params['w_tol']:.1e} "
                f"(best design: |sum(w) - 1| = {self.w_dev:.3e}). Widen the bounds, increase "
                f"n_particles / max_iter, or relax w_tol."
            )

        n_units = best['n_parallel'] * best['n_series']

        result = self.compute_score(
            self.HX, self.geom, n_units,
            self.Q_dot_constr, self.DP_h_constr, self.DP_c_constr,
        )

        self.score = result['score']
        self.penalty = result['penalty']
        self.pen_Q = result['pen_Q']
        self.pen_DP_h = result['pen_DP_h']
        self.pen_DP_c = result['pen_DP_c']
        self.V_tot = result['V_tot']
        self.A_tot = result['A_tot']
        self.A_source = result['A_source']

        self.best_design = best
        return best

    #%% ---------------------------------------------------------------------
    # Sorties

    def print_results(self):
        p, g = self.params, self.geom

        print(f"\nObjective : min {'V (volume total)' if p['obj'] == 'V' else 'A (surface d échange totale)'}")
        print("\nBest Position")
        print("-------------")
        print(f"D_ext      : {p['D_ext']:.3f} [m]")
        print(f"H          : {p['H']:.3f} [m]")
        print(f"C_canal_t  : {p['C_canal_t']*1e3:.2f} [mm]")
        print(f"H_canal_t  : {p['H_canal_t']*1e3:.2f} [mm]")
        print(f"r0_in      : {p['r0_in']*1e3:.1f} [mm]")
        print(f"t          : {p['t']*1e3:.2f} [mm] (fixe)")
        print(f"n_parallel : {p['n_parallel']} [-]")
        print(f"n_series   : {p['n_series']} [-]")

        print("\nGéométrie")
        print("-------------")
        hp = self.HX.params
        print(f"n_turns    : {hp.get('N_turns', g['n_turns']):.2f} [-]   (modèle ; estimation {g['n_turns']:.2f})")
        print(f"D_outer    : {hp.get('D_outer', float('nan')):.3f} [m]")
        print(f"L_hot      : {hp.get('L_hot', float('nan')):.2f} [m]   L_cold : {hp.get('L_cold', float('nan')):.2f} [m]")
        print(f"A_unit     : {hp.get('A_unit', g['A_ht']):.3f} [m^2] (modèle ; estimation Archimède {g['A_ht']:.2f})")
        print(f"A_tot      : {self.A_tot:.3f} [m^2] (source : {self.A_source})")
        print(f"V_unit     : {g['V_unit']*1e3:.2f} [L]")
        print(f"V_tot      : {self.V_tot*1e3:.2f} [L]")

        print("\nResults")
        print("-------------")
        print(f"Q_dot      : {self.HX.Q.Q_dot:,.1f} W   (cible {self.Q_dot_constr}; pincé : {getattr(self.HX, 'pinched', None)})")
        print(f"DP_h       : {self.HX.DP_h:,.1f} Pa  (max {self.DP_h_constr})")
        print(f"DP_c       : {self.HX.DP_c:,.1f} Pa  (max {self.DP_c_constr})")
        print(f"|sum(w)-1| : {self.w_dev:.2e}  (max {self.params['w_tol']:.1e})")
        print(f"penalty    : {self.penalty:.4g}  "
              f"(viol. rel. Q: {self.pen_Q}, DP_h: {self.pen_DP_h}, DP_c: {self.pen_DP_c})")

    def export_params_dict(self):
        return {
            "type": "SPHE",
            "Flow_Type": self.params["Flow_Type"],
            "inner_channel": self.params["inner_channel"],

            "Q_dot": self.HX.Q.Q_dot,
            "DP_h": self.HX.DP_h,
            "DP_c": self.HX.DP_c,

            "fluid_H": self.inputs['fluid_H'],
            "m_dot_H": self.inputs['m_dot_H'],
            "T_su_H": self.inputs['T_su_H'],
            "P_su_H": self.inputs['P_su_H'],

            "fluid_C": self.inputs['fluid_C'],
            "m_dot_C": self.inputs['m_dot_C'],
            "T_su_C": self.inputs['T_su_C'],
            "P_su_C": self.inputs['P_su_C'],

            "D_ext": self.params["D_ext"],
            "H": self.params["H"],
            "C_canal_t": self.params["C_canal_t"],
            "H_canal_t": self.params["H_canal_t"],
            "r0_in": self.params["r0_in"],
            "t": self.params["t"],
            "n_parallel": self.params["n_parallel"],
            "n_series": self.params["n_series"],

            "obj": self.params["obj"],
            "V_unit": self.geom['V_unit'],
            "V_tot": self.V_tot,
            "A_tot": self.A_tot,
            "A_source": self.A_source,
            "w_dev": self.w_dev,
        }


#%%

import CoolProp.CoolProp as CP

if __name__ == "__main__":
    import time
    start_time = time.time()
    
    HX_opt = SPHESizingOpt()
    case_study = "Evaporator"
    
    if case_study == "Condenser":

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
        
        # Contraintes à adapter au cahier des charges
        Q_dot_cstr = 6e3      # W
        DP_h_cstr  = 40e3     # Pa
        DP_c_cstr  = 40e3     # Pa
        
        HX_opt.set_bounds(
            D_ext     = [0.2, 0.3],
            H         = [0.05, 0.15],
            C_canal_t = [1e-3, 4e-3],
            H_canal_t = [1e-3, 2e-3],
            r0_in     = [2.49e-2, 2.51e-2],
        )
    
        # Seuls les paramètres différents des DEFAULT_PARAMETERS sont nécessaires
        HX_opt.set_parameters(
            t = 0.002,
            inner_channel = 'cold',
            n_disc = 100,
            
            H_Corr = {"1P": "Picon_Nunez", "2P": "Zakeralhoseini"},
            C_Corr = {"1P": "Picon_Nunez", "2P": None},
            H_DP   = {"1P": "Picon_Nunez", "2P": None},
            C_DP   = {"1P": "Picon_Nunez", "2P": None},
            
            Q_dot = Q_dot_cstr,
            DP_h  = DP_h_cstr,
            DP_c  = DP_c_cstr,
            
            charge_compute = False,
            
            k_wall = 50 # W/(m*K)
        )
        
    elif case_study == "Evaporator":
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
        
        # Contraintes à adapter au cahier des charges
        Q_dot_cstr = 4.4e3      # W
        DP_h_cstr  = 40e3     # Pa
        DP_c_cstr  = 40e3     # Pa
        
        HX_opt.set_bounds(
            D_ext     = [0.15, 0.3],
            H         = [0.05, 0.3],
            C_canal_t = [1e-3, 5e-3],
            H_canal_t = [1e-3, 5e-3],
            r0_in     = [2.49e-2, 2.51e-2],
        )
        
        # Seuls les paramètres différents des DEFAULT_PARAMETERS sont nécessaires
        HX_opt.set_parameters(
            t = 0.002,
            inner_channel = 'hot',
            n_disc = 100,
            
            H_Corr = {"1P": "Picon_Nunez", "2P": "Zakeralhoseini"},
            C_Corr = {"1P": "Picon_Nunez", "2P": "Zhao"},
            H_DP   = {"1P": "Picon_Nunez", "2P": None},
            C_DP   = {"1P": "Picon_Nunez", "2P": None},
            
            Q_dot = Q_dot_cstr,
            DP_h  = DP_h_cstr,
            DP_c  = DP_c_cstr,
            
            charge_compute = False,
            
            k_wall = 50 # W/(m*K)
        )
        
    # obj = 'V' : minimiser le volume total | obj = 'A' : minimiser la surface d'échange
    best = HX_opt.sizing(n_jobs=-1, n_particles=100, patience=10, obj='A')
    HX_opt.print_results()
    
    print(f"\nTemps total : {time.time() - start_time:.1f} s")
    
    from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.plot_SPHE import plot_spiral
    HX_opt.HX.plot_cells()
    plot_spiral(HX_opt.HX)