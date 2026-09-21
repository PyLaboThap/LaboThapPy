"""
Author : Basile Chaudoir
"""

import os
# os.environ.setdefault("OMP_NUM_THREADS", "1")
# os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
# os.environ.setdefault("MKL_NUM_THREADS", "1")
# os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
# os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

from labothappy.component.base_component import BaseComponent
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive

from labothappy.toolbox.heat_exchangers.PCHE.thicknesses import PCHE_thicknesses
from labothappy.toolbox.economics.cpi_data import actualize_price

import pyswarms as ps
import numpy as np
import warnings

#%%

# ---- joblib worker (process-based) ----
import os, numpy as np
from joblib import Parallel, delayed

from contextlib import contextmanager
from tqdm import tqdm
import joblib
from joblib.parallel import BatchCompletionCallBack

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

_SOLVER = None  # cached per-process

#%%

class PCHESizingOpt(BaseComponent):

    _CORR_KEYS = {
        'H_Corr', 'C_Corr', 'H_DP', 'C_DP',
        'htc_type', 'DP_type',
        'UD_H_HTC', 'UD_C_HTC', 'UD_H_DP', 'UD_C_DP',
    }
    _CONSTRAINT_KEYS = {'Q_dot', 'DP_h', 'DP_c'}

    # Bornes par défaut, typiques d'un canal PCHE zigzag CO2 supercritique.
    # set_bounds() ne fait qu'un update incrémental : ces valeurs restent
    # actives pour toute clé non explicitement surchargée.
    DEFAULT_BOUNDS = {
        'alpha':      [10, 40],
        'D_c':        [5e-4, 3e-3],
        'L_x':        [1, 1.5],
        'L_y':        [1.4, 2.3],
        'L_z':        [0.2, 0.6],
        'n_parallel': [1, 20],
        'n_series':   [1, 10],
    }

    # Paramètres par défaut, communs à la plupart des dimensionnements PCHE.
    # set_parameters() (hérité de BaseComponent) fait aussi un update
    # incrémental sur self.params, donc ces valeurs restent actives pour
    # toute clé non explicitement surchargée par l'appelant.
    DEFAULT_PARAMETERS = {
        'k_cond':    60,
        'R_p':       1,
        'n_disc':    30,
        'Flow_Type': 'CounterFlow',
        'H_DP_ON':   True,
        'C_DP_ON':   True,
        'H_Corr' : {"1P" : "Gnielinski", "SC" : "Gnielinski"},
        'C_Corr' : {"1P" : "Gnielinski", "SC" : "Gnielinski"},
        'H_DP'   : {"1P" : "Gnielinski_DP", "SC" : "Gnielinski_DP"},
        'C_DP'   : {"1P" : "Gnielinski_DP", "SC" : "Gnielinski_DP"},
        'obj': 'cost',
    }

    def __init__(self):
        super().__init__()

        self.HX = HexMBChargeSensitive('PCHE')
        self.bounds = dict(self.DEFAULT_BOUNDS)     # copie modifiable des defaults

        # self.params est initialisé par BaseComponent (probablement {}) ;
        # on pose les defaults sans écraser ce qui y serait déjà présent.
        self.params = {**self.DEFAULT_PARAMETERS, **getattr(self, 'params', {})}

        return

    #%%

    def set_corr(self, H_Corr=None, C_Corr=None, H_DP=None, C_DP=None, htc_type="Correlation", DP_type="Correlation", UD_H_HTC=None, UD_C_HTC=None, UD_H_DP=None, UD_C_DP=None):

        # Sauvegarder pour réutilisation dans simulate_HX
        self.corr_params = {
            'H_Corr': H_Corr, 'C_Corr': C_Corr,
            'H_DP': H_DP, 'C_DP': C_DP,
            'htc_type': htc_type, 'DP_type': DP_type,
            'UD_H_HTC': UD_H_HTC, 'UD_C_HTC': UD_C_HTC,
            'UD_H_DP': UD_H_DP, 'UD_C_DP': UD_C_DP,
        }

        self._apply_corr(self.HX)

    def _apply_corr(self, HX):
        """
        Applique les corrélations sauvegardées sur n'importe quel objet HX
        passé en argument. Lit self.corr_params/self.params['H_DP_ON']/
        ['C_DP_ON'] (lecture seule, jamais mutés pendant une optimisation) ;
        ne mute QUE le HX local reçu en argument -- thread-safe tant que
        chaque appelant lui passe son propre HX.
        """
        c = self.corr_params
    
        if c['htc_type'] == "Correlation":
            HX.set_htc(Corr_H=c['H_Corr'], Corr_C=c['C_Corr'])
        else:
            HX.set_htc(UD_H_HTC=c['UD_H_HTC'], UD_C_HTC=c['UD_C_HTC'])
    
        HX.params['H_DP_ON'] = self.params['H_DP_ON']
        HX.params['C_DP_ON'] = self.params['C_DP_ON']
    
        if c['DP_type'] == "Correlation":
            HX.set_DP(DP_type='Correlation_Disc', Corr_H=c['H_DP'], Corr_C=c['C_DP'])
        else:
            HX.set_DP(DP_type='User-Defined', UD_H_DP=c['UD_H_DP'], UD_C_DP=c['UD_C_DP'])


    def set_bounds(self, **bounds):
        for key, value in bounds.items():
            self.bounds[key] = value

    def _apply_deferred_parameters(self):
        """
        set_parameters() reste inchangé (hérité de BaseComponent) et stocke
        tout dans self.params, y compris les clés de corrélation et de
        contraintes. Cette méthode les en extrait et les dispatche vers
        set_corr()/set_constraints() — appelée juste avant l'optimisation.
        """
        corr_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CORR_KEYS}
        constraint_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CONSTRAINT_KEYS}

        if corr_kwargs:
            self.set_corr(**corr_kwargs)
        if constraint_kwargs:
            self.set_constraints(**constraint_kwargs)

    def set_constraints(self, Q_dot = None, DP_h = None, DP_c = None):

        self.Q_dot_constr = Q_dot
        self.DP_h_constr = DP_h
        self.DP_c_constr = DP_c

        return

    #%%
    def export_params_dict(self):

        return {
            "type": "PCHE",
            "type_channel": "Zigzag",
            "Flow_Type": self.params["Flow_Type"],

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

            "alpha": self.params["alpha"],
            "D_c": self.params["D_c"],
            "k_cond": self.params["k_cond"],
            "L_c": self.params["L_c"],
            "N_c": self.params["N_c"],
            "N_p": self.params["N_p"],
            "R_p": self.params["R_p"],
            "t_2": self.params["t_2"],
            "t_3": self.params["t_3"],
            "L_x": self.params["L_x"],
            "L_y": self.params["L_y"],
            "L_z": self.params["L_z"],

            "CAPEX": self.CAPEX,
        }

    #%%

    def compute_geom(self):

        P_max = max(self.inputs['P_su_H'], self.inputs['P_su_C'])*1.5
        T_max = max(self.inputs['T_su_H'], self.inputs['T_su_C'])*1.2

        self.params['t_2'], self.params['t_3'] = t_2, t_3 = PCHE_thicknesses(self.params['D_c'], P_max, T_max)
        self.params['L_c'] = L_c = self.params['L_x']/np.cos(np.pi*self.params['alpha']/180)

        self.params['N_p'] = N_p = np.floor((self.params['L_y']-t_3)/(self.params['D_c']/2 + t_3))
        self.params['N_c'] = N_c = np.floor((self.params['L_z']-t_2)/(self.params['D_c'] + t_2))

        V_channel = np.pi*self.params['D_c']**2 / 8 * L_c

        self.params['C_V_tot'] = 1/(1+self.params['R_p']) * N_p * N_c * V_channel
        self.params['H_V_tot'] = self.params['R_p']/(1+self.params['R_p']) * N_p * N_c * V_channel

        return

    #%%

    def compute_score(self, HX, m_HX, Q_dot_constr, DP_h_constr, DP_c_constr, obj):
        """
        UNIQUE fonction de calcul de score, dans tout le fichier. PURE :
        ne lit ni n'écrit aucun attribut mutable de self pendant le calcul
        (pas de self.HX, self.m_HX, self.score, self.penalty, self.CAPEX,
        etc.) -- tout lui est passé en argument, elle retourne tout dans un
        dict. Thread-safe par construction : peut être appelée
        concurremment par plusieurs threads tant que chacun lui passe son
        propre HX/m_HX locaux.

        L'assignation sur self (quand elle est utile, typiquement pour le
        rebuild final mono-thread) est laissée À L'APPELANT -- voir
        `sizing()`, seul endroit du fichier qui fait
        `self.score = result['score']` etc.

        Retourne un dict :
            {'score', 'penalty', 'pen_Q', 'pen_DP_h', 'pen_DP_c',
             'm_HX', 'UA', 'CAPEX_alt', 'CAPEX'}
        où 'CAPEX' est lui-même un dict {'HX','Currency','Install','Total'}
        ou None si obj != 'cost'.
        """
        PF_Q_dot = 10
        PF_DP = 1

        try:
            # Penalties
            pen_Q = max(Q_dot_constr - HX.Q.Q_dot, 0) if Q_dot_constr else 0
            pen_DP_h = max(HX.DP_h - DP_h_constr, 0) if DP_h_constr else 0
            pen_DP_c = max(HX.DP_c - DP_c_constr, 0) if DP_c_constr else 0

            penalty = PF_DP * abs(pen_DP_c) + PF_DP * abs(pen_DP_h) + PF_Q_dot * abs(pen_Q)

            # --- Objectif : masse ou coût (Weiland et al.) ---
            UA = CAPEX_alt = None
            capex = None
            if obj == 'cost':
                alpha_h_mean = sum(HX.alpha_h * HX.w) / HX.w_sum
                alpha_c_mean = sum(HX.alpha_c * HX.w) / HX.w_sum

                UA = HX.Q.Q_dot / np.sum(HX.LMTD * HX.w)

                # sCO2 Power Cycle Component Cost Correlations From DOE Data
                # Spanning Multiple Scales and Applications (2019), Eq.(10)
                CAPEX_alt = 49.45 * UA ** 0.7544

                HX_cost = actualize_price(CAPEX_alt, 2017, "USD")
                install = HX_cost * 0.35
                capex = {
                    "HX": HX_cost,
                    "Currency": "USD",
                    "Install": install,
                    "Total": HX_cost + install,
                }
                score = capex["Total"] + penalty
            else:
                score = m_HX + penalty

        except Exception:
            return {
                'score': 1e12, 'penalty': 1e12,
                'pen_Q': None, 'pen_DP_h': None, 'pen_DP_c': None,
                'm_HX': m_HX, 'UA': None, 'CAPEX_alt': None, 'CAPEX': None,
            }

        return {
            'score': score, 'penalty': penalty,
            'pen_Q': pen_Q, 'pen_DP_h': pen_DP_h, 'pen_DP_c': pen_DP_c,
            'm_HX': m_HX, 'UA': UA, 'CAPEX_alt': CAPEX_alt, 'CAPEX': capex,
        }

    #%%

    def simulate_HX(self, x):
        """
        THREAD-SAFE : construit un HX local, calcule un score local via
        compute_score(), et ne lit QUE des attributs self en lecture seule
        (self.inputs, self.params, self.corr_params, self.Q_dot_constr,
        self.DP_h_constr, self.DP_c_constr) -- aucun ne doit être muté par
        un autre thread pendant l'optimisation. N'écrit JAMAIS sur self
        (plus de self.HX = HX, self.set_parameters(**local_params), etc.).
        """
        warnings.filterwarnings('ignore')

        alpha      = x[0]
        D_c        = x[1]
        L_x        = x[2]
        L_y        = x[3]
        L_z        = x[4]
        n_parallel = np.round(x[5])
        n_series   = max(1, np.round(x[6]))

        # ---- local geometry ----
        P_max = max(self.inputs['P_su_H'], self.inputs['P_su_C']) * 1.5
        T_max = max(self.inputs['T_su_H'], self.inputs['T_su_C']) * 1.2

        t_2, t_3 = PCHE_thicknesses(D_c, P_max, T_max)
        L_c = L_x / np.cos(np.pi * alpha / 180)
        N_p = int(np.floor((L_y - t_3) / (D_c / 2 + t_3)))
        N_c = int(np.floor((L_z - t_2) / (D_c + t_2)))

        V_channel = np.pi * D_c**2 / 8 * L_c
        R_p = self.params['R_p']  # lecture seule
        C_V_tot = 1 / (1 + R_p) * N_p * N_c * V_channel
        H_V_tot = R_p / (1 + R_p) * N_p * N_c * V_channel

        # ---- feasibility ----
        V_block = L_x * L_y * L_z
        if (C_V_tot + H_V_tot) >= V_block:
            return 1e6

        # ---- build LOCAL param dict (ne touche jamais self.params) ----
        local_params = {
            **self.params,
            'alpha': alpha, 'D_c': D_c, 'L_x': L_x, 'L_y': L_y, 'L_z': L_z,
            'n_parallel': n_parallel, 'n_series': n_series,
            't_2': t_2, 't_3': t_3, 'L_c': L_c,
            'N_p': N_p, 'N_c': N_c,
            'C_V_tot': C_V_tot, 'H_V_tot': H_V_tot,
        }

        # ---- simulate on a LOCAL HX (propre à ce thread/cet appel) ----
        HX = HexMBChargeSensitive('PCHE')
        HX.set_inputs(**self.inputs)
        self._apply_corr(HX)              # ne mute que ce HX local
        HX.set_parameters(**local_params)

        try:
            HX.solve()
        except Exception:
            return 1e12

        rho_mat = 7850  # kg/m^3
        m_HX = n_parallel * rho_mat * (L_x * L_y * L_z - (C_V_tot + H_V_tot))

        result = self.compute_score(
            HX, m_HX, self.Q_dot_constr, self.DP_h_constr, self.DP_c_constr,
            self.params.get('obj', 'mass'),
        )

        return result['score']

    #%%

    def sizing(self, n_jobs=-1, backend="threading", chunksize="auto", n_particles = 30, max_iter = None, patience = 10, obj = "cost"):
        self._apply_deferred_parameters()
        
        self.params["obj"] = obj
        
        # ---- fixed order + bounds ----
        ORDER = ['alpha', 'D_c', 'L_x', 'L_y', 'L_z', 'n_parallel', 'n_series']
        def bounds_dict_to_arrays(bounds_dict, order=ORDER):
            lb = np.array([bounds_dict[k][0] for k in order], dtype=float)
            ub = np.array([bounds_dict[k][1] for k in order], dtype=float)
            if np.any(lb > ub):
                raise ValueError("Lower bound > upper bound for at least one variable.")
            return lb, ub

        lb, ub = bounds_dict_to_arrays(self.bounds, ORDER)   # <-- pass ORDER
        D = lb.size

        optimizer = ps.single.GlobalBestPSO(
            n_particles=n_particles, dimensions=D,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=(lb, ub)
        )

        tol = 1e-3

        if max_iter is None:
            max_iter = patience * 10

        # after creating `optimizer`
        try:
            n_particles = int(optimizer.swarm.n_particles)
        except AttributeError:
            # fallback if API differs
            n_particles = int(getattr(optimizer, "n_particles", len(optimizer.swarm.position)))

        total_evals = n_particles * int(max_iter)
        # pbar = tqdm(total=total_evals, desc="Function evaluations", unit="eval")

        def objective_wrapper(X):
            Xp = [np.clip(xi, lb, ub) for xi in np.asarray(X)]
            costs = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
                delayed(self.simulate_HX)(xi) for xi in Xp
            )
            # one update per PSO iteration
            pbar.update(len(Xp))
            return np.asarray(costs, dtype=float)

        no_improve, best_cost = 0, np.inf

        pbar = tqdm(total=total_evals, unit="eval", dynamic_ncols=True, desc='PCHE')

        for i in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            curr = optimizer.swarm.best_cost

            pbar.clear()   # remove previous bar
            pbar.set_description(f"Iter {i+1} | Best: {curr:.3f}")
            pbar.refresh()

            if curr < best_cost - tol:
                best_cost = curr
                no_improve = 0
            else:
                no_improve += 1

            if no_improve >= patience:
                break

        best_pos = optimizer.swarm.best_pos

        n_series_best = max(1, np.round(best_pos[6]))

        # ---- Safe single-threaded final rebuild ----
        self.params.update({
            'alpha': best_pos[0], 'D_c': best_pos[1],
            'L_x': best_pos[2],   'L_y': best_pos[3], 'L_z': best_pos[4],
            'n_parallel': np.round(best_pos[5]),
            'n_series': n_series_best,
        })
        self.compute_geom()

        self.HX = HexMBChargeSensitive('PCHE')
        self.HX.set_inputs(**self.inputs)
        self._apply_corr(self.HX)
        self.HX.set_parameters(**self.params)

        try:
            self.HX.solve()
        except Exception as e:
            raise RuntimeError(
                f"PCHE sizing did not converge to a feasible geometry "
                f"(best_cost={best_cost:.3e}); final solve failed with: {e}"
            ) from e

        self.m_HX = np.round(best_pos[5]) * n_series_best * 7850 * (
            self.params['L_x'] * self.params['L_y'] * self.params['L_z']
            - self.params['C_V_tot'] - self.params['H_V_tot']
        )

        # ---- SEUL endroit du fichier où compute_score() est appelée puis
        # son résultat assigné sur self -- mono-thread, donc sûr. On force
        # obj='cost' ici pour que self.CAPEX soit toujours peuplé au final,
        # même si l'optimisation a tourné avec obj='mass' (remplace l'appel
        # séparé à l'ancien cost_estimation() en fin de sizing()). ----
        result = self.compute_score(
            self.HX, self.m_HX, self.Q_dot_constr, self.DP_h_constr, self.DP_c_constr,
            obj='cost',
        )

        self.score = result['score']
        self.penalty = result['penalty']
        self.pen_Q = result['pen_Q']
        self.pen_DP_h = result['pen_DP_h']
        self.pen_DP_c = result['pen_DP_c']
        self.UA = result['UA']
        self.CAPEX_alt = result['CAPEX_alt']
        self.CAPEX = result['CAPEX']
        self.capex_score = result['CAPEX']['Total'] if result['CAPEX'] else None

        pbar.close()

        if __name__ == "__main__":
            print("\nBest Position")
            print("-------------")
            print(f"alpha : {round(self.params['alpha'],2)} [°]")
            print(f"D_c : {round(self.params['D_c']*1e3,2)} [mm]")
            print(f"L_x : {round(self.params['L_x'],3)} [m]")
            print(f"L_y : {round(self.params['L_y'],3)} [m]")
            print(f"L_z : {round(self.params['L_z'],3)} [m]")
            print(f"n_parallel : {self.params['n_parallel']} [-]")
            print(f"n_series : {self.params['n_series']} [-]")

            print("\nResults")
            print("-------------")
            print(f"A_h : {round(self.HX.A_h,2)} [m^2]")
            print(f"A_c : {round(self.HX.A_c,2)} [m^2]")
            print(f"Q_dot : {round(self.HX.Q.Q_dot,1)} [W]")
            print(f"DP_c : {round(self.HX.DP_c,1)} [Pa]")
            print(f"DP_h : {round(self.HX.DP_h,1)} [Pa]")
            print(f"m_HX : {round(self.m_HX,1)} [kg]")
            print(f"penalty : {round(self.penalty,1)}")
            print(f"CAPEX est. : {round(self.CAPEX['Total'],1)} [€ (2026)]")
            
            print(f"\n--- Diagnostic ---")
            print(f"Q_dot atteint  : {HX_opt.HX.Q.Q_dot:,.0f} W")
            print(f"Q_dot cible    : {HX_opt.Q_dot_constr:,.0f} W")
            print(f"pen_Q          : {HX_opt.pen_Q:,.1f}")
            print(f"pen_DP_h       : {HX_opt.pen_DP_h:,.1f}")
            print(f"pen_DP_c       : {HX_opt.pen_DP_c:,.1f}")
            print(f"penalty totale : {HX_opt.penalty:,.1f}")
            print(f"L_z retenu     : {HX_opt.params['L_z']:.4f}  (bornes {HX_opt.bounds['L_z']})")
            print(f"N_c retenu     : {HX_opt.params['N_c']}")
            print(f"N_p retenu     : {HX_opt.params['N_p']}")
            print(f"t_2/t_3        : {HX_opt.params['t_2']:.6f} / {HX_opt.params['t_3']:.6f}")

        return best_pos

#%%

if __name__ == "__main__":
    HX_opt = PCHESizingOpt()

    case_study = 'REC3'

    if case_study == "Reference":
        HX_opt.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 249 + 273.15, # K
            P_su_H = 96.4*1e5, # Pa
            m_dot_H = 5.35, # kg/s

            # Second fluid
            fluid_C = 'CO2',
            T_su_C = 52.77 + 273.15, # K
            P_su_C = 165.4*1e5, # Pa
            m_dot_C = 5.35, # kg/s  # Make sure to include fluid information
          )
        
        
        Q_dot_cstr = 6e5
        DP_h_cstr = 400000.0
        DP_c_cstr = 200000.0

    elif case_study == "REC":
        HX_opt.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 342.22, # K
            P_su_H = 6047124.7, # Pa
            m_dot_H = 373.62, # kg/s

            # Second fluid
            fluid_C = 'CO2',
            T_su_C = 303.078, # K
            P_su_C = 13260336.85, # Pa
            m_dot_C = 373.62, # kg/s  # Make sure to include fluid information
          )

    elif case_study == "REC2":
        HX_opt.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 329.3, # K
            P_su_H = 6196195.1, # Pa
            m_dot_H = 381.3, # kg/s

            # Second fluid
            fluid_C = 'CO2',
            T_su_C = 306.8, # K
            P_su_C = 15364581.8, # Pa
            m_dot_C = 381.3, # kg/s  # Make sure to include fluid information
          )

    elif case_study == "REC3":
        # Cas issu du print_resume() lors d'un échec de sizing en conditions réelles
        # (via TCO2_rec_comp_sizing) — servant de reproduction isolée pour le debug.
        HX_opt.set_inputs(
            # First fluid (hot, low pressure side)
            fluid_H = 'CO2',
            T_su_H = 316.1313920137559, # K
            P_su_H = 5336196.6651766505,  # Pa
            m_dot_H = 2819.158278285053, # kg/s

            # Second fluid (cold, high pressure side)
            fluid_C = 'CO2',
            T_su_C = 302.5694587811514, # K
            P_su_C = 15590873.653394502,  # Pa
            m_dot_C = 2819.158278285053, # kg/s  # Make sure to include fluid information
          )

        Q_dot_cstr = 53786439.826693796
        DP_h_cstr = 28635.0
        DP_c_cstr = 25732.5

    # Grâce à DEFAULT_PARAMETERS, seuls les paramètres qui diffèrent des
    # valeurs par défaut (n_disc, corrélations, contraintes) doivent être
    # passés ici. k_cond, R_p, Flow_Type, H_DP_ON, C_DP_ON sont déjà couverts.
    HX_opt.set_parameters(
        # n_disc = 30,

        H_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"},
        C_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"},
        H_DP   = {"1P" : "Gnielinski_DP", "SC" : "Gnielinski_DP"},
        C_DP   = {"1P" : "Gnielinski_DP", "SC" : "Gnielinski_DP"},

        Q_dot = Q_dot_cstr,
        DP_h  = DP_h_cstr,
        DP_c  = DP_c_cstr,
        )

    # Idem pour les bounds : DEFAULT_BOUNDS couvre déjà les 7 clés ; ici on
    # ne surcharge que D_c pour ce cas d'étude.
    # HX_opt.set_bounds(
    #     D_c = [2.5*1e-3, 3*1e-3],
    #     )

    best_pos = HX_opt.sizing(n_jobs=1, n_particles=50, patience = 10)