# -*- coding: utf-8 -*-
"""
co2_rc_pso_optimizer.py
CO2 Transcritical Rankine Cycle Optimizer — brique d'optimisation (PSO)
- Parallel PSO evaluation via joblib
- Support de 4 architectures : 'basic', 'REC', 'Recomp', 'Recomp_1_recup'
- Warm start

Ce module est destiné à être importé (ex: par co2_rc_full_design_optimizer.py)
qui ajoute le dimensionnement des composants par-dessus.

CORRECTIONS PAR RAPPORT A LA VERSION PRECEDENTE
------------------------------------------------
1. `_evaluate_final()` calcule maintenant explicitement l'effectivité
   (epsilon) de chaque échangeur après le solve de la meilleure position
   PSO, en déduit un NTU équivalent (NTU = -ln(1-eps), cohérent avec la
   fonction objectif utilisée dans `system_RC_parallel`), et stocke tout
   cela en attributs sur l'instance (`self.NTU_gh`, `self.NTU_rec`,
   `self.NTU_rec_LT`, `self.NTU_rec_HT`, `self.NTU_cd`,
   `self.NTU_weighted_sum`). Avant, ces valeurs n'étaient jamais calculées
   après le solve final, donc `extract_run_data` (dans le bloc __main__)
   les lisait toujours comme `None`.
2. Les facteurs de débit (`m_dot_HS_factor`, `m_dot_CS_factor`) sont
   maintenant conservés (`self.m_dot_HS_factor`, `self.m_dot_CS_factor`)
   au lieu d'être perdus après le calcul de `mdot_HS` / `mdot_CS`.
3. `P_high`, `m_dot`, `eta_gh`, `PP_gh`, `eta_rec`, `PP_cd`, `mdot_HS`,
   `mdot_CS`, `spliter_frac` sont aussi exposés en attributs plats
   (`self.P_high`, `self.m_dot`, ...) en plus du dict `self.it_var`, pour
   que `getattr(Optimizer, "P_high", None)` etc. fonctionnent réellement.
4. `W_net` est maintenant disponible sous le nom `self.W_dot_net` (déjà
   correct avant) ET sous l'alias `self.W_net` pour compatibilité avec
   `extract_run_data`.
5. `extract_run_data` (bloc __main__) a été réécrit pour être cohérent
   avec ces attributs et ne plus deviner des noms qui n'existaient pas.
"""

#%% Imports

from labothappy.machine.examples.ORC.fpi_TC_orc_example import (
    REC_CO2_TC, basic_CO2_TC, Recomp_CO2_TC, Recomp_CO2_TC_1_recup
)
from labothappy.connector.mass_connector import MassConnector

import numpy as np
from CoolProp.CoolProp import PropsSI
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm
from joblib import Parallel, delayed
import multiprocessing

import warnings
warnings.filterwarnings('ignore')

#%% Helper : effectivité -> NTU équivalent

def _epsilon_to_ntu(eps, eps_cap=1e-6):
    """NTU équivalent à partir d'une effectivité, cohérent avec la fonction
    objectif utilisée dans system_RC_parallel (-ln(1-eps))."""
    if eps is None:
        return None
    eps = float(np.clip(eps, 0.0, 1.0 - eps_cap))
    return -np.log(1.0 - eps)


def _lmtd_counterflow(T_h_in, T_h_out, T_c_in, T_c_out):
    """LMTD à contre-courant. Retourne None si croisement de température ou
    valeur non physique (protège contre les logs de nombres négatifs)."""
    dT1 = T_h_in - T_c_out
    dT2 = T_h_out - T_c_in
    if dT1 <= 0 or dT2 <= 0:
        return None
    if abs(dT1 - dT2) < 1e-6:
        return dT1
    return (dT1 - dT2) / np.log(dT1 / dT2)


def _estimate_UA_lmtd(model):
    """
    Estime UA = Q / LMTD à partir des températures terminales du modèle
    thermodynamique DÉJÀ résolu — ne nécessite aucun sizing détaillé.
    Approximatif : suppose des propriétés constantes le long de l'échangeur.
    Weiland & Lance (2019) rapportent jusqu'à ~80% d'erreur avec cette
    hypothèse pour un LTR proche du point critique du CO2 — d'où la
    calibration empirique (voir calibrate_ua_correction dans
    co2_rc_full_design_optimizer.py) qui corrige ce biais avec de vraies UA
    issues du sizing, avant d'utiliser ce proxy pour l'objectif du PSO.

    À VÉRIFIER : suppose que le modèle expose `.ex_H`/`.ex_C` (sorties
    chaude/froide), par analogie avec `.su_H`/`.su_C` déjà utilisés ailleurs
    dans ce module. Si le nom réel diffère dans tes modèles, adapter cette
    fonction — elle est conçue pour échouer proprement (retourne None) plutôt
    que de lever une exception si un attribut est absent.
    """
    try:
        Q = model.Q.Q_dot
        T_h_in  = model.su_H.T
        T_h_out = model.ex_H.T
        T_c_in  = model.su_C.T
        T_c_out = model.ex_C.T
        lmtd = _lmtd_counterflow(T_h_in, T_h_out, T_c_in, T_c_out)
        if lmtd is None or lmtd <= 0 or Q is None or Q <= 0:
            return None
        return Q / lmtd
    except Exception:
        return None


def _weiland_lance_cost(UA, T_max_C, a=49.45, b=0.7544, T_bp=550.0, c=0.02141):
    """
    Corrélation de coût récupérateur de Weiland & Lance (2019), réutilisée
    pour les 3 échangeurs liquide/CO2 de ce cycle (GasHeater, Recuperator,
    Condenser) — cf. discussion : le papier lui-même utilise ce modèle comme
    proxy pour un refroidisseur eau/sCO2 en l'absence de corrélation dédiée,
    et ses corrélations "primary heater" (charbon/gaz) ne s'appliquent pas à
    une source chaude liquide comme ici.
    UA en W/K, T_max_C = température max de fonctionnement en °C.
    Retourne le coût en $ (base 2017$, cf. papier).
    """
    if UA is None or UA <= 0:
        return None
    fT = 1.0
    if T_max_C >= T_bp:
        fT = 1.0 + c * (T_max_C - T_bp)
    return a * (UA ** b) * fT


#%% Top-level parallel evaluation function
# Doit rester au niveau module pour être picklable par joblib/loky.

def system_RC_parallel(x, input_data):
    """
    Évalue une particule x. Retourne (cost, penalty, eta).
    cost élevé positif = infaisable.
    """
    warnings.filterwarnings('ignore')

    fluid    = input_data['fluid']
    params   = input_data['params']
    obj      = input_data['obj']
    hs_props = input_data['HSource']
    cs_props = input_data['CSource']
    arch     = input_data['RC_ARCH']

    if arch == "Recomp":
        P_high, m_dot, m_dot_HS_fact, spliter_frac, eta_gh, PP_gh, eta_rec_LT, eta_rec_HT, PP_cd, m_dot_CS_fact = x
    elif arch == "Recomp_1_recup":
        P_high, m_dot, m_dot_HS_fact, spliter_frac, eta_gh, PP_gh, eta_rec, PP_cd, m_dot_CS_fact = x
    elif arch == 'REC':
        P_high, m_dot, m_dot_HS_fact, eta_gh, PP_gh, eta_rec, PP_cd, m_dot_CS_fact = x
    elif arch == 'basic':
        P_high, m_dot, m_dot_HS_fact, eta_gh, PP_gh, PP_cd, m_dot_CS_fact = x
    else:
        return 1000.0, np.inf, np.nan

    m_dot_HS = m_dot * m_dot_HS_fact
    m_dot_CS = m_dot * m_dot_CS_fact

    HSource = MassConnector()
    HSource.set_properties(T=hs_props['T'], P=hs_props['P'],
                            fluid=hs_props['fluid'], m_dot=m_dot_HS)

    CSource = MassConnector()
    CSource.set_properties(T=cs_props['T'], P=cs_props['P'],
                            fluid=cs_props['fluid'], m_dot=m_dot_CS)

    P_sat_CS    = input_data.get('P_sat_CS')
    P_crit      = input_data.get('P_crit')
    if P_sat_CS is None:
        P_sat_CS = PropsSI('P', 'T', cs_props['T'], 'Q', 0.5, fluid)
    if P_crit is None:
        P_crit = PropsSI('PCRIT', fluid)
    P_low_guess = min(1.3 * P_sat_CS, 0.8 * P_crit)

    if arch == 'REC':
        RC = REC_CO2_TC(
            HSource, CSource,
            PP_gh, params['PP_rec'], params['eta_pp'],
            params['eta_exp'], eta_gh, eta_rec,
            PP_cd, params['SC_cd'],
            P_low_guess, P_high, m_dot,
            DP_h_rec  = params.get('DP_h_rec',  0e5),
            DP_c_rec  = params.get('DP_c_rec',  0e5),
            DP_h_gh   = params.get('DP_h_gh',   0e5),
            DP_c_gh   = params.get('DP_c_gh',   0e5),
            DP_h_cond = params.get('DP_h_cond', 0e5),   
            DP_c_cond = params.get('DP_c_cond', 0e5),  
            mute_print_flag=1)

    elif arch == 'basic':
        RC = basic_CO2_TC(
            HSource, CSource,
            PP_gh, params['eta_pp'], params['eta_exp'],
            eta_gh, PP_cd, params['SC_cd'],
            P_low_guess, P_high, m_dot,
            DP_h_gh   = params.get('DP_h_gh',   0e5),
            DP_c_gh   = params.get('DP_c_gh',   0e5),
            DP_h_cond = params.get('DP_h_cond', 0e5),
            DP_c_cond = params.get('DP_c_cond', 0e5),  
            mute_print_flag=1)

    elif arch == "Recomp":
        RC = Recomp_CO2_TC(
            HSource, CSource,
            PP_gh, params['PP_rec'], params['eta_pp'],
            params['eta_exp'], params['eta_cp'],
            eta_gh, eta_rec_LT, eta_rec_HT,
            PP_cd, params['SC_cd'],
            P_low_guess, P_high, m_dot, spliter_frac,
            DP_h_rec  = params.get('DP_h_rec',  0e5),
            DP_c_rec  = params.get('DP_c_rec',  0e5),
            DP_h_gh   = params.get('DP_h_gh',   0e5),
            DP_c_gh   = params.get('DP_c_gh',   0e5),
            DP_h_cond = params.get('DP_h_cond', 0e5),   
            DP_c_cond = params.get('DP_c_cond', 0e5),  
            mute_print_flag=1)

    elif arch == "Recomp_1_recup":
        RC = Recomp_CO2_TC_1_recup(
            HSource, CSource,
            PP_gh, params['PP_rec'],
            params['eta_pp'], params['eta_exp'], params['eta_cp'],
            eta_rec, eta_gh,
            PP_cd, params['SC_cd'],
            P_low_guess, P_high, m_dot, spliter_frac,
            DP_h_rec  = params.get('DP_h_rec',  0e5),
            DP_c_rec  = params.get('DP_c_rec',  0e5),
            DP_h_gh   = params.get('DP_h_gh',   0e5),
            DP_c_gh   = params.get('DP_c_gh',   0e5),
            DP_h_cond = params.get('DP_h_cond', 0e5),   
            DP_c_cond = params.get('DP_c_cond', 0e5),  
            mute_print_flag=1)

    try:
        RC.solve()

        if not getattr(RC, 'converged', True):
            return 10000, np.inf, np.nan

        if arch in ("Recomp", "Recomp_1_recup"):
            W_cp = RC.components['Compressor'].model.W.W_dot
        else:
            W_cp = 0

        W_exp  = RC.components['Expander'].model.W.W_dot
        W_pump = RC.components['Pump'].model.W.W_dot
        Q_gh   = RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = RC.components['GasHeater'].model.su_H.D
        m_HS_act   = RC.components['GasHeater'].model.su_H.m_dot

        rho_CS     = RC.components['Condenser'].model.su_C.D
        m_CS_act   = RC.components['Condenser'].model.su_C.m_dot

        W_pump_aux_HS = params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * params.get('eta_pp_aux', 0.8))

        W_pump_aux_CS = params.get('DP_c_cond', 0.5e5) * m_CS_act / \
                     (rho_CS * params.get('eta_pp_aux', 0.8))

        W_dot_net = W_exp - W_pump - W_pump_aux_HS - W_pump_aux_CS - W_cp
        eta       = W_dot_net / Q_gh if Q_gh > 0 else 0.0

        penalty_W_dot = 0
        if abs(obj['W_dot'] - W_dot_net) / obj['W_dot'] > 2e-2:
            penalty_W_dot = abs(obj['W_dot'] - W_dot_net) / obj['W_dot']

        penalty_eta = 0
        if abs(obj['eta'] - eta) / obj['eta'] > 1e-2:
            penalty_eta = abs(obj['eta'] - eta) / obj['eta']

        RC.eta = eta
        RC.W_dot_net = W_dot_net

        eps = 1e-6

        if arch == 'REC':
            Q_cond = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon

            Q_rec = RC.components['Recuperator'].model.Q.Q_dot
            eta_rec = RC.components['Recuperator'].model.epsilon

            Q_gh = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh = RC.components['GasHeater'].model.epsilon

            # Pénalité si le récupérateur est quasi bypassé (Q_rec / eta_rec
            # réel proche de 0) : sinon le PSO peut "tricher" en désactivant
            # le récupérateur pour faire baisser artificiellement le NTU
            # pondéré (son terme disparaît du numérateur ET du dénominateur).
            eps_rec_min = params.get('eps_rec_min', 0.3)
            penalty_rec = max(0.0, eps_rec_min - eta_rec)

            eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
            eta_rec = np.clip(eta_rec, 0.0, 1.0 - eps)
            eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)

            # Poids de coût relatif par technologie d'échangeur (PCHE pour le
            # récupérateur vs Shell&Tube pour GasHeater/Condenser). Défaut=1
            # -> comportement d'origine (NTU pondéré par Q_dot, en assumant
            # un $/NTU identique entre technologies). À calibrer depuis des
            # données CAPEX réelles si les coûts par NTU diffèrent nettement
            # (cf. calibrate_cost_weights dans co2_rc_full_design_optimizer.py).
            c_gh   = params.get('cost_w_gh', 1.0)
            c_rec  = params.get('cost_w_rec', 1.0)
            c_cond = params.get('cost_w_cond', 1.0)
            
            # ------------------------------------------------------------
            # STAGE 1 vs STAGE 2 : voir discussion "optimiser les NTU d'abord,
            # remplacer par le vrai coût une fois les UA calculées".
            #
            # Stage 1 (cost_model_calibrated=False, défaut) : proxy NTU
            # pondéré par Q_dot (comportement d'origine, éventuellement
            # repondéré par cost_w_*).
            #
            # Stage 2 (cost_model_calibrated=True, activé par
            # cycle_design() après un premier size_components() réussi) :
            # coût réel via Weiland & Lance, C = a·UA^b·f_T(T_max), où UA est
            # estimée par LMTD (thermo seule, pas de sizing dans la boucle)
            # puis corrigée par un facteur k_x calibré sur les vraies UA
            # issues du sizing (voir calibrate_ua_correction).
            # ------------------------------------------------------------
            cost_model_calibrated = params.get('cost_model_calibrated', False)

            objective = None
            if cost_model_calibrated:
                k_gh   = params.get('ua_correction_gh', 1.0)
                k_rec  = params.get('ua_correction_rec', 1.0)
                k_cond = params.get('ua_correction_cond', 1.0)

                UA_gh   = _estimate_UA_lmtd(RC.components['GasHeater'].model)
                UA_rec  = _estimate_UA_lmtd(RC.components['Recuperator'].model)
                UA_cond = _estimate_UA_lmtd(RC.components['Condenser'].model)

                if UA_gh is not None and UA_rec is not None and UA_cond is not None:
                    T_max_gh   = RC.components['GasHeater'].model.su_H.T - 273.15
                    T_max_rec  = RC.components['Recuperator'].model.su_H.T - 273.15
                    T_max_cond = RC.components['Condenser'].model.su_H.T - 273.15

                    cost_gh   = _weiland_lance_cost(UA_gh * k_gh, T_max_gh)
                    cost_rec  = _weiland_lance_cost(UA_rec * k_rec, T_max_rec)
                    cost_cond = _weiland_lance_cost(UA_cond * k_cond, T_max_cond)

                    if None not in (cost_gh, cost_rec, cost_cond):
                        scale = params.get('cost_obj_scale', 1e6)  # $ -> échelle comparable au NTU pondéré / à PF
                        objective = (cost_gh + cost_rec + cost_cond) / scale

            if objective is None:
                # Stage 1 (ou repli si UA/LMTD non estimable ce coup-ci,
                # ex. croisement de température) : NTU pondéré d'origine.
                objective = (c_gh*Q_gh*(-np.log(1-eta_gh)) + c_rec*Q_rec*(-np.log(1-eta_rec)) + c_cond*Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta + penalty_rec)*PF
            cost = objective + penalty

        elif arch == 'basic':
            Q_cond = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon

            Q_gh = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh = RC.components['GasHeater'].model.epsilon

            eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
            eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)

            c_gh   = params.get('cost_w_gh', 1.0)
            c_cond = params.get('cost_w_cond', 1.0)

            objective = (c_gh*Q_gh*(-np.log(1-eta_gh)) + c_cond*Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta)*PF
            cost = objective + penalty

        elif arch == "Recomp":
            Q_cond = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon

            Q_rec_LT = RC.components['RecupLT'].model.Q.Q_dot
            eta_rec_LT = RC.components['RecupLT'].model.epsilon

            Q_rec_HT = RC.components['RecupHT'].model.Q.Q_dot
            eta_rec_HT = RC.components['RecupHT'].model.epsilon

            Q_gh = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh = RC.components['GasHeater'].model.epsilon

            eps_rec_min = params.get('eps_rec_min', 0.3)
            penalty_rec = (max(0.0, eps_rec_min - eta_rec_LT)
                           + max(0.0, eps_rec_min - eta_rec_HT))

            eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
            eta_rec_LT = np.clip(eta_rec_LT, 0.0, 1.0 - eps)
            eta_rec_HT = np.clip(eta_rec_HT, 0.0, 1.0 - eps)
            eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)

            c_gh     = params.get('cost_w_gh', 1.0)
            c_rec_lt = params.get('cost_w_rec_LT', 1.0)
            c_rec_ht = params.get('cost_w_rec_HT', 1.0)
            c_cond   = params.get('cost_w_cond', 1.0)

            objective = (c_gh*Q_gh*(-np.log(1-eta_gh)) + c_rec_lt*Q_rec_LT*(-np.log(1-eta_rec_LT)) + c_rec_ht*Q_rec_HT*(-np.log(1-eta_rec_HT)) + c_cond*Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec_LT + Q_rec_HT + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta + penalty_rec)*PF
            cost = objective + penalty

        elif arch == "Recomp_1_recup":
            Q_cond = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon

            Q_rec_LT = RC.components['RecupLT'].model.Q.Q_dot
            eta_rec_LT = RC.components['RecupLT'].model.epsilon

            Q_gh = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh = RC.components['GasHeater'].model.epsilon

            eps_rec_min = params.get('eps_rec_min', 0.3)
            penalty_rec = max(0.0, eps_rec_min - eta_rec_LT)

            eta_gh     = np.clip(eta_gh,     0.0, 1.0 - eps)
            eta_rec_LT = np.clip(eta_rec_LT, 0.0, 1.0 - eps)
            eta_cond   = np.clip(eta_cond,   0.0, 1.0 - eps)

            c_gh     = params.get('cost_w_gh', 1.0)
            c_rec_lt = params.get('cost_w_rec_LT', 1.0)
            c_cond   = params.get('cost_w_cond', 1.0)

            objective = (c_gh*Q_gh*(-np.log(1-eta_gh)) + c_rec_lt*Q_rec_LT*(-np.log(1-eta_rec_LT)) + c_cond*Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec_LT + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta + penalty_rec)*PF
            cost = objective + penalty

    except Exception:
        return 20000.0, np.inf, np.nan

    return cost, penalty, eta

#%% Optimizer Class (base)

class CO2RC_Cost_optimizer:

    def __init__(self, fluid):
        self.fluid  = fluid
        self.RC     = None

        self.inputs  = {}
        self.params  = {}
        self.it_var  = {}
        self.obj     = {}

        self._HSource_props = {}
        self._CSource_props = {}

        self.eta         = None
        self.W_dot_net   = None
        self.W_net       = None  # alias de W_dot_net, pour compatibilité
        self.penalty_log = {}
        self.allowable_positions = []
        self.top_positions = []

        # --- Attributs plats exposés après _evaluate_final() ---
        # (en plus de it_var, pour que getattr(Optimizer, "P_high", ...)
        # fonctionne vraiment, cf. extract_run_data dans le bloc __main__)
        self.P_high         = None
        self.m_dot           = None
        self.mdot            = None
        self.m_dot_HS        = None
        self.m_dot_CS        = None
        self.m_dot_HS_factor = None
        self.m_dot_CS_factor = None
        self.spliter_frac    = None
        self.eta_gh          = None
        self.PP_gh           = None
        self.eta_rec         = None
        self.eta_rec_LT      = None
        self.eta_rec_HT      = None
        self.PP_cd           = None

        # --- NTU calculés après le solve final (voir _evaluate_final) ---
        self.NTU_gh           = None
        self.NTU_rec           = None
        self.NTU_rec_LT        = None
        self.NTU_rec_HT        = None
        self.NTU_cd             = None
        self.NTU_weighted_sum  = None

        # --- Effectivités (epsilon) associées, calculées en même temps ---
        self.eps_gh      = None
        self.eps_rec     = None
        self.eps_rec_LT  = None
        self.eps_rec_HT  = None
        self.eps_cd      = None

    # ------------------------------------------------------------------ setters

    def set_inputs(self, **parameters):
        self.inputs.update(parameters)

    def set_parameters(self, **parameters):
        self.params.update(parameters)

    def set_it_var(self, **parameters):
        self.it_var.update(parameters)

    def set_obj(self, **parameters):
        self.obj.update(parameters)

    def set_HSource(self, T, P, fluid, m_dot=1.0):
        self._HSource_props = dict(T=T, P=P, fluid=fluid, m_dot=m_dot)

    def set_CSource(self, T, P, fluid, m_dot=1000.0):
        self._CSource_props = dict(T=T, P=P, fluid=fluid, m_dot=m_dot)

    # ------------------------------------------------------------------ helper : déballage d'une position PSO

    def _unpack_position(self, x):
        """Traduit un vecteur PSO x en dict de variables selon l'architecture.
        Conserve aussi les facteurs bruts (m_dot_HS_fact, m_dot_CS_fact) qui
        étaient auparavant perdus après le calcul de mdot_HS / mdot_CS."""
        arch = self.params.get('RC_ARCH', 'REC')
        x = np.asarray(x, dtype=float)

        if arch == "Recomp":
            P_high, m_dot, m_dot_HS_fact, spliter_frac, eta_gh, PP_gh, eta_rec_LT, eta_rec_HT, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        spliter_frac=spliter_frac, eta_gh=eta_gh, PP_gh=PP_gh,
                        eta_rec_LT=eta_rec_LT, eta_rec_HT=eta_rec_HT, PP_cd=PP_cd,
                        mdot_CS=m_dot * m_dot_CS_fact,
                        m_dot_HS_factor=m_dot_HS_fact, m_dot_CS_factor=m_dot_CS_fact)

        elif arch == "Recomp_1_recup":
            P_high, m_dot, m_dot_HS_fact, spliter_frac, eta_gh, PP_gh, eta_rec, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        spliter_frac=spliter_frac, eta_gh=eta_gh, PP_gh=PP_gh,
                        eta_rec=eta_rec, PP_cd=PP_cd, mdot_CS=m_dot * m_dot_CS_fact,
                        m_dot_HS_factor=m_dot_HS_fact, m_dot_CS_factor=m_dot_CS_fact)

        elif arch == "REC":
            P_high, m_dot, m_dot_HS_fact, eta_gh, PP_gh, eta_rec, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        eta_gh=eta_gh, PP_gh=PP_gh, eta_rec=eta_rec, PP_cd=PP_cd,
                        mdot_CS=m_dot * m_dot_CS_fact,
                        m_dot_HS_factor=m_dot_HS_fact, m_dot_CS_factor=m_dot_CS_fact)

        elif arch == "basic":
            P_high, m_dot, m_dot_HS_fact, eta_gh, PP_gh, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        eta_gh=eta_gh, PP_gh=PP_gh, PP_cd=PP_cd,
                        mdot_CS=m_dot * m_dot_CS_fact,
                        m_dot_HS_factor=m_dot_HS_fact, m_dot_CS_factor=m_dot_CS_fact)

        else:
            raise ValueError("'RC_ARCH' shall be 'basic', 'REC', 'Recomp' or 'Recomp_1_recup'")

    # ------------------------------------------------------------------ RC build

    def set_RC(self):
        """Builds self.RC from current it_var and source props."""
        HSource = MassConnector()
        HSource.set_properties(
            T     = self._HSource_props['T'],
            P     = self._HSource_props['P'],
            fluid = self._HSource_props['fluid'],
            m_dot = self.it_var['mdot_HS'],
        )

        CSource = MassConnector()
        CSource.set_properties(
            T     = self._CSource_props['T'],
            P     = self._CSource_props['P'],
            fluid = self._CSource_props['fluid'],
            m_dot = self.it_var['mdot_CS'],
        )

        P_sat_CS    = PropsSI('P', 'T', self._CSource_props['T'], 'Q', 0.5, self.fluid)
        P_crit      = PropsSI('PCRIT', self.fluid)
        P_low_guess = min(1.3 * P_sat_CS, 0.8 * P_crit)

        arch = self.params.get('RC_ARCH', 'REC')

        if arch == 'REC':
            self.RC = REC_CO2_TC(
                HSource, CSource,
                self.it_var['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'],
                self.it_var['eta_gh'], self.it_var['eta_rec'],
                self.it_var['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'],
                DP_h_rec  = self.params.get('DP_h_rec',  0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond',   0e5),
                DP_c_cond = self.params.get('DP_c_cond',   0e5),
                mute_print_flag=1,
            )
        elif arch == 'basic':
            self.RC = basic_CO2_TC(
                HSource, CSource,
                self.it_var['PP_gh'], self.params['eta_pp'],
                self.params['eta_exp'], self.it_var['eta_gh'],
                self.it_var['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'],
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond',   0e5),
                DP_c_cond = self.params.get('DP_c_cond',   0e5),
                mute_print_flag=1,
            )
        elif arch == 'Recomp':
            self.RC = Recomp_CO2_TC(
                HSource, CSource,
                self.it_var['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'], self.params['eta_cp'],
                self.it_var['eta_gh'], self.it_var['eta_rec_LT'], self.it_var['eta_rec_HT'],
                self.it_var['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'], self.it_var['spliter_frac'],
                DP_h_rec  = self.params.get('DP_h_rec',  0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond',   0e5),
                DP_c_cond = self.params.get('DP_c_cond',   0e5),
                mute_print_flag=1)
        elif arch == 'Recomp_1_recup':
            self.RC = Recomp_CO2_TC_1_recup(
                HSource, CSource,
                self.it_var['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'], self.params['eta_cp'],
                self.it_var['eta_rec'], self.it_var['eta_gh'],
                self.it_var['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'], self.it_var['spliter_frac'],
                DP_h_rec  = self.params.get('DP_h_rec',  0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond',   0e5),
                DP_c_cond = self.params.get('DP_c_cond',   0e5),
                mute_print_flag=1)
        else:
            raise ValueError("'RC_ARCH' parameter shall be either 'basic', 'REC', 'Recomp', 'Recomp_1_recup'")

    def _log_penalty(self, reason):
        self.penalty_log[reason] = self.penalty_log.get(reason, 0) + 1

    # ------------------------------------------------------------------ NTU (nouveau)

    def _compute_final_ntu(self):
        """Recalcule l'effectivité de chaque échangeur sur self.RC (déjà
        résolu avec la meilleure position) et en déduit un NTU équivalent,
        cohérent avec la fonction objectif de system_RC_parallel
        (NTU = -ln(1-epsilon)). Stocke tout en attributs sur l'instance.

        À appeler uniquement après un solve() réussi sur self.RC."""
        RC = self.RC
        arch = self.params.get('RC_ARCH', 'REC')

        # reset
        self.NTU_gh = self.NTU_rec = self.NTU_rec_LT = self.NTU_rec_HT = None
        self.NTU_cd = self.NTU_weighted_sum = None
        self.eps_gh = self.eps_rec = self.eps_rec_LT = self.eps_rec_HT = None
        self.eps_cd = None

        try:
            Q_gh    = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh  = RC.components['GasHeater'].model.epsilon
            ntu_gh  = _epsilon_to_ntu(eta_gh)

            Q_cond   = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon
            ntu_cd   = _epsilon_to_ntu(eta_cond)

            self.NTU_gh = ntu_gh
            self.NTU_cd = ntu_cd
            self.eps_gh = eta_gh
            self.eps_cd = eta_cond

            if arch == 'REC':
                Q_rec   = RC.components['Recuperator'].model.Q.Q_dot
                eta_rec = RC.components['Recuperator'].model.epsilon
                ntu_rec = _epsilon_to_ntu(eta_rec)
                self.NTU_rec = ntu_rec
                self.eps_rec = eta_rec
                denom = Q_gh + Q_rec + Q_cond
                self.NTU_weighted_sum = (Q_gh*ntu_gh + Q_rec*ntu_rec + Q_cond*ntu_cd) / denom

            elif arch == 'basic':
                denom = Q_gh + Q_cond
                self.NTU_weighted_sum = (Q_gh*ntu_gh + Q_cond*ntu_cd) / denom

            elif arch == 'Recomp':
                Q_rec_LT   = RC.components['RecupLT'].model.Q.Q_dot
                eta_rec_LT = RC.components['RecupLT'].model.epsilon
                ntu_rec_LT = _epsilon_to_ntu(eta_rec_LT)

                Q_rec_HT   = RC.components['RecupHT'].model.Q.Q_dot
                eta_rec_HT = RC.components['RecupHT'].model.epsilon
                ntu_rec_HT = _epsilon_to_ntu(eta_rec_HT)

                self.NTU_rec_LT = ntu_rec_LT
                self.NTU_rec_HT = ntu_rec_HT
                self.eps_rec_LT = eta_rec_LT
                self.eps_rec_HT = eta_rec_HT
                denom = Q_gh + Q_rec_LT + Q_rec_HT + Q_cond
                self.NTU_weighted_sum = (
                    Q_gh*ntu_gh + Q_rec_LT*ntu_rec_LT + Q_rec_HT*ntu_rec_HT + Q_cond*ntu_cd
                ) / denom

            elif arch == 'Recomp_1_recup':
                Q_rec_LT   = RC.components['RecupLT'].model.Q.Q_dot
                eta_rec_LT = RC.components['RecupLT'].model.epsilon
                ntu_rec_LT = _epsilon_to_ntu(eta_rec_LT)

                self.NTU_rec_LT = ntu_rec_LT
                self.eps_rec_LT = eta_rec_LT
                denom = Q_gh + Q_rec_LT + Q_cond
                self.NTU_weighted_sum = (Q_gh*ntu_gh + Q_rec_LT*ntu_rec_LT + Q_cond*ntu_cd) / denom

        except Exception as e:
            self._log_penalty(f"Calcul NTU final impossible : {e}")
            self.NTU_gh = self.NTU_rec = self.NTU_rec_LT = self.NTU_rec_HT = None
            self.NTU_cd = self.NTU_weighted_sum = None
            self.eps_gh = self.eps_rec = self.eps_rec_LT = self.eps_rec_HT = None
            self.eps_cd = None

    # ------------------------------------------------------------------ final eval

    def _evaluate_final(self, best_pos):
        """Re-evaluates the best PSO position with full diagnostics."""
        unpacked = self._unpack_position(best_pos)
        self.it_var.update(unpacked)

        self._HSource_props['m_dot'] = unpacked['mdot_HS']
        self._CSource_props['m_dot'] = unpacked['mdot_CS']

        # --- expose les variables de décision en attributs plats ---
        self.P_high         = unpacked.get('P_high')
        self.m_dot           = unpacked.get('mdot')
        self.mdot            = unpacked.get('mdot')
        self.m_dot_HS        = unpacked.get('mdot_HS')
        self.m_dot_CS        = unpacked.get('mdot_CS')
        self.m_dot_HS_factor = unpacked.get('m_dot_HS_factor')
        self.m_dot_CS_factor = unpacked.get('m_dot_CS_factor')
        self.spliter_frac    = unpacked.get('spliter_frac')
        self.eta_gh          = unpacked.get('eta_gh')
        self.PP_gh           = unpacked.get('PP_gh')
        self.eta_rec         = unpacked.get('eta_rec')
        self.eta_rec_LT      = unpacked.get('eta_rec_LT')
        self.eta_rec_HT      = unpacked.get('eta_rec_HT')
        self.PP_cd           = unpacked.get('PP_cd')

        self.set_RC()
        RC = self.RC

        try:
            RC.solve()
        except Exception as e:
            self._log_penalty(f"Final solve exception: {e}")
            self.eta = self.W_dot_net = self.W_net = None
            return

        if not getattr(RC, 'converged', True):
            self._log_penalty("Final solve did not converge")
            self.eta = self.W_dot_net = self.W_net = None
            return

        try:
            T_exp_ex  = RC.components['Expander'].model.ex.T
            P_exp_ex  = RC.components['Expander'].model.ex.p
            T_sat_exp = PropsSI('T', 'P', P_exp_ex, 'Q', 1, 'CO2')
            SH_exp    = T_exp_ex - T_sat_exp
        except Exception:
            SH_exp = 50.0

        if SH_exp < 0:
            self._log_penalty(f"Drops in expansion (SH = {SH_exp:.1f} K)")
            self.eta = self.W_dot_net = self.W_net = None
            return

        arch = self.params.get('RC_ARCH', 'REC')
        if arch in ("Recomp", "Recomp_1_recup"):
            W_cp = RC.components['Compressor'].model.W.W_dot
        else:
            W_cp = 0

        W_exp  = RC.components['Expander'].model.W.W_dot
        W_pump = RC.components['Pump'].model.W.W_dot
        Q_gh   = RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = RC.components['GasHeater'].model.su_H.D
        m_HS_act   = RC.components['GasHeater'].model.su_H.m_dot

        rho_CS     = RC.components['Condenser'].model.su_C.D
        m_CS_act   = RC.components['Condenser'].model.su_C.m_dot

        self.W_pump_aux_HS = W_pump_aux_HS = self.params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * self.params.get('eta_pp_aux', 0.8))

        self.W_pump_aux_CS = W_pump_aux_CS = self.params.get('DP_c_cond', 0.5e5) * m_CS_act / \
                     (rho_CS * self.params.get('eta_pp_aux', 0.8))

        self.W_dot_net = W_exp - W_pump - W_pump_aux_HS - W_pump_aux_CS - W_cp
        self.W_net     = self.W_dot_net  # alias, pour extract_run_data
        self.eta       = self.W_dot_net / Q_gh if Q_gh > 0 else 0.0

        # --- calcul du NTU final (corrige le bug : c'était toujours None) ---
        self._compute_final_ntu()

    # ------------------------------------------------------------------ optimise

    def opt_RC(self, n_jobs=1, n_particles=100, max_iter=30, patience=None,
               init_pos=None, warm_spread=0.05, warm_fraction=0.5, ntop=None):
        """
        PSO optimisation with parallel particle evaluation via joblib.
        Si `ntop` est fourni, self.top_positions est peuplé avec les `ntop`
        meilleures positions uniques (utile pour un dimensionnement ultérieur).
        """
        if patience is None:
            patience = max(1, max_iter // 5)

        arch = self.params.get('RC_ARCH', 'REC')

        if arch == "Recomp":
            eta_gh_disc     = self.params['eta_gh_disc']
            PP_gh_disc      = self.params['PP_gh_disc']
            eta_rec_disc    = self.params['eta_rec_disc']
            eta_rec_HT_disc = self.params['eta_rec_HT_disc']
            PP_cd_disc      = self.params['PP_cd_disc']

            lb = np.array([
                    self.params['P_high_bounds'][0], self.params['m_dot_bounds'][0],
                    self.params['m_dot_HS_fact_bounds'][0], self.params['spliter_frac_bounds'][0],
                    eta_gh_disc[0], PP_gh_disc[0], eta_rec_disc[0], eta_rec_HT_disc[0],
                    PP_cd_disc[0], self.params['m_dot_CS_fact_bounds'][0]
                ])
            ub = np.array([
                    self.params['P_high_bounds'][1], self.params['m_dot_bounds'][1],
                    self.params['m_dot_HS_fact_bounds'][1], self.params['spliter_frac_bounds'][1],
                    eta_gh_disc[-1], PP_gh_disc[-1], eta_rec_disc[-1], eta_rec_HT_disc[-1],
                    PP_cd_disc[-1], self.params['m_dot_CS_fact_bounds'][1]
                ])
            discrete_vars = {4: eta_gh_disc, 5: PP_gh_disc, 6: eta_rec_disc,
                              7: eta_rec_HT_disc, 8: PP_cd_disc}

        elif arch == "Recomp_1_recup":
            eta_gh_disc   = self.params['eta_gh_disc']
            PP_gh_disc    = self.params['PP_gh_disc']
            eta_rec_disc  = self.params['eta_rec_disc']
            PP_cd_disc    = self.params['PP_cd_disc']

            lb = np.array([
                    self.params['P_high_bounds'][0], self.params['m_dot_bounds'][0],
                    self.params['m_dot_HS_fact_bounds'][0], self.params['spliter_frac_bounds'][0],
                    eta_gh_disc[0], PP_gh_disc[0], eta_rec_disc[0], PP_cd_disc[0],
                    self.params['m_dot_CS_fact_bounds'][0]
                ])
            ub = np.array([
                    self.params['P_high_bounds'][1], self.params['m_dot_bounds'][1],
                    self.params['m_dot_HS_fact_bounds'][1], self.params['spliter_frac_bounds'][1],
                    eta_gh_disc[-1], PP_gh_disc[-1], eta_rec_disc[-1], PP_cd_disc[-1],
                    self.params['m_dot_CS_fact_bounds'][1]
                ])
            discrete_vars = {4: eta_gh_disc, 5: PP_gh_disc, 6: eta_rec_disc, 7: PP_cd_disc}

        elif arch == "REC":
            eta_gh_disc   = self.params['eta_gh_disc']
            PP_gh_disc    = self.params['PP_gh_disc']
            eta_rec_disc  = self.params['eta_rec_disc']
            PP_cd_disc    = self.params['PP_cd_disc']

            lb = np.array([
                    self.params['P_high_bounds'][0], self.params['m_dot_bounds'][0],
                    self.params['m_dot_HS_fact_bounds'][0], eta_gh_disc[0], PP_gh_disc[0],
                    eta_rec_disc[0], PP_cd_disc[0], self.params['m_dot_CS_fact_bounds'][0]
                ])
            ub = np.array([
                    self.params['P_high_bounds'][1], self.params['m_dot_bounds'][1],
                    self.params['m_dot_HS_fact_bounds'][1], eta_gh_disc[-1], PP_gh_disc[-1],
                    eta_rec_disc[-1], PP_cd_disc[-1], self.params['m_dot_CS_fact_bounds'][1]
                ])
            discrete_vars = {3: eta_gh_disc, 4: PP_gh_disc, 5: eta_rec_disc, 6: PP_cd_disc}

        elif arch == "basic":
            eta_gh_disc   = self.params['eta_gh_disc']
            PP_gh_disc    = self.params['PP_gh_disc']
            PP_cd_disc    = self.params['PP_cd_disc']

            lb = np.array([
                    self.params['P_high_bounds'][0], self.params['m_dot_bounds'][0],
                    self.params['m_dot_HS_fact_bounds'][0], eta_gh_disc[0], PP_gh_disc[0],
                    PP_cd_disc[0], self.params['m_dot_CS_fact_bounds'][0]
                ])
            ub = np.array([
                    self.params['P_high_bounds'][1], self.params['m_dot_bounds'][1],
                    self.params['m_dot_HS_fact_bounds'][1], eta_gh_disc[-1], PP_gh_disc[-1],
                    PP_cd_disc[-1], self.params['m_dot_CS_fact_bounds'][1]
                ])
            discrete_vars = {3: eta_gh_disc, 4: PP_gh_disc, 5: PP_cd_disc}

        else:
            raise ValueError()

        # --- warm start ---
        pso_init_pos = None
        if init_pos is not None:
            seed = np.asarray(init_pos, dtype=float)
            if seed.ndim == 1:
                seed    = np.clip(seed, lb, ub)
                n_warm  = max(1, int(round(warm_fraction * n_particles)))
                n_rand  = n_particles - n_warm
                noise   = np.random.uniform(-warm_spread, warm_spread,
                                            size=(n_warm, len(lb)))
                warm    = np.clip(seed[None, :] * (1.0 + noise), lb, ub)
                rand    = np.random.uniform(lb, ub, size=(n_rand, len(lb)))
                pso_init_pos = np.vstack([warm, rand])
                print(f"  → Warm start: {n_warm}/{n_particles} particles around "
                      f"P={seed[0]/1e5:.1f} bar, ṁ={seed[1]:.1f}")
            else:
                pso_init_pos = np.clip(seed, lb, ub)

        # Constantes indépendantes de la particule évaluée : calculées une
        # seule fois ici plutôt qu'à chaque appel de system_RC_parallel
        # (potentiellement des dizaines de milliers d'appels PropsSI économisés
        # sur une campagne complète).
        _P_sat_CS_const = PropsSI('P', 'T', self._CSource_props['T'], 'Q', 0.5, self.fluid)
        _P_crit_const   = PropsSI('PCRIT', self.fluid)

        input_data = {
            'fluid'   : self.fluid,
            'params'  : self.params,
            'obj'     : self.obj,
            'HSource' : {
                'T'     : self._HSource_props['T'],
                'P'     : self._HSource_props['P'],
                'fluid' : self._HSource_props['fluid'],
            },
            'CSource' : {
                'T'     : self._CSource_props['T'],
                'P'     : self._CSource_props['P'],
                'fluid' : self._CSource_props['fluid'],
            },
            'RC_ARCH' : arch,
            'discrete_vars' : discrete_vars,
            'P_sat_CS' : _P_sat_CS_const,
            'P_crit'   : _P_crit_const,
        }

        def discretize(x):
            x = np.array(x, dtype=float)
            for idx, allowed_vals in discrete_vars.items():
                allowed_vals = np.array(allowed_vals, dtype=float)
                x[idx] = allowed_vals[np.argmin(np.abs(allowed_vals - x[idx]))]
            return x

        def objective_wrapper(X):
            results = np.array(
                Parallel(n_jobs=n_jobs, backend='loky')(
                    delayed(system_RC_parallel)(x, input_data) for x in X
                )
            )
            costs     = results[:, 0]
            penalties = results[:, 1]

            for x_i, pen_i, cost_i in zip(X, penalties, costs):
                if pen_i == 0 and np.isfinite(cost_i):
                    x_disc = discretize(x_i)
                    self.allowable_positions.append({
                        'x'     : x_disc.copy(),
                        'score' : float(cost_i),
                    })

            return costs

        optimizer = GlobalBestPSO(
            n_particles = n_particles,
            dimensions  = len(ub),
            options     = {'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds      = (lb, ub),
            init_pos    = pso_init_pos,
        )

        best_cost  = np.inf
        no_improve = 0

        pbar = tqdm(range(max_iter), desc="PSO Optimizing", ncols=80)
        for i in pbar:
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current = optimizer.swarm.best_cost

            if current < best_cost - 1e-3:
                best_cost  = current
                no_improve = 0
            else:
                no_improve += 1

            pbar.set_postfix(best_cost=f"{best_cost:.6f}")

            if no_improve >= patience:
                pbar.set_description("Stopped (no improvement)")
                break
        pbar.close()

        self._evaluate_final(optimizer.swarm.best_pos)

        bp = optimizer.swarm.best_pos
        print("\n" + "="*55)
        print("  OPTIMAL RESULT")
        print("="*55)
        print(f"  P_high            : {bp[0]/1e5:.2f}  bar")
        print(f"  m_dot (CO2)       : {bp[1]:.4f}  kg/s")
        if self.W_dot_net is not None:
            print(f"  W_net             : {self.W_dot_net/1e3:.3f}  kW")
            print(f"  Thermal η         : {self.eta*100:.3f}  %")
        else:
            print("  ⚠️  Final solve failed — see penalty log.")
        print("="*55)

        # --- variables de décision brutes (pour affiner les bornes ensuite) ---
        if self.W_dot_net is not None:
            print("\n" + "-"*55)
            print("  VARIABLES DE DÉCISION (valeurs discrétisées PSO)")
            print("-"*55)
            print(f"  m_dot_HS          : {self.m_dot_HS:.3f}  kg/s"
                  + (f"   (facteur = {self.m_dot_HS_factor:.3f})" if self.m_dot_HS_factor is not None else ""))
            print(f"  m_dot_CS          : {self.m_dot_CS:.3f}  kg/s"
                  + (f"   (facteur = {self.m_dot_CS_factor:.3f})" if self.m_dot_CS_factor is not None else ""))
            print(f"  eta_gh (target)   : {self.eta_gh:.3f}")
            print(f"  PP_gh             : {self.PP_gh:.3f}  K")
            if arch == 'REC':
                print(f"  eta_rec (target)  : {self.eta_rec:.3f}")
            elif arch == 'Recomp':
                print(f"  eta_rec_LT (target): {self.eta_rec_LT:.3f}")
                print(f"  eta_rec_HT (target): {self.eta_rec_HT:.3f}")
                print(f"  spliter_frac      : {self.spliter_frac:.3f}")
            elif arch == 'Recomp_1_recup':
                print(f"  eta_rec (target)  : {self.eta_rec:.3f}")
                print(f"  spliter_frac      : {self.spliter_frac:.3f}")
            print(f"  PP_cd             : {self.PP_cd:.3f}  K")
            print("-"*55)

        # --- NTU / effectivité par échangeur ---
        if self.W_dot_net is not None and self.NTU_weighted_sum is not None:
            print("\n" + "-"*55)
            print("  NTU & EFFICACITÉ PAR ÉCHANGEUR")
            print("-"*55)
            print(f"  Gas Heater        : NTU = {self.NTU_gh:6.3f}   |   ε = {self.eps_gh*100:6.2f} %")

            if arch == 'REC':
                print(f"  Recuperator       : NTU = {self.NTU_rec:6.3f}   |   ε = {self.eps_rec*100:6.2f} %")
            elif arch == 'Recomp':
                print(f"  Recuperator LT    : NTU = {self.NTU_rec_LT:6.3f}   |   ε = {self.eps_rec_LT*100:6.2f} %")
                print(f"  Recuperator HT    : NTU = {self.NTU_rec_HT:6.3f}   |   ε = {self.eps_rec_HT*100:6.2f} %")
            elif arch == 'Recomp_1_recup':
                print(f"  Recuperator LT    : NTU = {self.NTU_rec_LT:6.3f}   |   ε = {self.eps_rec_LT*100:6.2f} %")
            # 'basic' n'a pas de récupérateur

            print(f"  Condenser         : NTU = {self.NTU_cd:6.3f}   |   ε = {self.eps_cd*100:6.2f} %")
            print("-"*55)
            print(f"  NTU pondéré total : {self.NTU_weighted_sum:.4f}")
            print("-"*55)

        total = sum(self.penalty_log.values())
        if total:
            print("\n" + "="*55)
            print("  PENALTY SUMMARY")
            print("="*55)
            for reason, count in sorted(self.penalty_log.items(),
                                        key=lambda kv: kv[1], reverse=True):
                print(f"  [{count:4d} | {count/total*100:5.1f}%] : {reason}")
            print("="*55)
        self.penalty_log = {}

        # --- top positions (utilisées par un éventuel dimensionnement en aval) ---
        if ntop is not None:
            unique_positions = {}
            for entry in self.allowable_positions:
                x = np.array(entry['x'], dtype=float)
                score = float(entry['score'])
                key = tuple(np.round(x, 8))
                if key not in unique_positions or score < unique_positions[key]['score']:
                    unique_positions[key] = {'x': x, 'score': score}

            unique_list = list(unique_positions.values())
            unique_list.sort(key=lambda e: e['score'])
            self.top_positions = unique_list[:ntop]

        return optimizer

#%% Main

if __name__ == "__main__":

    # -*- coding: utf-8 -*-
    """
    Optimisation de la répartition des NTU des échangeurs (CO2RC_HX_optimizer),
    à efficacité et puissance nette fixées, pour l'architecture REC.

    Pour chaque condition (architecture, T_H, eta_obj), on relance
    l'optimisation PSO N_RUNS fois et on conserve les BEST_KEEP runs
    présentant la plus faible valeur de l'objectif (somme pondérée des NTU) :
      - si moins de BEST_KEEP runs valides ont été trouvés, on les ajoute tous ;
      - une fois BEST_KEEP atteints, un nouveau run ne remplace le pire de la
        liste (NTU pondéré le plus élevé) que s'il fait mieux (NTU plus faible).

    Les runs sont enregistrés (et réécrits à chaque mise à jour, pour être
    robuste à une interruption) dans un fichier CSV.
    """

    import csv
    import multiprocessing
    import numpy as np
    import matplotlib.pyplot as plt

    # ---------------------------------------------------------------------
    # Gestion des "N meilleurs runs" par condition (architecture, T_H, eta_obj)
    # ---------------------------------------------------------------------
    BEST_KEEP = 5  # nombre de meilleurs runs conservés par condition
    CSV_PATH = "best_runs_HX.csv"

    FIELDNAMES = [
        "architecture", "T_H_C", "eta_obj", "rank",
        "P_high_Pa", "m_dot_CO2_kg_s", "m_dot_HS_factor", "m_dot_HS_kg_s",
        "m_dot_CS_factor", "m_dot_CS_kg_s", "spliter_frac",
        "eta_gh", "PP_gh", "eta_rec", "eta_rec_HT", "PP_cd",
        "NTU_gh", "NTU_rec", "NTU_cd", "NTU_weighted_sum",
        "eps_gh_actual", "eps_rec_actual", "eps_rec_HT_actual", "eps_cd_actual",
        "W_net_W", "eta_actual",
    ]

    # best_runs[(arch, T_H_C, eta_obj)] = liste des BEST_KEEP meilleurs runs
    # (dicts), triée par NTU_weighted_sum croissant (plus petit = meilleur)
    best_runs = {}

    # ---------------------------------------------------------------------
    # Seuil sous lequel le récupérateur est considéré comme "bypassé"
    # (effectivité réelle quasi nulle -> le cycle REC/Recomp dégénère en un
    # cycle 'basic' déguisé et doit être disqualifié, cf. discussion NTU).
    # ---------------------------------------------------------------------
    EPS_REC_MIN = 0.3

    def is_recuperator_valid(Optimizer, arch, eps_min=EPS_REC_MIN):
        """Retourne False si le(s) récupérateur(s) sont quasi bypassés
        (effectivité RÉELLE, post-solve, proche de 0), auquel cas la
        solution REC/Recomp dégénère en un cycle 'basic' déguisé."""
        if arch == 'basic':
            return True  # pas de récupérateur, rien à vérifier

        if arch == 'REC':
            eps = getattr(Optimizer, 'eps_rec', None)
            return eps is not None and eps >= eps_min

        if arch == 'Recomp_1_recup':
            eps = getattr(Optimizer, 'eps_rec_LT', None)
            return eps is not None and eps >= eps_min

        if arch == 'Recomp':
            eps_lt = getattr(Optimizer, 'eps_rec_LT', None)
            eps_ht = getattr(Optimizer, 'eps_rec_HT', None)
            return (eps_lt is not None and eps_lt >= eps_min and
                    eps_ht is not None and eps_ht >= eps_min)

        return True

    def extract_run_data(Optimizer, arch):
        """Récupère les grandeurs d'intérêt sur l'objet Optimizer après
        résolution.

        CORRIGE : les attributs lus ici existent désormais réellement sur
        l'instance (voir _evaluate_final / _compute_final_ntu dans la
        classe CO2RC_HX_optimizer). Avant cette correction, tous les
        NTU_* et P_high/m_dot/etc. étaient None car jamais calculés ou
        stockés sous ces noms — ce qui empêchait tout run d'être conservé
        par update_best_runs()."""

        return {
            "P_high_Pa"       : getattr(Optimizer, "P_high", None),
            "m_dot_CO2_kg_s"  : getattr(Optimizer, "m_dot", None),
            "m_dot_HS_factor" : getattr(Optimizer, "m_dot_HS_factor", None),
            "m_dot_HS_kg_s"   : getattr(Optimizer, "m_dot_HS", None),
            "m_dot_CS_factor" : getattr(Optimizer, "m_dot_CS_factor", None),
            "m_dot_CS_kg_s"   : getattr(Optimizer, "m_dot_CS", None),
            "spliter_frac"    : getattr(Optimizer, "spliter_frac", None) if arch in ("Recomp", "Recomp_1_recup") else None,
            "eta_gh"          : getattr(Optimizer, "eta_gh", None),
            "PP_gh"           : getattr(Optimizer, "PP_gh", None),
            "eta_rec"         : getattr(Optimizer, "eta_rec", None) if arch != "Recomp" else getattr(Optimizer, "eta_rec_LT", None),
            "eta_rec_HT"      : getattr(Optimizer, "eta_rec_HT", None) if arch == "Recomp" else None,
            "PP_cd"           : getattr(Optimizer, "PP_cd", None),
            "NTU_gh"          : getattr(Optimizer, "NTU_gh", None),
            "NTU_rec"         : getattr(Optimizer, "NTU_rec", None) if arch != "Recomp" else getattr(Optimizer, "NTU_rec_LT", None),
            "NTU_cd"          : getattr(Optimizer, "NTU_cd", None),
            "NTU_weighted_sum": getattr(Optimizer, "NTU_weighted_sum", None),
            "eps_gh_actual"    : getattr(Optimizer, "eps_gh", None),
            "eps_rec_actual"   : getattr(Optimizer, "eps_rec", None) if arch != "Recomp" else getattr(Optimizer, "eps_rec_LT", None),
            "eps_rec_HT_actual": getattr(Optimizer, "eps_rec_HT", None) if arch == "Recomp" else None,
            "eps_cd_actual"    : getattr(Optimizer, "eps_cd", None),
            "W_net_W"         : getattr(Optimizer, "W_net", None),
            "eta_actual"      : getattr(Optimizer, "eta", None),
        }

    def update_best_runs(key, new_run, max_keep=BEST_KEEP):
        """Ajoute new_run à la liste des meilleurs runs pour `key`. On minimise
        NTU_weighted_sum : si la liste dépasse max_keep après ajout, le run le
        plus mauvais (NTU pondéré le plus élevé) est écarté -- ce qui peut être
        new_run lui-même s'il n'est pas assez bon. Les runs sans valeur
        d'objectif exploitable (None) sont ignorés."""
        if new_run["NTU_weighted_sum"] is None:
            print("    [!] NTU_weighted_sum indisponible pour ce run -> non enregistré (solve final probablement échoué)")
            return

        lst = best_runs.setdefault(key, [])
        lst.append(new_run)
        lst.sort(key=lambda r: r["NTU_weighted_sum"])  # croissant : le meilleur (plus petit) en premier
        del lst[max_keep:]

    def write_best_runs_csv(path=CSV_PATH):
        rows = []
        for (arch, T_H_C, eta_obj), lst in best_runs.items():
            for rank, run in enumerate(lst, start=1):
                row = {"architecture": arch, "T_H_C": T_H_C, "eta_obj": eta_obj, "rank": rank}
                row.update(run)
                rows.append(row)
        rows.sort(key=lambda r: (r["architecture"], r["T_H_C"], r["eta_obj"], r["rank"]))

        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    def load_best_runs_csv(path=CSV_PATH):
        """Relit le CSV et renvoie les lignes de rang 1 (meilleur run par
        condition), pratique pour tracer sans avoir à relancer les
        optimisations (ex. après une interruption)."""
        rows = []
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for r in reader:
                if r["rank"] == "1":
                    rows.append(r)
        return rows

    def plot_results(path=CSV_PATH):
        """Trace :
          1) le NTU pondéré (meilleur run) vs T_H, une courbe par eta_obj ;
          2) la répartition du NTU par composant (gh / rec / cd) pour le
             meilleur run de chaque condition (barres empilées).
        """
        rows = load_best_runs_csv(path)
        if not rows:
            print("[!] Aucune donnée de rang 1 trouvée dans le CSV -> pas de graphique.")
            return

        archs = sorted(set(r["architecture"] for r in rows))
        eta_objs = sorted(set(float(r["eta_obj"]) for r in rows))

        # ---- 1) NTU pondéré (meilleur) vs T_H, une courbe par eta_obj ----
        for arch in archs:
            plt.figure(figsize=(8, 5))
            for eta_obj in eta_objs:
                sub = [r for r in rows if r["architecture"] == arch and float(r["eta_obj"]) == eta_obj]
                sub = sorted(sub, key=lambda r: float(r["T_H_C"]))
                T_vals = [float(r["T_H_C"]) for r in sub]
                ntu_vals = [
                    float(r["NTU_weighted_sum"]) if r["NTU_weighted_sum"] not in ("", "None") else np.nan
                    for r in sub
                ]
                plt.plot(T_vals, ntu_vals, marker='o', linewidth=2, label=f"eta_obj = {eta_obj:.2f}")

            plt.title(f"NTU pondéré minimal vs T_H -- architecture {arch}")
            plt.xlabel("Température source chaude [°C]")
            plt.ylabel("NTU pondéré (meilleur run)")
            plt.grid(True)
            plt.legend(title="Efficacité cible")
            plt.tight_layout()
            plt.show()

        # ---- 2) Répartition du NTU par composant (barres empilées) ----
        for arch in archs:
            sub = [r for r in rows if r["architecture"] == arch]
            sub = sorted(sub, key=lambda r: (float(r["T_H_C"]), float(r["eta_obj"])))

            labels = [f"{float(r['T_H_C']):.0f}°C\nη={float(r['eta_obj']):.2f}" for r in sub]

            def _val(r, key):
                v = r.get(key)
                return float(v) if v not in (None, "", "None") else 0.0

            ntu_gh = [_val(r, "NTU_gh") for r in sub]
            ntu_rec = [_val(r, "NTU_rec") for r in sub]
            ntu_cd = [_val(r, "NTU_cd") for r in sub]

            x = np.arange(len(sub))

            plt.figure(figsize=(max(6, len(sub) * 1.2), 5))
            plt.bar(x, ntu_gh, label="NTU_gh (gas heater)")
            plt.bar(x, ntu_rec, bottom=ntu_gh, label="NTU_rec (récupérateur)")
            bottom_2 = [a + b for a, b in zip(ntu_gh, ntu_rec)]
            plt.bar(x, ntu_cd, bottom=bottom_2, label="NTU_cd (condenseur)")

            plt.xticks(x, labels)
            plt.title(f"Répartition du NTU par échangeur -- architecture {arch}")
            plt.ylabel("NTU")
            plt.legend()
            plt.tight_layout()
            plt.show()

    if __name__ == "__main__":

        n_cores = multiprocessing.cpu_count()

        # ---- sweep ----
        # T_vec = np.linspace(150, 350, 5) + 273.15  # 150, 200, 250, 300, 350 °C
        T_vec = np.linspace(150, 150, 1) + 273.15  # 150, 200, 250, 300, 350 °C

        n_MW = 1  # W
        W_dot_obj = n_MW * 1e6  # W

        # Niveaux d'efficacité cible par T_H (issus de campagnes précédentes :
        # le max trouvé à chaque T_H, et deux valeurs légèrement dégradées).
        # -> à ajuster si T_vec est modifié.
        ETA_OBJ_BY_T_H_C = {
            150.0: [0.13, 0.12, 0.11],
            200.0: [0.17, 0.16, 0.15],
            250.0: [0.21, 0.20, 0.19],
            300.0: [0.24, 0.23, 0.22],
            350.0: [0.26, 0.25, 0.24],
        }

        ARCH_LIST = ['REC']  # seule architecture demandée pour cette campagne
        N_RUNS = 1  # nombre d'optimisations par condition (arch, T_H, eta_obj)

        # ---------------------------------------------------------------
        # Bornes ADAPTÉES À LA TEMPÉRATURE, d'après l'analyse de 45 runs
        # REC (T_H = 150→350°C, eta_obj = 0.11→0.26) :
        #
        #   T_H [°C] : 150          200          250          300          350
        #   P_high   : [133,198]    [147,196]    [154,199]    [154,195]    [138,197]  bar
        #   m_dot    : [325,427]    [234,295]    [171,212]    [141,189]    [117,156]  kg/s
        #
        # -> P_high reste dans la même plage quel que soit T_H (talonne
        #    systématiquement ~195-200 bar) : bornes FIXES, pas de
        #    dépendance en T retenue (élargies par rapport à l'observé,
        #    cf. campagne précédente).
        # -> m_dot chute d'un facteur ~3.6x entre 150°C et 350°C (meilleur
        #    rendement -> moins de débit nécessaire pour la même puissance
        #    visée) : bornes INTERPOLÉES linéairement en T_H (à plat hors
        #    de la plage [150,350]°C), avec une marge de ±15 % autour des
        #    valeurs observées.
        # ---------------------------------------------------------------
        m_dot_HS_fact_bounds = [0.1, 3]
        m_dot_CS_fact_bounds = [5, 15]
        P_high_bounds = np.array([110, 220]) * 1e5  # fixe (cf. justification ci-dessus)
        spliter_frac_bounds = np.array([0.01, 0.99])

        _T_ANCHORS_C = np.array([150.0, 200.0, 250.0, 300.0, 350.0])
        # bornes observées (min/max de m_dot par T_H, cf. tableau ci-dessus),
        # exprimées en multiple de n_MW (=10 dans les runs analysés), avec
        # marge -15 % / +15 % pour laisser de la place à l'exploration PSO.
        _M_DOT_MIN_MULT = np.array([325, 234, 171, 141, 117]) / (10*n_MW) * 0.85
        _M_DOT_MAX_MULT = np.array([427, 295, 212, 189, 156]) / (10*n_MW) * 1.15

        def m_dot_bounds_for_T(T_H_K, n_MW_local=n_MW):
            """Bornes de m_dot (CO2, kg/s) interpolées linéairement en T_H
            (extrapolation à plat hors de [150, 350]°C) depuis les runs
            observés. Retourne un np.array([lb, ub]) en kg/s."""
            T_C = T_H_K - 273.15
            lb_mult = np.interp(T_C, _T_ANCHORS_C, _M_DOT_MIN_MULT)
            ub_mult = np.interp(T_C, _T_ANCHORS_C, _M_DOT_MAX_MULT)
            return np.array([lb_mult, ub_mult]) * n_MW_local

        # Discrete Variable choices
        eta_gh_disc = np.arange(0.80, 1.00, 0.02)
        PP_gh_disc = np.arange(1, 10, 1)
        eta_rec_disc = np.arange(0.60, 0.98, 0.02)
        eta_rec_HT_disc = np.arange(0.60, 0.98, 0.02)
        PP_cd_disc = np.arange(1, 10, 1)

        for arch in ARCH_LIST:
            print(f"\n{'='*60}")
            print(f"  Architecture : {arch}")
            print(f"{'='*60}")

            for T in T_vec:
                T_H_C = round(T - 273.15, 1)

                # --- bornes m_dot adaptées à ce T_H ---
                m_dot_bounds = m_dot_bounds_for_T(T)
                print(f"\n[Bornes @ T_H={T_H_C:.1f}°C] "
                      f"m_dot ∈ [{m_dot_bounds[0]:.1f}, {m_dot_bounds[1]:.1f}] kg/s   "
                      f"P_high ∈ [{P_high_bounds[0]/1e5:.0f}, {P_high_bounds[1]/1e5:.0f}] bar")

                eta_obj_list = ETA_OBJ_BY_T_H_C.get(T_H_C, [0.15])  # défaut si T_H hors table

                for eta_obj in eta_obj_list:
                    print(f"\n--- T_H = {T_H_C:.1f} °C | eta_obj = {eta_obj:.3f} ---")

                    condition_key = (arch, T_H_C, eta_obj)

                    for run_idx in range(N_RUNS):
                        print(f"  Run {run_idx+1}/{N_RUNS}...")

                        try:
                            Optimizer = CO2RC_Cost_optimizer('CO2')

                            Optimizer.set_parameters(
                                RC_ARCH=arch,  # 'basic', 'REC', 'Recomp_1_recup', 'Recomp'

                                # Pump
                                eta_pp=0.85,
                                eta_pp_aux=0.8,

                                # Compressor (for recompression layouts)
                                eta_cp=0.8,

                                # GasHeater
                                DP_h_gh=50*1e3,
                                DP_c_gh=50*1e3,

                                # Recuperator
                                PP_rec=0,
                                DP_h_rec=50*1e3,
                                DP_c_rec=50*1e3,

                                # Expander
                                eta_exp=0.94,

                                # Condenser
                                SC_cd=0.1,
                                DP_h_cond=50*1e3,
                                DP_c_cond=50*1e3,

                                # Bounds
                                P_high_bounds=P_high_bounds,
                                m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
                                m_dot_CS_fact_bounds=m_dot_CS_fact_bounds,
                                m_dot_bounds=m_dot_bounds,
                                spliter_frac_bounds=spliter_frac_bounds,

                                # Discrete Values
                                eta_gh_disc=eta_gh_disc,
                                PP_gh_disc=PP_gh_disc,
                                eta_rec_disc=eta_rec_disc,
                                eta_rec_HT_disc=eta_rec_HT_disc,
                                PP_cd_disc=PP_cd_disc,
                            )

                            if Optimizer.params['RC_ARCH'] == "Recomp":
                                Optimizer.set_it_var(P_high=140e5, mdot=20.0 * n_MW, mdot_HS=15.0 * n_MW,
                                                      spliter_frac=0.9, eta_gh=0.95, PP_gh=5,
                                                      eta_rec_LT=0.8, eta_rec_HT=0.8, PP_cd=5, mdot_CS=200 * n_MW)
                            elif Optimizer.params['RC_ARCH'] == "Recomp_1_recup":
                                Optimizer.set_it_var(P_high=100e5, mdot=20.0 * n_MW, mdot_HS=15.0 * n_MW,
                                                      spliter_frac=1, eta_gh=0.95, PP_gh=5,
                                                      eta_rec=0.8, PP_cd=5, mdot_CS=200 * n_MW)
                            elif Optimizer.params['RC_ARCH'] == "REC":
                                Optimizer.set_it_var(P_high=100e5, mdot=20.0 * n_MW, mdot_HS=15.0 * n_MW,
                                                      eta_gh=0.95, PP_gh=5, eta_rec=0.8, PP_cd=5, mdot_CS=200 * n_MW)
                            elif Optimizer.params['RC_ARCH'] == "basic":
                                Optimizer.set_it_var(P_high=100e5, mdot=20.0 * n_MW, mdot_HS=15.0 * n_MW,
                                                      eta_gh=0.95, PP_gh=5, PP_cd=5, mdot_CS=200 * n_MW)

                            Optimizer.set_obj(W_dot=W_dot_obj, eta=eta_obj)

                            Optimizer.set_CSource(T=15 + 273.15, P=5e5, fluid='Water', m_dot=1000.0)
                            Optimizer.set_HSource(T=T, P=10e5, fluid='INCOMP::TVP1', m_dot=50.0)

                            Optimizer.set_RC()
                            Optimizer.opt_RC(n_jobs=n_cores - 1, n_particles=50, max_iter=50, patience=20)

                            if not is_recuperator_valid(Optimizer, arch, EPS_REC_MIN):
                                print(f"    [!] Récupérateur quasi bypassé "
                                      f"(ε_rec < {EPS_REC_MIN:.0%}) -> run disqualifié")
                            else:
                                run_data = extract_run_data(Optimizer, arch)
                                print(f"    -> NTU_weighted_sum = {run_data['NTU_weighted_sum']}")

                                update_best_runs(condition_key, run_data)
                                write_best_runs_csv()  # réécriture à chaque run -> robuste à une interruption

                        except Exception as e:
                            print(f"    ⚠️ Échec : {e}")

        print(f"\n{len(best_runs)} conditions (architecture, T_H, eta_obj) traitées.")
        print(f"Top {BEST_KEEP} runs par condition enregistrés dans {CSV_PATH}")

        # ---- Graphiques ----
        plot_results(CSV_PATH)