#%%

# -*- coding: utf-8 -*-
"""
co2_rc_full_design_optimizer.py
Étend CO2RCOptimizer (importé de co2_rc_pso_optimizer.py) avec le
dimensionnement des composants + CAPEX + boucle cycle_design + log CSV
des résultats (CAPEX, puissance nette, efficacité, efficacités d'échangeurs).

Architecture : les objets de sizing (Recuperator, GasHeater, Condenser,
Pump, Expander_Axial, Expander_Radial) sont créés et configurés une seule
fois, dans le __main__, avec tous leurs paramètres/bornes/corrélations
statiques + un attribut RUN_KWARGS (kwargs à passer à .sizing()).
size_all_components() boucle dessus et n'y injecte, à chaque appel, que ce
qui dépend du point de fonctionnement courant : les inputs thermo (T/P/mdot)
et, pour les échangeurs, les contraintes Q_dot/DP_h/DP_c.

PERFORMANCE (voir discussion) :
- cycle_design() warm-démarre désormais chaque itération PSO avec la
  meilleure position trouvée à l'itération précédente (au lieu de repartir
  du même init_pos, ou de rien, à chaque fois) : entre deux itérations,
  seuls eta_pp/eta_exp/DP_* changent légèrement, donc l'optimum ne se
  déplace que peu -> convergence bien plus rapide.
- Voir le __main__ pour la correction patience/max_iter et le ntop réduit.
"""

#%% Imports

import csv
import os
import time
from datetime import datetime

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.sizing.turbomachinery.turbine.axial.sizing_1D.mean_line_axial_turbine_loss_model_sizing import AxialTurbineMeanLineSizing
from labothappy.sizing.turbomachinery.turbine.radial.mean_line_radial_turbine_loss_model_sizing import RadialTurbineMeanLineSizing
from labothappy.sizing.heat_exchanger.shell_and_tube.shell_and_tube_sizing import ShellAndTubeSizingOpt
from labothappy.sizing.heat_exchanger.PCHE.PCHE_sizing import PCHESizingOpt
from labothappy.sizing.turbomachinery.pump.radial.radial_pump_0D_sizing import RadialPumpODSizing

# --- Import de la brique d'optimisation (fichier 1) ---
from labothappy.machine.optimization.thermodynamic.CO2_RC_HX_presize_Optimization import CO2RC_HX_optimizer

import warnings
warnings.filterwarnings('ignore')

#%% Extraction des entrées dynamiques (dépendent du point de fonctionnement courant)

def _hx_inputs(model):
    return dict(
        fluid_H=model.su_H.fluid, T_su_H=model.su_H.T, P_su_H=model.su_H.p, m_dot_H=model.su_H.m_dot,
        fluid_C=model.su_C.fluid, T_su_C=model.su_C.T, P_su_C=model.su_C.p, m_dot_C=model.su_C.m_dot,
    )

def _pump_inputs(model):
    return dict(P_su=model.su.p, P_ex=model.ex.p, T_su=model.su.T,
                H1=0, H2=0, v1=0, v2=0, m_dot=model.su.m_dot)

def _turbine_inputs(model):
    return dict(mdot=model.su.m_dot, W_dot=model.W.W_dot,
                p0_su=model.su.p, T0_su=model.su.T, p_ex=model.ex.p)

# Registre structurel : quel extracteur utiliser pour chaque clé de composant.
# 'Expander' est traité à part (choix axial/radial, voir size_all_components).
DYNAMIC_INPUT_EXTRACTORS = {
    'Recuperator': _hx_inputs,
    'GasHeater': _hx_inputs,
    'Condenser': _hx_inputs,
    'Pump': _pump_inputs,
}

# Seuil plancher pour DP_h/DP_c (mêmes valeurs que l'original max(DP_h, 1e3)/(1e4)).
# Recuperator n'avait pas de plancher dans l'original (0 = pas de max()).
HX_DP_FLOOR = {'Recuperator': 0.0, 'GasHeater': 1e3, 'Condenser': 1e4}

# GasHeater/Condenser ont aussi besoin de T_max_cycle/p_max_cycle
# (routés via set_parameters -> _apply_deferred_parameters -> set_max_cycle_prop).
# Recuperator (PCHE) ne les utilise pas.
HX_SOURCE_KEY = {'GasHeater': 'GH_Water', 'Condenser': 'CD_Water'}


def _set_dynamic_hx_constraints(sizing_obj, model, RC, key):
    dp_floor = HX_DP_FLOOR[key]
    sizing_obj.set_parameters(
        Q_dot=model.Q.Q_dot,
        DP_h=max(model.DP_h, dp_floor),
        DP_c=max(model.DP_c, dp_floor),
    )
    if key in HX_SOURCE_KEY:
        p_max_cycle = RC.components['Pump'].model.ex.p
        T_max_cycle = RC.sources[HX_SOURCE_KEY[key]].properties.T

        if p_max_cycle is None or T_max_cycle is None:
            raise ValueError(
                f"{key}: p_max_cycle/T_max_cycle indisponible "
                f"(Pump du cycle non convergé — p_max_cycle={p_max_cycle}, T_max_cycle={T_max_cycle})"
            )

        sizing_obj.set_parameters(T_max_cycle=T_max_cycle, p_max_cycle=p_max_cycle)


def size_all_components(RC, sizing_models, turb_choice="None"):
    """
    Boucle sur `sizing_models` (préconfigurés dans le main). N'y injecte,
    à chaque appel, que ce qui dépend du point de fonctionnement courant.
    Chaque sizing_obj porte un attribut `.RUN_KWARGS` (posé dans le main)
    avec les kwargs à passer à `.sizing()`.

    Retourne (ok, results, turb_choice) où results = {key: sizing_obj}.
    En cas d'échec, ok=False et le premier composant en échec est signalé.
    """
    results = {}

    # --- Composants "simples" : un seul sizing_obj par clé ---
    for key, sizing_obj in sizing_models.items():
        if key.startswith('Expander'):
            continue  # traité à part, plus bas

        model = RC.components[key].model
        RC.components[key].sizing = sizing_obj

        try:
            sizing_obj.set_inputs(**DYNAMIC_INPUT_EXTRACTORS[key](model))

            if key in HX_DP_FLOOR:
                _set_dynamic_hx_constraints(sizing_obj, model, RC, key)

            sizing_obj.sizing(**sizing_obj.RUN_KWARGS)

        except Exception as e:
            print(f"⚠️ Failed to design {key}: {e}")
            if hasattr(model, 'su_H'):
                model.su_H.print_resume()
                model.su_C.print_resume()
                print(f"Q_dot_cstr : {model.Q.Q_dot}")
                print(f"DP_h_cstr : {model.DP_h}")
                print(f"DP_c_cstr : {model.DP_c}")
            return False, results, "Fail"

        results[key] = sizing_obj

    # --- Turbine : choix axial vs radial (deux sizing_obj candidats) ---
    Turb_model = RC.components['Expander'].model
    turb_inputs = _turbine_inputs(Turb_model)
    eta_axial = eta_radial = 0
    Turb_axial_sizing = Turb_radial_sizing = None

    if turb_choice != 'Radial':
        try:
            Turb_axial_sizing = sizing_models['Expander_Axial']
            Turb_axial_sizing.set_inputs(**turb_inputs)
            Turb_axial_sizing.sizing(**Turb_axial_sizing.RUN_KWARGS)
            eta_axial = Turb_axial_sizing.eta_is
        except Exception as e:
            print(f"⚠️ Failed to design the axial Turbine: {e}")

    if turb_choice != 'Axial':
        try:
            Turb_radial_sizing = sizing_models['Expander_Radial']
            Turb_radial_sizing.set_inputs(**turb_inputs)
            Turb_radial_sizing.sizing(**Turb_radial_sizing.RUN_KWARGS)
            eta_radial = Turb_radial_sizing.eta_is
        except Exception as e:
            print(f"⚠️ Failed to design the radial Turbine: {e}")

    if eta_axial == 0 and eta_radial == 0:
        return False, results, "Fail"

    if eta_axial > eta_radial:
        RC.components['Expander'].sizing = results['Expander'] = Turb_axial_sizing
        turb_choice = "Axial"
    else:
        RC.components['Expander'].sizing = results['Expander'] = Turb_radial_sizing
        turb_choice = "Radial"

    print(f"eta_axial : {eta_axial}")
    print(f"eta_radial : {eta_radial}")

    return True, results, turb_choice

#%% Logging des résultats

def _hx_effectivenesses(RC, arch):
    """
    Recalcule/relit les epsilon ET les Q_dot des échangeurs, exactement
    comme dans system_RC_parallel (source de vérité pour ces valeurs).
    Renvoie un dict {nom_lisible: valeur}, avec NaN pour ce qui n'existe
    pas / échoue. Les Q_dot sont nécessaires pour calibrer ultérieurement
    les coefficients cost_w_* (voir calibrate_cost_weights_from_sizing ci-dessous).
    """
    out = {}

    def safe_epsilon(component_key, label):
        try:
            model = RC.components[component_key].model
            out[label] = model.epsilon
        except Exception:
            out[label] = float("nan")

    def safe_Q(component_key, label):
        try:
            model = RC.components[component_key].model
            out[label] = model.Q.Q_dot
        except Exception:
            out[label] = float("nan")

    # Condenser : epsilon n'est peuplé qu'après cet appel explicite
    try:
        RC.components['Condenser'].model.equivalent_effectiveness()
    except Exception:
        pass
    safe_epsilon('Condenser', 'eps_cond')
    safe_Q('Condenser', 'Q_cond')

    safe_epsilon('GasHeater', 'eps_gh')
    safe_Q('GasHeater', 'Q_gh')

    if arch == 'REC':
        safe_epsilon('Recuperator', 'eps_rec')
        safe_Q('Recuperator', 'Q_rec')
    elif arch == 'Recomp':
        safe_epsilon('RecupLT', 'eps_rec_LT')
        safe_Q('RecupLT', 'Q_rec_LT')
        safe_epsilon('RecupHT', 'eps_rec_HT')
        safe_Q('RecupHT', 'Q_rec_HT')
    elif arch == 'Recomp_1_recup':
        safe_epsilon('RecupLT', 'eps_rec_LT')
        safe_Q('RecupLT', 'Q_rec_LT')
    # arch == 'basic' : pas de récupérateur

    return out


def log_cycle_result(log_path, T_hot, T_cold, W_dot_obj, eta_obj, RC, arch,
                      Optimizer=None, duration_s=None, run_id=None):
    """
    Ajoute une ligne à un fichier CSV de log : CAPEX (total + détail par
    composant), puissance nette et efficacité atteintes, efficacités des
    échangeurs, en fonction de T_hot, T_cold, de l'objectif de puissance et
    de l'efficacité cible.

    `RC` est le cycle dimensionné (typiquement Optimizer.best_RC).

    Puissance nette et efficacité sont recalculées depuis les composants
    dimensionnés (pas les valeurs cibles) :
      W_dot_net = W_dot_expander - W_dot_pump
      eta       = W_dot_net / Q_dot_GasHeater
    """

    def safe_get(fn, default=float("nan")):
        try:
            return fn()
        except Exception:
            return default

    W_dot_exp = safe_get(lambda: RC.components['Expander'].sizing.W_dot)
    W_dot_pp = safe_get(lambda: RC.components['Pump'].sizing.W_dot)
    Q_dot_gh = safe_get(lambda: RC.components['GasHeater'].sizing.best_particle.Q)

    W_dot_net = safe_get(lambda: W_dot_exp - W_dot_pp)
    eta_achieved = safe_get(lambda: W_dot_net / Q_dot_gh if Q_dot_gh else float("nan"))

    eta_carnot = safe_get(lambda: 1 - T_cold / T_hot)
    eta_vs_carnot = safe_get(lambda: eta_achieved / eta_carnot if eta_carnot else float("nan"))

    capex_total = RC.CAPEX.get("Total", float("nan"))
    capex_specific = safe_get(lambda: capex_total / (W_dot_net / 1e3) if W_dot_net else float("nan"))  # $/kW

    turb_choice = safe_get(lambda: RC.components['Expander'].sizing.__class__.__name__)

    row = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "run_id": run_id,
        "duration_s": duration_s,

        "T_hot_C": T_hot - 273.15,
        "T_cold_C": T_cold - 273.15,
        "W_dot_obj_MW": W_dot_obj / 1e6,
        "eta_obj": eta_obj,

        "P_high_Pa": safe_get(lambda: RC.it_var.get('P_high')),
        "mdot_kg_s": safe_get(lambda: RC.it_var.get('mdot')),
        "mdot_HS_kg_s": safe_get(lambda: RC.it_var.get('mdot_HS')),
        "mdot_CS_kg_s": safe_get(lambda: RC.it_var.get('mdot_CS')),

        "W_dot_achieved_MW": safe_get(lambda: W_dot_net / 1e6),
        "eta_achieved": eta_achieved,
        "eta_carnot": eta_carnot,
        "eta_vs_carnot": eta_vs_carnot,

        "CAPEX_total": capex_total,
        "CAPEX_specific_USD_per_kW": capex_specific,

        "turbine_type": turb_choice,
        "eta_is_pump": safe_get(lambda: RC.components['Pump'].sizing.eta_is),
        "eta_is_expander": safe_get(lambda: RC.components['Expander'].sizing.eta_is),

        "DP_h_rec_Pa": safe_get(lambda: RC.components['Recuperator'].sizing.HX.DP_h),
        "DP_c_rec_Pa": safe_get(lambda: RC.components['Recuperator'].sizing.HX.DP_c),
        "DP_h_gh_Pa": safe_get(lambda: RC.components['GasHeater'].sizing.best_particle.DP_h),
        "DP_c_gh_Pa": safe_get(lambda: RC.components['GasHeater'].sizing.best_particle.DP_c),
        "DP_h_cond_Pa": safe_get(lambda: RC.components['Condenser'].sizing.best_particle.DP_h),
        "DP_c_cond_Pa": safe_get(lambda: RC.components['Condenser'].sizing.best_particle.DP_c),

        "n_iter_cycle_design": safe_get(lambda: Optimizer.criterion) if Optimizer is not None else None,
    }

    # Efficacités des échangeurs, telles qu'utilisées dans l'objectif du PSO
    row.update(_hx_effectivenesses(RC, arch))

    # Détail CAPEX par composant
    for key, value in RC.CAPEX.items():
        if key != "Total":
            row[f"CAPEX_{key}"] = value

    file_exists = os.path.isfile(log_path)
    with open(log_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

def calibrate_cost_weights_from_sizing(Optimizer, arch='REC', RC_list=None, reference='cond'):
    """
    Calibre cost_w_gh / cost_w_rec / cost_w_cond DIRECTEMENT depuis le CAPEX
    réel des échangeurs déjà sizés dans ce run -- pas de CSV, pas de
    régression multi-runs, pas de modèle UA/LMTD intermédiaire.

    Principe, par échangeur x : cost_w_x = (CAPEX_x / (Q_x * NTU_x)),
    normalisé par rapport à `reference` (donc cost_w_<reference> = 1.0). Ce
    ratio est directement injectable dans l'objectif du PSO thermodynamique
    (system_RC_parallel), qui calcule déjà Q_x*NTU_x pour chaque échangeur
    -- cost_w_x agit comme un simple facteur multiplicatif dessus.

    `RC_list` : liste de cycles déjà sizés (chacun avec RC.CAPEX[key] et
    RC.components[key].model peuplés). Par défaut, Optimizer.potential_RC
    (les candidats sizés à l'itération courante de cycle_design). Peut aussi
    recevoir la liste retournée par quick_calibration_points() pour calibrer
    depuis plusieurs points choisis à la main, en amont d'un cycle_design.

    Retourne un dict {'cost_w_gh':..., 'cost_w_rec':..., 'cost_w_cond':...}
    prêt à passer à Optimizer.set_parameters(**dict), ou None si aucun
    échangeur n'a pu être calibré.
    """
    if arch != 'REC':
        print(f"[calibrate_cost_weights_from_sizing] Architecture '{arch}' "
              f"non supportée (implémenté pour REC uniquement).")
        return None

    if RC_list is None:
        RC_list = Optimizer.potential_RC

    key_map = {'gh': 'GasHeater', 'rec': 'Recuperator', 'cond': 'Condenser'}
    ratios = {'gh': [], 'rec': [], 'cond': []}

    for RC in RC_list:
        for short, key in key_map.items():
            try:
                capex = RC.CAPEX.get(key)
                model = RC.components[key].model
                Q = model.Q.Q_dot
                if key == 'Condenser':
                    model.equivalent_effectiveness()
                eps = model.epsilon
                ntu = -np.log(1.0 - min(max(eps, 0.0), 1.0 - 1e-6))
                if capex and Q and ntu > 0:
                    ratios[short].append(capex / (Q * ntu))
            except Exception:
                continue

    ks = {}
    for short, vals in ratios.items():
        if vals:
            ks[short] = float(np.mean(vals))
            flag = "  (n=1, aucune robustesse)" if len(vals) == 1 else ""
            print(f"[calibrate_cost_weights_from_sizing] k_{short} = {ks[short]:.4g} "
                  f"(n={len(vals)}){flag}")
        else:
            print(f"[calibrate_cost_weights_from_sizing] Aucun point valide pour "
                  f"'{short}' -- cost_w_{short} laissé inchangé.")

    if reference not in ks:
        print(f"[calibrate_cost_weights_from_sizing] Référence '{reference}' "
              f"non calibrable -- abandon.")
        return None

    k_ref = ks[reference]
    weights = {f'cost_w_{short}': v / k_ref for short, v in ks.items()}
    print(f"[calibrate_cost_weights_from_sizing] cost_w_* (réf={reference}) : {weights}")
    return weights


def quick_calibration_points(Optimizer, sizing_models, points):
    """
    Size directement quelques positions thermodynamiques choisies à la main
    (PAS de PSO) pour obtenir de vrais CAPEX_gh/CAPEX_rec/CAPEX_cond en
    quelques minutes plutôt qu'en lançant une campagne cycle_design complète.

    `Optimizer` doit déjà avoir set_parameters/set_it_var/set_obj/
    set_CSource/set_HSource/set_RC appelés (comme dans le __main__).

    `points` : liste de dicts, ex. pour l'architecture REC :
        [{'P_high': 150e5, 'mdot': 200.0, 'mdot_HS': 150.0, 'mdot_CS': 2000.0,
          'eta_gh': 0.92, 'PP_gh': 5, 'eta_rec': 0.85, 'PP_cd': 5}, ...]
    Idéalement 3-5 points couvrant des Q_dot/NTU variés par échangeur, pour
    que calibrate_cost_weights_from_sizing() ait de quoi moyenner plutôt
    qu'un seul point isolé.

    Retourne la liste des cycles sizés avec succès (à passer directement à
    calibrate_cost_weights_from_sizing(Optimizer, RC_list=...)).
    """
    sized_RCs = []
    for i, point in enumerate(points):
        print(f"\n[quick_calibration] Point {i+1}/{len(points)} : {point}")

        Optimizer.it_var.update(point)
        Optimizer._HSource_props['m_dot'] = point['mdot_HS']
        Optimizer._CSource_props['m_dot'] = point['mdot_CS']

        try:
            Optimizer.set_RC()
            Optimizer.current_RC = Optimizer.RC
            Optimizer.current_RC.solve()
        except Exception as e:
            print(f"  ⚠️ Solve échoué pour ce point : {e}")
            continue

        ok, results, turb_choice = size_all_components(
            Optimizer.current_RC, sizing_models, Optimizer.turb_choice
        )
        if not ok:
            print("  ⚠️ Sizing échoué pour ce point, ignoré.")
            continue

        Optimizer.current_RC.CAPEX = {key: np.round(obj.CAPEX['Total']) for key, obj in results.items()}
        Optimizer.current_RC.CAPEX["Total"] = sum(Optimizer.current_RC.CAPEX.values())
        sized_RCs.append(Optimizer.current_RC)

    return sized_RCs



#%% Classe étendue : hérite de la brique d'optimisation importée

class CO2RCOptimizer(CO2RC_HX_optimizer):
    """
    Étend CO2RCOptimizer (co2_rc_pso_optimizer.py) avec :
    - dimensionnement des composants (size_all_components, via self.sizing_models)
    - calcul CAPEX
    - boucle itérative cycle_design (opt → size → ré-estime params → repeat)

    L'optimisation PSO (system_RC_parallel, set_RC, _evaluate_final, opt_RC)
    est intégralement héritée du module importé — non réécrite ici.
    """

    def __init__(self, fluid):
        super().__init__(fluid)
        self.CAPEX = {}
        self.turb_choice = "None"
        self.potential_RC = []
        self.best_RC = None
        self.sizing_models = {}   # <-- rempli depuis le main avant cycle_design()
        # Snapshot des DP_* ENCODÉS AVANT le lancement de l'optimisation (pris
        # une seule fois, au début de cycle_design). Sert de référence FIXE
        # pour le score de cohérence dans evaluate_systems() -- contrairement
        # à self.params['DP_h_gh'] etc. qui dérivent à chaque itération
        # (moyenne avec le DP réel), ce snapshot ne bouge jamais, donc un DP
        # réel qui reste sous la valeur d'origine n'est jamais pénalisé,
        # même après plusieurs itérations de mise à jour de self.params.
        self.original_DP_params = {}
        # Meilleur design (CAPEX le plus bas) jamais vu sur TOUT le run de
        # cycle_design, toutes itérations confondues -- voir discussion
        # "que se passe-t-il si le design redevient plus cher après une
        # itération". self.best_RC (seul, sans "_overall") reste le design
        # de la DERNIÈRE itération traitée par evaluate_systems() -- utile
        # pour la mise à jour de self.params, qui doit continuer à se baser
        # sur l'itération courante même si elle régresse en coût. C'est
        # self.best_RC_overall qui est restauré dans self.best_RC à la fin
        # de cycle_design(), pour que le résultat final ne soit jamais pire
        # que ce qui a déjà été trouvé en cours de route.
        self.best_RC_overall = None
        self.best_capex_overall = float('inf')
        self.capex_history = []  # [(iteration, CAPEX_Total), ...] sur tout le run

    def evaluate_systems(self):
        """
        Score chaque candidat sizé sur deux critères combinés :
          1) Cohérence :
             - eta_exp/eta_pp : écart entre la valeur ASSUMÉE dans le PSO
               thermodynamique (self.params, mise à jour à chaque itération)
               et celle RÉELLEMENT ATTEINTE après sizing.
             - DP_* : écart entre le DP RÉEL et le DP ORIGINAL encodé avant
               le lancement de l'optimisation (self.original_DP_params,
               figé une fois pour toutes -- PAS self.params qui dérive
               d'itération en itération). Le score n'est non nul QUE si le
               DP réel DÉPASSE cette valeur d'origine (voir le np.max ci-
               dessous) : un DP réel plus bas que ce qui a été encodé au
               départ n'est jamais pénalisé, même après que self.params ait
               dérivé vers une cible plus basse au fil des itérations. Ce
               score reste calculé et tracé (delta_dict, logs) dans tous les
               cas, il ne bloque simplement plus la convergence quand il est
               favorable.
          2) CAPEX : coût total du candidat, normalisé par le CAPEX minimum
             du lot (donc le moins cher a une contribution nulle, les autres
             sont pénalisés proportionnellement à leur surcoût relatif).

        `capex_weight` (dans self.params, défaut 1.0) règle l'importance
        relative du coût : 0 = comportement d'origine (cohérence seule),
        valeurs plus élevées = priorité croissante au coût. À 1.0, un
        candidat 20 % plus cher qu'un autre à cohérence égale perd si son
        écart de cohérence est inférieur à 0.2 (échelle des deltas au carré
        ci-dessous, typiquement de l'ordre de 1e-3 à 1e-1).
        """
        RC_scores = []
        delta_dicts = []
        capex_totals = []

        # Référence DP : self.original_DP_params si déjà initialisé par
        # cycle_design(), sinon repli sur self.params (permet d'appeler
        # evaluate_systems() de façon autonome, hors cycle_design).
        dp_ref = self.original_DP_params if self.original_DP_params else self.params

        for RC in self.potential_RC:
            delta_dicts.append({})

            eta_exp = RC.components['Expander'].sizing.eta_is
            eta_pp = RC.components['Pump'].sizing.eta_is

            DP_h_gh = RC.components['GasHeater'].sizing.best_particle.DP_h
            DP_c_gh = RC.components['GasHeater'].sizing.best_particle.DP_c

            DP_h_cond = RC.components['Condenser'].sizing.best_particle.DP_h
            DP_c_cond = RC.components['Condenser'].sizing.best_particle.DP_c

            DP_h_rec = RC.components['Recuperator'].sizing.HX.DP_h
            DP_c_rec = RC.components['Recuperator'].sizing.HX.DP_c

            delta_dicts[-1]['eta_exp'] = delta_exp = ((eta_exp - self.params['eta_exp']) / self.params['eta_exp']) ** 2
            delta_dicts[-1]['eta_pp'] = delta_pp = ((eta_pp - self.params['eta_pp']) / self.params['eta_pp']) ** 2

            delta_dicts[-1]['DP_h_gh'] = delta_h_gh = ((np.max([DP_h_gh, dp_ref['DP_h_gh']]) - dp_ref['DP_h_gh']) / dp_ref['DP_h_gh']) ** 2
            delta_dicts[-1]['DP_c_gh'] = delta_c_gh = ((np.max([DP_c_gh, dp_ref['DP_c_gh']]) - dp_ref['DP_c_gh']) / dp_ref['DP_c_gh']) ** 2

            delta_dicts[-1]['DP_h_cond'] = delta_h_cond = ((np.max([DP_h_cond, dp_ref['DP_h_cond']]) - dp_ref['DP_h_cond']) / dp_ref['DP_h_cond']) ** 2
            delta_dicts[-1]['DP_c_cond'] = delta_c_cond = ((np.max([DP_c_cond, dp_ref['DP_c_cond']]) - dp_ref['DP_c_cond']) / dp_ref['DP_c_cond']) ** 2

            delta_dicts[-1]['DP_h_rec'] = delta_h_rec = ((np.max([DP_h_rec, dp_ref['DP_h_rec']]) - dp_ref['DP_h_rec']) / dp_ref['DP_h_rec']) ** 2
            delta_dicts[-1]['DP_c_rec'] = delta_c_rec = ((np.max([DP_c_rec, dp_ref['DP_c_rec']]) - dp_ref['DP_c_rec']) / dp_ref['DP_c_rec']) ** 2

            score_current = delta_exp + delta_pp + delta_h_gh + delta_c_gh + delta_h_cond + delta_c_cond + delta_h_rec + delta_c_rec
            RC_scores.append(score_current)
            capex_totals.append(RC.CAPEX.get("Total", float("nan")))

        # --- pénalité CAPEX, ajoutée au score de cohérence ---
        capex_weight = self.params.get('capex_weight', 1.0)
        valid_capex = [c for c in capex_totals if np.isfinite(c)]
        if capex_weight > 0 and valid_capex:
            capex_min = min(valid_capex)
            for i, capex in enumerate(capex_totals):
                if np.isfinite(capex) and capex_min > 0:
                    RC_scores[i] += capex_weight * (capex / capex_min - 1.0)

        index_of_min = RC_scores.index(np.min(RC_scores))
        self.best_RC = best_RC = self.potential_RC[index_of_min]
        delta_dict = delta_dicts[index_of_min]

        new_params_dict = {
            'eta_exp': np.round(best_RC.components['Expander'].sizing.eta_is, 3),
            'eta_pp': np.round(best_RC.components['Pump'].sizing.eta_is, 3),
            'DP_h_gh': np.round(best_RC.components['GasHeater'].sizing.best_particle.DP_h),
            'DP_c_gh': np.round(best_RC.components['GasHeater'].sizing.best_particle.DP_c),
            'DP_h_cond': np.round(best_RC.components['Condenser'].sizing.best_particle.DP_h),
            'DP_c_cond': np.round(best_RC.components['Condenser'].sizing.best_particle.DP_c),
            'DP_h_rec': np.round(best_RC.components['Recuperator'].sizing.HX.DP_h),
            'DP_c_rec': np.round(best_RC.components['Recuperator'].sizing.HX.DP_c),
        }

        return new_params_dict, np.min(RC_scores), delta_dict

    def size_components(self):
        i = 0
        n_pos = len(self.top_positions)
        self.potential_RC = []
        turb_choices = []

        if self.obj['W_dot'] >= 9e6:
            self.turb_choice = 'Axial'

        for allowable_position in self.top_positions:
            print(f"Component Optimization for top position : {i+1}/{n_pos}")
            i += 1

            # Réutilise le helper hérité de co2_rc_pso_optimizer.py
            unpacked = self._unpack_position(allowable_position['x'])
            self.it_var.update(unpacked)

            self._HSource_props['m_dot'] = unpacked['mdot_HS']
            self._CSource_props['m_dot'] = unpacked['mdot_CS']

            try:
                self.set_RC()  # hérité : construit self.RC selon l'architecture
                self.current_RC = self.RC
                self.current_RC.solve()
            except Exception as e:
                print(f"⚠️ Failed to solve final RC circuit: {e}")
                continue

            ok, results, turb_choice = size_all_components(
                self.current_RC, self.sizing_models, self.turb_choice
            )

            if ok:
                self.current_RC.CAPEX = {key: np.round(obj.CAPEX['Total']) for key, obj in results.items()}
                self.current_RC.CAPEX["Total"] = sum(self.current_RC.CAPEX.values())
                self.potential_RC.append(self.current_RC)

            turb_choices.append(turb_choice)

        filtered = [c for c in turb_choices if c in ("Axial", "Radial")]
        if filtered:
            axial_count = filtered.count("Axial")
            radial_count = filtered.count("Radial")
            self.turb_choice = "Axial" if axial_count >= radial_count else "Radial"
            print("Most frequent choice:", self.turb_choice)
        else:
            print("No valid choices (Axial or Radial) found.")

    def cycle_design(self, n_jobs=None, n_particles=50, max_iter=30, patience=10,
                      ntop=5, init_pos=None):
        import multiprocessing as mp
        n_cores = mp.cpu_count()
        if n_jobs is None:
            n_jobs = n_cores - 1

        self.criterion = 0
        n_it_max = 10
        it = 0

        # Snapshot fixe des DP_* tels qu'encodés avant le lancement (voir
        # docstring d'evaluate_systems). Pris une seule fois -- si
        # cycle_design() est rappelée sur un Optimizer déjà utilisé, ne
        # réinitialise pas le snapshot (self.original_DP_params déjà non
        # vide), pour ne pas effacer la référence d'un run précédent.
        if not self.original_DP_params:
            self.original_DP_params = {
                'DP_h_gh': self.params['DP_h_gh'],
                'DP_c_gh': self.params['DP_c_gh'],
                'DP_h_cond': self.params['DP_h_cond'],
                'DP_c_cond': self.params['DP_c_cond'],
                'DP_h_rec': self.params['DP_h_rec'],
                'DP_c_rec': self.params['DP_c_rec'],
            }
            print(f"[cycle_design] DP de référence (fixes, pour tout le run) : "
                  f"{self.original_DP_params}")

        while self.criterion == 0 and it < n_it_max:

            self.allowable_positions = []  # reset à chaque itération

            # --- Optimisation PSO héritée du fichier 1 ---
            # PERFORMANCE : on capture l'objet `optimizer` retourné pour
            # réutiliser sa meilleure position comme graine (warm start) de
            # l'itération suivante. D'une itération de cycle_design à
            # l'autre, seuls eta_pp/eta_exp/DP_* bougent légèrement (moyenne
            # avec les anciennes valeurs, cf. plus bas) -> le nouvel optimum
            # est proche de l'ancien, donc le PSO reconverge en bien moins
            # d'itérations si on le démarre depuis là plutôt qu'à froid.
            pso_optimizer = self.opt_RC(n_jobs=n_jobs, n_particles=n_particles, max_iter=max_iter,
                        patience=patience, ntop=ntop, init_pos=init_pos)
            init_pos = pso_optimizer.swarm.best_pos

            # --- Dimensionnement des composants ---
            self.size_components()

            # ------------------------------------------------------------
            # Calibre cost_w_gh/rec/cond directement depuis le CAPEX réel des
            # échangeurs sizés à cette itération, et met à jour self.params
            # pour que le PSO thermodynamique des itérations suivantes en
            # tienne compte. Se recalibre à CHAQUE itération (les candidats
            # sizés changent d'une itération à l'autre, donc autant garder
            # les poids à jour plutôt que figés après la 1ère calibration).
            # ------------------------------------------------------------
            if self.potential_RC:
                cost_weights = calibrate_cost_weights_from_sizing(
                    self, self.params.get('RC_ARCH', 'REC')
                )
                if cost_weights is not None:
                    self.set_parameters(**cost_weights)
                    print(f"[cycle_design] cost_w_* mis à jour : {cost_weights}")
                else:
                    print("[cycle_design] Calibration cost_w_* impossible ce coup-ci "
                          "-- valeurs précédentes conservées.")

            new_params, best_score, delta_dict = self.evaluate_systems()
            self.new_params = new_params
            self.delta_dict = delta_dict

            # ------------------------------------------------------------
            # Suivi du meilleur design jamais vu (CAPEX), toutes itérations
            # confondues -- voir discussion "que faire si le design devient
            # plus cher après une itération". Silencieux : pas de message,
            # juste la mise à jour de self.best_RC_overall/self.capex_history.
            # ------------------------------------------------------------
            current_capex = self.best_RC.CAPEX.get('Total', float('inf'))
            self.capex_history.append((it + 1, current_capex))
            if current_capex < self.best_capex_overall:
                self.best_RC_overall = self.best_RC
                self.best_capex_overall = current_capex

            # ------------------------------------------------------------
            # Bloc d'affichage UNIQUE (fusionne score/CAPEX, cohérence
            # turbomachines, DP réel vs référence, et variables d'échangeurs
            # cible-PSO vs effectivité réellement atteinte sur best_RC).
            # ------------------------------------------------------------
            arch = self.params.get('RC_ARCH', 'REC')

            def _ntu(eps):
                if eps is None or not np.isfinite(eps):
                    return float('nan')
                eps_c = min(max(eps, 0.0), 1.0 - 1e-6)
                return -np.log(1.0 - eps_c)

            hx_eff = _hx_effectivenesses(self.best_RC, arch)

            print("\n" + "="*60)
            print(f"  cycle_design -- itération {it + 1}")
            print("="*60)
            print(f"  Best Score (cohérence + CAPEX) : {best_score:.6g}")
            print(f"  CAPEX Total (best_RC)          : {self.best_RC.CAPEX.get('Total'):,.0f}")
            print(f"  CAPEX candidats sizés           : "
                  f"{[rc.CAPEX.get('Total') for rc in self.potential_RC]}")

            print("-"*60)
            print("  Cohérence turbomachines (assumé -> réel, delta %) :")
            print(f"    eta_exp : {self.params['eta_exp']:.3f} -> {new_params['eta_exp']:.3f}  "
                  f"(delta={delta_dict['eta_exp']*100:.4f}%)")
            print(f"    eta_pp  : {self.params['eta_pp']:.3f} -> {new_params['eta_pp']:.3f}  "
                  f"(delta={delta_dict['eta_pp']*100:.4f}%)")

            print("-"*60)
            print("  DP réel (best_RC) vs référence fixe (pré-optimisation) :")
            for dp_key in ('DP_h_gh', 'DP_c_gh', 'DP_h_cond', 'DP_c_cond', 'DP_h_rec', 'DP_c_rec'):
                real_val = new_params.get(dp_key)
                ref_val = self.original_DP_params.get(dp_key)
                if real_val is not None and ref_val:
                    ratio = real_val / ref_val
                    tag = "OK (sous la référence)" if real_val <= ref_val else "⚠️ AU-DESSUS de la référence"
                    print(f"    {dp_key:10s} : réel={real_val:10.0f} Pa  |  "
                          f"référence={ref_val:10.0f} Pa  |  ratio={ratio:.2f}  |  {tag}")

            print("-"*60)
            print("  Échangeurs -- variable cible du PSO vs effectivité réelle (best_RC) :")
            eta_gh_cible = self.it_var.get('eta_gh')
            PP_gh_cible  = self.it_var.get('PP_gh')
            PP_cd_cible  = self.it_var.get('PP_cd')
            print(f"    GasHeater   : cible eta_gh={eta_gh_cible}, PP_gh={PP_gh_cible} K  |  "
                  f"eps réel={hx_eff.get('eps_gh'):.4f}  |  NTU réel={_ntu(hx_eff.get('eps_gh')):.3f}  |  "
                  f"Q={hx_eff.get('Q_gh', float('nan')):,.0f} W")

            if arch == 'REC':
                eta_rec_cible = self.it_var.get('eta_rec')
                print(f"    Recuperator : cible eta_rec={eta_rec_cible}  |  "
                      f"eps réel={hx_eff.get('eps_rec'):.4f}  |  NTU réel={_ntu(hx_eff.get('eps_rec')):.3f}  |  "
                      f"Q={hx_eff.get('Q_rec', float('nan')):,.0f} W")
            elif arch == 'Recomp':
                print(f"    RecupLT     : cible eta_rec_LT={self.it_var.get('eta_rec_LT')}  |  "
                      f"eps réel={hx_eff.get('eps_rec_LT'):.4f}  |  NTU réel={_ntu(hx_eff.get('eps_rec_LT')):.3f}  |  "
                      f"Q={hx_eff.get('Q_rec_LT', float('nan')):,.0f} W")
                print(f"    RecupHT     : cible eta_rec_HT={self.it_var.get('eta_rec_HT')}  |  "
                      f"eps réel={hx_eff.get('eps_rec_HT'):.4f}  |  NTU réel={_ntu(hx_eff.get('eps_rec_HT')):.3f}  |  "
                      f"Q={hx_eff.get('Q_rec_HT', float('nan')):,.0f} W")
            elif arch == 'Recomp_1_recup':
                print(f"    RecupLT     : cible eta_rec={self.it_var.get('eta_rec')}  |  "
                      f"eps réel={hx_eff.get('eps_rec_LT'):.4f}  |  NTU réel={_ntu(hx_eff.get('eps_rec_LT')):.3f}  |  "
                      f"Q={hx_eff.get('Q_rec_LT', float('nan')):,.0f} W")

            print(f"    Condenser   : cible PP_cd={PP_cd_cible} K  |  "
                  f"eps réel={hx_eff.get('eps_cond'):.4f}  |  NTU réel={_ntu(hx_eff.get('eps_cond')):.3f}  |  "
                  f"Q={hx_eff.get('Q_cond', float('nan')):,.0f} W")
            print("="*60)

            self.set_parameters(
                eta_exp=new_params['eta_exp'],
                eta_pp=new_params['eta_pp'],
                DP_h_gh=(new_params['DP_h_gh'] + self.params['DP_h_gh']) / 2,
                DP_c_gh=(new_params['DP_c_gh'] + self.params['DP_c_gh']) / 2,
                DP_h_cond=(new_params['DP_h_cond'] + self.params['DP_h_cond']) / 2,
                DP_c_cond=(new_params['DP_c_cond'] + self.params['DP_c_cond']) / 2,
                DP_h_rec=(new_params['DP_h_rec'] + self.params['DP_h_rec']) / 2,
                DP_c_rec=(new_params['DP_c_rec'] + self.params['DP_c_rec']) / 2,
            )

            self.criterion = 1
            for key in self.delta_dict:
                if self.delta_dict[key] > 1e-3:
                    self.criterion = 0
                    break
            it += 1

        # Restauration finale : si la dernière itération traitée n'est pas
        # celle qui a produit le meilleur CAPEX du run, self.best_RC (celui
        # de la dernière itération) est remplacé par self.best_RC_overall
        # (le meilleur jamais vu) -- pour que le design retourné par
        # cycle_design() ne soit jamais pire que ce qui a déjà été trouvé
        # en cours de route. self.new_params/self.delta_dict, en revanche,
        # restent ceux de la dernière itération (ils décrivent le processus
        # de convergence, pas le design final retenu).
        if self.best_RC_overall is not None:
            self.best_RC = self.best_RC_overall

        if self.params.get('save_file_path') is not None:
            import json, os

            class NumpyEncoder(json.JSONEncoder):
                def default(self, obj):
                    if isinstance(obj, np.integer):
                        return int(obj)
                    if isinstance(obj, np.floating):
                        return float(obj)
                    if isinstance(obj, np.ndarray):
                        return obj.tolist()
                    return super().default(obj)

            n_MW = int(self.obj["W_dot"] * 1e-6)
            eta = int(self.obj["eta"] * 100)
            T_hot = int(self._HSource_props['T'] - 273.15)
            T_cold = int(self._CSource_props['T'] - 273.15)

            folder_name = f"W{n_MW}_eta{eta}_TH{T_hot}_TC{T_cold}"
            save_folder = os.path.join(self.params['save_file_path'], folder_name)
            os.makedirs(save_folder, exist_ok=True)

            for component in self.best_RC.components:
                data = self.best_RC.components[component].sizing.export_params_dict()
                filepath = os.path.join(save_folder, f"{component}.json")
                with open(filepath, "w") as f:
                    json.dump(data, f, indent=4, cls=NumpyEncoder)

        return self

#%% Main

if __name__ == "__main__":

    # Cycle sizing parameters

    T_hot = 150 + 273.15
    T_cold = 10 + 273.15
    n_MW = 10
    W_dot_obj = n_MW * 1e6
    eta_obj = 0.12

    Optimizer = CO2RCOptimizer('CO2')

    m_dot_HS_fact_bounds = [0.1, 3]
    m_dot_CS_fact_bounds = [5, 15]
    P_high_bounds = np.array([110, 180]) * 1e5
    m_dot_bounds = np.array([10, 80]) * n_MW

    eta_gh_disc = np.arange(0.9, 0.98, 0.02)
    PP_gh_disc = np.arange(1, 10, 1)
    eta_rec_disc = np.arange(0.6, 0.96, 0.02)
    PP_cd_disc = np.arange(1, 10, 1)

    Optimizer.set_parameters(
        save_file_path=None,   # ou un chemin, comme dans le fichier 2 d'origine
        RC_ARCH='REC',          # seule architecture compatible avec le sizing actuel
        eta_pp=0.8,
        eta_pp_aux=0.8,
        DP_h_gh=100e3, DP_c_gh=4e5,
        PP_rec=0, DP_h_rec=4e5, DP_c_rec=2e5,
        eta_exp=0.9,
        SC_cd=0.1, DP_h_cond=2e5, DP_c_cond=50e3,
        P_high_bounds=P_high_bounds,
        m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
        m_dot_CS_fact_bounds=m_dot_CS_fact_bounds,
        m_dot_bounds=m_dot_bounds,
        eta_gh_disc=eta_gh_disc, PP_gh_disc=PP_gh_disc,
        eta_rec_disc=eta_rec_disc, PP_cd_disc=PP_cd_disc,
        # Poids du CAPEX dans evaluate_systems() (voir docstring) :
        #   0.0 = comportement d'origine (cohérence eta/DP seule)
        #   1.0 = coût et cohérence pèsent à peu près à égalité (défaut)
        #   >1.0 = priorité croissante au coût le plus bas
        capex_weight=1.0,
        # Poids de coût relatif par technologie d'échangeur, utilisés DANS
        # le PSO thermodynamique (system_RC_parallel). Valeur initiale =1
        # partout (NTU pondéré par Q_dot, coût $/NTU supposé identique entre
        # PCHE et Shell&Tube) -- recalibré automatiquement à CHAQUE itération
        # de cycle_design() depuis le CAPEX réel des échangeurs sizés (voir
        # calibrate_cost_weights_from_sizing, appelée dans cycle_design()).
        # Pour pré-calibrer AVANT même de lancer cycle_design (utile avec
        # ntop=1, où un seul candidat par itération donne peu de robustesse
        # au début) :
        #   points = [...]  # voir docstring de quick_calibration_points
        #   sized = quick_calibration_points(Optimizer, sizing_models, points)
        #   weights = calibrate_cost_weights_from_sizing(Optimizer, RC_list=sized)
        #   if weights: Optimizer.set_parameters(**weights)
        cost_w_gh=1.0, cost_w_rec=1.0, cost_w_cond=1.0,
    )
    
    if Optimizer.params['RC_ARCH'] == "Recomp":
        Optimizer.set_it_var(P_high=140e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, spliter_frac = 0.9, eta_gh=0.95, PP_gh=5, eta_rec_LT=0.8, eta_rec_HT=0.8, PP_cd=5, mdot_CS=450*n_MW)
    elif Optimizer.params['RC_ARCH'] == "Recomp_1_recup":
        Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, spliter_frac = 1, eta_gh=0.95, PP_gh=5, eta_rec=0.8, PP_cd=5, mdot_CS=450*n_MW)
    elif Optimizer.params['RC_ARCH'] == "REC":
        Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, eta_gh=0.95, PP_gh=5, eta_rec=0.8, PP_cd=5, mdot_CS=450*n_MW)
    elif Optimizer.params['RC_ARCH'] == "basic":
        Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, eta_gh=0.95, PP_gh=5, PP_cd=5, mdot_CS=450*n_MW)

    Optimizer.set_obj(W_dot=W_dot_obj, eta=eta_obj)

    Optimizer.set_CSource(T=T_cold, P=5e5,  fluid='Water', m_dot=450*n_MW)
    Optimizer.set_HSource(T=T_hot,      P=100e5, fluid='Water', m_dot=50.0)

    Optimizer.set_RC()
    
    #%% Composants — configuration statique (paramètres, bornes, corrélations, RUN_KWARGS)

    sizing_models = {}

    # --- Recuperator (PCHE) ---
    REC = sizing_models["Recuperator"] = PCHESizingOpt()
    REC.set_parameters(
        H_Corr={"1P": "Gnielinski", "SC": "Gnielinski", "2P": "Thome_Condensation"},
        C_Corr={"1P": "Gnielinski", "SC": "Gnielinski", "2P": "Flow_boiling"},
        H_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
        C_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
    )
    REC.RUN_KWARGS = dict(n_jobs=-1, n_particles=50, max_iter=50, patience=10)

    # --- GasHeater / Condenser (Shell&Tube) : géométrie + paramètres communs ---

    # PERFORMANCE : n_jobs=-1 ajouté ici (absent dans le script d'origine) —
    # par cohérence avec Recuperator/Expander_Axial/Expander_Radial qui
    # utilisent tous n_jobs=-1, ShellAndTubeSizingOpt.sizing() accepte très
    # probablement ce kwarg (même famille d'optimiseur PSO interne). Sans
    # lui, GasHeater/Condenser tournaient vraisemblablement en série sur un
    # seul coeur pendant que tout le reste utilisait tous les coeurs.
    # -> À VÉRIFIER : si .sizing() ne connaît pas n_jobs, retirer cette ligne.
    shell_tube_run_kwargs = dict(n_particles=100, max_iterations=50, obj='mass', print_flag=0, n_jobs=-1)

    GH = sizing_models["GasHeater"] = ShellAndTubeSizingOpt()
    GH.set_parameters(
        Shell_Side='H',
        H_Corr={"SC": "Shell_Kern_HTC", "1P": "Shell_Kern_HTC", "2P": "Shell_Kern_HTC"},
        C_Corr={"SC": "Gnielinski", "1P": "Gnielinski", "2P": "Flow_boiling"},
        H_DP={"SC": "Shell_Kern_DP", "1P": "Shell_Kern_DP", "2P": "Shell_Kern_DP"},
        C_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Gnielinski_DP"},
    )
    GH.RUN_KWARGS = shell_tube_run_kwargs

    CD = sizing_models["Condenser"] = ShellAndTubeSizingOpt()
    CD.set_parameters(
        Shell_Side='C',
        H_Corr={"SC": "Gnielinski", "1P": "Gnielinski", "2P": "Thome_Condensation"},
        C_Corr={"SC": "Shell_Kern_HTC", "1P": "Shell_Kern_HTC", "2P": "Shell_Kern_HTC"},
        H_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
        C_DP={"SC": "Shell_Kern_DP", "1P": "Shell_Kern_DP", "2P": "Shell_Kern_DP"},
    )
    CD.RUN_KWARGS = shell_tube_run_kwargs

    # --- Pump ---
    PP = sizing_models["Pump"] = RadialPumpODSizing(Optimizer.fluid)
    PP.RUN_KWARGS = dict()

    # --- Turbine : deux candidats, axial et radial ---
    TA = sizing_models["Expander_Axial"] = AxialTurbineMeanLineSizing(Optimizer.fluid)
    TA.RUN_KWARGS = dict(n_jobs=-1, n_particles=30, max_iter=50)

    TR = sizing_models["Expander_Radial"] = RadialTurbineMeanLineSizing(Optimizer.fluid)
    TR.RUN_KWARGS = dict(max_iter=3, n_jobs=-1)

    Optimizer.sizing_models = sizing_models
    
    #%%
    t0 = time.perf_counter()

    # PERFORMANCE — corrections apportées par rapport à l'appel d'origine
    # `cycle_design(ntop=5, n_particles=100, n_jobs=-1, patience=30)` :
    #
    #   1) `patience=30` avec `max_iter` par défaut (=30) : patience >= max_iter
    #      empêche TOUT arrêt anticipé (le critère `no_improve >= patience` ne
    #      peut jamais se déclencher avant la fin des 30 itérations). Le PSO
    #      tournait donc systématiquement à fond même quand il convergeait
    #      bien avant (cf. logs précédents : convergence en 10-40 itérations
    #      sur 50 dans la plupart des cas). Corrigé : max_iter=50, patience=15.
    #
    #   2) `ntop=5` : chaque position dimensionnée déclenche l'intégralité du
    #      pipeline de sizing (Recuperator PSO + GasHeater PSO + Condenser PSO
    #      + turbine). Réduit à 3 : les 2 positions les moins bonnes du top 5
    #      changent rarement le résultat final (cf. evaluate_systems, qui ne
    #      garde que la meilleure), donc ce coût est le plus souvent perdu.
    #      Remettre à 5 si la robustesse du choix final est prioritaire sur
    #      la vitesse.
    Optimizer.cycle_design(ntop=3, n_particles=100, max_iter=50, n_jobs=-1, patience=15)

    elapsed = time.perf_counter() - t0
    
    if Optimizer.best_RC is not None:
        log_cycle_result(
            log_path="co2_rc_results_log.csv",
            T_hot=T_hot, T_cold=T_cold,
            W_dot_obj=W_dot_obj, eta_obj=eta_obj,
            RC=Optimizer.best_RC, arch=Optimizer.params['RC_ARCH'],
            Optimizer=Optimizer, duration_s=round(elapsed, 1),
        )
    else:
        print("⚠️ Aucun RC valide trouvé — rien à logger.")
        
#     #%% 
    
# if __name__ == "__main__":
        
#     #!/usr/bin/env python3
#     # -*- coding: utf-8 -*-
#     """
#     plot_co2_rc_sweep_results.py
    
#     Génère DEUX scatter plots à partir du CSV de résultats du sweep
#     (co2_rc_sweep_results/co2_rc_sweep_results_log.csv, généré par
#     log_cycle_result() dans co2_rc_full_design_optimizer.py) :
    
#       1) CAPEX total (EUR2026)         vs W_dot_net
#       2) CAPEX / kW_net (EUR2026/kW)   vs W_dot_net   -- CAPEX spécifique
    
#     Dans les deux cas :
#       - COULEUR -> T_hot_C (température de source chaude)
#       - FORME   -> efficacité du second principe, càd eta_vs_carnot
#                    (= eta_achieved / eta_carnot), regroupée par paliers de
#                    --eta-bin.
    
#     NOTE UNITÉ : la colonne CAPEX_total du CSV est en EUR2026, pas en USD
#     (malgré son nom hérité de actualize_price(..., currency="USD") côté
#     sizing) -- seul le libellé des axes est corrigé ici, aucune conversion de
#     devise n'est appliquée.
    
#     W_dot_net utilisé (pour les deux graphes ET pour le calcul du CAPEX
#     spécifique) est la valeur BRUTE (W_dot_achieved_MW = W_dot_turbine -
#     W_dot_pompe, sans W_dot_pompe_aux non récupérable) remplacée par un palier
#     fixe via bucket_W_dot_net_MW() -- voir sa docstring. eta_achieved et
#     eta_vs_carnot sont recalculés à partir de cette valeur remplacée plutôt que
#     lus tels quels depuis le CSV (voir load_rows).
    
#     Pourquoi regrouper eta_vs_carnot par paliers plutôt que d'utiliser la
#     valeur brute : c'est une grandeur quasi-continue, donc l'utiliser telle
#     quelle donnerait presque une forme de marqueur par point -- illisible et
#     ça épuiserait vite les formes disponibles. On la regroupe donc par paliers
#     (0.05 par défaut, réglable via --eta-bin).
    
#     Usage :
#         python plot_co2_rc_sweep_results.py
#         python plot_co2_rc_sweep_results.py --log-path autre_dossier/log.csv --eta-bin 0.02 \\
#             --out capex_total.png --out-specific capex_specific.png
#     """
    
#     import argparse
#     import csv
#     import math
#     from pathlib import Path
    
#     import matplotlib.pyplot as plt
#     import numpy as np
#     from matplotlib.lines import Line2D
    
#     def candidate_log_paths():
#         """
#         Liste, dans l'ordre de priorité, les emplacements où le CSV de résultats
#         est susceptible de se trouver :
#           1) Répertoire de travail courant (CWD) -- c'est là que le sweep crée
#              co2_rc_sweep_results/, puisque save_root="co2_rc_sweep_results" est
#              un chemin RELATIF dans co2_rc_full_design_optimizer.py : le dossier
#              atterrit là où tu étais au moment de lancer `python ...`, pas
#              nécessairement à côté du script.
#           2) Dossier du package labothappy.machine.optimization (au cas où le
#              sweep aurait été lancé depuis ce dossier précis).
#           3) Dossier de CE script (au cas où plot_co2_rc_sweep_results.py serait
#              copié à côté du CSV).
#         """
#         candidates = [Path.cwd() / "co2_rc_sweep_results" / "co2_rc_sweep_results_log.csv"]
    
#         try:
#             import labothappy.machine.optimization as opt_pkg
#             candidates.append(
#                 Path(opt_pkg.__file__).resolve().parent
#                 / "co2_rc_sweep_results" / "co2_rc_sweep_results_log.csv"
#             )
#         except ImportError:
#             pass
    
#         candidates.append(
#             Path(__file__).resolve().parent
#             / "co2_rc_sweep_results" / "co2_rc_sweep_results_log.csv"
#         )
    
#         return candidates
    
    
#     def resolve_default_log_path():
#         """
#         Renvoie le premier chemin candidat qui existe réellement sur disque, ou
#         le premier candidat (CWD) si aucun n'existe -- dans ce dernier cas,
#         main() affichera la liste complète des chemins testés pour permettre un
#         diagnostic rapide, plutôt qu'un "fichier introuvable" sans contexte.
#         """
#         candidates = candidate_log_paths()
#         for c in candidates:
#             if c.is_file():
#                 return c
#         return candidates[0]
    
    
#     DEFAULT_LOG_PATH = resolve_default_log_path()
    
#     # Palette de formes de marqueurs matplotlib, cyclée si plus de catégories
#     # d'eta_vs_carnot (après binning) que de formes listées ici.
#     MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<', '>', 'p', '8']
    
    
#     def bucket_W_dot_net_MW(w_raw_MW):
#         """
#         Remplace ENTIÈREMENT la valeur brute de W_dot_net -- actuellement
#         W_dot_turbine - W_dot_pompe, colonne W_dot_achieved_MW du CSV, sans le
#         terme W_dot_pump_aux qui n'est pas récupérable depuis le cycle résolu
#         (voir discussion) -- par un palier fixe :
#             W_dot_net_brut <  10 MW -> 1  MW
#             W_dot_net_brut <  30 MW -> 10 MW
#             W_dot_net_brut <  50 MW -> 30 MW
#             W_dot_net_brut >= 50 MW -> 50 MW
#         Ce n'est PAS une soustraction d'un terme estimé : la valeur brute sert
#         uniquement à choisir le palier, puis est jetée au profit de ce palier.
#         """
#         if w_raw_MW < 10:
#             return 1.0
#         elif w_raw_MW < 30:
#             return 10.0
#         elif w_raw_MW < 50:
#             return 30.0
#         else:
#             return 50.0
    
    
#     def load_rows(log_path):
#         """
#         Lit le CSV et ne garde que les lignes exploitables : W_dot_achieved_MW,
#         CAPEX_total, T_hot_C, Q_gh, eta_carnot doivent être des nombres valides
#         (présents et non-NaN). Les échecs n'ont de toute façon pas de
#         CAPEX/W_dot_achieved -- ils sont dans co2_rc_sweep_fails_log.csv, pas
#         dans ce fichier de succès.
    
#         W_dot_achieved_MW (brut) est remplacé par bucket_W_dot_net_MW() -- voir
#         sa docstring -- pour tenir lieu de "W_dot_net" dans les deux plots.
#         eta_achieved et eta_vs_carnot sont ALORS RECALCULÉS à partir de ce
#         W_dot_net remplacé (et non plus lus tels quels depuis le CSV, qui les
#         avait calculés à partir du W_dot_net brut) :
#             eta_achieved_new  = (W_dot_net_bucketed * 1e6) / Q_gh
#             eta_vs_carnot_new = eta_achieved_new / eta_carnot
#         Q_gh (W) et eta_carnot viennent directement du CSV (colonnes déjà
#         présentes, écrites respectivement par _hx_effectivenesses et
#         log_cycle_result) -- pas besoin de les recalculer.
    
#         Le CAPEX spécifique (EUR2026/kW) est calculé à partir de ce même
#         W_dot_net remplacé :
#             CAPEX_specific = CAPEX_total / (W_dot_net_bucketed * 1000)
    
#         CAPEX_Recuperator (colonne écrite par log_cycle_result pour chaque
#         composant du cycle, architecture REC uniquement) est lue si présente,
#         pour calculer capex_recup_pct = CAPEX_Recuperator / CAPEX_total * 100 --
#         utilisée par plot_recuperator_capex_share_heatmaps(). Optionnelle : si
#         absente ou invalide sur une ligne, capex_recup_pct vaut NaN pour cette
#         ligne (elle reste utilisable pour les deux scatter plots, juste exclue
#         de la heatmap).
#         """
#         rows = []
#         with open(log_path, "r", newline="") as f:
#             reader = csv.DictReader(f)
#             for r in reader:
#                 try:
#                     w_raw = float(r["W_dot_achieved_MW"])
#                     capex = float(r["CAPEX_total"])
#                     t_hot = float(r["T_hot_C"])
#                     Q_gh = float(r["Q_gh"])
#                     eta_carnot = float(r["eta_carnot"])
#                 except (KeyError, ValueError, TypeError):
#                     continue
#                 if any(math.isnan(v) for v in (w_raw, capex, t_hot, Q_gh, eta_carnot)):
#                     continue
#                 if Q_gh == 0 or eta_carnot == 0:
#                     continue
    
#                 w = bucket_W_dot_net_MW(w_raw)
#                 eta_achieved = (w * 1e6) / Q_gh
#                 eta = eta_achieved / eta_carnot
#                 capex_specific = capex / (w * 1000.0)  # EUR2026 / kW_net
    
#                 try:
#                     capex_recup = float(r["CAPEX_Recuperator"])
#                     capex_recup_pct = (capex_recup / capex) * 100.0 if capex else float("nan")
#                     if math.isnan(capex_recup_pct):
#                         capex_recup_pct = float("nan")
#                 except (KeyError, ValueError, TypeError):
#                     capex_recup_pct = float("nan")
    
#                 # Optionnelles, pour plot_Phigh_vs_Thot() et
#                 # plot_turbomachinery_efficiency_vs_scale() -- NaN si absentes/
#                 # invalides sur une ligne donnée, sans affecter le reste.
#                 try:
#                     p_high_bar = float(r["P_high_Pa"]) / 1e5
#                 except (KeyError, ValueError, TypeError):
#                     p_high_bar = float("nan")
    
#                 try:
#                     eta_is_pump = float(r["eta_is_pump"])
#                 except (KeyError, ValueError, TypeError):
#                     eta_is_pump = float("nan")
    
#                 try:
#                     eta_is_expander = float(r["eta_is_expander"])
#                 except (KeyError, ValueError, TypeError):
#                     eta_is_expander = float("nan")
    
#                 rows.append({
#                     "W_dot": w, "CAPEX": capex, "CAPEX_specific": capex_specific,
#                     "eta_vs_carnot": eta, "T_hot_C": t_hot,
#                     "capex_recup_pct": capex_recup_pct,
#                     "P_high_bar": p_high_bar,
#                     "eta_is_pump": eta_is_pump,
#                     "eta_is_expander": eta_is_expander,
#                 })
    
#         return rows
    
    
#     def bin_eta(eta, bin_size):
#         """
#         Regroupe eta_vs_carnot par paliers de `bin_size` (ex. 0.05), pour
#         assigner un nombre raisonnable de formes de marqueurs distinctes plutôt
#         qu'une par valeur flottante quasi-unique.
#         """
#         return round(round(eta / bin_size) * bin_size, 6)
    
    
#     def build_color_map(values):
#         """
#         Associe une couleur distincte à chaque valeur de `values` (ici les
#         T_hot_C uniques), via un colormap continu échantillonné en autant de
#         points que de valeurs -- fonctionne quel que soit le nombre de
#         températures balayées dans le sweep (pas limité à une palette fixe,
#         contrairement aux formes de marqueurs).
#         """
#         n = max(len(values), 1)
#         try:
#             # matplotlib >= 3.7 : API non dépréciée
#             cmap = plt.colormaps["viridis"].resampled(n)
#         except AttributeError:
#             # matplotlib < 3.7 : repli sur l'ancienne API
#             cmap = plt.cm.get_cmap("viridis", n)
#         return {v: cmap(i) for i, v in enumerate(values)}
    
    
#     def _pow10_formatter(x, pos):
#         """
#         Formatte un tick d'échelle log sous la forme "1e6", "1e8", etc. --
#         plus lisible que la notation scientifique par défaut de matplotlib
#         ("$\\mathdefault{10^{6}}$") pour ce genre de graphe.
#         """
#         if x <= 0:
#             return ""
#         exponent = int(round(np.log10(x)))
#         return f"1e{exponent}"
    
    
#     def _plain_formatter(x, pos):
#         """
#         Affiche un tick sous forme de nombre simple ("1", "10", "100"...) --
#         utilisé une fois l'axe déjà exprimé en millions (y_scale=1e6), où la
#         notation "1eN" façon _pow10_formatter n'a plus lieu d'être.
#         """
#         if x <= 0:
#             return ""
#         if float(x).is_integer():
#             return f"{int(x)}"
#         return f"{x:g}"
    
    
#     def make_scatter(rows, T_hot_values, eta_bins, color_of_T, marker_of_eta,
#                       y_key, y_label, log_x, log_y, y_scale=1.0, y_ticks=None,
#                       label_fontsize=18, legend_fontsize=13):
#         """
#         Construit UNE figure scatter (y_key vs W_dot) avec les deux légendes
#         (couleur = T_hot, forme = palier d'eta_vs_carnot). Factorisé pour être
#         appelé une fois par grandeur tracée (CAPEX total, puis CAPEX
#         spécifique) sans dupliquer toute la logique de légendes.
    
#         y_scale : diviseur appliqué aux valeurs de y_key avant traçage (ex. 1e6
#         pour passer d'EUR2026 à M€2026). N'affecte que l'affichage : les données
#         sources (rows) restent inchangées.
    
#         y_ticks : liste optionnelle de valeurs -- DANS L'UNITÉ AFFICHÉE (donc
#         après division par y_scale) -- pour fixer explicitement les ticks de
#         l'axe y en échelle log, plutôt que de laisser matplotlib choisir. Sans
#         effet si log_y=False. Formatés en "1eN" si y_scale=1 (valeurs brutes,
#         grandes), ou en nombre simple ("1", "10", "100") si y_scale != 1
#         (valeurs déjà réduites, ex. en millions).
#         """
#         fig, ax = plt.subplots(figsize=(9, 6))
    
#         for t_hot in T_hot_values:
#             for eta_bin in eta_bins:
#                 subset = [r for r in rows if r["T_hot_C"] == t_hot and r["eta_bin"] == eta_bin]
#                 if not subset:
#                     continue
#                 ax.scatter(
#                     [r["W_dot"] for r in subset],
#                     [r[y_key] / y_scale for r in subset],
#                     color=color_of_T[t_hot],
#                     marker=marker_of_eta[eta_bin],
#                     edgecolors="black", linewidths=0.4,
#                     s=60, alpha=0.85,
#                 )
    
#         # Titre du graphe retiré, labels d'axes agrandis, W_dot_net en notation
#         # mathématique (accent point = dérivée temporelle, comme dans le code).
#         ax.set_xlabel(r"$\dot{W}_{net}$ (MW)", fontsize=label_fontsize)
#         ax.set_ylabel(y_label, fontsize=label_fontsize)
#         ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#         ax.grid(True, which="both", alpha=0.3)
    
#         if log_x:
#             ax.set_xscale("log")
#         if log_y:
#             ax.set_yscale("log")
#             if y_ticks is not None:
#                 ax.set_yticks(y_ticks)
#                 formatter = _plain_formatter if y_scale != 1.0 else _pow10_formatter
#                 ax.yaxis.set_major_formatter(plt.FuncFormatter(formatter))
#                 ax.yaxis.set_minor_formatter(plt.NullFormatter())
    
#         # Deux légendes séparées (une par dimension encodée), à L'INTÉRIEUR du
#         # graphe (haut-gauche / bas-droite) -- et sans titre, pour rester compact.
#         color_handles = [
#             Line2D([0], [0], marker='o', linestyle='', markerfacecolor=color_of_T[t],
#                    markeredgecolor='black', markersize=8, label=f"T_hot = {t:.0f} °C")
#             for t in T_hot_values
#         ]
#         marker_handles = [
#             Line2D([0], [0], marker=marker_of_eta[e], linestyle='', markerfacecolor='grey',
#                    markeredgecolor='black', markersize=8, label=f"η/η_carnot ≈ {e:.2f}")
#             for e in eta_bins
#         ]
    
#         legend1 = ax.legend(handles=color_handles, loc="upper left", framealpha=0.9,
#                              fontsize=legend_fontsize)
#         ax.add_artist(legend1)
#         ax.legend(handles=marker_handles, loc="lower right", framealpha=0.9,
#                   fontsize=legend_fontsize)
    
#         fig.tight_layout()
#         return fig
    
    
#     def plot_recuperator_capex_share_heatmaps(rows, label_fontsize=15, legend_fontsize=11):
#         """
#         Pour chaque palier d'efficacité de second principe (eta_bin, càd
#         η/η_carnot -- voir bin_eta), trace une heatmap T_hot vs W_dot_net où la
#         couleur code la part (%) du CAPEX du récupérateur dans le CAPEX total
#         du cycle (capex_recup_pct = CAPEX_Recuperator / CAPEX_total * 100, voir
#         load_rows).
    
#         Une figure à sous-graphes, un par palier eta_bin, avec une échelle de
#         couleur COMMUNE à tous les sous-graphes (même vmin/vmax) pour rester
#         comparable d'un palier à l'autre -- une colorbar par sous-graphe
#         donnerait sinon l'impression trompeuse que chaque palier a la même
#         plage de valeurs alors que les échelles différeraient.
    
#         Si plusieurs lignes tombent dans la même case (même T_hot, même W_dot
#         bucketé, même eta_bin -- possible puisque W_dot est déjà regroupé en 4
#         valeurs fixes), la valeur affichée est leur MOYENNE, annotée dans la
#         case avec le nombre de points moyennés entre parenthèses si > 1.
#         """
#         eta_bins = sorted({r["eta_bin"] for r in rows if not math.isnan(r["capex_recup_pct"])})
#         if not eta_bins:
#             raise ValueError(
#                 "Aucune ligne exploitable pour la heatmap (capex_recup_pct manquant partout -- "
#                 "vérifie que la colonne CAPEX_Recuperator existe dans le CSV, architecture REC)."
#             )
    
#         T_hot_values = sorted({r["T_hot_C"] for r in rows})
#         W_dot_values = sorted({r["W_dot"] for r in rows})
    
#         # Grille (eta_bin, T_hot, W_dot) -> (moyenne, effectif), construite une
#         # seule fois pour calculer aussi le vmin/vmax communs avant de tracer.
#         grids = {}
#         all_values = []
#         for eta_bin in eta_bins:
#             grid = np.full((len(T_hot_values), len(W_dot_values)), np.nan)
#             counts = np.zeros((len(T_hot_values), len(W_dot_values)), dtype=int)
#             for i, t in enumerate(T_hot_values):
#                 for j, w in enumerate(W_dot_values):
#                     vals = [
#                         r["capex_recup_pct"] for r in rows
#                         if r["eta_bin"] == eta_bin and r["T_hot_C"] == t and r["W_dot"] == w
#                         and not math.isnan(r["capex_recup_pct"])
#                     ]
#                     if vals:
#                         grid[i, j] = sum(vals) / len(vals)
#                         counts[i, j] = len(vals)
#                         all_values.append(grid[i, j])
#             grids[eta_bin] = (grid, counts)
    
#         vmin, vmax = min(all_values), max(all_values)
    
#         ncols = min(3, len(eta_bins))
#         nrows = math.ceil(len(eta_bins) / ncols)
#         fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.0 * nrows), squeeze=False)
    
#         im = None
#         for idx, eta_bin in enumerate(eta_bins):
#             ax = axes[idx // ncols][idx % ncols]
#             grid, counts = grids[eta_bin]
    
#             im = ax.imshow(grid, aspect="auto", origin="lower", cmap="viridis",
#                             vmin=vmin, vmax=vmax)
    
#             ax.set_xticks(range(len(W_dot_values)))
#             ax.set_xticklabels([f"{w:g}" for w in W_dot_values], fontsize=legend_fontsize)
#             ax.set_yticks(range(len(T_hot_values)))
#             ax.set_yticklabels([f"{t:.0f}" for t in T_hot_values], fontsize=legend_fontsize)
#             ax.set_title(f"η/η_carnot ≈ {eta_bin:.2f}", fontsize=label_fontsize * 0.75)
    
#             for i in range(len(T_hot_values)):
#                 for j in range(len(W_dot_values)):
#                     v = grid[i, j]
#                     if math.isnan(v):
#                         continue
#                     label = f"{v:.0f}%" if counts[i, j] <= 1 else f"{v:.0f}%\n(n={counts[i, j]})"
#                     color = "white" if v > (vmin + vmax) / 2 else "black"
#                     ax.text(j, i, label, ha="center", va="center",
#                             color=color, fontsize=legend_fontsize * 0.85)
    
#         # Sous-graphes inutilisés (si len(eta_bins) ne remplit pas exactement
#         # la grille nrows x ncols) masqués plutôt que laissés vides.
#         for idx in range(len(eta_bins), nrows * ncols):
#             axes[idx // ncols][idx % ncols].axis("off")
    
#         fig.supxlabel(r"$\dot{W}_{net}$ (MW)", fontsize=label_fontsize)
#         fig.supylabel(r"$T_{hot}$ (°C)", fontsize=label_fontsize)
    
#         if im is not None:
#             cbar = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
#             cbar.set_label("CAPEX récupérateur / CAPEX total (%)", fontsize=label_fontsize * 0.85)
#             cbar.ax.tick_params(labelsize=legend_fontsize)
    
#         return fig
    
    
#     def plot_eta_envelope_vs_Thot(rows, label_fontsize=18, legend_fontsize=13):
#         """
#         Enveloppe du rendement de second principe RÉELLEMENT ATTEINT (pas la
#         cible visée) en fonction de T_hot : pour chaque puissance W_dot_net
#         (palier), la meilleure valeur de eta_vs_carnot obtenue à chaque T_hot,
#         plus l'enveloppe globale (meilleure valeur tous W_dot confondus) en
#         trait épais noir en pointillés.
    
#         Contrairement aux scatter CAPEX, ce graphe utilise eta_vs_carnot comme
#         variable CONTINUE en y (pas binnée par --eta-bin) : l'objectif n'est
#         pas de distinguer plein de paliers cibles, mais de voir la limite
#         pratique réellement atteignable par le pipeline, indépendamment de ce
#         qui était visé (eta_obj).
#         """
#         T_hot_values = sorted({r["T_hot_C"] for r in rows})
#         W_dot_values = sorted({r["W_dot"] for r in rows})
#         color_of_W = build_color_map(W_dot_values)
    
#         fig, ax = plt.subplots(figsize=(9, 6))
    
#         # Nuage de tous les points en fond (contexte, sans légende dédiée) --
#         # montre la dispersion autour de l'enveloppe, pas seulement son maximum.
#         for w in W_dot_values:
#             subset = [r for r in rows if r["W_dot"] == w]
#             ax.scatter([r["T_hot_C"] for r in subset], [r["eta_vs_carnot"] for r in subset],
#                        color=color_of_W[w], alpha=0.25, s=25, edgecolors="none")
    
#         # Enveloppe par W_dot : meilleur eta_vs_carnot atteint à chaque T_hot.
#         for w in W_dot_values:
#             env = []
#             for t in T_hot_values:
#                 vals = [r["eta_vs_carnot"] for r in rows if r["T_hot_C"] == t and r["W_dot"] == w]
#                 env.append(max(vals) if vals else np.nan)
#             ax.plot(T_hot_values, env, color=color_of_W[w], marker="o", linewidth=2,
#                     markersize=7, label=f"{w:g} MW")
    
#         # Enveloppe globale, tous W_dot confondus.
#         global_env = []
#         for t in T_hot_values:
#             vals = [r["eta_vs_carnot"] for r in rows if r["T_hot_C"] == t]
#             global_env.append(max(vals) if vals else np.nan)
#         ax.plot(T_hot_values, global_env, color="black", linewidth=2.5, linestyle="--",
#                 marker="*", markersize=13, label="Enveloppe globale")
    
#         ax.set_xlabel(r"$T_{hot}$ (°C)", fontsize=label_fontsize)
#         ax.set_ylabel(r"$\eta / \eta_{carnot}$ (meilleur atteint)", fontsize=label_fontsize)
#         ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#         ax.grid(True, alpha=0.3)
#         ax.legend(fontsize=legend_fontsize, loc="best", framealpha=0.9, title="$\\dot{W}_{net}$")
    
#         fig.tight_layout()
#         return fig
    
    
#     def plot_Phigh_vs_Thot(rows, eta_bins, marker_of_eta, label_fontsize=18, legend_fontsize=13):
#         """
#         Pression haute optimale (P_high, en bar) choisie par le PSO
#         thermodynamique en fonction de T_hot -- résultat classique de la
#         littérature sCO2 (la pression haute optimale dépend fortement de
#         T_hot, typiquement croissante). Sert de test de cohérence physique du
#         pipeline : si la tendance connue n'apparaît pas dans ces données, ça
#         pointe vers un problème dans le PSO thermodynamique plutôt que dans le
#         sizing des composants.
    
#         COULEUR -> W_dot_net (palier), FORME -> palier eta_bin -- inversé par
#         rapport aux scatter CAPEX (où COULEUR = T_hot), puisque T_hot est ici
#         l'axe des x et ne peut plus servir de couleur.
#         """
#         valid = [r for r in rows if not math.isnan(r["P_high_bar"])]
#         if not valid:
#             raise ValueError(
#                 "Colonne P_high_Pa absente/inexploitable dans le CSV -- rien à tracer."
#             )
    
#         W_dot_values = sorted({r["W_dot"] for r in valid})
#         color_of_W = build_color_map(W_dot_values)
    
#         fig, ax = plt.subplots(figsize=(9, 6))
    
#         for w in W_dot_values:
#             for eta_bin in eta_bins:
#                 subset = [r for r in valid if r["W_dot"] == w and r["eta_bin"] == eta_bin]
#                 if not subset:
#                     continue
#                 ax.scatter(
#                     [r["T_hot_C"] for r in subset],
#                     [r["P_high_bar"] for r in subset],
#                     color=color_of_W[w], marker=marker_of_eta[eta_bin],
#                     edgecolors="black", linewidths=0.4, s=60, alpha=0.85,
#                 )
    
#         ax.set_xlabel(r"$T_{hot}$ (°C)", fontsize=label_fontsize)
#         ax.set_ylabel(r"$P_{high}$ (bar)", fontsize=label_fontsize)
#         ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#         ax.grid(True, alpha=0.3)
    
#         color_handles = [
#             Line2D([0], [0], marker='o', linestyle='', markerfacecolor=color_of_W[w],
#                    markeredgecolor='black', markersize=8, label=f"{w:g} MW")
#             for w in W_dot_values
#         ]
#         marker_handles = [
#             Line2D([0], [0], marker=marker_of_eta[e], linestyle='', markerfacecolor='grey',
#                    markeredgecolor='black', markersize=8, label=f"η/η_carnot ≈ {e:.2f}")
#             for e in eta_bins
#             if any(r["eta_bin"] == e for r in valid)
#         ]
    
#         legend1 = ax.legend(handles=color_handles, loc="upper left", framealpha=0.9,
#                              fontsize=legend_fontsize)
#         ax.add_artist(legend1)
#         ax.legend(handles=marker_handles, loc="lower right", framealpha=0.9,
#                   fontsize=legend_fontsize)
    
#         fig.tight_layout()
#         return fig
    
    
#     def plot_turbomachinery_efficiency_vs_scale(rows, T_hot_values, eta_bins, color_of_T,
#                                                  marker_of_eta, log_x=False,
#                                                  label_fontsize=18, legend_fontsize=13):
#         """
#         Rendement isentropique de la pompe et de la turbine (eta_is_pump,
#         eta_is_expander) en fonction de l'échelle (W_dot_net) -- pour vérifier
#         empiriquement, sur les vraies données du sweep, la dégradation de
#         rendement aux petites puissances (turbomachines miniatures) discutée
#         précédemment de façon qualitative.
    
#         Deux sous-graphes côte à côte (pompe, turbine), même style que les
#         scatter CAPEX (COULEUR = T_hot, FORME = palier eta_bin), axe y partagé
#         pour comparer directement les deux composants.
#         """
#         has_pump = any(not math.isnan(r["eta_is_pump"]) for r in rows)
#         has_exp = any(not math.isnan(r["eta_is_expander"]) for r in rows)
#         if not (has_pump or has_exp):
#             raise ValueError(
#                 "Colonnes eta_is_pump/eta_is_expander absentes/inexploitables -- rien à tracer."
#             )
    
#         fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
#         specs = [("eta_is_pump", "Pump", axes[0]), ("eta_is_expander", "Turbine", axes[1])]
    
#         for key, title, ax in specs:
#             for t_hot in T_hot_values:
#                 for eta_bin in eta_bins:
#                     subset = [
#                         r for r in rows
#                         if r["T_hot_C"] == t_hot and r["eta_bin"] == eta_bin
#                         and not math.isnan(r[key])
#                     ]
#                     if not subset:
#                         continue
#                     ax.scatter(
#                         [r["W_dot"] for r in subset],
#                         [r[key] for r in subset],
#                         color=color_of_T[t_hot], marker=marker_of_eta[eta_bin],
#                         edgecolors="black", linewidths=0.4, s=60, alpha=0.85,
#                     )
#             ax.set_xlabel(r"$\dot{W}_{net}$ (MW)", fontsize=label_fontsize)
#             ax.set_title(title, fontsize=label_fontsize)
#             ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#             ax.grid(True, which="both", alpha=0.3)
#             if log_x:
#                 ax.set_xscale("log")
    
#         axes[0].set_ylabel(r"$\eta_{is}$", fontsize=label_fontsize)
#         axes[0].set_ylim(0.8, 1.0)  # axe y partagé (sharey=True) -- s'applique aux deux panneaux
    
#         # Zones indicatives de préférence radiale/axiale sur le panneau Turbine
#         # uniquement (voir discussion : en dessous du seuil, le débit volumique
#         # CO2 est trop faible pour qu'une turbine axiale reste dans sa plage de
#         # vitesse spécifique efficace -- la radiale devient préférable ;
#         # au-dessus, c'est l'inverse). Seuil purement indicatif (15 MW), pas
#         # issu d'un calcul de vitesse spécifique sur CES données.
#         ax_turbine = axes[1]
#         xlims = ax_turbine.get_xlim()
#         boundary = 15.0
#         ax_turbine.axvspan(xlims[0], boundary, color="tab:blue", alpha=0.08, zorder=0)
#         ax_turbine.axvspan(boundary, xlims[1], color="tab:orange", alpha=0.08, zorder=0)
#         ax_turbine.axvline(boundary, color="grey", linestyle="--", linewidth=1, zorder=0)
#         ax_turbine.set_xlim(xlims)  # axvspan/axvline ne doivent pas élargir la vue
    
#         # Étiquettes texte plutôt qu'une légende séparée -- transform mixte
#         # (x en coordonnées data, y en fraction d'axe) pour rester bien placées
#         # en haut du graphe quel que soit log_x et quelles que soient les
#         # limites y (fixées à [0.8, 1.0] ci-dessus).
#         trans = ax_turbine.get_xaxis_transform()
#         if log_x:
#             x_left = math.sqrt(max(xlims[0], 1e-9) * boundary)
#             x_right = math.sqrt(boundary * xlims[1])
#         else:
#             x_left = (xlims[0] + boundary) / 2
#             x_right = (boundary + xlims[1]) / 2
#         ax_turbine.text(x_left, 0.95, f"Radial turbine\npreferred (< {boundary:g} MW)",
#                          transform=trans, ha="center", va="top",
#                          fontsize=legend_fontsize, color="tab:blue")
#         ax_turbine.text(x_right, 0.95, f"Axial turbine\npreferred (≥ {boundary:g} MW)",
#                          transform=trans, ha="center", va="top",
#                          fontsize=legend_fontsize, color="tab:orange")
    
#         color_handles = [
#             Line2D([0], [0], marker='o', linestyle='', markerfacecolor=color_of_T[t],
#                    markeredgecolor='black', markersize=8, label=f"T_hot = {t:.0f} °C")
#             for t in T_hot_values
#         ]
#         marker_handles = [
#             Line2D([0], [0], marker=marker_of_eta[e], linestyle='', markerfacecolor='grey',
#                    markeredgecolor='black', markersize=8, label=f"η/η_carnot ≈ {e:.2f}")
#             for e in eta_bins
#         ]
    
#         axes[0].legend(handles=color_handles, loc="upper left", framealpha=0.9,
#                        fontsize=legend_fontsize)
#         ax_turbine.legend(handles=marker_handles, loc="lower right", framealpha=0.9,
#                            fontsize=legend_fontsize)
    
#         fig.tight_layout()
#         return fig
    
    
#     def main():
#         parser = argparse.ArgumentParser(
#             description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
#         )
#         parser.add_argument("--log-path", default=DEFAULT_LOG_PATH,
#                              help=f"Chemin du CSV de résultats (défaut : {DEFAULT_LOG_PATH})")
#         parser.add_argument("--eta-bin", type=float, default=0.05,
#                              help="Taille des paliers pour regrouper eta_vs_carnot (défaut : 0.05)")
#         parser.add_argument("--log-x", action="store_true",
#                              help="Échelle log sur l'axe W_dot_net (les deux graphes)")
#         parser.add_argument("--linear-y", dest="log_y", action="store_false", default=True,
#                              help="Désactive l'échelle log sur l'axe y (activée par défaut, les "
#                                   "deux graphes) : le CAPEX (total ou spécifique) s'étale "
#                                   "typiquement sur plusieurs ordres de grandeur vu la plage de "
#                                   "W_dot du sweep, 1 à 100 MW")
#         parser.add_argument("--out", default=None,
#                              help="Chemin de sauvegarde du graphe CAPEX total (PNG). "
#                                   "Si omis, affichage interactif.")
#         parser.add_argument("--out-specific", default=None,
#                              help="Chemin de sauvegarde du graphe CAPEX spécifique "
#                                   "(EUR2026/kW_net) (PNG). Si omis, affichage interactif.")
#         parser.add_argument("--out-recuperator-heatmap", default=None,
#                              help="Chemin de sauvegarde de la heatmap part CAPEX récupérateur "
#                                   "(PNG, un sous-graphe par palier eta_bin). Si omis, affichage interactif.")
#         parser.add_argument("--out-eta-envelope", default=None,
#                              help="Chemin de sauvegarde de l'enveloppe eta_vs_carnot réellement "
#                                   "atteint vs T_hot (PNG). Si omis, affichage interactif.")
#         parser.add_argument("--out-phigh", default=None,
#                              help="Chemin de sauvegarde du graphe P_high vs T_hot (PNG). "
#                                   "Si omis, affichage interactif.")
#         parser.add_argument("--out-turbomachinery", default=None,
#                              help="Chemin de sauvegarde du graphe rendement pompe/turbine vs "
#                                   "échelle (PNG). Si omis, affichage interactif.")
#         args = parser.parse_args()
    
#         log_path = Path(args.log_path)
#         if not log_path.is_file():
#             tried = "\n".join(f"  - {c}" for c in candidate_log_paths())
#             raise SystemExit(
#                 f"Fichier introuvable : {log_path}\n\n"
#                 f"Emplacements testés automatiquement (aucun trouvé) :\n{tried}\n\n"
#                 f"Précise le bon chemin avec --log-path, par exemple :\n"
#                 f"  python plot_co2_rc_sweep_results.py --log-path \"chemin/vers/co2_rc_sweep_results_log.csv\""
#             )
    
#         rows = load_rows(log_path)
#         if not rows:
#             raise SystemExit(
#                 f"Aucune ligne exploitable dans {log_path} (colonnes attendues : "
#                 "W_dot_achieved_MW, CAPEX_total, T_hot_C, Q_gh, eta_carnot)."
#             )
    
#         for r in rows:
#             r["eta_bin"] = bin_eta(r["eta_vs_carnot"], args.eta_bin)
    
#         T_hot_values = sorted({r["T_hot_C"] for r in rows})
#         eta_bins = sorted({r["eta_bin"] for r in rows})
    
#         if len(eta_bins) > len(MARKERS):
#             print(f"⚠️  {len(eta_bins)} paliers d'eta_vs_carnot pour seulement {len(MARKERS)} "
#                   f"formes de marqueurs disponibles -- certains paliers partageront la même forme. "
#                   f"Augmente --eta-bin pour réduire le nombre de paliers si besoin.")
    
#         color_of_T = build_color_map(T_hot_values)
#         marker_of_eta = {e: MARKERS[i % len(MARKERS)] for i, e in enumerate(eta_bins)}
    
#         fig_total = make_scatter(
#             rows, T_hot_values, eta_bins, color_of_T, marker_of_eta,
#             y_key="CAPEX", y_label="CAPEX (M€2026)",
#             log_x=args.log_x, log_y=args.log_y,
#             y_scale=1e6, y_ticks=[1, 10, 100],
#         )
#         fig_specific = make_scatter(
#             rows, T_hot_values, eta_bins, color_of_T, marker_of_eta,
#             y_key="CAPEX_specific", y_label="CAPEX / kW_net (EUR2026/kW)",
#             log_x=args.log_x, log_y=args.log_y,
#         )
    
#         fig_heatmap = None
#         n_with_recup = sum(1 for r in rows if not math.isnan(r["capex_recup_pct"]))
#         if n_with_recup:
#             fig_heatmap = plot_recuperator_capex_share_heatmaps(rows)
#         else:
#             print("⚠️  Colonne CAPEX_Recuperator absente/inexploitable dans le CSV -- "
#                   "heatmap part récupérateur non générée.")
    
#         fig_envelope = plot_eta_envelope_vs_Thot(rows)
        
#         fig_phigh = None
#         n_with_phigh = sum(1 for r in rows if not math.isnan(r["P_high_bar"]))
#         if n_with_phigh:
#             fig_phigh = plot_Phigh_vs_Thot(rows, eta_bins, marker_of_eta)
#         else:
#             print("⚠️  Colonne P_high_Pa absente/inexploitable dans le CSV -- "
#                   "graphe P_high vs T_hot non généré.")
    
#         fig_turbo = None
#         n_with_turbo = sum(
#             1 for r in rows
#             if not math.isnan(r["eta_is_pump"]) or not math.isnan(r["eta_is_expander"])
#         )
#         if n_with_turbo:
#             fig_turbo = plot_turbomachinery_efficiency_vs_scale(
#                 rows, T_hot_values, eta_bins, color_of_T, marker_of_eta, log_x=args.log_x,
#             )
#         else:
#             print("⚠️  Colonnes eta_is_pump/eta_is_expander absentes/inexploitables dans le CSV -- "
#                   "graphe rendement pompe/turbine non généré.")
    
#         if args.out:
#             fig_total.savefig(args.out, dpi=150, bbox_inches="tight")
#             print(f"Graphe CAPEX total sauvegardé : {args.out}")
#         if args.out_specific:
#             fig_specific.savefig(args.out_specific, dpi=150, bbox_inches="tight")
#             print(f"Graphe CAPEX spécifique sauvegardé : {args.out_specific}")
#         if args.out_recuperator_heatmap and fig_heatmap is not None:
#             fig_heatmap.savefig(args.out_recuperator_heatmap, dpi=150, bbox_inches="tight")
#             print(f"Heatmap part récupérateur sauvegardée : {args.out_recuperator_heatmap}")
#         if args.out_eta_envelope:
#             fig_envelope.savefig(args.out_eta_envelope, dpi=150, bbox_inches="tight")
#             print(f"Enveloppe eta_vs_carnot sauvegardée : {args.out_eta_envelope}")
#         if args.out_phigh and fig_phigh is not None:
#             fig_phigh.savefig(args.out_phigh, dpi=150, bbox_inches="tight")
#             print(f"Graphe P_high vs T_hot sauvegardé : {args.out_phigh}")
#         if args.out_turbomachinery and fig_turbo is not None:
#             fig_turbo.savefig(args.out_turbomachinery, dpi=150, bbox_inches="tight")
#             print(f"Graphe rendement pompe/turb
    
# if __name__ == "__main__":
        
#     #!/usr/bin/env python3
#     # -*- coding: utf-8 -*-
#     """
#     plot_co2_rc_sweep_results.py
    
#     Génère DEUX scatter plots à partir du CSV de résultats du sweep
#     (co2_rc_sweep_results/co2_rc_sweep_results_log.csv, généré par
#     log_cycle_result() dans co2_rc_full_design_optimizer.py) :
    
#       1) CAPEX total (EUR2026)         vs W_dot_net
#       2) CAPEX / kW_net (EUR2026/kW)   vs W_dot_net   -- CAPEX spécifique
    
#     Dans les deux cas :
#       - COULEUR -> T_hot_C (température de source chaude)
#       - FORME   -> efficacité du second principe, càd eta_vs_carnot
#                    (= eta_achieved / eta_carnot), regroupée par paliers de
#                    --eta-bin.
    
#     NOTE UNITÉ : la colonne CAPEX_total du CSV est en EUR2026, pas en USD
#     (malgré son nom hérité de actualize_price(..., currency="USD") côté
#     sizing) -- seul le libellé des axes est corrigé ici, aucune conversion de
#     devise n'est appliquée.
    
#     W_dot_net utilisé (pour les deux graphes ET pour le calcul du CAPEX
#     spécifique) est la valeur BRUTE (W_dot_achieved_MW = W_dot_turbine -
#     W_dot_pompe, sans W_dot_pompe_aux non récupérable) remplacée par un palier
#     fixe via bucket_W_dot_net_MW() -- voir sa docstring. eta_achieved et
#     eta_vs_carnot sont recalculés à partir de cette valeur remplacée plutôt que
#     lus tels quels depuis le CSV (voir load_rows).
    
#     Pourquoi regrouper eta_vs_carnot par paliers plutôt que d'utiliser la
#     valeur brute : c'est une grandeur quasi-continue, donc l'utiliser telle
#     quelle donnerait presque une forme de marqueur par point -- illisible et
#     ça épuiserait vite les formes disponibles. On la regroupe donc par paliers
#     (0.05 par défaut, réglable via --eta-bin).
    
#     Usage :
#         python plot_co2_rc_sweep_results.py
#         python plot_co2_rc_sweep_results.py --log-path autre_dossier/log.csv --eta-bin 0.02 \\
#             --out capex_total.png --out-specific capex_specific.png
#     """
    
#     import argparse
#     import csv
#     import math
#     from pathlib import Path
    
#     import matplotlib.pyplot as plt
#     import numpy as np
#     from matplotlib.lines import Line2D
    
#     def candidate_log_paths():
#         """
#         Liste, dans l'ordre de priorité, les emplacements où le CSV de résultats
#         est susceptible de se trouver :
#           1) Répertoire de travail courant (CWD) -- c'est là que le sweep crée
#              co2_rc_sweep_results/, puisque save_root="co2_rc_sweep_results" est
#              un chemin RELATIF dans co2_rc_full_design_optimizer.py : le dossier
#              atterrit là où tu étais au moment de lancer `python ...`, pas
#              nécessairement à côté du script.
#           2) Dossier du package labothappy.machine.optimization (au cas où le
#              sweep aurait été lancé depuis ce dossier précis).
#           3) Dossier de CE script (au cas où plot_co2_rc_sweep_results.py serait
#              copié à côté du CSV).
#         """
#         candidates = [Path.cwd() / "co2_rc_sweep_results" / "co2_rc_sweep_results_log.csv"]
    
#         try:
#             import labothappy.machine.optimization as opt_pkg
#             candidates.append(
#                 Path(opt_pkg.__file__).resolve().parent
#                 / "co2_rc_sweep_results" / "co2_rc_sweep_results_log.csv"
#             )
#         except ImportError:
#             pass
    
#         candidates.append(
#             Path(__file__).resolve().parent
#             / "co2_rc_sweep_results" / "co2_rc_sweep_results_log.csv"
#         )
    
#         return candidates
    
    
#     def resolve_default_log_path():
#         """
#         Renvoie le premier chemin candidat qui existe réellement sur disque, ou
#         le premier candidat (CWD) si aucun n'existe -- dans ce dernier cas,
#         main() affichera la liste complète des chemins testés pour permettre un
#         diagnostic rapide, plutôt qu'un "fichier introuvable" sans contexte.
#         """
#         candidates = candidate_log_paths()
#         for c in candidates:
#             if c.is_file():
#                 return c
#         return candidates[0]
    
    
#     DEFAULT_LOG_PATH = resolve_default_log_path()
    
#     # Palette de formes de marqueurs matplotlib, cyclée si plus de catégories
#     # d'eta_vs_carnot (après binning) que de formes listées ici.
#     MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<', '>', 'p', '8']
    
    
#     def bucket_W_dot_net_MW(w_raw_MW):
#         """
#         Remplace ENTIÈREMENT la valeur brute de W_dot_net -- actuellement
#         W_dot_turbine - W_dot_pompe, colonne W_dot_achieved_MW du CSV, sans le
#         terme W_dot_pump_aux qui n'est pas récupérable depuis le cycle résolu
#         (voir discussion) -- par un palier fixe :
#             W_dot_net_brut <  10 MW -> 1  MW
#             W_dot_net_brut <  30 MW -> 10 MW
#             W_dot_net_brut <  50 MW -> 30 MW
#             W_dot_net_brut >= 50 MW -> 50 MW
#         Ce n'est PAS une soustraction d'un terme estimé : la valeur brute sert
#         uniquement à choisir le palier, puis est jetée au profit de ce palier.
#         """
#         if w_raw_MW < 10:
#             return 1.0
#         elif w_raw_MW < 30:
#             return 10.0
#         elif w_raw_MW < 50:
#             return 30.0
#         else:
#             return 50.0
    
    
#     def load_rows(log_path):
#         """
#         Lit le CSV et ne garde que les lignes exploitables : W_dot_achieved_MW,
#         CAPEX_total, T_hot_C, Q_gh, eta_carnot doivent être des nombres valides
#         (présents et non-NaN). Les échecs n'ont de toute façon pas de
#         CAPEX/W_dot_achieved -- ils sont dans co2_rc_sweep_fails_log.csv, pas
#         dans ce fichier de succès.
    
#         W_dot_achieved_MW (brut) est remplacé par bucket_W_dot_net_MW() -- voir
#         sa docstring -- pour tenir lieu de "W_dot_net" dans les deux plots.
#         eta_achieved et eta_vs_carnot sont ALORS RECALCULÉS à partir de ce
#         W_dot_net remplacé (et non plus lus tels quels depuis le CSV, qui les
#         avait calculés à partir du W_dot_net brut) :
#             eta_achieved_new  = (W_dot_net_bucketed * 1e6) / Q_gh
#             eta_vs_carnot_new = eta_achieved_new / eta_carnot
#         Q_gh (W) et eta_carnot viennent directement du CSV (colonnes déjà
#         présentes, écrites respectivement par _hx_effectivenesses et
#         log_cycle_result) -- pas besoin de les recalculer.
    
#         Le CAPEX spécifique (EUR2026/kW) est calculé à partir de ce même
#         W_dot_net remplacé :
#             CAPEX_specific = CAPEX_total / (W_dot_net_bucketed * 1000)
    
#         CAPEX_Recuperator (colonne écrite par log_cycle_result pour chaque
#         composant du cycle, architecture REC uniquement) est lue si présente,
#         pour calculer capex_recup_pct = CAPEX_Recuperator / CAPEX_total * 100 --
#         utilisée par plot_recuperator_capex_share_heatmaps(). Optionnelle : si
#         absente ou invalide sur une ligne, capex_recup_pct vaut NaN pour cette
#         ligne (elle reste utilisable pour les deux scatter plots, juste exclue
#         de la heatmap).
#         """
#         rows = []
#         with open(log_path, "r", newline="") as f:
#             reader = csv.DictReader(f)
#             for r in reader:
#                 try:
#                     w_raw = float(r["W_dot_achieved_MW"])
#                     capex = float(r["CAPEX_total"])
#                     t_hot = float(r["T_hot_C"])
#                     Q_gh = float(r["Q_gh"])
#                     eta_carnot = float(r["eta_carnot"])
#                 except (KeyError, ValueError, TypeError):
#                     continue
#                 if any(math.isnan(v) for v in (w_raw, capex, t_hot, Q_gh, eta_carnot)):
#                     continue
#                 if Q_gh == 0 or eta_carnot == 0:
#                     continue
    
#                 w = bucket_W_dot_net_MW(w_raw)
#                 eta_achieved = (w * 1e6) / Q_gh
#                 eta = eta_achieved / eta_carnot
#                 capex_specific = capex / (w * 1000.0)  # EUR2026 / kW_net
    
#                 try:
#                     capex_recup = float(r["CAPEX_Recuperator"])
#                     capex_recup_pct = (capex_recup / capex) * 100.0 if capex else float("nan")
#                     if math.isnan(capex_recup_pct):
#                         capex_recup_pct = float("nan")
#                 except (KeyError, ValueError, TypeError):
#                     capex_recup_pct = float("nan")
    
#                 # Optionnelles, pour plot_Phigh_vs_Thot() et
#                 # plot_turbomachinery_efficiency_vs_scale() -- NaN si absentes/
#                 # invalides sur une ligne donnée, sans affecter le reste.
#                 try:
#                     p_high_bar = float(r["P_high_Pa"]) / 1e5
#                 except (KeyError, ValueError, TypeError):
#                     p_high_bar = float("nan")
    
#                 try:
#                     eta_is_pump = float(r["eta_is_pump"])
#                 except (KeyError, ValueError, TypeError):
#                     eta_is_pump = float("nan")
    
#                 try:
#                     eta_is_expander = float(r["eta_is_expander"])
#                 except (KeyError, ValueError, TypeError):
#                     eta_is_expander = float("nan")
    
#                 rows.append({
#                     "W_dot": w, "CAPEX": capex, "CAPEX_specific": capex_specific,
#                     "eta_vs_carnot": eta, "T_hot_C": t_hot,
#                     "capex_recup_pct": capex_recup_pct,
#                     "P_high_bar": p_high_bar,
#                     "eta_is_pump": eta_is_pump,
#                     "eta_is_expander": eta_is_expander,
#                 })
    
#         return rows
    
    
#     def bin_eta(eta, bin_size):
#         """
#         Regroupe eta_vs_carnot par paliers de `bin_size` (ex. 0.05), pour
#         assigner un nombre raisonnable de formes de marqueurs distinctes plutôt
#         qu'une par valeur flottante quasi-unique.
#         """
#         return round(round(eta / bin_size) * bin_size, 6)
    
    
#     def build_color_map(values):
#         """
#         Associe une couleur distincte à chaque valeur de `values` (ici les
#         T_hot_C uniques), via un colormap continu échantillonné en autant de
#         points que de valeurs -- fonctionne quel que soit le nombre de
#         températures balayées dans le sweep (pas limité à une palette fixe,
#         contrairement aux formes de marqueurs).
#         """
#         n = max(len(values), 1)
#         try:
#             # matplotlib >= 3.7 : API non dépréciée
#             cmap = plt.colormaps["viridis"].resampled(n)
#         except AttributeError:
#             # matplotlib < 3.7 : repli sur l'ancienne API
#             cmap = plt.cm.get_cmap("viridis", n)
#         return {v: cmap(i) for i, v in enumerate(values)}
    
    
#     def _pow10_formatter(x, pos):
#         """
#         Formatte un tick d'échelle log sous la forme "1e6", "1e8", etc. --
#         plus lisible que la notation scientifique par défaut de matplotlib
#         ("$\\mathdefault{10^{6}}$") pour ce genre de graphe.
#         """
#         if x <= 0:
#             return ""
#         exponent = int(round(np.log10(x)))
#         return f"1e{exponent}"
    
    
#     def _plain_formatter(x, pos):
#         """
#         Affiche un tick sous forme de nombre simple ("1", "10", "100"...) --
#         utilisé une fois l'axe déjà exprimé en millions (y_scale=1e6), où la
#         notation "1eN" façon _pow10_formatter n'a plus lieu d'être.
#         """
#         if x <= 0:
#             return ""
#         if float(x).is_integer():
#             return f"{int(x)}"
#         return f"{x:g}"
    
    
#     def make_scatter(rows, T_hot_values, eta_bins, color_of_T, marker_of_eta,
#                       y_key, y_label, log_x, log_y, y_scale=1.0, y_ticks=None,
#                       label_fontsize=18, legend_fontsize=13):
#         """
#         Construit UNE figure scatter (y_key vs W_dot) avec les deux légendes
#         (couleur = T_hot, forme = palier d'eta_vs_carnot). Factorisé pour être
#         appelé une fois par grandeur tracée (CAPEX total, puis CAPEX
#         spécifique) sans dupliquer toute la logique de légendes.
    
#         y_scale : diviseur appliqué aux valeurs de y_key avant traçage (ex. 1e6
#         pour passer d'EUR2026 à M€2026). N'affecte que l'affichage : les données
#         sources (rows) restent inchangées.
    
#         y_ticks : liste optionnelle de valeurs -- DANS L'UNITÉ AFFICHÉE (donc
#         après division par y_scale) -- pour fixer explicitement les ticks de
#         l'axe y en échelle log, plutôt que de laisser matplotlib choisir. Sans
#         effet si log_y=False. Formatés en "1eN" si y_scale=1 (valeurs brutes,
#         grandes), ou en nombre simple ("1", "10", "100") si y_scale != 1
#         (valeurs déjà réduites, ex. en millions).
#         """
#         fig, ax = plt.subplots(figsize=(9, 6))
    
#         for t_hot in T_hot_values:
#             for eta_bin in eta_bins:
#                 subset = [r for r in rows if r["T_hot_C"] == t_hot and r["eta_bin"] == eta_bin]
#                 if not subset:
#                     continue
#                 ax.scatter(
#                     [r["W_dot"] for r in subset],
#                     [r[y_key] / y_scale for r in subset],
#                     color=color_of_T[t_hot],
#                     marker=marker_of_eta[eta_bin],
#                     edgecolors="black", linewidths=0.4,
#                     s=60, alpha=0.85,
#                 )
    
#         # Titre du graphe retiré, labels d'axes agrandis, W_dot_net en notation
#         # mathématique (accent point = dérivée temporelle, comme dans le code).
#         ax.set_xlabel(r"$\dot{W}_{net}$ (MW)", fontsize=label_fontsize)
#         ax.set_ylabel(y_label, fontsize=label_fontsize)
#         ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#         ax.grid(True, which="both", alpha=0.3)
    
#         if log_x:
#             ax.set_xscale("log")
#         if log_y:
#             ax.set_yscale("log")
#             if y_ticks is not None:
#                 ax.set_yticks(y_ticks)
#                 formatter = _plain_formatter if y_scale != 1.0 else _pow10_formatter
#                 ax.yaxis.set_major_formatter(plt.FuncFormatter(formatter))
#                 ax.yaxis.set_minor_formatter(plt.NullFormatter())
    
#         # Deux légendes séparées (une par dimension encodée), à L'INTÉRIEUR du
#         # graphe (haut-gauche / bas-droite) -- et sans titre, pour rester compact.
#         color_handles = [
#             Line2D([0], [0], marker='o', linestyle='', markerfacecolor=color_of_T[t],
#                    markeredgecolor='black', markersize=8, label=f"T_hot = {t:.0f} °C")
#             for t in T_hot_values
#         ]
#         marker_handles = [
#             Line2D([0], [0], marker=marker_of_eta[e], linestyle='', markerfacecolor='grey',
#                    markeredgecolor='black', markersize=8, label=f"η/η_carnot ≈ {e:.2f}")
#             for e in eta_bins
#         ]
    
#         legend1 = ax.legend(handles=color_handles, loc="upper left", framealpha=0.9,
#                              fontsize=legend_fontsize)
#         ax.add_artist(legend1)
#         ax.legend(handles=marker_handles, loc="lower right", framealpha=0.9,
#                   fontsize=legend_fontsize)
    
#         fig.tight_layout()
#         return fig
    
    
#     def plot_recuperator_capex_share_heatmaps(rows, label_fontsize=15, legend_fontsize=11):
#         """
#         Pour chaque palier d'efficacité de second principe (eta_bin, càd
#         η/η_carnot -- voir bin_eta), trace une heatmap T_hot vs W_dot_net où la
#         couleur code la part (%) du CAPEX du récupérateur dans le CAPEX total
#         du cycle (capex_recup_pct = CAPEX_Recuperator / CAPEX_total * 100, voir
#         load_rows).
    
#         Une figure à sous-graphes, un par palier eta_bin, avec une échelle de
#         couleur COMMUNE à tous les sous-graphes (même vmin/vmax) pour rester
#         comparable d'un palier à l'autre -- une colorbar par sous-graphe
#         donnerait sinon l'impression trompeuse que chaque palier a la même
#         plage de valeurs alors que les échelles différeraient.
    
#         Si plusieurs lignes tombent dans la même case (même T_hot, même W_dot
#         bucketé, même eta_bin -- possible puisque W_dot est déjà regroupé en 4
#         valeurs fixes), la valeur affichée est leur MOYENNE, annotée dans la
#         case avec le nombre de points moyennés entre parenthèses si > 1.
#         """
#         eta_bins = sorted({r["eta_bin"] for r in rows if not math.isnan(r["capex_recup_pct"])})
#         if not eta_bins:
#             raise ValueError(
#                 "Aucune ligne exploitable pour la heatmap (capex_recup_pct manquant partout -- "
#                 "vérifie que la colonne CAPEX_Recuperator existe dans le CSV, architecture REC)."
#             )
    
#         T_hot_values = sorted({r["T_hot_C"] for r in rows})
#         W_dot_values = sorted({r["W_dot"] for r in rows})
    
#         # Grille (eta_bin, T_hot, W_dot) -> (moyenne, effectif), construite une
#         # seule fois pour calculer aussi le vmin/vmax communs avant de tracer.
#         grids = {}
#         all_values = []
#         for eta_bin in eta_bins:
#             grid = np.full((len(T_hot_values), len(W_dot_values)), np.nan)
#             counts = np.zeros((len(T_hot_values), len(W_dot_values)), dtype=int)
#             for i, t in enumerate(T_hot_values):
#                 for j, w in enumerate(W_dot_values):
#                     vals = [
#                         r["capex_recup_pct"] for r in rows
#                         if r["eta_bin"] == eta_bin and r["T_hot_C"] == t and r["W_dot"] == w
#                         and not math.isnan(r["capex_recup_pct"])
#                     ]
#                     if vals:
#                         grid[i, j] = sum(vals) / len(vals)
#                         counts[i, j] = len(vals)
#                         all_values.append(grid[i, j])
#             grids[eta_bin] = (grid, counts)
    
#         vmin, vmax = min(all_values), max(all_values)
    
#         ncols = min(3, len(eta_bins))
#         nrows = math.ceil(len(eta_bins) / ncols)
#         fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.0 * nrows), squeeze=False)
    
#         im = None
#         for idx, eta_bin in enumerate(eta_bins):
#             ax = axes[idx // ncols][idx % ncols]
#             grid, counts = grids[eta_bin]
    
#             im = ax.imshow(grid, aspect="auto", origin="lower", cmap="viridis",
#                             vmin=vmin, vmax=vmax)
    
#             ax.set_xticks(range(len(W_dot_values)))
#             ax.set_xticklabels([f"{w:g}" for w in W_dot_values], fontsize=legend_fontsize)
#             ax.set_yticks(range(len(T_hot_values)))
#             ax.set_yticklabels([f"{t:.0f}" for t in T_hot_values], fontsize=legend_fontsize)
#             ax.set_title(f"η/η_carnot ≈ {eta_bin:.2f}", fontsize=label_fontsize * 0.75)
    
#             for i in range(len(T_hot_values)):
#                 for j in range(len(W_dot_values)):
#                     v = grid[i, j]
#                     if math.isnan(v):
#                         continue
#                     label = f"{v:.0f}%" if counts[i, j] <= 1 else f"{v:.0f}%\n(n={counts[i, j]})"
#                     color = "white" if v > (vmin + vmax) / 2 else "black"
#                     ax.text(j, i, label, ha="center", va="center",
#                             color=color, fontsize=legend_fontsize * 0.85)
    
#         # Sous-graphes inutilisés (si len(eta_bins) ne remplit pas exactement
#         # la grille nrows x ncols) masqués plutôt que laissés vides.
#         for idx in range(len(eta_bins), nrows * ncols):
#             axes[idx // ncols][idx % ncols].axis("off")
    
#         fig.supxlabel(r"$\dot{W}_{net}$ (MW)", fontsize=label_fontsize)
#         fig.supylabel(r"$T_{hot}$ (°C)", fontsize=label_fontsize)
    
#         if im is not None:
#             cbar = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
#             cbar.set_label("CAPEX récupérateur / CAPEX total (%)", fontsize=label_fontsize * 0.85)
#             cbar.ax.tick_params(labelsize=legend_fontsize)
    
#         return fig
    
    
#     def plot_eta_envelope_vs_Thot(rows, label_fontsize=18, legend_fontsize=13):
#         """
#         Enveloppe du rendement de second principe RÉELLEMENT ATTEINT (pas la
#         cible visée) en fonction de T_hot : pour chaque puissance W_dot_net
#         (palier), la meilleure valeur de eta_vs_carnot obtenue à chaque T_hot,
#         plus l'enveloppe globale (meilleure valeur tous W_dot confondus) en
#         trait épais noir en pointillés.
    
#         Contrairement aux scatter CAPEX, ce graphe utilise eta_vs_carnot comme
#         variable CONTINUE en y (pas binnée par --eta-bin) : l'objectif n'est
#         pas de distinguer plein de paliers cibles, mais de voir la limite
#         pratique réellement atteignable par le pipeline, indépendamment de ce
#         qui était visé (eta_obj).
#         """
#         T_hot_values = sorted({r["T_hot_C"] for r in rows})
#         W_dot_values = sorted({r["W_dot"] for r in rows})
#         color_of_W = build_color_map(W_dot_values)
    
#         fig, ax = plt.subplots(figsize=(9, 6))
    
#         # Nuage de tous les points en fond (contexte, sans légende dédiée) --
#         # montre la dispersion autour de l'enveloppe, pas seulement son maximum.
#         for w in W_dot_values:
#             subset = [r for r in rows if r["W_dot"] == w]
#             ax.scatter([r["T_hot_C"] for r in subset], [r["eta_vs_carnot"] for r in subset],
#                        color=color_of_W[w], alpha=0.25, s=25, edgecolors="none")
    
#         # Enveloppe par W_dot : meilleur eta_vs_carnot atteint à chaque T_hot.
#         for w in W_dot_values:
#             env = []
#             for t in T_hot_values:
#                 vals = [r["eta_vs_carnot"] for r in rows if r["T_hot_C"] == t and r["W_dot"] == w]
#                 env.append(max(vals) if vals else np.nan)
#             ax.plot(T_hot_values, env, color=color_of_W[w], marker="o", linewidth=2,
#                     markersize=7, label=f"{w:g} MW")
    
#         # Enveloppe globale, tous W_dot confondus.
#         global_env = []
#         for t in T_hot_values:
#             vals = [r["eta_vs_carnot"] for r in rows if r["T_hot_C"] == t]
#             global_env.append(max(vals) if vals else np.nan)
#         ax.plot(T_hot_values, global_env, color="black", linewidth=2.5, linestyle="--",
#                 marker="*", markersize=13, label="Enveloppe globale")
    
#         ax.set_xlabel(r"$T_{hot}$ (°C)", fontsize=label_fontsize)
#         ax.set_ylabel(r"$\eta / \eta_{carnot}$ (meilleur atteint)", fontsize=label_fontsize)
#         ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#         ax.grid(True, alpha=0.3)
#         ax.legend(fontsize=legend_fontsize, loc="best", framealpha=0.9, title="$\\dot{W}_{net}$")
    
#         fig.tight_layout()
#         return fig
    
    
#     def plot_Phigh_vs_Thot(rows, eta_bins, marker_of_eta, label_fontsize=18, legend_fontsize=13):
#         """
#         Pression haute optimale (P_high, en bar) choisie par le PSO
#         thermodynamique en fonction de T_hot -- résultat classique de la
#         littérature sCO2 (la pression haute optimale dépend fortement de
#         T_hot, typiquement croissante). Sert de test de cohérence physique du
#         pipeline : si la tendance connue n'apparaît pas dans ces données, ça
#         pointe vers un problème dans le PSO thermodynamique plutôt que dans le
#         sizing des composants.
    
#         COULEUR -> W_dot_net (palier), FORME -> palier eta_bin -- inversé par
#         rapport aux scatter CAPEX (où COULEUR = T_hot), puisque T_hot est ici
#         l'axe des x et ne peut plus servir de couleur.
#         """
#         valid = [r for r in rows if not math.isnan(r["P_high_bar"])]
#         if not valid:
#             raise ValueError(
#                 "Colonne P_high_Pa absente/inexploitable dans le CSV -- rien à tracer."
#             )
    
#         W_dot_values = sorted({r["W_dot"] for r in valid})
#         color_of_W = build_color_map(W_dot_values)
    
#         fig, ax = plt.subplots(figsize=(9, 6))
    
#         for w in W_dot_values:
#             for eta_bin in eta_bins:
#                 subset = [r for r in valid if r["W_dot"] == w and r["eta_bin"] == eta_bin]
#                 if not subset:
#                     continue
#                 ax.scatter(
#                     [r["T_hot_C"] for r in subset],
#                     [r["P_high_bar"] for r in subset],
#                     color=color_of_W[w], marker=marker_of_eta[eta_bin],
#                     edgecolors="black", linewidths=0.4, s=60, alpha=0.85,
#                 )
    
#         ax.set_xlabel(r"$T_{hot}$ (°C)", fontsize=label_fontsize)
#         ax.set_ylabel(r"$P_{high}$ (bar)", fontsize=label_fontsize)
#         ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#         ax.grid(True, alpha=0.3)
    
#         color_handles = [
#             Line2D([0], [0], marker='o', linestyle='', markerfacecolor=color_of_W[w],
#                    markeredgecolor='black', markersize=8, label=f"{w:g} MW")
#             for w in W_dot_values
#         ]
#         marker_handles = [
#             Line2D([0], [0], marker=marker_of_eta[e], linestyle='', markerfacecolor='grey',
#                    markeredgecolor='black', markersize=8, label=f"η/η_carnot ≈ {e:.2f}")
#             for e in eta_bins
#             if any(r["eta_bin"] == e for r in valid)
#         ]
    
#         legend1 = ax.legend(handles=color_handles, loc="upper left", framealpha=0.9,
#                              fontsize=legend_fontsize)
#         ax.add_artist(legend1)
#         ax.legend(handles=marker_handles, loc="lower right", framealpha=0.9,
#                   fontsize=legend_fontsize)
    
#         fig.tight_layout()
#         return fig
    
    
#     def plot_turbomachinery_efficiency_vs_scale(rows, T_hot_values, eta_bins, color_of_T,
#                                                  marker_of_eta, log_x=False,
#                                                  label_fontsize=18, legend_fontsize=13):
#         """
#         Rendement isentropique de la pompe et de la turbine (eta_is_pump,
#         eta_is_expander) en fonction de l'échelle (W_dot_net) -- pour vérifier
#         empiriquement, sur les vraies données du sweep, la dégradation de
#         rendement aux petites puissances (turbomachines miniatures) discutée
#         précédemment de façon qualitative.
    
#         Deux sous-graphes côte à côte (pompe, turbine), même style que les
#         scatter CAPEX (COULEUR = T_hot, FORME = palier eta_bin), axe y partagé
#         pour comparer directement les deux composants.
#         """
#         has_pump = any(not math.isnan(r["eta_is_pump"]) for r in rows)
#         has_exp = any(not math.isnan(r["eta_is_expander"]) for r in rows)
#         if not (has_pump or has_exp):
#             raise ValueError(
#                 "Colonnes eta_is_pump/eta_is_expander absentes/inexploitables -- rien à tracer."
#             )
    
#         fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
#         specs = [("eta_is_pump", "Pump", axes[0]), ("eta_is_expander", "Turbine", axes[1])]
    
#         for key, title, ax in specs:
#             for t_hot in T_hot_values:
#                 for eta_bin in eta_bins:
#                     subset = [
#                         r for r in rows
#                         if r["T_hot_C"] == t_hot and r["eta_bin"] == eta_bin
#                         and not math.isnan(r[key])
#                     ]
#                     if not subset:
#                         continue
#                     ax.scatter(
#                         [r["W_dot"] for r in subset],
#                         [r[key] for r in subset],
#                         color=color_of_T[t_hot], marker=marker_of_eta[eta_bin],
#                         edgecolors="black", linewidths=0.4, s=60, alpha=0.85,
#                     )
#             ax.set_xlabel(r"$\dot{W}_{net}$ (MW)", fontsize=label_fontsize)
#             ax.set_title(title, fontsize=label_fontsize)
#             ax.tick_params(axis="both", labelsize=label_fontsize * 0.75)
#             ax.grid(True, which="both", alpha=0.3)
#             if log_x:
#                 ax.set_xscale("log")
    
#         axes[0].set_ylabel(r"$\eta_{is}$", fontsize=label_fontsize)
#         axes[0].set_ylim(0.8, 1.0)  # axe y partagé (sharey=True) -- s'applique aux deux panneaux
    
#         # Zones indicatives de préférence radiale/axiale sur le panneau Turbine
#         # uniquement (voir discussion : en dessous du seuil, le débit volumique
#         # CO2 est trop faible pour qu'une turbine axiale reste dans sa plage de
#         # vitesse spécifique efficace -- la radiale devient préférable ;
#         # au-dessus, c'est l'inverse). Seuil purement indicatif (15 MW), pas
#         # issu d'un calcul de vitesse spécifique sur CES données.
#         ax_turbine = axes[1]
#         xlims = ax_turbine.get_xlim()
#         boundary = 15.0
#         ax_turbine.axvspan(xlims[0], boundary, color="tab:blue", alpha=0.08, zorder=0)
#         ax_turbine.axvspan(boundary, xlims[1], color="tab:orange", alpha=0.08, zorder=0)
#         ax_turbine.axvline(boundary, color="grey", linestyle="--", linewidth=1, zorder=0)
#         ax_turbine.set_xlim(xlims)  # axvspan/axvline ne doivent pas élargir la vue
    
#         # Étiquettes texte plutôt qu'une légende séparée -- transform mixte
#         # (x en coordonnées data, y en fraction d'axe) pour rester bien placées
#         # en haut du graphe quel que soit log_x et quelles que soient les
#         # limites y (fixées à [0.8, 1.0] ci-dessus).
#         trans = ax_turbine.get_xaxis_transform()
#         if log_x:
#             x_left = math.sqrt(max(xlims[0], 1e-9) * boundary)
#             x_right = math.sqrt(boundary * xlims[1])
#         else:
#             x_left = (xlims[0] + boundary) / 2
#             x_right = (boundary + xlims[1]) / 2
#         ax_turbine.text(x_left, 0.95, f"Radial turbine\npreferred (< {boundary:g} MW)",
#                          transform=trans, ha="center", va="top",
#                          fontsize=legend_fontsize, color="tab:blue")
#         ax_turbine.text(x_right, 0.95, f"Axial turbine\npreferred (≥ {boundary:g} MW)",
#                          transform=trans, ha="center", va="top",
#                          fontsize=legend_fontsize, color="tab:orange")
    
#         color_handles = [
#             Line2D([0], [0], marker='o', linestyle='', markerfacecolor=color_of_T[t],
#                    markeredgecolor='black', markersize=8, label=f"T_hot = {t:.0f} °C")
#             for t in T_hot_values
#         ]
#         marker_handles = [
#             Line2D([0], [0], marker=marker_of_eta[e], linestyle='', markerfacecolor='grey',
#                    markeredgecolor='black', markersize=8, label=f"η/η_carnot ≈ {e:.2f}")
#             for e in eta_bins
#         ]
    
#         axes[0].legend(handles=color_handles, loc="upper left", framealpha=0.9,
#                        fontsize=legend_fontsize)
#         ax_turbine.legend(handles=marker_handles, loc="lower right", framealpha=0.9,
#                            fontsize=legend_fontsize)
    
#         fig.tight_layout()
#         return fig
    
    
#     def main():
#         parser = argparse.ArgumentParser(
#             description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
#         )
#         parser.add_argument("--log-path", default=DEFAULT_LOG_PATH,
#                              help=f"Chemin du CSV de résultats (défaut : {DEFAULT_LOG_PATH})")
#         parser.add_argument("--eta-bin", type=float, default=0.05,
#                              help="Taille des paliers pour regrouper eta_vs_carnot (défaut : 0.05)")
#         parser.add_argument("--log-x", action="store_true",
#                              help="Échelle log sur l'axe W_dot_net (les deux graphes)")
#         parser.add_argument("--linear-y", dest="log_y", action="store_false", default=True,
#                              help="Désactive l'échelle log sur l'axe y (activée par défaut, les "
#                                   "deux graphes) : le CAPEX (total ou spécifique) s'étale "
#                                   "typiquement sur plusieurs ordres de grandeur vu la plage de "
#                                   "W_dot du sweep, 1 à 100 MW")
#         parser.add_argument("--out", default=None,
#                              help="Chemin de sauvegarde du graphe CAPEX total (PNG). "
#                                   "Si omis, affichage interactif.")
#         parser.add_argument("--out-specific", default=None,
#                              help="Chemin de sauvegarde du graphe CAPEX spécifique "
#                                   "(EUR2026/kW_net) (PNG). Si omis, affichage interactif.")
#         parser.add_argument("--out-recuperator-heatmap", default=None,
#                              help="Chemin de sauvegarde de la heatmap part CAPEX récupérateur "
#                                   "(PNG, un sous-graphe par palier eta_bin). Si omis, affichage interactif.")
#         parser.add_argument("--out-eta-envelope", default=None,
#                              help="Chemin de sauvegarde de l'enveloppe eta_vs_carnot réellement "
#                                   "atteint vs T_hot (PNG). Si omis, affichage interactif.")
#         parser.add_argument("--out-phigh", default=None,
#                              help="Chemin de sauvegarde du graphe P_high vs T_hot (PNG). "
#                                   "Si omis, affichage interactif.")
#         parser.add_argument("--out-turbomachinery", default=None,
#                              help="Chemin de sauvegarde du graphe rendement pompe/turbine vs "
#                                   "échelle (PNG). Si omis, affichage interactif.")
#         args = parser.parse_args()
    
#         log_path = Path(args.log_path)
#         if not log_path.is_file():
#             tried = "\n".join(f"  - {c}" for c in candidate_log_paths())
#             raise SystemExit(
#                 f"Fichier introuvable : {log_path}\n\n"
#                 f"Emplacements testés automatiquement (aucun trouvé) :\n{tried}\n\n"
#                 f"Précise le bon chemin avec --log-path, par exemple :\n"
#                 f"  python plot_co2_rc_sweep_results.py --log-path \"chemin/vers/co2_rc_sweep_results_log.csv\""
#             )
    
#         rows = load_rows(log_path)
#         if not rows:
#             raise SystemExit(
#                 f"Aucune ligne exploitable dans {log_path} (colonnes attendues : "
#                 "W_dot_achieved_MW, CAPEX_total, T_hot_C, Q_gh, eta_carnot)."
#             )
    
#         for r in rows:
#             r["eta_bin"] = bin_eta(r["eta_vs_carnot"], args.eta_bin)
    
#         T_hot_values = sorted({r["T_hot_C"] for r in rows})
#         eta_bins = sorted({r["eta_bin"] for r in rows})
    
#         if len(eta_bins) > len(MARKERS):
#             print(f"⚠️  {len(eta_bins)} paliers d'eta_vs_carnot pour seulement {len(MARKERS)} "
#                   f"formes de marqueurs disponibles -- certains paliers partageront la même forme. "
#                   f"Augmente --eta-bin pour réduire le nombre de paliers si besoin.")
    
#         color_of_T = build_color_map(T_hot_values)
#         marker_of_eta = {e: MARKERS[i % len(MARKERS)] for i, e in enumerate(eta_bins)}
    
#         fig_total = make_scatter(
#             rows, T_hot_values, eta_bins, color_of_T, marker_of_eta,
#             y_key="CAPEX", y_label="CAPEX (M€2026)",
#             log_x=args.log_x, log_y=args.log_y,
#             y_scale=1e6, y_ticks=[1, 10, 100],
#         )
#         fig_specific = make_scatter(
#             rows, T_hot_values, eta_bins, color_of_T, marker_of_eta,
#             y_key="CAPEX_specific", y_label="CAPEX / kW_net (EUR2026/kW)",
#             log_x=args.log_x, log_y=args.log_y,
#         )
    
#         fig_heatmap = None
#         n_with_recup = sum(1 for r in rows if not math.isnan(r["capex_recup_pct"]))
#         if n_with_recup:
#             fig_heatmap = plot_recuperator_capex_share_heatmaps(rows)
#         else:
#             print("⚠️  Colonne CAPEX_Recuperator absente/inexploitable dans le CSV -- "
#                   "heatmap part récupérateur non générée.")
    
#         fig_envelope = plot_eta_envelope_vs_Thot(rows)
        
#         fig_phigh = None
#         n_with_phigh = sum(1 for r in rows if not math.isnan(r["P_high_bar"]))
#         if n_with_phigh:
#             fig_phigh = plot_Phigh_vs_Thot(rows, eta_bins, marker_of_eta)
#         else:
#             print("⚠️  Colonne P_high_Pa absente/inexploitable dans le CSV -- "
#                   "graphe P_high vs T_hot non généré.")
    
#         fig_turbo = None
#         n_with_turbo = sum(
#             1 for r in rows
#             if not math.isnan(r["eta_is_pump"]) or not math.isnan(r["eta_is_expander"])
#         )
#         if n_with_turbo:
#             fig_turbo = plot_turbomachinery_efficiency_vs_scale(
#                 rows, T_hot_values, eta_bins, color_of_T, marker_of_eta, log_x=args.log_x,
#             )
#         else:
#             print("⚠️  Colonnes eta_is_pump/eta_is_expander absentes/inexploitables dans le CSV -- "
#                   "graphe rendement pompe/turbine non généré.")
    
#         if args.out:
#             fig_total.savefig(args.out, dpi=150, bbox_inches="tight")
#             print(f"Graphe CAPEX total sauvegardé : {args.out}")
#         if args.out_specific:
#             fig_specific.savefig(args.out_specific, dpi=150, bbox_inches="tight")
#             print(f"Graphe CAPEX spécifique sauvegardé : {args.out_specific}")
#         if args.out_recuperator_heatmap and fig_heatmap is not None:
#             fig_heatmap.savefig(args.out_recuperator_heatmap, dpi=150, bbox_inches="tight")
#             print(f"Heatmap part récupérateur sauvegardée : {args.out_recuperator_heatmap}")
#         if args.out_eta_envelope:
#             fig_envelope.savefig(args.out_eta_envelope, dpi=150, bbox_inches="tight")
#             print(f"Enveloppe eta_vs_carnot sauvegardée : {args.out_eta_envelope}")
#         if args.out_phigh and fig_phigh is not None:
#             fig_phigh.savefig(args.out_phigh, dpi=150, bbox_inches="tight")
#             print(f"Graphe P_high vs T_hot sauvegardé : {args.out_phigh}")
#         if args.out_turbomachinery and fig_turbo is not None:
#             fig_turbo.savefig(args.out_turbomachinery, dpi=150, bbox_inches="tight")
#             print(f"Graphe rendement pompe/turbine sauvegardé : {args.out_turbomachinery}")
    
#         all_out_args = [args.out, args.out_specific, args.out_eta_envelope]
#         all_figs = [(fig_heatmap, args.out_recuperator_heatmap),
#                     (fig_phigh, args.out_phigh), (fig_turbo, args.out_turbomachinery)]
#         nothing_saved_for_existing_fig = any(fig is not None and not out for fig, out in all_figs)
#         if not all(all_out_args) or nothing_saved_for_existing_fig:
#             plt.show()
    
    
#     if __name__ == "__main__":
#         main()ine sauvegardé : {args.out_turbomachinery}")
    
#         all_out_args = [args.out, args.out_specific, args.out_eta_envelope]
#         all_figs = [(fig_heatmap, args.out_recuperator_heatmap),
#                     (fig_phigh, args.out_phigh), (fig_turbo, args.out_turbomachinery)]
#         nothing_saved_for_existing_fig = any(fig is not None and not out for fig, out in all_figs)
#         if not all(all_out_args) or nothing_saved_for_existing_fig:
#             plt.show()
    
    
#     if __name__ == "__main__":
#         main()