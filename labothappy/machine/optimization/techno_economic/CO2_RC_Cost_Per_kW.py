#%%

# -*- coding: utf-8 -*-
"""
co2_rc_full_design_optimizer.py

Étend CO2RC_HX_optimizer avec le dimensionnement des composants + CAPEX
+ boucle cycle_design + log CSV.

Architectures supportées (RC_ARCH) :
    'basic'           : pas de récupérateur
    'REC'             : Recuperator
    'Recomp_1_recup'  : RecupLT + Compressor (recompresseur)
    'Recomp'          : RecupLT + RecupHT + Compressor

Les objets de sizing sont créés UNE fois dans le __main__ (paramètres, bornes,
corrélations + attribut RUN_KWARGS). size_all_components() boucle dessus et
n'injecte à chaque appel que les entrées dépendant du point de fonctionnement.

Corrections par rapport à la version précédente :
- clés de récupérateurs dépendantes de l'architecture (RecupLT / RecupHT / Recuperator)
  dans TOUS les registres, evaluate_systems, log, calibration ;
- ajout du compresseur (wrapper RadialCompressorSizing : échec explicite via
  `penalty`, W_dot, CAPEX, reset de L_z) ;
- W_dot_net = W_turbine - W_pompe - W_compresseur dans le log ;
- eta_cp intégré à la boucle de cohérence ;
- affichage de l'exception complète (traceback) quand un sizing échoue ;
- les composants absents de RC.components sont signalés (plus d'échec muet).
"""

#%% Imports

import csv
import os
import time
import traceback
from datetime import datetime

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.sizing.turbomachinery.turbine.axial.sizing_1D.mean_line_axial_turbine_loss_model_sizing import AxialTurbineMeanLineSizing
from labothappy.sizing.turbomachinery.turbine.radial.mean_line_radial_turbine_loss_model_sizing import RadialTurbineMeanLineSizing
from labothappy.sizing.turbomachinery.compressor.radial.sizing_1D.mean_line_radial_compressor_sizing import RadialCPMLDesign
from labothappy.sizing.heat_exchanger.shell_and_tube.shell_and_tube_sizing import ShellAndTubeSizingOpt
from labothappy.sizing.heat_exchanger.PCHE.PCHE_sizing import PCHESizingOpt
from labothappy.sizing.turbomachinery.pump.radial.radial_pump_0D_sizing import RadialPumpODSizing

from labothappy.machine.optimization.thermodynamic.CO2_RC_HX_presize_Optimization import CO2RC_HX_optimizer

import warnings
warnings.filterwarnings('ignore')

#%% Clés de composants selon l'architecture

# Récupérateurs présents dans RC.components pour chaque architecture.
REC_KEYS_BY_ARCH = {
    'basic': [],
    'REC': ['Recuperator'],
    'Recomp_1_recup': ['RecupLT'],
    'Recomp': ['RecupLT', 'RecupHT'],
}

# Architectures avec recompresseur.
COMPRESSOR_ARCHS = ('Recomp_1_recup', 'Recomp')


def rec_keys(arch):
    return REC_KEYS_BY_ARCH.get(arch, [])


def _rec_dp(RC, arch):
    """(DP_h, DP_c) de référence pour les récupérateurs : max sur les récupérateurs présents."""
    dh, dc = [], []
    for k in rec_keys(arch):
        if k in RC.components:
            hx = RC.components[k].sizing.HX
            dh.append(hx.DP_h)
            dc.append(hx.DP_c)
    if not dh:
        return None, None
    return max(dh), max(dc)

#%% Extraction des entrées dynamiques (dépendent du point de fonctionnement courant)

def _hx_inputs(model):
    return dict(
        fluid_H=model.su_H.fluid, T_su_H=model.su_H.T, P_su_H=model.su_H.p, m_dot_H=model.su_H.m_dot,
        fluid_C=model.su_C.fluid, T_su_C=model.su_C.T, P_su_C=model.su_C.p, m_dot_C=model.su_C.m_dot,
    )

def _pump_inputs(model):
    return dict(P_su=model.su.p, P_ex=model.ex.p, T_su=model.su.T,
                H1=0, H2=0, v1=0, v2=0, m_dot=model.su.m_dot)

def _compressor_inputs(model):
    return dict(mdot=model.su.m_dot, p0_su=model.su.p, T0_su=model.su.T, p_ex=model.ex.p)

def _turbine_inputs(model):
    return dict(mdot=model.su.m_dot, W_dot=model.W.W_dot,
                p0_su=model.su.p, T0_su=model.su.T, p_ex=model.ex.p)

# 'Expander' est traité à part (choix axial/radial).
DYNAMIC_INPUT_EXTRACTORS = {
    'Recuperator': _hx_inputs,
    'RecupLT': _hx_inputs,
    'RecupHT': _hx_inputs,
    'GasHeater': _hx_inputs,
    'Condenser': _hx_inputs,
    'Pump': _pump_inputs,
    'Compressor': _compressor_inputs,
}

# Plancher DP_h/DP_c (0 = pas de max()).
HX_DP_FLOOR = {'Recuperator': 0.0, 'RecupLT': 0.0, 'RecupHT': 0.0,
               'GasHeater': 1e3, 'Condenser': 1e4}

# GasHeater/Condenser ont besoin de T_max_cycle/p_max_cycle (PCHE : non).
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
                f"(Pump non convergé — p_max_cycle={p_max_cycle}, T_max_cycle={T_max_cycle})"
            )
        sizing_obj.set_parameters(T_max_cycle=T_max_cycle, p_max_cycle=p_max_cycle)


def size_all_components(RC, sizing_models, turb_choice="None"):
    """
    Boucle sur `sizing_models`. Retourne (ok, results, turb_choice), results = {key: sizing_obj}.
    En cas d'échec, ok=False.
    """
    results = {}

    for key, sizing_obj in sizing_models.items():
        if key.startswith('Expander'):
            continue

        if key not in RC.components:
            print(f"⚠️ '{key}' est dans sizing_models mais absent de RC.components "
                  f"({list(RC.components.keys())}) -- ignoré.")
            continue

        model = RC.components[key].model
        RC.components[key].sizing = sizing_obj

        try:
            sizing_obj.set_inputs(**DYNAMIC_INPUT_EXTRACTORS[key](model))

            if key in HX_DP_FLOOR:
                _set_dynamic_hx_constraints(sizing_obj, model, RC, key)

            sizing_obj.sizing(**sizing_obj.RUN_KWARGS)

        except Exception as e:
            print(f"⚠️ Failed to design {key}: {type(e).__name__}: {e}")
            traceback.print_exc()
            if hasattr(model, 'su_H'):
                model.su_H.print_resume()
                model.su_C.print_resume()
                print(f"Q_dot_cstr : {model.Q.Q_dot}")
                print(f"DP_h_cstr : {model.DP_h}")
                print(f"DP_c_cstr : {model.DP_c}")
            return False, results, "Fail"

        results[key] = sizing_obj

    # --- Turbine : axial vs radial ---
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
            print(f"⚠️ Failed to design the axial Turbine: {type(e).__name__}: {e}")

    if turb_choice != 'Axial':
        try:
            Turb_radial_sizing = sizing_models['Expander_Radial']
            Turb_radial_sizing.set_inputs(**turb_inputs)
            Turb_radial_sizing.sizing(**Turb_radial_sizing.RUN_KWARGS)
            eta_radial = Turb_radial_sizing.eta_is
        except Exception as e:
            print(f"⚠️ Failed to design the radial Turbine: {type(e).__name__}: {e}")

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
    """Epsilon et Q_dot des échangeurs (NaN si indisponible)."""
    out = {}

    def safe_epsilon(component_key, label):
        try:
            out[label] = RC.components[component_key].model.epsilon
        except Exception:
            out[label] = float("nan")

    def safe_Q(component_key, label):
        try:
            out[label] = RC.components[component_key].model.Q.Q_dot
        except Exception:
            out[label] = float("nan")

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

    return out


def log_cycle_result(log_path, T_hot, T_cold, W_dot_obj, eta_obj, RC, arch,
                     Optimizer=None, duration_s=None, run_id=None):
    """
    Ajoute une ligne au CSV : CAPEX, puissance nette, efficacité, efficacités d'échangeurs.
      W_dot_net = W_expander - W_pump - W_compressor (si présent)
      eta       = W_dot_net / Q_dot_GasHeater
    """

    def safe_get(fn, default=float("nan")):
        try:
            return fn()
        except Exception:
            return default

    W_dot_exp = safe_get(lambda: RC.components['Expander'].sizing.W_dot)
    W_dot_pp = safe_get(lambda: RC.components['Pump'].sizing.W_dot)
    has_comp = 'Compressor' in RC.components and getattr(RC.components['Compressor'], 'sizing', None) is not None
    W_dot_cp = safe_get(lambda: RC.components['Compressor'].sizing.W_dot) if has_comp else 0.0
    Q_dot_gh = safe_get(lambda: RC.components['GasHeater'].sizing.best_particle.Q)

    W_dot_net = safe_get(lambda: W_dot_exp - W_dot_pp - W_dot_cp)
    eta_achieved = safe_get(lambda: W_dot_net / Q_dot_gh if Q_dot_gh else float("nan"))

    eta_carnot = safe_get(lambda: 1 - T_cold / T_hot)
    eta_vs_carnot = safe_get(lambda: eta_achieved / eta_carnot if eta_carnot else float("nan"))

    capex_total = RC.CAPEX.get("Total", float("nan"))
    capex_specific = safe_get(lambda: capex_total / (W_dot_net / 1e3) if W_dot_net else float("nan"))  # /kW

    turb_choice = safe_get(lambda: RC.components['Expander'].sizing.__class__.__name__)

    dp_h_rec, dp_c_rec = safe_get(lambda: _rec_dp(RC, arch), (None, None))

    row = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "run_id": run_id,
        "duration_s": duration_s,
        "arch": arch,

        "T_hot_C": T_hot - 273.15,
        "T_cold_C": T_cold - 273.15,
        "W_dot_obj_MW": W_dot_obj / 1e6,
        "eta_obj": eta_obj,

        "P_high_Pa": safe_get(lambda: RC.it_var.get('P_high')),
        "mdot_kg_s": safe_get(lambda: RC.it_var.get('mdot')),
        "mdot_HS_kg_s": safe_get(lambda: RC.it_var.get('mdot_HS')),
        "mdot_CS_kg_s": safe_get(lambda: RC.it_var.get('mdot_CS')),

        "W_dot_exp_MW": safe_get(lambda: W_dot_exp / 1e6),
        "W_dot_pump_MW": safe_get(lambda: W_dot_pp / 1e6),
        "W_dot_comp_MW": safe_get(lambda: W_dot_cp / 1e6),
        "W_dot_achieved_MW": safe_get(lambda: W_dot_net / 1e6),
        "eta_achieved": eta_achieved,
        "eta_carnot": eta_carnot,
        "eta_vs_carnot": eta_vs_carnot,

        "CAPEX_total": capex_total,
        "CAPEX_specific_USD_per_kW": capex_specific,

        "turbine_type": turb_choice,
        "eta_is_pump": safe_get(lambda: RC.components['Pump'].sizing.eta_is),
        "eta_is_expander": safe_get(lambda: RC.components['Expander'].sizing.eta_is),
        "eta_is_compressor": safe_get(lambda: RC.components['Compressor'].sizing.eta_is) if has_comp else float("nan"),

        "DP_h_rec_Pa": dp_h_rec,
        "DP_c_rec_Pa": dp_c_rec,
        "DP_h_gh_Pa": safe_get(lambda: RC.components['GasHeater'].sizing.best_particle.DP_h),
        "DP_c_gh_Pa": safe_get(lambda: RC.components['GasHeater'].sizing.best_particle.DP_c),
        "DP_h_cond_Pa": safe_get(lambda: RC.components['Condenser'].sizing.best_particle.DP_h),
        "DP_c_cond_Pa": safe_get(lambda: RC.components['Condenser'].sizing.best_particle.DP_c),

        "n_iter_cycle_design": safe_get(lambda: Optimizer.criterion) if Optimizer is not None else None,
    }

    row.update(_hx_effectivenesses(RC, arch))

    for key, value in RC.CAPEX.items():
        if key != "Total":
            row[f"CAPEX_{key}"] = value

    file_exists = os.path.isfile(log_path)
    with open(log_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

#%% Calibration des poids de coût

def calibrate_cost_weights_from_sizing(Optimizer, arch='REC', RC_list=None, reference='cond'):
    """
    Calibre cost_w_gh / cost_w_rec / cost_w_cond depuis le CAPEX réel des échangeurs
    sizés : cost_w_x = CAPEX_x / (Q_x * NTU_x), normalisé par `reference`.

    Supporté : 'REC' (Recuperator) et 'Recomp_1_recup' (RecupLT).
    (Pour 'Recomp' à deux récupérateurs, non supporté : ambigu.)
    Pour Recomp_1_recup, cost_w_rec n'a d'effet que si le PSO thermo l'utilise pour cette architecture.

    Retourne un dict {'cost_w_gh','cost_w_rec','cost_w_cond'} ou None.
    """
    if arch not in ('REC', 'Recomp_1_recup'):
        print(f"[calibrate_cost_weights_from_sizing] Architecture '{arch}' non supportée.")
        return None

    if RC_list is None:
        RC_list = Optimizer.potential_RC

    key_map = {'gh': 'GasHeater', 'rec': rec_keys(arch)[0], 'cond': 'Condenser'}
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
            print(f"[calibrate_cost_weights_from_sizing] k_{short} = {ks[short]:.4g} (n={len(vals)}){flag}")
        else:
            print(f"[calibrate_cost_weights_from_sizing] Aucun point valide pour '{short}' "
                  f"-- cost_w_{short} laissé inchangé.")

    if reference not in ks:
        print(f"[calibrate_cost_weights_from_sizing] Référence '{reference}' non calibrable -- abandon.")
        return None

    k_ref = ks[reference]
    weights = {f'cost_w_{short}': v / k_ref for short, v in ks.items()}
    print(f"[calibrate_cost_weights_from_sizing] cost_w_* (réf={reference}) : {weights}")
    return weights


def quick_calibration_points(Optimizer, sizing_models, points):
    """
    Size directement quelques positions choisies à la main (sans PSO) pour obtenir de vrais
    CAPEX. `points` : liste de dicts de variables (P_high, mdot, mdot_HS, mdot_CS, eta_gh, ...).
    Retourne la liste des cycles sizés avec succès.
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

#%% Classe étendue

class CO2RCOptimizer(CO2RC_HX_optimizer):
    """
    Étend CO2RC_HX_optimizer avec sizing des composants, CAPEX et boucle cycle_design.
    L'optimisation PSO (system_RC_parallel, set_RC, _evaluate_final, opt_RC) est héritée.
    """

    def __init__(self, fluid):
        super().__init__(fluid)
        self.CAPEX = {}
        self.turb_choice = "None"
        self.potential_RC = []
        self.best_RC = None
        self.sizing_models = {}
        # Snapshot FIXE des DP_* avant optimisation (référence du score de cohérence).
        self.original_DP_params = {}
        # Meilleur design (CAPEX le plus bas) vu sur tout le run, restauré à la fin.
        self.best_RC_overall = None
        self.best_capex_overall = float('inf')
        self.capex_history = []

    @property
    def arch(self):
        return self.params.get('RC_ARCH', 'REC')

    def evaluate_systems(self):
        """
        Score de chaque candidat sizé :
          1) Cohérence : écart eta_exp/eta_pp/(eta_cp) assumé vs réel ; DP réel vs DP de
             référence figé (pénalisé seulement s'il DÉPASSE la référence).
          2) CAPEX normalisé par le CAPEX minimum du lot, pondéré par `capex_weight`.
        """
        arch = self.arch
        has_comp = arch in COMPRESSOR_ARCHS
        has_rec = bool(rec_keys(arch))

        RC_scores = []
        delta_dicts = []
        capex_totals = []

        dp_ref = self.original_DP_params if self.original_DP_params else self.params

        def dp_delta(real, ref):
            return ((np.max([real, ref]) - ref) / ref) ** 2

        for RC in self.potential_RC:
            delta_dicts.append({})
            d = delta_dicts[-1]

            eta_exp = RC.components['Expander'].sizing.eta_is
            eta_pp = RC.components['Pump'].sizing.eta_is

            DP_h_gh = RC.components['GasHeater'].sizing.best_particle.DP_h
            DP_c_gh = RC.components['GasHeater'].sizing.best_particle.DP_c
            DP_h_cond = RC.components['Condenser'].sizing.best_particle.DP_h
            DP_c_cond = RC.components['Condenser'].sizing.best_particle.DP_c

            d['eta_exp'] = ((eta_exp - self.params['eta_exp']) / self.params['eta_exp']) ** 2
            d['eta_pp'] = ((eta_pp - self.params['eta_pp']) / self.params['eta_pp']) ** 2
            d['DP_h_gh'] = dp_delta(DP_h_gh, dp_ref['DP_h_gh'])
            d['DP_c_gh'] = dp_delta(DP_c_gh, dp_ref['DP_c_gh'])
            d['DP_h_cond'] = dp_delta(DP_h_cond, dp_ref['DP_h_cond'])
            d['DP_c_cond'] = dp_delta(DP_c_cond, dp_ref['DP_c_cond'])

            if has_comp and 'Compressor' in RC.components:
                eta_cp = RC.components['Compressor'].sizing.eta_is
                d['eta_cp'] = ((eta_cp - self.params['eta_cp']) / self.params['eta_cp']) ** 2

            if has_rec:
                DP_h_rec, DP_c_rec = _rec_dp(RC, arch)
                d['DP_h_rec'] = dp_delta(DP_h_rec, dp_ref['DP_h_rec'])
                d['DP_c_rec'] = dp_delta(DP_c_rec, dp_ref['DP_c_rec'])

            RC_scores.append(sum(d.values()))
            capex_totals.append(RC.CAPEX.get("Total", float("nan")))

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
        }
        if has_comp and 'Compressor' in best_RC.components:
            new_params_dict['eta_cp'] = np.round(best_RC.components['Compressor'].sizing.eta_is, 3)
        if has_rec:
            dh, dc = _rec_dp(best_RC, arch)
            new_params_dict['DP_h_rec'] = np.round(dh)
            new_params_dict['DP_c_rec'] = np.round(dc)

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

            unpacked = self._unpack_position(allowable_position['x'])
            self.it_var.update(unpacked)

            self._HSource_props['m_dot'] = unpacked['mdot_HS']
            self._CSource_props['m_dot'] = unpacked['mdot_CS']

            try:
                self.set_RC()
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

        arch = self.arch
        self.criterion = 0
        n_it_max = 10
        it = 0

        if not self.original_DP_params:
            self.original_DP_params = {
                'DP_h_gh': self.params['DP_h_gh'],
                'DP_c_gh': self.params['DP_c_gh'],
                'DP_h_cond': self.params['DP_h_cond'],
                'DP_c_cond': self.params['DP_c_cond'],
                'DP_h_rec': self.params['DP_h_rec'],
                'DP_c_rec': self.params['DP_c_rec'],
            }
            print(f"[cycle_design] DP de référence (fixes, pour tout le run) : {self.original_DP_params}")

        def _ntu(eps):
            if eps is None or not np.isfinite(eps):
                return float('nan')
            eps_c = min(max(eps, 0.0), 1.0 - 1e-6)
            return -np.log(1.0 - eps_c)

        def _f(v, fmt):
            try:
                return format(v, fmt)
            except Exception:
                return "n/a"

        while self.criterion == 0 and it < n_it_max:

            self.allowable_positions = []

            # Warm start : l'optimum bouge peu d'une itération à l'autre.
            pso_optimizer = self.opt_RC(n_jobs=n_jobs, n_particles=n_particles, max_iter=max_iter,
                                        patience=patience, ntop=ntop, init_pos=init_pos)
            init_pos = getattr(getattr(pso_optimizer, 'swarm', None), 'best_pos', init_pos)

            self.size_components()

            if not self.potential_RC:
                print("⚠️ Aucun candidat dimensionné à cette itération -- arrêt de cycle_design.")
                break

            cost_weights = calibrate_cost_weights_from_sizing(self, arch)
            if cost_weights is not None:
                self.set_parameters(**cost_weights)
                print(f"[cycle_design] cost_w_* mis à jour : {cost_weights}")
            else:
                print("[cycle_design] Calibration cost_w_* impossible ce coup-ci -- valeurs précédentes conservées.")

            new_params, best_score, delta_dict = self.evaluate_systems()
            self.new_params = new_params
            self.delta_dict = delta_dict

            current_capex = self.best_RC.CAPEX.get('Total', float('inf'))
            self.capex_history.append((it + 1, current_capex))
            if current_capex < self.best_capex_overall:
                self.best_RC_overall = self.best_RC
                self.best_capex_overall = current_capex

            hx_eff = _hx_effectivenesses(self.best_RC, arch)

            print("\n" + "=" * 60)
            print(f"  cycle_design -- itération {it + 1}")
            print("=" * 60)
            print(f"  Best Score (cohérence + CAPEX) : {best_score:.6g}")
            print(f"  CAPEX Total (best_RC)          : {self.best_RC.CAPEX.get('Total'):,.0f}")
            print(f"  CAPEX candidats sizés           : {[rc.CAPEX.get('Total') for rc in self.potential_RC]}")

            print("-" * 60)
            print("  Cohérence turbomachines (assumé -> réel, delta %) :")
            print(f"    eta_exp : {self.params['eta_exp']:.3f} -> {new_params['eta_exp']:.3f}  "
                  f"(delta={delta_dict['eta_exp']*100:.4f}%)")
            print(f"    eta_pp  : {self.params['eta_pp']:.3f} -> {new_params['eta_pp']:.3f}  "
                  f"(delta={delta_dict['eta_pp']*100:.4f}%)")
            if 'eta_cp' in new_params:
                print(f"    eta_cp  : {self.params['eta_cp']:.3f} -> {new_params['eta_cp']:.3f}  "
                      f"(delta={delta_dict.get('eta_cp', float('nan'))*100:.4f}%)")

            print("-" * 60)
            print("  DP réel (best_RC) vs référence fixe (pré-optimisation) :")
            for dp_key in ('DP_h_gh', 'DP_c_gh', 'DP_h_cond', 'DP_c_cond', 'DP_h_rec', 'DP_c_rec'):
                real_val = new_params.get(dp_key)
                ref_val = self.original_DP_params.get(dp_key)
                if real_val is not None and ref_val:
                    ratio = real_val / ref_val
                    tag = "OK (sous la référence)" if real_val <= ref_val else "⚠️ AU-DESSUS de la référence"
                    print(f"    {dp_key:10s} : réel={real_val:10.0f} Pa  |  "
                          f"référence={ref_val:10.0f} Pa  |  ratio={ratio:.2f}  |  {tag}")

            print("-" * 60)
            print("  Échangeurs -- variable cible du PSO vs effectivité réelle (best_RC) :")
            print(f"    GasHeater   : cible eta_gh={self.it_var.get('eta_gh')}, PP_gh={self.it_var.get('PP_gh')} K  |  "
                  f"eps réel={_f(hx_eff.get('eps_gh'), '.4f')}  |  NTU réel={_f(_ntu(hx_eff.get('eps_gh')), '.3f')}  |  "
                  f"Q={_f(hx_eff.get('Q_gh'), ',.0f')} W")

            if arch == 'REC':
                print(f"    Recuperator : cible eta_rec={self.it_var.get('eta_rec')}  |  "
                      f"eps réel={_f(hx_eff.get('eps_rec'), '.4f')}  |  NTU réel={_f(_ntu(hx_eff.get('eps_rec')), '.3f')}  |  "
                      f"Q={_f(hx_eff.get('Q_rec'), ',.0f')} W")
            elif arch == 'Recomp':
                for lab, k in (('LT', 'eps_rec_LT'), ('HT', 'eps_rec_HT')):
                    print(f"    Recup{lab}     : cible eta_rec_{lab}={self.it_var.get('eta_rec_' + lab)}  |  "
                          f"eps réel={_f(hx_eff.get(k), '.4f')}  |  NTU réel={_f(_ntu(hx_eff.get(k)), '.3f')}  |  "
                          f"Q={_f(hx_eff.get('Q_rec_' + lab), ',.0f')} W")
            elif arch == 'Recomp_1_recup':
                print(f"    RecupLT     : cible eta_rec={self.it_var.get('eta_rec')}  |  "
                      f"eps réel={_f(hx_eff.get('eps_rec_LT'), '.4f')}  |  NTU réel={_f(_ntu(hx_eff.get('eps_rec_LT')), '.3f')}  |  "
                      f"Q={_f(hx_eff.get('Q_rec_LT'), ',.0f')} W")

            print(f"    Condenser   : cible PP_cd={self.it_var.get('PP_cd')} K  |  "
                  f"eps réel={_f(hx_eff.get('eps_cond'), '.4f')}  |  NTU réel={_f(_ntu(hx_eff.get('eps_cond')), '.3f')}  |  "
                  f"Q={_f(hx_eff.get('Q_cond'), ',.0f')} W")
            print("=" * 60)

            updates = dict(
                eta_exp=new_params['eta_exp'],
                eta_pp=new_params['eta_pp'],
                DP_h_gh=(new_params['DP_h_gh'] + self.params['DP_h_gh']) / 2,
                DP_c_gh=(new_params['DP_c_gh'] + self.params['DP_c_gh']) / 2,
                DP_h_cond=(new_params['DP_h_cond'] + self.params['DP_h_cond']) / 2,
                DP_c_cond=(new_params['DP_c_cond'] + self.params['DP_c_cond']) / 2,
            )
            if 'eta_cp' in new_params:
                updates['eta_cp'] = new_params['eta_cp']
            if 'DP_h_rec' in new_params:
                updates['DP_h_rec'] = (new_params['DP_h_rec'] + self.params['DP_h_rec']) / 2
                updates['DP_c_rec'] = (new_params['DP_c_rec'] + self.params['DP_c_rec']) / 2
            self.set_parameters(**updates)

            self.criterion = 1
            for key in self.delta_dict:
                if self.delta_dict[key] > 1e-3:
                    self.criterion = 0
                    break
            it += 1

        # Le design retourné n'est jamais pire que le meilleur CAPEX vu pendant le run.
        if self.best_RC_overall is not None:
            self.best_RC = self.best_RC_overall

        if self.params.get('save_file_path') is not None and self.best_RC is not None:
            import json

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
                sizing_obj = getattr(self.best_RC.components[component], 'sizing', None)
                if sizing_obj is None or not hasattr(sizing_obj, 'export_params_dict'):
                    continue
                try:
                    data = sizing_obj.export_params_dict()
                    with open(os.path.join(save_folder, f"{component}.json"), "w") as f:
                        json.dump(data, f, indent=4, cls=NumpyEncoder)
                except Exception as e:
                    print(f"⚠️ Export JSON de {component} impossible : {e}")

        return self

#%% Main

if __name__ == "__main__":

    fluid = 'CO2'

    arch = "Recomp"  # 'basic', 'REC', 'Recomp_1_recup', 'Recomp'
    
    # T_hot_vec = np.array([150,200,250,300,350]) + 273.15
    T_hot_vec = np.array([350]) + 273.15
    # n_MW_vec = np.array([1,10,30,50]) 
    n_MW_vec = np.array([30]) 
    eta_obj_carnot = np.array([0.5]) 
    
    n_cases = len(T_hot_vec)*len(n_MW_vec)*len(eta_obj_carnot)
    
    case = 0
    
    for T_hot in T_hot_vec:
        for n_MW in n_MW_vec:
            for eta_obj_car in eta_obj_carnot:
                
                case += 1
                
                print("="*50)
                print(f"T_hot : {T_hot} | n_MW : {n_MW} | eta_obj_car : {eta_obj_car}")
                print(f"Cas {case} / {n_cases}")
                print("="*50)

                T_cold = 10 + 273.15
                W_dot_obj = n_MW * 1e6
                eta_obj = 0.29 # eta_obj_car*(1 - T_cold/T_hot)
                
                Optimizer = CO2RCOptimizer(fluid)
            
                m_dot_HS_fact_bounds = [0.05, 3]
                m_dot_CS_fact_bounds = [5, 30]
                P_high_bounds = np.array([100, 180]) * 1e5
                m_dot_bounds = np.array([5, 80]) * n_MW
                # Borne basse 0.01 => recompresseur à débit quasi nul, sizing quasi toujours en échec.
                spliter_frac_bounds = np.array([0.01, 0.99])
            
                eta_gh_disc = np.arange(0.8, 0.98, 0.02)
                PP_gh_disc = np.arange(1, 10, 1)
                eta_rec_disc = np.arange(0.6, 0.96, 0.02)
                PP_cd_disc = np.arange(1, 10, 1)
            
                Optimizer.set_parameters(
                    save_file_path=None,
                    RC_ARCH=arch,
            
                    eta_pp=0.85,
                    eta_cp=0.85,
                    eta_pp_aux=0.8,
            
                    DP_h_gh=50e3, DP_c_gh=50e3,
                    DP_h_rec=50e3, DP_c_rec=50e3,
                    DP_h_cond=50e3, DP_c_cond=50e3,
            
                    PP_rec=0,
                    eta_exp=0.94,
                    SC_cd=0.1,
            
                    P_high_bounds=P_high_bounds,
                    m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
                    m_dot_CS_fact_bounds=m_dot_CS_fact_bounds,
                    m_dot_bounds=m_dot_bounds,
                    spliter_frac_bounds=spliter_frac_bounds,
            
                    eta_gh_disc=eta_gh_disc, PP_gh_disc=PP_gh_disc,
                    eta_rec_disc=eta_rec_disc, PP_cd_disc=PP_cd_disc,
            
                    # 0.0 = cohérence seule ; 1.0 = coût et cohérence à égalité ; >1 = priorité au coût
                    capex_weight=1.0,
                    cost_w_gh=1.0, cost_w_rec=1.0, cost_w_cond=1.0,
                )
                
                if arch == "Recomp":
                    Optimizer.set_it_var(P_high=140e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, spliter_frac=0.9, eta_gh=0.95, PP_gh=5,
                                         eta_rec_LT=0.8, eta_rec_HT=0.8, PP_cd=5, mdot_CS=450*n_MW)
                elif arch == "Recomp_1_recup":
                    Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, spliter_frac=1, eta_gh=0.95, PP_gh=5,
                                         eta_rec=0.8, PP_cd=5, mdot_CS=450*n_MW)
                elif arch == "REC":
                    Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, eta_gh=0.95, PP_gh=5,
                                         eta_rec=0.8, PP_cd=5, mdot_CS=450*n_MW)
                elif arch == "basic":
                    Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, eta_gh=0.95, PP_gh=5,
                                         PP_cd=5, mdot_CS=450*n_MW)
            
                Optimizer.set_obj(W_dot=W_dot_obj, eta=eta_obj)
                Optimizer.set_CSource(T=T_cold, P=5e5, fluid='Water', m_dot=450*n_MW)
                Optimizer.set_HSource(T=T_hot, P=10e5, fluid='INCOMP::TVP1', m_dot=50.0*n_MW)
                
                Optimizer.set_RC()
                
                # Vérification des noms de composants (les clés de sizing_models doivent y figurer)
                print("Composants du cycle :", list(Optimizer.RC.components.keys()))
            
                #%% Composants -- configuration statique
            
                sizing_models = {}
                
                # n_jobs=-1 : à retirer si ShellAndTubeSizingOpt.sizing() ne connaît pas ce kwarg.
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
            
                PP = sizing_models["Pump"] = RadialPumpODSizing(Optimizer.fluid)
                PP.RUN_KWARGS = dict()
            
                TA = sizing_models["Expander_Axial"] = AxialTurbineMeanLineSizing(Optimizer.fluid)
                TA.RUN_KWARGS = dict(n_jobs=-1, n_particles=50, max_iter=50)
            
                TR = sizing_models["Expander_Radial"] = RadialTurbineMeanLineSizing(Optimizer.fluid)
                TR.RUN_KWARGS = dict(max_iter=10, n_jobs=-1, patience = 5)
            
                # --- Récupérateur(s) PCHE selon l'architecture ---
                def make_pche():
                    r = PCHESizingOpt()
                    r.set_parameters(
                        H_Corr={"1P": "Gnielinski", "SC": "Gnielinski", "2P": "Thome_Condensation"},
                        C_Corr={"1P": "Gnielinski", "SC": "Gnielinski", "2P": "Flow_boiling"},
                        H_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
                        C_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
                    )
                    r.RUN_KWARGS = dict(n_jobs=-1, n_particles=50, max_iter=50, patience=10)
                    return r
            
                for k in rec_keys(arch):
                    sizing_models[k] = make_pche()
            
                # --- Recompresseur ---
                if arch in COMPRESSOR_ARCHS:
                    # cost_fn : f(sizing_obj) -> coût. À fournir (même forme que la corrélation de la pompe,
                    # ex. loi en puissance de COMP.W_dot). None => CAPEX compresseur = 0 (avertissement affiché).
                    COMP = sizing_models["Compressor"] = RadialCPMLDesign(Optimizer.fluid)
                    COMP.set_parameters(
                        t_b=0.762e-3,
                        eps_imp=0.254e-3,
                        eps_bf_imp=0.254e-3,
                        k_imp=0.01e-3,
                    )
                    COMP.set_bounds(
                        Omega_bounds=[1000, 200000],   # Omega optimisé par le PSO
                        b2_r2_bounds=[0.02, 0.1],
                    )
                    COMP.RUN_KWARGS = dict(n_jobs=-1, n_particles=50, max_iter=50, patience=10)
            
                Optimizer.sizing_models = sizing_models
            
                #%% Lancement 
                # try: 
                t0 = time.perf_counter()
            
                # patience < max_iter, sinon aucun arrêt anticipé possible ; ntop réduit à 3 pour la vitesse.
                Optimizer.cycle_design(ntop=1, n_particles=100, max_iter=50, n_jobs=-1, patience=15)
            
                elapsed = time.perf_counter() - t0
            
                # except:
                #     Optimizer.best_RC = None
                
                if Optimizer.best_RC is not None:
                    log_cycle_result(
                        log_path="co2_rc_results_log.csv",
                        T_hot=T_hot, T_cold=T_cold,
                        W_dot_obj=W_dot_obj, eta_obj=eta_obj,
                        RC=Optimizer.best_RC, arch=arch,
                        Optimizer=Optimizer, duration_s=round(elapsed, 1),
                    )
                else:
                    print("⚠️ Aucun RC valide trouvé — rien à logger.")
        