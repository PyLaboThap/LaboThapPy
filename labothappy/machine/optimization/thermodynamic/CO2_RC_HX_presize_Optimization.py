# -*- coding: utf-8 -*-
"""
co2_rc_pso_optimizer.py
CO2 Transcritical Rankine Cycle Optimizer — brique d'optimisation (PSO)
- Parallel PSO evaluation via joblib
- Support de 4 architectures : 'basic', 'REC', 'Recomp', 'Recomp_1_recup'
- Warm start

Ce module est destiné à être importé (ex: par co2_rc_full_design_optimizer.py)
qui ajoute le dimensionnement des composants par-dessus.
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

    P_sat_CS    = PropsSI('P', 'T', cs_props['T'], 'Q', 0.5, fluid)
    P_crit      = PropsSI('PCRIT', fluid)
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
        W_pump_aux = params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * params.get('eta_pp', 0.8))

        W_dot_net = W_exp - W_pump - W_pump_aux - W_cp
        eta       = W_dot_net / Q_gh if Q_gh > 0 else 0.0

        penalty_W_dot = 0
        if abs(obj['W_dot'] - W_dot_net) / obj['W_dot'] > 2e-2:
            penalty_W_dot = abs(obj['W_dot'] - W_dot_net) / obj['W_dot']

        penalty_eta = 0
        if abs(obj['eta'] - eta) / obj['eta'] > 2e-2:
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

            eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
            eta_rec = np.clip(eta_rec, 0.0, 1.0 - eps)
            eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)

            objective = (Q_gh*(-np.log(1-eta_gh)) + Q_rec*(-np.log(1-eta_rec)) + Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta)*PF
            cost = objective + penalty

        elif arch == 'basic':
            Q_cond = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon

            Q_gh = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh = RC.components['GasHeater'].model.epsilon

            eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
            eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)

            objective = (Q_gh*(-np.log(1-eta_gh)) + Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_gh)

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

            eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
            eta_rec_LT = np.clip(eta_rec_LT, 0.0, 1.0 - eps)
            eta_rec_HT = np.clip(eta_rec_HT, 0.0, 1.0 - eps)
            eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)

            objective = (Q_gh*(-np.log(1-eta_gh)) + Q_rec_LT*(-np.log(1-eta_rec_LT)) + Q_rec_HT*(-np.log(1-eta_rec_HT)) + Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec_LT + Q_rec_HT + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta)*PF
            cost = objective + penalty

        elif arch == "Recomp_1_recup":
            Q_cond = RC.components['Condenser'].model.Q.Q_dot
            RC.components['Condenser'].model.equivalent_effectiveness()
            eta_cond = RC.components['Condenser'].model.epsilon

            Q_rec_LT = RC.components['RecupLT'].model.Q.Q_dot
            eta_rec_LT = RC.components['RecupLT'].model.epsilon

            Q_gh = RC.components['GasHeater'].model.Q.Q_dot
            eta_gh = RC.components['GasHeater'].model.epsilon

            eta_gh     = np.clip(eta_gh,     0.0, 1.0 - eps)
            eta_rec_LT = np.clip(eta_rec_LT, 0.0, 1.0 - eps)
            eta_cond   = np.clip(eta_cond,   0.0, 1.0 - eps)

            objective = (Q_gh*(-np.log(1-eta_gh)) + Q_rec_LT*(-np.log(1-eta_rec_LT)) + Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec_LT + Q_gh)

            PF = 1000
            penalty = (penalty_W_dot + penalty_eta)*PF
            cost = objective + penalty

    except Exception:
        return 20000.0, np.inf, np.nan

    return cost, penalty, eta

#%% Optimizer Class (base)

class CO2RC_HX_optimizer:

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
        self.penalty_log = {}
        self.allowable_positions = []
        self.top_positions = []

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
        """Traduit un vecteur PSO x en dict de variables selon l'architecture."""
        arch = self.params.get('RC_ARCH', 'REC')
        x = np.asarray(x, dtype=float)

        if arch == "Recomp":
            P_high, m_dot, m_dot_HS_fact, spliter_frac, eta_gh, PP_gh, eta_rec_LT, eta_rec_HT, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        spliter_frac=spliter_frac, eta_gh=eta_gh, PP_gh=PP_gh,
                        eta_rec_LT=eta_rec_LT, eta_rec_HT=eta_rec_HT, PP_cd=PP_cd,
                        mdot_CS=m_dot * m_dot_CS_fact)

        elif arch == "Recomp_1_recup":
            P_high, m_dot, m_dot_HS_fact, spliter_frac, eta_gh, PP_gh, eta_rec, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        spliter_frac=spliter_frac, eta_gh=eta_gh, PP_gh=PP_gh,
                        eta_rec=eta_rec, PP_cd=PP_cd, mdot_CS=m_dot * m_dot_CS_fact)

        elif arch == "REC":
            P_high, m_dot, m_dot_HS_fact, eta_gh, PP_gh, eta_rec, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        eta_gh=eta_gh, PP_gh=PP_gh, eta_rec=eta_rec, PP_cd=PP_cd,
                        mdot_CS=m_dot * m_dot_CS_fact)

        elif arch == "basic":
            P_high, m_dot, m_dot_HS_fact, eta_gh, PP_gh, PP_cd, m_dot_CS_fact = x
            return dict(P_high=P_high, mdot=m_dot, mdot_HS=m_dot * m_dot_HS_fact,
                        eta_gh=eta_gh, PP_gh=PP_gh, PP_cd=PP_cd,
                        mdot_CS=m_dot * m_dot_CS_fact)

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

    # ------------------------------------------------------------------ final eval

    def _evaluate_final(self, best_pos):
        """Re-evaluates the best PSO position with full diagnostics."""
        unpacked = self._unpack_position(best_pos)
        self.it_var.update(unpacked)

        self._HSource_props['m_dot'] = unpacked['mdot_HS']
        self._CSource_props['m_dot'] = unpacked['mdot_CS']

        self.set_RC()
        RC = self.RC

        try:
            RC.solve()
        except Exception as e:
            self._log_penalty(f"Final solve exception: {e}")
            self.eta = self.W_dot_net = None
            return

        if not getattr(RC, 'converged', True):
            self._log_penalty("Final solve did not converge")
            self.eta = self.W_dot_net = None
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
            self.eta = self.W_dot_net = None
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
        W_pump_aux = self.params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * self.params.get('eta_pp', 0.8))

        self.W_dot_net = W_exp - W_pump - W_pump_aux - W_cp
        self.eta       = self.W_dot_net / Q_gh if Q_gh > 0 else 0.0

        self.Q_dot_waste = RC.components['GasHeater'].model.ex_H.m_dot * (
            RC.components['GasHeater'].model.ex_H.h
            - PropsSI('H', 'T', 273.15 + 15, 'P',
                      RC.components['GasHeater'].model.ex_H.p,
                      RC.components['GasHeater'].model.ex_H.fluid)
        )

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

    n_cores = multiprocessing.cpu_count()

    # ---- sweep ----
    T_vec = np.linspace(150, 150, 1) + 273.15

    n_MW = 10 # W
    W_dot_obj = n_MW*1e6 # W
    
    eta_obj = 0.12
    
    # Create optimizer instance
    Optimizer = CO2RC_HX_optimizer('CO2')
    
    # Sweep parameters
    m_dot_HS_fact_bounds = [0.1,3]
    m_dot_CS_fact_bounds = [5,15]
    P_high_bounds = np.array([90, 200]) * 1e5
    m_dot_bounds = np.array([10,80])*n_MW
    spliter_frac_bounds = np.array([0.01,0.99])

    # Discrete Variable choices
    eta_gh_disc = np.arange(0.8,0.98,0.02)
    PP_gh_disc = np.arange(1,10,1)
    eta_rec_disc= np.arange(0.6,0.96,0.02)
    eta_rec_HT_disc= np.arange(0.6,0.96,0.02)
    PP_cd_disc = np.arange(1,10,1)

    for T in T_vec:

        Optimizer.set_parameters(            
            RC_ARCH= 'REC', # 'basic', 'REC', 'Recomp_1_recup', 'Recomp'
            
            # Pump
            eta_pp=0.85,
            
            # Compressor (for recompression layouts)
            eta_cp = 0.8,
            
            # GasHeater
            DP_h_gh = 1*1e5,
            DP_c_gh = 4*1e5,
    
            # Recuperator
            PP_rec=0,
            DP_h_rec = 4*1e5,
            DP_c_rec = 2*1e5,
            
            # Expander
            eta_exp=0.93,
            
            # Condenser
            SC_cd=0.1,
            DP_h_cond = 2*1e5,
            DP_c_cond = 100*1e3,
            
            # Bounds
            P_high_bounds=P_high_bounds,
            m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
            m_dot_CS_fact_bounds=m_dot_CS_fact_bounds,
            m_dot_bounds = m_dot_bounds,
            spliter_frac_bounds = spliter_frac_bounds, 
            
            # Discrete Values
            eta_gh_disc=eta_gh_disc,
            PP_gh_disc=PP_gh_disc,
            eta_rec_disc=eta_rec_disc,
            eta_rec_HT_disc=eta_rec_HT_disc,
            PP_cd_disc=PP_cd_disc,
        )
        
        if Optimizer.params['RC_ARCH'] == "Recomp":
            Optimizer.set_it_var(P_high=140e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, spliter_frac = 0.9, eta_gh=0.95, PP_gh=5, eta_rec_LT=0.8, eta_rec_HT=0.8, PP_cd=5, mdot_CS=200*n_MW)
        elif Optimizer.params['RC_ARCH'] == "Recomp_1_recup":
            Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, spliter_frac = 1, eta_gh=0.95, PP_gh=5, eta_rec=0.8, PP_cd=5, mdot_CS=200*n_MW)
        elif Optimizer.params['RC_ARCH'] == "REC":
            Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, eta_gh=0.95, PP_gh=5, eta_rec=0.8, PP_cd=5, mdot_CS=200*n_MW)
        elif Optimizer.params['RC_ARCH'] == "basic":
            Optimizer.set_it_var(P_high=100e5, mdot=20.0*n_MW, mdot_HS=15.0*n_MW, eta_gh=0.95, PP_gh=5, PP_cd=5, mdot_CS=200*n_MW)

        Optimizer.set_obj(W_dot=W_dot_obj, eta=eta_obj)

        Optimizer.set_CSource(T=15 + 273.15, P=5e5,  fluid='Water', m_dot=1000.0)
        Optimizer.set_HSource(T=T,           P=100e5, fluid='Water', m_dot=50.0)

        Optimizer.set_RC()
        Optimizer.opt_RC(n_jobs = n_cores - 1, n_particles=50, max_iter=50, patience = 20)
        # Optimizer.opt_RC(n_jobs = 1, n_particles=100, max_iter=50, patience = 10)

    Optimizer.RC.plot_cycle_Ts()

