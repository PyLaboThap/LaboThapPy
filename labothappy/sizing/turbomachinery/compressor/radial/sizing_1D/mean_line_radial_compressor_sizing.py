# --- loading libraries 

from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize, brentq, least_squares
from joblib import Parallel, delayed
from tqdm import tqdm

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pyswarms as ps

from labothappy.toolbox.economics.cpi_data import actualize_price
from labothappy.correlations.turbomachinery.radial_compressor_losses import radial_compressor_rotor_losses, radial_compressor_stator_losses

import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Helper: create a blank state dict with keys H, S, P, D, A, V
# indexed by integer positions 1–5 (matching the old DataFrame indices).
# ---------------------------------------------------------------------------

class RadialCPMLDesign(object):

    # -----------------------------------------------------------------
    # Valeurs par défaut, appliquées à l'instanciation puis surchargeables
    # au cas par cas via set_parameters()/set_bounds(). Reprennent les
    # valeurs communes à tous les cas d'usage observés (CO2, R134a, Air).
    #
    # ⚠️ Volontairement ABSENTS des défauts (doivent être fournis explicitement,
    # sous peine d'erreur explicite dans sizing()) :
    #   - t_b, eps_imp, eps_bf_imp, k_imp : tolérances de fabrication,
    #     dépendantes de la machine réelle -- un défaut silencieux serait
    #     trompeur (l'écart CO2 0.762mm vs Air 2.11mm est significatif).
    #   - Omega (params) / Omega_bounds ou omega_bounds (bounds) : la vitesse
    #     de rotation (fixée ou optimisée) doit être un choix explicite.
    # -----------------------------------------------------------------
    DEFAULT_PARAMS = {
        # Géométrie / corrélations, constantes dans tous les cas observés
        'b3_b2_ratio': 1,
        'b5_b3': 1,
        'CP': 0.44,

        # Contraintes de conception, identiques dans tous les cas observés
        'M1s_rel_max': 1.4,
        'M1_rel_max': 0.9,
        'W2_W1s_min': 0.25,
        'alpha2_max': 85.0,
        'o1_min': 1e-3,
        'o1_max': 50e-3,
        'DR_min': -0.1,
        'DR_max': 0.9,
        'U2_max': 400.0,
        'r3_r2_min': 1.05,
        'r3_r2_max': 2.0,
    }

    # Clés de self.params qui n'ont pas de défaut et doivent être fournies
    # explicitement par set_parameters() avant sizing().
    REQUIRED_PARAMS_NO_DEFAULT = ('t_b', 'eps_imp', 'eps_bf_imp', 'k_imp')

    DEFAULT_BOUNDS = {
        'psi_is_bounds':  [0.3, 1.1],
        'r1s_r2_bounds':  [0.4, 0.7],
        'r1h_r1s_bounds': [0.25, 0.4],
        'b2_r2_bounds':   [0.02, 0.3],
        'r5_r3_bounds':   [1.01, 1.5],
        'r3_r2_bounds':   [1.05, 2.0],
        'xhi1_bounds':    [40, 70],
        'xhi2_bounds':    [20, 55],
        # Pas de défaut pour Omega_bounds/omega_bounds -- voir ci-dessus.
    }

    def __init__(self, fluid):
        # Inputs (point de fonctionnement thermo : mdot, p0_su, T0_su, p_ex)
        self.inputs = {}

        # Params (valeurs fixes/scalaires : t_b, eps_imp, CP, Omega si fixé, contraintes max/min, ...)
        # Initialisé avec les défauts de classe -- set_parameters() les surcharge ensuite au besoin.
        self.params = dict(self.DEFAULT_PARAMS)

        # Bounds (bornes d'optimisation *_bounds : psi_is_bounds, r1s_r2_bounds, ...)
        # Idem : initialisé avec les défauts de classe, surchargeable via set_bounds().
        self.bounds = dict(self.DEFAULT_BOUNDS)

        # RUN_KWARGS : kwargs passés à .sizing() par size_all_components(), comme pour les
        # autres sizing_models (REC.RUN_KWARGS, GH.RUN_KWARGS, etc.). Vide par défaut -> les
        # valeurs par défaut de sizing() (n_particles=100, max_iter=100, patience=15) s'appliquent.
        self.RUN_KWARGS = dict()

        # Sentinelle de pénalité, lue par size_all_components() via
        # `if hasattr(sizing_obj, "penalty"): if sizing_obj.penalty >= 1e6: raise ...`.
        # 1e6 = pas encore dimensionné / échec structurel. Mise à jour à chaque appel de
        # designSystem() (voir plus bas) et reflète, après sizing(), le meilleur design trouvé.
        self.penalty = 1e6

        # Flag résolu une seule fois dans sizing() à partir de self.bounds (et non plus
        # self.params comme dans l'original -- voir designSystem).
        self.optimize_omega = False

        # Abstract State 
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Blade Dictionary
        self.stages = []

        # Velocity Triangle Data
        self.Vel_Tri_R = {}
        self.Vel_Tri_S = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None
        
        self._STATE_KEYS = ('H', 'S', 'P', 'D', 'A', 'V')
        
        # State dicts – replaces pd.DataFrame(columns=[…], index=[1,2,3,4,5])
        self.total_states  = {k: {i: np.nan for i in range(1, 6)} for k in self._STATE_KEYS}
        self.static_states = {k: {i: np.nan for i in range(1, 6)} for k in self._STATE_KEYS}
        
        self.AS = CP.AbstractState('HEOS', fluid)
            
        # Nozzle and rotor losses initiated to 0
        self.losses = { 
            'DP0_S_volute' : 0,
            'Dh_S_nozzle' : 0,
        }
        
        self.CAPEX = {}
        
    def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.total_states['H'][position] = self.AS.hmass()            
        self.total_states['S'][position] = self.AS.smass()            
        self.total_states['P'][position] = self.AS.p()            
        self.total_states['D'][position] = self.AS.rhomass()            

        try:        
            self.total_states['A'][position] = self.AS.speed_sound()            
        except:
            self.total_states['A'][position] = -1  
            
        self.total_states['V'][position] = self.AS.viscosity()            
        
        return
    
    def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.static_states['H'][position] = self.AS.hmass()            
        self.static_states['S'][position] = self.AS.smass()            
        self.static_states['P'][position] = self.AS.p()            
        self.static_states['D'][position] = self.AS.rhomass()    
        
        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1            
            
        self.static_states['V'][position] = self.AS.viscosity()            

        return
    
    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_inputs(self, **parameters):
        """Point de fonctionnement thermo courant : mdot, p0_su, T0_su, p_ex."""
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
        """
        Paramètres fixes/scalaires : géométrie de référence (t_b, eps_imp, eps_bf_imp,
        k_imp, b3_b2_ratio, b5_b3, L_z, CP), Omega si imposé (sinon voir set_bounds avec
        Omega_bounds/omega_bounds), et les contraintes de conception
        (M1s_rel_max, M1_rel_max, W2_W1s_min, alpha2_max, o1_min/max, DR_min/max,
        U2_max, r3_r2_min/max).

        Un jeu de défauts raisonnable est déjà chargé à l'instanciation
        (voir RadialCPMLDesign.DEFAULT_PARAMS) pour les contraintes de
        conception et la géométrie constante entre cas -- set_parameters()
        n'a donc besoin de fournir QUE ce qui diffère du défaut, plus
        obligatoirement t_b/eps_imp/eps_bf_imp/k_imp (pas de défaut,
        dépendants de la machine réelle) et Omega si fixe.

        NE PREND PLUS les clés *_bounds -- utiliser set_bounds() pour celles-ci.
        """
        for key, value in parameters.items():
            if key.endswith('_bounds'):
                raise ValueError(
                    f"'{key}' est une borne d'optimisation : utiliser set_bounds({key}=...) "
                    f"plutôt que set_parameters({key}=...)."
                )
            self.params[key] = value

    def set_bounds(self, **bounds):
        """
        Bornes d'optimisation PSO : psi_is_bounds, r1s_r2_bounds, r1h_r1s_bounds,
        b2_r2_bounds, r5_r3_bounds, r3_r2_bounds, xhi1_bounds, xhi2_bounds, et
        optionnellement Omega_bounds (ou omega_bounds) si la vitesse de rotation
        doit être optimisée plutôt que fixée via set_parameters(Omega=...).

        Un jeu de bornes géométriques par défaut est déjà chargé à
        l'instanciation (voir RadialCPMLDesign.DEFAULT_BOUNDS) -- set_bounds()
        n'a besoin de fournir que ce qui diffère du défaut, plus
        Omega_bounds/omega_bounds si Omega n'est pas fixé via set_parameters.
        """
        for key, value in bounds.items():
            if not key.endswith('_bounds'):
                raise ValueError(
                    f"'{key}' n'est pas une clé de borne (doit finir par '_bounds')."
                )
            self.bounds[key] = value

    def _check_required_params(self):
        """
        Vérifie que les paramètres sans défaut (tolérances de fabrication,
        voir REQUIRED_PARAMS_NO_DEFAULT) ont bien été fournis via
        set_parameters() avant de lancer sizing(). Lève une erreur explicite
        sinon, plutôt que de planter plus loin avec un KeyError peu clair.
        """
        missing = [k for k in self.REQUIRED_PARAMS_NO_DEFAULT if k not in self.params]
        if missing:
            raise ValueError(
                f"Paramètre(s) requis manquant(s) (pas de valeur par défaut, "
                f"dépendants de la machine réelle) : {missing}. "
                f"À fournir via set_parameters(...)."
            )
    
    def _clone_light(self):
        """
        Clone léger pour évaluation parallèle d'une particule : recrée un objet
        RadialCPMLDesign avec les mêmes inputs/params/bounds/optimize_omega,
        mais un AbstractState CoolProp propre à lui. On évite deepcopy() car
        self.AS (CP.AbstractState) n'est en général pas deep-copiable proprement.
        """
        clone = RadialCPMLDesign(self.fluid)
        clone.inputs  = dict(self.inputs)
        clone.params  = dict(self.params)
        clone.bounds  = dict(self.bounds)
        clone.optimize_omega = self.optimize_omega
        return clone

    def _evaluate_particle_isolated(self, xi):
        """
        Évalue une particule (vecteur xi) sur un clone isolé de l'objet, pour
        permettre l'évaluation parallèle des particules d'une même itération
        PSO sans que les appels concurrents ne s'écrasent mutuellement l'état
        (self.total_states, self.omega, self.n_blade_R, ...).
        Retourne uniquement le coût scalaire -- l'état détaillé du clone est
        jeté après l'appel (seul le résultat sur self.designSystem(best_pos),
        exécuté en série à la fin de sizing(), fait foi pour l'état final).
        """
        clone = self._clone_light()
        try:
            return clone.designSystem(xi)
        except Exception:
            return 1e6
        
    # ---------------- Flow Computations ------------------------------------------------------------------

    def designRotor(self):
        
        "R0) Initiate Rotor Design"
        
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.total_states['S'][1])
        
        h0s = self.AS.hmass()
        self.Dh0s = h0s - self.total_states['H'][1]
        
        self.Vel_Tri_R['u2'] = u2 = np.sqrt(self.Dh0s / self.inputs['psi_is'])
        
        self.params['r2'] = r2 = u2 / self.omega
        
        self.Q   = self.inputs['mdot'] / self.total_states['D'][1]
        self.phi = 4 * self.Q / (np.pi * u2 * (2*r2)**2)
        
        "R1) Rotor Inlet"
        
        self.params['r1s'] = r1s = r2 * self.inputs['r1s_r2']
        self.params['r1h'] = r1h = r1s * self.inputs['r1h_r1s']
        self.params['b1']  = b1  = r1s - r1h
        
        self.params['r1']    = r1     = np.sqrt((r1s**2 + r1h**2) / 2)
        self.params['pitch1']= pitch1 = 2*np.pi*r1 / self.n_blade_R
        self.params['b2']    = b2     = r2 * self.inputs['b2_r2']

        self.params['r1']       = r1  = (r1s + r1h) / 2
        self.Vel_Tri_R['u1']    = u1  = self.omega * r1
        
        u1h = self.omega * r1h
        u1s = self.omega * r1s
        
        self.A1 = b1 * pitch1 * self.n_blade_R
        self.o1 = pitch1 * np.cos(np.pi/180 * self.inputs['xhi1']) - self.params['t_b']
        
        self.impeller_inlet_blockage = (
            1 - self.n_blade_R * self.params['t_b']
            / (np.pi * r1 * 2 * np.sin(np.pi * abs(self.inputs['xhi1']) / 180))
        )
        self.A1_th = np.pi * (r1s**2 - r1h**2) * self.impeller_inlet_blockage

        def compute_h1_new(h1): 
            self.update_static_AS(CP.HmassSmass_INPUTS, h1, self.total_states['S'][1], 1)
                        
            self.Vel_Tri_R['vm1'] = vm1 = self.inputs['mdot'] / (self.static_states['D'][1] * self.A1_th)
            
            self.Vel_Tri_R['beta1'] = beta1 = 0
            self.Vel_Tri_R['w1']    = w1    = vm1 / np.cos(beta1)
            self.Vel_Tri_R['wu1']   = w1 * np.sin(beta1)
            
            self.Vel_Tri_R['vu1']   = vu1 = self.Vel_Tri_R['wu1'] + self.Vel_Tri_R['u1']
            self.Vel_Tri_R['v1']    = v1  = np.sqrt(vm1**2 + vu1**2)
            self.Vel_Tri_R['alpha1']= np.arccos(vm1 / v1)
            
            self.beta1h = np.arctan(u1h / vm1)
            self.beta1s = np.arctan(u1s / vm1)
            
            self.w1s = vm1 / np.cos(self.beta1s)
            self.w1h = vm1 / np.cos(self.beta1h)
                        
            h1_new = self.total_states['H'][1] - v1**2 / 2
                        
            return h1_new

        def residual_h1(h1):
            return h1 - compute_h1_new(h1)
        
        h1_min = self.total_states['H'][1] * 0.95
        h1_max = self.total_states['H'][1]
        
        h1_solution = brentq(residual_h1, h1_min, h1_max, xtol=1e-6)
                
        "R2) Rotor Outlet"
        
        self.params['pitch2'] = pitch2 = 2*np.pi*r2 / self.n_blade_R
        self.A2_th = b2 * pitch2 * self.n_blade_R
        
        self.A2_th = (
            (2*np.pi*r2*b2 - self.n_blade_R*b2*self.params['t_b'])
            * np.cos(np.pi * abs(self.inputs['xhi2']) / 180)
        )
        
        # Slip Factor
        self.sigma = 1 - np.sqrt(np.cos(self.inputs['xhi2']*np.pi/180)) / self.n_blade_R**0.7
        angle_star_deg = 19 + 0.2*(90 - self.inputs['xhi2'])
        sigma_star = np.sin(np.pi/180 * angle_star_deg)
        
        r1_r2_lim = (self.sigma - sigma_star) / (1 - sigma_star)
        
        if r1/r2 > r1_r2_lim:
            num  = (r1/r2) - r1_r2_lim**(np.sqrt((90 - self.inputs['xhi2'])/10))
            fact = 1 - (num / (1 - sigma_star))
            self.sigma = self.sigma * fact
        
        def system_rotor(vm2):

            self.Vel_Tri_R['vm2'] = vm2
        
            rho2 = self.inputs['mdot'] / (vm2 * self.A2_th)
            
            self.Vel_Tri_R['vu2'] = vu2 = (
                self.sigma * self.Vel_Tri_R['u2']
                - vm2 * np.tan(self.inputs['xhi2'] * np.pi/180)
            )
            
            self.dh0 = vu2*self.Vel_Tri_R['u2'] - self.Vel_Tri_R['vu1']*self.Vel_Tri_R['u1']
            
            h02 = self.total_states['H'][1] + self.dh0
            
            self.Vel_Tri_R['v2'] = v2 = np.sqrt(vu2**2 + vm2**2)
            h2 = h02 - v2**2 / 2
            
            try:
                self.AS.update(CP.DmassHmass_INPUTS, h2, rho2)
                s2 = self.AS.smass()
            except:   
                try:
                    s2 = PropsSI('S', 'H', h2, 'D', rho2, self.AS.fluid_names()[0]) 
                except:
                    return 1e6

            self.update_static_AS(CP.HmassSmass_INPUTS, h2, s2, 2)
            self.update_total_AS(CP.HmassSmass_INPUTS, h02, self.static_states['S'][2], 2)

            self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][2], self.total_states['S'][1])
            h_is = self.AS.hmass()
            
            self.Vel_Tri_R['alpha2'] = alpha2 = np.arccos(vm2 / v2)
            
            self.params['L_z'] = L_z = self.params['r2'] * (0.1 + 2*self.phi)

            self.Vel_Tri_R['wu2'] = wu2 = self.Vel_Tri_R['vu2'] - self.Vel_Tri_R['u2']
            self.Vel_Tri_R['w2']  = w2  = np.sqrt(wu2**2 + vm2**2)
            self.Vel_Tri_R['beta2'] = beta2 = np.arccos(vm2 / w2)
            
            self.rotor_losses = radial_compressor_rotor_losses(
                A1=self.A1, A1_th=self.A1_th, alpha2=alpha2,
                beta1=self.Vel_Tri_R['beta1'], beta1h=self.beta1h, beta1s=self.beta1s,
                beta2=self.Vel_Tri_R['beta2'], b2=self.params['b2'], C_df=0.004, C_fi=0.004,
                Dh0=self.dh0, eps_a=self.params['eps_imp'], eps_b=self.params['eps_bf_imp'],
                eps_r=self.params['eps_imp'], k_roughness=self.params['k_imp'], L_z=L_z,
                mdot=self.inputs['mdot'], mu1=self.static_states['V'][1],
                mu2=self.static_states['V'][2], n_bl_r=self.n_blade_R,
                rho1=self.static_states['D'][1], rho2=self.static_states['D'][2],
                r1h=self.params['r1h'], r1s=self.params['r1s'], r2=self.params['r2'],
                u2=self.Vel_Tri_R['u2'], vu2=vu2, v1m=self.Vel_Tri_R['vm1'], v2=v2,
                w1=self.Vel_Tri_R['w1'], w1_th=self.Vel_Tri_R['w1'], w1h=self.w1h,
                w1s=self.w1s, w2=w2, xhi1=self.inputs['xhi1']*np.pi/180,
                xhi2=self.inputs['xhi2']*np.pi/180,
            )
            
            h02_new = self.rotor_losses['tot'] + h_is
            h2_new  = h02_new - v2**2 / 2
            
            self.update_static_AS(CP.HmassSmass_INPUTS, h2_new, self.static_states['S'][2], 2)
            
            try:
                s2 = self.static_states['S'][2]
                self.update_total_AS(CP.HmassSmass_INPUTS, h02, s2, 2)
            except Exception:
                return 1e6
            
            res = (h2 - h2_new) / h2_new
            self.res_rotor_ex = res
            
            return res
        
        rho2_max_phys = 3.0 * self.static_states['D'][1]
        vm2_min = self.inputs['mdot'] / (rho2_max_phys * self.A2_th)
        vm2_max = self.inputs['mdot'] / (self.static_states['D'][1] * self.A2_th)
        
        res_min = system_rotor(vm2_min)
        res_max = system_rotor(vm2_max)
        
        if res_min * res_max > 0:
            raise ValueError()
        
        self.sol_rotor_ex = brentq(system_rotor, vm2_min, vm2_max, xtol=1e-6)
        
        if self.rotor_losses['tot'] > self.dh0:
            raise ValueError()
        
        "Compute constraint terms"
        
        self.M1s_rel = self.w1s / self.static_states['A'][1]
        self.M1_rel  = self.Vel_Tri_R['w1'] / self.static_states['A'][1]
        self.W2_W1s  = self.Vel_Tri_R['w2'] / self.w1s
        
        vu1 = self.Vel_Tri_R['vu1']
        vu2 = self.Vel_Tri_R['vu2']
        u2  = self.Vel_Tri_R['u2']
        self.DR = 1.0 - (vu2 + vu1) / (2.0 * u2)
        
        return

#%%

    def designStator(self):

        "S3) Vaneless Space Exhaust"
        self.static_states['P'][3] = p3 = (
            self.static_states['P'][2]
            + self.params["CP"] * (self.total_states['P'][2] - self.static_states['P'][2])
        )
        self.total_states['H'][3] = h03 = self.total_states['H'][2]
        
        def vaneless_system(x):
        
            s3 = x
        
            self.update_total_AS(CP.HmassSmass_INPUTS, self.total_states['H'][3], s3, 3)
            self.update_static_AS(CP.PSmass_INPUTS, self.static_states['P'][3], s3, 3)
        
            p03 = self.total_states['P'][3]
            K   = (self.total_states['P'][2] - p03) / (self.total_states['P'][2] - self.static_states['P'][2])
            CP_id = self.params['CP'] + K 
        
            self.Vel_Tri_S['v3']    = v3  = np.sqrt(2*(self.total_states['H'][3] - self.static_states['H'][3]))
            self.Vel_Tri_S['alpha3'] = self.Vel_Tri_R['alpha2']
            self.Vel_Tri_S['vm3']   = vm3 = v3 * np.cos(self.Vel_Tri_S['alpha3'])
            self.Vel_Tri_S['vu3']   = vu3 = v3 * np.sin(self.Vel_Tri_S['alpha3'])
            self.Vel_Tri_S['u3']    = u3  = self.Vel_Tri_R['u2'] * self.inputs['r3_r2']
            self.Vel_Tri_S['wu3']   = wu3 = vu3 - u3
            self.Vel_Tri_S['wm3']   = wm3 = vm3
            self.Vel_Tri_S['w3']    = np.sqrt(wm3**2 + wu3**2)
            self.Vel_Tri_S['beta3'] = beta3 = np.arctan(wu3 / wm3)
            
            res = (self.inputs['mdot'] - (vm3 * self.static_states['D'][3] * self.A2_th) * (1-CP_id)**(-0.5))**2
            
            return res
            
        self.AS.update(CP.HmassP_INPUTS, h03, p3)
        s3_max = self.AS.smass()
        s3_min = self.static_states['S'][2]
        
        sol = minimize(vaneless_system, s3_min, method='L-BFGS-B', bounds=[(s3_min, s3_max)],
                       options={'ftol': 1e-8, 'gtol': 1e-8})
    
        "S4) Vaned Diffuser Inlet"
        s4  = self.static_states['S'][3]
        h04 = self.total_states['H'][3]
        
        self.update_total_AS(CP.HmassSmass_INPUTS, h04, s4, 4)
        self.params['r3']  = self.params['r2'] * self.inputs['r3_r2']
        self.Vel_Tri_R['u3'] = self.Vel_Tri_R['u2'] * self.inputs['r3_r2']
        self.params['b3']  = self.params['b2'] * self.params['b3_b2_ratio']

        self.params['pitch3'] = pitch3 = 2*np.pi*self.params['r3'] / self.n_blade_R
        self.A3 = pitch3 * self.params['b2'] * self.n_blade_R
        
        def vaned_diffuser_inlet_system(x):
            h4 = x
            self.update_static_AS(CP.HmassSmass_INPUTS, h4, s4, 4)

            self.Vel_Tri_S['v4']  = v4  = np.sqrt(2*(self.total_states['H'][4] - h4))
            self.Vel_Tri_S['vu4'] = vu4 = self.Vel_Tri_S['vu3']
            self.Vel_Tri_S['vm4'] = vm4 = np.sqrt(max(v4**2 - vu4**2, 0))
            self.Vel_Tri_S['alpha4'] = alpha4 = np.arctan(vu4 / vm4)
            
            self.Vel_Tri_S['u4']   = u4  = self.Vel_Tri_R['u3']
            self.Vel_Tri_S['wu4']  = wu4 = vu4 - u4
            self.Vel_Tri_S['wm4']  = wm4 = vm4
            self.Vel_Tri_S['w4']   = np.sqrt(wm4**2 + wu4**2)
            self.Vel_Tri_S['beta4']= beta4 = np.arctan(wu4 / wm4)
            
            res = vm4 - self.inputs['mdot'] / (self.static_states['D'][4] * self.A3)
            self.res_diff_in = res
            return res
        
        h4_min = self.static_states['H'][3] * 0.9
        h4_max = self.total_states['H'][4]
        
        self.sol_vaned_diff = brentq(vaned_diffuser_inlet_system, h4_min, h4_max, xtol=1e-6)
        
        self.o4    = pitch3 * np.cos(self.Vel_Tri_S['alpha4']) - self.params['t_b']
        self.A3_th = self.o4 * self.params['b2'] * self.n_blade_R
        
        self.params['xhi4'] = self.Vel_Tri_S['alpha4']

        "S5) Vaned Diffuser Outlet"
        
        self.params['r5'] = self.params['r3'] * self.inputs['r5_r3']
        self.A5 = self.A3_th * self.params['b5_b3'] * self.inputs['r5_r3']
        self.params['b5'] = self.params['b3'] * self.params['b5_b3']
        
        def system_stator(x):
            alpha5 = x[0]
            v5     = x[1]
        
            vu5 = v5 * np.sin(alpha5)
            vm5 = v5 * np.cos(alpha5)
        
            self.Vel_Tri_S['alpha5'] = alpha5
            self.Vel_Tri_S['vm5']    = vm5
            self.Vel_Tri_S['v5']     = v5
            self.Vel_Tri_S['vu5']    = vu5
        
            h05 = self.total_states['H'][4]
            h5  = h05 - v5**2 / 2
        
            self.stator_losses = radial_compressor_stator_losses(
                A4=self.A3, A4_th=self.A3_th, beta4=self.Vel_Tri_S['beta4'], C_f=0.004,
                r3=self.params['r3'], r4=self.params['r3'], r5=self.params['r5'],
                vm=vm5, w4=self.Vel_Tri_S['w4'], xhi3=self.Vel_Tri_S['alpha3'],
                xhi4=self.Vel_Tri_S['alpha4'], xhi5=alpha5,
            )
        
            self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.static_states['S'][4])
            h5is = self.AS.hmass()
        
            h5_new = h5is + self.stator_losses['tot']
            self.AS.update(CP.HmassSmass_INPUTS, h5_new, self.static_states['S'][4])
            p5 = self.AS.p()
            self.update_static_AS(CP.HmassP_INPUTS, h5_new, p5, 5)
            self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)
        
            res1 = (h5 - h5_new) / h05
            res2 = (self.A5*self.Vel_Tri_S['vm5']*self.static_states['D'][5] - self.inputs['mdot'])/self.inputs['mdot']
        
            self.params['xhi5'] = alpha5
        
            return [res1, res2]
        
        alpha5_guess = self.Vel_Tri_S['alpha4'] * 0.5
        v5_guess     = self.inputs['mdot'] / (self.static_states['D'][4] * self.A5)
        
        v5_upper    = self.inputs['mdot'] / (self.static_states['D'][4] * self.A5)
        alpha5_lb   = 0 * np.pi/180
        alpha5_ub   = 85 * np.pi/180
        
        self.sol_sys_stator = least_squares(
            system_stator,
            x0     = [alpha5_guess, v5_guess],
            bounds = ([alpha5_lb, 1.0], [alpha5_ub, v5_upper]),
            method = 'trf',
        )
        
        alpha5_sol, v5_sol = self.sol_sys_stator.x
        
        return

#%%

    def cost_estimation(self):
        
        if self.fluid == 'CO2' or self.fluid == 'CarbonDioxide' or self.fluid == 'R744':
            """
            SCO2 POWER CYCLE COMPONENT COST CORRELATIONS FROM DOE DATA
            SPANNING MULTIPLE SCALES AND APPLICATIONS (2019)
            
            Nathan T. Weiland,  Blake W. Lance, Sandeep R. Pidaparti
            
            Especially good for 1.5 - 200 MW 
            Based on 2017 CEPCI (chemical plant cost index) for dollars
            """
        
            W_dot_MW = self.W_dot/1e6
            CAPEX_compressor = actualize_price(1230000 * W_dot_MW**0.3392, 2017, "USD")

            self.CAPEX['Compressor'] = CAPEX_compressor
        
        else:
            """
            Key components for Carnot Battery: Technology review, technical barriers and selection criteria
            
            Ting Liang, Andrea Vecchi, Kai Knobloch, Adriano Sciacovelli, Kurt Engelbrecht, Yongliang Li, Yulong Ding
            
            Based on 2017 CEPCI (chemical plant cost index) for dollars
            """
            
            CAPEX_fun = 39.5 * self.inputs["mdot"] * (self.PR * np.log10(self.PR))/(0.9-self.eta_is)
            CAPEX_compressor = actualize_price(CAPEX_fun, 2017, "USD")
            self.CAPEX['Compressor'] = CAPEX_compressor

        
        # Generator Costs
        
        """
        SCO2 POWER CYCLE COMPONENT COST CORRELATIONS FROM DOE DATA
        SPANNING MULTIPLE SCALES AND APPLICATIONS (2019)
        
        Nathan T. Weiland,  Blake W. Lance, Sandeep R. Pidaparti
        
        Especially good for 1.5 - 200 MW 
        Based on 2017 CEPCI (chemical plant cost index) for dollars
        """
        
        self.eta_alt = 0.97
        self.W_dot_el = self.W_dot*self.eta_alt
        
        W_dot_el_MW = self.W_dot_el/1e6
        
        CAPEX_alternator = actualize_price(108900 * W_dot_el_MW**0.5463, 2017, "EUR")
        self.CAPEX['Alternator'] = CAPEX_alternator
        
        self.f_install = 0.35 
        self.CAPEX['Installation'] = self.f_install*(self.CAPEX['Alternator'] + self.CAPEX['Compressor'])
        
        self.CAPEX['Total'] = self.CAPEX['Compressor'] + self.CAPEX['Alternator'] + self.CAPEX['Installation']
            
        return
    
#%%

    def designSystem(self, x):
        
        self.penalty_factor = 100

        # Sentinelle : tant qu'on n'a pas atteint un design complet, self.penalty
        # reste à 1e6 -- lu par size_all_components() via hasattr(sizing_obj, "penalty").
        self.penalty = 1e6

        if self.optimize_omega:
            self.inputs['psi_is']   = x[0]
            self.inputs['r1s_r2']   = x[1]
            self.inputs['r1h_r1s']  = x[2]
            self.inputs['b2_r2']    = x[3]
            self.inputs['r5_r3']    = x[4]
            self.inputs['r3_r2']    = x[5]
            self.inputs['xhi1']     = x[6]
            self.inputs['xhi2']     = x[7]
            self.inputs['Omega']    = self.params['Omega'] = x[8]
        else:
            self.inputs['psi_is']   = x[0]
            self.inputs['r1s_r2']   = x[1]
            self.inputs['r1h_r1s']  = x[2]
            self.inputs['b2_r2']    = x[3]
            self.inputs['r5_r3']    = x[4]
            self.inputs['r3_r2']    = x[5]
            self.inputs['xhi1']     = x[6]
            self.inputs['xhi2']     = x[7]
        
        self.update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        self.update_static_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
    
        self.PR       = self.inputs['p_ex'] / self.inputs['p0_su']
        self.n_blade_R = np.floor(12.03 + 2.544*self.PR)
        self.omega    = self.params['Omega'] * (2*np.pi) / 60
        
        try:            
            self.designRotor()
            
            eta_diff_min = 0.95
            p_req = self.inputs['p_ex'] / eta_diff_min
            
            if self.total_states['P'][2] < p_req:
                res = abs((p_req - self.total_states['P'][2]) / p_req)
                # échec structurel : self.penalty reste à 1e6 (déjà posé plus haut)
                return res * 500

        except:
            return 2222
        
        try:            
            self.designStator()
        except:
            return 1100
        
        p = self.params
    
        self.error_log = []
        penalty = 0.0
        
        if self.M1s_rel > p['M1s_rel_max']:
            delta = (self.M1s_rel - p['M1s_rel_max']) / p['M1s_rel_max']
            penalty += delta
            self.error_log.append(('M1s_rel', delta))
        
        if self.M1_rel > p['M1_rel_max']:
            delta = (self.M1_rel - p['M1_rel_max']) / p['M1_rel_max']
            penalty += delta
            self.error_log.append(('M1_rel', delta))
        
        if self.W2_W1s < p['W2_W1s_min']:
            delta = (p['W2_W1s_min'] - self.W2_W1s) / p['W2_W1s_min']
            penalty += delta
            self.error_log.append(('W2_W1s', delta))
        
        alpha2_deg = abs(self.Vel_Tri_R['alpha2']) * 180.0 / np.pi
        if alpha2_deg > p['alpha2_max']:
            delta = (alpha2_deg - p['alpha2_max']) / p['alpha2_max']
            penalty += delta
            self.error_log.append(('alpha2_deg', delta))
        
        if self.o1 < p['o1_min']:
            delta = (p['o1_min'] - self.o1) / p['o1_min']
            penalty += delta
            self.error_log.append(('o1_low', delta))
        elif self.o1 > p['o1_max']:
            delta = (self.o1 - p['o1_max']) / p['o1_max']
            penalty += delta
            self.error_log.append(('o1_high', delta))
        
        if self.DR < p['DR_min']:
            delta = (p['DR_min'] - self.DR) / (p['DR_max'] - p['DR_min'])
            penalty += delta
            self.error_log.append(('DR_low', delta))
        elif self.DR > p['DR_max']:
            delta = (self.DR - p['DR_max']) / (p['DR_max'] - p['DR_min'])
            penalty += delta
            self.error_log.append(('DR_high', delta))
        
        if self.Vel_Tri_R['u2'] > p['U2_max']:
            delta = (self.Vel_Tri_R['u2'] - p['U2_max']) / p['U2_max']
            penalty += delta
            self.error_log.append(('u2', delta))
        
        if self.inputs['r3_r2'] < p['r3_r2_min']:
            delta = (p['r3_r2_min'] - self.inputs['r3_r2']) / p['r3_r2_min']
            penalty += delta
            self.error_log.append(('r3_r2_low', delta))
        elif self.inputs['r3_r2'] > p['r3_r2_max']:
            delta = (self.inputs['r3_r2'] - p['r3_r2_max']) / p['r3_r2_max']
            penalty += delta
            self.error_log.append(('r3_r2_high', delta))
        
        # Total-static efficiency
        self.AS.update(CP.PSmass_INPUTS, self.static_states['P'][5], self.total_states['S'][1])
        hout_is = self.AS.hmass()
        self.eta_is = (hout_is - self.total_states['H'][1]) / \
                      (self.static_states['H'][5] - self.total_states['H'][1])
        
        # Total-total efficiency
        self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][5], self.total_states['S'][1])
        hout_is = self.AS.hmass()
        self.eta_is_tt = (hout_is - self.total_states['H'][1]) / \
                         (self.total_states['H'][5] - self.total_states['H'][1])
        
        if abs(self.res_rotor_ex) > 1e-3:
            penalty += self.res_rotor_ex
            
        if penalty > 0:
            # Design faisable mais contraintes violées : self.penalty reflète
            # l'ampleur de la violation (généralement << 1e6, donc pas rejeté
            # par size_all_components, mais dégradé dans l'objectif).
            self.penalty = penalty
            return -self.eta_is + penalty * self.penalty_factor
        
        self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][5], self.total_states['S'][1])
        h05_is = self.AS.hmass()
        self.Dh0s_2 = h05_is - self.total_states['H'][1]

        # Design pleinement faisable : pénalité nulle.
        self.penalty = 0.0
        
        self.W_dot = (self.total_states['H'][5] - self.total_states['H'][1])*self.inputs["mdot"]
        
        return -self.eta_is
    
    def sizing(self, n_particles=100, max_iter=100, patience=15, n_jobs=1):
        """
        n_jobs : nombre de workers joblib pour l'évaluation parallèle des
        particules à chaque itération PSO (1 = série, -1 = tous les coeurs).
        Backend "threading" utilisé (et non "loky"/process) car
        CP.AbstractState n'est en général pas picklable pour le multiprocessing ;
        le gain vient du fait que les appels CoolProp/scipy relâchent le GIL
        pendant une bonne partie du calcul.
        """
        # Vérifie que t_b/eps_imp/eps_bf_imp/k_imp ont bien été fournis
        # (pas de défaut -- voir REQUIRED_PARAMS_NO_DEFAULT).
        self._check_required_params()

        self.optimize_omega = ("Omega_bounds" in self.bounds) or ("omega_bounds" in self.bounds)
    
        if self.optimize_omega and "Omega" not in self.params:
            omega_b = self.bounds.get("Omega_bounds", self.bounds.get("omega_bounds"))
            self.params['Omega'] = 0.5 * (omega_b[0] + omega_b[1])
        elif not self.optimize_omega and "Omega" not in self.params:
            raise ValueError(
                "Omega n'est ni fixé (set_parameters(Omega=...)) ni optimisé "
                "(set_bounds(Omega_bounds=...) ou set_bounds(omega_bounds=...))."
            )
    
        if self.optimize_omega:
            omega_b = self.bounds.get("Omega_bounds", self.bounds.get("omega_bounds"))
            bounds = (np.array([
                self.bounds['psi_is_bounds'][0], self.bounds['r1s_r2_bounds'][0],
                self.bounds['r1h_r1s_bounds'][0], self.bounds['b2_r2_bounds'][0],
                self.bounds['r5_r3_bounds'][0], self.bounds['r3_r2_bounds'][0],
                self.bounds['xhi1_bounds'][0], self.bounds['xhi2_bounds'][0],
                omega_b[0],
            ]),
            np.array([
                self.bounds['psi_is_bounds'][1], self.bounds['r1s_r2_bounds'][1],
                self.bounds['r1h_r1s_bounds'][1], self.bounds['b2_r2_bounds'][1],
                self.bounds['r5_r3_bounds'][1], self.bounds['r3_r2_bounds'][1],
                self.bounds['xhi1_bounds'][1], self.bounds['xhi2_bounds'][1],
                omega_b[1],
            ]))
        else:
            bounds = (np.array([
                self.bounds['psi_is_bounds'][0], self.bounds['r1s_r2_bounds'][0],
                self.bounds['r1h_r1s_bounds'][0], self.bounds['b2_r2_bounds'][0],
                self.bounds['r5_r3_bounds'][0], self.bounds['r3_r2_bounds'][0],
                self.bounds['xhi1_bounds'][0], self.bounds['xhi2_bounds'][0],
            ]),
            np.array([
                self.bounds['psi_is_bounds'][1], self.bounds['r1s_r2_bounds'][1],
                self.bounds['r1h_r1s_bounds'][1], self.bounds['b2_r2_bounds'][1],
                self.bounds['r5_r3_bounds'][1], self.bounds['r3_r2_bounds'][1],
                self.bounds['xhi1_bounds'][1], self.bounds['xhi2_bounds'][1],
            ]))
    
        def objective_wrapper(x):
            if n_jobs == 1:
                costs = [self.designSystem(xi) for xi in x]
            else:
                costs = Parallel(n_jobs=n_jobs, prefer="threads")(
                    delayed(self._evaluate_particle_isolated)(xi) for xi in x
                )
            return np.asarray(costs, dtype=float)
    
        optimizer = ps.single.GlobalBestPSO(
            n_particles=n_particles,
            dimensions=len(bounds[0]),
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds,
        )
    
        tol                 = 1e-3
        no_improve_counter  = 0
        best_cost           = np.inf
    
        pbar = tqdm(range(max_iter), desc="Radial compressor sizing", unit="it")
        for i in pbar:
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = optimizer.swarm.best_cost
    
            batch_best = getattr(self, "_last_batch_max_wdot", self.inputs.get("W_dot", 0.0))
            if batch_best > self.inputs.get("W_dot", 0.0):
                self.inputs["W_dot"] = batch_best
    
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
    
            pbar.set_postfix({
                "best_cost": f"{best_cost:.4g}",
                "stagnation": f"{no_improve_counter}/{patience}",
            })
    
            if no_improve_counter >= patience:
                pbar.set_postfix({
                    "best_cost": f"{best_cost:.4g}",
                    "status": "stopped (stagnation)",
                })
                pbar.close()
                break
        else:
            pbar.close()
    
        best_pos = optimizer.swarm.best_pos
        # Recalcule en série sur self (pas un clone) pour que l'état final
        # (self.eta_is, self.penalty, self.total_states, ...) reflète bien
        # le meilleur design trouvé.
        self.designSystem(best_pos)
        
        self.cost_estimation()
        
        return

    # Alias de rétrocompatibilité si du code appelle encore .design()
    def design(self, n_particles=100, max_iter=100, patience=15, n_jobs=1):
        return self.sizing(n_particles=n_particles, max_iter=max_iter, patience=patience, n_jobs=n_jobs)

    # ---------------- Compatibilité pipeline (export JSON dans cycle_design) --------------------

    def export_params_dict(self):
        """
        ⚠️ TODO CAPEX : contrairement aux autres sizing_models (REC, GH, CD, PP, TA, TR),
        cet objet n'a pas encore de corrélation de coût -> pas de self.CAPEX['Total'].
        Il ne peut donc PAS être inséré tel quel dans sizing_models du pipeline principal
        (size_all_components/size_components planteront sur obj.CAPEX['Total']) tant
        qu'une corrélation de coût n'a pas été ajoutée et que self.CAPEX = {'Total': ...}
        n'est pas peuplé par sizing().
        Stub minimal en attendant, pour permettre un export JSON individuel du compresseur.
        """
        return {
            'inputs': dict(self.inputs),
            'params': {k: v for k, v in self.params.items() if not k.endswith('_bounds')},
            'bounds': dict(self.bounds),
            'eta_is': getattr(self, 'eta_is', None),
            'eta_is_tt': getattr(self, 'eta_is_tt', None),
            'penalty': getattr(self, 'penalty', None),
            'n_blade_R': getattr(self, 'n_blade_R', None),
            'Omega': self.params.get('Omega'),
            'r2': self.params.get('r2'),
        }


if __name__ == "__main__":

    fluid = "CO2_MW" # CO2 / CO2_MW / R134a / Air_1 / Air_2 / Air_3    
    
    eta_is_vec = []
    
    for i in range(1):
        
        if fluid == "CO2":
            Comp_des = RadialCPMLDesign('CO2')
            
            Comp_des.set_inputs(
                mdot  = 2.15,
                p0_su = 76.9*1e5,
                T0_su = 305.97,
                p_ex  = 96.894*1e5,
            )
            
            # Seuls t_b/eps_imp/eps_bf_imp/k_imp (pas de défaut) et Omega/L_z
            # (cas particulier) sont fournis ici -- tout le reste (CP,
            # b3_b2_ratio, b5_b3, contraintes M1s_rel_max, etc.) vient des
            # défauts de classe (DEFAULT_PARAMS).
            Comp_des.set_parameters(
                t_b          = 0.762*1e-3,
                eps_imp      = 0.254*1e-3,
                eps_bf_imp   = 0.254*1e-3,
                k_imp        = 0.01*1e-3,
                L_z          = 0.1137,
                Omega        = 50000,
            )

            # Toutes les bornes géométriques ici correspondent au défaut
            # (DEFAULT_BOUNDS) -- cet appel est donc facultatif dans ce cas
            # précis, laissé pour l'exemple/la clarté.
            Comp_des.set_bounds(
                psi_is_bounds  = [0.3, 1.1],
                r1s_r2_bounds  = [0.4, 0.7],
                r1h_r1s_bounds = [0.25, 0.4],
                b2_r2_bounds   = [0.02, 0.1],
                r5_r3_bounds   = [1.01, 1.5],
                r3_r2_bounds   = [1.05, 2],
                xhi1_bounds    = [40, 70],
                xhi2_bounds    = [20, 55],
            )
            
            Comp_des.sizing()
        
        elif fluid == "CO2_MW":
            Comp_des = RadialCPMLDesign('CO2')
            
            Comp_des.set_inputs(
                mdot  = 5*10.84,
                p0_su = 4069717,
                T0_su = 299.13,
                p_ex  = 17364328,
            )
            
            Comp_des.set_parameters(
                t_b          = 0.762*1e-3,
                eps_imp      = 0.254*1e-3,
                eps_bf_imp   = 0.254*1e-3,
                k_imp        = 0.01*1e-3,
            )

            Comp_des.set_bounds(
                Omega_bounds   = [1000, 200000],   # Omega non fixé ici -> optimisé par le PSO
                b2_r2_bounds   = [0.02, 0.1],       # diffère du défaut [0.02, 0.3]
            )
                        
            Comp_des.sizing()
        
        elif fluid == "R134a":
            Comp_des = RadialCPMLDesign('R134a')
            
            Comp_des.set_inputs(
                mdot  = 0.039,
                p0_su = 1.65*1e5,
                T0_su = 265,
                p_ex  = 3.8775*1e5,
            )
            
            Comp_des.set_parameters(
                t_b          = 0.1*1e-3,
                eps_imp      = 0.15*1e-3,
                eps_bf_imp   = 1*1e-3,
                k_imp        = 0.01*1e-3,
                L_z          = 0.007693,
                Omega        = 180000,
            )

            Comp_des.set_bounds(
                psi_is_bounds  = [0.3, 0.5],
                r1s_r2_bounds  = [0.6, 0.7],
                r1h_r1s_bounds = [0.3, 0.4],
                b2_r2_bounds   = [0.02, 0.3],
                r5_r3_bounds   = [1.01, 1.3],
                r3_r2_bounds   = [1.05, 1.2],
                xhi1_bounds    = [40, 60],
                xhi2_bounds    = [40, 55],
            )
            
            Comp_des.sizing()
        
        elif fluid == "Air_1":
            Comp_des = RadialCPMLDesign('Air')
            
            Comp_des.set_inputs(
                mdot  = 5.32,
                p0_su = 1.01*1e5,
                T0_su = 288.15,
                p_ex  = 2.02*1e5,
            )
            
            Comp_des.set_parameters(
                t_b          = 2.11*1e-3,
                eps_imp      = 0.372*1e-3,
                eps_bf_imp   = 0.372*1e-3,
                k_imp        = 0.002*1e-3,
                L_z          = 0.13,
                Omega        = 14000,
            )

            Comp_des.set_bounds(
                psi_is_bounds  = [0.3, 0.7],
                b2_r2_bounds   = [0.02, 0.8],
            )
            
            Comp_des.sizing()
        
        elif fluid == "Air_2":
            Comp_des = RadialCPMLDesign('Air')
            
            Comp_des.set_inputs(
                mdot  = 4.54,
                p0_su = 1.01*1e5,
                T0_su = 288.15,
                p_ex  = 1.8281*1e5,
            )
            
            Comp_des.set_parameters(
                t_b          = 2.11*1e-3,
                eps_imp      = 0.235*1e-3,
                eps_bf_imp   = 0.235*1e-3,
                k_imp        = 0.002*1e-3,
                L_z          = 0.13,
                Omega        = 14000,
            )

            Comp_des.set_bounds(
                b2_r2_bounds   = [0.02, 0.8],
            )
            
            Comp_des.sizing()
            
        elif fluid == "Air_3":
            Comp_des = RadialCPMLDesign('Air')
            
            Comp_des.set_inputs(
                mdot  = 4.54,
                p0_su = 1.01*1e5,
                T0_su = 288.15,
                p_ex  = 1.6867*1e5,
            )
            
            Comp_des.set_parameters(
                t_b          = 2.11*1e-3,
                eps_imp      = 0.372*1e-3,
                eps_bf_imp   = 0.372*1e-3,
                k_imp        = 0.002*1e-3,
                L_z          = 0.13,
                Omega        = 14000,
            )

            Comp_des.set_bounds(
                b2_r2_bounds   = [0.02, 0.8],
            )
            
            Comp_des.sizing()
        
        try:
            eta_is_vec.append(Comp_des.eta_is)    
        except:
            eta_is_vec.append(-1)