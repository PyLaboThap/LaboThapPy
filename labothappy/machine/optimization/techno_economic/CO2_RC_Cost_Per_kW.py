#%%

# -*- coding: utf-8 -*-
"""
co2_rc_full_design_optimizer.py
Étend CO2RCOptimizer (importé de co2_rc_pso_optimizer.py) avec le
dimensionnement des composants + CAPEX + boucle cycle_design.
"""

#%% Imports

import numpy as np
from CoolProp.CoolProp import PropsSI

from labothappy.sizing.turbomachinery.turbine.axial.sizing_1D.mean_line_axial_turbine_loss_model_sizing import AxialTurbineMeanLineSizing
from labothappy.sizing.turbomachinery.turbine.radial.mean_line_radial_turbine_loss_model_sizing import RadialTurbineMeanLineSizing
from labothappy.sizing.heat_exchanger.shell_and_tube.shell_and_tube_sizing import ShellAndTubeSizingOpt
from labothappy.sizing.heat_exchanger.PCHE.PCHE_sizing import PCHESizingOpt
from labothappy.sizing.turbomachinery.pump.radial.radial_pump_0D_design import RadialPumpODDesign

# --- Import de la brique d'optimisation (fichier 1) ---
from labothappy.machine.optimization.thermodynamic.CO2_RC_HX_presize_Optimization import CO2RC_HX_optimizer

import warnings
warnings.filterwarnings('ignore')

#%% Dimensionnement des composants (inchangé par rapport au fichier 2 d'origine)

def TCO2_rec_comp_sizing(RC, turb_choice):

    # --- Recuperator ---
    REC_model = RC.components['Recuperator'].model
    REC_sizing = RC.components['Recuperator'].sizing = PCHESizingOpt()

    REC_sizing.set_inputs(
        fluid_H=REC_model.su_H.fluid, T_su_H=REC_model.su_H.T,
        P_su_H=REC_model.su_H.p, m_dot_H=REC_model.su_H.m_dot,
        fluid_C=REC_model.su_C.fluid, T_su_C=REC_model.su_C.T,
        P_su_C=REC_model.su_C.p, m_dot_C=REC_model.su_C.m_dot,
    )
    
    REC_sizing.set_parameters(
        k_cond=20, R_p=1, n_disc=50,
        Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True,

        H_Corr={"1P": "Gnielinski", "SC": "Gnielinski"},
        C_Corr={"1P": "Gnielinski", "SC": "Gnielinski"},
        H_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP"},
        C_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP"},

        Q_dot=REC_model.Q.Q_dot,
        DP_h=REC_model.DP_h,
        DP_c=REC_model.DP_c,
    )
    
    REC_sizing.set_bounds(
        alpha=[10, 40], D_c=[1e-3, 3e-3],
        L_x=[0.2, 1.5], L_y=[0.2, 2.3], L_z=[0.2, 0.6],
        n_parallel=[1, 20], n_series=[1, 20],
    )
        
    try:
        REC_sizing.sizing_parallel(n_jobs=-1, n_particles=100, max_iter=50, patience=20)
    
        if REC_sizing.score >= 1000000:
            raise ValueError("Recuperator Sizing did not Converge")
    
    except Exception as e:
        print(f"⚠️ Failed to design Recuperator: {e}")
        # REC_model.su_H.print_resume()
        # REC_model.su_C.print_resume()
        # print(f"Q_dot_constr : {REC_sizing.Q_dot_constr}")
        # print(f"DP_h_constr : {REC_sizing.DP_h_constr}")
        # print(f"DP_c_constr : {REC_sizing.DP_c_constr}")
        
        return RC, 0, "Fail"
    
    # REC_model.su_H.print_resume()
    # REC_model.su_C.print_resume()
    # print(f"Q_dot_constr : {REC_sizing.Q_dot_constr}")
    # print(f"DP_h_constr : {REC_sizing.DP_h_constr}")
    # print(f"DP_c_constr : {REC_sizing.DP_c_constr}")

    # --- GasHeater ---
    try:
        GH_model = RC.components['GasHeater'].model
        GH_sizing = RC.components['GasHeater'].sizing = ShellAndTubeSizingOpt()
    
        choice_vectors = {
            'D_o_inch': [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
            'Shell_ID_inch': [25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
            'Tube_pass': [1, 2, 4],
            'tube_layout': [0, 45, 60],
            'n_parallel': [1, 2, 3, 4, 5, 6],
        }
        
        GH_sizing.set_inputs(
            fluid_H=GH_model.su_H.fluid, T_su_H=GH_model.su_H.T,
            P_su_H=GH_model.su_H.p, m_dot_H=GH_model.su_H.m_dot,
            fluid_C=GH_model.su_C.fluid, T_su_C=GH_model.su_C.T,
            P_su_C=GH_model.su_C.p, m_dot_C=GH_model.su_C.m_dot,
        )
        
        GH_sizing.set_parameters(
            n_series=1, foul_t=0.000176, foul_s=0.000176, tube_cond=20,
            Overdesign=0, Shell_Side='H', Flow_Type='Shell&Tube',
            H_DP_ON=True, C_DP_ON=True, n_disc=50,
    
            opt_vars=['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'],
            T_max_cycle=RC.sources['GH_Water'].properties.T,
            p_max_cycle=RC.components['Pump'].model.ex.p,
    
            H_Corr={"SC": "Shell_Kern_HTC", "1P": "Shell_Kern_HTC", "2P": "Shell_Kern_HTC"},
            C_Corr={"SC": "Gnielinski", "1P": "Gnielinski", "2P": "Flow_boiling"},
            H_DP={"SC": "Shell_Kern_DP", "1P": "Shell_Kern_DP", "2P": "Shell_Kern_DP"},
            C_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Gnielinski_DP"},
    
            Q_dot=GH_model.Q.Q_dot,
            DP_h=max(GH_model.DP_h, 1e3),
            DP_c=max(GH_model.DP_c, 1e3),
        )
        
        bounds = {
            "L_shell": [1, 15],
            "D_o_inch": [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
            "Shell_ID_inch": [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]],
            "Tube_pass": [choice_vectors['Tube_pass'][0], choice_vectors['Tube_pass'][-1]],
            "tube_layout": [choice_vectors['tube_layout'][0], choice_vectors['tube_layout'][-1]],
            "Baffle_cut": [15, 45],
        }
        GH_sizing.set_bounds(bounds, choice_vectors=choice_vectors)
        
        GH_sizing.opt_size()
        
    except Exception as e:
        print(f"⚠️ Failed to design GasHeater: {e}")
        return RC, 0, "Fail"

    # --- Condenser ---
    try:
        CD_model = RC.components['Condenser'].model
        CD_sizing = RC.components['Condenser'].sizing = ShellAndTubeSizingOpt()
    
        choice_vectors = {
            'D_o_inch': [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
            'Shell_ID_inch': [25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
            'Tube_pass': [1, 2, 4],
            'tube_layout': [0, 45, 60],
            'n_parallel': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        }
    
        CD_sizing.set_inputs(
            fluid_H=CD_model.su_H.fluid, T_su_H=CD_model.su_H.T,
            P_su_H=CD_model.su_H.p, m_dot_H=CD_model.su_H.m_dot,
            fluid_C=CD_model.su_C.fluid, T_su_C=CD_model.su_C.T,
            P_su_C=CD_model.su_C.p, m_dot_C=CD_model.su_C.m_dot,
        )
    
        CD_sizing.set_parameters(
            n_series=1, foul_t=0.000176, foul_s=0.000176, tube_cond=20,
            Overdesign=0, Shell_Side='C', Flow_Type='Shell&Tube',
            H_DP_ON=True, C_DP_ON=True, n_disc=50,
    
            opt_vars=['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'],
            T_max_cycle=RC.sources['CD_Water'].properties.T,
            p_max_cycle=RC.components['Pump'].model.ex.p,
    
            H_Corr={"SC": "Gnielinski", "1P": "Gnielinski", "2P": "Thome_Condensation"},
            C_Corr={"SC": "Shell_Kern_HTC", "1P": "Shell_Kern_HTC", "2P": "Shell_Kern_HTC"},
            H_DP={"SC": "Gnielinski_DP", "1P": "Gnielinski_DP", "2P": "Choi_DP"},
            C_DP={"SC": "Shell_Kern_DP", "1P": "Shell_Kern_DP", "2P": "Shell_Kern_DP"},
    
            Q_dot=CD_model.Q.Q_dot,
            DP_h=max(CD_model.DP_h, 1e4),
            DP_c=max(CD_model.DP_c, 1e4),
        )
    
        bounds = {
            "L_shell": [1, 15],
            "D_o_inch": [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
            "Shell_ID_inch": [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]],
            "Tube_pass": [choice_vectors['Tube_pass'][0], choice_vectors['Tube_pass'][-1]],
            "tube_layout": [choice_vectors['tube_layout'][0], choice_vectors['tube_layout'][-1]],
            "Baffle_cut": [15, 45],
        }
        CD_sizing.set_bounds(bounds, choice_vectors=choice_vectors)
    
        CD_sizing.opt_size()
    
    except Exception as e:
        print(f"⚠️ Failed to design Condenser: {e}")
        return RC, 0, "Fail"

    # --- Pump ---
    try:
        Pump_model = RC.components['Pump'].model
        Pump_sizing = RC.components['Pump'].sizing = RadialPumpODDesign(RC.fluid)

        Pump_sizing.set_inputs(
            P_su=Pump_model.su.p, P_ex=Pump_model.ex.p, T_su=Pump_model.su.T,
            H1=0, H2=0, v1=0, v2=0, m_dot=Pump_model.su.m_dot,
        )

        Pump_sizing.set_parameters(
            Omega_choices=np.array([750, 1000, 1500, 3000]),
            n_parallel_choices=np.array([1, 2, 3, 4, 5, 6, 7, 8]),
        )
        Pump_sizing.design()

    except Exception as e:
        print(f"⚠️ Failed to design Pump: {e}")
        return RC, 0, "Fail"

    # --- Turbine ---
    eta_axial = 0
    eta_radial = 0

    try:
        if turb_choice != 'Radial':
            Turb_model = RC.components['Expander'].model
            Turb_axial_sizing = AxialTurbineMeanLineSizing(RC.fluid)

            Turb_axial_sizing.set_inputs(
                mdot=Turb_model.su.m_dot, W_dot=Turb_model.W.W_dot,
                p0_su=Turb_model.su.p, T0_su=Turb_model.su.T, p_ex=Turb_model.ex.p,
            )

            Turb_axial_sizing.set_parameters(
                Zweifel=0.8,
                damping=0.3, p_rel_tol=0.05, delta_tip=0.4e-3, N_lw=0, D_lw=0,
                e_blade=0.002e-3, t_TE_o=0.05, t_TE_min=5e-4,
            )
            
            Turb_axial_sizing.set_bounds(
                AR_min=0.8, r_hub_tip_max=0.95, r_hub_tip_min=0.6,
                Re_bounds=[3e6, 8e6], psi_bounds=[1, 1.9], phi_bounds=[0.5, 0.8],
                R_bounds=[0.45, 0.55], M_1st_bounds=[0.3, 0.6], r_m_bounds=[0.1, 0.6],
            )

            Turb_axial_sizing.sizing_parallel(n_jobs=-1, n_particles=30, max_iter=50)
            eta_axial = Turb_axial_sizing.eta_is

    except Exception as e:
        print(f"⚠️ Failed to design the axial Turbine: {e}")

    try:
        if turb_choice != 'Axial':
            Turb_model = RC.components['Expander'].model
            Turb_radial_sizing = RadialTurbineMeanLineSizing(RC.fluid)

            Turb_radial_sizing.set_inputs(
                mdot=Turb_model.su.m_dot, W_dot=Turb_model.W.W_dot,
                p0_su=Turb_model.su.p, T0_su=Turb_model.su.T, p_ex=Turb_model.ex.p,
            )

            Turb_radial_sizing.set_parameters(
                S_b4_ratio=1.05, t_TE_c_S_max=0.02, t_TE_S=5e-4,
                cl_a=0.4e-3, cl_r=0.4e-3, damping=0.5,
                Mth_target=0.4, r5t_guess=0.15, r4_guess=0.22,
            )
            
            Turb_radial_sizing.set_bounds(
                r5_r4_bounds=[0.3, 0.7], psi_bounds=[0.5, 1.5], phi_bounds=[0.3, 0.6],
                xhi_bounds=[0.3, 0.6], r5h_r5t_bounds=[0.3, 0.4],
            )

            Turb_radial_sizing.sizing_parallel(n_jobs=-1, n_particles=20, max_iter=50)
            eta_radial = Turb_radial_sizing.eta_is

    except Exception as e:
        print(f"⚠️ Failed to design the radial Turbine: {e}")

    if eta_axial == 0 and eta_radial == 0:
        return RC, 0, "Fail"
    else:
        if eta_axial > eta_radial:
            RC.components['Expander'].sizing = Turb_sizing = Turb_axial_sizing
            turb_choice = "Axial"
        else:
            RC.components['Expander'].sizing = Turb_sizing = Turb_radial_sizing
            turb_choice = "Radial"

    print(f"eta_axial : {eta_axial}")
    print(f"eta_radial : {eta_radial}")

    RC.CAPEX = {
        "Pump": np.round(Pump_sizing.CAPEX['Total']),
        "GasHeater": np.round(GH_sizing.CAPEX['Total']),
        "Recuperator": np.round(REC_sizing.CAPEX['Total']),
        "Expander": np.round(Turb_sizing.CAPEX['Total']),
        "Condenser": np.round(CD_sizing.CAPEX['Total']),
    }
    RC.CAPEX["Total"] = sum(v for k, v in RC.CAPEX.items() if k != "Total")

    return RC, 1, turb_choice

#%% Classe étendue : hérite de la brique d'optimisation importée

class CO2RCOptimizer(CO2RC_HX_optimizer):
    """
    Étend CO2RCOptimizer (co2_rc_pso_optimizer.py) avec :
    - dimensionnement des composants (TCO2_rec_comp_sizing)
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

    def evaluate_systems(self):
        RC_scores = []
        delta_dicts = []

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

            delta_dicts[-1]['DP_h_gh'] = delta_h_gh = ((np.max([DP_h_gh, self.params['DP_h_gh']]) - self.params['DP_h_gh']) / self.params['DP_h_gh']) ** 2
            delta_dicts[-1]['DP_c_gh'] = delta_c_gh = ((np.max([DP_c_gh, self.params['DP_c_gh']]) - self.params['DP_c_gh']) / self.params['DP_c_gh']) ** 2

            delta_dicts[-1]['DP_h_cond'] = delta_h_cond = ((np.max([DP_h_cond, self.params['DP_h_cond']]) - self.params['DP_h_cond']) / self.params['DP_h_cond']) ** 2
            delta_dicts[-1]['DP_c_cond'] = delta_c_cond = ((np.max([DP_c_cond, self.params['DP_c_cond']]) - self.params['DP_c_cond']) / self.params['DP_c_cond']) ** 2

            delta_dicts[-1]['DP_h_rec'] = delta_h_rec = ((np.max([DP_h_rec, self.params['DP_h_rec']]) - self.params['DP_h_rec']) / self.params['DP_h_rec']) ** 2
            delta_dicts[-1]['DP_c_rec'] = delta_c_rec = ((np.max([DP_c_rec, self.params['DP_c_rec']]) - self.params['DP_c_rec']) / self.params['DP_c_rec']) ** 2

            score_current = delta_exp + delta_pp + delta_h_gh + delta_c_gh + delta_h_cond + delta_c_cond + delta_h_rec + delta_c_rec
            RC_scores.append(score_current)

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

            self.current_RC, flag, turb_choice = TCO2_rec_comp_sizing(self.current_RC, self.turb_choice)

            if flag == 1:
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

        while self.criterion == 0 and it < n_it_max:

            self.allowable_positions = []  # reset à chaque itération

            # --- Optimisation PSO héritée du fichier 1 ---
            self.opt_RC(n_jobs=n_jobs, n_particles=n_particles, max_iter=max_iter,
                        patience=patience, ntop=ntop, init_pos=init_pos)

            # --- Dimensionnement des composants (fichier 2) ---
            self.size_components()

            new_params, best_score, delta_dict = self.evaluate_systems()
            self.new_params = new_params
            self.delta_dict = delta_dict

            print("\n----------------------------------------")
            print(f"New Values - Best Score : {best_score}")
            print("----------------------------------------")
            for k in new_params:
                print(f"{k} : {new_params[k]} - {delta_dict[k]*100}")

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

    T_test = 200 + 273.15
    n_MW = 10
    W_dot_obj = n_MW * 1e6
    eta_obj = 0.14

    Optimizer = CO2RCOptimizer('CO2')

    m_dot_HS_fact_bounds = [0.1, 3]
    m_dot_CS_fact_bounds = [2, 15]
    P_high_bounds = np.array([9, 200]) * 1e5
    m_dot_bounds = np.array([10, 80]) * n_MW

    eta_gh_disc = np.arange(0.9, 0.98, 0.02)
    PP_gh_disc = np.arange(1, 10, 1)
    eta_rec_disc = np.arange(0.6, 0.96, 0.02)
    PP_cd_disc = np.arange(1, 10, 1)

    Optimizer.set_parameters(
        save_file_path=None,   # ou un chemin, comme dans le fichier 2 d'origine
        RC_ARCH='REC',          # seule architecture compatible avec le sizing actuel
        eta_pp=0.85,
        DP_h_gh=100e3, DP_c_gh=4e5,
        
        PP_rec=0, DP_h_rec=4e5, DP_c_rec=2e5,
        eta_exp=0.92,
        
        SC_cd=0.1, DP_h_cond=2e5, DP_c_cond=100e3,
        P_high_bounds=P_high_bounds,
        m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
        m_dot_CS_fact_bounds=m_dot_CS_fact_bounds,
        m_dot_bounds=m_dot_bounds,
        eta_gh_disc=eta_gh_disc, PP_gh_disc=PP_gh_disc,
        eta_rec_disc=eta_rec_disc, PP_cd_disc=PP_cd_disc,
    )

    Optimizer.set_it_var(P_high=120e5, mdot=17, mdot_HS=20,
                          eta_gh=0.95, PP_gh=5, eta_rec=0.8, PP_cd=5, mdot_CS=200)

    Optimizer.set_obj(W_dot=W_dot_obj, eta=eta_obj)

    Optimizer.set_CSource(T=15 + 273.15, P=5e5, fluid='Water', m_dot=1000.0)
    Optimizer.set_HSource(T=T_test, P=100e5, fluid='Water', m_dot=50.0)

    Optimizer.cycle_design(ntop=5, n_particles=100, n_jobs=-1)

