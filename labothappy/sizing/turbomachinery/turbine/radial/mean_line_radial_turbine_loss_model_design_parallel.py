# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 13:31:47 2025

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, root, least_squares
from correlations.turbomachinery.radial_turbine_losses import nozzle_losses, rotor_losses

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyswarms as ps
import time

import warnings
warnings.filterwarnings("ignore")

#%%

# ---- joblib worker (process-based) ----
import os, numpy as np
from joblib import Parallel, delayed

from contextlib import contextmanager
from tqdm import tqdm
import joblib
from joblib.parallel import BatchCompletionCallBack

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

# --- worker for joblib ---
def _eval_particle(x, cls, fluid, params, inputs):
    """Evaluate one particle using a per-process cached solver."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")
    warnings.filterwarnings("ignore")

    global _SOLVER
    if _SOLVER is None:
        s = cls(fluid)
        s.set_parameters(**params)
        s.set_inputs(**inputs)
        _SOLVER = s

    # Re-apply inputs every call to avoid cross-particle contamination
    _SOLVER.set_inputs(**inputs)
    _SOLVER.W_dot = 0  # start clean for this evaluation

    x = np.asarray(x, dtype=float)
    cost = float(_SOLVER.design_system(x))
    return cost, float(_SOLVER.W_dot)

#%%

class RadialTurbineMeanLineDesign(object):

    def __init__(self, fluid):
        # Inputs
        self.inputs = {}
        
        # Params
        self.params = {}  

        # Abstract State 
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Blade Dictionnary
        self.stages = []

        # Velocity Triangle Data
        self.Vel_Tri_R = {}
        self.Vel_Tri_S = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None
        
        self.total_states  = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3,4,5])
        self.static_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3,4,5])
        self.AS = CP.AbstractState('HEOS', fluid)
            
        # Nozzle and rotor losses initiated to 0
        self.losses = { 
            'DP0_S_volute' : 0,
            'Dh_S_nozzle' : 0,
        }
        
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
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value
                
    # ---------------- Blade row ----------------------------------------------------------------------
     
    def designRotor(self):

        "1) -------- Velocity Triangle ------------------------"
        
        "1.1) -------- (4) Rotor Inlet ------------------------"
        # Rotor + Meridional velocities
        self.Vel_Tri_R['u4'] = u4 = np.sqrt(self.Dh0/self.inputs['psi'])
        self.Vel_Tri_R['vm5'] = vm5 = u4*self.inputs['phi']
        self.Vel_Tri_R['vm4'] = vm4 = vm5/self.inputs['xhi']

        # Absolute velocities        
        self.Vel_Tri_R['vu4'] = vu4 = self.inputs['psi']*self.Vel_Tri_R['u4']
        self.Vel_Tri_R['v4'] = np.sqrt(vu4**2 + vm4**2)
        self.Vel_Tri_R['alpha4'] = alpha4 = np.arctan(vu4/vm4)

        # Relative velocities        
        self.Vel_Tri_R['wu4'] = wu4 = vu4 - u4
        self.Vel_Tri_R['w4'] = w4 = np.sqrt(wu4**2 + vm4**2)
        self.Vel_Tri_R['beta4'] = beta4 = np.arctan(wu4/vm4)

        "1.2) -------- (5) Rotor Outlet ------------------------"
        
        # Rotor Velocity
        self.Vel_Tri_R['u5'] = u5 = u4*self.params['r5_r4_ratio']
        
        # Absolute velocities        
        self.Vel_Tri_R['alpha5'] = alpha5 = 0 # No outlet swirl
        self.Vel_Tri_R['vu5'] = vu5 = vm5*np.tan(alpha5)
        self.Vel_Tri_R['v5'] = np.sqrt(vm5**2 + vu5**2)

        # Relative velocities        
        self.Vel_Tri_R['wu5'] = wu5 = vu5 - u5
        self.Vel_Tri_R['w5'] = w5 = np.sqrt(wu5**2 + vm5**2)
        self.Vel_Tri_R['beta5'] = beta5 = np.arctan(wu5/vm5)
        
        "2) -------- Rotor Computation -----------------------------"

        "2.1) -------- (4) Rotor Inlet ------------------------"
        
        self.update_total_AS(CP.HmassSmass_INPUTS, self.total_states['H'][3], self.total_states['S'][3], 4)
        
        h4 = self.total_states['H'][4] - self.Vel_Tri_R['v4']**2 / 2
        
        self.update_static_AS(CP.HmassSmass_INPUTS, h4, self.total_states['S'][3], 4)
        
        self.M4 = self.Vel_Tri_R['v4']/self.static_states['A'][4]
        self.M4_rel = self.Vel_Tri_R['w4']/self.static_states['A'][4]
        
        self.A4 = self.inputs['mdot']/(self.Vel_Tri_R['vm4']*self.static_states['D'][4])
        
        alpha4_deg = self.Vel_Tri_R['alpha4']*180/np.pi
        self.n_blades_R = np.ceil((np.pi*(110-alpha4_deg)*np.tan(self.Vel_Tri_R['alpha4']))/30) # Jamieson-Glassman
        
        "2.2) -------- (5) Rotor Outlet ------------------------"

        def system_rotor(x): # Mass Balance solving
            # guess outlet state (imposing input p_ex) and radii   
            h5 = x[0]*1e5
            r5t = x[1]
            r4 = x[2]
            
            self.update_static_AS(CP.HmassP_INPUTS, h5, self.inputs['p_ex'], 5)
            
            # Compute flow area
            self.M5 = self.Vel_Tri_R['v5']/self.static_states['A'][5]
            self.M5_rel = M5_rel = self.Vel_Tri_R['w5']/self.static_states['A'][5]
            
            self.A5 = self.inputs['mdot']/(self.Vel_Tri_R['vm5']*self.static_states['D'][5])
            
            # Compute geomtry satisfying the mass balance
            r5h = r5t*self.params['r5h_r5t_ratio']
            
            self.params['r5'] = r5 = r4*self.params['r5_r4_ratio']
            self.params['r4'] = r4
            
            self.params['r5h'] = r5h
            self.params['r5t'] = r5t
            
            # Blockage computation
            beta_5h = np.arctan(r5h*np.tan(self.Vel_Tri_R['beta5'])/r5)
            beta_5t = np.arctan(r5t*np.tan(self.Vel_Tri_R['beta5'])/r5)
            
            # From Aungier rules for preliminary design 
            self.params['t4'] = t4 = 0.04*r4
            self.params['t5h'] = t5h = 0.04*r4
            self.params['t5t'] = t5t = 0.04*r4
            
            self.params['b4'] = b4 = self.A4/(2*np.pi*r4 - self.n_blades_R*t4)
            self.params['b5'] = b5 = (r5t - r5h)
            self.params['Lz'] = L_z = 1.5*self.params['b5']
            
            teh = t5h/np.cos(beta_5h)
            tet = t5t/np.cos(beta_5t)
            
            Abb = (r5t - r5h)*(tet+teh)/2
            BK5 = self.n_blades_R*Abb/(np.pi*(r5t**2 - r5h**2))
            
            self.AS.update(CP.HmassP_INPUTS, h5, self.inputs['p_ex'])
            gamma5 = self.AS.cpmass()/self.AS.cvmass()
            
            # Evaluate losses
            self.rotor_losses = rotor_losses(alpha4, beta4, beta5, b4, b5, self.params['cl_a'], self.params['cl_r'],
                                             gamma5, L_z, M5_rel, self.n_blades_R, r4, r5, r5h, r5t, t5t, u4, vm4, vm5, w4, w5)
            
            self.eta_blade_rotor = (self.total_states['H'][4] - h5)/(self.total_states['H'][4]  - (h5 - self.rotor_losses['Dh_tot']))

            self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.static_states['S'][4])         
            h5_is = self.AS.hmass()
            
            h5_new = h5_is + self.rotor_losses['Dh_tot']
            
            h05 = h5_new + self.Vel_Tri_R['w5']**2 / 2 
                        
            # Reconciliate with enthalpy guess
            self.update_static_AS(CP.HmassP_INPUTS, h5_new, self.inputs['p_ex'], 5)
                        
            self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)  

            f1 = (h5_new - h5)/h5
            f2 = (self.A5 - np.pi*(r5t**2 - r5h**2)*(1.0 - BK5))/self.A5
            f3 = ((r5**2) - 0.5*(r5t**2 + r5h**2))/r5**2
            
            return np.sum(np.array([f1, f2, f3])**2)

        x0 = [self.static_states['H'][4]*1e-5, self.params['r5t_guess'], self.params['r4_guess']]

        self.AS.update(CP.PQ_INPUTS, self.inputs['p_ex'], 0.5)

        h_sat = self.AS.hmass()

        bounds = [
            (h_sat*1.01*1e-5, self.static_states['H'][4] * 1e-5),
            (0.01, 1),
            (0.011, 1),
        ]

        self.sol_rotor = minimize(system_rotor, x0, method='L-BFGS-B', bounds=bounds,
                        options={'ftol': 1e-8, 'gtol': 1e-8})
        
        self.losses['Dh_R_incidence'] = self.rotor_losses['Dh_inc']
        self.losses['Dh_R_passage'] = self.rotor_losses['Dh_p']
        self.losses['Dh_R_clearance'] = self.rotor_losses['Dh_cl']
        self.losses['Dh_R_TE'] = self.rotor_losses['Dh_TE']
        self.losses['Dh_R_tot'] = self.rotor_losses['Dh_tot']
                
        return
        
    def designStator(self):

        "1) -------- (2) Volute Outlet ------------------------"
  
        p0loss_volute = self.losses['DP0_S_volute'] # asssumption
        
        p0_2 = self.total_states['P'][1] - p0loss_volute
        T0_2 = self.inputs['T0_su']

        self.update_total_AS(CP.PT_INPUTS, p0_2, T0_2, 2)        

        "2) -------- (3) Stator Outlet ------------------------"
        
        def system_MB_stator(x):
            h3 = x[0]*1e5
            s3 = x[1]*1e3
            rho_th = x[2]*1e2
            n_s = x[3]
                          
            # Params : From Aungier's rule
                     
            S3_c_ratio = 0.6
            theta_n = 0*np.pi/180 # deflection
            # d_c_ratio = 0.4 
            t2_c_ratio = 0.012
            t3_c_ratio = 0.025
            tmax_c_ratio = 0.06
            
            self.n_blades_S = n_s = round(n_s)
            
            self.params['pitch_S'] = S3 = 2*np.pi*self.params['r3']/n_s
            self.params['chord_S'] = c = S3/S3_c_ratio
            # d = d_c_ratio*c
            self.params['t2'] = t2_c_ratio*c
            self.params['t3'] = t3_c_ratio*c
            self.params['tmaxS'] = tmax = tmax_c_ratio*c
            
            self.update_static_AS(CP.HmassSmass_INPUTS, h3, s3, 3)

            p3 = self.static_states['P'][3]
            
            self.Vel_Tri_S['vu3'] = self.Vel_Tri_R['vu4']*self.params['r4']/self.params['r3']
            self.Vel_Tri_S['vm3'] = self.inputs['mdot']/(2*np.pi*self.params['r3']*self.params['b3']*self.static_states['D'][3])
            self.Vel_Tri_S['v3'] = v3 = np.sqrt(self.Vel_Tri_S['vm3']**2 + self.Vel_Tri_S['vu3']**2)
            self.Vel_Tri_S['alpha3'] = alpha3 = np.arctan(self.Vel_Tri_S['vu3']/self.Vel_Tri_S['vm3'])
    
            # Rodgers correlation for stator loss (1967)
            self.Re_3 = Re_3 = v3*self.params['b3']/self.static_states['V'][3]
            self.losses['Dh_S_nozzle'] = nozzle_losses(v3, Re_3, alpha3, S3, c, self.params['b3'])
            
            self.AS.update(CP.PSmass_INPUTS, p3, self.total_states['S'][2])

            h3_new = self.AS.hmass() + self.losses['Dh_S_nozzle']

            self.update_static_AS(CP.HmassP_INPUTS, h3_new, p3, 3)
            
            h03 = h3_new + (self.Vel_Tri_S['v3']**2)/2 
            
            self.update_total_AS(CP.HmassSmass_INPUTS, h03, self.static_states['S'][3], 3)
                
            "3) -------- (2-3) Stator Throat ------------------------"

            r_th = self.params['r3'] + c/2*np.sin(self.Vel_Tri_S['alpha3']-theta_n/2)          

            # r_th, rho_th
            alpha_th = np.arctan((self.params['r3']/r_th)*(rho_th/self.static_states['D'][3])*np.tan(self.Vel_Tri_S['alpha3']))
            # o_th = S3*np.cos(alpha_th)
            
            BK = n_s*tmax
            A_th = (2*np.pi*r_th - BK) * self.params['b3']
            v_th = self.inputs['mdot']/(rho_th*A_th)
            
            h_th = self.total_states['H'][3] - v_th**2 / 2

            self.AS.update(CP.HmassSmass_INPUTS, h_th, self.total_states['S'][3])

            a_th = self.AS.speed_sound()
            
            self.M_th = v_th/a_th
            
            f1 = ((h3 - h3_new)/h3_new)**2
            f2 = ((s3 - self.static_states['S'][3])/self.static_states['S'][3])**2
            f3 = ((rho_th - self.AS.rhomass())/self.AS.rhomass())**2
            f4 = ((self.M_th - self.params['Mth_target'])/self.params['Mth_target'])**2   # Mach at throat = target
            
            return np.sum(np.array([f1, f2, f3, 100*f4]))

        # Initial guess
        # x0 = [self.static_states['H'][4] * 1e-5, self.static_states['P'][4] * 1e-5 * 1.5,  self.static_states['D'][4] * 1e-2, self.n_blades_R + 3]
        x0 = [self.static_states['H'][4] * 1e-5, self.total_states['S'][2] * 1e-3,  self.static_states['D'][4] * 1e-2, self.n_blades_R + 3]
        
        # Bounds (in minimize, you need a sequence of (low, high) tuples)
        
        if self.total_states['S'][2] >= self.total_states['S'][4]:
            bounds = [
                (self.static_states['H'][4] * 1e-5, self.total_states['H'][2] * 1e-5),
                (self.total_states['S'][2] * 1e-3, self.total_states['S'][2] * 1e-3* 1.05),
                (self.static_states['D'][4] * 1e-2 * 0.5, self.static_states['D'][1] * 1e-2 * 2),
                (self.n_blades_R, self.n_blades_R * 2)
            ]
        else:
            # print(f"P4 : {self.static_states['P'][4]}")
            
            self.AS.update(CP.PQ_INPUTS, self.static_states['P'][4], 1)
            h_sat = self.AS.hmass() + 100
            
            bounds = [
                (h_sat*1e-5, self.total_states['H'][2] * 1e-5),
                (self.total_states['S'][2] * 1e-3, self.total_states['S'][4] * 1e-3),
                (self.static_states['D'][4] * 1e-2 * 0.5, self.static_states['D'][1] * 1e-2 * 2),
                (self.n_blades_R, self.n_blades_R * 2)
            ]
            
        # Call minimize (trust-constr works well with bounds, but L-BFGS-B is simpler)
        self.sol_stator1 = minimize(system_MB_stator, x0, method='L-BFGS-B', bounds=bounds,
                       options={'ftol': 1e-8, 'gtol': 1e-8})
        
        # # Check result
        # if sol.success:
        #     print("Solver succeeded!")
        # else:
        #     print("Solver failed.")
        
        # print("Message:", sol.message)
        # print("Solution:", sol.x)
        
        "4) -------- (2) Stator Inlet ------------------------"
        
        # blade angles based on Aungier's function for the blade shape
        a = 0.5*self.params['chord_S'] # Assumption as throat assumed at the middle of the blade
        b = self.params['tmaxS']
        c = self.params['chord_S']
        
        xhi2 = (180/np.pi) * np.arctan(4*b/(4*a - c))
        xhi3 = (180/np.pi) * np.arctan(4*b/(3*c - 4*a))

        # Minimize incidence
        fact = 3.6*np.sqrt(10*self.params['t2']/self.params['chord_S']) + abs(xhi3 - xhi2)/3.4
        i_opt = fact*np.sqrt(self.params['chord_S']/self.params['pitch_S']) - abs(xhi3 - xhi2)/2
        alpha_opt = (xhi2 - i_opt * np.sign(xhi3 - xhi2))*np.pi/180
        
        self.Vel_Tri_S['alpha2'] = alpha_opt
        self.params['r2'] = self.params['r3'] + self.params['chord_S']*np.sin((self.Vel_Tri_S['alpha3']+self.Vel_Tri_S['alpha2'])/2) 
        self.params['b2'] = self.params['b3']

        def stator_inlet_calc(x):
            rho2 = x[0]
            self.Vel_Tri_S['v2'] = v2 = self.inputs['mdot']/(2*np.pi*self.params['r2']*self.params['b2']*rho2)
            h2 = self.total_states['H'][1] - v2**2 / 2
            
            self.AS.update(CP.HmassSmass_INPUTS, h2, self.total_states['S'][1])
            rho2_calc = self.AS.rhomass()

            return np.array([rho2 - rho2_calc])

        x0 = [self.total_states['D'][1]]
        
        [rho_lo] = [self.static_states['D'][3]]
        [rho_hi] = [self.total_states['D'][1]]
        
        self.sol_stator2 = least_squares(stator_inlet_calc, x0,
                    bounds=([rho_lo],[rho_hi]),
                    method='trf', xtol=1e-8, ftol=1e-8, gtol=1e-8)

        h2 = self.total_states['H'][1] - self.Vel_Tri_S['v2']**2 / 2
        s2 = PropsSI('S', 'D', self.sol_stator2.x[0], 'H', h2, self.fluid)

        self.update_static_AS(CP.HmassSmass_INPUTS, h2, s2, 2)

        return
     
    def design_system(self, x):    
        
        self.inputs['psi'] = x[0]
        self.inputs['phi'] = x[1]
        self.inputs['xhi'] = x[2]
        self.params['r5_r4_ratio'] = x[3]
        self.params['r5h_r5t_ratio'] = x[4]
        
        self.update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        self.update_static_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        
        "------------- 1) Isentropic Expansion Calculation -----------------------------------------------" 
        s_in = self.total_states['S'][1]
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s_in)
        
        h_is_ex = self.AS.hmass()
        Dh0s = self.total_states['H'][1] - h_is_ex
                
        self.Dh0 = self.inputs['W_dot']/self.inputs['mdot']
        self.eta_is = self.Dh0/Dh0s
        
        def determine_stator_inlet(x):
            h03, p03 = x*1e5
            
            # print(f"old : {h03}, {p03}")
            
            "------------- 2) Rotor Design -------------------------------------"       
            self.update_total_AS(CP.HmassP_INPUTS, h03, p03, 3)
        
            self.designRotor()

            "------------- 3) Stator Design  ------------------------------------"         
            # Stator - rotor interspace
            self.params['r3'] = self.params['r4'] + self.params['S_b4_ratio'] * self.params['b4'] * np.cos(self.Vel_Tri_R['alpha4'])
            self.params['b3'] = self.params['b4']
            
            self.designStator()

            p03_new = self.total_states['P'][3]
            h03_new = self.total_states['H'][3]
            
            f1 = ((h03 - h03_new)/h03)**2
            f2 = ((p03 - p03_new)/p03)**2

            return f1 +  f2
            # return np.array([h03_new, p03_new])*1e-5    

        # res = 1
        # x_in = x0_disc
        
        # c = 0
        
        try:
            
            start_time = time.time()
            max_seconds = 20  # limit
            
            def time_limited_callback(xk, *args):
                if time.time() - start_time > max_seconds:
                    raise TimeoutError("Optimization exceeded time limit")
            
            # Initial guess vector
            x0_disc = np.concatenate(([self.total_states['H'][1]], [self.total_states['P'][1]]))*1e-5
            bounds_arr= np.array([(self.total_states['H'][1]-self.Dh0, 
                         self.total_states['H'][1]), 
                        (self.inputs['p_ex'], 
                         self.total_states['P'][1])])*1e-5


            sol = minimize(determine_stator_inlet, x0=x0_disc, bounds = bounds_arr, callback=time_limited_callback)     
            
            self.designRotor()
        
        except TimeoutError:
            obj = 1000
            # print(f"Fail time: {obj}")
            return obj
    
        except:
            obj = 1000
            # print(f"Fail err: {obj}")
            return obj

        hin = self.total_states['H'][1]
        hout = self.static_states['H'][5]
        
        self.AS.update(CP.PSmass_INPUTS, self.static_states['P'][5], self.static_states['S'][1])

        self.hout_s = self.AS.hmass()
        
        self.Omega = self.Vel_Tri_R['u4']/self.params['r4']
        self.W_dot = self.inputs['mdot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - self.hout_s)

        self.exit_loss = self.inputs['mdot']*(self.Vel_Tri_R['v5']**2)/2        
        
        if self.sol_rotor.success and self.sol_stator1.success and self.sol_stator2.success:
            if self.static_states['S'][5] <= self.static_states['S'][1] or self.static_states['S'][5] <= self.static_states['S'][3] or self.static_states['S'][3] <= self.static_states['S'][1]:
                obj = 1000
                # print(f"Fail entrop: {obj}")
                return obj  
            else:
                obj = -self.eta_is
        else:
            obj = 1000
            # print(f"Fail sol: {obj}")
            return obj        
        
        # print(f"obj : {obj}")
        
        return obj
    
#%%

    def design(self):
        bounds = (np.array([
            self.params['psi_bounds'][0],
            self.params['phi_bounds'][0],
            self.params['xhi_bounds'][0],
            self.params['r5_r4_bounds'][0],
            self.params['r5h_r5t_bounds'][0],
        ]),
        np.array([
            self.params['psi_bounds'][1],
            self.params['phi_bounds'][1],
            self.params['xhi_bounds'][1],
            self.params['r5_r4_bounds'][1],
            self.params['r5h_r5t_bounds'][1],
        ]))
    
        def objective_wrapper(x):
            rounded_x = np.copy(x)
            costs, wdots = [], []
            for xi in rounded_x:
                c = self.design_system(xi)
                costs.append(c)
            return np.asarray(costs, dtype=float)
    
        optimizer = ps.single.GlobalBestPSO(
            n_particles=10,
            dimensions=5,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        patience = 5
        tol = 1e-3
        max_iter = 2
        no_improve_counter = 0
        best_cost = np.inf
    
        for i in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = optimizer.swarm.best_cost
    
            print(f"--------------------------")
            print(f"Iteration: {i+1}/{max_iter}")
            print(f"Current best: {current_best}")
            print(f"--------------------------")

            # between-iteration W_dot raise
            batch_best = getattr(self, "_last_batch_max_wdot", self.inputs.get("W_dot", 0.0))
            if batch_best > self.inputs.get("W_dot", 0.0):
                self.inputs["W_dot"] = batch_best
                # print(f"[iter {i+1}] raised target W_dot to {self.inputs['W_dot']:.3f} W")
    
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
            # print(f"[{i+1}] Best cost: {best_cost:.6f}")
            if no_improve_counter >= patience:
                print("Stopping early due to stagnation.")
                break
    
        best_pos = optimizer.swarm.best_pos
        self.design_system(best_pos)
    
        print(f"Parameters : {self.inputs['psi'], self.inputs['phi'], self.inputs['xhi'], self.params['r5_r4_ratio'], self.params['r5h_r5t_ratio']}")
        print(f"Turbine rotation speed: {self.Omega} [RPM]")
        print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")
        print(f"Turbine Generation : {self.W_dot} [W]")
        return best_pos

    def design_parallel(self, n_jobs=-1, backend="loky", chunksize="auto"):
        os.environ["PYTHONWARNINGS"] = "ignore" 
        
        bounds = (np.array([
            self.params['psi_bounds'][0],
            self.params['phi_bounds'][0],
            self.params['xhi_bounds'][0],
            self.params['r5_r4_bounds'][0],
            self.params['r5h_r5t_bounds'][0],
        ]),
        np.array([
            self.params['psi_bounds'][1],
            self.params['phi_bounds'][1],
            self.params['xhi_bounds'][1],
            self.params['r5_r4_bounds'][1],
            self.params['r5h_r5t_bounds'][1],
        ]))
        
        dimensions = 5
        
        # snapshot of class + parameters
        inp = dict(self.inputs)
    
        def pick(*names, default=None):
            for n in names:
                if n in inp and inp[n] is not None:
                    return inp[n]
            return default
    
        # this dict is updated each iteration to carry the latest W_dot target
        inputs_snapshot = {
            "p0_su": pick("p0_su", "P0_su", "P_su", "p_su"),
            "T0_su": pick("T0_su", "t0_su", "T_su", "t_su"),
            "p_ex" : pick("p_ex", "P_ex"),
            "mdot" : pick("mdot", "m_dot"),
            "W_dot": pick("W_dot", "W"),
        }
    
        snapshot = {
            "cls": type(self),
            "fluid": self.fluid,
            "params": dict(self.params),
            "inputs": inputs_snapshot,
        }
        
        def objective_wrapper(X):
            X = np.asarray(X, dtype=float)
            with tqdm_joblib(tqdm(total=len(X), desc="Particles", unit="pt")):
                results = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
                    delayed(_eval_particle)(
                        xi, snapshot["cls"], snapshot["fluid"],
                        snapshot["params"], snapshot["inputs"]
                    ) for xi in X
                )
            # results: list of (cost, W_dot)
            costs = [r[0] for r in results]
            wdots = [r[1] for r in results]
            self._last_batch_max_wdot = max(wdots) if wdots else self.inputs.get("W_dot", 0.0)
            return np.asarray(costs, dtype=float)
    
        # --- PSO optimizer ---
        optimizer = ps.single.GlobalBestPSO(
            n_particles=20, dimensions=dimensions,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        patience, tol, max_iter = 5, 1e-3, 3
        no_improve, best_cost = 0, float("inf")
    
        for i in range(max_iter):
            # ensure workers see the current target THIS iteration
            snapshot["inputs"]["W_dot"] = self.inputs.get("W_dot", snapshot["inputs"]["W_dot"])
    
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            cur = optimizer.swarm.best_cost
    
            # --- between-iteration W_dot raise ---
            batch_best = getattr(self, "_last_batch_max_wdot", self.inputs.get("W_dot", 0.0))
            if batch_best > self.inputs.get("W_dot", 0.0):
                self.inputs["W_dot"] = batch_best
                # optional trace:
                # print(f"[iter {i+1}] raised target W_dot to {self.inputs['W_dot']:.3f} W")
    
            if cur < best_cost - tol:
                best_cost, no_improve = cur, 0
            else:
                no_improve += 1
            print(f"[{i+1}] Best cost: {best_cost:.6f}")
            if no_improve >= patience:
                print("Stopping early due to stagnation.")
                break
    
        best_pos = optimizer.swarm.best_pos
    
        # Finalize
        self.design_system(best_pos)
        # self.cost_estimation()
    
        print(f"Work Coef : {self.inputs['psi']}")
        print(f"Flow Coef : {self.inputs['phi']}")
        print(f"Xhi : {self.inputs['xhi']}")
        print(f"r5_r4_ratio : {self.params['r5_r4_ratio']}")
        print(f"r5h_r5t_ratio  : {self.params['r5h_r5t_ratio']}")
    
        print(f"P_in : {self.total_states['P'][0]} [Pa]")
        print(f"P_out: {self.static_states['P'][-1]} [Pa]")
        print(f"Omega: {self.params['Omega']} [RPM]")
        print(f"eta_is: {self.eta_is}")
        print(f"W_dot : {self.W_dot} [W]")
        
        print(f"r4 : {self.params['r4']} [m]")
        print(f"r5 : {self.params['r5']} [m]")
        
        return best_pos

Turb = RadialTurbineMeanLineDesign('CO2')

Turb.set_inputs(
    mdot = 100, # kg/s
    W_dot = 4.69*1e6, # W
    p0_su = 140*1e5, # Pa
    T0_su = 273.15 + 121, # K
    p_ex = 39.8*1e5, # Pa
    )

Turb.set_parameters(
    r5_r4_bounds = [0.3,0.7], # [-] : r5/r4 ratio
    psi_bounds = [0.5, 1.5],
    phi_bounds = [0.3, 0.6],
    xhi_bounds = [0.3, 0.6],
    r5h_r5t_bounds = [0.3, 0.4], # [-] : hub_tip ratio at the exit
    
    S_b4_ratio = 1.05, # flow path length to blade height ratio -> from 1 to 2 depending on the app, 1.05 max for CO2
    t_TE_c_S_max = 0.02, # [-]
    t_TE_S = 5*1e-4, # [m]
    cl_a = 0.4*1e-3, # [m] : Axial clearance
    cl_r = 0.4*1e-3, # [m] : Radial clearance
    
    damping = 0.5, # [-]

    Mth_target = 0.3, # [-]    
    r5t_guess = 0.15, # [m]
    r4_guess = 0.22, # [m]
    )
    
Turb.design_parallel()

