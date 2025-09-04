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

import warnings
warnings.filterwarnings("ignore")

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
            
            h05 = self.total_states['H'][4] - self.rotor_losses['Dh_tot']   
            
            # Reconciliate with enthalpy guess
            h5_new = self.AS.hmass() - self.Vel_Tri_R['w5']**2 / 2 
            self.update_static_AS(CP.HmassP_INPUTS, h5_new, self.inputs['p_ex'], 5)
    
            self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)  
            
            f1 = (h5_new - h5)/h5
            f2 = (self.A5 - np.pi*(r5t**2 - r5h**2)*(1.0 - BK5))/self.A5
            f3 = ((r5**2) - 0.5*(r5t**2 + r5h**2))/r5**2
            
            return np.sum(np.array([f1, f2, f3])**2)

        x0 = [self.static_states['H'][4]*1e-5, 0.15, 0.22]

        self.AS.update(CP.PQ_INPUTS, self.inputs['p_ex'], 0.5)

        h_sat = self.AS.hmass()

        bounds = [
            (h_sat*1.01*1e-5, self.static_states['H'][4] * 1e-5),
            (0.01, 1),
            (0.011, 1),
        ]

        self.sol_rotor = minimize(system_rotor, x0, method='L-BFGS-B', bounds=bounds,
                       options={'ftol': 1e-12, 'gtol': 1e-12})
        
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
            p3 = x[1]*1e5
            rho_th = x[2]*1e2
            n_s = x[3]
                        
            # print(x)
            
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
            
            self.update_static_AS(CP.HmassP_INPUTS, h3, p3, 3)
            
            self.Vel_Tri_S['vu3'] = self.Vel_Tri_R['vu4']*self.params['r4']/self.params['r3']
            self.Vel_Tri_S['vm3'] = self.inputs['mdot']/(2*np.pi*self.params['r3']*self.params['b3']*self.static_states['D'][3])
            self.Vel_Tri_S['v3'] = v3 = np.sqrt(self.Vel_Tri_S['vm3']**2 + self.Vel_Tri_S['vu3']**2)
            self.Vel_Tri_S['alpha3'] = alpha3 = np.arctan(self.Vel_Tri_S['vu3']/self.Vel_Tri_S['vm3'])
    
            # Rodgers correlation for stator loss (1967)
            self.Re_3 = Re_3 = v3*self.params['b3']/self.static_states['V'][3]
            self.losses['Dh_S_nozzle'] = nozzle_losses(v3, Re_3, alpha3, S3, c, self.params['b3'])
            
            h03 = self.total_states['H'][2] - self.losses['Dh_S_nozzle']
                        
            self.AS.update(CP.HmassSmass_INPUTS, h03, self.total_states['S'][2])
            p03 = self.AS.p()
            
            self.update_total_AS(CP.HmassP_INPUTS, self.total_states['H'][2], p03, 3)
            h3_new = h03 - (self.Vel_Tri_S['v3']**2)/2   
                
            self.update_static_AS(CP.HmassSmass_INPUTS, h3_new, self.total_states['S'][3], 3)
        
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
            f2 = ((p3 - self.static_states['P'][3])/self.static_states['P'][3])**2
            f3 = ((rho_th - self.AS.rhomass())/self.AS.rhomass())**2
            f4 = ((self.M_th - self.params['Mth_target'])/self.params['Mth_target'])**2   # Mach at throat = target
            
            return np.sum(np.array([f1, f2, f3, 100*f4]))

        # Initial guess
        x0 = [self.static_states['H'][4] * 1e-5, self.static_states['P'][4] * 1e-5 * 1.5,  self.static_states['D'][4] * 1e-2, self.n_blades_R + 3]
        
        # Bounds (in minimize, you need a sequence of (low, high) tuples)
        bounds = [
            (self.static_states['H'][4] * 1e-5, self.total_states['H'][2] * 1e-5),
            (self.static_states['P'][4] * 1e-5, self.static_states['P'][2] * 1e-5),
            (self.static_states['D'][4] * 1e-2 * 0.5, self.static_states['D'][1] * 1e-2 * 2),
            (self.n_blades_R, self.n_blades_R * 2)
        ]
            
        # Call minimize (trust-constr works well with bounds, but L-BFGS-B is simpler)
        self.sol_stator1 = minimize(system_MB_stator, x0, method='L-BFGS-B', bounds=bounds,
                       options={'ftol': 1e-10, 'gtol': 1e-10})
        
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
                    method='trf', xtol=1e-10, ftol=1e-10, gtol=1e-10)

        h2 = self.total_states['H'][1] - self.Vel_Tri_S['v2']**2 / 2
        s2 = PropsSI('S', 'D', self.sol_stator2.x[0], 'H', h2, self.fluid)

        self.update_static_AS(CP.HmassSmass_INPUTS, h2, s2, 2)

        return
     
    def design(self):    
        
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
                        
            "------------- 2) Rotor Design -------------------------------------"       
            # !!! : Guess on stator outlet
            self.update_total_AS(CP.HmassP_INPUTS, h03, p03,3)
        
            self.designRotor()

            "------------- 3) Stator Design  ------------------------------------"         
            # Stator - rotor interspace
            self.params['r3'] = self.params['r4'] + self.params['S_b4_ratio'] * self.params['b4'] * np.cos(self.Vel_Tri_R['alpha4'])
            self.params['b3'] = self.params['b4']
            
            self.designStator()

            p03_new = self.total_states['P'][3]
            h03_new = self.total_states['H'][3]

            return np.array([h03_new, p03_new])*1e-5

        # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])         
        
        # Initial guess vector
        x0_disc = np.concatenate(([self.total_states['H'][1]], [self.total_states['P'][1]]))*1e-5
        
        res = 1
        x_in = x0_disc
        
        c = 0
        
        while res > 1e-6:
            
            print(f"iteration {c+1}")
            
            if c > 100:
                exit()
            
            # print(f"x_in : {x_in}")
            
            x_out = determine_stator_inlet(x_in)

            # print(f"x_out : {x_out}")
            
            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            
            x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                          
            # print(f"new x_in : {x_in}")

            c += 1
            
            print(f"res : {res}")
                        
        determine_stator_inlet(x_out)

        hin = self.total_states['H'][1]
        hout = self.static_states['H'][5]
        
        self.AS.update(CP.PSmass_INPUTS, self.static_states['P'][5], self.static_states['S'][1])

        hout_s = self.AS.hmass()
        
        self.W_dot = self.inputs['mdot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - hout_s)

        self.exit_loss = self.inputs['mdot']*(self.Vel_Tri_R['v5']**2)/2        
        
        return

Turb = RadialTurbineMeanLineDesign('CO2')

Turb.set_inputs(
    mdot = 100, # kg/s
    W_dot = 4.69*1e6, # W
    p0_su = 140*1e5, # Pa
    T0_su = 273.15 + 121, # K
    p_ex = 39.8*1e5, # Pa
    psi = 1, # [-] : Iterate
    phi = 0.4, # [-] : Iterate
    xhi = 0.4, # [-] : Iterate
    )

Turb.set_parameters(
    r5_r4_ratio = 0.5, # [-] : Iterate
    r5h_r5t_ratio = 0.3, # [-] : Iterate
    S_b4_ratio = 1.05, # flow path length to blade height ratio -> from 1 to 2 depending on the app, 1.05 max for CO2
    Mth_target = 0.4, # [-]
    t_TE_c_S_max = 0.02, # [-]
    t_TE_S = 5*1e-4, # [m]
    cl_a = 0.4*1e-3, # [m] : Axial clearance
    cl_r = 0.4*1e-3, # [m] : Radial clearance
    damping = 0.5, # [-]
    )
    
Turb.design()

