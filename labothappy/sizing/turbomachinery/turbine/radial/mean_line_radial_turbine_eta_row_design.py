# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 13:31:47 2025

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, root, least_squares

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
        self.Dh_n = 0
        self.Dh_r = 0
        
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
        self.Vel_Tri_R['u4'] = np.sqrt(self.Dh0/self.inputs['psi'])
        self.Vel_Tri_R['vm5'] = self.Vel_Tri_R['u4']*self.inputs['phi']
        self.Vel_Tri_R['vm4'] = self.Vel_Tri_R['vm5']/self.inputs['xhi']

        # Absolute velocities        
        self.Vel_Tri_R['vu4'] = self.inputs['psi']*self.Vel_Tri_R['u4']
        self.Vel_Tri_R['v4'] = np.sqrt(self.Vel_Tri_R['vu4']**2 + self.Vel_Tri_R['vm4']**2)
        self.Vel_Tri_R['alpha4'] = np.arctan(self.Vel_Tri_R['vu4']/self.Vel_Tri_R['vm4'])

        # Relative velocities        
        self.Vel_Tri_R['wu4'] = self.Vel_Tri_R['vu4'] - self.Vel_Tri_R['u4']
        self.Vel_Tri_R['w4'] = np.sqrt(self.Vel_Tri_R['wu4']**2 + self.Vel_Tri_R['vm4']**2)
        self.Vel_Tri_R['beta4'] = np.arctan(self.Vel_Tri_R['wu4']/self.Vel_Tri_R['vm4'])

        "1.2) -------- (5) Rotor Outlet ------------------------"
        
        # Rotor Velocity
        self.Vel_Tri_R['u5'] = self.Vel_Tri_R['u4']*self.params['r5_r4_ratio']
        
        # Absolute velocities        
        self.Vel_Tri_R['alpha5'] = 0 # No outlet swirl
        self.Vel_Tri_R['vu5'] = self.Vel_Tri_R['vm5']*np.tan(self.Vel_Tri_R['alpha5'])
        self.Vel_Tri_R['v5'] = np.sqrt(self.Vel_Tri_R['vm5']**2 + self.Vel_Tri_R['vu5']**2)

        # Relative velocities        
        self.Vel_Tri_R['wu5'] = self.Vel_Tri_R['vu5'] - self.Vel_Tri_R['u5']
        self.Vel_Tri_R['w5'] = np.sqrt(self.Vel_Tri_R['wu5']**2 + self.Vel_Tri_R['vm5']**2)
        self.Vel_Tri_R['beta5'] = np.arctan(self.Vel_Tri_R['wu5']/self.Vel_Tri_R['vm5'])
        
        "2) -------- Rotor Computation -----------------------------"

        "2.1) -------- (4) Rotor Inlet ------------------------"
        self.update_total_AS(CP.HmassP_INPUTS, self.total_states['H'][1], self.total_states['P'][1], 4) # !!! Assumption on p
        
        h4 = self.total_states['H'][4] - self.Vel_Tri_R['v4']**2 / 2
        h4s = h4 - self.Dh_n
         
        self.update_static_AS(CP.HmassSmass_INPUTS, h4s, self.total_states['S'][1], 4)
        
        self.M4 = self.Vel_Tri_R['v4']/self.static_states['A'][4]
        self.M4_rel = self.Vel_Tri_R['w4']/self.static_states['A'][4]
        
        self.A4 = self.inputs['mdot']/(self.Vel_Tri_R['vm4']*self.static_states['D'][4])
        
        alpha4_deg = self.Vel_Tri_R['alpha4']*180/np.pi
        self.n_blades_R = np.ceil((np.pi*(110-alpha4_deg)*np.tan(self.Vel_Tri_R['alpha4']))/30) # Jamieson-Glassman
        
        "2.2) -------- (5) Rotor Outlet ------------------------"

        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.total_states['S'][1])
        h5 = self.AS.hmass() - self.Dh_n - self.Dh_r

        self.update_static_AS(CP.HmassP_INPUTS, h5, self.inputs['p_ex'], 5)

        h05 = self.static_states['H'][5] + self.Vel_Tri_R['v5']**2 / 2
        self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)        
        
        self.M5 = self.Vel_Tri_R['v5']/self.static_states['A'][5]
        self.M5_rel = self.Vel_Tri_R['w5']/self.static_states['A'][5]
        
        self.A5 = self.inputs['mdot']/(self.Vel_Tri_R['vm5']*self.static_states['D'][5])
        
        def system_MB_rotor(x): # Mass Balance solving
            r5t, r4 = x
            
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
            
            self.params['b4'] = self.A4/(2*np.pi*r4 - self.n_blades_R*t4)
            self.params['b5'] = (r5t - r5h)
            self.params['Lz'] = 1.5*self.params['b5']
            
            teh = t5h/np.cos(beta_5h)
            tet = t5t/np.cos(beta_5t)
            
            Abb = (r5t - r5h)*(tet+teh)/2
            BK5 = self.n_blades_R*Abb/(np.pi*(r5t**2 - r5h**2))
            
            f1 = self.A5 - np.pi*(r5t**2 - r5h**2)*(1.0 - BK5) 
            f2 = (r5**2) - 0.5*(r5t**2 + r5h**2)
            
            return np.array([f1, f2])

        x0 = [0.15, 0.22]

        sol = root(system_MB_rotor, x0)
        
        if not sol.success:
            raise RuntimeError(sol.message)
            
        return
        
    def designStator(self):

        "1) -------- (2) Volute Outlet ------------------------"
  
        p0loss_volute = 0 # asssumption
        
        p0_2 = self.total_states['P'][1] - p0loss_volute
        T0_2 = self.inputs['T0_su']

        self.update_total_AS(CP.PT_INPUTS, p0_2, T0_2, 2)        

        "2) -------- (3) Stator Outlet ------------------------"
        # !!!
        # ITERATION REQUIRED HERE TO INTERACT BETWEEN RHO3 AND VM3 (ADD LOSS CORR FOR S3)
        self.Vel_Tri_S['vu3'] = self.Vel_Tri_R['vu4']*self.params['r4']/self.params['r3']
        self.Vel_Tri_S['vm3'] = self.inputs['mdot']/(2*np.pi*self.params['r3']*self.params['b3']*self.total_states['D'][2]) # !!! Assumption for now
        self.Vel_Tri_S['v3'] = np.sqrt(self.Vel_Tri_S['vm3']**2 + self.Vel_Tri_S['vu3']**2)
        self.Vel_Tri_S['alpha3'] = np.arctan(self.Vel_Tri_S['vu3']/self.Vel_Tri_S['vm3'])

        h03 = self.total_states['H'][2]
        h3 = h03 - (self.Vel_Tri_S['v3']**2)/2   
    
        s3 = self.total_states['S'][2] # !!! Assumption for now

        self.update_total_AS(CP.HmassSmass_INPUTS, h03, s3, 3)
        self.update_static_AS(CP.HmassSmass_INPUTS, h3, s3, 3)
        
        "3) -------- (2-3) Stator Throat ------------------------"
        # Params : From Aungier's rule

        def system_MB_stator(x): # Mass Balance solving
            rho_th, n_s = x
                        
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
            
            r_th = self.params['r3'] + c/2*np.sin(self.Vel_Tri_S['alpha3']-theta_n/2)          

            # r_th, rho_th
            alpha_th = np.arctan((self.params['r3']/r_th)*(rho_th/self.static_states['D'][3])*np.tan(self.Vel_Tri_S['alpha3']))
            # o_th = S3*np.cos(alpha_th)
            
            BK = n_s*tmax
            A_th = (2*np.pi*r_th - BK) * self.params['b3']
            v_th = self.inputs['mdot']/(rho_th*A_th)
            
            h_th = h03 - v_th**2 / 2

            self.AS.update(CP.HmassSmass_INPUTS, h_th, s3)

            a_th = self.AS.speed_sound()
                        
            f1 = rho_th - self.AS.rhomass()
            f2 = v_th/a_th - self.params['Mth_target']        # Mach at throat = target
                        
            return np.array([f1, f2])

        x0 = [self.static_states['D'][3], self.n_blades_R + 3]
        
        [rho_lo, r_lo] = [self.static_states['D'][3], self.n_blades_R]
        [rho_hi, r_hi] = [self.static_states['D'][1], self.n_blades_R*2]
        
        sol = least_squares(system_MB_stator, x0,
                    bounds=([rho_lo, r_lo],[rho_hi, r_hi]),
                    method='trf', xtol=1e-10, ftol=1e-10, gtol=1e-10)
        
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
        
        sol = least_squares(stator_inlet_calc, x0,
                    bounds=([rho_lo],[rho_hi]),
                    method='trf', xtol=1e-10, ftol=1e-10, gtol=1e-10)

        h2 = self.total_states['H'][1] - self.Vel_Tri_S['v2']**2 / 2
        s2 = PropsSI('S', 'D', sol.x[0], 'H', h2, self.fluid)

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
        
        "------------- 2) Rotor Design -------------------------------------"         
        self.designRotor()

        "------------- 3) Stator Design  ------------------------------------"         
        # Stator - rotor interspace
        self.params['r3'] = self.params['r4'] + self.params['S_b4_ratio'] * self.params['b4'] * np.cos(self.Vel_Tri_R['alpha4'])
        self.params['b3'] = self.params['b4']
        
        self.designStator()

        self.exit_loss = self.inputs['mdot']*(self.Vel_Tri_R['v5']**2)/2        
        
        return

Turb = RadialTurbineMeanLineDesign('CO2')

Turb.set_inputs(
    mdot = 100, # kg/s
    W_dot = 4.69*1e6, # W
    p0_su = 140*1e5, # Pa
    T0_su = 273.15 + 121, # K
    p_ex = 39.8*1e5, # Pa
    psi = 2, # [-] : Iterate
    phi = 0.4, # [-] : Iterate
    xhi = 0.4, # [-] : Iterate
    )

Turb.set_parameters(
    r5_r4_ratio = 0.5, # [-] : Iterate
    r5h_r5t_ratio = 0.3, # [-] : Iterate
    S_b4_ratio = 1.05, # flow path length to blade height ratio -> from 1 to 2 depending on the app, 1.05 max for CO2
    Mth_target = 0.5, # [-]
    )
    
Turb.design()

