# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 13:31:47 2025

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, root_scalar

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
            
    def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.total_states['H'][position] = self.AS.hmass()            
        self.total_states['S'][position] = self.AS.smass()            
        self.total_states['P'][position] = self.AS.p()            
        self.total_states['D'][position] = self.AS.rhomass()            

        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1  
            
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
    
    def computeStator(self):
   
        p0loss_volute = 0
        
        p0_2 = self.total_states['P'][1] - p0loss_volute
        T0_2 = self.inputs['T0_su']

        self.update_total_AS(CP.PT_INPUTS, p0_2, T0_2, 2)
  
        h02 = self.total_states['H'][2]
        h2 = h02 - self.Vel_Tri_S['v2']**2 / 2
        s2 = self.total_states['S'][2]
        
        self.update_static_AS(CP.HmassSmass_INPUTS, h2, s2, 2)
        
        h03 = h02
        
        h3 = h03 - (self.Vel_Tri_S['v3']**2)/2            
        h3_s = h2 - (h2-h3)/self.eta_blade_row
        
        self.AS.update(CP.HmassSmass_INPUTS, h3_s, self.total_states['S'][2])
        p3 = self.AS.p()
        
        self.update_static_AS(CP.HmassP_INPUTS, h3, p3, 3)    
        self.update_total_AS(CP.HmassSmass_INPUTS, h03, self.static_states['S'][3], 3)
        
        return
     
    def computeRotor(self):

        self.update_total_AS(CP.HmassP_INPUTS, self.total_states['H'][3], self.total_states['P'][3], 4)
        self.update_static_AS(CP.HmassP_INPUTS, self.static_states['H'][3], self.static_states['P'][3], 4)
        
        h4 = self.static_states['H'][4]
        
        h05 = self.total_states['H'][4] - self.Dh0
        
        h5 = h05 - self.Vel_Tri_R['v5']**2 /2             
        h5_s = h4 - (h4-h5)/self.eta_blade_row
        
        self.AS.update(CP.HmassSmass_INPUTS, h5_s, self.total_states['S'][4])
        p5 = self.AS.p()
        
        self.update_static_AS(CP.HmassP_INPUTS, h5, p5, 5)        
        self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)        
        
        return          

    def computeVelTriangles(self):
        
        "1) -------- Rotor ------------------------"
        self.Vel_Tri_R['u4'] = np.sqrt(self.Dh0/self.inputs['psi'])
        self.Vel_Tri_R['vm5'] = self.Vel_Tri_R['u4']*self.inputs['phi']
        self.Vel_Tri_R['vm4'] = self.Vel_Tri_R['vm5']/self.inputs['xhi']

        # Rotor Geometry
        self.params['r4'] = self.Vel_Tri_R['u4']/self.omega_rads
        self.params['r5'] = self.params['r4']*self.params['r5_r4_ratio']

        self.Vel_Tri_R['alpha5'] = 0 # No outlet swirl
        self.Vel_Tri_R['vu5'] = self.Vel_Tri_R['vm5']*np.tan(self.Vel_Tri_R['alpha5'])
        self.Vel_Tri_R['v5'] = np.sqrt(self.Vel_Tri_R['vm5']**2 + self.Vel_Tri_R['vu5']**2)

        # Rotor Outlet
        self.Vel_Tri_R['u5'] = self.omega_rads*self.params['r5']
        self.Vel_Tri_R['wu5'] = self.Vel_Tri_R['vu5'] - self.Vel_Tri_R['u5']
        self.Vel_Tri_R['w5'] = np.sqrt(self.Vel_Tri_R['wu5']**2 + self.Vel_Tri_R['vm5']**2)
        self.Vel_Tri_R['beta5'] = np.arctan(self.Vel_Tri_R['wu5']/self.Vel_Tri_R['vm5'])
        
        # Euler Equation determines vt4
        self.Vel_Tri_R['vu4'] = (self.Dh0 + self.Vel_Tri_R['u5']*self.Vel_Tri_R['vm5'])/self.Vel_Tri_R['u4']
        self.Vel_Tri_R['alpha4'] = np.arctan(self.Vel_Tri_R['vu4']/self.Vel_Tri_R['vm4'])
        self.Vel_Tri_R['v4'] = np.sqrt(self.Vel_Tri_R['vu4']**2 + self.Vel_Tri_R['vm4']**2)

        self.Vel_Tri_R['wu4'] = self.Vel_Tri_R['vu4'] - self.Vel_Tri_R['u4']
        self.Vel_Tri_R['w4'] = np.sqrt(self.Vel_Tri_R['wu4']**2 + self.Vel_Tri_R['vm4']**2)
        self.Vel_Tri_R['beta4'] = np.arctan(self.Vel_Tri_R['wu4']/self.Vel_Tri_R['vm4'])

        "2) -------- Stator ------------------------"
        # Outlet
        self.Vel_Tri_S['alpha3'] = self.Vel_Tri_R['alpha4']
        self.Vel_Tri_S['beta3'] = self.Vel_Tri_R['beta4']
        self.Vel_Tri_S['vm3'] = self.Vel_Tri_R['vm4']
        self.Vel_Tri_S['u3'] = self.Vel_Tri_R['u4']
        self.Vel_Tri_S['vu3'] = self.Vel_Tri_R['vu4']
        self.Vel_Tri_S['wu3'] = self.Vel_Tri_R['wu4']
        self.Vel_Tri_S['v3'] = self.Vel_Tri_R['v4']
        self.Vel_Tri_S['w3'] = self.Vel_Tri_R['w4']
        
        # Inlet
        self.Vel_Tri_S['alpha2'] = 0 # design choice
        self.Vel_Tri_S['vm2'] = self.Vel_Tri_S['vm3']
        self.Vel_Tri_S['vu2'] = self.Vel_Tri_S['vm2']*np.tan(self.Vel_Tri_S['alpha2'])
        self.Vel_Tri_S['v2'] = np.sqrt(self.Vel_Tri_S['vm2']**2 + self.Vel_Tri_S['vu2']**2)
        
        self.Vel_Tri_S['wu2'] = self.Vel_Tri_S['vu2'] - self.Vel_Tri_S['u3']
        self.Vel_Tri_R['beta2'] = np.arctan(self.Vel_Tri_S['vu2']/self.Vel_Tri_S['vm2'])

        self.Vel_Tri_S['w2'] = np.sqrt(self.Vel_Tri_S['vm2']**2 + self.Vel_Tri_S['wu2']**2)
        
        return
     
    def design_system(self):    
                
        self.update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        self.update_static_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        
        "------------- 1) Isentropic Expansion Calculation -----------------------------------------------" 
        s_in = self.total_states['S'][1]
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s_in)
        
        h_is_ex = self.AS.hmass()
        Dh0s = self.total_states['H'][1] - h_is_ex
                
        self.Dh0 = self.inputs['W_dot']/self.inputs['mdot']
        self.eta_is = self.Dh0/Dh0s

        hout = self.total_states['H'][1] - Dh0s * self.eta_is
        
        self.AS.update(CP.HmassP_INPUTS, hout, self.inputs['p_ex'])
        rho_out = self.AS.rhomass()
        
        "------------- 2) Preliminary design parameters -----------------------------------------------" 
        
        if self.eta_is >= 0.87:
            self.omega_s = 0.55
        else:   
            def eta_is_func(omega_s):
                return 0.87 - 1.07*(omega_s - 0.55)**2 - 0.5*(omega_s - 0.55)**3
            
            def find_omega_s(eta_target, bounds=(0.1, 0.55)): # bounds can be inverted to 0.1 to 0.55 
                func = lambda omega: eta_is_func(omega) - eta_target
                sol = root_scalar(func, bracket=bounds, method='brentq')
                if sol.converged:
                    return sol.root
                else:
                    raise ValueError("No solution found in given bounds")
        
            self.omega_s = find_omega_s(self.eta_is)
        
        V_dot = self.inputs['mdot']/rho_out

        self.omega_rads = self.omega_s * (Dh0s)**0.75 /np.sqrt(V_dot)    
        self.Omega = self.omega_rads*60/(2*np.pi)
        
        self.v_s = 0.737*self.omega_s**0.2
        
        "------------- 3) Rotor Sizing -------------------------------------"
        
        # 3.1 Rotor Tip Sizing
        
        v_0s = np.sqrt(2*Dh0s) # discharge spouting velocity
        self.Vel_Tri_R['u4'] = self.v_s*v_0s
        self.Vel_Tri_R['vu4'] = self.Vel_Tri_R['u4']*self.eta_is/(2*self.v_s**2)
        
        h5 = self.static_states['H'][1] - self.Dh0
        self.params['r4'] = self.Vel_Tri_R['u4']/self.omega_rads
        
        Pt4_guess = self.total_states['P'][1]-self.static_states['D'][1]*Dh0s*(1-self.eta_is)/4

        # 3.2 Rotor Design Specifications
        self.Vel_Tri_R['alpha4'] = np.pi*(10.8+14.2*self.omega_s**2)/180
        
        self.params['tb4'] = 0.04*self.params['r4']
        self.params['tb5'] = 0.02*self.params['r4']
        self.params['rh5'] = 0.185*self.params['r4']
        
        self.Vel_Tri_R['vu5'] = 0 # assumed

        self.params['n_blade_R'] = np.round(12+0.03*(33 - self.Vel_Tri_R['alpha4']*180/np.pi),0)       


        "------------- 4) Velocity Triangle Computation -------------------------------------"         

        self.computeVelTriangles()
        
        "------------- 5) Find eta_blade_row and stator sizing ------------------------------------"         
        
        def find_eta_blade(x):
            print(x)
            self.eta_blade_row = x[0]
            self.computeStator()
    
            self.computeStator()
            self.computeRotor()
    
            pn_comp = self.static_states['P'][5]

            return (self.inputs["p_ex"] - pn_comp)**2
        
        sol = minimize(find_eta_blade, 1, bounds=[(self.eta_is-0.1, 1)], tol = 1e-4)
        
        self.exit_loss = self.inputs['mdot']*(self.Vel_Tri_R['v5']**2)/2        
        
        self.params['b4'] = self.inputs['mdot']/(2*np.pi*self.static_states['D'][4]*self.params['r4']*self.Vel_Tri_R['vm4'])
        self.params['b5'] = self.inputs['mdot']/(2*np.pi*self.static_states['D'][4]*self.params['r5']*self.Vel_Tri_R['vm5'])
        self.params['b3'] = self.params['b4']
        
        return

Turb = RadialTurbineMeanLineDesign('CO2')

Turb.set_inputs(
    mdot = 100, # kg/s
    W_dot = 4.69*1e6, # W
    p0_su = 140*1e5, # Pa
    T0_su = 273.15 + 121, # K
    p_ex = 39.8*1e5, # Pa
    psi = 1, # [-]
    phi = 0.4, # [-]
    xhi = 0.4, # [-]
    )

Turb.set_parameters(
    r5_r4_ratio = 0.5, # [-] : Maybe impose r5/r4 ratio instead - Iterate on it
    )
    
Turb.design_system()

