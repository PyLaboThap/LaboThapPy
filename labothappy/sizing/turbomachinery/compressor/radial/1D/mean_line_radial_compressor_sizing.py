#!/usr/bin/python3

# --- loading libraries 

from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize, brentq

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyswarms as ps

from labothappy.correlations.turbomachinery.radial_compressor_losses import radial_compressor_rotor_losses, radial_compressor_stator_losses

import warnings
warnings.filterwarnings("ignore")

class RadialCPMLDesign(object):

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
    
    # ---------------- Result Plot Methods ----------------------------------------------------------------

    def plot_geometry(self, fontsize = 16, ticksize = 12):
        
        r_m_line = np.ones(len(self.r_tip))*self.r_m
        
        x = np.linspace(0,len(self.r_tip)-1, len(self.r_tip))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1
        
        plt.figure()
        plt.plot(self.r_tip)
        plt.plot(self.r_hub)
        plt.plot(r_m_line)
                
        plt.axis([-0.5, len(self.r_tip)-0.5, 0, max(self.r_tip)*1.2])
        plt.legend(["$r_{tip}$", "$r_{hub}$", "$r_{m}$"])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.grid()
        plt.ylabel("Length or radius [m]", fontsize= fontsize)
        plt.show()

    def plot_n_blade(self, fontsize = 16, ticksize = 12):
        n_blade_plot = self.n_blade.flatten()

        x = np.linspace(0,len(n_blade_plot)-1, len(n_blade_plot))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1

        for i in range(len(n_blade_plot)*2):
              if np.mod(i,4) == 1 or np.mod(i,4) == 2: # Stator
                    n_blade_plot = np.insert(n_blade_plot,i,None)      

        n_blade_plot = n_blade_plot.reshape(int(len(n_blade_plot)/2),2)

        plt.figure()
        plt.plot(n_blade_plot[:, 0], 'o', label="Stator Blades")  # Plot first column with label
        plt.plot(n_blade_plot[:, 1], 'o', label="Rotor Blades")  # Plot second column with label
        plt.axis([-0.5, len(self.r_tip)-0.5, 0, max(n_blade_plot.flatten())*1.2])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.legend()
        plt.grid()
        plt.ylabel("Blade number [-]", fontsize= fontsize)
        plt.show()

    def plot_radius_verif(self, fontsize = 16, ticksize = 12):
        
        x = np.linspace(0,len(self.r_ratio2)-1, len(self.r_ratio2))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1
            
        plt.figure()
        plt.plot(self.r_ratio2)
        plt.axis([-0.5, len(self.r_ratio2)-0.5, 0, max(self.r_ratio2)*1.2])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.grid()
        plt.ylabel("$\\left[ r_{ext}/r_{hub} \\right]^2$ [-]", fontsize= fontsize)
        plt.show()

        plt.figure()
        plt.plot(self.r_hub_tip)
        plt.axis([-0.5, len(self.r_hub_tip)-0.5, 0, 1])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.grid()
        plt.ylabel("$\\left[ r_{hub}/r_{tip} \\right]$ [-]", fontsize= fontsize)
        plt.show()

    def plot_Mollier(self, fontsize = 16, ticksize = 12):
        # Thermo Prop
        x = np.linspace(0,len(self.r_tip)-1, len(self.r_tip))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1
        
        x2 = np.linspace(0,len(self.r_tip), len(self.r_tip)+1)
        labels2 = ['0'] + labels
        
        p = [self.stages[0].static_states['P'][1]]
        s = [self.stages[0].static_states['S'][1]]
        h = [self.stages[0].static_states['H'][1]]
        
        for i in range(self.nStages):
            p.append(self.stages[i].static_states['P'][2])
            p.append(self.stages[i].static_states['P'][3])
        
            s.append(self.stages[i].static_states['S'][2])
            s.append(self.stages[i].static_states['S'][3])
            
            h.append(self.stages[i].static_states['H'][2])
            h.append(self.stages[i].static_states['H'][3])
        
        plt.figure()
        plt.plot(np.array(p)*1e-3)
        plt.axis([-0.5, len(self.r_tip)+0.5, 0, max(np.array(p)*1e-3)*1.2])
        plt.xticks(ticks=x2, labels=labels2, size=ticksize)
        plt.grid()
        plt.ylabel("Oulet Pressure [kPa]", fontsize= fontsize)
        plt.show()
        
        plt.figure()
        plt.plot(s, h)
        plt.plot([s[0], s[0]], [h[0], h[-1]])
        
        # Define entropy range (in J/(kg·K))
        entropy_range = np.linspace(s[0], s[-1], 100)  # Adjust range for your fluid
        
        for P in p:
            enthalpy = [PropsSI('H', 'S', s, 'P', P, self.fluid) for s in entropy_range]  # Enthalpy in kJ/kg
            entropy = entropy_range  # Entropy in kJ/(kg·K)
            plt.plot(entropy, enthalpy, color = 'grey', alpha=0.3, label=f'P = {P/1e5} bar')
        
        plt.ylabel("$Enthalpy$ [J/kg]", fontsize= fontsize)
        plt.xlabel("$Entropy$ [J/(kg x K)]", fontsize= fontsize)
        plt.legend(["real", "isentropic"])
        plt.show()


    # ---------------- Flow Computations ------------------------------------------------------------------

    def computeVelTriangle(self):

        # Velocities over u
        self.Vel_Tri['vu2OverU'] = (2*(1-self.inputs['R']) + self.inputs['psi'])/2
        self.Vel_Tri['vu3OverU'] = (2*(1-self.inputs['R']) - self.inputs['psi'])/2
        self.Vel_Tri['vmOverU']  = self.inputs['phi']
        
        self.Vel_Tri['wu2OverU']  = self.Vel_Tri['vu2OverU'] - 1
        self.Vel_Tri['wu3OverU']  = self.Vel_Tri['vu3OverU'] - 1

        self.Vel_Tri['v2OverU']  = np.sqrt(self.Vel_Tri['vu2OverU']*self.Vel_Tri['vu2OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['w2OverU']  = np.sqrt(self.Vel_Tri['wu2OverU']*self.Vel_Tri['wu2OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['v3OverU']  = np.sqrt(self.Vel_Tri['vu3OverU']*self.Vel_Tri['vu3OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['w3OverU']  = np.sqrt(self.Vel_Tri['wu3OverU']*self.Vel_Tri['wu3OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])

        # Angles in radians
        self.Vel_Tri['alpha1'] = self.Vel_Tri['alpha3'] = np.arctan(self.Vel_Tri['vu3OverU']/self.Vel_Tri['vmOverU'])
        self.Vel_Tri['alpha2'] = np.arctan(self.Vel_Tri['vu2OverU']/self.Vel_Tri['vmOverU'])

        self.Vel_Tri['beta1'] = self.Vel_Tri['beta3'] = np.arctan(self.Vel_Tri['wu3OverU']/self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta2'] = np.arctan(self.Vel_Tri['wu2OverU']/self.Vel_Tri['vmOverU'])
        
        return 
    
    def computeBladeRow(self,stage,row_type):
        if row_type == 'S': # Stator
            hin = stage.static_states['H'][2]
            h0in = hin + (self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2  

            stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
            
            hout = h0in - (self.Vel_Tri['vu3']**2 + self.Vel_Tri['vm']**2)/2            
            hout_s = hin + (hout-hin)*self.eta_blade_row
            
            self.AS.update(CP.HmassSmass_INPUTS, hout_s, stage.total_states['S'][2])
            pout = self.AS.p()
            
            stage.update_static_AS(CP.HmassP_INPUTS, hout, pout, 3)            
                        
        else: # Rotor
            hin = stage.static_states['H'][1]
            h0in = hin + (self.Vel_Tri['wu1']**2 + self.Vel_Tri['vm']**2)/2  
                        
            stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
            
            hout = h0in - (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2            
            hout_s = hin + (hout-hin)*self.eta_blade_row
            
            self.AS.update(CP.HmassSmass_INPUTS, hout_s, stage.total_states['S'][1])
            pout = self.AS.p()
            
            stage.update_static_AS(CP.HmassP_INPUTS, hout, pout, 2)        
        
        return
            
    def computeRepeatingStages(self):
        
        for i in range(int(self.nStages)):
            if i == 0:
                self.computeBladeRow(self.stages[i], 'R')
                self.computeBladeRow(self.stages[i], 'S')
            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
                self.computeBladeRow(self.stages[i], 'R')
                self.computeBladeRow(self.stages[i], 'S')
            
        return
    
    def designRotor(self):
        
        "R0) Initiate Rotor Design"
        
        # Determine u2 from psi and Dh0s
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.total_states['S'][1])
        
        h0s = self.AS.hmass()
        self.Dh0s = h0s - self.total_states['H'][1]
        
        self.Vel_Tri_R['u2'] = u2 = np.sqrt(self.Dh0s/self.inputs['psi_is'])
        
        # Determine r2 from u2 and Omega
        self.params['r2'] = r2 = u2/self.omega
        
        # Flow coefficient : THE DESIGN OF A FAMILY OF PROCESS COMPRESSOR STAGES - H. Hazby – M. Casey – C. Robinson – R. Spataro
        self.Q = self.inputs['mdot']/self.total_states['D'][1] # volume flowrate based on inlet toal prop
        self.phi = 4*self.Q/(np.pi*u2*(2*r2)**2)
        
        "R1) Rotor Inlet"
        
        # Determine r1s, r1h, pitch1 and b2 from r2
        self.params['r1s'] = r1s = r2*self.inputs['r1s_r2']
        self.params['r1h'] = r1h = r1s*self.inputs['r1h_r1s']
        self.params['b1'] = b1 = r1s - r1h
        
        self.params['r1'] = r1 = np.sqrt((r1s**2 + r1h**2) / 2)

        self.params['pitch1'] = pitch1 = 2*np.pi*r1/(self.n_blade_R)
        self.params['b2'] = b2 = r2*self.inputs['b2_r2']

        self.params['r1'] = r1 = (r1s + r1h)/2
        self.Vel_Tri_R['u1'] = u1 = self.omega*r1
        
        u1h = self.omega*r1h
        u1s = self.omega*r1s
        
        # Rotor inlet velocity triangle
        self.A1 = b1*pitch1*self.n_blade_R
        self.o1 = pitch1*np.cos(np.pi/180*self.inputs['xhi1']) - self.params['t_b']
        
        self.impeller_inlet_blockage = 1 - self.n_blade_R*self.params['t_b']/(np.pi*r1*2*np.sin(np.pi*abs(self.inputs['xhi1'])/180)) 
        self.A1_th = np.pi*(r1s**2 - r1h**2)*self.impeller_inlet_blockage

        def compute_h1_new(h1): 
            
            self.update_static_AS(CP.HmassSmass_INPUTS, h1, self.total_states['S'][1], 1) 
                        
            self.Vel_Tri_R['vm1'] = vm1 = self.inputs['mdot']/(self.static_states['D'][1]*self.A1_th)
            
            self.Vel_Tri_R['beta1'] = beta1 = 0 # np.arctan(u1/vm1)
            self.Vel_Tri_R['w1'] = w1 = vm1/np.cos(beta1)
            self.Vel_Tri_R['wu1'] = w1*np.sin(beta1)
            
            self.Vel_Tri_R['vu1'] = vu1 = self.Vel_Tri_R['wu1'] + self.Vel_Tri_R['u1']
            self.Vel_Tri_R['v1'] = v1 = np.sqrt(vm1**2 + vu1**2)
            self.Vel_Tri_R['alpha1'] = np.arccos(vm1/v1)
            
            self.beta1h = np.arctan(u1h/vm1)
            self.beta1s = np.arctan(u1s/vm1)
            
            self.w1s = vm1/np.cos(self.beta1s)
            self.w1h = vm1/np.cos(self.beta1h)
                        
            h1_new = (self.total_states['H'][1] - v1**2 / 2) 
                        
            return h1_new


        def residual_h1(h1):
            
            # print(f"x : {h1}")
            res = h1 - compute_h1_new(h1)
            
            # print(f"res : {res}")

            return res
        
        # Bounds for brentq
        h1_min = self.total_states['H'][1] * 0.95   # large kinetic energy → low static h
        h1_max = self.total_states['H'][1]           # zero velocity → h1 = h01
        
        h1_solution = brentq(residual_h1, h1_min, h1_max, xtol=1e-6)
                
        "R2) Rotor Outlet"
        
        # Rotor Outlet Area
        self.params['pitch2'] = pitch2 = 2*np.pi*r2/self.n_blade_R
        self.A2_th = b2*pitch2*self.n_blade_R
        
        self.A2_th = (2*np.pi*r2*b2-self.n_blade_R*b2*self.params['t_b'])*np.cos(np.pi*abs(self.inputs['xhi2'])/180)
        
        # Slip Factor
        self.sigma = 1 - np.sqrt(np.cos(self.inputs['xhi2']*np.pi/180))/(self.n_blade_R)**0.7 # Aungier
        angle_star_deg = 19 + 0.2*(90 - self.inputs['xhi2'])
        sigma_star = np.sin(np.pi/180 * angle_star_deg)
        
        r1_r2_lim = (self.sigma-sigma_star) / (1-sigma_star)
        
        if r1/r2 > r1_r2_lim:
            num = (r1/r2) - r1_r2_lim**(np.sqrt((90 - self.inputs['xhi2'])/10))
            fact = 1 - (num/(1-sigma_star))
            self.sigma = self.sigma*fact
        
        def system_rotor(vm2): # Mass Balance solving

            self.Vel_Tri_R['vm2'] = vm2
        
            rho2 = self.inputs['mdot']/(vm2*self.A2_th) # vm2 <-> rho2
            
            self.Vel_Tri_R['vu2'] = vu2 = self.sigma*self.Vel_Tri_R['u2'] - vm2*np.tan(self.inputs['xhi2']*np.pi/180) # vu2 <-> vm2
            
            self.dh0 = vu2*self.Vel_Tri_R['u2'] - self.Vel_Tri_R['vu1']*self.Vel_Tri_R['u1']
            
            h02 = self.total_states['H'][1] + self.dh0 # vu2 <-> h02
            
            self.Vel_Tri_R['v2'] = v2 = np.sqrt(vu2**2 + vm2**2)

            h2 = h02 - v2**2 /2
            
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

            self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][2], self.total_states['S'][1]) # get p2 ? 
            h_is = self.AS.hmass()
            
            self.Vel_Tri_R['alpha2'] = alpha2 = np.arccos(vm2/v2)
            
            # self.params['L_z'] = L_z = self.params['r2']*(0.08 + 3.16*self.phi) 
            # From Aungier, R. H., (2000). Centrifugal compressors – a strategy for aerodynamic design and analysis, ASME Press, New York, USA 
            
            if "L_z" not in self.params:
                self.params['L_z'] = L_z = self.params['r2']*(0.1 + 2*self.phi) 
                # From THE DESIGN OF A FAMILY OF PROCESS COMPRESSOR STAGES - H. Hazby – M. Casey – C. Robinson – R. Spataro (2017)  
            else:
                L_z = self.params['L_z'] 
                
            self.Vel_Tri_R['wu2'] = wu2 =  self.Vel_Tri_R['vu2'] - self.Vel_Tri_R['u2']
            
            self.Vel_Tri_R['w2'] = w2 = np.sqrt(wu2**2 + vm2**2)
            
            self.Vel_Tri_R['beta2'] = beta2 = np.arccos(vm2/w2)
            
            self.rotor_losses = radial_compressor_rotor_losses(
                A1 = self.A1, A1_th = self.A1_th, alpha2 = alpha2, beta1 = self.Vel_Tri_R['beta1'], beta1h = self.beta1h, 
                beta1s = self.beta1s, beta2 = self.Vel_Tri_R['beta2'],  b2 = self.params['b2'], C_df = 0.004, C_fi = 0.004, Dh0 = self.dh0, eps_a = self.params['eps_imp'],
                eps_b = self.params['eps_bf_imp'], eps_r = self.params['eps_imp'], k_roughness = self.params['k_imp'], L_z = L_z, mdot = self.inputs['mdot'], 
                mu1 = self.static_states['V'][1], mu2 = self.static_states['V'][2], n_bl_r = self.n_blade_R, rho1 = self.static_states['D'][1], 
                rho2 = self.static_states['D'][2], r1h = self.params['r1h'], r1s = self.params['r1s'], r2 = self.params['r2'], u2 = self.Vel_Tri_R['u2'], 
                vu2 = vu2, v1m = self.Vel_Tri_R['vm1'], v2 = v2, w1 = self.Vel_Tri_R['w1'], w1_th = self.Vel_Tri_R['w1'], w1h = self.w1h,
                w1s = self.w1s, w2 = w2, xhi1 = self.inputs['xhi1']*np.pi/180, xhi2 = self.inputs['xhi2']*np.pi/180)
            
            h02_new = self.rotor_losses['tot'] + h_is
            h2_new = h02_new - v2**2 /2
            
            self.update_static_AS(CP.HmassSmass_INPUTS, h2_new, self.static_states['S'][2], 2)
            
            try:
                s2 = self.static_states['S'][2]
                self.update_total_AS(CP.HmassSmass_INPUTS, h02, s2, 2)
            except Exception:
                return 1e6
            
            res = (h2 - h2_new)/h2_new
            self.res_rotor_ex = res
            
            return res
        
        # Lower bound: vm2 that gives rho2 = 3x inlet density (physically motivated)
        rho2_max_phys = 3.0 * self.static_states['D'][1]
        vm2_min = self.inputs['mdot'] / (rho2_max_phys * self.A2_th)
        
        # Upper bound: vm2 that gives rho2 = inlet density (minimum compression)
        vm2_max = self.inputs['mdot'] / (self.static_states['D'][1] * self.A2_th)
        
        # Safety check: ensure sign change exists before calling brentq
        res_min = system_rotor(vm2_min)
        res_max = system_rotor(vm2_max)
        
        if res_min * res_max > 0:
            raise ValueError()  # no root in interval — penalize this design point
        
        self.sol_rotor_ex = brentq(system_rotor, vm2_min, vm2_max, xtol=1e-6)
        
        if self.rotor_losses['tot'] > self.dh0:
            raise ValueError()
        
        "Compute constraint terms"
        
        # Mach relatif au shroud de l'inducteur (contrainte M1s_rel ≤ 1.4)
        self.M1s_rel = self.w1s / self.static_states['A'][1]
        
        # Mach relatif moyen de l'inducteur (contrainte M1_rel ≤ 0.9)
        w1_mean = self.Vel_Tri_R['w1']
        self.M1_rel = w1_mean / self.static_states['A'][1]
        
        # Ratio de diffusion W2/W1s (contrainte ≥ 0.25)
        self.W2_W1s = self.Vel_Tri_R['w2'] / self.w1s
        
        vu1 = self.Vel_Tri_R['vu1']
        vu2 = self.Vel_Tri_R['vu2']
        u2  = self.Vel_Tri_R['u2']
        self.DR = 1.0 - (vu2 + vu1) / (2.0 * u2)
        
        return
    
    def designStator(self):

        "S3) Vaneless Space Exhaust"
        self.static_states['P'][3] = p3 = self.static_states['P'][2] + self.params["CP"]*(self.total_states['P'][2]-self.static_states['P'][2])
        self.total_states['H'][3] = h03 = self.total_states['H'][2]
        
        def vaneless_system(x):
        
            s3 = x
        
            self.update_total_AS(CP.HmassSmass_INPUTS, self.total_states['H'][3], s3, 3)
            self.update_static_AS(CP.PSmass_INPUTS, self.static_states['P'][3], s3, 3)
        
            p03 = self.total_states['P'][3]
            K = (self.total_states['P'][2] - p03)/(self.total_states['P'][2] - self.static_states['P'][2])

            CP_id = self.params['CP'] + K 
        
            self.Vel_Tri_S['v3'] = v3 = np.sqrt(2*(self.total_states['H'][3] - self.static_states['H'][3]))
            self.Vel_Tri_S['alpha3'] = self.Vel_Tri_R['alpha2']

            self.Vel_Tri_S['vm3'] = vm3 = v3*np.cos(self.Vel_Tri_S['alpha3'])
            self.Vel_Tri_S['vu3'] = vu3 = v3*np.sin(self.Vel_Tri_S['alpha3'])

            self.Vel_Tri_S['u3'] = u3 = self.Vel_Tri_R['u2']*self.inputs['r3_r2']

            self.Vel_Tri_S['wu3'] = wu3 = vu3-u3
            self.Vel_Tri_S['wm3'] = wm3 = vm3
            self.Vel_Tri_S['w3'] = np.sqrt(wm3**2 + wu3**2)
            self.Vel_Tri_S['beta3'] = beta3 = np.arctan(wu3/wm3)
            
            res = (self.inputs['mdot'] - (vm3*self.static_states['D'][3]*self.A2_th)*(1-CP_id)**(-0.5))**2
            
            return res
            
        self.AS.update(CP.HmassP_INPUTS, h03, p3)
        s3_max = self.AS.smass()  # entropy where h3 = h03, v3 = 0
        s3_min = self.static_states['S'][2]        
        
        sol = minimize(vaneless_system, s3_min, method='L-BFGS-B', bounds=[(s3_min,s3_max)],
                       options={'ftol': 1e-8, 'gtol': 1e-8})
    
        "S4) Vaned Diffuser Inlet"
        s4 = self.static_states['S'][3]
        h04 = self.total_states['H'][3]
        
        self.update_total_AS(CP.HmassSmass_INPUTS, h04, s4, 4)
        self.params['r3'] = self.params['r2']*self.inputs['r3_r2']
        self.Vel_Tri_R['u3'] = self.Vel_Tri_R['u2']*self.inputs['r3_r2']
        self.params['b3'] = self.params['b2'] *  self.params['b3_b2_ratio']

        self.params['pitch3'] = pitch3 = 2*np.pi*self.params['r3'] / self.n_blade_R  # n_blade_S = n_blade_R
        self.A3 = pitch3 * self.params['b2'] *  self.n_blade_R 
        
        def vaned_diffuser_inlet_system(x):

            h4 = x
            self.update_static_AS(CP.HmassSmass_INPUTS, h4, s4, 4)

            self.Vel_Tri_S['v4'] = v4 = np.sqrt(2*(self.total_states['H'][4] - h4))
            self.Vel_Tri_S['vu4'] = vu4 = self.Vel_Tri_S['vu3'] # Since r4 = r3
            self.Vel_Tri_S['vm4'] = vm4 = np.sqrt(max(v4**2 - vu4**2,0))
            self.Vel_Tri_S['alpha4'] = alpha4 = np.arctan(vu4/vm4)
            
            self.Vel_Tri_S['u4'] = u4 = self.Vel_Tri_R['u3']

            self.Vel_Tri_S['wu4'] = wu4 = vu4-u4
            self.Vel_Tri_S['wm4'] = wm4 = vm4
            self.Vel_Tri_S['w4'] = np.sqrt(wm4**2 + wu4**2)
            self.Vel_Tri_S['beta4'] = beta4 = np.arctan(wu4/wm4)
            
            res = (vm4 - self.inputs['mdot']/(self.static_states['D'][4]*self.A3))
            
            self.res_diff_in = res
            
            return res
        
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s4)
        h4_min = self.static_states['H'][3]*0.9  # entropy where h3 = h03, v3 = 0
        h4_max = self.total_states['H'][4]      
        
        self.sol_vaned_diff = brentq(vaned_diffuser_inlet_system, h4_min, h4_max, xtol=1e-6)
        
        self.o4 = pitch3 * np.cos(self.Vel_Tri_S['alpha4']) - self.params['t_b']  # throat opening
        self.A3_th = self.o4 * self.params['b2'] * self.n_blade_R
        
        "S5) Vaned Diffuser Outlet"
        
        # alpha 5 so that p5 = p_ex
        self.params['r5'] = self.params['r3']*self.inputs['r5_r3']
        self.A5 = self.A3_th*self.params['b5_b3']*self.inputs['r5_r3']
        self.params['b5'] = self.params['b3']*self.params['b5_b3']
        
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
                xhi4=self.Vel_Tri_S['alpha4'], xhi5=alpha5)
        
            self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.static_states['S'][4])
            h5is = self.AS.hmass()
        
            h5_new = h5is + self.stator_losses['tot']
            self.AS.update(CP.HmassSmass_INPUTS, h5_new, self.static_states['S'][4])
            p5     = self.AS.p()
            self.update_static_AS(CP.HmassP_INPUTS, h5_new, p5, 5)
            self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)
        
            # Residual 1 — enthalpy consistency (normalized)
            res1 = (h5 - h5_new) / h05
        
            # Residual 2 — exit pressure matches imposed back-pressure (normalized)
            res2 = (p5 - self.inputs['p_ex']) / self.inputs['p_ex']
        
            self.params['xhi5'] = alpha5
        
            return [res1, res2]
        
        from scipy.optimize import least_squares

        # Initial guess
        alpha5_guess = self.Vel_Tri_S['alpha4'] * 0.5  # some residual swirl
        v5_guess     = self.inputs['mdot'] / (self.static_states['D'][4] * self.A5)
        
        x0 = [alpha5_guess, v5_guess]
        
        v5_upper = self.inputs['mdot'] / (self.static_states['D'][4] * self.A5)
        
        alpha5_lb = 0*np.pi/180
        alpha5_ub = 85*np.pi/180
        
        self.sol_sys_stator = least_squares(
            system_stator,
            x0     = [alpha5_guess, v5_guess],
            bounds = ([alpha5_lb, 1.0],          # lower bounds [alpha5, v5]
                      [alpha5_ub, v5_upper]),     # upper bounds [alpha5, v5]
            method = 'trf'                        # trust region reflective, handles bounds
        )
        
        alpha5_sol, v5_sol = self.sol_sys_stator.x
        
        return
    
    def designSystem(self, x):
        
        self.penalty_factor = 100
        
        if "Omega_bounds" in self.params or "omega_bounds" in self.params:   
            self.inputs['psi_is']   = x[0]
            self.inputs['r1s_r2']   = x[1]
            self.inputs['r1h_r1s']  = x[2]
            self.inputs['b2_r2']    = x[3]
            self.inputs['r5_r3']    = x[4]
            self.inputs['r3_r2']    = x[5]
            self.inputs['xhi1']     = x[6]
            self.inputs['xhi2']     = x[7]
            
            # omega_step = self.params.get('omega_step_rpm', 1000)
            # x[8] = np.round(x[8] / omega_step) * omega_step
            
            self.inputs['Omega']  = self.params['Omega']   = x[8]
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
    
        self.PR = self.inputs['p_ex']/self.inputs['p0_su']
        self.n_blade_R = np.floor(12.03 + 2.544*self.PR) # Without splitter : PR between 1 and 5
        # self.n_blade_R = np.floor(-4.527*np.exp(1.865/self.PR) + 32.22) # With splitters : PR between 1.8 and 8
        self.omega = self.params['Omega'] * (2*np.pi) / 60 # rad/s
        
        try:            
            self.designRotor()
            
            eta_diff_min = 0.95
            # eta_diff_min = 0.8
            
            p_req = self.inputs['p_ex'] / eta_diff_min
            
            if self.total_states['P'][2] < p_req:
                res = abs((p_req - self.total_states['P'][2])/p_req)
                
                return res*500

        except:
            return 2222
        
        try:            
            self.designStator()
        except:
            return 1100
        
        # ------------------------------------------------------------------ #
        # Contraintes physiques — Meroni et al. (2018), Additional constraints
        # Retourne une pénalité si une contrainte est violée.
        # ------------------------------------------------------------------ #
    
        p = self.params   # raccourci
    
        # # Pas de prérotation sans IGV
        # if self.Vel_Tri_R['vu1'] > 0.1 * self.Vel_Tri_R['u1']:
        #     return 1000
        
        self.error_log = []
        penalty = 0.0
        
        # C1 — Mach relatif shroud inducteur : M1s_rel ≤ M1s_rel_max
        if self.M1s_rel > p['M1s_rel_max']:
            delta = (self.M1s_rel - p['M1s_rel_max']) / p['M1s_rel_max']
            penalty += delta
            self.error_log.append(('M1s_rel', delta))
        
        # C2 — Mach relatif moyen inducteur : M1_rel ≤ M1_rel_max
        if self.M1_rel > p['M1_rel_max']:
            delta = (self.M1_rel - p['M1_rel_max']) / p['M1_rel_max']
            penalty += delta
            self.error_log.append(('M1_rel', delta))
        
        # C3 — Ratio de diffusion : W2/W1s ≥ W2_W1s_min
        if self.W2_W1s < p['W2_W1s_min']:
            delta = (p['W2_W1s_min'] - self.W2_W1s) / p['W2_W1s_min']
            penalty += delta
            self.error_log.append(('W2_W1s', delta))
        
        # C4 — Angle absolu sortie rotor : alpha2 ≤ alpha2_max
        alpha2_deg = abs(self.Vel_Tri_R['alpha2']) * 180.0 / np.pi
        if alpha2_deg > p['alpha2_max']:
            delta = (alpha2_deg - p['alpha2_max']) / p['alpha2_max']
            penalty += delta
            self.error_log.append(('alpha2_deg', delta))
        
        # C5 — Ouverture gorge inducteur : o1_min ≤ o1 ≤ o1_max
        if self.o1 < p['o1_min']:
            delta = (p['o1_min'] - self.o1) / p['o1_min']
            penalty += delta
            self.error_log.append(('o1_low', delta))
        elif self.o1 > p['o1_max']:
            delta = (self.o1 - p['o1_max']) / p['o1_max']
            penalty += delta
            self.error_log.append(('o1_high', delta))
        
        # C6 — Degré de réaction : DR_min ≤ DR ≤ DR_max
        if self.DR < p['DR_min']:
            delta = (p['DR_min'] - self.DR) / (p['DR_max'] - p['DR_min'])
            penalty += delta
            self.error_log.append(('DR_low', delta))
        elif self.DR > p['DR_max']:
            delta = (self.DR - p['DR_max']) / (p['DR_max'] - p['DR_min'])
            penalty += delta
            self.error_log.append(('DR_high', delta))
        
        # C7 — Tip speed : U2 ≤ U2_max
        if self.Vel_Tri_R['u2'] > p['U2_max']:
            delta = (self.Vel_Tri_R['u2'] - p['U2_max']) / p['U2_max']
            penalty += delta
            self.error_log.append(('u2', delta))
        
        # C8 — Ratio diffuseur vaneless : r3_r2_min ≤ r3/r2 ≤ r3_r2_max
        if self.inputs['r3_r2'] < p['r3_r2_min']:
            delta = (p['r3_r2_min'] - self.inputs['r3_r2']) / p['r3_r2_min']
            penalty += delta
            self.error_log.append(('r3_r2_low', delta))
        elif self.inputs['r3_r2'] > p['r3_r2_max']:
            delta = (self.inputs['r3_r2'] - p['r3_r2_max']) / p['r3_r2_max']
            penalty += delta
            self.error_log.append(('r3_r2_high', delta))
        
        # ------------------------------------------------------------------ #
        # Efficacité totale-statique
        # ------------------------------------------------------------------ #
        self.AS.update(CP.PSmass_INPUTS, self.static_states['P'][5], self.total_states['S'][1])
        hout_is = self.AS.hmass()
        self.eta_is = (hout_is - self.total_states['H'][1]) / \
                      (self.static_states['H'][5] - self.total_states['H'][1])
        
        # ------------------------------------------------------------------ #
        # Efficacité totale-totale
        # ------------------------------------------------------------------ #
        self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][5], self.total_states['S'][1])
        hout_is = self.AS.hmass()
        self.eta_is_tt = (hout_is - self.total_states['H'][1]) / \
                         (self.total_states['H'][5] - self.total_states['H'][1])
        
        if self.res_rotor_ex > 1e-3:
            penalty += self.res_rotor_ex
            
        if penalty > 0:
            return -self.eta_is + penalty*self.penalty_factor
        
        self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][5], self.total_states['S'][1])
        h05_is = self.AS.hmass()
        self.Dh0s_2 = h05_is - self.total_states['H'][1]
               
        return -self.eta_is
    
    def design(self):
        
        if "omega_bounds" in self.params:
            bounds = (np.array([
                self.params['psi_is_bounds'][0],
                self.params['r1s_r2_bounds'][0],
                self.params['r1h_r1s_bounds'][0],
                self.params['b2_r2_bounds'][0],
                self.params['r5_r3_bounds'][0],
                self.params['r3_r2_bounds'][0],
                self.params['xhi1_bounds'][0],
                self.params['xhi2_bounds'][0],
                self.params['omega_bounds'][0],
            ]),
            np.array([
                self.params['psi_is_bounds'][1],
                self.params['r1s_r2_bounds'][1],
                self.params['r1h_r1s_bounds'][1],
                self.params['b2_r2_bounds'][1],
                self.params['r5_r3_bounds'][1],
                self.params['r3_r2_bounds'][1],
                self.params['xhi1_bounds'][1],
                self.params['xhi2_bounds'][1],
                self.params['omega_bounds'][1],
            ]))
            
        elif "Omega_bounds" in self.params:
            bounds = (np.array([
                self.params['psi_is_bounds'][0],
                self.params['r1s_r2_bounds'][0],
                self.params['r1h_r1s_bounds'][0],
                self.params['b2_r2_bounds'][0],
                self.params['r5_r3_bounds'][0],
                self.params['r3_r2_bounds'][0],
                self.params['xhi1_bounds'][0],
                self.params['xhi2_bounds'][0],
                self.params['Omega_bounds'][0],
            ]),
            np.array([
                self.params['psi_is_bounds'][1],
                self.params['r1s_r2_bounds'][1],
                self.params['r1h_r1s_bounds'][1],
                self.params['b2_r2_bounds'][1],
                self.params['r5_r3_bounds'][1],
                self.params['r3_r2_bounds'][1],
                self.params['xhi1_bounds'][1],
                self.params['xhi2_bounds'][1],
                self.params['Omega_bounds'][1],
            ]))
            
        else:
            bounds = (np.array([
                self.params['psi_is_bounds'][0],
                self.params['r1s_r2_bounds'][0],
                self.params['r1h_r1s_bounds'][0],
                self.params['b2_r2_bounds'][0],
                self.params['r5_r3_bounds'][0],
                self.params['r3_r2_bounds'][0],
                self.params['xhi1_bounds'][0],
                self.params['xhi2_bounds'][0],
            ]),
            np.array([
                self.params['psi_is_bounds'][1],
                self.params['r1s_r2_bounds'][1],
                self.params['r1h_r1s_bounds'][1],
                self.params['b2_r2_bounds'][1],
                self.params['r5_r3_bounds'][1],
                self.params['r3_r2_bounds'][1],
                self.params['xhi1_bounds'][1],
                self.params['xhi2_bounds'][1],
            ]))
    
        def objective_wrapper(x):
            costs, wdots = [], []
            for xi in x:
                c = self.designSystem(xi)
                costs.append(c)
            return np.asarray(costs, dtype=float)
    
        optimizer = ps.single.GlobalBestPSO(
            n_particles=100,
            dimensions=len(bounds[0]),
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        patience = 15
        tol = 1e-3
        max_iter = 100
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
        self.designSystem(best_pos)
        
        return

if __name__ == "__main__":

    fluid = "Air_1"    
    
    eta_is_vec = []
    
    for i in range(1):
        
        if fluid == "CO2":
            Comp = RadialCPMLDesign('CO2')
            
            # CO2 Case
            # Comp.set_inputs(
            #     mdot = 0.8, # kg/s
            #     p0_su = 76.9*1e5, # Pa
            #     T0_su = 305.97, # K
            #     p_ex = 127*1e5, # Pa
            #     )
            
            # Comp.set_inputs(
            #     mdot = 2.15, # kg/s
            #     p0_su = 76.9*1e5, # Pa
            #     T0_su = 305.97, # K
            #     p_ex = 127*1e5, # Pa
            #     )
                    
            Comp.set_inputs(
                mdot = 2.15, # kg/s
                p0_su = 76.9*1e5, # Pa
                T0_su = 305.97, # K
                p_ex = 96.894*1e5, # Pa
                )
            
            Comp.set_parameters(
                # Fixed Values from source
                t_b = 0.762*1e-3, # [m] : blade thickness
                
                eps_imp = 0.254*1e-3, # [m] : Impeller clearance
                eps_bf_imp = 0.254*1e-3, # [m] : Impeller back face clearance
                k_imp = 0.01*1e-3, # [m] : Impeller surface roughness
                
                b3_b2_ratio = 1, # [-] : Diffuser width ratio
                b5_b3 = 1, 
                
                L_z = 0.1137,
                
                CP = 0.44, # [-] : Diffuser pressure coefficient 0.44
                Omega = 50000, 
                
                psi_is_bounds = [0.3, 1.1],
                r1s_r2_bounds = [0.4, 0.7], 
                r1h_r1s_bounds = [0.25, 0.4],
                b2_r2_bounds = [0.02, 0.1], 
                r5_r3_bounds = [1.01, 1.5],
                r3_r2_bounds = [1.05, 2],
                xhi1_bounds = [40, 70],
                xhi2_bounds = [20, 55],
                # omega_bounds = np.array([30, 80])*1e3,
                
                M1s_rel_max  = 1.4,    # Mach relatif shroud inducteur  ≤ 1.4
                M1_rel_max   = 0.9,    # Mach relatif moyen inducteur   ≤ 0.9
                W2_W1s_min   = 0.25,   # Ratio diffusion W2/W1s         ≥ 0.25
                alpha2_max   = 85.0,   # Angle absolu sortie rotor [°]  ≤ 85°
                o1_min       = 1.49e-3,# Ouverture gorge inducteur [m]  ≥ 1.49 mm
                o1_max       = 50e-3,  # Ouverture gorge inducteur [m]  ≤ 50 mm
                DR_min       = -0.1,   # Degré de réaction              ≥ -0.1
                DR_max       = 0.9,    # Degré de réaction              ≤ 0.9
                U2_max       = 400.0,  # Tip speed [m/s]                ≤ 400
                r3_r2_min    = 1.05,   # Vaneless diffuser radius ratio ≥ 1.05
                r3_r2_max    = 2.0,    # Vaneless diffuser radius ratio ≤ 2.0  
                )
            
            Comp.design()
        
        elif fluid == "R134a":
            Comp = RadialCPMLDesign('R134a')
            
            # Zorlu Case
            
            Comp.set_inputs(
                mdot = 0.039, # kg/s
                p0_su = 1.65*1e5, # Pa
                T0_su = 265, # K
                p_ex = 3.8775*1e5, # Pa
                )
            
            Comp.set_parameters(
                # Fixed Values from source
                t_b = 0.1*1e-3, # [m] : blade thickness
                
                eps_imp = 0.15*1e-3, # [m] : Impeller clearance
                eps_bf_imp = 1*1e-3, # [m] : Impeller back face clearance
                k_imp = 0.01*1e-3, # [m] : Impeller surface roughness
                
                b3_b2_ratio = 1, # [-] : Diffuser width ratio
                b5_b3 = 1, 
            
                L_z = 0.007693,
                CP = 0.44, # [-] : Diffuser pressure coefficient 0.44
                Omega = 180000, 
                
                psi_is_bounds = [0.3, 0.5],
                r1s_r2_bounds = [0.6, 0.7], 
                r1h_r1s_bounds = [0.3, 0.4],
                b2_r2_bounds = [0.02, 0.3], 
                r5_r3_bounds = [1.01, 1.3],
                r3_r2_bounds = [1.05, 1.2],
                xhi1_bounds = [40, 60],
                xhi2_bounds = [40, 55],
                # omega_bounds = np.array([150, 210])*1e3,
    
                M1s_rel_max  = 1.4,    # Mach relatif shroud inducteur  ≤ 1.4
                M1_rel_max   = 0.9,    # Mach relatif moyen inducteur   ≤ 0.9
                W2_W1s_min   = 0.25,   # Ratio diffusion W2/W1s         ≥ 0.25
                alpha2_max   = 85.0,   # Angle absolu sortie rotor [°]  ≤ 85°
                o1_min       = 1e-3,# Ouverture gorge inducteur [m]  ≥ 1.49 mm
                o1_max       = 50e-3,  # Ouverture gorge inducteur [m]  ≤ 50 mm
                DR_min       = -0.1,   # Degré de réaction              ≥ -0.1
                DR_max       = 0.9,    # Degré de réaction              ≤ 0.9
                U2_max       = 400.0,  # Tip speed [m/s]                ≤ 400
                r3_r2_min    = 1.05,   # Vaneless diffuser radius ratio ≥ 1.05
                r3_r2_max    = 2.0,    # Vaneless diffuser radius ratio ≤ 2.0  
                )
            
            Comp.design()
        
        elif fluid == "Air_1":
            Comp = RadialCPMLDesign('Air')
            
            # Zorlu Case
            
            Comp.set_inputs(
                mdot = 5.32, # kg/s
                p0_su = 1.01*1e5, # Pa
                T0_su = 288.15, # K
                p_ex = 2.02*1e5, # Pa
                )
            
            Comp.set_parameters(
                # Fixed Values from source
                t_b = 2.11*1e-3, # [m] : blade thickness
                
                eps_imp = 0.372*1e-3, # [m] : Impeller clearance
                eps_bf_imp = 0.372*1e-3, # [m] : Impeller back face clearance
                k_imp = 0.002*1e-3, # [m] : Impeller surface roughness
                
                b3_b2_ratio = 1, # [-] : Diffuser width ratio
                b5_b3 = 1, 
            
                L_z = 0.13,
                CP = 0.44, # [-] : Diffuser pressure coefficient 0.44
                Omega = 14000, 
                
                # psi_is_bounds = [0.3, 1.1],
                psi_is_bounds = [0.3, 0.7],
                r1s_r2_bounds = [0.4, 0.7], 
                r1h_r1s_bounds = [0.25, 0.4],
                b2_r2_bounds = [0.02, 0.8], 
                r5_r3_bounds = [1.01, 1.5],
                r3_r2_bounds = [1.05, 2],
                xhi1_bounds = [40, 70],
                xhi2_bounds = [20, 55],
                # omega_bounds = np.array([10, 20])*1e3,
    
                M1s_rel_max  = 1.4,    # Mach relatif shroud inducteur  ≤ 1.4
                M1_rel_max   = 0.9,    # Mach relatif moyen inducteur   ≤ 0.9
                W2_W1s_min   = 0.25,   # Ratio diffusion W2/W1s         ≥ 0.25
                alpha2_max   = 85.0,   # Angle absolu sortie rotor [°]  ≤ 85°
                o1_min       = 1e-3,# Ouverture gorge inducteur [m]  ≥ 1.49 mm
                o1_max       = 50e-3,  # Ouverture gorge inducteur [m]  ≤ 50 mm
                DR_min       = -0.1,   # Degré de réaction              ≥ -0.1
                DR_max       = 0.9,    # Degré de réaction              ≤ 0.9
                U2_max       = 400.0,  # Tip speed [m/s]                ≤ 400
                r3_r2_min    = 1.05,   # Vaneless diffuser radius ratio ≥ 1.05
                r3_r2_max    = 2.0,    # Vaneless diffuser radius ratio ≤ 2.0  
                )
            
            Comp.design()
        
        elif fluid == "Air_2":
            Comp = RadialCPMLDesign('Air')
            
            # Zorlu Case
            
            Comp.set_inputs(
                mdot = 4.54, # kg/s
                p0_su = 1.01*1e5, # Pa
                T0_su = 288.15, # K
                p_ex = 1.8281*1e5, # Pa
                )
            
            Comp.set_parameters(
                # Fixed Values from source
                t_b = 2.11*1e-3, # [m] : blade thickness
                
                eps_imp = 0.235*1e-3, # [m] : Impeller clearance
                eps_bf_imp = 0.235*1e-3, # [m] : Impeller back face clearance
                k_imp = 0.002*1e-3, # [m] : Impeller surface roughness
                
                b3_b2_ratio = 1, # [-] : Diffuser width ratio
                b5_b3 = 1, 
            
                L_z = 0.13,
                CP = 0.44, # [-] : Diffuser pressure coefficient 0.44
                Omega = 14000, 
                
                psi_is_bounds = [0.3, 1.1],
                r1s_r2_bounds = [0.4, 0.7], 
                r1h_r1s_bounds = [0.25, 0.4],
                b2_r2_bounds = [0.02, 0.8], 
                r5_r3_bounds = [1.01, 1.5],
                r3_r2_bounds = [1.05, 2],
                xhi1_bounds = [40, 70],
                xhi2_bounds = [20, 55],
                # omega_bounds = np.array([10, 20])*1e3,
    
                M1s_rel_max  = 1.4,    # Mach relatif shroud inducteur  ≤ 1.4
                M1_rel_max   = 0.9,    # Mach relatif moyen inducteur   ≤ 0.9
                W2_W1s_min   = 0.25,   # Ratio diffusion W2/W1s         ≥ 0.25
                alpha2_max   = 85.0,   # Angle absolu sortie rotor [°]  ≤ 85°
                o1_min       = 1e-3,# Ouverture gorge inducteur [m]  ≥ 1.49 mm
                o1_max       = 50e-3,  # Ouverture gorge inducteur [m]  ≤ 50 mm
                DR_min       = -0.1,   # Degré de réaction              ≥ -0.1
                DR_max       = 0.9,    # Degré de réaction              ≤ 0.9
                U2_max       = 400.0,  # Tip speed [m/s]                ≤ 400
                r3_r2_min    = 1.05,   # Vaneless diffuser radius ratio ≥ 1.05
                r3_r2_max    = 2.0,    # Vaneless diffuser radius ratio ≤ 2.0  
                )
            
            Comp.design()
            
        elif fluid == "Air_3":
            Comp = RadialCPMLDesign('Air')
            
            # Zorlu Case
            
            Comp.set_inputs(
                mdot = 4.54, # kg/s
                p0_su = 1.01*1e5, # Pa
                T0_su = 288.15, # K
                p_ex = 1.6867*1e5, # Pa
                )
            
            Comp.set_parameters(
                # Fixed Values from source
                t_b = 2.11*1e-3, # [m] : blade thickness
                
                eps_imp = 0.372*1e-3, # [m] : Impeller clearance
                eps_bf_imp = 0.372*1e-3, # [m] : Impeller back face clearance
                k_imp = 0.002*1e-3, # [m] : Impeller surface roughness
                
                b3_b2_ratio = 1, # [-] : Diffuser width ratio
                b5_b3 = 1, 
            
                L_z = 0.13,
                CP = 0.44, # [-] : Diffuser pressure coefficient 0.44
                Omega = 14000, 
                
                psi_is_bounds = [0.3, 1.1],
                r1s_r2_bounds = [0.4, 0.7], 
                r1h_r1s_bounds = [0.25, 0.4],
                b2_r2_bounds = [0.02, 0.8], 
                r5_r3_bounds = [1.01, 1.5],
                r3_r2_bounds = [1.05, 2],
                xhi1_bounds = [40, 70],
                xhi2_bounds = [20, 55],
                # omega_bounds = np.array([10, 20])*1e3,
    
                M1s_rel_max  = 1.4,    # Mach relatif shroud inducteur  ≤ 1.4
                M1_rel_max   = 0.9,    # Mach relatif moyen inducteur   ≤ 0.9
                W2_W1s_min   = 0.25,   # Ratio diffusion W2/W1s         ≥ 0.25
                alpha2_max   = 85.0,   # Angle absolu sortie rotor [°]  ≤ 85°
                o1_min       = 1e-3,# Ouverture gorge inducteur [m]  ≥ 1.49 mm
                o1_max       = 50e-3,  # Ouverture gorge inducteur [m]  ≤ 50 mm
                DR_min       = -0.1,   # Degré de réaction              ≥ -0.1
                DR_max       = 0.9,    # Degré de réaction              ≤ 0.9
                U2_max       = 400.0,  # Tip speed [m/s]                ≤ 400
                r3_r2_min    = 1.05,   # Vaneless diffuser radius ratio ≥ 1.05
                r3_r2_max    = 2.0,    # Vaneless diffuser radius ratio ≤ 2.0  
                )
            
            Comp.design()
        
        try:
            eta_is_vec.append(Comp.eta_is)    
        except:
            eta_is_vec.append(-1)    

