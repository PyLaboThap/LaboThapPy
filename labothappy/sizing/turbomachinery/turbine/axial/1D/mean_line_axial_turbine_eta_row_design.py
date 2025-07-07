#!/usr/bin/python3

# --- loading libraries 

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

class AxialTurbineMeanLineDesign(object):

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
        self.Vel_Tri = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None

    # ---------------- Stage Sub Class ----------------------------------------------------------------------
    
    class stage(object):
        
        def __init__(self, fluid):
            self.total_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3])
            self.static_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3])
            self.AS = CP.AbstractState('HEOS', fluid)
            
        def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.total_states['H'][position] = self.AS.hmass()            
            self.total_states['S'][position] = self.AS.smass()            
            self.total_states['P'][position] = self.AS.p()            
            self.total_states['D'][position] = self.AS.rhomass()            
            self.total_states['A'][position] = self.AS.speed_sound()            
            self.total_states['V'][position] = self.AS.viscosity()            
            
            return
        
        def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.static_states['H'][position] = self.AS.hmass()            
            self.static_states['S'][position] = self.AS.smass()            
            self.static_states['P'][position] = self.AS.p()            
            self.static_states['D'][position] = self.AS.rhomass()            
            self.static_states['A'][position] = self.AS.speed_sound()            
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
            hin = stage.static_states['H'][1]
            h0in = hin + (self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)/2  

            stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
            
            hout = h0in - (self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2            
            hout_s = hin - (hin-hout)/self.eta_blade_row
            
            self.AS.update(CP.HmassSmass_INPUTS, hout_s, stage.total_states['S'][1])
            pout = self.AS.p()
            
            stage.update_static_AS(CP.HmassP_INPUTS, hout, pout, 2)            
                        
        else: # Rotor
            hin = stage.static_states['H'][2]
            h0in = hin + (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2  
            
            stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
            
            hout = h0in - (self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)/2            
            hout_s = hin - (hin-hout)/self.eta_blade_row
            
            self.AS.update(CP.HmassSmass_INPUTS, hout_s, stage.total_states['S'][2])
            pout = self.AS.p()
            
            stage.update_static_AS(CP.HmassP_INPUTS, hout, pout, 3)        
        
        return
            
    def computeRepeatingStages(self):
        
        for i in range(int(self.nStages)):
            if i == 0:
                self.computeBladeRow(self.stages[i], 'S')
                self.computeBladeRow(self.stages[i], 'R')
            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
                self.computeBladeRow(self.stages[i], 'S')
                self.computeBladeRow(self.stages[i], 'R')
            
        return
    
    # ---------------- Main Method ------------------------------------------------------------------------
    
    def design(self):
        
        # First Stator Instanciation
        self.stages.append(self.stage(self.fluid))
        self.stages[0].update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        
        "------------- 1) Isentropic Expansion Calculation -----------------------------------------------" 
        s_in = self.stages[0].total_states['S'][1]
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s_in)
        
        h_is_ex = self.AS.hmass()
        Dh0s = self.stages[0].total_states['H'][1] - h_is_ex
        
        Dh0 = self.inputs['W_dot']/self.inputs['mdot']
        
        self.eta_is = Dh0/Dh0s
        
        "------------- 2) Velocity Triangle Computation (+ Solodity) -------------------------------------" 
        self.computeVelTriangle()
        
        self.solidityStator = 2*np.cos(self.Vel_Tri['alpha2'])/np.cos(self.Vel_Tri['alpha1'])*np.sin(abs(self.Vel_Tri['alpha2']-self.Vel_Tri['alpha1']))/self.params['Zweifel']
        self.solidityRotor  = 2*np.cos(self.Vel_Tri['beta3'])/np.cos(self.Vel_Tri['beta2'])*np.sin(abs(self.Vel_Tri['beta3']-self.Vel_Tri['beta2']))/self.params['Zweifel']
        
        "------------- 3) Guess u from vMax (subsonic flow)  ---------------------------------------------" 
        
        vMax = self.AS.speed_sound() * self.inputs['Mmax']
        
        # Assume u based on the maximum speed
        self.Vel_Tri['u'] = vMax / max([self.Vel_Tri['v2OverU'],self.Vel_Tri['w3OverU']])
        
        "------------- 4) Compute number of stage + recompute u  -----------------------------------------" 
        
        # Compute required number of stages based on assumed u
        Dh0Stage = self.inputs['psi'] * self.Vel_Tri['u']**2
        self.nStages = int(round(Dh0/Dh0Stage))
        
        for i in range(self.nStages-1):
            self.stages.append(self.stage(self.fluid))
        
        # Recompute u based on the number of stages to satisfy the work. As r_m is constant, u is contant accross stages
        Dh0Stage = Dh0/self.nStages
        self.Vel_Tri['u'] = np.sqrt(Dh0Stage/self.inputs['psi'])

        "------------- 5) Compute complete velocity triangles and exit losses ----------------------------" 

        # Compute velocity triangle with the value of u
        self.Vel_Tri['vm'] = self.Vel_Tri['vmOverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu2'] = self.Vel_Tri['vu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu3'] = self.Vel_Tri['vu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu2'] = self.Vel_Tri['wu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu3'] = self.Vel_Tri['wu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu1'] = self.Vel_Tri['vu3']

        self.exit_loss = (self.Vel_Tri['vm']**2+self.Vel_Tri['vu3']**2)/2

        "------------- 6) Find eta_blade_row by iterating on the repeating stages ------------------------" 

        h_in = self.stages[0].total_states['H'][1] - (self.Vel_Tri['vm']**2)/2
        self.stages[0].update_static_AS(CP.HmassSmass_INPUTS, h_in, s_in, 1)

        def find_eta_blade(x):
            self.eta_blade_row = x[0]
            self.computeRepeatingStages()
    
            pn_comp = self.stages[-1].static_states['P'][3]


            return (self.inputs["p_ex"] - pn_comp)**2

        sol = minimize(find_eta_blade, 1, bounds=[(self.eta_is, 1)], tol = 1e-4)
        
        "------------- 7) Iterate on r_m to satisfy hub to tip ratio -------------------------------------" 
        self.cord = np.zeros([self.nStages,2])
        self.h_blade = np.zeros([self.nStages,2])
        self.AR = np.zeros([self.nStages,2])

        self.A_flow = np.zeros([self.nStages,2])
        
        self.pitch = np.zeros([self.nStages,2])
        self.n_blade = np.zeros([self.nStages,2])

        def find_r_m(x):
            self.r_m = x[0]
    
            self.r_tip = []
            self.r_hub = []
            self.r_hub_tip = []
            self.r_ratio2 = []
    
            for i in range(self.nStages):
                self.A_flow[i][0] = self.inputs['mdot']/(self.stages[i].static_states['D'][2]*self.Vel_Tri['vm'])
                self.A_flow[i][1] = self.inputs['mdot']/(self.stages[i].static_states['D'][3]*self.Vel_Tri['vm'])
    
                # Determine minimum chord to satisfy minimum Reynolds
                # by using velocity, density and viscosity at the blade outlet
    
                self.h_blade[i][0] = self.A_flow[i][0]/(4*np.pi*self.r_m)
                self.h_blade[i][1] = self.A_flow[i][1]/(4*np.pi*self.r_m)
    
                self.cord[i][0] = (self.params['Re_min']*self.stages[i].static_states['V'][2])/(self.stages[i].static_states['D'][2]*self.Vel_Tri['vm'])
                self.cord[i][1] = (self.params['Re_min']*self.stages[i].static_states['V'][3])/(self.stages[i].static_states['D'][3]*self.Vel_Tri['vm'])
        
                self.AR[i][0] = self.h_blade[i][0]/self.cord[i][0]
                self.AR[i][1] = self.h_blade[i][1]/self.cord[i][1]
            
                self.r_tip.append(self.r_m + self.h_blade[i][0]/2)
                self.r_hub.append(self.r_m - self.h_blade[i][0]/2)
                self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
                self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
            
                self.r_tip.append(self.r_m + self.h_blade[i][1]/2)
                self.r_hub.append(self.r_m - self.h_blade[i][1]/2)
                self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
                self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)

            if self.r_hub_tip[-1] > 0: # Penalty to prevent converging to values not satisfying conditions on r_hub_tip
                penalty_1 = max(self.r_hub_tip[0] - self.params['r_hub_tip_max'],0)*1000
                penalty_2 = max(self.params['r_hub_tip_min'] - self.r_hub_tip[-1],0)*1000
                
                return self.r_m + penalty_1 + penalty_2
            
            else: # A very high penalty prevents converging to r_m values very close to 0,  
                return self.r_m + 100000

        sol = minimize(find_r_m, bounds=[(0, 10)], x0=0.2, tol = 1e-4)        

        "------------- 8) Compute rotation speed and number of blades per stage ---------------------------" 

        self.omega_rads = self.Vel_Tri['u']/self.r_m # rad/s
        self.omega_RPM = self.omega_rads*60/(2*np.pi) 

        for i in range(self.nStages):
              self.pitch[i][0] = self.solidityStator*self.cord[i][0]
              self.pitch[i][1] = self.solidityRotor*self.cord[i][1]

              self.n_blade[i][0] = round(2*np.pi*self.r_m/self.pitch[i][0])
              self.n_blade[i][1] = round(2*np.pi*self.r_m/self.pitch[i][1])

        "------------- 9) Print Main Results -------------------------------------------------------------" 
        
        print(f"Turbine mean radius: {self.r_m} [m]")
        print(f"Turbine rotation speed: {self.omega_RPM} [RPM]")
        print(f"Turbine number of stage : {self.nStages} [-]")
        print(f"Turbine static-to-static blade efficiency : {self.eta_blade_row} [-]")

        return

case_study = "TCO2_ORC"

if case_study == 'Cuerva':

    Turb = AxialTurbineMeanLineDesign('Cyclopentane')
    
    Turb.set_inputs(
        mdot = 46.18, # kg/s
        W_dot_req = 4257000, # W
        p0_su = 1230000, # Pa
        T0_su = 273.15 + 158, # K
        p_ex = 78300, # Pa
        psi = 1, # [-]
        phi = 0.6, # [-]
        R = 0.5, # [-]
        Mmax = 0.8 # [-]
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        Re_min = 5e5, # [-]
        AR_min = 1, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        r_m_bounds = [0.4,0.5], # [m]
        psi_bounds = [0.9,1.4], # [-]
        phi_bounds = [0.4,0.6], # [-] 
        R_bounds = [0.49,0.51], # [-]
        n_disc = 10
        )

elif case_study == 'Zorlu':
    
    Turb = AxialTurbineMeanLineDesign('Cyclopentane')

    Turb.set_inputs(
        mdot = 34.51, # kg/s
        W_dot_req = 2506000, # W
        p0_su = 767800, # Pa
        T0_su = 273.15 + 131, # K
        p_ex = 82000, # Pa
        psi = 1, # [-]
        phi = 0.6, # [-]
        R = 0.5, # [-]
        Mmax = 0.8 # [-]
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        Re_min = 5e5, # [-]
        AR_min = 1, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        r_m_bounds = [0.16,0.4], # [m]
        psi_bounds = [0.8,1.4], # [-]
        phi_bounds = [0.5,0.7], # [-] 
        R_bounds = [0.49,0.51], # [-]
        n_disc = 15
        )
    
elif case_study == 'TCO2_ORC':

    Turb = AxialTurbineMeanLineDesign('CO2')

    Turb.set_inputs(
        mdot = 100, # kg/s
        W_dot = 4.69*1e6, # W
        p0_su = 140*1e5, # Pa
        T0_su = 273.15 + 121, # K
        p_ex = 39.8*1e5, # Pa
        psi = 1, # [-]
        phi = 0.6, # [-]
        R = 0.5, # [-]
        Mmax = 0.8 # [-]
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        Re_min = 5e5, # [-]
        AR_min = 1, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        r_m_bounds = [0.01,1], # [m]
        psi_bounds = [0.8,1.4], # [-]
        phi_bounds = [0.5,0.7], # [-] 
        R_bounds = [0.49,0.51], # [-]
        n_disc = 3
        )


Turb.design()

Turb.plot_geometry()
Turb.plot_n_blade()
Turb.plot_radius_verif()
Turb.plot_Mollier()
