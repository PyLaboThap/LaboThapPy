#!/usr/bin/python3

# --- loading libraries 

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution, least_squares, root
import pyswarms as ps

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
        self.Vel_Tri_Last_Stage = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None

    def reset(self):

        self.AS = CP.AbstractState('HEOS', self.fluid)
        
        # Blade Dictionnary
        self.stages = []
        
        # Velocity Triangle Data
        self.Vel_Tri = {}
        self.Vel_Tri_Last_Stage = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None
        
        return
    # ---------------- Stage Sub Class ----------------------------------------------------------------------
    
    class stage(object):
        
        def __init__(self, fluid, n_disc):
            # Define your column and index layout
            cols = ['H', 'S', 'P', 'D', 'A', 'V']
            index = [1, 2, 3]
            
            # Initialize all DataFrames with dtype=object to store arrays
            self.total_states = pd.DataFrame(columns=cols, index=index, dtype=object)
            self.static_states = pd.DataFrame(columns=cols, index=index, dtype=object)
            self.total_states_disc = pd.DataFrame(np.empty((len(index), len(cols)), dtype=object), index=index, columns=cols)
            self.static_states_disc = pd.DataFrame(np.empty((len(index), len(cols)), dtype=object), index=index, columns=cols)
            
            # Fill each cell with the correct array
            for row in self.total_states_disc.index:
                for col in self.total_states_disc.columns:
                    if n_disc == 1:
                        self.total_states_disc.at[row, col] = [0] 
                        self.static_states_disc.at[row, col] = [0] 
                    else:
                        self.total_states_disc.at[row, col] = np.zeros(n_disc) 
                        self.static_states_disc.at[row, col] = np.zeros(n_disc)  
            
            self.AS = CP.AbstractState('HEOS', fluid)
            
            self.eta_is_R = None
            self.eta_is_S = None
            
            self.A_flow_S = None
            self.A_flow_R = None
            
            self.h_blade_S = None
            self.h_blade_R = None
            
            self.chord_S = None
            self.chord_R = None

            self.disc_chord_S = np.zeros([n_disc])
            self.disc_chord_R = np.zeros([n_disc])
            
            self.AR = None
            
            self.disc_Vel_Tri_S = {}
            self.disc_Vel_Tri_R = {}
            
            Vel_Tri_args = ['u', 'psi', 'phi', 'R', 'r', 'vu1', 'v1', 'wu1', 'w1', 'vu2', 'vm', 'v2', 'wu2', 'w2', 'alpha1', 'alpha2', 'beta1', 'beta2']
            
            for arg in Vel_Tri_args:
                self.disc_Vel_Tri_S[arg] = np.zeros(n_disc)
                self.disc_Vel_Tri_R[arg] = np.zeros(n_disc)

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

        def update_total_AS_disc(self, CP_INPUTS, input_1, input_2, position, i):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.total_states_disc['H'][position][i] = self.AS.hmass()            
            self.total_states_disc['S'][position][i] = self.AS.smass()            
            self.total_states_disc['P'][position][i] = self.AS.p()            
            self.total_states_disc['D'][position][i] = self.AS.rhomass()       
            
            try:
                self.total_states_disc['A'][position][i] = self.AS.speed_sound()            
            except:
                self.total_states_disc['A'][position][i] = -1            
                
            self.total_states_disc['V'][position][i] = self.AS.viscosity()            
            
            return
        
        def update_static_AS_disc(self, CP_INPUTS, input_1, input_2, position, i):
            self.AS.update(CP_INPUTS, input_1, input_2)
                        
            self.static_states_disc['H'][position][i] = self.AS.hmass()            
            self.static_states_disc['S'][position][i] = self.AS.smass()            
            self.static_states_disc['P'][position][i] = self.AS.p()            
            self.static_states_disc['D'][position][i] = self.AS.rhomass()            

            try:
                self.static_states_disc['A'][position][i] = self.AS.speed_sound()            
            except:
                self.static_states_disc['A'][position][i] = -1
    
            self.static_states_disc['V'][position][i] = self.AS.viscosity()            

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
        n_blade_plot = np.array(self.n_blade)

        x = np.linspace(0,len(n_blade_plot)-1, len(n_blade_plot))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1     

        plt.figure()
        plt.plot(x[::2], n_blade_plot[::2], 'o', label="Stator Blades")  # even indices
        plt.plot(x[1::2], n_blade_plot[1::2], 'o', label="Rotor Blades")  # odd indices
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

    # ---------------- Speed Profile ------------------------------------------------------------------------

    def compute_speed_profile(self, stage, blade_row_type):
        
        a = self.Vel_Tri['u']*(1-self.R)
        
        if blade_row_type == 'R':
            
            r_h = self.r_m - stage.h_blade_R/2
            
            for i in range(self.params['n_disc']):
                r_h = self.r_m - stage.h_blade_R/2
                stage.disc_Vel_Tri_R['r'][i] = r = r_h + (i+1/2)*stage.h_blade_R/self.params['n_disc']
                
                stage.disc_Vel_Tri_R['u'][i] = r*self.Vel_Tri['u']/self.r_m
                
                stage.disc_Vel_Tri_R['psi'][i] = self.psi * self.Vel_Tri['u']**2 / stage.disc_Vel_Tri_R['u'][i]**2
                stage.disc_Vel_Tri_R['R'][i] = 1-a/stage.disc_Vel_Tri_R['u'][i]
                
                if i == 0:
                    fact_1 = 2*(1-stage.disc_Vel_Tri_R['R'][0])**2 * (np.log(self.r_m/r_h) + self.psi/(2*(1-stage.disc_Vel_Tri_R['R'][0]))*(self.r_m/r_h)*(self.r_m/r_h - 1))
                    stage.disc_Vel_Tri_R['phi'][i] = np.sqrt((self.phi*self.r_m/r_h)**2 + fact_1)
                else:
                    fact_1 = 2*(1-stage.disc_Vel_Tri_R['R'][0])**2 * (np.log(r/r_h) + stage.disc_Vel_Tri_R['psi'][i]/(2*(1-stage.disc_Vel_Tri_R['R'][0]))*(r/r_h)*(r/r_h - 1))
                    stage.disc_Vel_Tri_R['phi'][i] = (r_h/r)*np.sqrt(stage.disc_Vel_Tri_R['phi'][0]**2 - fact_1)
        else:
            
            for i in range(self.params['n_disc']):
                r_h = self.r_m - stage.h_blade_S/2
                stage.disc_Vel_Tri_S['r'][i] = r = r_h + (i+1/2)*stage.h_blade_S/self.params['n_disc']
                
                stage.disc_Vel_Tri_S['u'][i] = r*self.Vel_Tri['u']/self.r_m
                
                stage.disc_Vel_Tri_S['psi'][i] = self.Vel_Tri['u']**2*self.psi/stage.disc_Vel_Tri_S['u'][i]**2
                stage.disc_Vel_Tri_S['R'][i] = 1-a/stage.disc_Vel_Tri_S['u'][i]

                if i == 0:
                    fact_1 = 2*(1-stage.disc_Vel_Tri_S['R'][0])**2 * (np.log(self.r_m/r_h) + self.psi/(2*(1-stage.disc_Vel_Tri_S['R'][0]))*(self.r_m/r_h)*(self.r_m/r_h - 1))
                    stage.disc_Vel_Tri_S['phi'][i] = np.sqrt((self.phi*self.r_m/r_h)**2 + fact_1)
                else:
                    fact_1 = 2*(1-stage.disc_Vel_Tri_S['R'][0])**2 * (np.log(r/r_h) + stage.disc_Vel_Tri_S['psi'][i]/(2*(1-stage.disc_Vel_Tri_S['R'][0]))*(r/r_h)*(r/r_h - 1))
                    stage.disc_Vel_Tri_S['phi'][i] = (r_h/r)*np.sqrt(stage.disc_Vel_Tri_S['phi'][0]**2 - fact_1)
                    
        return

    def computeStageVelTriangle(self, stage_index, blade_row_type):
        
        if blade_row_type == 'R':
            
            for i in range(self.params['n_disc']):
            
                Vel_Tri = self.stages[stage_index].disc_Vel_Tri_R
                
                u = Vel_Tri['u'][i]
                R = Vel_Tri['R'][i]
                psi = Vel_Tri['psi'][i]
                phi = Vel_Tri['phi'][i]
        
                # Inlet Absolute velocities 
                Vel_Tri['vu1'][i] = self.stages[stage_index].disc_Vel_Tri_S['vu2'][i]
                Vel_Tri['vm'][i]  = self.stages[stage_index].disc_Vel_Tri_S['vm'][i]
                Vel_Tri['v1'][i]  = self.stages[stage_index].disc_Vel_Tri_S['v2'][i]

                # Inlet Relative velocities 
                Vel_Tri['wu1'][i] = self.stages[stage_index].disc_Vel_Tri_S['wu2'][i]
                Vel_Tri['w1'][i]  = self.stages[stage_index].disc_Vel_Tri_S['w2'][i]
        
                # Outlet Absolute velocities 
                Vel_Tri['vu2'][i] = vu2 = u*(2*(1-R) - psi)/2
                Vel_Tri['vm'][i]  = vm  = u*phi
                Vel_Tri['v2'][i]  = np.sqrt(vu2**2+vm**2)
    
                # Outlet Relative velocities 
                Vel_Tri['wu2'][i] = wu2 = vu2 - u
                Vel_Tri['w2'][i]  = np.sqrt(wu2**2+vm**2)
    
                # Angles in radians
                Vel_Tri['alpha1'][i] = self.stages[stage_index].disc_Vel_Tri_S['alpha2'][i]
                Vel_Tri['alpha2'][i] = np.arctan(vu2/vm)
        
                Vel_Tri['beta1'][i] = self.stages[stage_index].disc_Vel_Tri_S['beta2'][i]
                Vel_Tri['beta2'][i] = np.arctan(wu2/vm)
        
        else:
            
            for i in range(self.params['n_disc']):
          
                Vel_Tri = self.stages[stage_index].disc_Vel_Tri_S
                
                u = Vel_Tri['u'][i]
                R = Vel_Tri['R'][i]
                psi = Vel_Tri['psi'][i]
                phi = Vel_Tri['phi'][i]
        
                if stage_index == 0: # First Stage

                    # Inlet Absolute velocities 
                    Vel_Tri['vu1'][i] = self.Vel_Tri['vu1']
                    Vel_Tri['vm'][i]  = vm = self.Vel_Tri['vm']
                    Vel_Tri['v1'][i]  = self.Vel_Tri['v1']
    
                    # Inlet Relative velocities 
                    Vel_Tri['wu1'][i] = self.Vel_Tri['wu1']
                    Vel_Tri['w1'][i]  = self.Vel_Tri['w1']
                
                else:
                    
                    # Inlet Absolute velocities 
                    Vel_Tri['vu1'][i] = self.stages[stage_index-1].disc_Vel_Tri_R['vu2'][i]
                    Vel_Tri['vm'][i]  = vm = self.stages[stage_index-1].disc_Vel_Tri_R['vm'][i]
                    Vel_Tri['v1'][i]  = self.stages[stage_index-1].disc_Vel_Tri_R['v2'][i]
    
                    # Inlet Relative velocities 
                    Vel_Tri['wu1'][i] = self.stages[stage_index-1].disc_Vel_Tri_R['wu2'][i]
                    Vel_Tri['w1'][i]  = self.stages[stage_index-1].disc_Vel_Tri_R['w2'][i]
                
                # Outlet Absolute velocities 
                Vel_Tri['vu2'][i] = vu2 = u*(2*(1-R) + psi)/2
                Vel_Tri['v2'][i]  = np.sqrt(vu2**2+vm**2)
    
                # Outlet Relative velocities 
                Vel_Tri['wu2'][i] = wu2 = vu2 - u
                Vel_Tri['w2'][i]  = np.sqrt(wu2**2+vm**2)
    
                # Angles in radians
                if stage_index == 0:   
                    Vel_Tri['alpha1'][i] = self.Vel_Tri['alpha1']
                    Vel_Tri['alpha2'][i] = np.arctan(vu2/vm)
            
                    Vel_Tri['beta1'][i] = self.Vel_Tri['beta1']
                    Vel_Tri['beta2'][i] = np.arctan(wu2/vm)                
                else:
                    Vel_Tri['alpha1'][i] = self.stages[stage_index-1].disc_Vel_Tri_R['alpha2'][i]
                    Vel_Tri['alpha2'][i] = np.arctan(vu2/vm)
            
                    Vel_Tri['beta1'][i] = self.stages[stage_index-1].disc_Vel_Tri_R['beta2'][i]
                    Vel_Tri['beta2'][i] = np.arctan(wu2/vm)
            
        return 

    # ---------------- Loss Models ------------------------------------------------------------------------

    def stator_blade_row_system(self, x):    
        # print(x)
        stage_index = self.curr_stage_index
                
        stage = self.stages[stage_index]
        
        # 1) Guess outlet state
        n = int((len(x) - 2) / 2)
        h_static_out_mean = x[0]*1e5
        p_static_out_mean = x[1]*1e5
        h_static_out = x[2:2 + n]*1e5
        p_static_out = x[2 + n:]*1e5      
                
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu2']**2)
        stage.M_S = max(v,w)/stage.static_states['A'][2]
        
        pout_calc = []
        hout = []
        
        # 2) Solve Mean Line Blade system
        # 2.1) Compute A_flow and h_blade based on r_m guess

        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out_mean, p_static_out_mean, 2)

        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(4*np.pi*self.r_m)

        # 2.2) Compute cord and aspect ratio
            
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])    
        stage.AR_S = stage.h_blade_S/stage.chord_S

        # 2.3) Compute total inlet state
        hin = stage.static_states['H'][1]

        h0in = hin + (self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)/2 

        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)

        # 2.4) Estimate pressure losses 
        # 2.4.1) Balje-Binsley
        
        H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
        t_TE = 5e-4 # m Tailing-edge blade thickness design estimate  - # !!! Limite de fabrication, voir 
        theta = 0.036*stage.chord_S/self.params['Re_min']**0.2 # Boundary layer momentum thickness : c
        t_blade = 0.12*stage.chord_S # Blade thickness estimate : Assumption for NACA 0012 airfoil
        lambda_2_rad = (self.Vel_Tri['beta2']+self.Vel_Tri['beta1'])/2
        
        A = 1-(1+H_TE)*theta-t_TE/t_blade
        B = 1-H_TE*theta-t_TE/t_blade
        
        num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
        den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
        Yp = 1- num_Yp/den_Yp

        # 2.4.2) Secondary loss : Kacker-Okaapu
        Z = self.solidityStator*(self.Vel_Tri['beta1']-self.Vel_Tri['beta2'])/np.cos(self.Vel_Tri['beta2']) # Loading Factor
        Ys = abs(0.0334*1/stage.AR_S*(np.cos(self.Vel_Tri['alpha2'])/np.cos(self.Vel_Tri['beta1']))*Z)

        # 2.5) Pressure loss 
        DP_loss = (Yp+Ys)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)*stage.static_states['D'][2]/2
        p0_out = stage.total_states['P'][1]-DP_loss

        # 2.6) Computation of static outlet state
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]               

        hout_m = h0in-(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2

        stage.update_static_AS(CP.HmassSmass_INPUTS, hout_m, sout, 2)
        
        pout_m = stage.static_states['P'][2]

        # 2.7) Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_m, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        # 3) Solve the ISRE system

        for i in range(self.params['n_disc']):

            stage.update_static_AS_disc(CP.HmassP_INPUTS, h_static_out[i], p_static_out[i], 2, i)
        
            # 3.1) Speed Profile along blade
            self.compute_speed_profile(stage, 'S')
            self.computeStageVelTriangle(stage_index, 'S')
                
            # 3.2) Compute total inlet state
            hin = stage.static_states_disc['H'][1][i]
            h0in = hin + (stage.disc_Vel_Tri_S['vu1'][i]**2 + stage.disc_Vel_Tri_S['vm'][i]**2)/2  
        
            stage.update_total_AS_disc(CP.HmassSmass_INPUTS, h0in, stage.static_states_disc['S'][1][i], 1, i)            
            
            stage.disc_chord_S[i] = (self.params['Re_min']*stage.static_states_disc['V'][2][i])/(stage.static_states_disc['D'][2][i]*stage.disc_Vel_Tri_S['vm'][i])    
            
            # 3.3) Estimate pressure losses 
            # 3.3.1) Balje-Binsley
            H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
            t_TE = 5e-4 # m Tailing-edge blade thickness design estimate  - # !!! Limite de fabrication, voir 
            theta = 0.036*stage.disc_chord_S[i]/self.params['Re_min']**0.2 # Boundary layer momentum thickness : c
            t_blade = 0.12*stage.disc_chord_S[i] # Blade thickness estimate : Assumption for NACA 0012 airfoil
            lambda_2_rad = (stage.disc_Vel_Tri_S['beta2'][i]+stage.disc_Vel_Tri_S['beta1'][i])/2
            
            A = 1-(1+H_TE)*theta-t_TE/t_blade
            B = 1-H_TE*theta-t_TE/t_blade
            
            num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
            den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
            Yp = 1 - num_Yp/den_Yp
    
            # 3.3.2) Secondary loss : Kacker-Okaapu
            Z = self.solidityStator*(stage.disc_Vel_Tri_S['beta1'][i]-stage.disc_Vel_Tri_S['beta2'][i])/np.cos(stage.disc_Vel_Tri_S['beta2'][i]) # Loading Factor
            Ys = abs(0.0334*1/stage.AR_S*(np.cos(stage.disc_Vel_Tri_S['alpha2'][i])/np.cos(stage.disc_Vel_Tri_S['beta1'][i]))*Z)
    
            # 3.4) Pressure loss 
            DP_loss_i = (Yp+Ys)*(stage.disc_Vel_Tri_S['vm'][i]**2 + stage.disc_Vel_Tri_S['vu2'][i]**2)*stage.static_states_disc['D'][2][i]/2
            p0_out_i = stage.total_states_disc['P'][1][i]-DP_loss_i
            
            # 3.5) Computation of static outlet pressure
            # print("----------------------------------------")
            # print(f"h0in : {h0in}")
            # print(f"p0_out_i : {p0_out_i}")
            # print(f"DP_loss_i : {DP_loss_i}")
            
            # print(f"chord : {stage.chord_S}")
            
            # print(f"num_Yp : {num_Yp}")
            # print(f"den_Yp : {den_Yp}")
            # print(f"Yp : {Yp}")
            # print(f"Ys : {Ys}")
            # print(f"rho2 : stage.static_states_disc['D'][2][i]")
            
            stage.update_total_AS_disc(CP.HmassP_INPUTS, h0in, p0_out_i, 2, i)
            sout = stage.total_states_disc['S'][2][i]
            
            # 3.6) Outlet state
            hout_i = h0in-(stage.disc_Vel_Tri_S['vu2'][i]**2 + stage.disc_Vel_Tri_S['vm'][i]**2)/2

            hout.append(hout_i)
            stage.update_static_AS_disc(CP.HmassSmass_INPUTS, hout_i, sout, 2, i)
            
            pout_calc.append(stage.static_states_disc['P'][2][i])
        
        return np.array([hout_m] + [pout_m] + hout + pout_calc)*1e-5 # [(p_static_out - pout_calc)**2 , (h_static_out - hout)**2] # np.concatenate([(p_static_out - pout_calc)/pout_calc, 0.5*(h_static_out - hout)/hout])

    def rotor_blade_row_system(self, x):

        stage_index = self.curr_stage_index

        stage = self.stages[stage_index]

        # print(f"x: {x}")

        # 1) Guess outlet state
        n = int((len(x) - 2) / 2)
        h_static_out_mean = x[0]*1e5
        p_static_out_mean = x[1]*1e5
        h_static_out = x[2:2 + n]*1e5
        p_static_out = x[2 + n:]*1e5  
        
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu3']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)
        stage.M_R = max(v,w)/stage.static_states['A'][3]
        
        pout_calc = []
        hout = []
        
        # 2) Solve Mean Line Blade system
        # 2.1) Compute A_flow and h_blade based on r_m guess
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out_mean, p_static_out_mean, 3)
        
        stage.A_flow_R = self.inputs['mdot']/(stage.static_states['D'][3]*self.Vel_Tri['vm'])
        stage.h_blade_R = stage.A_flow_R/(4*np.pi*self.r_m)
        
        # 2.2) Compute cord and aspect ratio

        stage.chord_R = (self.params['Re_min']*stage.static_states['V'][3])/(stage.static_states['D'][3]*self.Vel_Tri['vm'])            
        stage.AR_R = stage.h_blade_R/stage.chord_R
            
        # 2.3) Compute total inlet state
        hin = stage.static_states['H'][2]
        h0in = hin + (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        # 2.4) Estimate pressure losses 
        # 2.4.1) Balje-Binsley : Profile pressure losses         
        H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
        t_TE = 5e-4 #  Tailing-edge blade thickness design estimate 
        theta = 0.036*stage.chord_R/self.params['Re_min']**0.2 # Boundary layer momentum thickness : Empirical equation for turbulent plate
        t_blade = 0.12*stage.chord_R # Blade thickness estimate : Assumption for NACA 0012 airfoil
        lambda_2_rad = abs((self.Vel_Tri['beta3']+self.Vel_Tri['beta2'])/2)
        
        A = 1-(1+H_TE)*theta-t_TE/t_blade
        B = 1-H_TE*theta-t_TE/t_blade
        
        num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
        den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
        Yp = abs(1- num_Yp/den_Yp)

        # 2.4.2) Kacker-Okaapu : Secondary pressure losses
        Z = self.solidityRotor*(self.Vel_Tri['beta2']-self.Vel_Tri['beta3'])/np.cos(self.Vel_Tri['beta3']) # Loading Factor
        Ys = abs(0.0334*1/stage.AR_R*(np.cos(self.Vel_Tri['alpha3'])/np.cos(self.Vel_Tri['beta2']))*Z)

        # 2.5) Pressure loss 
        DP_loss = (Yp+Ys)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)*stage.static_states['D'][3]/2
        p0_out = stage.total_states['P'][2]-DP_loss
        
        # 2.6) Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.total_states['S'][3]
        
        hout_m = h0in-(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout_m, sout, 3)
        
        pout_m = stage.static_states['P'][3]

        # 2.7) Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_m, stage.static_states['S'][2])
        hout_s = self.AS.hmass()

        stage.eta_is_R = (stage.static_states['H'][2]-stage.static_states['H'][3])/(stage.static_states['H'][2]-hout_s)

        # 3) Solve the ISRE system

        for i in range(self.params['n_disc']):

            stage.update_static_AS_disc(CP.HmassP_INPUTS, h_static_out[i], p_static_out[i], 3, i)
            
            # 3.1) Speed profiles along blade
            
            self.compute_speed_profile(stage, 'R')
            self.computeStageVelTriangle(stage_index, 'R')
                        
            # 3.2) Compute total inlet state
            hin = stage.static_states_disc['H'][2][i]
            h0in = hin + (stage.disc_Vel_Tri_R['wu1'][i]**2 + stage.disc_Vel_Tri_R['vm'][i]**2)/2  
            
            stage.update_total_AS_disc(CP.HmassSmass_INPUTS, h0in, stage.static_states_disc['S'][2][i], 2, i)   
                 
            stage.disc_chord_R[i] = (self.params['Re_min']*stage.static_states_disc['V'][3][i])/(stage.static_states_disc['D'][3][i]*stage.disc_Vel_Tri_R['vm'][i])    
                           
            # 3.3) Estimate pressure losses 
            # 3.3.1) Balje-Binsley : Profile pressure losses         
            H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
            t_TE = 5e-4 #  Tailing-edge blade thickness design estimate 
            theta = 0.036*stage.disc_chord_R[i]/self.params['Re_min']**0.2 # Boundary layer momentum thickness : Empirical equation for turbulent plate
            t_blade = 0.12*stage.disc_chord_R[i] # Blade thickness estimate : Assumption for NACA 0012 airfoil
            lambda_2_rad = (stage.disc_Vel_Tri_R['beta2'][i]+stage.disc_Vel_Tri_R['beta1'][i])/2
            
            A = 1-(1+H_TE)*theta-t_TE/t_blade
            B = 1-H_TE*theta-t_TE/t_blade
            
            num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
            den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
            Yp = abs(1- num_Yp/den_Yp)
    
            # 3.3.2) Kacker-Okaapu : Secondary pressure losses
            Z = self.solidityRotor*(stage.disc_Vel_Tri_R['beta1'][i]-stage.disc_Vel_Tri_R['beta2'][i])/np.cos(stage.disc_Vel_Tri_R['beta2'][i]) # Loading Factor
            Ys = abs(0.0334*1/stage.AR_R*(np.cos(stage.disc_Vel_Tri_R['alpha2'][i])/np.cos(stage.disc_Vel_Tri_R['beta1'][i]))*Z)
    
            # 3.4) Pressure loss 
            DP_loss = (Yp+Ys)*(stage.disc_Vel_Tri_R['vm'][i]**2 + stage.disc_Vel_Tri_R['wu2'][i]**2)*stage.static_states_disc['D'][3][i]/2
            p0_out = stage.total_states_disc['P'][2][i]-DP_loss
            
            # 3.5) Computation of static outlet pressure
            
            stage.update_total_AS_disc(CP.HmassP_INPUTS, h0in, p0_out, 3, i)
            sout = stage.total_states_disc['S'][3][i]
            
            # 3.6) Outlet State
            hout.append(h0in-(stage.disc_Vel_Tri_R['wu2'][i]**2 + stage.disc_Vel_Tri_R['vm'][i]**2)/2)
            stage.update_static_AS_disc(CP.HmassSmass_INPUTS, hout[-1], sout, 3, i)
            
            pout_calc.append(stage.static_states_disc['P'][3][i])
            
        return np.array([hout_m] + [pout_m] + hout + pout_calc)*1e-5 # np.concatenate([(p_static_out - pout_calc)/pout_calc, (h_static_out - hout)/hout])

    def last_blade_row_system(self, x):
        
        stage = self.stages[-1]
        
        # 1) Guess outlet state
        n = int((len(x) - 2) / 2)
        h_static_out_mean = x[0]*1e5
        p_static_out_mean = x[1]*1e5
        h_static_out = x[2:2 + n]*1e5
        p_static_out = x[2 + n:]*1e5      
                
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu2']**2)
        stage.M_S = max(v,w)/stage.static_states['A'][2]
        
        pout_calc = []
        hout = []
        
        # 2) Solve Mean Line Blade system
        # 2.1) Compute A_flow and h_blade based on r_m guess

        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out_mean, p_static_out_mean, 2)

        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(4*np.pi*self.r_m)       
        
        # 2.2) Compute cord and aspect ratio
            
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        # 2.3) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)/2 
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)

        # 2.4) Estimate pressure losses 
        # 2.4.1) Aungier : Profile pressure losses           

        v_1 = np.sqrt(self.Vel_Tri_Last_Stage["vm"]**2 + self.Vel_Tri_Last_Stage["vu1"]**2)
        v_2 = np.sqrt(self.Vel_Tri_Last_Stage["vm"]**2 + self.Vel_Tri_Last_Stage["vu2"]**2)
        
        a = 0.0117 # NACA blade - 0.007 : C.4 circular-arc blade
        
        alpha = 0
        alpha_star = 0
        
        D_e = (np.cos(self.Vel_Tri_Last_Stage['alpha2'])/np.cos(self.Vel_Tri_Last_Stage['alpha1']))*(1.12+a*(alpha - alpha_star)+0.61*np.cos(self.Vel_Tri_Last_Stage['alpha1'])**2 / self.solidityStator * (np.tan(self.Vel_Tri_Last_Stage['alpha1'])-np.tan(self.Vel_Tri_Last_Stage['alpha2'])))
        P_cst = np.cos(self.Vel_Tri_Last_Stage["alpha2"])/2 * self.solidityStator * (v_1/v_2)**2 # Profile Constant
        Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
        # 2.4.2) Cohen : Endwall losses
        EW_Cst = np.cos((self.Vel_Tri_Last_Stage["alpha1"]+self.Vel_Tri_Last_Stage["alpha2"])/2)**3 / np.cos(self.Vel_Tri_Last_Stage["alpha1"])**2  # Endwall Constant
        Yew = 0.02*(self.solidityStator/stage.AR_S)/EW_Cst

        # Pressure loss 
        DP_loss = (Yp+Yew)*(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu1']**2)*stage.static_states['D'][1]/2
        p0_out = stage.total_states['P'][1]-DP_loss

        # 2.6) Computation of static outlet state
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout_m = h0in-(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout_m, sout, 2)
        
        pout_m = stage.static_states['P'][2]

        # 2.7) Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_m, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        # 3) Solve the ISRE system
            
        for i in range(self.params['n_disc']):
    
            stage.update_static_AS_disc(CP.HmassP_INPUTS, h_static_out[i], p_static_out[i], 2, i)
        
            # 3.1) Speed Profile along blade
            self.compute_speed_profile(stage, 'S')
            self.computeStageVelTriangle(-1, 'S')
                
            # 3.2) Compute total inlet state
            hin = stage.static_states_disc['H'][1][i]
            h0in = hin + (stage.disc_Vel_Tri_S['vu1'][i]**2 + stage.disc_Vel_Tri_S['vm'][i]**2)/2  
        
            stage.update_total_AS_disc(CP.HmassSmass_INPUTS, h0in, stage.static_states_disc['S'][1][i], 1, i)            
                        
            # 3.3) Estimate pressure losses 
            # 3.3.1) Aungier : Profile pressure losses           
    
            v_1 = np.sqrt(stage.disc_Vel_Tri_S['vm'][i]**2 + stage.disc_Vel_Tri_S['vu1'][i]**2)
            v_2 = np.sqrt(stage.disc_Vel_Tri_S['vm'][i]**2 + stage.disc_Vel_Tri_S['vu2'][i]**2)
            
            a = 0.0117 # NACA blade - 0.007 : C.4 circular-arc blade
            
            alpha = 0
            alpha_star = 0
            
            D_e = (np.cos(stage.disc_Vel_Tri_S['alpha2'][i])/np.cos(stage.disc_Vel_Tri_S['alpha1'][i]))*(1.12+a*(alpha - alpha_star)+0.61*np.cos(stage.disc_Vel_Tri_S['alpha1'][i])**2 / self.solidityStator * (np.tan(stage.disc_Vel_Tri_S['alpha1'][i])-np.tan(stage.disc_Vel_Tri_S['alpha2'][i])))
            P_cst = np.cos(stage.disc_Vel_Tri_S["alpha2"][i])/2 * self.solidityStator * (v_1/v_2)**2 # Profile Constant
            Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
        
            # 3.3.2) Cohen : Endwall losses
            EW_Cst = np.cos((stage.disc_Vel_Tri_S["alpha1"][i]+stage.disc_Vel_Tri_S["alpha2"][i])/2)**3 / np.cos(stage.disc_Vel_Tri_S["alpha1"][i])**2  # Endwall Constant
            Yew = 0.02*(self.solidityStator/stage.AR_S)/EW_Cst
    
            # Pressure loss 
            DP_loss = (Yp+Yew)*(stage.disc_Vel_Tri_S['vm'][i]**2 + stage.disc_Vel_Tri_S['vu1'][i]**2)*stage.static_states_disc['D'][1][i]/2
            p0_out = stage.total_states_disc['P'][1][i]-DP_loss
            
            # 3.5) Computation of static outlet pressure
            stage.update_total_AS_disc(CP.HmassP_INPUTS, h0in, p0_out, 2, i)
            sout = stage.total_states_disc['S'][2][i]
            
            # 3.6) Outlet state
            hout.append(h0in-(stage.disc_Vel_Tri_S['vu2'][i]**2 + stage.disc_Vel_Tri_S['vm'][i]**2)/2)
            stage.update_static_AS_disc(CP.HmassSmass_INPUTS, hout[-1], sout, 2, i)
            
            pout_calc.append(stage.static_states_disc['P'][2][i])
        
        return np.array([hout_m] + [pout_m] + hout + pout_calc)*1e-5

    # ---------------- Flow Computations ------------------------------------------------------------------

    def computeVelTriangle(self):

        # Velocities over u
        self.Vel_Tri['vu2OverU'] = (2*(1-self.R) + self.psi)/2
        self.Vel_Tri['vu3OverU'] = (2*(1-self.R) - self.psi)/2
        self.Vel_Tri['vmOverU']  = self.phi
        
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
    
    def computeVelTriangleLastStage(self):

        self.Vel_Tri_Last_Stage['u'] = self.Vel_Tri['u']
        self.Vel_Tri_Last_Stage['vu2'] = 0
        self.Vel_Tri_Last_Stage['vu1'] = self.Vel_Tri['vu3']
        self.Vel_Tri_Last_Stage['vm']  = self.Vel_Tri['vm']
        
        self.Vel_Tri_Last_Stage['wu2'] = self.Vel_Tri_Last_Stage['vu2'] - self.Vel_Tri_Last_Stage['u']
        self.Vel_Tri['v2'] = np.sqrt(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w2'] = np.sqrt(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)

        # Angles in radians
        self.Vel_Tri_Last_Stage['alpha1'] = self.Vel_Tri['alpha3'] 
        self.Vel_Tri_Last_Stage['alpha2'] = 0

        self.Vel_Tri_Last_Stage['beta1'] = self.Vel_Tri['beta3']
        self.Vel_Tri_Last_Stage['beta2'] = np.arctan(self.Vel_Tri['u']/self.Vel_Tri['vm'])
        
        a = self.Vel_Tri['u']*(1-self.R)
        
        for i in range(self.params['n_disc']):
      
            Vel_Tri = self.stages[-1].disc_Vel_Tri_S
            
            u = self.stages[-2].disc_Vel_Tri_R['u'][i]
            R = Vel_Tri['R'][i]
            psi = Vel_Tri['psi'][i]
            phi = Vel_Tri['phi'][i]
                
            # Inlet Absolute velocities 
            Vel_Tri['vu1'][i] = self.stages[-2].disc_Vel_Tri_R['vu2'][i]
            Vel_Tri['vm'][i]  = vm = self.stages[-2].disc_Vel_Tri_R['vm'][i]
            Vel_Tri['v1'][i]  = self.stages[-2].disc_Vel_Tri_R['v2'][i]

            # Inlet Relative velocities 
            Vel_Tri['wu1'][i] = self.stages[-2].disc_Vel_Tri_R['wu2'][i]
            Vel_Tri['w1'][i]  = self.stages[-2].disc_Vel_Tri_R['w2'][i]
            
            # Outlet Absolute velocities 
            Vel_Tri['vu2'][i] = vu2 = 0
            Vel_Tri['v2'][i]  = np.sqrt(vu2**2+vm**2)

            # Outlet Relative velocities 
            Vel_Tri['wu2'][i] = wu2 = vu2 - u
            Vel_Tri['w2'][i]  = np.sqrt(wu2**2+vm**2)

            Vel_Tri['alpha1'][i] = self.stages[-2].disc_Vel_Tri_R['alpha2'][i]
            Vel_Tri['alpha2'][i] = np.arctan(vu2/vm)
    
            Vel_Tri['beta1'][i] = self.stages[-2].disc_Vel_Tri_R['beta2'][i]
            Vel_Tri['beta2'][i] = np.arctan(wu2/vm)
    
        return 
    
    def computeBladeRow(self,stage_index,row_type):
        stage = self.stages[stage_index]
        
        self.curr_stage_index = stage_index

        if row_type == 'S': # Stator 
                        
            h_out_guess = stage.static_states['H'][1] - self.Dh0Stage/2
            h_out_guess_disc = np.array(stage.static_states_disc['H'][1]) - self.Dh0Stage/2
            
            self.AS.update(CP.HmassSmass_INPUTS, h_out_guess, stage.static_states['S'][1])
            pout_guess = self.AS.p()
            pout_guess_disc = np.full(len(h_out_guess_disc), pout_guess)
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess], h_out_guess_disc,pout_guess_disc))*1e-5
            
            res = 1
            x_in = x0_disc

            while res > 1e-8:

                x_out = self.stator_blade_row_system(x_in)
                
                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = x_out
            
            self.stator_blade_row_system(x_out)
            
        else: # Rotor 
                    
            h_out_guess = stage.static_states['H'][2] - self.Dh0Stage/2
            h_out_guess_disc = np.array(stage.static_states_disc['H'][2]) - self.Dh0Stage/2
            
            self.AS.update(CP.HmassSmass_INPUTS, h_out_guess, stage.static_states['S'][2])
            pout_guess = self.AS.p()
            pout_guess_disc = np.full(len(h_out_guess_disc), pout_guess)

            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess], h_out_guess_disc,pout_guess_disc))*1e-5

            res = 1
            x_in = x0_disc
            
            while res > 1e-8:

                x_out = self.rotor_blade_row_system(x_in) 

                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = x_out 
            
            self.rotor_blade_row_system(x_out)
            
        return
            
    def computeRepeatingStages(self):
        
        for i in range(int(self.nStages)):
            if i == 0:
                self.computeBladeRow(i, 'S')
                self.computeBladeRow(i, 'R')
            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                self.stages[i].static_states_disc.loc[1] = self.stages[i-1].static_states_disc.loc[3]
                
                self.computeBladeRow(i, 'S')
                self.computeBladeRow(i, 'R')
            
            self.r_tip.append(self.r_m + self.stages[i].h_blade_S/2)
            self.r_hub.append(self.r_m - self.stages[i].h_blade_S/2)
            self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
            self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
        
            self.r_tip.append(self.r_m + self.stages[i].h_blade_R/2)
            self.r_hub.append(self.r_m - self.stages[i].h_blade_R/2)
            self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
            self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
            
        return
    
    def computeLastStage(self):
        self.stages[-1].static_states.loc[1] = self.stages[-2].static_states.loc[3]
        self.stages[-1].static_states_disc.loc[1] = self.stages[-2].static_states_disc.loc[3]
                
        stage = self.stages[-1]
        
        h_out_guess = stage.static_states['H'][1] - self.Dh0Stage/2
        h_out_guess_disc = np.array(stage.static_states_disc['H'][1]) - self.Dh0Stage/2
        
        self.AS.update(CP.HmassSmass_INPUTS, h_out_guess, stage.static_states['S'][1])
        pout_guess = self.AS.p()
        pout_guess_disc = np.full(len(h_out_guess_disc), pout_guess)
        
        # Initial guess vector
        x0_disc = np.concatenate(([h_out_guess], [pout_guess], h_out_guess_disc,pout_guess_disc))*1e-5
            
        res = 1
        x_in = x0_disc

        while res > 1e-6:

            x_out = self.last_blade_row_system(x_in) 
            
            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            
            x_in = x_out
        
        return
    
    def design_system(self, x):
        # try:
            print(x)
            
            self.reset()
            
            self.psi = x[0]
            self.phi = x[1]
            self.R = x[2]
            self.r_m = x[3]
                                    
            self.r_tip = []
            self.r_hub = []
            self.r_hub_tip = []
            self.r_ratio2 = []
            
            # First Stator Instanciation
            self.stages.append(self.stage(self.fluid, self.params['n_disc']))
            self.stages[0].update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
            
            for i in range(self.params['n_disc']):
                self.stages[0].update_total_AS_disc(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1, i)                
            
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
            Dh0Stage = self.psi * self.Vel_Tri['u']**2
            self.nStages = int(round(Dh0/Dh0Stage))
            
            for i in range(self.nStages-1):
                self.stages.append(self.stage(self.fluid, self.params['n_disc']))
            
            # Recompute u based on the number of stages to satisfy the work. As r_m is constant, u is contant accross stages
            self.Dh0Stage = Dh0/self.nStages
            self.Vel_Tri['u'] = np.sqrt(self.Dh0Stage/self.psi)
    
            self.omega_rads = self.Vel_Tri['u']/self.r_m # rad/s
    
            for stage in self.stages:
                stage.disc_Vel_Tri_S['u'] = np.zeros(self.params['n_disc'])
                stage.disc_Vel_Tri_R['u'] = np.zeros(self.params['n_disc'])

                stage.disc_Vel_Tri_S['psi'] = np.zeros(self.params['n_disc'])
                stage.disc_Vel_Tri_R['psi'] = np.zeros(self.params['n_disc'])
                
                stage.disc_Vel_Tri_S['phi'] = np.zeros(self.params['n_disc'])
                stage.disc_Vel_Tri_R['phi'] = np.zeros(self.params['n_disc'])
                
                stage.disc_Vel_Tri_S['R'] = np.zeros(self.params['n_disc'])
                stage.disc_Vel_Tri_R['R'] = np.zeros(self.params['n_disc'])

                stage.disc_Vel_Tri_S['r'] = np.zeros(self.params['n_disc'])
                stage.disc_Vel_Tri_R['r'] = np.zeros(self.params['n_disc'])
        
            "------------- 5) Compute complete velocity triangles and exit losses ----------------------------" 
    
            # Compute velocity triangle with the value of u
            self.Vel_Tri['vm'] = self.Vel_Tri['vmOverU'] * self.Vel_Tri['u']
            self.Vel_Tri['vu2'] = self.Vel_Tri['vu2OverU'] * self.Vel_Tri['u']
            self.Vel_Tri['vu3'] = self.Vel_Tri['vu3OverU'] * self.Vel_Tri['u']
            self.Vel_Tri['wu2'] = self.Vel_Tri['wu2OverU'] * self.Vel_Tri['u']
            self.Vel_Tri['wu3'] = self.Vel_Tri['wu3OverU'] * self.Vel_Tri['u']
            
            self.Vel_Tri['vu1'] = self.Vel_Tri['vu3']
            self.Vel_Tri['v1'] = np.sqrt(self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)
            self.Vel_Tri['wu1'] = self.Vel_Tri['wu3']
            self.Vel_Tri['w1'] = np.sqrt(self.Vel_Tri['wu1']**2 + self.Vel_Tri['vm']**2)
    
            self.exit_loss = (self.Vel_Tri['vm']**2+self.Vel_Tri['vu3']**2)/2
            self.exit_loss_W = self.exit_loss*self.inputs['mdot']
    
            "------------- 6) Iterate on r_m and compute repeating stages ------------------------" 
    
            h_in = self.stages[0].total_states['H'][1] - (self.Vel_Tri['vm']**2)/2
            self.stages[0].update_static_AS(CP.HmassSmass_INPUTS, h_in, s_in, 1)
    
            for i in range(self.params['n_disc']):
                self.stages[0].update_static_AS_disc(CP.HmassSmass_INPUTS, h_in, s_in, 1, i)   
    
            try:
                self.computeRepeatingStages()        
            except:
                print("Error in stages")
                obj = 10000
                return obj
            
            "------------- 7) Compute rotation speed and number of blades per stage ---------------------------" 
    
            self.omega_RPM = self.omega_rads*60/(2*np.pi)
    
            self.n_blade = []
    
            for stage in self.stages:
                  stage.pitch_S = self.solidityStator*stage.chord_S
                  stage.pitch_R = self.solidityRotor*stage.chord_R
    
                  stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
                  self.n_blade.append(stage.n_blade_S)
    
                  stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
                  self.n_blade.append(stage.n_blade_R)
    
            "------------- 8) Add last redirecting stator stage -------------------------------------------------------------" 
        
            self.stages.append(self.stage(self.fluid, self.params['n_disc']))
            
            self.computeVelTriangleLastStage()
            
            try:
                self.computeLastStage()
            except:
                print("Error in last stage")
                obj = 10000
                return obj 

            "------------- 9) Compute main outputs -------------------------------------------------------------" 
            
            hin = self.stages[0].total_states['H'][1]
            hout = self.stages[-1].static_states['H'][2]
            
            self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])
    
            hout_s = self.AS.hmass()
            
            self.W_dot = self.inputs['mdot']*(hin-hout)
            self.eta_is = (hin - hout)/(hin - hout_s)
    
            self.penalty_1 = max(self.r_hub_tip[0] - self.params['r_hub_tip_max'],0)*1000
            self.penalty_2 = max(self.params['r_hub_tip_min'] - self.r_hub_tip[-1],0)*1000
                     
            rel_p_m_dev = abs((self.inputs["p_ex"] - self.stages[-1].static_states['P'][2])/self.inputs["p_ex"])
            rel_p_dev = sum(abs((self.inputs["p_ex"] - np.array(self.stages[-1].static_states_disc['P'][2]))/self.inputs["p_ex"]))
            
            if rel_p_m_dev + rel_p_dev >= 5e-2:
                self.penalty_3 = (rel_p_m_dev + rel_p_dev)*10
                self.Pressure_Deviation_m = self.inputs["p_ex"] - self.stages[-1].static_states['P'][2]
                self.Pressure_Deviations = self.inputs["p_ex"] - np.array(self.stages[-1].static_states_disc['P'][2])
            else:
                self.penalty_3 = 0
                self.Pressure_Deviation_m = self.inputs["p_ex"] - self.stages[-1].static_states['P'][2]
                self.Pressure_Deviations = self.inputs["p_ex"] - np.array(self.stages[-1].static_states_disc['P'][2])

            # print(f"penalty_1, 2, 3 : {penalty_1, penalty_2, penalty_3}")
    
            self.penalty = self.penalty_1 + self.penalty_2 + self.penalty_3
    
            if self.eta_is > 0 and self.eta_is <= 1:
                if self.penalty == 0:
                    if self.W_dot > self.inputs['W_dot']:
                        self.inputs['W_dot'] = self.W_dot
                        
                obj = -self.eta_is + self.penalty
                print(f"opt 'success' : {obj}")
            else:
                print("Bad eta_is")
                obj = 10000
    
            return obj

    def design(self):
        bounds = (np.array([
            self.params['psi_bounds'][0],
            self.params['phi_bounds'][0],
            self.params['R_bounds'][0],
            self.params['r_m_bounds'][0],
        ]),
        np.array([
            self.params['psi_bounds'][1],
            self.params['phi_bounds'][1],
            self.params['R_bounds'][1],
            self.params['r_m_bounds'][1],
        ]))
    
        def objective_wrapper(x):
            # Round inputs if needed for discrete steps
            rounded_x = np.copy(x)
            # (Optional rounding logic goes here)
            return np.array([self.design_system(xi) for xi in rounded_x])
    
        # Initialize the optimizer
        optimizer = ps.single.GlobalBestPSO(
            n_particles=10,
            dimensions=4,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        # Custom stopping logic
        patience = 5
        tol = 1e-3
        max_iter = 100
        no_improve_counter = 0
        best_cost = np.inf
    
        for i in range(max_iter):
            # One iteration step
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
    
            current_best = optimizer.swarm.best_cost
    
            if current_best < best_cost - tol:
                best_cost = current_best
                # best_penalties = [self.penalty_1, self.penalty_2, self.penalty_3]
                no_improve_counter = 0
            else:
                no_improve_counter += 1
    
            print(f"[{i+1}] Best cost: {best_cost:.6f}")
    
            if no_improve_counter >= patience:
                print("Stopping early due to stagnation.")
                break
             
        self.best_pos = optimizer.swarm.best_pos
             
        self.design_system(self.best_pos)
        
        "------------- Print Main Results -------------------------------------------------------------" 
        
        print(f"Parameters : {self.psi, self.phi, self.R}")
        # print(f"Associated penalties : {best_penalties}")
        
        print(f"Turbine mean radius: {self.r_m} [m]")
        print(f"Turbine rotation speed: {self.omega_RPM} [RPM]")
        print(f"Turbine number of stage : {self.nStages} [-]")
        print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")        
        print(f"Turbine Generation : {self.W_dot} [W]")        
        
        return 

case_study = "TCO2_ORC"

if case_study == 'Cuerva':

    Turb = AxialTurbineMeanLineDesign('Cyclopentane')
    
    Turb.set_inputs(
        mdot = 46.18, # kg/s
        W_dot = 4500*1e3, # W
        p0_su = 1230*1e3, # Pa
        T0_su = 273.15 + 158, # K
        p_ex = 78300, # Pa
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
        W_dot = 2506000, # W
        p0_su = 767800, # Pa
        T0_su = 273.15 + 131, # K
        p_ex = 82000, # Pa
        Mmax = 0.8 # [-]
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        Re_min = 1e6, # [-] # 5e5, # [-]
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
        Mmax = 0.8 # [-]
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        Re_min = 6e5, # [-]
        AR_min = 1, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        r_m_bounds = [0.03,0.07], # [m]
        psi_bounds = [0.8,1.4], # [-]
        phi_bounds = [0.4,0.7], # [-] 
        R_bounds = [0.49,0.51], # [-]
        n_disc = 10
        )

Turb.design()

Turb.plot_geometry()
Turb.plot_n_blade()
Turb.plot_radius_verif()
Turb.plot_Mollier()
