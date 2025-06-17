#!/usr/bin/python3

# --- loading libraries 

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

class AxialCPMLLossCorrDesign(object):

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
            
            self.eta_is_R = None
            self.eta_is_S = None
            
            self.A_flow_S = None
            self.A_flow_R = None
            
            self.h_blade_S = None
            self.h_blade_R = None
            
            self.chord_S = None
            self.chord_R = None
            
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
                self.total_states['A'][position] = self.AS.speed_sound()            
            except:
                self.total_states['A'][position] = -1
    
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

    # ---------------- Loss Models ------------------------------------------------------------------------

    def stator_blade_row_system(self, x, stage):
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(4*np.pi*self.r_m)

        # 4) If first stage : Compute max AR allowed in the turbine to estimate good AR for all stages
        
        if stage == self.stages[0]:
            chord_min = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
            self.AR_max = stage.h_blade_S/chord_min
            
            self.AR = np.linspace(self.params['AR_min'],self.AR_max,self.nStages*2)

            c = 0
            
            for i in range(len(self.AR)):
                if np.mod(i,2): # Rotor
                
                    self.stages[c].AR_R = self.AR[i]
                    c = c+1
                else: # Stator
                    self.stages[c].AR_S = self.AR[i]
        
        # 5) Compute cord based on h_blade and AR
        stage.chord_S = stage.h_blade_S/stage.AR_S
        
        # 6) Estimate pressure losses 
        # 6.1) Balje-Binsley
        H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
        t_TE = stage.chord_S*0.03 # Tailing-edge blade thickness design estimate 
        theta = 0.036*stage.chord_S/self.params['Re_min']**0.2 # Boundary layer momentum thickness : c
        t_blade = 0.12*stage.chord_S # Blade thickness estimate : Assumption for NACA 0012 airfoil
        lambda_2_rad = (self.Vel_Tri['beta2']+self.Vel_Tri['beta1'])/2
        
        A = 1-(1+H_TE)*theta-t_TE/t_blade
        B = 1-H_TE*theta-t_TE/t_blade
        
        num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
        den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
        Yp = 1- num_Yp/den_Yp

        # Secondary loss : Kacker-Okaapu
        Z = self.solidityStator*(self.Vel_Tri['beta1']-self.Vel_Tri['beta2'])/np.cos(self.Vel_Tri['beta2']) # Loading Factor
        Ys = abs(0.0334*1/self.params['AR_min']*(np.cos(self.Vel_Tri['alpha2'])/np.cos(self.Vel_Tri['beta1']))*Z)

        # Pressure loss 
        DP_loss = (Yp+Ys)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)*stage.static_states['D'][2]/2
        p0_out = stage.total_states['P'][1]-DP_loss
        
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
        
        pout_calc = stage.static_states['P'][2]

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        return (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def rotor_blade_row_system(self, x, stage):
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][2]
        h0in = hin + (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_R = self.inputs['mdot']/(stage.static_states['D'][3]*self.Vel_Tri['vm'])
        stage.h_blade_R = stage.A_flow_R/(4*np.pi*self.r_m)
        
        # 4) Compute cord based on h_blade and AR
        stage.chord_R = stage.h_blade_R/stage.AR_R
        
        # 5) Estimate pressure losses 
        # 5.1) Balje-Binsley : Profile pressure losses         
        H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
        t_TE = stage.chord_R*0.03 #  Tailing-edge blade thickness design estimate 
        theta = 0.036*stage.chord_R/self.params['Re_min']**0.2 # Boundary layer momentum thickness : Empirical equation for turbulent plate
        t_blade = 0.12*stage.chord_R # Blade thickness estimate : Assumption for NACA 0012 airfoil
        lambda_2_rad = abs((self.Vel_Tri['beta3']+self.Vel_Tri['beta2'])/2)
        
        A = 1-(1+H_TE)*theta-t_TE/t_blade
        B = 1-H_TE*theta-t_TE/t_blade
        
        num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
        den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
        Yp = abs(1- num_Yp/den_Yp)

        # 5.2) Kacker-Okaapu : Secondary pressure losses
        Z = self.solidityStator*(self.Vel_Tri['beta2']-self.Vel_Tri['beta3'])/np.cos(self.Vel_Tri['beta3']) # Loading Factor
        Ys = abs(0.0334*1/self.params['AR_min']*(np.cos(self.Vel_Tri['alpha3'])/np.cos(self.Vel_Tri['beta2']))*Z)

        # Pressure loss 
        DP_loss = (Yp+Ys)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)*stage.static_states['D'][3]/2
        p0_out = stage.total_states['P'][2]-DP_loss
        
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.total_states['S'][3]
        
        hout = h0in-(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 3)
        
        pout_calc = stage.static_states['P'][3]

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][2])
        hout_s = self.AS.hmass()

        stage.eta_is_R = (stage.static_states['H'][2]-stage.static_states['H'][3])/(stage.static_states['H'][2]-hout_s)

        return (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

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
            
            # 1) Compute total inlet state
            hin = stage.static_states['H'][2]
            h0in = hin + (self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2  
            
            stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
            
            # 2) Compute A_flow and h_blade based on r_m guess
            stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
            stage.h_blade_S = stage.A_flow_S/(4*np.pi*self.r_m)
            
            # 3) Compute cord based on h_blade and AR
            stage.chord_S = stage.h_blade_S/stage.AR_S
            stage.pitch_S = self.solidityStator*stage.chord_S

            stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
            self.n_blade.append(stage.n_blade_S)
            
            # 5) Estimate pressure losses 
            # 5.1) Aungier : Profile pressure losses         
            D_e = 1.30*(stage.h_blade_S*2 * stage.pitch_S)**0.625 / (stage.h_blade_S*2 + stage.pitch_S)**0.25 # Equivalent diameter : Huebscher
            
            v_2 = np.sqrt(self.Vel_Tri["vm"]**2 + self.Vel_Tri["vu2"]**2)
            v_3 = np.sqrt(self.Vel_Tri["vm"]**2 + self.Vel_Tri["vu3"]**2)
            
            P_cst = np.cos(self.Vel_Tri["alpha3"])/2 * self.solidityStator * (v_2/v_3)**2 # Profile Constant
            
            Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
            # 5.2) Cohen : Endwall losses
            EW_Cst = np.cos((self.Vel_Tri["alpha2"]+self.Vel_Tri["alpha3"])/2)**3 / np.cos(self.Vel_Tri["alpha2"])**2  # Endwall Constant
    
            Yew = 0.02*(self.solidityStator/stage.AR_S)/EW_Cst
    
            # Pressure loss 
            DP_loss = (Yp+Yew)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)*stage.static_states['D'][2]/2
            p0_out = stage.total_states['P'][2]-DP_loss
            
            # Computation of static outlet pressure
            stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
            sout = stage.total_states['S'][3]
            
            hout = h0in-(self.Vel_Tri['vu3']**2 + self.Vel_Tri['vm']**2)/2
            stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 3)
            
            pout_calc = stage.static_states['P'][3]
    
            # Isentropic efficiency of the blade
            self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][2])
            hout_s = self.AS.hmass()
    
            stage.eta_is_S = (hout_s - stage.static_states['H'][2])/(stage.static_states['H'][3] - stage.static_states['H'][2])     
        else: # Rotor
                        
            # 1) Compute total inlet state
            hin = stage.static_states['H'][1]
            h0in = hin + (self.Vel_Tri['wu1']**2 + self.Vel_Tri['vm']**2)/2  
            
            stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
            
            # 2) Compute A_flow and h_blade based on r_m guess
            stage.A_flow_R = self.inputs['mdot']/(stage.static_states['D'][1]*self.Vel_Tri['vm'])
            stage.h_blade_R = stage.A_flow_R/(4*np.pi*self.r_m)
            
            # 3) If first stage : Compute max AR allowed in the compressor to estimate good AR for all stages
            
            if stage == self.stages[0]:
                chord_min = (self.params['Re_min']*stage.static_states['V'][1])/(stage.static_states['D'][1]*self.Vel_Tri['vm'])
                self.AR_max = stage.h_blade_R/chord_min
                
                self.AR = np.linspace(self.AR_max, self.params['AR_min'],self.nStages*2)
    
                c = 0
                
                for i in range(len(self.AR)):
                    if np.mod(i,2): # Stator
                    
                        self.stages[c].AR_S = self.AR[i]
                        c = c+1
                    else: # Rotor
                        self.stages[c].AR_R = self.AR[i]
            
            # 4) Compute cord based on h_blade and AR
            stage.chord_R = stage.h_blade_R/stage.AR_R
            stage.pitch_R = self.solidityRotor*stage.chord_R

            stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
            self.n_blade.append(stage.n_blade_R)
            
            # 5) Estimate pressure losses 
            # 5.1) Aungier : Profile pressure losses         
            D_e = 1.30*(stage.h_blade_R*2 * stage.pitch_R)**0.625 / (stage.h_blade_R*2 + stage.pitch_R)**0.25 # Equivalent diameter : Huebscher
            
            w_1 = np.sqrt(self.Vel_Tri["vm"]**2 + self.Vel_Tri["wu1"]**2)
            w_2 = np.sqrt(self.Vel_Tri["vm"]**2 + self.Vel_Tri["wu2"]**2)
            
            P_cst = np.cos(self.Vel_Tri["beta2"])/2 * self.solidityRotor * (w_1/w_2)**2 # Profile Constant
            
            Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
            # 5.2) Cohen : Endwall losses
            EW_Cst = np.cos((self.Vel_Tri["beta1"]+self.Vel_Tri["beta2"])/2)**3 / np.cos(self.Vel_Tri["beta1"])**2  # Endwall Constant
    
            Yew = 0.02*(self.solidityRotor/stage.AR_R)/EW_Cst
    
            # Pressure loss 
            DP_loss = (Yp+Yew)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu1']**2)*stage.static_states['D'][1]/2
            p0_out = stage.total_states['P'][1]-DP_loss
            
            # Computation of static outlet pressure
            stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
            sout = stage.total_states['S'][2]
            
            hout = h0in-(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2
            stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
            
            pout_calc = stage.static_states['P'][2]
    
            # Isentropic efficiency of the blade
            self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
            hout_s = self.AS.hmass()
    
            stage.eta_is_R = (hout_s - stage.static_states['H'][1])/(stage.static_states['H'][2] - stage.static_states['H'][1])
    
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
            
            self.r_tip.append(self.r_m + self.stages[i].h_blade_S/2)
            self.r_hub.append(self.r_m - self.stages[i].h_blade_S/2)
            self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
            self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
        
            self.r_tip.append(self.r_m + self.stages[i].h_blade_R/2)
            self.r_hub.append(self.r_m - self.stages[i].h_blade_R/2)
            self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
            self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
            
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
        
        Dh0 = self.inputs['W_dot_req']/self.inputs['mdot']
        
        self.eta_is = Dh0/Dh0s
        
        "------------- 2) Velocity Triangle Computation (+ Solidity) -------------------------------------" 
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
        self.Dh0Stage = Dh0/self.nStages
        self.Vel_Tri['u'] = np.sqrt(self.Dh0Stage/self.inputs['psi'])

        "------------- 5) Compute complete velocity triangles and exit losses ----------------------------" 

        # Compute velocity triangle with the value of u
        self.Vel_Tri['vm'] = self.Vel_Tri['vmOverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu2'] = self.Vel_Tri['vu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu3'] = self.Vel_Tri['vu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu2'] = self.Vel_Tri['wu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu3'] = self.Vel_Tri['wu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu1'] = self.Vel_Tri['vu3']
        self.Vel_Tri['wu1'] = self.Vel_Tri['wu3']

        self.exit_loss = (self.Vel_Tri['vm']**2+self.Vel_Tri['vu3']**2)/2

        # "------------- 6) Iterate on r_m and compute repeating stages ------------------------" 

        h_in = self.stages[0].total_states['H'][1] - (self.Vel_Tri['vm']**2)/2
        self.stages[0].update_static_AS(CP.HmassSmass_INPUTS, h_in, s_in, 1)

        def find_r_m(x):
            self.r_m = np.round(x[0],3)
            
            self.n_blade = []
            self.r_tip = []
            self.r_hub = []
            self.r_hub_tip = []
            self.r_ratio2 = []
    
            self.computeRepeatingStages()
        
            if self.r_hub_tip[-1] > 0: # Penalty to prevent converging to values not satisfying conditions on r_hub_tip
                penalty_1 = max(self.r_hub_tip[-1] - self.params['r_hub_tip_max'],0)*1000
                penalty_2 = max(self.params['r_hub_tip_min'] - self.r_hub_tip[0],0)*1000
                                
                return self.r_m + penalty_1 + penalty_2
            
            else: # A very high penalty prevents converging to r_m values very close to 0,  
                return self.r_m + 100
        
        # # sol = minimize(find_r_m, bounds=[(0.01, 10)], x0=0.1, tol = 1e-2)        
        sol = differential_evolution(find_r_m, bounds = self.params['r_m_bounds'], popsize = 5, strategy='best1bin', tol=1e-3, maxiter=100)        
        
        # # sol = differential_evolution(find_r_m,
        # #             bounds= self.params['r_m_bounds'],        # replace with your bounds
        # #             strategy='best1bin',    # fast, exploitative
        # #             popsize=8,              # smaller population
        # #             mutation=0.5,           # smaller step sizes
        # #             recombination=0.9,      # aggressive crossover
        # #             tol=1e-3,               # looser tolerance
        # #             maxiter=100,            # early stopping
        # #             workers=-1              # parallel execution
        # #         )
        
        # "------------- 7) Compute rotation speed and number of blades per stage ---------------------------" 

        # self.omega_rads = self.Vel_Tri['u']/self.r_m # rad/s
        # self.omega_RPM = self.omega_rads*60/(2*np.pi) 

        # self.n_blade = []

        # for stage in self.stages:
        #       stage.pitch_S = self.solidityStator*stage.chord_S
        #       stage.pitch_R = self.solidityRotor*stage.chord_R

        #       stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
        #       self.n_blade.append(stage.n_blade_S)

        #       stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
        #       self.n_blade.append(stage.n_blade_R)

        # "------------- 8) Compute main outputs -------------------------------------------------------------" 
        
        # hin = self.stages[0].total_states['H'][1]
        # hout = self.stages[-1].static_states['H'][3]
        
        # self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][3], self.stages[0].static_states['S'][1])

        # hout_s = self.AS.hmass()
        
        # self.W_dot = self.inputs['mdot']*(hin-hout)
                
        # self.eta_is = (hin - hout)/(hin - hout_s)

        # "------------- 9) Print Main Results -------------------------------------------------------------" 
        
        # print(f"Turbine mean diameter: {self.r_m} [m]")
        # print(f"Turbine rotation speed: {self.omega_RPM} [RPM]")
        # print(f"Turbine number of stage : {self.nStages} [-]")
        # print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")

        return

Comp = AxialCPMLLossCorrDesign('Cyclopentane')

# Cuerva Case

# Comp.set_inputs(
#     mdot = 53.52, # kg/s
#     W_dot_req = 3150*1e3, # W
#     p0_su = 1009*1e3, # Pa
#     T0_su = 273.15 + 182.3, # K
#     p_ex = 2493*1e3, # Pa
#     psi = 1, # [-] # 0.25
#     phi = 0.35, # [-]
#     R = 0.5, # [-] # 0.875
#     Mmax = 0.8 # [-]
#     )

# Comp.set_parameters(
#     Zweifel = 0.8, # [-]
#     Re_min = 5e5, # [-]
#     AR_min = 1, # [-]
#     r_hub_tip_max = 0.95, # [-]
#     r_hub_tip_min = 0.6, # [-]
#     )

# Zorlu Case

Comp.set_inputs(
    mdot = 19.24, # kg/s
    W_dot_req = 1187*1e3, # W
    p0_su = 468.4*1e3, # Pa
    T0_su = 273.15 + 138.3, # K
    p_ex = 1149*1e3, # Pa
    psi = 1, # [-]
    phi = 0.35, # [-]
    R = 0.5, # [-]
    Mmax = 0.8 # [-]
    )

Comp.set_parameters(
    Zweifel = 0.8, # [-]
    Re_min = 5e5, # [-]
    AR_min = 1, # [-]
    r_hub_tip_max = 0.95, # [-]
    r_hub_tip_min = 0.6, # [-]
    r_m_bounds = [(0.05,0.2)]
    )

Comp.design()

Comp.plot_geometry()
Comp.plot_n_blade()
Comp.plot_radius_verif()
Comp.plot_Mollier()