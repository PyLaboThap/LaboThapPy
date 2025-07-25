#!/usr/bin/python3

# --- loading libraries 

from connector.mass_connector import MassConnector
from correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution
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
            
            self.stage = None
            self.AR = None
            
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
    def compute_rotor_t_max(self, stage):
                
        def z_Saravanamuttoo(xhi, t, c):
            
            B_values = np.array([951.083, 904.889, 829.035, 753.251, 665.223, 568.022, 461.654, 359.362, 273.379])
            n_values = np.array([1.694, 1.605, 1.520, 1.436, 1.347, 1.268, 1.179, 1.094, 1])
            xhi_values = np.array([40,50,60,70,80,90,100,110,120])
            
            B = np.interp(xhi, xhi_values, B_values)
            n = np.interp(xhi, xhi_values, n_values)
            
            z = 1/B * (10*t/c)**n 
            
            return z
        
        def system(x):

            stage.t_blade_R = x[0]
            # Maximum Bending Stress : from Saravanamuttoo et. al. - Gas Turbine Theory
            xhi = abs(180*(self.Vel_Tri["alpha2"] - self.Vel_Tri["alpha3"])/np.pi) # Rotor blade camber angle
             
            z = z_Saravanamuttoo(xhi, stage.t_blade_R, stage.chord_R)
            stage.sigma_bmax = (1/(z*stage.chord_R**3)) * (stage.h_blade_R/2) * self.inputs['mdot']*self.Vel_Tri['vm']*(np.tan(self.Vel_Tri['alpha2'])+np.tan(self.Vel_Tri['alpha3']))/stage.n_blade_R
    
            return stage.sigma_bmax - (130*1e6*2)

        # Initial guess (e.g., 5 mm thickness)
        x0 = [0.005]  # in meters
        
        # Solve
        solution = fsolve(system, x0)

    def compute_stator_t_max(self, stage):
                
        def z_Saravanamuttoo(xhi, t, c):
            
            B_values = np.array([951.083, 904.889, 829.035, 753.251, 665.223, 568.022, 461.654, 359.362, 273.379])
            n_values = np.array([1.694, 1.605, 1.520, 1.436, 1.347, 1.268, 1.179, 1.094, 1])
            xhi_values = np.array([40,50,60,70,80,90,100,110,120])
            
            B = np.interp(xhi, xhi_values, B_values)
            n = np.interp(xhi, xhi_values, n_values)
            
            z = 1/B * (10*t/c)**n 
            
            return z
        
        def system(x):

            stage.t_blade_S = x[0]
            # Maximum Bending Stress : from Saravanamuttoo et. al. - Gas Turbine Theory
            xhi = abs(180*(self.Vel_Tri["alpha1"] - self.Vel_Tri["alpha2"])/np.pi) # Rotor blade camber angle
             
            z = z_Saravanamuttoo(xhi, stage.t_blade_S, stage.chord_S)
            stage.sigma_bmax = (1/(z*stage.chord_S**3)) * (stage.h_blade_S/2) * self.inputs['mdot']*self.Vel_Tri['vm']*(np.tan(self.Vel_Tri['alpha1'])+np.tan(self.Vel_Tri['alpha2']))/stage.n_blade_S
    
            return stage.sigma_bmax - (130*1e6*2)

        # Initial guess (e.g., 5 mm thickness)
        x0 = [0.005]  # in meters
        
        # Solve
        solution = fsolve(system, x0)
        return 

    def stator_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]

        # 1) Guess outlet state
        h_static_out = x[0]*1e5
        p_static_out = x[1]*1e5
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu2']**2)
        stage.M_S = max(v,w)/stage.static_states['A'][2]
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(4*np.pi*self.r_m)

        # 4) Compute cord, aspect ratio, blade pitch and blade number
                
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        stage.pitch_S = self.solidityStator*stage.chord_S
        stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
        
        self.compute_stator_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha2'] - self.Vel_Tri['alpha1']
        
        R_c = 2*np.sin(abs(camber)/2)*stage.chord_S # Geometrical estimation of blade curvature radius
                
        stage.M1_S = self.Vel_Tri['v1']/stage.static_states['A'][1]
        stage.M2_S = self.Vel_Tri['v2']/stage.static_states['A'][2]
        
        stage.Y_vec_S = aungier_loss_model(self.Vel_Tri['alpha1'], self.Vel_Tri['alpha2'], self.Vel_Tri['alpha1'], stage.chord_S, 
                               0, self.params['D_lw'], self.params['e_blade'], stage.h_blade_S, stage.static_states['V'][2], 
                               stage.M1_S, stage.M2_S, self.params['N_lw'], R_c, stage.static_states['D'][2], stage.pitch_S, stage.t_blade_S, self.params['t_TE'],
                               self.Vel_Tri['vm'], self.Vel_Tri['v2'],1)
        
        # stage.Y_vec_S = aungier_loss_model(self.Vel_Tri['beta1'], self.Vel_Tri['beta2'], self.Vel_Tri['beta1'], stage.chord_S, stage.delta_tip, stage.e_blade)    
        Y = stage.Y_vec_S['Y_tot']
                        
        # p0_out = stage.total_states['P'][1] - (Y*(stage.total_states['P'][1] - p_static_out))
        p0_out = (stage.total_states['P'][1] + Y * p_static_out)/(1+Y)

        
        # print(f"Rotor DP_loss : {DP_loss}")
        
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

        return np.array([hout, pout_calc])*1e-5 # (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def rotor_blade_row_system(self, x):
        
        stage = self.stages[self.curr_stage_index]
        
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x*1e5
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)
        
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu3']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)
        stage.M_R = max(v,w)/stage.static_states['A'][3]
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][2]
        h0in = hin + (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_R = self.inputs['mdot']/(stage.static_states['D'][3]*self.Vel_Tri['vm'])
        stage.h_blade_R = stage.A_flow_R/(4*np.pi*self.r_m)
        
        # 4) Compute cord and aspect ratio
        stage.chord_R = (self.params['Re_min']*stage.static_states['V'][3])/(stage.static_states['D'][3]*self.Vel_Tri['vm'])            
        stage.AR_R = stage.h_blade_R/stage.chord_R
        
        stage.pitch_R = self.solidityRotor*stage.chord_R
        stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
        
        self.compute_rotor_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha3'] - self.Vel_Tri['alpha2']
        
        R_c = 2*np.sin(abs(camber)/2)*stage.chord_R # Geometrical estimation of blade curvature radius
                
        stage.M2_R = self.Vel_Tri['w2']/stage.static_states['A'][2]
        stage.M3_R = self.Vel_Tri['w3']/stage.static_states['A'][3]
        
        stage.Y_vec_R = aungier_loss_model(self.Vel_Tri['beta2'], self.Vel_Tri['beta3'], self.Vel_Tri['beta2'], stage.chord_R, 
                               self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.static_states['V'][3], 
                               stage.M2_R, stage.M3_R, self.params['N_lw'], R_c, stage.static_states['D'][3], stage.pitch_R, stage.t_blade_R, self.params['t_TE'],
                               self.Vel_Tri['vm'], self.Vel_Tri['w3'],1)

        Y = stage.Y_vec_R['Y_tot']
                
        p0_out = (stage.total_states['P'][2] + Y * p_static_out)/(1+Y)

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

        return np.array([hout, pout_calc])*1e-5 # (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def last_blade_row_system(self, x):
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x*1e5
        
        stage = self.stages[-1]
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        v = np.sqrt(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu2']**2)
        w = np.sqrt(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['wu2']**2)
        stage.M_S = max(v,w)/stage.static_states['A'][2]
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (self.Vel_Tri_Last_Stage['vu1']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri_Last_Stage['vm'])
        stage.h_blade_S = stage.A_flow_S/(4*np.pi*self.r_m)

        # 4) Compute cord and aspect ratio
                
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri_Last_Stage['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        # 5) Estimate pressure losses 
        # 5.1) Aungier : Profile pressure losses                     
        v_1 = np.sqrt(self.Vel_Tri_Last_Stage["vm"]**2 + self.Vel_Tri_Last_Stage["vu1"]**2)
        v_2 = np.sqrt(self.Vel_Tri_Last_Stage["vm"]**2 + self.Vel_Tri_Last_Stage["vu2"]**2)
        
        a = 0.0117 # NACA blade - 0.007 : C.4 circular-arc blade
        
        alpha = 0
        alpha_star = 0
        
        D_e = (np.cos(self.Vel_Tri_Last_Stage['alpha2'])/np.cos(self.Vel_Tri_Last_Stage['alpha1']))*(1.12+a*(alpha - alpha_star)+0.61*np.cos(self.Vel_Tri_Last_Stage['alpha1'])**2 / self.solidityStator * (np.tan(self.Vel_Tri_Last_Stage['alpha1'])-np.tan(self.Vel_Tri_Last_Stage['alpha2'])))
        
        P_cst = np.cos(self.Vel_Tri_Last_Stage["alpha2"])/2 * self.solidityStator * (v_1/v_2)**2 # Profile Constant
        
        Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
        # 5.2) Cohen : Endwall losses
        EW_Cst = np.cos((self.Vel_Tri_Last_Stage["alpha1"]+self.Vel_Tri_Last_Stage["alpha2"])/2)**3 / np.cos(self.Vel_Tri_Last_Stage["alpha1"])**2  # Endwall Constant

        Yew = 0.02*(self.solidityStator/stage.AR_S)/EW_Cst

        # Pressure loss 
        DP_loss = (Yp+Yew)*(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu1']**2)*stage.static_states['D'][1]/2
        p0_out = stage.total_states['P'][1]-DP_loss
                
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(self.Vel_Tri_Last_Stage['vu2']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
        
        pout_calc = stage.static_states['P'][2]

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        # print(f"h0in: {h0in}")
        # print(f"h1: {stage.static_states['H'][1]}")
        # print(f"kinetic1: {(self.Vel_Tri_Last_Stage['vu1']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2}")
        # print(f"h2: {stage.static_states['H'][2]}")
        # print(f"kinetic2: {(self.Vel_Tri_Last_Stage['vu2']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2}")

        return np.array([hout, pout_calc])*1e-5 # return (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

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
        self.Vel_Tri['w3'] = np.sqrt(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)

        # Angles in radians
        self.Vel_Tri_Last_Stage['alpha1'] = self.Vel_Tri['alpha3'] 
        self.Vel_Tri_Last_Stage['alpha2'] = 0

        self.Vel_Tri_Last_Stage['beta1'] = self.Vel_Tri['beta3']
        self.Vel_Tri_Last_Stage['beta2'] = np.arctan(self.Vel_Tri['u']/self.Vel_Tri['vm'])
        
        return 
    
    def computeBladeRow(self,stage_index,row_type):
        stage = self.stages[stage_index]

        self.curr_stage_index = stage_index
               
        if row_type == 'S': # Stator
            
            # print("Stator")
        
            RP_1_row = (self.inputs['p0_su']/self.inputs['p_ex'])**(1/(2*self.nStages))
            h_out_guess = stage.static_states['H'][1] - self.Dh0Stage/2
            pout_guess = stage.static_states['P'][1]/RP_1_row
            # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])         
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
            
            res = 1
            x_in = x0_disc
            
            c = 0
            
            while res > 1e-6:
                
                if c > 100:
                    exit()
                
                # print(f"x_in : {x_in}")
                
                x_out = self.stator_blade_row_system(x_in)

                # print(f"x_out : {x_out}")
                
                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                              
                # print(f"new x_in : {x_in}")

                c += 1
                
            self.stator_blade_row_system(x_out)
            
        else: # Rotor

            # print("Rotor")

            RP_1_row = (self.inputs['p0_su']/self.inputs['p_ex'])**(1/(2*self.nStages))
            h_out_guess = stage.static_states['H'][2] - self.Dh0Stage/2
            pout_guess = stage.static_states['P'][2]/RP_1_row
            # sol = minimize(self.rotor_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])    
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

            res = 1
            x_in = x0_disc
            
            c = 0
            
            while res > 1e-6:

                if c > 100:
                    exit()                

                # print(f"x_in : {x_in}")

                x_out = self.rotor_blade_row_system(x_in) 

                # print(f"x_out : {x_out}")

                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                       
                # print(f"new x_in : {x_in}")

                c += 1
            
            self.rotor_blade_row_system(x_out)
                        
        return
            
    def computeRepeatingStages(self):
        
        for i in range(int(self.nStages)):
            if i == 0:
                self.computeBladeRow(i, 'S')
                self.computeBladeRow(i, 'R')

            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
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
        
        h_out_guess = self.stages[-2].total_states['H'][3]
        pout_guess = self.stages[-2].total_states['P'][3]
        # sol = minimize(self.last_blade_row_system, x0=(h_out_guess,pout_guess), bounds=[(self.stages[-1].static_states['H'][1], h_out_guess), (self.stages[-1].static_states['P'][1], pout_guess)])         
        
        # Initial guess vector
        x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

        res = 1
        x_in = x0_disc
        
        c = 0
            
        while res > 1e-6:

            if c > 100:
                exit()                

            x_out = self.last_blade_row_system(x_in) 

            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            
            x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                        
            c += 1
        
        self.last_blade_row_system(x_out)
        
        return
    
    def design_system(self, x):
        self.penalty = -1
        # try:
            # print(x)
        self.reset()
        
        self.psi = x[0]
        self.phi = x[1]
        self.R = x[2]
        self.r_m = x[3]
        self.M_max = x[4]
        
        # print(f"pos : [{self.psi, self.phi, self.R, self.r_m}]")
        
        self.r_tip = []
        self.r_hub = []
        self.r_hub_tip = []
        self.r_ratio2 = []
        
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
        
        vMax = self.AS.speed_sound() * self.M_max
        
        # Assume u based on the maximum speed
        self.Vel_Tri['u'] = vMax / max([self.Vel_Tri['v2OverU'],self.Vel_Tri['w3OverU']])
        
        "------------- 4) Compute number of stage + recompute u  -----------------------------------------" 
        
        # Compute required number of stages based on assumed u
        Dh0Stage = self.psi * self.Vel_Tri['u']**2
        self.nStages = int(round(Dh0/Dh0Stage))
        
        for i in range(self.nStages-1):
            self.stages.append(self.stage(self.fluid))
        
        # Recompute u based on the number of stages to satisfy the work. As r_m is constant, u is contant accross stages
        self.Dh0Stage = Dh0/self.nStages
        self.Vel_Tri['u'] = np.sqrt(self.Dh0Stage/self.psi)

        "------------- 5) Compute complete velocity triangles and exit losses ----------------------------" 

        # Compute velocity triangle with the value of u
        self.Vel_Tri['vm'] = self.Vel_Tri['vmOverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu2'] = self.Vel_Tri['vu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu3'] = self.Vel_Tri['vu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu2'] = self.Vel_Tri['wu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu3'] = self.Vel_Tri['wu1'] = self.Vel_Tri['wu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu1'] = self.Vel_Tri['vu3']

        self.Vel_Tri['v1'] = np.sqrt(self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['v2'] = np.sqrt(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['v3'] = np.sqrt(self.Vel_Tri['vu3']**2 + self.Vel_Tri['vm']**2)

        self.Vel_Tri['w1'] = np.sqrt(self.Vel_Tri['wu1']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w2'] = np.sqrt(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w3'] = np.sqrt(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)

        self.exit_loss = (self.Vel_Tri['vm']**2+self.Vel_Tri['vu3']**2)/2
        self.exit_loss_W = self.exit_loss*self.inputs['mdot']

        "------------- 6) Iterate on r_m and compute repeating stages ------------------------" 

        h_in = self.stages[0].total_states['H'][1] - (self.Vel_Tri['vm']**2)/2
        self.stages[0].update_static_AS(CP.HmassSmass_INPUTS, h_in, s_in, 1)
        
        try:
            self.computeRepeatingStages()        
        except:
            print("Error in stages")
            obj = 10000
            return obj
        
        "------------- 7) Compute rotation speed and number of blades per stage ---------------------------" 

        self.omega_rads = self.Vel_Tri['u']/self.r_m # rad/s
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

        self.stages.append(self.stage(self.fluid))
        
        self.computeVelTriangleLastStage()
        
        # try:
        self.computeLastStage()
        # except:
        #     print("Error in last stage")
        #     obj = 10000
        #     return obj

        "------------- 9) Compute main outputs -------------------------------------------------------------" 
        
        hin = self.stages[0].total_states['H'][1]
        hout = self.stages[-1].static_states['H'][2]
        
        self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])

        hout_s = self.AS.hmass()
        
        self.W_dot = self.inputs['mdot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - hout_s)


        penalty_1 = max(self.r_hub_tip[0] - self.params['r_hub_tip_max'],0)*100
        penalty_2 = max(self.params['r_hub_tip_min'] - self.r_hub_tip[-1],0)*100
        
        if abs((self.inputs["p_ex"] - self.stages[-1].static_states['P'][2])/self.inputs["p_ex"]) >= 5e-2:
            penalty_3 = abs((self.inputs["p_ex"] - self.stages[-1].static_states['P'][2])/self.inputs["p_ex"])*10
            self.Pressure_Deviation = self.inputs["p_ex"] - self.stages[-1].static_states['P'][2]
        else:
            penalty_3 = 0
            self.Pressure_Deviation = self.inputs["p_ex"] - self.stages[-1].static_states['P'][2]

        self.penalty = penalty_1 + penalty_2 + penalty_3

        if self.eta_is > 0 and self.eta_is <= 1:
            if self.penalty == 0:
                if self.W_dot > self.inputs['W_dot']:
                    self.inputs['W_dot'] = self.W_dot
                    
            obj = -self.eta_is + self.penalty
            # print(f"opt 'success' : {obj}")
        else:
            # print("Bad eta_is")
            obj = 10000

        
        print(obj)
        
        if obj < 10000:
            self.allowable_positions.append([self.psi, self.phi, self.R, self.r_m])
        
        return obj

        # except:
        #     # print("FAIL")
        #     return 10000

    # def design(self):
    #     bounds = (np.array([
    #         self.params['psi_bounds'][0],
    #         self.params['phi_bounds'][0],
    #         self.params['R_bounds'][0],
    #         self.params['r_m_bounds'][0],
    #     ]),
    #     np.array([
    #         self.params['psi_bounds'][1],
    #         self.params['phi_bounds'][1],
    #         self.params['R_bounds'][1],
    #         self.params['r_m_bounds'][1],
    #     ]))
    
    #     # Wrap your function to return a scalar
    #     def objective_wrapper(x):
    #         return np.array([self.design_system(xi) for xi in x])
    
    #     optimizer = ps.single.GlobalBestPSO(
    #         n_particles=10,
    #         dimensions=4,
    #         options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
    #         bounds=bounds
    #     )
    
    #     cost, pos = optimizer.optimize(objective_wrapper, iters=30, verbose=True)

    def design(self):
        bounds = (np.array([
            self.params['psi_bounds'][0],
            self.params['phi_bounds'][0],
            self.params['R_bounds'][0],
            self.params['r_m_bounds'][0],
            self.params['M_1_st_bounds'][0],
        ]),
        np.array([
            self.params['psi_bounds'][1],
            self.params['phi_bounds'][1],
            self.params['R_bounds'][1],
            self.params['r_m_bounds'][1],
            self.params['M_1_st_bounds'][1],
        ]))
    
        def objective_wrapper(x):
            # Round inputs if needed for discrete steps
            rounded_x = np.copy(x)
            # (Optional rounding logic goes here)
            return np.array([self.design_system(xi) for xi in rounded_x])
    
        # Initialize the optimizer
        optimizer = ps.single.GlobalBestPSO(
            n_particles=20,
            dimensions=5,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
        
        self.allowable_positions = []
        
        # Custom stopping logic
        patience = 5
        tol = 1e-3
        max_iter = 20
        no_improve_counter = 0
        best_cost = np.inf
    
        for i in range(max_iter):
            # One iteration step
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
    
            current_best = optimizer.swarm.best_cost
    
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
    
            print(f"[{i+1}] Best cost: {best_cost:.6f}")
    
            if no_improve_counter >= patience:
                print("Stopping early due to stagnation.")
                break
             
        best_pos = optimizer.swarm.best_pos
             
        self.design_system(best_pos)
             
        "------------- Print Main Results -------------------------------------------------------------" 
        
        print(f"Parameters : {self.psi, self.phi, self.R, self.r_m, self.M_max}")
        
        print(f"Turbine mean radius: {self.r_m} [m]")
        print(f"Turbine rotation speed: {self.omega_RPM} [RPM]")
        print(f"Turbine number of stage : {self.nStages} [-]")
        print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")        
        print(f"Turbine Generation : {self.W_dot} [W]")        
        
        return best_pos

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
        damping = 0.5 # [-]
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
        Re_min = 5e5, # [-]
        AR_min = 1, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        r_m_bounds = [0.16,0.4], # [m]
        psi_bounds = [0.8,1.4], # [-]
        phi_bounds = [0.5,0.7], # [-]
        R_bounds = [0.49,0.51], # [-]
        damping = 0.5 # [-]
        )
    
elif case_study == 'TCO2_ORC':

    Turb = AxialTurbineMeanLineDesign('CO2')

    Turb.set_inputs(
        mdot = 5*100, # kg/s
        W_dot = 5*4.69*1e6, # W
        p0_su = 140*1e5, # Pa
        T0_su = 273.15 + 121, # K
        p_ex = 39.8*1e5, # Pa
        Mmax = 0.5 # [-]
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        Re_min = 5e6, # [-]
        AR_min = 1, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        r_m_bounds = [0.1,0.3], # [m]
        psi_bounds = [1,2.5], # [-]
        phi_bounds = [0.4,0.7], # [-] 
        R_bounds = [0.49,0.51], # [-]
        M_1_st_bounds = [0.4, 0.7], # [-]
        damping = 0.2, # [-]
        delta_tip = 0.2*1e-3, # [m] : tip clearance
        N_lw = 0, # [-] : Number of lashing wires
        D_lw = 0, # [m] : Diameter of lashing wires
        e_blade = 0.002*1e-3, # [m] : blade roughness
        t_TE = 0.5*1e-3 # [m] : trailing edge thickness
        )
    

best_pos = Turb.design()

Turb.plot_geometry()
Turb.plot_n_blade()
Turb.plot_radius_verif()
Turb.plot_Mollier()
