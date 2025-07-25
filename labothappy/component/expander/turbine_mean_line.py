# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 13:52:34 2025

@author: Basile
"""

"External modules"

import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import pandas as pd
import matplotlib.pyplot as plt

"Internal modules"
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model

class TurbMeanLine(BaseComponent):
    #%%
    
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
    
    
    
    
    #%%
    def __init__(self, fluid):
        super().__init__()

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

        # Define connectors: suction (su), exhaust (ex), and work (W_turb)
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_turb = WorkConnector()

    def get_required_inputs(self):
        """
        Determines and sets the required input values from connected components.
        Returns the list of all required input labels.
        """

        return ["su_p", "su_T", "ex_p", "N_rot", "su_fluid"]

    def get_required_parameters(self):
        """
        Returns a list of parameter names that must be defined to perform calculations.
        """
        return [
            "D_inlet",  # Inlet hydraulic diameter
            "N_turb_rated",  # Rated turbine speed
            "turb_voltage",  # Turbine electrical voltage
            "turb_phases",  # Number of turbine phases
            "eta_max_motor",  # Max motor (generator) efficiency
            "W_dot_el_rated",  # Rated electrical output
            "eta_m",  # Mechanical efficiency
            "eta_is_coefs",  # Polynomial coefficients for isentropic efficiency
            "eta_is_coefs_red",  # (Unused here but listed)
            "A_th",  # Turbine throat area
        ]

    def print_setup(self):
        """
        Prints the current configuration of the turbine component, including connectors, inputs, and parameters.
        """
        print("=== Turbine Setup ===")
        print("Connectors:")
        print(
            f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}"
        )
        print(
            f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}"
        )
        print(f"  - W_dot: speed={self.W_turb.N}")

        print("\nInputs:")
        for input in self.get_required_inputs():
            if input in self.inputs:
                print(f"  - {input}: {self.inputs[input]}")
            else:
                print(f"  - {input}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")

        print("======================")

    #%%

    def eta_el_fun(self, P_ratio):
        """
        Computes the electrical efficiency as a function of normalized power output (P_ratio).
        Uses two separate polynomial fits depending on the magnitude of P_ratio.
        """

        # Polynomial coefficients for high and low power ranges
        coefs_20p = [
            78.74503721229180,
            1.54402269709448,
            -0.04662069008665,
            0.00069559243591,
            -0.00000499382422,
            0.00000001349770,
        ]
        coefs_m20 = [
            0.82025554862776,
            4.78234707015054,
            0.73842411551209,
            -0.06398392686793,
            0.00134594665523,
        ]

        eta_ratio = 0

        # Use appropriate polynomial based on power ratio
        if P_ratio >= 20:
            for i in range(len(coefs_20p)):
                eta_ratio += coefs_20p[i] * P_ratio**i
        else:
            for i in range(len(coefs_m20)):
                eta_ratio += coefs_m20[i] * P_ratio**i

        return eta_ratio * self.params["eta_max_motor"]

    #%%
    
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

    #%%
    
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

        # 4) Compute cord and aspect ratio
                
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        stage.Y_vec_S = aungier_loss_model(self.Vel_Tri['beta1'], self.Vel_Tri['beta2'], self.Vel_Tri['beta1'], self.chord_S, stage.delta_tip, stage.e_blade)
        
        Y = stage.Y_vec_S[0]
        
        # aungier_loss_model(alpha1_pr, alpha2_pr, beta1, c, delta, D_lw, e, h, mu2, M1_pr, M2_pr, N_lw, R_c, rho2, s, t_max, t_TE, vm2, w2)
        
        # # 5) Estimate pressure losses 
        # # 5.1) Balje-Binsley
        # H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
        # t_TE = 5e-4 # m Tailing-edge blade thickness design estimate  - # !!! Limite de fabrication, voir 
        # theta = 0.036*stage.chord_S/self.params['Re_min']**0.2 # Boundary layer momentum thickness 
        # t_blade = 0.12*stage.chord_S # Blade thickness estimate : Assumption for NACA 0012 airfoil
        # lambda_2_rad = (self.Vel_Tri['beta2']+self.Vel_Tri['beta1'])/2
        
        # A = 1-(1+H_TE)*theta-t_TE/t_blade
        # B = 1-H_TE*theta-t_TE/t_blade
        
        # num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
        # den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
        # Yp = 1- num_Yp/den_Yp

        # # Secondary loss : Kacker-Okaapu
        # Z = self.solidityStator*(self.Vel_Tri['beta1']-self.Vel_Tri['beta2'])/np.cos(self.Vel_Tri['beta2']) # Loading Factor
        # Ys = abs(0.0334*1/stage.AR_S*(np.cos(self.Vel_Tri['alpha2'])/np.cos(self.Vel_Tri['beta1']))*Z)

        # # Pressure loss 
        # DP_loss = (Yp+Ys)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)*stage.static_states['D'][2]/2
        
        # p0_out = max(stage.total_states['P'][1]-DP_loss, self.inputs['p_ex']*0.6)
        
        # print(f"Rotor DP_loss : {DP_loss}")
                
        p0_out = stage.total_states['P'][1] - (Y*(stage.total_states['P'][1] - p_static_out))
        
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
        
        stage.Y_vec_R = aungier_loss_model(self.Vel_Tri['beta1'], self.Vel_Tri['beta2'], self.Vel_Tri['beta1'], self.chord_S, stage.delta_tip, stage.e_blade)
        
        Y = stage.Y_vec_R[0]
        
        # # 5) Estimate pressure losses 
        # # 5.1) Balje-Binsley : Profile pressure losses         
        # H_TE = 1.4 + 300/self.params['Re_min']**0.5 # Trailing-edge boundary layer shape factor : Aungier's Correlation for fully turbulent flow
        # t_TE = 5e-4 #  Tailing-edge blade thickness design estimate 
        # theta = 0.036*stage.chord_R/self.params['Re_min']**0.2 # Boundary layer momentum thickness : Empirical equation for turbulent plate
        # t_blade = 0.12*stage.chord_R # Blade thickness estimate : Assumption for NACA 0012 airfoil
        # lambda_2_rad = abs((self.Vel_Tri['beta3']+self.Vel_Tri['beta2'])/2)
        
        # A = 1-(1+H_TE)*theta-t_TE/t_blade
        # B = 1-H_TE*theta-t_TE/t_blade
        
        # num_Yp = (np.cos(lambda_2_rad)**2 * A**2) / B**2 + (np.sin(lambda_2_rad)**2) * B**2
        # den_Yp = 1 + 2 * (np.sin(lambda_2_rad)**2) * lambda_2_rad * (B**2 - A)
        # Yp = abs(1- num_Yp/den_Yp)

        # # 5.2) Kacker-Okaapu : Secondary pressure losses
        # Z = self.solidityRotor*(self.Vel_Tri['beta2']-self.Vel_Tri['beta3'])/np.cos(self.Vel_Tri['beta3']) # Loading Factor
        # Ys = abs(0.0334*1/stage.AR_R*(np.cos(self.Vel_Tri['alpha3'])/np.cos(self.Vel_Tri['beta2']))*Z)

        # # Pressure loss 
        # DP_loss = (Yp+Ys)*(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)*stage.static_states['D'][3]/2
        # p0_out = max(stage.total_states['P'][2]-DP_loss, self.inputs['p_ex']*0.6)
        
        # print(f"Rotor DP_loss : {DP_loss}")

        p0_out = stage.total_states['P'][2] - (Y*(stage.total_states['P'][2] - p_static_out))

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

    #%%
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
        
        return 

    #%%

    def solve(self):
        """
        Solves the turbine model: calculates outlet state, flow rate, mechanical and electrical work output.
        """

        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:

            # Calculate speed of sound at inlet using ideal gas approximation
            R = 8.3144
            R_M = R / PropsSI("M", self.su.fluid)
            cp_in, cv_in = PropsSI(
                ("C", "CVMASS"), "T", self.su.T, "P", self.su.p, self.su.fluid
            )
            gamma = cp_in / cv_in
            self.a = np.sqrt(gamma * R_M * self.su.T)

            # Mass flow rate using choked flow assumption
            self.m_dot = self.a * self.su.D * self.params["A_th"]
            self.su.set_m_dot(self.m_dot)

            # Calculate isentropic efficiency from the empirical map
            self.eta_is = self.eta_is_turb()

            # Compute isentropic enthalpy at exhaust pressure and entropy
            h_ex_is = PropsSI("H", "P", self.ex.p, "S", self.su.s, self.su.fluid)

            # Actual enthalpy based on isentropic efficiency
            h_ex = self.su.h - (self.eta_is / 100) * (self.su.h - h_ex_is)

            # Set exhaust state
            self.ex.set_h(h_ex)
            self.ex.set_m_dot(self.su.m_dot)

            # Compute shaft work
            self.W_dot = (self.su.h - h_ex) * self.su.m_dot

            # Electrical efficiency and electrical output
            self.eta_el = (
                self.eta_el_fun(100 * self.W_dot / self.params["W_dot_el_rated"]) / 100
            )
            self.W_el = self.W_dot * self.eta_el

            # Set mechanical output in the work connector
            self.W_turb.set_W_dot(self.W_dot)

            self.defined = True

        else:
            if not self.calculable:
                print("Input of the component not completely known")
            if not self.parametrized:
                print("Parameters of the component not completely known")
