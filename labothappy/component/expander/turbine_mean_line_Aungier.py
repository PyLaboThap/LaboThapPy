from connector.mass_connector import MassConnector
from correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector

from toolbox.turbomachinery.mean_line_axial_turbine_mapping import map_plot, plot_power_eta_vs_mdot
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution
import pyswarms as ps

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

class AxialTurbineMeanLine(BaseComponent):

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
        self.Vel_Tri_Last_Stage = {}

        self.su = MassConnector()
        self.ex = MassConnector()
        
        self.Dh0_stage_guess = 0

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
            
            self.xhi_S1 = None
            self.xhi_S2 = None
            
            self.xhi_R1 = None
            self.xhi_R2 = None
            
            self.Vel_Tri_R = {}
            self.Vel_Tri_S = {}
            
        def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.total_states.loc[position, 'H'] = self.AS.hmass()
            self.total_states.loc[position, 'S'] = self.AS.smass()
            self.total_states.loc[position, 'P'] = self.AS.p()
            self.total_states.loc[position, 'D'] = self.AS.rhomass()
            try:
                self.total_states.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                self.total_states.loc[position, 'A'] = -1
            self.total_states.loc[position, 'V'] = self.AS.viscosity()        
            
            return
        
        def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.static_states.loc[position, 'H'] = self.AS.hmass()
            self.static_states.loc[position, 'S'] = self.AS.smass()
            self.static_states.loc[position, 'P'] = self.AS.p()
            self.static_states.loc[position, 'D'] = self.AS.rhomass()
            try:
                self.static_states.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                self.static_states.loc[position, 'A'] = -1
            self.static_states.loc[position, 'V'] = self.AS.viscosity()         

            return

    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_stage_parameters(self, **parameters):
        """
        Assign stage parameters from arrays.
        If a stage doesn't exist, instantiate it.
        """
        # number of stages required = longest array among all parameters
        n_stages = max(len(arr) for arr in parameters.values())
    
        # ensure self.stages exists and has enough elements
        if not hasattr(self, "stages"):
            self.stages = []
        while len(self.stages) < n_stages:
            # ⬇️ replace Stage() with your actual Stage class constructor
            self.stages.append(self.stage(self.fluid))
    
        # assign parameters
        for key, array in parameters.items():
            for i in range(len(array)):
                setattr(self.stages[i], key, array[i])
                    
    # ---------------- Loss Models ------------------------------------------------------------------------

    def stator_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]

        # 1) Guess outlet state
        h_static_out = x[0]*1e5
        p_static_out = x[1]*1e5
                
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        stage.Vel_Tri_S['u'] = u = self.u
        
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][2]*A_flow)
        
        if self.curr_stage_index == 0:
            stage.Vel_Tri_S['alpha1'] = alpha1 = stage.xhi_S1
            stage.Vel_Tri_S['vu1'] = vu1 = vm*np.tan(alpha1)
        else:
            stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
            stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u 
            stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
        
        
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][1]

        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        h02 = h0in
        
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][2]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)

        # 4) Compute cord, aspect ratio, blade pitch and blade number
        
        stage.Re_s = stage.chord_S*(stage.static_states['D'][2]*vm)/(stage.static_states['V'][2])
        stage.AR_S = stage.h_blade_S/stage.chord_S
                
        # 5) Loss model
        
        stage.beta_g_S = np.arcsin(stage.o_S/stage.pitch_S)
        
        stage.Y_vec_S = aungier_loss_model(stage.Vel_Tri_S['alpha1'], stage.Vel_Tri_S['alpha2'], stage.beta_g_S*180/np.pi, stage.xhi_S1, stage.chord_S, 
                               0, self.params['D_lw'], self.params['e_blade'], stage.h_blade_S, stage.static_states['V'][2], 
                               stage.M1_S, stage.M2_S, self.params['N_lw'], stage.R_c_S, stage.static_states['D'][2], stage.pitch_S, stage.t_blade_S, stage.t_TE_S,
                               vm, v2,1)
        
        self.compute_deviation_stator(stage)
        alpha2_calc = stage.xhi_S2 + stage.delta_S
        
        v2_new = vm/np.cos(alpha2_calc)
        
        Y = stage.Y_vec_S['Y_tot']
        p0_out = (stage.total_states['P'][1] + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = stage.total_states['H'][2]-v2_new**2/2
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
        
        stage.Vel_Tri_R['u'] = u = self.u
        
        A_flow = stage.h_blade_R*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_R['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][3]*A_flow)
        stage.Vel_Tri_R['vu2'] = vu2 = vm*np.tan(stage.Vel_Tri_R['alpha2'])    
        
        stage.Vel_Tri_R['wu2'] = wu2 = vu2 - u
        stage.Vel_Tri_R['w2'] = w2 = np.sqrt(wu2**2 + vm**2)
        stage.M2_R = w2/stage.static_states['A'][2]
        stage.Vel_Tri_R['beta2'] = np.arctan(wu2/vm)
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][2]
        h0in = hin + (w2**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        h03 = stage.total_states['H'][2]
        stage.Vel_Tri_R['w3'] = w3 = np.sqrt(2*(h03 - h_static_out))
        stage.M3_R = w3/stage.static_states['A'][3]
        stage.Vel_Tri_R['beta3'] = -np.arccos(vm/w3)
                
        # 4) Compute cord, aspect ratio, pitch and blade number
        stage.Re_r = stage.chord_R*(stage.static_states['D'][3]*vm)/(stage.static_states['V'][3])
        stage.AR_R = stage.h_blade_R/stage.chord_R
                
        # 5) Loss model
        
        stage.beta_g_R =  np.arcsin(stage.o_R/stage.pitch_R)  # mid-passage metal angle        
                
        stage.Y_vec_R = aungier_loss_model(-stage.Vel_Tri_R['beta2'], -stage.Vel_Tri_R['beta3'], stage.beta_g_R*180/np.pi, -stage.xhi_R1, stage.chord_R, 
                               self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.static_states['V'][3], 
                               stage.M2_R, stage.M3_R, self.params['N_lw'], stage.R_c_R, stage.static_states['D'][3], stage.pitch_R, stage.t_blade_R, stage.t_TE_R,
                               vm, w3,1)

        self.compute_deviation_rotor(stage)
        beta3_calc = stage.xhi_R2 + stage.delta_R

        w3_new = vm/np.cos(beta3_calc)

        Y = stage.Y_vec_R['Y_tot']
                
        p0_out = (stage.total_states['P'][2] + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.total_states['S'][3]
        
        hout = h0in-(w3_new**2)/2
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
        
        stage.Vel_Tri_S['u'] = u = self.u
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][2]*A_flow)
        
        stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
        stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u 
        stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
        
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][1]

        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)    
                
        h02 = h0in
        
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][2]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)
        
        # 5) Estimate pressure losses 
        AR_S = stage.h_blade_S/stage.chord_S
        solidity = (stage.chord_S/stage.pitch_S)
        
        a = 0.0117 # NACA blade - 0.007 : C.4 circular-arc blade
        
        D_e = (np.cos(alpha2)/np.cos(alpha1))*(1.12+a*(alpha1 - stage.xhi_S1)+0.61*np.cos(alpha1)**2 / solidity * (np.tan(alpha1)-np.tan(alpha2)))
        
        P_cst = np.cos(alpha2)/2 * solidity * (v1/v2)**2 # Profile Constant
        
        Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
        # 5.2) Cohen : Endwall losses
        EW_Cst = np.cos((alpha1+alpha2)/2)**3 / np.cos(alpha1)**2  # Endwall Constant

        Yew = 0.02*(solidity/AR_S)/EW_Cst

        # Pressure loss 
        DP_loss = (Yp+Yew)*(v1**2)*stage.static_states['D'][1]/2
        p0_out = stage.total_states['P'][1]-DP_loss
                
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(v2**2)/2
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

    def compute_deviation_stator(self, stage):
        
        delta_0S = np.arcsin((stage.o_S/stage.pitch_S)*(1+(1-stage.o_S/stage.pitch_S)*(2*stage.beta_g_S/np.pi)**2)) - stage.beta_g_S
        
        if stage.M2_S <= 0.5:
            stage.delta_S = delta_0S
        else:
            X = 2*stage.M2_S-1
            stage.delta_S = delta_0S*(1-10*X**3 + 15*X**4 - 6*X**5)
        
        return 

    def compute_deviation_rotor(self, stage):
                
        delta_0R = np.arcsin((stage.o_R/stage.pitch_R)*(1+(1-stage.o_R/stage.pitch_R)*(2*stage.beta_g_R/np.pi)**2)) - abs(stage.beta_g_R)
        
        if stage.M3_R <= 0.5:
            stage.delta_R = delta_0R
        else:
            X = 2*stage.M3_R-1
            stage.delta_R = delta_0R*(1-10*X**3 + 15*X**4 - 6*X**5)

        return 

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
        
            RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
            
            if self.Dh0_stage_guess !=0:
                h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2    
            else:
                h_out_guess = stage.static_states['H'][1]*0.99
                
            pout_guess = stage.static_states['P'][1]/RP_1_row
            # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])         
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
            
            res = 1
            x_in = x0_disc
            
            c = 0
            
            while res > 1e-8:
                
                if c > 1000:
                    raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")
                
                # print(f"x_in : {x_in}")
                
                x_out = self.stator_blade_row_system(x_in)

                # print(f"x_out : {x_out}")
                
                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                              
                # print(f"new x_in : {x_in}")

                c += 1
                
            self.stator_blade_row_system(x_out)
                        
            # print(f'Y_S : {stage.Y_vec_S}')

        else: # Rotor

            # print("Rotor")

            RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
            
            if self.Dh0_stage_guess !=0:
                h_out_guess = stage.static_states['H'][2] - self.Dh0_stage_guess/2    
            else:
                h_out_guess = stage.static_states['H'][2]*0.99
            
            pout_guess = stage.static_states['P'][2]/RP_1_row
            # sol = minimize(self.rotor_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])    
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

            res = 1
            x_in = x0_disc
            
            c = 0
            
            while res > 1e-8:

                if c > 1000:
                    raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")

                # print(f"x_in : {x_in}")

                x_out = self.rotor_blade_row_system(x_in) 

                # print(f"x_out : {x_out}")

                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                       
                # print(f"new x_in : {x_in}")

                c += 1
            
            self.rotor_blade_row_system(x_out)
            self.compute_deviation_rotor(stage)
                        
            # print(f'Y_R : {stage.Y_vec_R}')

        return
            
    def computeRepeatingStages(self):
                
        self.nStages = self.params['nStages']
        
        for i in range(int(self.nStages)):
                    
            if i == 0:
                self.computeBladeRow(i, 'S')
                
                self.compute_deviation_stator(self.stages[i])
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']
                
                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

                self.Dh0_stage_guess = self.stages[i].total_states['H'][1] - self.stages[i].total_states['H'][3]

            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
                self.computeBladeRow(i, 'S')
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']

                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

        return
    
    def computeLastStage(self):
        stage = self.stages[-1]
        
        stage.static_states.loc[1] = self.stages[-2].static_states.loc[3]
        
        RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))     
        
        h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2  
        pout_guess = stage.static_states['P'][1]/RP_1_row
        # sol = minimize(self.last_blade_row_system, x0=(h_out_guess,pout_guess), bounds=[(self.stages[-1].static_states['H'][1], h_out_guess), (self.stages[-1].static_states['P'][1], pout_guess)])         
        
        # Initial guess vector
        x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

        res = 1
        x_in = x0_disc
        
        c = 0
            
        while res > 1e-8:

            if c > 1000:
                raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")                

            x_out = self.last_blade_row_system(x_in) 

            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            
            x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                        
            c += 1
        
        self.last_blade_row_system(x_out)
        
        return
    
    def generate_map_m_dot_N_rot(
        self,
        m_dot_grid=None,
        N_rot_grid=None,
        *,
        m_dot_range=None,  # (min, max, n)
        N_rot_range=None,  # (min, max, n)
        fixed_P_su=None,   # override supply pressure for the map (Pa)
        fixed_T_su=None,   # override supply temperature for the map (K)
        fixed_P_ex=None,   # override exhaust/static outlet pressure (Pa)
        per_point_hook=None,  # callable(self) to tweak params before solve
        max_retries=2,
        mach_limit=1.2,
        pressure_tol=0.02,
        verbose=False
        ):
        """
        Build an operation map by sweeping mass flow and speed.
    
        Parameters
        ----------
        m_dot_grid : iterable of float
            Mass-flow values [kg/s]. If None, use m_dot_range.
        N_rot_grid : iterable of float
            Rotational speeds [rpm]. If None, use N_rot_range.
        m_dot_range : (mmin, mmax, n)
            Range spec if m_dot_grid is None.
        N_rot_range : (nmin, nmax, n)
            Range spec if N_rot_grid is None.
        fixed_P_su, fixed_T_su, fixed_P_ex : float or None
            If given, overrides the current inputs for the whole map (Pa, K, Pa).
        per_point_hook : callable(self) -> None
            Called right before `solve()` for each (ṁ, N) point to adjust geometry,
            clearances, or params based on current operating point.
        max_retries : int
            Retries for a point if convergence fails (e.g., tweak damping).
        mach_limit : float
            Warn if any stage exit Mach exceeds this value.
        pressure_tol : float
            Relative tolerance on outlet pressure: warn if |P_ex_calc - P_ex_target| / P_ex_target > pressure_tol.
        verbose : bool
            Print progress.
    
        Returns
        -------
        pandas.DataFrame
            One row per (ṁ, N). Columns include:
            ['m_dot','N_rot','P_su','T_su','P_ex_target','P_ex_calc',
             'PR','W_dot','eta_is','converged','mach_warn','pressure_warn','notes']
        """
        import numpy as _np
        import pandas as _pd
        
        # Build grids if ranges given
        if m_dot_grid is None:
            if m_dot_range is None:
                # default: ±30% around current
                m0 = self.inputs.get('m_dot', 1.0)
                m_dot_grid = _np.linspace(0.7*m0, 1.3*m0, 9)
            else:
                m_dot_grid = _np.linspace(*m_dot_range)
        else:
            m_dot_grid = _np.array(list(m_dot_grid), dtype=float)
    
        if N_rot_grid is None:
            if N_rot_range is None:
                N0 = self.inputs.get('N_rot', 1000.0)
                N_rot_grid = _np.linspace(0.6*N0, 1.2*N0, 9)
            else:
                N_rot_grid = _np.linspace(*N_rot_range)
        else:
            N_rot_grid = _np.array(list(N_rot_grid), dtype=float)
    
        # Cache original boundary conditions to restore later
        _P_su0 = self.inputs.get('P_su', None)
        _T_su0 = self.inputs.get('T_su', None)
        _P_ex0 = self.inputs.get('P_ex', None)
    
        rows = []
        total_pts = len(m_dot_grid) * len(N_rot_grid)
        idx = 0
    
        for N in N_rot_grid:
            for m in m_dot_grid:
                idx += 1
                if verbose:
                    print(f"[{idx}/{total_pts}] N={N:.2f} rpm, ṁ={m:.3f} kg/s")
    
                # Set operating point
                self.set_inputs(
                    m_dot = float(m),
                    P_su  = float(fixed_P_su if fixed_P_su is not None else _P_su0),
                    T_su  = float(fixed_T_su if fixed_T_su is not None else _T_su0),
                    N_rot = float(N),
                    fluid = self.fluid,
                    P_ex  = float(fixed_P_ex if fixed_P_ex is not None else _P_ex0),
                )
    
                # Optionally let user tweak geometry/params per point
                if per_point_hook is not None:
                    try:
                        per_point_hook(self)
                    except Exception as e:
                        # Non-fatal: record the note and proceed
                        hook_note = f"per_point_hook failed: {e}"
                    else:
                        hook_note = ""
    
                # Try to solve with limited retries (e.g., adjust damping on the fly)
                converged = False
                notes = []
                for attempt in range(max_retries + 1):
                    try:
                        # Example adaptive damping: increase damping on retries
                        if attempt > 0 and 'damping' in self.params:
                            self.params['damping'] = min(0.8, self.params['damping'] * 1.5)
    
                        # Re-init the *first* stage inlet each attempt to reduce drift
                        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
    
                        self.solve()
                        converged = True
                        break
                    except Exception as e:
                        notes.append(f"attempt {attempt}: {e}")
    
                # Collect metrics (even if failed; fill NaNs)
                if converged:
                    try:
                        P_su = self.inputs['P_su']
                        T_su = self.inputs['T_su']
                        P_ex_target = self.inputs['P_ex']
                        P_ex_calc = float(self.stages[-1].static_states['P'][2])
                
                        W_dot = getattr(self, 'W_dot', _np.nan)
                        eta   = getattr(self, 'eta_is', _np.nan)
                
                        # >>> INSERT THESE THREE LINES HERE <<<
                        RP_target = P_su / P_ex_target if P_ex_target else _np.nan
                        RP_calc   = P_su / P_ex_calc   if _np.isfinite(P_ex_calc) and P_ex_calc else _np.nan
                        PR = RP_target  # keep old name for compatibility
                        # >>> END INSERT <<<
                
                        # ---- health checks (unchanged) ----
                        machs = []
                        for st in self.stages:
                            for label in ('M2_S','M3_R'):
                                if hasattr(st, label):
                                    val = getattr(st, label)
                                    if val is not None:
                                        machs.append(val)
                        mach_warn = bool(len([x for x in machs if _np.isfinite(x) and x > mach_limit]) > 0)
                
                        pressure_warn = False
                        if P_ex_target and _np.isfinite(P_ex_calc):
                            rel_err = abs(P_ex_calc - P_ex_target)/P_ex_target
                            pressure_warn = rel_err > pressure_tol
                            if pressure_warn:
                                notes.append(f"P_ex mismatch {rel_err:.2%}")
    
                    except Exception as e:
                        # If something unexpected happens while reading results
                        P_su = self.inputs.get('P_su', _np.nan)
                        T_su = self.inputs.get('T_su', _np.nan)
                        P_ex_target = self.inputs.get('P_ex', _np.nan)
                        P_ex_calc = _np.nan
                        PR = _np.nan
                        W_dot = _np.nan
                        eta = _np.nan
                        mach_warn = False
                        pressure_warn = True
                        notes.append(f"post-process error: {e}")
    
                else:
                    # Not converged
                    P_su = self.inputs.get('P_su', _np.nan)
                    T_su = self.inputs.get('T_su', _np.nan)
                    P_ex_target = self.inputs.get('P_ex', _np.nan)
                    P_ex_calc = _np.nan
                    PR = _np.nan
                    W_dot = _np.nan
                    eta = _np.nan
                    mach_warn = False
                    pressure_warn = True
    
                rows.append(dict(
                    m_dot=float(m),
                    N_rot=float(N),
                    P_su=float(P_su) if P_su is not None else _np.nan,
                    T_su=float(T_su) if T_su is not None else _np.nan,
                    P_ex_target=float(P_ex_target) if P_ex_target is not None else _np.nan,
                    P_ex_calc=float(P_ex_calc) if _np.isfinite(P_ex_calc) else _np.nan,
                    PR=float(PR) if _np.isfinite(PR) else _np.nan,
                    W_dot=float(W_dot) if _np.isfinite(W_dot) else _np.nan,
                    eta_is=float(eta) if _np.isfinite(eta) else _np.nan,
                    converged=bool(converged),
                    mach_warn=bool(mach_warn),
                    pressure_warn=bool(pressure_warn),
                    notes="; ".join([hook_note] + notes) if 'hook_note' in locals() else "; ".join(notes)
                ))
    
        # Restore original boundary conditions
        self.set_inputs(
            m_dot=self.inputs['m_dot'],
            P_su=_P_su0,
            T_su=_T_su0,
            N_rot=self.inputs['N_rot'],
            fluid=self.fluid,
            P_ex=_P_ex0
        )
    
        df = _pd.DataFrame(rows)
    
        # Useful sorted ordering
        df.sort_values(by=['N_rot','m_dot'], inplace=True, ignore_index=True)
        return df

    
    def solve(self):
        
        self.omega_rads = 2*np.pi*self.inputs['N_rot']/60
        self.u = self.omega_rads*self.params['r_m']*2
        
        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
        # self.stages.append(self.stage(self.fluid))
        
        self.computeRepeatingStages()
        
        self.computeLastStage()
        
        hin = self.stages[0].total_states['H'][1]
        hout = self.stages[-1].static_states['H'][2]
        
        self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])

        hout_s = self.AS.hmass()
        
        self.W_dot = self.inputs['m_dot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - hout_s)
        
        return 


#%%

if __name__ == "__main__":
    
    # # --- 1) Faire un snapshot picklable de la machine ---
    def _snapshot_from_machine(machine):
        base_inputs = dict(
            P_su = float(machine.inputs['P_su']),
            T_su = float(machine.inputs['T_su']),
            P_ex = float(machine.inputs['P_ex']),
            m_dot = float(machine.inputs['m_dot']),
            N_rot = float(machine.inputs['N_rot']),
            fluid = machine.fluid
        )
        base_params = dict(machine.params)  # copie légère

        stage_params = {}
        import numbers, numpy as np
        def _ok(v): return (v is None) or isinstance(v, (numbers.Number, np.floating, np.integer))
        if machine.stages:
            keys = set().union(*(vars(st).keys() for st in machine.stages))
            blacklist = {'AS','total_states','static_states','Vel_Tri_S','Vel_Tri_R',
                          'eta_is_R','eta_is_S','M1_S','M2_S','M2_R','M3_R',
                          'Y_vec_S','Y_vec_R','delta_S','delta_R','beta_g_S','beta_g_R'}
            for k in sorted(keys - blacklist):
                vals = [getattr(st, k, None) for st in machine.stages]
                if any((_ok(v) for v in vals)) and all((_ok(v) for v in vals)):
                    stage_params[k] = vals
        return base_inputs, base_params, stage_params

    def _eval_point_from_snapshot(m, N, base_inputs, base_params, stage_params):
        try:
            turb = AxialTurbineMeanLine(base_inputs['fluid'])
            turb.set_inputs(
                m_dot=float(m),
                P_su=base_inputs['P_su'],
                T_su=base_inputs['T_su'],
                N_rot=float(N),
                fluid=base_inputs['fluid'],
                P_ex=base_inputs['P_ex'],
            )
            turb.set_parameters(**base_params)
            
            if stage_params:
                turb.set_stage_parameters(**stage_params)

            # print(f"Computing for N = {float(N)} RPM and m_dot = {float(m)} kg/s \n ")

            turb.solve()

            P_ex_calc = float(turb.stages[-1].static_states['P'][2])
            RP_calc = turb.inputs['P_su'] / P_ex_calc if P_ex_calc else np.nan
            RP_target = turb.inputs['P_su'] / turb.inputs['P_ex'] if turb.inputs.get('P_ex') else np.nan

            # ✅ Ajout du print de debug
            # print(f"[DONE] m={m:.2f}, N={N:.1f}, W={turb.W_dot:.2e}, eta={turb.eta_is:.3f}, RP={RP_calc:.2f}")
            
            if turb.eta_is < 0.3:
                return dict(
                    m_dot=float(m), N_rot=float(N),
                    P_su=float(base_inputs.get('P_su', np.nan)),
                    T_su=float(base_inputs.get('T_su', np.nan)),
                    P_ex_target=float(base_inputs.get('P_ex', np.nan)),
                    P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
                    W_dot=np.nan, eta_is=np.nan, converged=False
                )
                
            else:
                return dict(
                    m_dot=float(m),
                    N_rot=float(N),
                    P_su=float(turb.inputs['P_su']),
                    T_su=float(turb.inputs['T_su']),
                    P_ex_target=float(turb.inputs['P_ex']),
                    P_ex_calc=P_ex_calc,
                    RP_target=RP_target,
                    RP_calc=RP_calc,
                    W_dot=float(getattr(turb, 'W_dot', np.nan)),
                    eta_is=float(getattr(turb, 'eta_is', np.nan)),
                    converged=True,
                    note=""
                )

        except Exception as e:
            return dict(
                m_dot=float(m), N_rot=float(N),
                P_su=float(base_inputs.get('P_su', np.nan)),
                T_su=float(base_inputs.get('T_su', np.nan)),
                P_ex_target=float(base_inputs.get('P_ex', np.nan)),
                P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
                W_dot=np.nan, eta_is=np.nan, converged=False, note=str(e)
            )

    import sys
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from tqdm import tqdm
    import pandas as pd

    def generate_map_threads(machine, m_grid, N_grid, max_workers=6, desc="Operation map"):
        base_inputs, base_params, stage_params = _snapshot_from_machine(machine)
        tasks = [(m, N) for N in N_grid for m in m_grid]
        results = []

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = [ex.submit(_eval_point_from_snapshot, m, N, base_inputs, base_params, stage_params)
                        for (m, N) in tasks]
            with tqdm(total=len(futures),
                      desc=desc,
                      unit="pt",
                      dynamic_ncols=True,
                      miniters=1,
                      mininterval=0,
                      ascii=True,
                      file=sys.stdout) as bar:
                for fut in as_completed(futures):
                    results.append(fut.result())
                    bar.update(1)
                    bar.refresh()

        return pd.DataFrame(results).sort_values(['N_rot','m_dot'], ignore_index=True)
    
#%%
    
    case_study = "TCO2_ORC"
    
    if case_study == "Salah_Case":
        Turb_OD = AxialTurbineMeanLine('CO2')
        
        Turb_OD.set_inputs(
              m_dot = 655.18,
              P_su = 25000000.0,
              T_su = 923,
              N_rot = 1506.9946780513785, # 996.4330963327212,
              fluid = 'CO2', 
              P_ex = 100*1e5
            )
        
        Turb_OD.set_parameters(
            r_m = 0.261423771889,
            nStages = 12,
            mdot_rated = 655.18,
            DP_rated = 2.5,
            damping = 0.5,
            delta_tip = 0.0004,
            N_lw = 0,
            D_lw = 0,
            e_blade = 2e-06
            )
        
        Turb_OD.set_stage_parameters(
            h_blade_S = [0.05893535333, 0.06254061127, 0.06644686421, 0.07068270152, 0.07527979324, 0.08027326565, 0.08570212832, 0.09160976073, 0.09804446776, 0.105060115, 0.1127168568, 0.1210819728, 0.1254195395],
            chord_S = [0.008645525688, 0.009066896935, 0.009520522543, 0.01000898242, 0.01053510451, 0.01110199118, 0.0117130487, 0.0123720203, 0.0130830232, 0.01385059029, 0.01467971707, 0.01557591465, 0.01603830059],
            xhi_S1 = [-0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816],
            xhi_S2 = [1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265],
            pitch_S = [0.006966193979, 0.007305716867, 0.007671228935, 0.008064808966, 0.008488735591, 0.008945508564, 0.009437872521, 0.009968843587, 0.01054173924, 0.0111602119, 0.01182828672, 0.01255040431, 0.01292297508],
            o_S = [0.001744270233, 0.001829283609, 0.001920804434, 0.002019353216, 0.002125500504, 0.002239872211, 0.00236315557, 0.002496105791, 0.002639553539, 0.002794413345, 0.00296169307, 0.003142504605, 0.003235792862],
            t_TE_S = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
            t_blade_S = [0.003025933991, 0.003173413927, 0.00333218289, 0.003503143847, 0.003687286579, 0.003885696912, 0.004099567045, 0.004330207105, 0.00457905812, 0.004847706601, 0.005137900974, 0.005451570129, 0.005613405206],
            n_blade_S = [236, 225, 214, 204, 194, 184, 174, 165, 156, 147, 139, 131, 127],
            R_c_S = [0.01352981312, 0.01418923794, 0.01489913922, 0.01566355437, 0.01648690899, 0.01737405813, 0.01833033244, 0.01936158987, 0.02047427367, 0.02167547758, 0.02297301932, 0.02437552349, 0.02509913423],
        
            h_blade_R = [0.06085456338, 0.06461845758, 0.06869830659, 0.07312430505, 0.07792992342, 0.08315231023, 0.08883274986, 0.09501718485, 0.1017568126, 0.1091087686, 0.1171369103, 0.1259127187, None],
            chord_R = [0.008870626618, 0.009309055087, 0.009781089003, 0.01028943789, 0.01083707309, 0.01142725573, 0.01206356808, 0.01274994866, 0.01349073167, 0.01429069144, 0.01515509243, 0.01608974592, None],
            xhi_R1 = [0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, None],
            xhi_R2 = [-1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, None],
            pitch_R = [0.007763880051, 0.008147607851, 0.008560748305, 0.009005672883, 0.009484982197, 0.01000153051, 0.01055845315, 0.01115919723, 0.01180755622, 0.01250770875, 0.01326426248, 0.01408230362, None],
            o_R = [0.00229560867, 0.002409068546, 0.002531225097, 0.002662779514, 0.002804500743, 0.002957232725, 0.003121902507, 0.003299529328, 0.003491234828, 0.00369825454, 0.003921950849, 0.004163827632, None],
            t_TE_R = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
            t_blade_R = [0.003104719316, 0.00325816928, 0.003423381151, 0.003601303262, 0.00379297558, 0.003999539504, 0.004222248827, 0.00446248203, 0.004721756085, 0.005001742002, 0.005304282352, 0.00563141107, None],
            n_blade_R = [212, 202, 192, 182, 173, 164, 156, 147, 139, 131, 124, 117, None],
            R_c_R = [0.01388208476, 0.01456820328, 0.01530691263, 0.01610245309, 0.01695947464, 0.0178830808, 0.01887887765, 0.01995302875, 0.02111231693, 0.02236421375, 0.02371695786, 0.02517964358, None],
            )
        
        df_map = generate_map_threads(
            Turb_OD,
            m_grid=np.linspace(0.5*Turb_OD.inputs['m_dot'], 1.5*Turb_OD.inputs['m_dot'], 10),
            N_grid=np.linspace(0.5*Turb_OD.inputs['N_rot'], 1.5*Turb_OD.inputs['N_rot'], 10),
            max_workers=6  # 8 threads logiques - 2 pour laisser de la marge
        )
        
        fig, ax = map_plot(
            df_map, levels=24, dpi=300,
            use_grid=True, nx=500, ny=500, smooth_sigma=0.8
        )
        
        plt.show()
        
        _ = plot_power_eta_vs_mdot(df_map, speeds=None, max_lines=5)  # auto-picks up to 5 speeds
        plt.show()
        
    elif case_study == "TCO2_ORC":
        Turb_OD = AxialTurbineMeanLine('CO2')
        
        Turb_OD.set_inputs(
            m_dot = 10*100, # kg/s
            P_su = 140*1e5, # Pa
            T_su = 273.15 + 121, # K
            fluid = 'CO2', 
            N_rot = 542.2770417507232, # RPM
            P_ex = 39.8*1e5, # Pa
            )
        
        Turb_OD.set_inputs(
              m_dot = 10*100,
              P_su = 140*1e5,
              T_su = 273.15 + 121,
              N_rot = 542.277, 
              fluid = 'CO2', 
              P_ex = 39.8*1e5
            )
        
        Turb_OD.set_parameters(
            r_m = 0.429443878686,
            nStages = 17,
            mdot_rated = 1000,
            DP_rated = 3.5175879397,
            damping = 0.2,
            delta_tip = 0.0004,
            N_lw = 0,
            D_lw = 0,
            e_blade = 2e-06
            )
        
        Turb_OD.set_stage_parameters(
            h_blade_S = [0.06186052046, 0.06497930968, 0.06837023781, 0.07205830785, 0.0760708064, 0.08043755356, 0.08519118879, 0.09036750015, 0.09600581217, 0.1021494566, 0.1088463398, 0.1161495644, 0.1241179972, 0.1328167815, 0.1423181134, 0.1527026414, 0.1640611996, 0.1691870839],
            chord_S = [0.01297062444, 0.01322718946, 0.01351015762, 0.01382120246, 0.01416206854, 0.01453458074, 0.01494065578, 0.01538231639, 0.01586170955, 0.01638113189, 0.01694306329, 0.01755019954, 0.01820546352, 0.01891199241, 0.01967315282, 0.02049264073, 0.0213746114, 0.02176946983],
            xhi_S1 = [-0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107, -0.3044754107],
            xhi_S2 = [1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293, 1.184134293],
            pitch_S = [0.01392769756, 0.01420319394, 0.01450704168, 0.01484103782, 0.01520705564, 0.01560705467, 0.01604309307, 0.01651734283, 0.01703210932, 0.01758985866, 0.01819325372, 0.01884518919, 0.01954880361, 0.02030746568, 0.0211247904, 0.02200474647, 0.02295179577, 0.02337578991],
            o_S = [0.006072654584, 0.006192774535, 0.006325256041, 0.006470882638, 0.006630471097, 0.006804875799, 0.00699499413, 0.007201773099, 0.007426217887, 0.007669403745, 0.007932491722, 0.008216743937, 0.008523528841, 0.008854315224, 0.009210679273, 0.009594351395, 0.01000727702, 0.0101921439],
            t_TE_S = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
            t_blade_S = [0.004539718553, 0.004629516312, 0.004728555166, 0.004837420861, 0.004956723989, 0.005087103258, 0.005229229523, 0.005383810737, 0.005551598342, 0.00573339616, 0.005930072153, 0.00614256984, 0.006371912231, 0.006619197345, 0.006885603488, 0.007172424256, 0.007481113989, 0.007619314439],
            n_blade_S = [194, 190, 186, 182, 177, 173, 168, 163, 158, 153, 148, 143, 138, 133, 128, 123, 118, 115],
            R_c_S = [0.01778930871, 0.01814118957, 0.01852928251, 0.01895588285, 0.01942338324, 0.01993428651, 0.02049122147, 0.02109696232, 0.02175445363, 0.02246684526, 0.0232375384, 0.02407023032, 0.02496892978, 0.02593793946, 0.02698187667, 0.02810581048, 0.02931543987, 0.02985699116],

            h_blade_R = [0.06318950062, 0.06642241346, 0.06993785474, 0.07376176811, 0.07792248332, 0.08245098083, 0.08738119491, 0.09275036538, 0.09859945767, 0.1049736752, 0.1119230598, 0.1195031032, 0.1277752753, 0.1368075922, 0.1466756406, 0.1574642019, 0.1692696806, None],
            chord_R = [0.01308030048, 0.01334806021, 0.0136429212, 0.01396658553, 0.01432083011, 0.01470751692, 0.01512860555, 0.0155861688, 0.01608241361, 0.01661971019, 0.01720062716, 0.01782795751, 0.01850471708, 0.0192341344, 0.02001969956, 0.02086529172, 0.02177535744, None],
            xhi_R1 = [0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, 0.5507188137, None],
            xhi_R2 = [-1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, -1.173070172, None],
            pitch_R = [0.0114033135, 0.01163674454, 0.01189380226, 0.01217597054, 0.01248479846, 0.01282190929, 0.01318901138, 0.01358791179, 0.01402053451, 0.01448894587, 0.01499538517, 0.01554228734, 0.01613228155, 0.01676818241, 0.01745303256, 0.01819021383, 0.0189836027, None],
            o_R = [0.003445506951, 0.003516038054, 0.00359370795, 0.00367896499, 0.003772277232, 0.003874135144, 0.003985054905, 0.004105582516, 0.00423629931, 0.004377829628, 0.004530850077, 0.004696096365, 0.004874362896, 0.005066500106, 0.005273427326, 0.005496166374, 0.005735888526, None],
            t_TE_R = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
            t_blade_R = [0.004578105168, 0.004671821074, 0.004775022419, 0.004888304934, 0.005012290538, 0.005147630922, 0.005295011943, 0.005455159081, 0.005628844764, 0.005816898566, 0.006020219505, 0.006239785129, 0.006476650977, 0.006731947039, 0.007006894845, 0.007302852102, 0.007560651135, None],
            n_blade_R = [237, 232, 227, 222, 216, 210, 205, 199, 192, 186, 180, 174, 167, 161, 155, 148, 142, None],
            R_c_R = [0.01793973023, 0.01830696471, 0.01871136876, 0.01915527682, 0.01964112592, 0.02017146978, 0.020748996, 0.02137654743, 0.02205715091, 0.02279405719, 0.02359078917, 0.02445117746, 0.02537935828, 0.0263797596, 0.02745716811, 0.02861690411, 0.02986506607, None],
            )
        
        df_map = generate_map_threads(
            Turb_OD,
            m_grid=np.linspace(0.5*Turb_OD.inputs['m_dot'], 1.5*Turb_OD.inputs['m_dot'], 10),
            N_grid=np.linspace(0.5*Turb_OD.inputs['N_rot'], 1.5*Turb_OD.inputs['N_rot'], 10),
            max_workers=6  # 8 threads logiques - 2 pour laisser de la marge
        )
        
        fig, ax = map_plot(
            df_map, levels=24, dpi=300,
            use_grid=True, nx=500, ny=500, smooth_sigma=0.8
        )
        
        plt.show()
        
        _ = plot_power_eta_vs_mdot(df_map, speeds=None, max_lines=5)  # auto-picks up to 5 speeds
        plt.show()