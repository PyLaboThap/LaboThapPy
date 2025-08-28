from connector.mass_connector import MassConnector
from correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector

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
        while len(self.stages) < n_stages-1:
            # ⬇️ replace Stage() with your actual Stage class constructor
            self.stages.append(self.stage(self.fluid))
    
        # assign parameters
        for key, array in parameters.items():
            for i in range(len(array)-1):
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
        
        print(w3)
        print(w3_new)
        print("---------------")
        print(hout)
        print(pout_calc)
        print(stage.static_states['D'][3])
        print("---------------")
        
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
        stage.h_blade_S = stage.A_flow_S/(2*np.pi*self.r_m)

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
            self.compute_deviation_rotor(stage)
                        
            # print(f'Y_R : {stage.Y_vec_R}')

        return
            
    def computeRepeatingStages(self):
                
        self.nStages = self.params['nStages']
        
        for i in range(int(self.nStages)):
        
            print(i)
            
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
    
    def solve(self):
        
        self.omega_rads = 2*np.pi*self.inputs['N_rot']/60
        self.u = self.omega_rads*self.params['r_m']*2
        
        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
        
        self.stages.append(self.stage(self.fluid))
        
        self.computeRepeatingStages()
        
        
        
        return 

Turb_OD = AxialTurbineMeanLine('CO2')

Turb_OD.set_inputs(
     m_dot = 655.18,
     P_su = 25000000.0,
     T_su = 923,
     N_rot = 1859.8609382646837,
     fluid = 'CO2', 
     P_ex = 100*1e5
    )

Turb_OD.set_parameters(
    r_m = 0.244343625629,
    nStages = 14,
    mdot_rated = 655.18,
    DP_rated = 2.5,
    damping = 0.2,
    delta_tip = 0.0004,
    N_lw = 0,
    D_lw = 0,
    e_blade = 2e-06
    )

Turb_OD.set_stage_parameters(
    h_blade_S = [0.05432719733, 0.05722765621, 0.06033982489, 0.06368121979, 0.06727099697, 0.07113012869, 0.07528159678, 0.07975060915, 0.08456484736, 0.08975473956, 0.09535377108, 0.1013988309, 0.1079306086, 0.1149940376, 0.1182025127],
    chord_S = [0.008662625433, 0.0090311472, 0.009424359133, 0.009843981111, 0.01029187648, 0.01077006579, 0.01128074127, 0.01182628288, 0.01240927658, 0.01303253391, 0.01369911408, 0.0144123481, 0.01517586629, 0.01599362819, 0.01636401415],
    xhi_S1 = [-0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, -0.05819637398, None],
    xhi_S2 = [1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, 1.055766315, None],
    pitch_S = [0.008163086796, 0.008510357401, 0.008880894389, 0.009276318462, 0.009698385516, 0.01014899957, 0.01063022646, 0.01114430889, 0.0116936837, 0.01228100029, 0.01290914148, 0.01358124617, 0.01430073534, 0.01507134022, None],
    o_S = [0.004015711895, 0.004186546621, 0.004368826906, 0.004563350031, 0.004770979784, 0.004992652819, 0.005229385393, 0.00548228077, 0.005752537721, 0.006041459582, 0.00635046451, 0.006681096644, 0.007035038884, 0.007414126753, None],
    t_TE_S = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
    t_blade_S = [0.003031918902, 0.00316090152, 0.003298525696, 0.003445393389, 0.003602156769, 0.003769523026, 0.003948259445, 0.004139199008, 0.004343246803, 0.004561386867, 0.004794689929, 0.005044321836, 0.0053115532, 0.005597769868, None],
    n_blade_S = [188, 180, 173, 166, 158, 151, 144, 138, 131, 125, 119, 113, 107, 102, None],
    R_c_S = [0.009385810203, 0.009785097392, 0.01021113596, 0.01066578938, 0.01115107655, 0.01166918668, 0.01222249505, 0.01281358029, 0.01344524424, 0.01412053316, 0.01484276167, 0.01561553884, 0.01644279805, 0.01732882945, None],

    h_blade_R = [0.05555881627, 0.05854748153, 0.06175507837, 0.06519979435, 0.0689015314, 0.07288208398, 0.07716534482, 0.08177753586, 0.08674746124, 0.09210679875, 0.09789042501, 0.1041367804, 0.1108882875, 0.1181918185, None],
    chord_R = [0.008820158445, 0.009199013019, 0.009603267441, 0.0100347005, 0.01049524049, 0.01098697843, 0.01151218371, 0.0120733212, 0.01267306943, 0.01331434164, 0.01400030867, 0.01473442425, 0.0155204538, 0.01636250586, None],
    xhi_R1 = [-1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, -1.024956929, None],
    xhi_R2 = [0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, 0.2915700805, None],
    pitch_R = [0.006677156358, 0.006963962003, 0.007269996186, 0.007596605514, 0.007945249757, 0.008317511908, 0.008715109953, 0.009139909891, 0.009593939458, 0.01007940407, 0.01059870417, 0.01115445433, 0.01174950512, 0.01238696682, None],
    o_R = [0.002353831751, 0.002454936502, 0.002562819699, 0.002677956048, 0.002800860148, 0.002932090035, 0.003072251332, 0.003222001844, 0.003382056387, 0.003553192415, 0.003736256132, 0.003932169229, 0.004141936584, 0.004366654641, None],
    t_TE_R = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
    t_blade_R = [0.003087055456, 0.003219654557, 0.003361143604, 0.003512145176, 0.00367333417, 0.003845442451, 0.004029264298, 0.00422566242, 0.004435574301, 0.004660019572, 0.004900108035, 0.005157048488, 0.00543215883, 0.005726877051, None],
    n_blade_R = [230, 220, 211, 202, 193, 185, 176, 168, 160, 152, 145, 138, 131, 124, None],
    R_c_R = [0.009556494595, 0.009966977208, 0.01040498013, 0.01087243066, 0.01137141805, 0.01190420792, 0.01247325908, 0.01308124219, 0.0137310594, 0.01442586714, 0.01516910098, 0.01596450297, 0.01681615288, 0.01772850224, None],
    )


Turb_OD.solve()

