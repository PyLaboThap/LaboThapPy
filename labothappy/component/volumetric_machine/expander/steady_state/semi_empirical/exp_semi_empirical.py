# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 2024

@author: Elise Neven (elise.neven@uliege.be)
"""

"External modules"
from CoolProp.CoolProp import AbstractState
import CoolProp.CoolProp as CoolProp
from scipy.optimize import fsolve
import numpy as np
import time

"Internal modules"
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

class ExpanderSE(BaseComponent):
    """
        Component: Volumetric expander

        Model: The model is based on the thesis of V. Lemort (2008) and is a semi-empirical model.

        **Descritpion**:

            This model is used to simulate the performance of a volumetric expander. 
            The parameters of the model need to be calibrated with experimental datas to represent the real behavior of the expander.

        **Assumptions**:

            - Steady-state operation.

        **Connectors**:

            su (MassConnector): Mass connector for the suction side.

            ex (MassConnector): Mass connector for the exhaust side.

            W_exp (WorkConnector): Work connector.

            Q_amb (HeatConnector): Heat connector for the ambient heat transfer.

        **Parameters**:

            AU_amb: Heat transfer coefficient for the ambient heat transfer. [W/K]

            AU_su_n: Nominal heat transfer coefficient for the suction side heat transfer. [W/K]

            AU_ex_n: Nominal heat transfer coefficient for the exhaust side heat transfer. [W/K]

            d_su1: Pressure drop diameter. [m]

            m_dot_n: Nominal mass flow rate. [kg/s]

            A_leak: Leakage area. [m^2]

            W_dot_loss_0: Constant loss in the compressor. [W]

            alpha: Proportionality rate for mechanical losses. [-]

            C_loss: Torque losses. [N.m]

            rv_in: Inlet volume ratio. [-]

            V_s: Swept volume. [m^3]

        **Inputs**:

            su_p: Suction side pressure. [Pa]

            su_T: Suction side temperature. [K]

            ex_p: Exhaust side pressure. [Pa]

            su_fluid: Suction side fluid. [-]

            N_rot: Rotational speed. [rpm]

            T_amb: Ambient temperature. [K]

        **Ouputs**:

            eta_is: Isentropic efficiency. [-]

            ex_h: Exhaust side specific enthalpy. [J/kg]

            ex_T: Exhaust side temperature. [K]

            W_dot_exp: Compressor power. [W]

            m_dot: Mass flow rate. [kg/s]
            
            epsilon_v: Volumetric efficiency. [-]
    """
    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_exp = WorkConnector()
        self.Q_amb = HeatConnector()

    def get_required_inputs(self):
        if self.params['mode'] == 'N_rot':
            return ['P_su', 'T_su', 'P_ex', 'N_rot', 'T_amb', 'fluid']
        elif self.params['mode'] == 'm_dot':
            return ['P_su', 'T_su', 'P_ex', 'm_dot', 'T_amb', 'fluid']

    # def get_required_parameters(self):
    #     return [
    #         'AU_amb', 'AU_su_n', 'AU_ex_n', 'd_su1', 'm_dot_n', 
    #         'A_leak', 'W_dot_loss_0', 'alpha', 'C_loss', 'rv_in', 'V_s',
    #         'mode'
    #     ]
    
    # def get_required_parameters(self):
    #     required_params = [
    #         'AU_amb', 'AU_su_n', 'AU_ex_n', 'd_su1', 'm_dot_n', 
    #         'A_leak', 'W_dot_loss_0', 'alpha', 'C_loss', 'rv_in', 'V_s',
    #         'mode'
    #     ]

    #     # Set default values for missing parameters
    #     for param in required_params:
    #         self.params.setdefault(param, 0)  # Sets missing values to 0
        
    #     return required_params
    
    def get_required_parameters(self):
        default_values = {
            'AU_amb': 0, #  Heat transfer to the ambient neglected
            'AU_su_n': 0, # Heat transfer to the suction side neglected
            'AU_ex_n': 0, # Heat transfer to the exhaust side neglected
            'd_su1': 0, # Pressure drop at the suction side neglected
            'm_dot_n': 1,  # Nominal mass flow rate
            'A_leak': 1e-12, # Very small leakage area to neglect the leakage effect 
            'W_dot_loss_0': 0, # Idle losses neglected
            'alpha': 0,  # Proportionality rate for mechanical losses neglected
            'C_loss': 0, # Torque losses neglected
            'rv_in': 2.5,  # Default volume ratio
            'V_s': 0.001,  # Default swept volume
            'mode': 'N_rot'  # Default mode
        }

        # Ensure all required parameters are set, assigning default values if missing
        for param, default in default_values.items():
            self.params.setdefault(param, default)
        
        return list(default_values.keys())  # Return the list of required parameters


    def System(self, x):
        if self.params['mode'] == 'N_rot':
            self.m_dot, self.T_w = x
            self.N_rot = self.inputs['N_rot']
        elif self.params['mode'] == 'm_dot':
            self.m_dot = self.inputs['m_dot']
            self.T_w = x[0]
            
        #Boundary on the mass flow rate
        self.m_dot = max(self.m_dot, 1e-5)
        print('A_leak', self.params['A_leak'])
        "1. Supply conditions: su"
        T_su = self.su.T
        P_su = self.su.p
        h_su = self.su.h
        s_su = self.su.s
        rho_su = self.su.D
        P_ex = self.ex.p
        self.P_ex = P_ex

        self.AS.update(CoolProp.PSmass_INPUTS, P_ex, s_su)
        h_ex_is = self.AS.hmass()
        self.AS.update(CoolProp.PT_INPUTS, 4e6, 500)
        h_max = self.AS.hmass()
        self.AS.update(CoolProp.PQ_INPUTS, P_su, 1)  
        T_sat_su = self.AS.T()
        if T_su<T_sat_su:
            print('----Warning the expander inlet stream is not in vapor phase---')
        
        "2. Supply pressure drop: su->su1"
        h_su1 = h_su #Isenthalpic valve
        if self.params['d_su1'] == 0:
            P_su1 = P_su
            T_su1 = T_su
        else:
            # Assumption that the density doesn't change too much
            A_su = np.pi*(self.params['d_su1']/2)**2
            V_dot_su = self.m_dot/rho_su
            C_su = V_dot_su/A_su
            h_su_thr1 = h_su-(C_su**2)/2
            h_su_thr = max(h_su_thr1, h_ex_is)
            self.AS.update(CoolProp.HmassSmass_INPUTS, h_su_thr, s_su)
            P_su_thr = self.AS.p()
            P_su1 = max(P_su_thr, P_ex+1)
            self.DP_su = P_su-P_su1
            self.AS.update(CoolProp.HmassP_INPUTS, h_su1, P_su1)
            T_su1 = self.AS.T()

        "3. Cooling at the entrance: su1->su2"  
        try:
            self.AS.update(CoolProp.HmassP_INPUTS, h_su1, P_su1)
            cp_su1 = self.AS.cpmass()
        except:
            self.AS.update(CoolProp.PQ_INPUTS, P_su1, 0)
            cp_su1 = self.AS.cpmass()
        
        if self.params['AU_su_n'] == 0:
            h_su2 = h_su1
            Q_dot_su = 0
        else:          
            AU_su = self.params['AU_su_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
            C_dot_su = self.m_dot*cp_su1
            NTU_su = AU_su/C_dot_su
            epsilon_su1 = max(0,1-np.exp(-NTU_su))
            Q_dot_su = max(0, epsilon_su1*self.m_dot*cp_su1*(T_su1-self.T_w))
            
            self.AS.update(CoolProp.PQ_INPUTS, P_su1, 0.1)
            h_su2 = min(h_max, max(max(h_ex_is, self.AS.hmass()), h_su1 - Q_dot_su/self.m_dot))
        
        P_su2 = P_su1 #No pressure drop just heat transfer
        self.AS.update(CoolProp.HmassP_INPUTS, h_su2, P_su2)
        rho_su2 = self.AS.rhomass()
        s_su2 = self.AS.smass()
        
        "4. Leakage"
        try:
            self.AS.update(CoolProp.HmassP_INPUTS, h_su1, P_su1)
            cv_su1 = self.AS.cvmass()
        except:
            self.AS.update(CoolProp.PQ_INPUTS, P_su1, 0)
            cv_su1 = self.AS.cvmass()
        gamma = max(1e-2, cp_su1/cv_su1)
        P_crit = P_su2*(2/(gamma+1))**(gamma/(gamma-1))
        P_thr_leak = max(P_ex, P_crit)
        self.AS.update(CoolProp.PSmass_INPUTS, P_thr_leak, s_su2)
        rho_thr_leak = self.AS.rhomass()
        h_thr_leak = self.AS.hmass()
        C_thr_leak = np.sqrt(2*(h_su2-h_thr_leak))
        m_dot_leak = self.params['A_leak']*C_thr_leak*rho_thr_leak

        if self.params['mode'] == 'N_rot':
            m_dot_in = (self.N_rot/60)*self.params['V_s']*rho_su2
            m_dot_leak_bis = self.m_dot-m_dot_in
        elif self.params['mode'] == 'm_dot':
            m_dot_in = self.m_dot-m_dot_leak
            self.N_rot = (m_dot_in/(self.params['V_s']*rho_su2))*60
        
        "5. Internal expansion"
        "Isentropic expansion to the internal pressure: su2->in"
        rho_in = rho_su2/self.params['rv_in']
        #Trick to not have problems with CoolProp
        try:
            self.AS.update(CoolProp.DmassSmass_INPUTS, rho_in, s_su2)
            P_in = self.AS.p()
        except:
            delta = 0.0001
            self.AS.update(CoolProp.DmassSmass_INPUTS, rho_in*(1+delta), s_su2)
            P_in1 = self.AS.p()
            self.AS.update(CoolProp.DmassSmass_INPUTS, rho_in*(1-delta), s_su2)
            P_in2 = self.AS.p()
            P_in = 0.5*P_in1+0.5*P_in2
            
        self.AS.update(CoolProp.DmassP_INPUTS, rho_in, P_in)
        h_in = self.AS.hmass()
        w_in_s = h_su2-h_in
        "Expansion at constant volume: in->ex2"
        w_in_v = (P_in-P_ex)/rho_in
        h_ex2 = h_in-w_in_v
        "Total work"
        W_dot_in = m_dot_in*(w_in_s+w_in_v)

        "6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
        h_ex1 = max(min((m_dot_in*h_ex2 + m_dot_leak*h_su2)/self.m_dot, h_su2), h_ex2)
        P_ex1 = P_ex
        self.AS.update(CoolProp.HmassP_INPUTS, h_ex1, P_ex1)
        T_ex1 = self.AS.T()
        try:
            self.AS.update(CoolProp.HmassP_INPUTS, h_ex1, P_ex1)
            cp_ex2 = self.AS.cpmass()
        except:
            self.AS.update(CoolProp.PQ_INPUTS, P_ex1, 0)
            cp_ex2 = self.AS.cpmass()
        
        "7. Exhaust side heat transfer: ex1->ex"
        if self.params['AU_ex_n'] == 0:
            self.h_ex = h_ex1
            Q_dot_ex = 0
        else:
            AU_ex = self.params['AU_ex_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
            C_dot_ex = self.m_dot*cp_ex2
            NTU_ex=AU_ex/C_dot_ex
            epsilon_ex = max(0, 1-np.exp(-NTU_ex))
            Q_dot_ex = max(0, epsilon_ex*C_dot_ex*(self.T_w-T_ex1))
            self.h_ex = h_ex1+Q_dot_ex/self.m_dot

        "8. Energy balance"
        self.Q_dot_amb = self.params['AU_amb']*(self.T_w-self.inputs['T_amb'])
        W_dot_loss = self.params['alpha']*W_dot_in + self.params['W_dot_loss_0'] + self.params['C_loss']*(self.N_rot/60)*2*np.pi
        self.W_dot_exp = W_dot_in - W_dot_loss

        "9. Performances"
        W_dot_s = self.m_dot*(h_su-h_ex_is)
        self.epsilon_is = self.W_dot_exp/W_dot_s
        self.m_dot_th = (self.N_rot/60)*self.params['V_s']*rho_su
        self.epsilon_v = self.m_dot/self.m_dot_th
        
        "10. Residuals"
        self.res_E = abs((Q_dot_su + W_dot_loss - Q_dot_ex - self.Q_dot_amb)/(Q_dot_su + W_dot_loss))
        if self.params['mode'] == 'N_rot':
            self.res_m_leak = abs((m_dot_leak_bis-m_dot_leak)/m_dot_leak)
            self.res = [self.res_E, self.res_m_leak]
        elif self.params['mode'] == 'm_dot':
            self.res = [self.res_E]
        print('residuals', self.res)
        return self.res

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        fluid = self.su.fluid  # Extract fluid name
        self.AS = AbstractState("HEOS", fluid)  # Create a reusable state object

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("ExpanderSE could not be solved. It is not calculable and/or not parametrized")
            self.print_setup()
            return
        
        ff_guess = [0.7, 1.2, 0.8, 1.3, 0.4, 1.7]
        x_T_guess = [0.7, 0.95, 0.8, 0.9]
        stop = 0
        j = 0
        
        if self.params['mode'] == 'N_rot':
            try:  
                while not stop and j < len(x_T_guess):
                    k = 0
                    while not stop and k < len(ff_guess):
                        # Loop to permit multiple attempts to solve the implicit calculation 
                        self.AS.update(CoolProp.PT_INPUTS, self.inputs['P_su'], self.inputs['T_su'])
                        m_dot_guess = ff_guess[k] * self.params['V_s'] * self.inputs['N_rot'] / 60 * self.AS.rhomass()
                        T_w_guess = x_T_guess[j] * self.inputs['T_su'] + (1 - x_T_guess[j]) * self.inputs['T_amb']
                        #---------------------------------------------------------------------
                        args = ()
                        x = [m_dot_guess, T_w_guess]
                        #--------------------------------------------------------------------------
                        try: # !!!
                            fsolve(self.System, x, args=args)
                            res_norm = np.linalg.norm(self.res)
                        except: # !!!
                            res_norm = 1
                        if res_norm < 1e-4:
                            stop = 1
                        k += 1
                    j += 1
            except Exception as e: # !!!
                print(f"ExpanderSE could not be solved. Error: {e}")
                self.solved = False
        
        if self.params['mode'] == 'm_dot':
            try:  
                while not stop and j < len(x_T_guess):
                    # Loop to permit multiple attempts to solve the implicit calculation 
                    T_w_guess = x_T_guess[j] * self.inputs['T_su'] + (1 - x_T_guess[j]) * self.inputs['T_amb']
                    #---------------------------------------------------------------------
                    args = ()
                    x = [T_w_guess]
                    #--------------------------------------------------------------------------
                    try:
                        fsolve(self.System, x, args=args)
                        res_norm = np.linalg.norm(self.res)
                    except:
                        res_norm = 1
                    if res_norm < 1e-4:
                        stop = 1
                    j += 1

            except Exception as e:
                print(f"ExpanderSE could not be solved. Error: {e}")
                self.solved = False
    
        self.convergence = stop

        if self.convergence:
            self.update_connectors()
            self.solved = True

    def update_connectors(self):
        """Update the connectors with the calculated values."""
        self.su.set_m_dot(self.m_dot)
        
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.m_dot)
        self.ex.set_h(self.h_ex)
        self.ex.set_p(self.P_ex)

        self.W_exp.set_W_dot(self.W_dot_exp)
        self.W_exp.set_N(self.N_rot)
        self.Q_amb.set_Q_dot(self.Q_dot_amb)
        self.Q_amb.set_T_hot(self.T_w)

    def print_results(self):
        print("=== Expander Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_exp: {self.W_exp.W_dot} [W]")
        print(f"  - epsilon_is: {self.epsilon_is} [-]")
        print(f"  - m_dot: {self.m_dot} [kg/s]")
        print(f"  - epsilon_v: {self.epsilon_v} [-]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Expander Results ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")
        print("Work connector:")
        print(f"  - W_dot_exp: {self.W_exp.W_dot} [W]")
        print("=========================")
        print("Heat connector:")
        print(f"  - Q_dot_amb: {self.Q_amb.Q_dot} [W]")
        print(f"  - T_hot: {self.Q_amb.T_hot} [K]")
        print(f"  - T_cold: {self.Q_amb.T_cold} [K]")
        print("=========================")