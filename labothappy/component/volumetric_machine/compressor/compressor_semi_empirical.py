# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:38:08 2023

@author: Elise
"""

"External modules"
from CoolProp.CoolProp import AbstractState
import CoolProp.CoolProp as CoolProp
from scipy.optimize import fsolve
import numpy as np

"Internal modules"
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector


class CompressorSE(BaseComponent):
    """
        Component: Volumetric compressor

        Model: The model is based on the thesis of V. Lemort (2008) and is a semi-empirical model.

        **Descritpion**:

            This model is used to simulate the performance of a volumetric compressor. 
            The parameters of the model need to be calibrated with experimental datas to represent the real behavior of the compressor.

        **Assumptions**:

            - Steady-state operation.

        **Connectors**:

            su (MassConnector): Mass connector for the suction side.

            ex (MassConnector): Mass connector for the exhaust side.

            W_mec (WorkConnector): Work connector.

            Q_amb (HeatConnector): Heat connector for the ambient heat transfer.

        **Parameters**:

            AU_amb: Heat transfer coefficient for the ambient heat transfer. [W/K]

            AU_su_n: Nominal heat transfer coefficient for the suction side heat transfer. [W/K]

            AU_ex_n: Nominal heat transfer coefficient for the exhaust side heat transfer. [W/K]

            d_ex: Pressure drop diameter. [m]

            m_dot_n: Nominal mass flow rate. [kg/s]

            A_leak: Leakage area. [m^2]

            W_dot_loss_0: Constant loss in the compressor. [W]

            alpha: Loss coefficient. [-]

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

            W_dot_cp: Compressor power. [W]

            m_dot: Mass flow rate. [kg/s]
            
            epsilon_v: Volumetric efficiency. [-]
    """
    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Suction side mass connector
        self.ex = MassConnector() # Exhaust side mass connector
        self.W_mec = WorkConnector() # Work connector of the compressor
        self.Q_amb = HeatConnector() # Heat connector to the ambient

    def get_required_inputs(self):
        # Define the required inputs for the component
        # If the mode is 'N_rot', the rotational speed of the compressor is required while the mass flow rate is calculated
        if self.params['mode'] == 'N_rot':
            return ['P_su', 'h_su', 'P_ex', 'N_rot', 'T_amb', 'fluid']
        # If the mode is 'm_dot', the mass flow rate is required while the rotational speed of the compressor is calculated
        elif self.params['mode'] == 'm_dot':
            return ['P_su', 'h_su', 'P_ex', 'm_dot', 'T_amb', 'fluid']

    def get_required_parameters(self):
        default_values = {
            'AU_amb': 0, #  Heat transfer to the ambient neglected
            'AU_su_n': 0, # Heat transfer to the suction side neglected
            'AU_ex_n': 0, # Heat transfer to the exhaust side neglected
            'd_ex': 0, # Pressure drop at the suction side neglected
            'm_dot_n': 1,  # Nominal mass flow rate
            'A_leak': 1e-12, # Very small leakage area to neglect the leakage effect 
            'W_dot_loss_0': 0, # Idle losses neglected
            'alpha': 0,  # Proportionality rate for mechanical losses neglected
            'C_loss': 0, # Torque losses neglected
            'rv_in': 2.5,  # Default volume ratio
            'V_s': 0.001,  # Default swept volume
            'mode': 'N_rot'  # Default mode
        }

        # Ensure all required parameters are set, assigning default values if missing (neglecting the effect associated)
        for param, default in default_values.items():
            self.params.setdefault(param, default)
        
        return list(default_values.keys())  # Return the list of required parameters
    
    def solve(self):
        """Solve the compressor model."""

        # Check if the component is calculable and parametrized
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("CompressorSE could not be solved. It is not calculable and/or not parametrized")
            self.print_setup()
            return
        
        fluid = self.su.fluid  # Extract fluid name
        self.AS = AbstractState("HEOS", fluid)  # Create a reusable state object

        # Check if the fluid is in the two-phase region at the suction side
        self.AS.update(CoolProp.PQ_INPUTS, self.su.p, 1)
        self.T_sat_su = self.AS.T()
        if self.su.T < self.T_sat_su:
            print('----Warning the compressor inlet stream is not in vapor phase---')

        x_m_guess = [1.1, 0.99, 1.3, 1, 1.15] # guesses on the filling factor to provide suitable initial point for the iteration
        x_T_guess = [0.9, 0.8, 0.7, 0.2] # For the iteration on the T_w
        stop = 0 # Stop condition
        j = 0 # Index for the iteration on the T_w

        if self.params['mode'] == 'N_rot': # If the rotational speed is known
            try:
                # Loop to permit multiple attempts to solve the implicit calculation
                while not stop and j < len(x_T_guess):
                    k = 0
                    while not stop and k < len(x_m_guess):
                        # Guesses for the initial values
                        self.AS.update(CoolProp.PSmass_INPUTS, self.ex.p, self.su.s)
                        T_w_guess = x_T_guess[j]*self.AS.T()+5 # Guess on the wall temperature
                        m_dot_guess = x_m_guess[k]*self.params['V_s']*self.inputs['N_rot']/60*self.su.D # Guess on the mass flow rate
                        h_ex2_guess = self.AS.hmass() + 50000 # Guess on the enthalpy at the ex2
                        P_ex2_guess = 0.9*self.ex.p # Guess on the pressure at the ex2
                        #---------------------------------------------------------------------
                        args = ()
                        x = [T_w_guess, m_dot_guess, h_ex2_guess, P_ex2_guess]
                        #--------------------------------------------------------------------------
                        try:
                            fsolve(self.system, x, args = args)
                            res_norm = np.linalg.norm(self.res)
                        except:
                            res_norm = 1e6
                    
                        if res_norm < 1e-3:
                            stop = 1 # Stop the loop if the convergence is reached
                        k = k + 1
                    j = j + 1

            except Exception as e:
                print(f"CompressorSE could not be solved. Error: {e}")
                self.solved = False

        if self.params['mode'] == 'm_dot':
            try:
                # Loop to permit multiple attempts to solve the implicit calculation
                while not stop and j < len(x_T_guess):
                    # Guesses for the initial values
                    self.AS.update(CoolProp.PSmass_INPUTS, self.ex.p, self.su.s)
                    T_w_guess = x_T_guess[j]*self.AS.T()+5 # Guess on the wall temperature
                    h_ex2_guess = self.AS.hmass() + 50000 # Guess on the enthalpy at the ex2
                    P_ex2_guess = 0.9*self.ex.p # Guess on the pressure at the ex2
                    #---------------------------------------------------------------------
                    args = ()
                    x = [T_w_guess, h_ex2_guess, P_ex2_guess]
                    #--------------------------------------------------------------------------
                    try:
                        fsolve(self.system, x, args = args)
                        res_norm = np.linalg.norm(self.res)
                    except:
                        res_norm = 1e6
                
                    if res_norm < 1e-3:
                        stop = 1
                    j = j + 1

            except Exception as e:
                print(f"CompressorSE could not be solved. Error: {e}")
                self.solved = False

        self.convergence = stop
                
        if self.convergence: # If the calculation converged
            self.update_connectors() # Update the connectors with the calculated values
            self.solved = True


    def system(self, x):
        "Modelling section of the code"

        if self.params['mode'] == 'N_rot': # The rotational speed is given as an input
            self.T_w, self.m_dot, h_ex2_bis, P_ex2 = x # Values on which the system iterates
            self.N_rot = self.inputs['N_rot']
            #Boundary on the mass flow rate
            self.m_dot = max(self.m_dot, 1e-5)
        if self.params['mode'] == 'm_dot': # The mass flow rate is given as an input
            self.T_w, h_ex2_bis, P_ex2 = x # Values on which the system iterates
            self.m_dot = self.inputs['m_dot']
        
        #------------------------------------------------------------------------
        "1. Supply conditions: su"
        T_su = self.su.T
        P_su = self.su.p
        h_su = self.su.h
        s_su = self.su.s
        rho_su = self.su.D
        P_ex = self.ex.p

        #------------------------------------------------------------------------
        "2. Supply heating: su->su1"
        self.AS.update(CoolProp.HmassP_INPUTS, h_su, P_su)
        cp_su = self.AS.cpmass()
        if self.params['AU_su_n'] == 0: # If the heat transfer to the suction side is neglected
            Q_dot_su = 0
            h_su1 = h_su
            P_su1 = P_su
        else:
            C_dot_su = self.m_dot*cp_su
            AU_su = self.params['AU_su_n']*(self.m_dot/self.params['m_dot_n'])**0.8
            NTU_su = AU_su/C_dot_su
            epsilon_su = 1-np.exp(-NTU_su)
            Q_dot_su = epsilon_su*C_dot_su*(self.T_w-T_su)
    
            h_su1 = h_su + Q_dot_su/self.m_dot
            P_su1 = P_su

        #------------------------------------------------------------------------
        "3. Internal leakage: ex2->su1"
        self.AS.update(CoolProp.HmassP_INPUTS, h_ex2_bis, P_ex2)
        s_ex2_bis = self.AS.smass()
        x_ex2_bis = self.AS.Q()
        s_thr = s_ex2_bis

        if x_ex2_bis > 0 and x_ex2_bis < 1:
            # Two-phase
            self.AS.update(CoolProp.PQ_INPUTS, P_ex2, 0)
            T_ex2_bis = self.AS.T()
            cv_leak = self.AS.cvmass()
            cp_leak = self.AS.cpmass()
        else:
            # Single-phase
            T_ex2_bis = self.AS.T()
            cv_leak = self.AS.cvmass()
            cp_leak = self.AS.cpmass()
        
        gamma_ex = cp_leak/cv_leak
        P_thr_crit = P_ex2*(2/(gamma_ex+1))**(gamma_ex/(gamma_ex-1))
        P_thr = max(P_thr_crit, P_su1)
        self.AS.update(CoolProp.PSmass_INPUTS, P_thr, s_thr)
        rho_thr = self.AS.rhomass()
        h_thr = self.AS.hmass()
        C_thr = min(300, np.sqrt(2*(h_ex2_bis-h_thr)))
        V_dot_leak = self.params['A_leak']*C_thr
        m_dot_leak = V_dot_leak*rho_thr

        #------------------------------------------------------------------------
        "4. Flow rate calculation: su2"
        P_su2 = P_su1
        m_dot_in = m_dot_leak + self.m_dot
        h_su2 = max(min((self.m_dot*h_su1 + m_dot_leak*h_ex2_bis)/m_dot_in, h_ex2_bis), h_su1)
        self.AS.update(CoolProp.HmassP_INPUTS, h_su2, P_su2)
        rho_su2 = self.AS.rhomass()
        s_su2 = self.AS.smass()
        T_su2 = self.AS.T()
        self.N_rot_bis = m_dot_in/self.params['V_s']/rho_su2*60
        if self.params['mode'] == 'm_dot':
            self.N_rot = self.N_rot_bis

        #------------------------------------------------------------------------
        "5. Internal compression: su2->ex2"
        "Isentropic compression: su2->in"
        s_in = s_su2
        rho_in = rho_su2*self.params['rv_in']
        self.AS.update(CoolProp.DmassSmass_INPUTS, rho_in, s_in)
        P_in = self.AS.p()
        h_in = self.AS.hmass()
        w_in_is = h_in-h_su2
        
        "Isochoric compression: in->ex2"
        w_in_v = (P_ex2-P_in)/rho_in
        
        "Total internal work"
        w_in = w_in_is + w_in_v
        h_ex2 = h_su2 + w_in
        self.AS.update(CoolProp.HmassP_INPUTS, h_ex2, P_ex2)
        T_ex2 = self.AS.T()
        x_ex2 = self.AS.Q()
        
        #------------------------------------------------------------------------
        "6. Pressure drops: ex2->ex1"
        h_ex1 = h_ex2 #Isenthalpic valve
        if self.params['d_ex'] == 0: # If the pressure drop is neglected
            h_ex1_bis = h_ex2
            self.AS.update(CoolProp.HmassP_INPUTS, h_ex1_bis, P_ex)
            T_ex1 = self.AS.T()
        else:
            A_ex = np.pi*(self.params['d_ex']/2)**2
            P_crit_ex = P_ex*(2/(gamma_ex+1))**(gamma_ex/(gamma_ex-1))
            P_thr_ex = max(P_crit_ex, P_ex)
            self.AS.update(CoolProp.PSmass_INPUTS, P_thr_ex, s_ex2_bis)
            h_thr_ex = self.AS.hmass()
            v_thr_ex = 1./self.AS.rhomass()
            
            V_dot_ex = self.m_dot*v_thr_ex
            C_ex = V_dot_ex/A_ex
            h_ex1_bis = h_thr_ex + (C_ex**2)/2
            self.AS.update(CoolProp.HmassP_INPUTS, h_ex1_bis, P_ex)
            T_ex1 = self.AS.T()

        #------------------------------------------------------------------------           
        "7. Cooling at exit: ex1 -> ex "
        if self.params['AU_ex_n'] == 0: # If the heat transfer to the exhaust side is neglected
            Q_dot_ex = 0
            self.h_ex = h_ex1
        else:
            AU_ex = self.params['AU_ex_n']*(self.m_dot/self.params['m_dot_n'])**0.8
            if x_ex2 < 0 or x_ex2 >= 1:
                cp_ex = self.AS.cpmass()
                C_dot_ex = self.m_dot*cp_ex
            else:
                self.AS.update(CoolProp.PQ_INPUTS, P_ex, 1)
                cp_sat_v_ex = self.AS.cpmass()
                self.AS.update(CoolProp.PQ_INPUTS, P_ex, 0)
                cp_sat_l_ex = self.AS.cpmass()
                C_dot_ex = self.m_dot*((cp_sat_v_ex*x_ex2)+(1-x_ex2)*cp_sat_l_ex)
            NTU_ex = AU_ex/C_dot_ex
            epsilon_ex = 1 - np.exp(-NTU_ex)
            Q_dot_ex = epsilon_ex*C_dot_ex*(self.T_w-T_ex1)
            self.h_ex = h_ex1 + Q_dot_ex/self.m_dot
        
        #------------------------------------------------------------------------
        "8. Energy balance"
        # Fictious enveloppe heat balance
        self.Q_dot_amb = self.params['AU_amb']*(self.T_w-self.inputs['T_amb'])
        # Compression work and power
        W_dot_in = m_dot_in*w_in
        W_dot_loss = self.params['alpha']*W_dot_in + self.params['W_dot_loss_0'] + self.params['C_loss']*(self.N_rot/60)*2*np.pi
        self.W_dot_cp = W_dot_in + W_dot_loss
        
        "9. Performances"
        # Isentropic efficiency
        self.AS.update(CoolProp.PSmass_INPUTS, P_ex, s_su)
        h_ex_is = self.AS.hmass()
        w_s = h_ex_is-h_su
        W_dot_s = self.m_dot*w_s
        self.epsilon_is = W_dot_s/self.W_dot_cp
        
        #Volumetric efficiency
        #Theoretical flowrate
        V_s_dot = self.params['V_s']*self.N_rot/60
        m_dot_th = V_s_dot*rho_su
        
        m_dot_in_bis = V_s_dot*rho_su2
        
        # Volumetric efficiencies definitions
        self.epsilon_v = self.m_dot/m_dot_th
        self.epsilon_v_l = self.m_dot/m_dot_in
        self.epsilon_v_PT = m_dot_in/m_dot_th
        
        "10. Residue"
        self.resE = abs((W_dot_loss - Q_dot_ex - Q_dot_su - self.Q_dot_amb)/(W_dot_loss))
        self.res_h_ex1 = abs(h_ex1_bis-h_ex1)/h_ex1
        self.res_h_ex2 = abs((h_ex2_bis-h_ex2)/h_ex2)
        if self.params['mode'] == 'N_rot':
            self.res_Nrot = abs((self.inputs['N_rot']-self.N_rot_bis)/self.inputs['N_rot'])
            self.res = [self.resE, self.res_h_ex1, self.res_h_ex2, self.res_Nrot]
        if self.params['mode'] == 'm_dot':
            self.res = [self.resE, self.res_h_ex1, self.res_h_ex2]

        return self.res
    
    def update_connectors(self):
        """Update the connectors with the calculated values."""

        self.ex.set_fluid(self.su.fluid)
        self.ex.set_h(self.h_ex)
        self.ex.set_m_dot(self.m_dot)

        self.W_mec.set_W_dot(self.W_dot_cp)
        if self.params['mode'] == 'N_rot':
            self.su.set_m_dot(self.m_dot)
        if self.params['mode'] == 'm_dot':
            self.W_mec.set_N(self.N_rot)
            
        self.Q_amb.set_Q_dot(self.Q_dot_amb)
        self.Q_amb.set_T_hot(self.T_w)

    def print_results(self):
        print("=== Expander Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_cp: {self.W_mec.W_dot} [W]")
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
        print(f"  - W_dot_cp: {self.W_mec.W_dot} [W]")
        print("=========================")
        print("Heat connector:")
        print(f"  - Q_dot_amb: {self.Q_amb.Q_dot} [W]")
        print(f"  - T_hot: {self.Q_amb.T_hot} [K]")
        print(f"  - T_cold: {self.Q_amb.T_cold} [K]")
        print("=========================")
