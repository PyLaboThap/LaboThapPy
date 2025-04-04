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

"Internal modules"
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

class ExpanderSE(BaseComponent):
    """
        Component: Volumetric expander

        Model: The model is based on the thesis of V. Lemort (2008) and is a semi-empirical model.

        Reference: 

        **Descritpion**:

            This model is used to simulate the performance of a volumetric expander. 
            The parameters of the model need to be calibrated with experimental datas to represent the real behavior of the expander.

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

            d_su1: Pressure drop diameter. [m]

            m_dot_n: Nominal mass flow rate. [kg/s]

            A_leak: Leakage area. [m^2]

            W_dot_loss_0: Constant loss in the compressor. [W]

            alpha: Proportionality rate for mechanical losses. [-]

            C_loss: Torque losses. [N.m]

            rv_in: Inlet volume ratio. [-]

            V_s: Swept volume. [m^3]

            mode: Mode of operation ('N_rot' if N_rot is given in the inputs or 'm_dot' if m_dot is given in the inputs).

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
        self.su = MassConnector() # Suction side mass connector
        self.ex = MassConnector() # Exhaust side mass connector
        self.W_mec = WorkConnector() # Work connector of the expander
        self.Q_amb = HeatConnector() # Heat connector to the ambient

    def get_required_inputs(self):
        # Define the required inputs for the component
        # If the mode is 'N_rot', the rotational speed of the expander is required while the mass flow rate is calculated
        if self.params['mode'] == 'N_rot':
            return ['P_su', 'h_su', 'P_ex', 'N_rot', 'T_amb', 'fluid']
        # If the mode is 'm_dot', the mass flow rate is required while the rotational speed of the expander is calculated
        elif self.params['mode'] == 'm_dot':
            return ['P_su', 'h_su', 'P_ex', 'm_dot', 'T_amb', 'fluid']
    
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

        # Ensure all required parameters are set, assigning default values if missing (neglecting the effect associated)
        for param, default in default_values.items():
            self.params.setdefault(param, default)
        
        return list(default_values.keys())  # Return the list of required parameters

    def solve(self):
        """Solve the expander model."""

        # Check if the component is calculable and parametrized
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized): # If the component is not calculable and/or not parametrized
            self.solved = False
            print("ExpanderSE could not be solved. It is not calculable and/or not parametrized")
            self.print_setup()
            return

        fluid = self.su.fluid  # Extract fluid name
        self.AS = AbstractState("HEOS", fluid)  # Create a reusable state object

        # Check if the fluid is in the two-phase region at the suction side
        self.AS.update(CoolProp.PQ_INPUTS, self.su.p, 1)  
        T_sat_su = self.AS.T()
        if self.su.T<T_sat_su:
            print('----Warning the expander inlet stream is not in vapor phase---')

        ff_guess = [0.7, 1.2, 0.8, 1.3, 0.4, 1.7] # Guesses for the filling factor (ff)
        x_T_guess = [0.7, 0.95, 0.8, 0.9] # Guesses for the temperature ratio (x_T)
        stop = 0 # Stop flag
        j = 0 # Index for the temperature ratio

        if self.params['mode'] == 'N_rot': # The rotational speed is given as an input
            try: 
                # Loop to permit multiple attempts to solve the implicit calculation 
                # Loop stops when the residuals are below the threshold
                while not stop and j < len(x_T_guess):
                    k = 0 # Index for the filling factor
                    while not stop and k < len(ff_guess):
                        # Loop to permit multiple attempts to solve the implicit calculation 
                        self.AS.update(CoolProp.PT_INPUTS, self.inputs['P_su'], self.su.T)
                        m_dot_guess = ff_guess[k] * self.params['V_s'] * self.inputs['N_rot'] / 60 * self.AS.rhomass() # Guess for the mass flow rate
                        T_w_guess = x_T_guess[j] * self.su.T + (1 - x_T_guess[j]) * self.inputs['T_amb'] # Guess for the wall temperature
                        #---------------------------------------------------------------------
                        args = ()
                        x = [m_dot_guess, T_w_guess]
                        #--------------------------------------------------------------------------
                        try:
                            fsolve(self.System, x, args=args) # Solve the system of equations
                            res_norm = np.linalg.norm(self.res) # Calculate the norm of the residuals
                        except:
                            res_norm = 1
                        if res_norm < 1e-4:
                            stop = 1 # Stop the loop if the norm of the residuals are below the threshold
                        k += 1
                    j += 1
            except Exception as e:
                print(f"ExpanderSE could not be solved. Error: {e}")
                self.solved = False # Unable to solve the component
        
        if self.params['mode'] == 'm_dot': # The mass flow rate is given as an input
            try:  
                while not stop and j < len(x_T_guess):
                    T_w_guess = x_T_guess[j] * self.su.T + (1 - x_T_guess[j]) * self.inputs['T_amb'] # Guess for the wall temperature
                    #---------------------------------------------------------------------
                    args = ()
                    x = [T_w_guess]
                    #--------------------------------------------------------------------------
                    try:
                        fsolve(self.system, x, args=args) # Solve the system of equations
                        res_norm = np.linalg.norm(self.res) # Calculate the norm of the residuals
                    except:
                        res_norm = 1
                    if res_norm < 1e-4:
                        stop = 1 # Stop the loop if the norm of the residuals are below the threshold
                    j += 1

            except Exception as e:
                print(f"ExpanderSE could not be solved. Error: {e}")
                self.solved = False # Unable to solve the component
    
        self.convergence = stop

        if self.convergence: # If the component is solved
            self.update_connectors() # Update the connectors with the calculated values
            self.solved = True


    def system(self, x):
        """System of equations to solve the expander model."""

        if self.params['mode'] == 'N_rot': # The rotational speed is given as an input
            self.m_dot, self.T_w = x # Values on which the system iterates
            self.N_rot = self.inputs['N_rot']
            # Boundary on the mass flow rate guess
            self.m_dot = max(self.m_dot, 1e-5)
        elif self.params['mode'] == 'm_dot': # The mass flow rate is given as an input
            self.m_dot = self.inputs['m_dot']
            self.T_w = x[0] # Values on which the system iterates
 
        #------------------------------------------------------------------------------------------------
        "1. Supply conditions: su"
        T_su = self.su.T
        P_su = self.su.p
        h_su = self.su.h
        s_su = self.su.s
        rho_su = self.su.D
        P_ex = self.ex.p

        self.AS.update(CoolProp.PSmass_INPUTS, P_ex, s_su)
        h_ex_is = self.AS.hmass()
        self.AS.update(CoolProp.PT_INPUTS, 4e6, 500)
        h_max = self.AS.hmass()
        
        #------------------------------------------------------------------------------------------------
        "2. Supply pressure drop: su->su1"
        h_su1 = h_su #Isenthalpic valve

        if self.params['d_su1'] == 0: # No pressure drop
            P_su1 = P_su
            T_su1 = T_su
        else: # Pressure drop
            A_su = np.pi*(self.params['d_su1']/2)**2
            V_dot_su = self.m_dot/rho_su # Assumption that the density doesn't change too much
            C_su = V_dot_su/A_su
            h_su_thr1 = h_su-(C_su**2)/2
            h_su_thr = max(h_su_thr1, h_ex_is)
            self.AS.update(CoolProp.HmassSmass_INPUTS, h_su_thr, s_su)
            P_su_thr = self.AS.p()
            P_su1 = max(P_su_thr, P_ex+1)
            self.DP_su = P_su-P_su1
            self.AS.update(CoolProp.HmassP_INPUTS, h_su1, P_su1)
            T_su1 = self.AS.T()

        #------------------------------------------------------------------------------------------------
        "3. Cooling at the entrance: su1->su2"  
        try:
            self.AS.update(CoolProp.HmassP_INPUTS, h_su1, P_su1)
            cp_su1 = self.AS.cpmass()
        except:
            self.AS.update(CoolProp.PQ_INPUTS, P_su1, 0)
            cp_su1 = self.AS.cpmass()
        
        if self.params['AU_su_n'] == 0: # No heat transfer
            h_su2 = h_su1
            Q_dot_su = 0
        else: # Heat transfer          
            AU_su = self.params['AU_su_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
            C_dot_su = self.m_dot*cp_su1
            NTU_su = AU_su/C_dot_su
            epsilon_su1 = max(0,1-np.exp(-NTU_su))
            Q_dot_su = max(0, epsilon_su1*self.m_dot*cp_su1*(T_su1-self.T_w))
            self.AS.update(CoolProp.PQ_INPUTS, P_su1, 0.1)
            h_su2 = min(h_max, max(max(h_ex_is, self.AS.hmass()), h_su1 - Q_dot_su/self.m_dot))
        
        P_su2 = P_su1 # No pressure drop just heat transfer
        self.AS.update(CoolProp.HmassP_INPUTS, h_su2, P_su2)
        rho_su2 = self.AS.rhomass()
        s_su2 = self.AS.smass()
        
        #------------------------------------------------------------------------------------------------
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
            m_dot_leak_bis = self.m_dot-m_dot_in # residual on the leakage flow rate to get the right mass flow rate
        elif self.params['mode'] == 'm_dot':
            m_dot_in = self.m_dot-m_dot_leak
            self.N_rot = (m_dot_in/(self.params['V_s']*rho_su2))*60 # rotational speed calculation
        
        #------------------------------------------------------------------------------------------------
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
        "Total internal work"
        W_dot_in = m_dot_in*(w_in_s+w_in_v)

        #------------------------------------------------------------------------------------------------
        "6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
        h_ex1 = max(min((m_dot_in*h_ex2 + m_dot_leak*h_su2)/self.m_dot, h_su2), h_ex2) # Adiabatic mixing
        P_ex1 = P_ex
        self.AS.update(CoolProp.HmassP_INPUTS, h_ex1, P_ex1)
        T_ex1 = self.AS.T()
        try:
            self.AS.update(CoolProp.HmassP_INPUTS, h_ex1, P_ex1)
            cp_ex2 = self.AS.cpmass()
        except:
            self.AS.update(CoolProp.PQ_INPUTS, P_ex1, 0)
            cp_ex2 = self.AS.cpmass()
        
        #------------------------------------------------------------------------------------------------
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

        #------------------------------------------------------------------------------------------------
        "8. Energy balance"
        self.Q_dot_amb = self.params['AU_amb']*(self.T_w-self.inputs['T_amb'])
        W_dot_loss = self.params['alpha']*W_dot_in + self.params['W_dot_loss_0'] + self.params['C_loss']*(self.N_rot/60)*2*np.pi
        self.W_dot_exp = W_dot_in - W_dot_loss

        #------------------------------------------------------------------------------------------------
        "9. Performances"
        W_dot_s = self.m_dot*(h_su-h_ex_is)
        self.epsilon_is = self.W_dot_exp/W_dot_s
        self.m_dot_th = (self.N_rot/60)*self.params['V_s']*rho_su
        self.epsilon_v = self.m_dot/self.m_dot_th
        
        #------------------------------------------------------------------------------------------------
        "10. Residuals"
        self.res_E = abs((Q_dot_su + W_dot_loss - Q_dot_ex - self.Q_dot_amb)/(Q_dot_su + W_dot_loss))
        if self.params['mode'] == 'N_rot':
            self.res_m_leak = abs((m_dot_leak_bis-m_dot_leak)/m_dot_leak)
            self.res = [self.res_E, self.res_m_leak]
        elif self.params['mode'] == 'm_dot':
            self.res = [self.res_E]
        return self.res

    def update_connectors(self):
        """Update the connectors with the calculated values."""
        self.su.set_m_dot(self.m_dot)
        
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.m_dot)
        self.ex.set_h(self.h_ex)

        self.W_mec.set_W_dot(self.W_dot_exp)
        self.W_mec.set_N(self.N_rot)
        self.Q_amb.set_Q_dot(self.Q_dot_amb)
        self.Q_amb.set_T_hot(self.T_w)

    def print_results(self):
        print("=== Expander Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_exp: {self.W_mec.W_dot} [W]")
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
        print(f"  - W_dot_exp: {self.W_mec.W_dot} [W]")
        print("=========================")
        print("Heat connector:")
        print(f"  - Q_dot_amb: {self.Q_amb.Q_dot} [W]")
        print(f"  - T_hot: {self.Q_amb.T_hot} [K]")
        print(f"  - T_cold: {self.Q_amb.T_cold} [K]")
        print("=========================")