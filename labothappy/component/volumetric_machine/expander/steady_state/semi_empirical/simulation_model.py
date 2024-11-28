# -*- coding: utf-8 -*-

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np
import time

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
            self.sync_inputs()
            # Return a list of required inputs
            return ['su_p', 'su_T', 'ex_p', 'N_rot', 'T_amb', 'su_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p
        if self.W_exp.N is not None:
            self.inputs['N_rot'] = self.W_exp.N
        if self.Q_amb.T_cold is not None:
            self.inputs['T_amb'] = self.Q_amb.T_cold

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])
        if 'N_rot' is self.inputs:
            self.W_exp.set_N(self.inputs['N_rot'])
        if 'T_amb' is self.inputs:
            self.Q_amb.set_T_cold(self.inputs['T_amb'])


    def get_required_parameters(self):
        return [
            'AU_amb', 'AU_su_n', 'AU_ex_n', 'd_su1', 'm_dot_n', 
            'A_leak', 'W_dot_loss_0', 'alpha', 'C_loss', 'rv_in', 'V_s'
        ]
    
    def print_setup(self):
        print("=== Expander Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - W_exp: N={self.W_exp.N}, W={self.W_exp.W_dot}")
        print(f"  - Q_amb: T_cold={self.Q_amb.T_cold}, T_hot={self.Q_amb.T_hot}, Q={self.Q_amb.Q_dot}")

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

    def System(self, x):
        self.m_dot, self.T_w = x
        self.N = self.inputs['N_rot'] /60
            
        #Boundary on the mass flow rate
        self.m_dot = max(self.m_dot, 1e-5)
        
        "1. Supply conditions: su"
        T_su = self.su.T
        P_su = self.su.p
        h_su = self.su.h
        s_su = self.su.s
        rho_su = self.su.D
        Fluid = self.su.fluid
        P_ex = self.ex.p
        self.P_ex = P_ex

        h_ex_is = PropsSI('H', 'P', P_ex, 'S', s_su, Fluid)
        h_max = PropsSI('H', 'P', 4e6, 'T', 500, Fluid)

        T_sat_su = PropsSI('T', 'P', P_su, 'Q', 1, Fluid)
        if T_su<T_sat_su:
            print('----Warning the expander inlet stream is not in vapor phase---')

        "2. Supply pressure drop: su->su1"
        h_su1 = h_su #Isenthalpic valve
        # Assumption that the density doesn't change too much
        A_su = np.pi*(self.params['d_su1']/2)**2
        V_dot_su = self.m_dot/rho_su
        C_su = V_dot_su/A_su
        h_su_thr1 = h_su-(C_su**2)/2
        h_su_thr = max(h_su_thr1, h_ex_is)
        P_su_thr = PropsSI('P', 'S', s_su, 'H', h_su_thr, Fluid)
        P_su1 = max(P_su_thr, P_ex+1)
        self.DP_su = P_su-P_su1
        T_su1 = PropsSI('T', 'P', P_su1, 'H', h_su1, Fluid)

        "3. Cooling at the entrance: su1->su2"  
        try:
            cp_su1 = PropsSI('CPMASS', 'P', P_su1, 'H', h_su1, Fluid)
        except:
            cp_su1 = PropsSI('CPMASS', 'P', P_su1, 'Q', 0, Fluid)
        AU_su = self.params['AU_su_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
        C_dot_su = self.m_dot*cp_su1
        NTU_su = AU_su/C_dot_su
        epsilon_su1 = max(0,1-np.exp(-NTU_su))
        Q_dot_su = max(0, epsilon_su1*self.m_dot*cp_su1*(T_su1-self.T_w))
        
        h_su2 = min(h_max, max(max(h_ex_is, PropsSI('H','Q', 0.1, 'P', P_su1, Fluid)), h_su1 - Q_dot_su/self.m_dot))
        P_su2 = P_su1 #No pressure drop just heat transfer
        
        rho_su2 = PropsSI('D', 'P', P_su2, 'H', h_su2, Fluid)
        T_su2 = PropsSI('T', 'P', P_su2, 'H', h_su2, Fluid)
        s_su2 = PropsSI('S', 'P', P_su2, 'H', h_su2, Fluid)
        
        "4. Leakage"
        try:
            cv_su1 = PropsSI('CVMASS', 'P', P_su1, 'H', h_su1, Fluid)
        except:
            cv_su1 = PropsSI('CVMASS', 'P', P_su1, 'Q', 0, Fluid)
        gamma = max(1e-2, cp_su1/cv_su1)
        P_crit = P_su2*(2/(gamma+1))**(gamma/(gamma-1))
        P_thr_leak = max(P_ex, P_crit)
        rho_thr_leak = PropsSI('D', 'P', P_thr_leak, 'S', s_su2, Fluid)
        h_thr_leak = PropsSI('H', 'P', P_thr_leak, 'S', s_su2, Fluid)
        C_thr_leak = np.sqrt(2*(h_su2-h_thr_leak))
        m_dot_leak = self.params['A_leak']*C_thr_leak*rho_thr_leak
        if self.su.m_dot == None:
            m_dot_in = self.N*self.params['V_s']*rho_su2
            m_dot_leak_bis = self.m_dot-m_dot_in
        elif self.su.m_dot != None:
            m_dot_in = self.m_dot-m_dot_leak
            if self.params['N_rot'] == None:
                self.N = m_dot_in/(self.params['V_s']*rho_su2)
        

        "5. Internal expansion"
        "Isentropic expansion to the internal pressure: su2->in"
        rho_in = rho_su2/self.params['rv_in']
        #Trick to not have problems with CoolProp
        try:
            P_in = PropsSI('P', 'D', rho_in, 'S', s_su2, Fluid)
        except:
            delta = 0.0001
            P_in = 0.5*PropsSI('P', 'D', rho_in*(1+delta), 'S', s_su2, Fluid)+0.5*PropsSI('P', 'D', rho_in*(1-delta), 'S', s_su2, Fluid)
            
        h_in = PropsSI('H', 'D', rho_in, 'P', P_in, Fluid)
        w_in_s = h_su2-h_in
        "Expansion at constant volume: in->ex2"
        w_in_v = (P_in-P_ex)/rho_in
        h_ex2 = h_in-w_in_v
        "Total work"
        W_dot_in = m_dot_in*(w_in_s+w_in_v)
        
        "6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
        h_ex1 = max(min((m_dot_in*h_ex2 + m_dot_leak*h_su2)/self.m_dot, h_su2), h_ex2)
        P_ex1 = P_ex
        T_ex1 = PropsSI('T', 'P', P_ex1, 'H', h_ex1, Fluid)
        try:
            cp_ex2 = PropsSI('CVMASS', 'P', P_ex1, 'H', h_ex1, Fluid)
        except:
            cp_ex2 = PropsSI('CVMASS', 'P', P_ex1, 'Q', 0, Fluid)
        AU_ex = self.params['AU_ex_n']*(self.m_dot/self.params['m_dot_n'])**(0.8)
        C_dot_ex = self.m_dot*cp_ex2
        NTU_ex=AU_ex/C_dot_ex
        epsilon_ex = max(0, 1-np.exp(-NTU_ex))
        Q_dot_ex = max(0, epsilon_ex*C_dot_ex*(self.T_w-T_ex1))
        self.h_ex = h_ex1+Q_dot_ex/self.m_dot

        "8. Energy balance"
        self.Q_dot_amb = self.params['AU_amb']*(self.T_w-self.inputs['T_amb'])
        W_dot_loss = self.params['alpha']*W_dot_in + self.params['W_dot_loss_0'] + self.params['C_loss']*self.N*2*np.pi
        self.W_dot_exp = W_dot_in - W_dot_loss

        "9. Performances"
        W_dot_s = self.m_dot*(h_su-h_ex_is)
        self.epsilon_is = self.W_dot_exp/W_dot_s
        self.m_dot_th = self.N*self.params['V_s']*rho_su
        self.epsilon_v = self.m_dot/self.m_dot_th
        
        "10. Residuals"
        self.res_E = abs((Q_dot_su + W_dot_loss - Q_dot_ex - self.Q_dot_amb)/(Q_dot_su + W_dot_loss))
        self.res = self.res_E
        self.res_m_leak = abs((m_dot_leak_bis-m_dot_leak)/m_dot_leak)
        self.res = [self.res, self.res_m_leak]
        return self.res

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("ExpanderSE could not be solved. It is not calculable and/or not parametrized")
            return

        try:
            start_time = time.time()
            ff_guess = [0.7, 1.2, 0.8, 1.3, 0.4, 1.7, 3]
            x_T_guess = [0.7, 0.95, 0.8, 0.99, 0.9]
            stop = 0
            j = 0
            while not stop and j < len(x_T_guess):
                k = 0
                while not stop and k < len(ff_guess):
                    # Loop to permit multiple attempts to solve the implicit calculation 
                    m_dot_guess = ff_guess[k] * self.params['V_s'] * self.inputs['N_rot'] / 60 * PropsSI('D', 'P', self.inputs['su_p'], 'T', self.inputs['su_T'], self.inputs['su_fluid'])
                    T_w_guess = x_T_guess[j] * self.inputs['su_T'] + (1 - x_T_guess[j]) * self.inputs['T_amb']
                    #---------------------------------------------------------------------
                    args = ()
                    x = [m_dot_guess, T_w_guess]
                    #--------------------------------------------------------------------------
                    try:
                        fsolve(self.System, x, args=args)
                        res_norm = np.linalg.norm(self.res)
                    except:
                        res_norm = 1
                    if res_norm < 1e-4:
                        stop = 1
                    k += 1
                j += 1

            self.convergence = stop

            if self.convergence:
                self.update_connectors()
                self.solved = True

        except Exception as e:
            print(f"ExpanderSE could not be solved. Error: {e}")
            self.solved = False

    def update_connectors(self):
        """Update the connectors with the calculated values."""
        self.su.set_m_dot(self.m_dot)
        
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.m_dot)
        self.ex.set_h(self.h_ex)
        self.ex.set_p(self.P_ex)

        self.W_exp.set_W_dot(self.W_dot_exp)
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



