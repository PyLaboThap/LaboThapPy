# -*- coding: utf-8 -*-
"""
Aug 08 14:32:39 2024

@author: Elise Neven
@email: elise.neven@uliege.be
"""

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

from component.heat_exchanger.moving_boundary.simple_model.modules.U import U_Gnielinski_calibrated, U_DittusBoelter, U_Cooper_calibrater, U_Thonon

from CoolProp.CoolProp import PropsSI
import numpy as np
import math

class HXSimpleMB(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su_wf = MassConnector() # Working fluid supply
        self.su_sf = MassConnector() # Secondary fluid supply
        self.ex_wf = MassConnector()
        self.ex_sf = MassConnector()
        self.Q_dot = HeatConnector()

    def get_required_inputs(self):
        if self.inputs == {}:
            # Hot Fluid
            if self.su_wf.fluid is not None:
                self.inputs['fluid_wf'] = self.su_wf.fluid
            if self.su_wf.T is not None:
                self.inputs['su_wf_T'] = self.su_wf.T
            if self.su_wf.m_dot is not None:
                self.inputs['su_wf_m_dot'] = self.su_wf.m_dot
            if self.su_wf.x is not None:
                self.inputs['su_wf_x'] = self.su_wf.x
            if self.su_sf.fluid is not None:
                self.inputs['fluid_sf'] = self.su_sf.fluid
            if self.su_sf.T is not None:
                self.inputs['su_sf_T'] = self.su_sf.T
            if self.su_sf.p is not None:
                self.inputs['su_sf_p'] = self.su_sf.p
            if self.su_sf.m_dot is not None:
                self.inputs['su_sf_m_dot'] = self.su_sf.m_dot
            if self.ex_sf.T is not None:
                self.inputs['ex_sf_T'] = self.ex_sf.T
        
        print(self.inputs)
        if self.inputs != {}:
            # Working Fluid
            if 'fluid_wf' in self.inputs:
                self.su_wf.set_fluid(self.inputs['fluid_wf'])
            if 'su_wf_T' in self.inputs:
                self.su_wf.set_T(self.inputs['su_wf_T'])
            if 'su_wf_m_dot' in self.inputs:
                self.su_wf.set_m_dot(self.inputs['su_wf_m_dot'])
            if 'su_wf_x' in self.inputs:
                self.su_wf.set_x(self.inputs['su_wf_x'])
            # Secondary Fluid
            if 'fluid_sf' in self.inputs:
                self.su_sf.set_fluid(self.inputs['fluid_sf'])
            if 'su_sf_T' in self.inputs:
                self.su_sf.set_T(self.inputs['su_sf_T'])
            if 'su_sf_p' in self.inputs:
                self.su_sf.set_p(self.inputs['su_sf_p'])
            if 'su_sf_m_dot' in self.inputs:
                self.su_sf.set_m_dot(self.inputs['su_sf_m_dot'])
            if 'ex_sf_T' in self.inputs:
                self.ex_sf.set_T(self.inputs['ex_sf_T'])
            
        return['fluid_wf', 'su_wf_T', 'su_wf_m_dot', 'su_wf_x', 'fluid_sf', 'su_sf_T', 'su_sf_p', 'su_sf_m_dot', 'ex_sf_T']
    
    def get_required_parameters(self):
        return [
            'HX_type', # Type of heat exchanger
            'HX_D', # Diameter of the heat exchanger
            'HX_A', # Area of the heat exchanger
            'min_pinch', # Minimum pinch point
            'Delta_T_sup_or_sub', # Superheating or subcooling
        ]
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - su_wf: fluid={self.su_wf.fluid}, T={self.su_wf.T}, p={self.su_wf.p}, m_dot={self.su_wf.m_dot}")
        print(f"  - su_sf: fluid={self.su_sf.fluid}, T={self.su_sf.T}, p={self.su_sf.p}, m_dot={self.su_sf.m_dot}")
        print(f"  - ex_wf: fluid={self.ex_wf.fluid}, T={self.ex_wf.T}, p={self.ex_wf.p}, m_dot={self.ex_wf.m_dot}")
        print(f"  - ex_sf: fluid={self.ex_sf.fluid}, T={self.ex_sf.T}, p={self.ex_sf.p}, m_dot={self.ex_sf.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")

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

    def MB_model(self, P_wf):
        # TEMPERATURES AND THERMAL POWERS
        # evaporating/consensing phase - wf side
        T_wf_sat = PropsSI('T', 'P', P_wf, 'Q', 0, self.su_wf.fluid)

        # case depending temperatures
        if self.params['HX_type'] == 'evaporator':
            sf_ht_dir = 'cooled'
            wf_ht_dir = 'heated'
            T_sf_h = self.su_sf.T
            T_sf_c = self.ex_sf.T    #T_sf_h - self.glide_sf
            T_wf_h = T_wf_sat + self.params['Delta_T_sup_or_sub']
            T_wf_c = self.su_wf.T
            self.T_wf_out = T_wf_h
        elif self.params['HX_type'] == 'condenser':
            sf_ht_dir = 'heated'
            wf_ht_dir = 'cooled'
            T_sf_c = self.su_sf.T
            T_sf_h = self.ex_sf.T   #T_sf_c + self.glide_sf
            T_wf_c = T_wf_sat - self.params['Delta_T_sup_or_sub']
            T_wf_h = self.su_wf.T
            self.T_wf_out = T_wf_c

        # Energy balance - wf side
        if self.params['HX_type'] == 'evaporator':
            if self.su_wf.x >= 0 and self.su_wf.x <= 1:
                h_wf_sat_liq = PropsSI('H','P', P_wf,'Q', self.su_wf.x, self.su_wf.fluid)
                h_wf_c = h_wf_sat_liq
            else:
                h_wf_sat_liq = PropsSI('H','P', P_wf,'Q', 0, self.su_wf.fluid)
                h_wf_c = PropsSI('H','P', P_wf,'T', T_wf_c, self.su_wf.fluid)
            h_wf_sat_vap = PropsSI('H','P', P_wf,'Q', 1, self.su_wf.fluid)
            h_wf_h = PropsSI('H','P', P_wf,'T', T_wf_h, self.su_wf.fluid)

        if self.params['HX_type'] == 'condenser':
            if self.su_wf.x >= 0 and self.su_wf.x <= 1:
                h_wf_sat_vap = PropsSI('H','P', P_wf,'Q', self.su_wf.x, self.su_wf.fluid)
                h_wf_h = h_wf_sat_vap
            else:
                h_wf_sat_vap = PropsSI('H','P', P_wf,'Q', 1, self.su_wf.fluid)
                h_wf_h = PropsSI('H','P', P_wf,'T', T_wf_h, self.su_wf.fluid)
            h_wf_sat_liq = PropsSI('H','P', P_wf,'Q', 0, self.su_wf.fluid)
            h_wf_c = PropsSI('H','P', P_wf,'T', T_wf_c, self.su_wf.fluid)

        Q_sub = max(0, self.su_wf.m_dot*(h_wf_sat_liq-h_wf_c))
        Q_sat = max(0, self.su_wf.m_dot*(h_wf_sat_vap-h_wf_sat_liq))
        Q_sup = max(0, self.su_wf.m_dot*(h_wf_h-h_wf_sat_vap))

        self.Q = Q_sub + Q_sat + Q_sup

        # Energy balance - sf side
        T_sf_mean = (T_sf_h + T_sf_c)/2
        cp_sf = PropsSI('C','T', T_sf_mean,'P', self.su_sf.p, self.su_sf.fluid)
        self.m_dot_sf = self.Q/(cp_sf*(T_sf_h - T_sf_c))
        T_sf_ch = T_sf_c + Q_sub/(self.m_dot_sf*cp_sf)
        T_sf_hc = T_sf_h - Q_sup/(self.m_dot_sf*cp_sf)

        # FLUID MASS DISTRIBUTION
        # Heat tranfer coefficients 1P - sf side: Dittus_Boelter
        U_sf = U_DittusBoelter(self.m_dot_sf, sf_ht_dir, self.params['HX_D'], self.su_sf.fluid, self.su_sf.p, T_sf_mean)
        # Heat tranfer coefficients - wf side
        if self.params['HX_type'] == 'condenser':
            try:
                self.U_wf_h = U_Thonon(self.su_wf.m_dot, self.params['HX_D'], self.su_wf.fluid, P_wf, (T_wf_h + T_wf_sat)/2)
            except:
                self.U_wf_h = 0 #On arrive déjà en 2 phases dans le condenseur
            self.U_wf_c = U_Thonon(self.su_wf.m_dot, self.params['HX_D'], self.su_wf.fluid, P_wf, (T_wf_c + T_wf_sat)/2)
            self.U_wf_2P = U_Cooper_calibrater(self.Q, self.params['HX_A'], P_wf, self.su_wf.fluid)

        elif self.params['HX_type'] == 'evaporator':
            self.U_wf_h = U_Thonon(self.su_wf.m_dot, self.params['HX_D'], self.su_wf.fluid, P_wf, (T_wf_h + T_wf_sat)/2)
            try:
                self.U_wf_c = U_Thonon(self.su_wf.m_dot, self.params['HX_D'], self.su_wf.fluid, P_wf, (T_wf_c + T_wf_sat)/2)
            except:
                self.U_wf_c = 0 #On arrive déjà en 2 phases dans l'évaporateur
            self.U_wf_2P = U_Gnielinski_calibrated(self.su_wf.m_dot, self.params['HX_D'], self.su_wf.fluid, P_wf)

        try:
            self.U_sub = (1/U_sf + 1/self.U_wf_c)**(-1)
        except:
            self.U_sub = 0

        self.U_sat = (1/U_sf + 1/self.U_wf_2P)**(-1)

        try:
            self.U_sup = (1/U_sf + 1/self.U_wf_h)**(-1)
        except:
            self.U_sup = 0 #On arrive déjà en 2 phases

        # Mean logarithm temperature difference
        if self.params['HX_type'] == 'evaporator':
            tau_sub = [(T_sf_c-T_wf_c), (T_sf_ch-T_wf_sat)]
            tau_sat = [(T_sf_ch-T_wf_sat), (T_sf_hc-T_wf_sat)]
            tau_sup = [(T_sf_hc-T_wf_sat), (T_sf_hc-T_wf_h)]
        elif self.params['HX_type'] == 'condenser':
            tau_sub = [-(T_sf_c-T_wf_c), -(T_sf_ch-T_wf_sat)]
            tau_sat = [-(T_sf_ch-T_wf_sat), -(T_sf_hc-T_wf_sat)]
            tau_sup = [-(T_sf_hc-T_wf_sat), -(T_sf_h-T_wf_h)]

        if min(tau_sub) * max(tau_sub) < 0: #If are of opposite signs
            dT_ml_sub = 0
            self.flag_tau = 1
        else:
            dT_ml_sub = (max(tau_sub)-min(tau_sub))/np.log(max(tau_sub)/min(tau_sub))
            self.flag_tau = 0
        
        if min(tau_sat) * max(tau_sat) < 0: #If are of opposite signs
            dT_ml_sat = 0
            self.flag_tau = 1
        else:
            dT_ml_sat = (max(tau_sat)-min(tau_sat))/np.log(max(tau_sat)/min(tau_sat))
        
        if min(tau_sup) * max(tau_sup) < 0: #If are of opposite signs
            dT_ml_sup = 0
            self.flag_tau = 1
        else:
            dT_ml_sup = (max(tau_sup)-min(tau_sup))/np.log(max(tau_sup)/min(tau_sup))
            self.flag_tau = 0

        self.tau = [tau_sub, tau_sat, tau_sup]

        # Energy balance UA
        if Q_sub == 0:
            self.A_sub = 0
        else:
            self.A_sub = Q_sub/(self.U_sub*dT_ml_sub)
        if Q_sat == 0:
            self.A_sat = 0
        else:
            self.A_sat = Q_sat/(self.U_sat*dT_ml_sat)
        
        if Q_sup == 0: #Déja en 2 phases
            self.A_sup = 0
        else:
            self.A_sup = Q_sup/(self.U_sup*dT_ml_sup)

        self.A_sat = abs(self.params['HX_A'] - self.A_sub - self.A_sup)

        # Closure equations
        Q_sat_calc = self.A_sat*self.U_sat*dT_ml_sat
        self.res = abs(Q_sat-Q_sat_calc)

        "Results for TQ plot"
        self.Q_wf_plot = [0, Q_sub, Q_sub+Q_sat, self.Q]
        self.Q_sf_plot = [0, self.Q]
        self.T_wf_plot = [T_wf_c, T_wf_sat, T_wf_sat, T_wf_h]
        self.T_sf_plot = [T_sf_c, T_sf_h]

    def solve(self):
        self.check_calculable()
        if self.params['HX_type'] == 'evaporator':
            print(self.su_wf.T)
            lb = PropsSI('P', 'T', self.su_wf.T, 'Q', 0, self.su_wf.fluid)
            ub = PropsSI('P', 'T', self.su_sf.T-self.params['Delta_T_sup_or_sub']-self.params['min_pinch'], 'Q', 0, self.su_wf.fluid)

        if self.params['HX_type'] == 'condenser':
            lb = PropsSI('P', 'T', self.su_sf.T+self.params['Delta_T_sup_or_sub']+self.params['min_pinch'], 'Q', 0, self.su_wf.fluid)
            try:
                ub = PropsSI('P', 'T', self.su_wf.T, 'Q', 0, self.su_wf.fluid)
            except:
                ub = PropsSI('P_max', 'T', self.su_wf.T, 'Q', 0, self.su_wf.fluid)
            "We make the assumption that the pinch point is at the exit of the condenser BUT it's not always the case! BE CAREFUL"
        
        self.flag = 0
        # Minimization problem
        if abs(ub - lb) < 10 or (ub < lb):
            print('The temperature of the working fluid is too high or too low to reach the secondary fluid temperature.')
            P_wf_f = (lb + ub)/2
            self.flag = 0
        else:
            P_wf0 = np.linspace(lb, ub, 20) # Initial guess
            res_vector = []
            for P_wf in P_wf0:
                try:
                    results_temp = self.MB_model(P_wf)
                    res = self.res
                    if math.isnan(res):
                        res = 1e6
                    else:
                        res = self.res
                        # Append res to the list

                    res_vector.append(res)
                    
                except:
                    res_vector.append(1e6)

                min_indx = np.argmin(res_vector)
                self.res = res_vector[min_indx]
                P_wf_f = P_wf0[min_indx]

        results = self.MB_model(P_wf_f)

        self.ex_wf.set_fluid(self.su_wf.fluid)
        self.ex_wf.set_T(self.T_wf_out)
        self.ex_wf.set_p(P_wf_f)
        self.ex_wf.set_m_dot(self.su_wf.m_dot)
        self.ex_wf.set_fluid(self.su_wf.fluid)
        # self.Q_dot.set_Q_dot(self.Q)

        self.ex_sf.set_p(self.su_sf.p)
        self.ex_sf.set_m_dot(self.m_dot_sf)
        self.ex_sf.set_fluid(self.ex_sf.fluid)

    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q}")
        # print(f"Q_sub: {self.Q_sub}")
        # print(f"Q_sat: {self.Q_sat}")
        # print(f"Q_sup: {self.Q_sup}")
        print(f"m_dot_sf: {self.m_dot_sf}")
        print(f"T_wf_out: {self.T_wf_out}")
        print(f"U_sub: {self.U_sub}")
        print(f"U_sat: {self.U_sat}")
        print(f"U_sup: {self.U_sup}")
        print(f"tau: {self.tau}")
        # print(f"A_sub: {self.A_sub}")
        # print(f"A_sat: {self.A_sat}")
        # print(f"A_sup: {self.A_sup}")
        print(f"res: {self.res}")
        print(f"flag_tau: {self.flag_tau}")
        print("======================")





