"""
author: Elise Neven
email: elise.neven@uliege.be

"""


from component.base_component import BaseComponent

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np
import time

class HXSimpleMB(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su_hot = MassConnector()
        self.su_cold = MassConnector()
        self.ex_hot = MassConnector()
        self.su_cold = MassConnector()
        self.Q_dot = HeatConnector()

    def get_required_inputs(self):
        
        if self.inputs == {}:
            # Hot Fluid
            if self.su['H'].T is not None:
                self.inputs['Hsu_T'] = self.su['H'].T
            elif self.su['H'].h is not None:
                self.inputs['Hsu_h'] = self.su['H'].h
            if self.su['H'].p is not None:
                                self.inputs['Hsu_p'] = self.su['H'].p
            if self.su['H'].fluid is not None:
                self.inputs['Hsu_fluid'] = self.su['H'].fluid
            if self.su['H'].m_dot is not None:
                self.inputs['Hsu_m_dot'] = self.su['H'].m_dot
                
            # Cold Fluid                
            if self.su['C'].T is not None:
                self.inputs['Csu_T'] = self.su['C'].T
            elif self.su['C'].h is not None:
                self.inputs['Csu_h'] = self.su['C'].h
            if self.su['C'].p is not None:
                self.inputs['Csu_p'] = self.su['C'].p
            if self.su['C'].fluid is not None:
                self.inputs['Csu_fluid'] = self.su['C'].fluid
            if self.su['C'].m_dot is not None:
                self.inputs['Csu_m_dot'] = self.su['C'].m_dot
                
        if self.inputs != {}:
            # Hot Fluid
            self.su['H'].set_fluid(self.inputs['Hsu_fluid'])
            if 'Hsu_T' in self.inputs:
                self.su['H'].set_T(self.inputs['Hsu_T'])
            elif 'Hsu_h' in self.inputs:
                self.su['H'].set_h(self.inputs['Hsu_h'])
            if 'Hsu_p' in self.inputs:
                self.su['H'].set_p(self.inputs['Hsu_p'])
            if 'Hsu_m_dot' in self.inputs:
                self.su['H'].set_m_dot(self.inputs['Hsu_m_dot'])

            # Cold Fluid
            self.su['C'].set_fluid(self.inputs['Csu_fluid'])
            if 'Csu_T' in self.inputs:
                self.su['C'].set_T(self.inputs['Csu_T'])
            elif 'Csu_h' in self.inputs:
                self.su['C'].set_h(self.inputs['Csu_h'])
            if 'Csu_p' in self.inputs:
                self.su['C'].set_p(self.inputs['Csu_p'])
            if 'Csu_m_dot' in self.inputs:
                self.su['C'].set_m_dot(self.inputs['Csu_m_dot'])

        return ['Hsu_p', 'Hsu_T', 'Hsu_m_dot', 'Hsu_fluid', 'Csu_p', 'Csu_T', 'Csu_m_dot', 'Csu_fluid']

    def get_required_parameters(self):
        return ['A_htx', 'L_HTX', 'V_HTX', 'Flow_Type',
                'A_canal_h', 'A_canal_c', 'D_h',
                'k_plate', 't_plate', 'n_plates',
                'co_pitch', 'chevron_angle', 'fouling']
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - H_su: fluid={self.su['H'].fluid}, T={self.su['H'].T}, p={self.su['H'].p}, m_dot={self.su['H'].m_dot}")
        print(f"  - C_su: fluid={self.su['C'].fluid}, T={self.su['C'].T}, p={self.su['C'].p}, m_dot={self.su['C'].m_dot}")

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
    
    def e_NTU(self, NTU, C_r):
        
        if self.params['Flow_Type'] == "CounterFlow":
            eps = (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
                
        return eps

#%%
    
    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:
            
            # Detect Phase change
            # self.detect_phase_change()
            
            # Calcul de C_r
            cp_h = PropsSI('C', 'H', self.su['H'].h, 'P', self.su['H'].p, self.su['H'].fluid)
            cp_c = PropsSI('C', 'H', self.su['C'].h, 'P', self.su['C'].p, self.su['C'].fluid)
            
            C_h = cp_h*self.su['H'].m_dot
            C_c = cp_c*self.su['C'].m_dot
            
            C_min = min(C_h, C_c)
            C_max = max(C_h, C_c)
            C_r = C_min/C_max
                        
            # Calcul de NTU
            T_w = (self.su['H'].T + self.su['C'].T)/2
            
            mu_h, Pr_h, k_h = PropsSI(('V','PRANDTL','L'), 'H', self.su['H'].h, 'P', self.su['H'].p, self.su['H'].fluid)
            mu_c, Pr_c, k_c = PropsSI(('V','PRANDTL','L'), 'H', self.su['C'].h, 'P', self.su['C'].p, self.su['C'].fluid)
            
            G_h = self.su['H'].m_dot/self.params['A_canal_h']
            G_c = self.su['C'].m_dot/self.params['A_canal_c']
            
            h_h = Gnielinski_Pipe_HTC(mu_h, Pr_h, Pr_h, k_h, G_h, self.params['D_h'], self.params['L_HTX'])
            h_c = Gnielinski_Pipe_HTC(mu_c, Pr_c, Pr_c, k_c, G_c, self.params['D_h'], self.params['L_HTX'])
                        
            AU = (1/(self.params['A_htx']*h_h) + 1/(self.params['A_htx']*h_c) + self.params['t_plate']/(self.params['k_plate']*self.params['A_htx']) + self.params['fouling']/self.params['A_htx'])**(-1)         

            NTU = AU/C_min
                        
            # epsilon NTU 
            eps = self.e_NTU(NTU, C_r)
                        
            # Calcul de Q
            h_c_Th = PropsSI('H','T',self.su['H'].T,'P',self.su['C'].p,self.su['C'].fluid)
            h_h_Tc = PropsSI('H','T',self.su['C'].T,'P',self.su['H'].p,self.su['H'].fluid)
            
            DH_pc_c = PropsSI('H','Q',1,'P',self.su['C'].p,self.su['C'].fluid) - PropsSI('H','Q',0,'P',self.su['C'].p,self.su['C'].fluid)
            DH_pc_h = PropsSI('H','Q',1,'P',self.su['H'].p,self.su['H'].fluid) - PropsSI('H','Q',0,'P',self.su['H'].p,self.su['H'].fluid)
            
            Qmax_c = self.su['C'].m_dot*((h_c_Th - self.su['C'].h))
            Qmax_h = self.su['H'].m_dot*((self.su['H'].h - h_h_Tc))
                        
            Qmax = min(Qmax_c, Qmax_h)
            
            Q = eps*Qmax

            # Exhaust conditions 
            self.ex['H'].set_properties(H = self.su['H'].h - Q/self.su['H'].m_dot, fluid = self.su['H'].fluid, m_dot = self.su['H'].m_dot, P = self.su['H'].p)
            self.ex['C'].set_properties(H = self.su['C'].h + Q/self.su['C'].m_dot, fluid = self.su['C'].fluid, m_dot = self.su['C'].m_dot, P = self.su['C'].p)
                        
            self.defined = True
        
        else:
            if not self.calculable:
                print("Input of the component not completely known. Required inputs:")
                for input in self.get_required_inputs():
                    if input not in self.inputs:
                        print(f"  - {input}")
                        
            if not self.parametrized:
                print("Parameters of the component not completely known. Required parameters:")
                for param in self.get_required_parameters():
                    if param not in self.params:
                        print(f"  - {param}")
    
    def print_results(self):
        if self.defined:
            print("=== Heat Exchanger Results ===")
            print(f"  - H_ex: fluid={self.ex['H'].fluid}, T={self.ex['H'].T}, p={self.ex['H'].p}, m_dot={self.ex['H'].m_dot}")
            print(f"  - C_ex: fluid={self.ex['C'].fluid}, T={self.ex['C'].T}, p={self.ex['C'].p}, m_dot={self.ex['C'].m_dot}")
            print(f"  - Q_dot: temperature_in={self.Q_dot}")

        else:
            print("Heat Exchanger component is not defined. Ensure it is solved first.")

#-------------------------------------------------------------
if __name__ == '__main__':
    # Example usage
    HX = HXeNTU()
    # expander.print_setup()
    # Connectors (if not already set elsewhere)
    # expander.su = Mass_connector()
    # expander.ex = Mass_connector()

    # # Set fluid states for connectors
    # expander.su.set_fluid('R134a')
    # expander.ex.set_fluid('R134a')

    # # # Set properties for su connector
    # expander.su.set_p(1e5)
    # expander.su.set_T(300)  # You need to set su.h appropriately

    # # # Set properties for ex connector
    # expander.ex.set_p(1e4)
    
    # # Set rotational speed
    # expander.work_exp.set_speed(3000)

    # # Set ambient temperature
    # expander.heat_amb.set_temperature_in(298.15)

    # Setting inputs
    HX.set_inputs(
        # First fluid
        Hsu_fluid = 'Cyclopentane',
        Hsu_T = 205 + 273.15, # K
        Hsu_p = 1*1e5, # Pa
        Hsu_m_dot = 0.014, # kg/s

        # Second fluid
        Csu_fluid = 'Water',
        Csu_T = 12 + 273.15, # K
        Csu_p = 4*1e5, # Pa
        Csu_m_dot = 0.08, # kg/s  # Make sure to include fluid information
    )
    
    # #Either set as connectors or set as inputs

    # Setting parameters
    
    # Heat transfer area  
    A_htx = 0.752 # m^2
    
    # HTX dimensions
    L = 0.393 # [m] : length
    w = 0.243 # [m] : width
    h = 0.0446 # [m] : height
    l_v = 0.324 # [m] : length between ports
    casing_t = 0.005 # [m] : casing thickness # !!! arbitrary
    
    # Number and thickness of plates
    n_plates = 10
    t_plates = 0.0008 # [m] # !!! arbitrary
    
    # Fooling factor
    fooling = 98.51/1000 # (m^2*K)/W

    # Number of canals
    C_n_canals = 5
    H_n_canals = 4
    
    # Plate values 
    plate_cond = 45 # [W/(m*K)] : plate conduction
    plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
    chevron_angle = 20*np.pi/180 # !!! arbitrary

    # Total volume of each part
    V_tot = 0.9*1e-3 # [m^3]
    
    # Canal thickness
    C_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*C_n_canals)
    H_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*H_n_canals)

    # Total Canal Surface
    C_CS = C_canal_t*(w-2*casing_t)*C_n_canals
    H_CS = H_canal_t*(w-2*casing_t)*H_n_canals

    # Dh : hydraulic diameter
    C_Dh = (4*C_canal_t*w)/(2*C_canal_t+2*w)
    H_Dh = (4*H_canal_t*w)/(2*H_canal_t+2*w)
    
    HX.set_parameters(
        A_htx=0.752, L_HTX=0.393, V_HTX=0.9*1e-3, Flow_Type = 'CounterFlow',
        A_canal_h=H_CS, A_canal_c=C_CS, D_h=H_Dh, 
        k_plate=45, t_plate=0.0008, n_plates = 10,
        co_pitch=0.005, chevron_angle=20*np.pi/180, fouling=98.51/1000
    )

    # Solve the expander component
    HX.solve()
    HX.print_setup()
    HX.print_results()
    
    # print(expander.defined)  # Should print True if the component was successfully solved

