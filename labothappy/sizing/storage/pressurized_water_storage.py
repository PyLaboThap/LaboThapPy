# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 13:52:57 2025

@author: Basile
"""

import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize

#%%

class PressuzedWaterStorageDesign(object):

    def __init__(self, fluid):
        # Inputs
        self.inputs = {}
        
        # Params
        self.params = {}

        # Abstract State
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Gravity Constant
        self.g = 9.81 # m/s^2
        
        return

    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value

    # ---------------- Computation Related Methods ----------------------------------------------------------------------

    def req_thickness(self, P_des, R_i, sigma_max):
        """
        Inputs:
        -------------
        P_des     : Design pressure [Pa]
        sigma_max : Maximum Allowable tensile stress of the material [Pa]
        R_i         : For pipe, the inside radius R is determined by the nominal outside radius minus the nominal wall thickness [m]
        
        Reference 
        -------------
        BVPC 2007 
        """
        
        QF = 0.8 # Quality Factor oof values not to exceed as a security margin
        
        E_joint = 0.85 # Joint efficiency : Butt joints welded from both sides, or from one side only without a permanent
                       # backing strip, that are verified as having achieved full penetration and fusion as
                       # required by UW-35. Case b 
                    
        t_req_circ = (P_des*R_i)/(sigma_max*E_joint - 0.6*P_des) # Circumferential Stress (Longitudinal Joints)
        t_req_long = (P_des*R_i)/(2*sigma_max*E_joint + 0.4*P_des) # Longitudinal Stress (Circumferential Joints)
        
        t_head = (P_des*R_i)/(2*sigma_max*E_joint - 0.2*P_des) # Spherical Stress (Spherical Shell)
        
        print(f"t_req_P : {max(t_req_circ, t_req_long, t_head)}")
        print(f"t_req_min : {self.t_req_min}")
        
        return max(max(t_req_circ, t_req_long, t_head)/QF, self.t_req_min)
    
    def material_properties(self):
        
        if self.params['Material'] == 'Cast_Iron':
            self.rho_mat = 7300 # kg/m^3
            self.sigma_max = 80e6 # Pa     
            self.t_req_min = 10*1e-3 # m
            
        elif self.params['Material'] == 'Carbon_Steel':
            self.rho_mat = 7850 # kg/m^3
            self.sigma_max = 140e6 # Pa 
            self.t_req_min = 0.6*1e-3
        
        else: 
            raise ValueError(f"Unknown material: {self.params['Material']}")
        
        return    
    
    # ---------------- Design Methods ----------------------------------------------------------------------
    
    def design_system(self, x):
        
        self.material_properties()
        
        # self.T_des = x[0]
        self.D_i   = x[0]
        self.H     = x[1]
    
        self.T_des = self.inputs['T_des']
    
        if self.T_des >= 100+273.15:
            (P_sat, rho_l) = PropsSI(('P','D'), 'T', self.T_des, 'Q', 0, self.fluid)
        else:
            P_sat = 101325
            rho_l = PropsSI('D', 'T', self.T_des, 'P', P_sat, self.fluid)
    
        P_stat = P_sat + rho_l*self.g*self.H
        self.P_des = P_stat*(1+self.params['P_OD'])
    
        self.t_req = self.req_thickness(self.P_des, self.D_i/2, self.sigma_max)
    
        self.lateral_plate_vol = np.pi*(self.D_i/2 + self.t_req)**2 * self.t_req * 2
        self.shell_vol = (np.pi*(self.D_i/2 + self.t_req)**2 - np.pi*(self.D_i/2)**2)*self.H
        self.metal_vol = self.shell_vol + self.lateral_plate_vol
        
        self.m_mat = self.rho_mat*self.metal_vol  
    
        self.fluid_vol = np.pi*(self.D_i/2)**2 * self.H
        h_cap = (PropsSI('H', 'T', self.T_des, 'P', self.P_des, self.fluid) - PropsSI('H', 'T', self.inputs['T_in'], 'P', self.P_des, self.fluid))
        E_cap = h_cap*self.fluid_vol*rho_l
        
        self.E_cap_MWh = E_cap/(3600*1e6)
        
        self.Cost_cap_dollar = 4500*self.fluid_vol**0.6 + 10000
        
        return self.m_mat/self.E_cap_MWh # self.Cost_cap_dollar/self.E_cap_MWh

    def opt_design(self):

        # Bounds: [T_des, D_i, H]
        bounds = [
            # tuple(self.params['bounds_T_des']),
            tuple(self.params['bounds_D_i']),  # inner diameter in meters
            tuple(self.params['bounds_H']),
        ]
        
        # Objective function to minimize
        def objective(x):
            try:
                return self.design_system(x)
            except Exception as e:
                print(f"Error at x={x}: {e}")
                return np.inf
        
        # x0 = [self.params['bounds_T_des'][-1], self.params['bounds_D_i'][-1], self.params['bounds_H'][-1]]
        x0 = [self.params['bounds_D_i'][-1], self.params['bounds_H'][-1]]
        
        result = minimize(
            fun=objective,
            x0=x0,
            method='L-BFGS-B',
            bounds=bounds,
            options={'disp': False}
        )
        
        return result

#%% 

if __name__ == "__main__":
    
    T_des_vec = np.linspace(100,150,21) + 273.15
    
    # Storage
    fluid = 'Water'
    fun_vec = []
    cap_vec = []
    
    for T_des in T_des_vec:
        Tank = PressuzedWaterStorageDesign(fluid)
        
        T_in = 15 + 273.15 # K
        P_OD = 0.2
        
        Tank.set_inputs(
            T_in = T_in, # [K]
            T_des = T_des # [K]
            )
        
        Tank.set_parameters(
            P_OD = P_OD, # [-]
            Material = 'Carbon_Steel', # [-]
            
            bounds_T_des = [80 + 273.15, 150 + 273.15], # K
            bounds_D_i = [1, 100], # m
            bounds_H = [0.5, 100], # m
            )
        
        opt_result = Tank.opt_design()
        
        fun_vec.append(opt_result.fun)
        cap_vec.append(Tank.E_cap_MWh*1e3/Tank.fluid_vol)
    
    plt.figure()
    plt.plot(T_des_vec, fun_vec)
    plt.show()
    
    plt.figure()
    plt.plot(T_des_vec, cap_vec)
    plt.show()

    # if opt_result.success:
    #     print("\n✅ Optimisation réussie !")
    #     print(f"T_des optimal : {opt_result.x[0] - 273.15:.1f} °C")
    #     print(f"D_i optimal   : {opt_result.x[1]:.3f} m")
    #     print(f"H optimal     : {opt_result.x[2]:.3f} m")
    #     print(f"Coût [kg/MWh] : {opt_result.fun:.2f} kg/MWh")
    # else:
    #     print("\n❌ Optimisation échouée.")
    #     print(opt_result.message)

# print(f"E_cap [J] : {E_cap}")
# print(f"E_cap [Wh] : {E_cap/3600}")
# print(f"E_cap [kWh] : {E_cap/(3600*1000)}")
# print(f"E_cap [MWh] : {E_cap/(3600*1000*1000)}")

# print(f"mass_mat [kg] : {mass_mat}")
# print(f"mass_mat [t] : {mass_mat/1000}")
