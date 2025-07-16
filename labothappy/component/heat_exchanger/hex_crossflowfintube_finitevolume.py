# -*- coding: utf-8 -*-
"""
Created on Wed May  8 09:38:42 2024

@author: Basile & Sam
"""
import __init__

# External Toolbox 
import numpy as np
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP

# Connectors
from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

# HTC correlations
from correlations.convection.fins_htc import htc_tube_and_fins
from correlations.convection.pipe_htc import gnielinski_pipe_htc
from correlations.convection.pipe_htc import horizontal_tube_internal_condensation
from correlations.convection.pipe_htc import horizontal_flow_boiling

#Pressure drop correlations
from correlations.pressure_drop.fins_DP import DP_tube_and_fins
from correlations.pressure_drop.pipe_DP import gnielinski_pipe_DP, Muller_Steinhagen_Heck_DP


# Phase related import
from correlations.properties.void_fraction import void_fraction

# Component base frame
from component.base_component import BaseComponent

class CrossFlowTubeAndFinsHTX(BaseComponent):
    """
    Component: Cross-Flow Tube-and-Fins Heat Exchanger (HTX)
    
    Model: Discretized Cross-Flow Heat Exchanger with Fin Efficiency and Pressure Drop
    
    **Description**:
    
        This model simulates a cross-flow heat exchanger with tube-and-fins geometry, suitable for systems where one fluid flows through tubes and the other over an external finned surface. The model uses detailed correlations for heat transfer coefficients and pressure drops, accounting for phase change phenomena (condensation and evaporation). It is designed for steady-state, on-design simulations.
        The model uses a finite volume discretization to take into consideration the fluid properties evolution accross the heat exchanger.
    
    **Assumptions**:
    
        - Steady-state operation.
        - Uniform flow distribution across tubes and fins.
        - Correlations used for heat transfer and pressure drop calculations (e.g., Gnielinski, Müller-Steinhagen-Heck).
        - Two-phase flow effects are considered where applicable.
        - Fluid properties are retrieved from CoolProp.
    
    **Connectors**:
    
        su_H (MassConnector): Mass connector for the hot-side supply (tube or fin side depending on configuration).
        su_C (MassConnector): Mass connector for the cold-side supply (tube or fin side depending on configuration).
        
        ex_H (MassConnector): Mass connector for the hot-side exhaust.
        ex_C (MassConnector): Mass connector for the cold-side exhaust.
        
        Q_dot (HeatConnector): Heat connection (not actively used in current implementation).
    
    **Parameters**:
    
        H_DP_ON (bool): Flag to enable pressure drop calculation on hot side.
        C_DP_ON (bool): Flag to enable pressure drop calculation on cold side.
        n_disc (int): Number of discretization segments along the heat exchanger length.
        Fin_Side (str): Specifies which side has fins ('H' for hot side, 'C' for cold side).
        
        Geometry Parameters:
            
            Fin_OD : Fin diameter/length [m]
            
            Fin_per_m : Fin density per m [1/m]
            
            Fin_t : Fin Thickness : [m]
            
            Fin_type : 'Annular' or 'Square'
            
            fouling : Global fouling resistance [K/W] 
            
            h : Bundle height [m]
            
            k_fin : Fin thermal conductivity [W/(m*K)]
            
            pitch_V : Vertical Pitch [m]
            
            pitch_H : Horizontal Pitch [m]
            
            tube_arrang : Tube arrangement 'Inline' or 'Staggered'
            
            Tube_cond : Tube material conductivity [W/(m*K)]
            
            Tube_L : Tube length [m]
            
            Tube_OD : Tube outer diameter [m]
            
            Tube_t : Tube thickness [m]
            
            w : Bundle width [m]
            
            Fin_Side : Fin side 'H' or 'C'
    
    **Inputs**:
    
        fluid_H (str): Hot-side fluid.
        
        h_su_H (float): Hot-side inlet specific enthalpy [J/kg].
        
        P_su_H (float): Hot-side inlet pressure [Pa].
        
        m_dot_H (float): Hot-side mass flow rate [kg/s].
    
        fluid_C (str): Cold-side fluid.
        
        h_su_C (float): Cold-side inlet specific enthalpy [J/kg].
        
        P_su_C (float): Cold-side inlet pressure [Pa].
        
        m_dot_C (float): Cold-side mass flow rate [kg/s].
    
    **Outputs**:

        h_ex_C: Cold exhaust side enthalpy. [J/kg]

        P_ex_C: Cold exhaust side pressure. [Pa]

        h_ex_H: Hot exhaust side enthalpy. [J/kg]

        P_ex_H: Hot exhaust side pressure. [Pa]

        Q_dot: Heat Exchanger's heat duty. [W]
        
        M_tube: Total mass inventory in the tube-side [kg].
        
        M_bank: Total mass inventory in the fin-side (bundle side) [kg].
        
        T_matrix, H_matrix, P_matrix, x_matrix: Discretized distributions of temperature, enthalpy, pressure, and vapor quality within the heat exchanger.
    
    """

    def __init__(self):
        """
        su : Supply - 'H' : hot
                    - 'C' : cold

        ex : Exhaust - 'H' : hot
                     - 'C' : cold
                    
        Q_dot : Heat connection to the ambient
        
        HTX_Type : Type of HTX - Plate 
                               - Shell and Tube
                               - Tube and Fins
        """
        
        super().__init__()
        
        self.su_H = MassConnector()
        self.su_C = MassConnector()
        
        self.ex_H = MassConnector()
        self.ex_C = MassConnector() # Mass_connector
        
        self.Q_dot = HeatConnector()
        self.debug = 0

    #%%    
    
    def get_required_inputs(self): # Used in check_calculable to see if all of the required inputs are set
        """
        Hot side required inputs :
            
            - T_su_H or h_su_H : Hot supply temperature or enthalpy
            - P_su_H          : Hot supply pressure
            - Fluid_H     : Hot supply fluid
            - m_dot_H      : Hot supply flow rate
            
        Cold side required inputs : 
            
            - T_su_C or h_su_C  : Cold supply temperature or enthalpy
            - P_su_C          : Cold supply pressure
            - Fluid_C      : Cold supply fluid
            - m_dot_C      : Cold supply flow rate
        """
        self.sync_inputs()
        # Return a list of required inputs
        return['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']
    
#%%

    def get_required_parameters(self):
        """
        General Parameters : 
            
            - htc_type  : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - H_DP_ON   : Hot side pressure drop considered or not
            - C_DP_ON   : Cold side pressure drop considered or not
            - n_disc    : number of discretizations
        
        Geometry Parameters depend on specific geometry python files.
            
        """

        general_parameters = ['H_DP_ON', 'C_DP_ON','n_disc']
            
        geometry_parameters = ['Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type', 'fouling', 
                               'h', 'k_fin', 'pitch_V', 'pitch_H', 'tube_arrang',  
                               'Tube_cond', 'Tube_L', 'Tube_OD', 'Tube_t', 'w', 'Fin_Side']
        
        return general_parameters + geometry_parameters

#%%
            
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - H_su: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - C_su: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")

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
    def set_debug(self):
        self.debug = 1
        return

    def compute_detailed_geom(self):
        
        # Computation - n_tubes - n_row
        n_row_1 = np.floor(self.params['h']/self.params['pitch_V'])
        n_col = np.floor(self.params['w']/self.params['pitch_H'])

        if self.params['tube_arrang'] == 'Inline':
            n_row = n_row_1
            self.params['n_tubes'] = n_col*n_row
        elif self.params['tube_arrang'] == 'Staggered':
            if np.mod(int(n_col),2) == 0:
                n_row = (2*n_row_1 - 1)/2
                self.params['n_tubes'] = n_col*n_row
            else:
                n_row = (np.ceil(n_col/2)*n_row_1 + np.floor(n_col/2)*(n_row_1-1))/n_col
                self.params['n_tubes'] = np.floor(n_col*n_row)   
        else:
            print("Tube Arrangement is not 'Inline' or 'Staggered'.")
          
        # n_rows
        self.params['n_rows'] = int(np.ceil(n_row))
        
        # A_flow
        self.params['A_flow'] = self.params['Tube_L']*self.params['h']
        
        # Heat Transfer Area
        N_fins = self.params['Tube_L']*self.params['Fin_per_m'] - 1
        self.params['A_unfinned'] = 2*np.pi*(self.params['Tube_OD']/2)*self.params['n_tubes']*self.params['Tube_L']

        A_out_tube = 2*np.pi*(self.params['Tube_OD']/2)*self.params['n_tubes']*(self.params['Tube_L']- N_fins*self.params['Fin_t'])
        A_out_fin = 2*np.pi*(self.params['Fin_OD']/2)*self.params['Fin_t']*N_fins*self.params['n_tubes']
        A_out_plate_fin = 2*np.pi*((self.params['Fin_OD']/2)**2 - (self.params['Tube_OD']/2)**2)*N_fins*self.params['n_tubes']

        self.params['A_finned'] = A_out_tube + A_out_fin + A_out_plate_fin
        
        self.params['A_in_tot'] = self.params['Tube_L']*self.params['n_tubes']*2*np.pi*(self.params['Tube_OD']/2 - self.params['Tube_t'])
        
        # B_V_tot and T_V_tot computation
        
        self.params['B_V_tot'] = self.params['Tube_L']*self.params['w']*self.params['h']
        self.params['T_V_tot'] = self.params['A_in_tot']*self.params['Tube_L']
                
        return

    def compute_cell(self, T_b_in, p_b_in, h_b_in, m_dot_b_in_all, T_t_in, p_t_in, h_t_in, m_dot_1_tube_in,j):    
              
        """
        T_b_in, p_b_in, h_b_in: Inlet state of the bundle side (fin side).
        T_t_in, p_t_in, h_t_in: Inlet state of the tube side.
        m_dot_b_in_all: Total bundle-side mass flow.
        m_dot_1_tube_in: Flow rate in one tube.
        j: index for discretization segment.
        """

        if self.debug:
            print("-----------------")
            print("h_b_in",h_b_in)
            print("h_t_in",h_t_in)
            print("p_b_in",p_b_in)
            print("p_t_in",p_t_in)
            print("T_b_in",T_b_in)
            print("T_t_in",T_t_in)
            print("-----")
        
        "1) Heat Transfer Coefficients"
    
        T_wall = (T_b_in + T_t_in)/2
    
        # Initialize AbstractStates
        AS_b = CP.AbstractState("HEOS", self.B_su.fluid)
        AS_b.update(CP.HmassP_INPUTS, h_b_in, p_b_in)
    
        AS_t = CP.AbstractState("HEOS", self.T_su.fluid)
        AS_t.update(CP.HmassP_INPUTS, h_t_in, p_t_in)
        
        # Tube Bank
        self.alpha_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1])
        alpha_b = htc_tube_and_fins(self.B_su.fluid, self.params, p_b_in, h_b_in, m_dot_b_in_all, self.params['Fin_type'])[0]
        
        A_in_one_tube = (np.pi/4)*(self.params['Tube_OD'] - 2*self.params['Tube_t'])**2
        P_in_one_tube = (np.pi)*(self.params['Tube_OD'] - 2*self.params['Tube_t'])
        D_h_one_tube = 4*A_in_one_tube/P_in_one_tube
        G_1t = m_dot_1_tube_in/A_in_one_tube
    
        x = AS_t.Q()
        
        if x < 0: # 1 phase case
            # Tube
            mu = AS_t.viscosity()
            Pr = AS_t.Prandtl()
            k = AS_t.conductivity()
            rho = AS_t.rhomass()
    
            AS_t.update(CP.PT_INPUTS, p_t_in, T_wall)
            Pr_w = AS_t.Prandtl()     
    
            if self.debug:
                print("Pr",Pr)
                print("Pr_w",Pr_w)
                print("ratio",(Pr/Pr_w)**0.11)
                        
            alpha_t = gnielinski_pipe_htc(mu, Pr, Pr_w, k, G_1t, self.params['Tube_OD'] - self.params['Tube_t'], self.params['Tube_L']/self.params['n_disc'])[0]
            DP_t = gnielinski_pipe_DP(mu, rho, G_1t, D_h_one_tube, self.params['Tube_L']/self.params['n_disc'])
            
        else: # 2 phase flow
            AS_t.update(CP.PQ_INPUTS, p_t_in, 0)
            P_sat = p_t_in
            DP_t = Muller_Steinhagen_Heck_DP(self.T_su.fluid, G_1t, P_sat, x, self.params['Tube_L']/self.params['n_disc'], D_h_one_tube)           
            if T_wall <= T_t_in: # Condensation
                alpha_t = horizontal_tube_internal_condensation(self.T_su.fluid,self.T_su.m_dot,P_sat,h_t_in,T_wall,self.params['Tube_OD'] - self.params['Tube_t'])
            else: # Evaporation
                x_t = x
                def equation(q):
                    alpha_t_local = horizontal_flow_boiling(self.T_su.fluid, G_1t, P_sat, x_t, self.params['Tube_OD'] - self.params['Tube_t'], q)
                    AU = 1 / (1 / (alpha_t_local * A_in_one_tube) + 1 / (alpha_b * self.params['A_finned']/(self.params['n_tubes']*self.params['n_disc'])))
                    C_b = m_dot_b_in_all * AS_b.cpmass()
                    C_t = 20000
                    C_min = min(C_b, C_t)
                    C_max = max(C_b, C_t)
                    NTU = max(0, AU / C_min)
                    eps = 1 - np.exp((1 / (C_min/C_max)) * NTU**0.22 * (np.exp(-(C_min/C_max) * NTU**0.78) - 1))
                    Q_dot_max = C_min * abs(T_b_in - T_t_in)
                    Q_dot_1_tube = Q_dot_max * eps
                    return q - Q_dot_1_tube / (self.params['A_in_tot']/(self.params['n_tubes']*self.params['n_disc']))
    
                q_solution = fsolve(equation, 1.0)[0]
    
                # Use q_solution in place of q
                alpha_t = horizontal_flow_boiling(self.T_su.fluid, G_1t, P_sat, x_t, self.params['Tube_OD'] - self.params['Tube_t'], q_solution)
                
        A_out_1_tube = self.params['A_finned']/(self.params['n_tubes']*self.params['n_disc']) 
        A_in_1_tube = self.params['A_unfinned']/(self.params['n_tubes']*self.params['n_disc']) 
        
        AU = 1/(1/(alpha_t*A_in_1_tube) + 1/(alpha_b*A_out_1_tube))
        
        if self.debug == 1:
            print("alpha_b",alpha_b)
            print("alpha_t",alpha_t)
            print("-----")
            
        "2) Q_dot_cell : e-NTU method"
        
        m_dot_b_in = m_dot_b_in_all/(self.params['n_disc']*(self.params['n_tubes']/self.params['n_rows']))
        m_dot_t_in = m_dot_1_tube_in
        
        C_b = m_dot_b_in * AS_b.cpmass()
        
        if x < 0: # 1 phase case
            C_t = m_dot_t_in * AS_t.cpmass()
        else: # 2 phase case
            C_t = 20000
        
        C_min = min(C_b, C_t)
        C_max = max(C_b, C_t)
            
        C_r = C_min/C_max
            
        NTU = max(0,AU/C_min)
            
        # CrossFlow, both unmixed
        eps = 1 - np.exp((1/C_r)*NTU**0.22*(np.exp(-C_r*NTU**0.78)-1))
            
        Q_dot_max = C_min*(abs(T_b_in - T_t_in))
        Q_dot_1_tube = Q_dot_max*eps
        
        Q_dot_1_row = Q_dot_1_tube*(self.params['n_tubes']/self.params['n_rows'])
    
        if self.debug == 1:
            print("AU", AU)
            print("C_r", C_r)
            print("NTU", NTU)
            print("eps", eps)
            print("Q_dot_max", Q_dot_max)
            print("Q_dot_1_tube", Q_dot_1_tube)
            print("Q_dot_1_row", Q_dot_1_row)
        
        "3) Outlet Conditions"
        
        if T_b_in < T_t_in: 
            h_b_out = (h_b_in*m_dot_b_in + Q_dot_1_tube)/m_dot_b_in
            h_t_out = (h_t_in*m_dot_t_in - Q_dot_1_tube)/m_dot_t_in
        else:
            h_t_out = (h_t_in*m_dot_t_in + Q_dot_1_tube)/m_dot_t_in
            h_b_out = (h_b_in*m_dot_b_in - Q_dot_1_tube)/m_dot_b_in  
        
        DP_b = DP_tube_and_fins(self.B_su.fluid, self.params, p_b_in, T_b_in, m_dot_b_in)
        p_b_out = p_b_in-DP_b
        p_t_out = p_t_in-DP_t
    
        AS_b.update(CP.HmassP_INPUTS, h_b_out, p_b_out)
        AS_t.update(CP.HmassP_INPUTS, h_t_out, p_t_out)
    
        T_b_out = AS_b.T()
        T_t_out = AS_t.T()
    
        if self.debug == 1:
            print("-----------------")
            print("h_b_out", h_b_out)
            print("h_t_out", h_t_out)
            print("p_b_out", p_b_out)
            print("p_t_out", p_t_out)
            print("T_b_out", T_b_out)
            print("T_t_out", T_t_out)
            print("---------")
            
        return h_b_out, h_t_out, p_b_out, p_t_out, T_b_out, T_t_out, Q_dot_1_row
        
    def compute_mass(self, H_vec, P_vec, fluid, i):
        D_vec = np.zeros(len(H_vec))
        M_vec = np.zeros(len(H_vec))
        x_vec = np.zeros(len(H_vec))

        # Determine segment volume
        if np.mod(i, 2) == 0:
            Volume = self.params['B_V_tot'] / (self.params['n_disc'] * (self.params['n_rows'] + 1))
        else:
            N_tpr = self.params['n_tubes'] / self.params['n_rows']
            Tube_vol = self.params['Tube_L'] * np.pi * ((self.params['Tube_OD'] - self.params['Tube_t']) / 2) ** 2
            Volume = N_tpr * Tube_vol / self.params['n_disc']

        # Initialize AbstractState
        AS = CP.AbstractState("HEOS", fluid)

        for j in range(len(H_vec)):
            try:
                AS.update(CP.HmassP_INPUTS, H_vec[j], P_vec[j])
                x = AS.Q()
            except ValueError:
                x = -1  # Single-phase fallback if not two-phase

            x_vec[j] = x

            if x < 0 or x > 1:
                rho = AS.rhomass()
                D_vec[j] = rho
                M_vec[j] = rho * Volume
            else:
                AS.update(CP.PQ_INPUTS, P_vec[j], 1)
                rho_g = AS.rhomass()
                AS.update(CP.PQ_INPUTS, P_vec[j], 0)
                rho_l = AS.rhomass()

                eps_void, rho_avg = void_fraction(x, rho_g, rho_l)
                D_vec[j] = rho_avg
                M_vec[j] = rho_avg * Volume

        return [D_vec, M_vec, x_vec]


#%%

    def solve(self):
        self.check_calculable()
        self.check_parametrized()   
        self.compute_detailed_geom()

        if self.calculable and self.parametrized:        

            "------- 1) Create Finite Difference matrices ------------------------"
            
            self.T_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            self.P_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            self.H_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            self.D_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            self.M_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            self.x_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            self.Q_dot_matrix = np.zeros([2*self.params['n_rows'] + 1, self.params['n_disc'] + 1]) # self.create_matrix(self.n_disc, self.geom.n_rows)
            
            self.M_tube = 0
            self.M_bank = 0

            "------- 2) Determine Which Side is Fin Side ------------------------"

            if self.params['Fin_Side'] == 'C':
                self.B_su = self.su_C
                self.T_su = self.su_H
            elif self.params['Fin_Side'] == 'H':
                self.B_su = self.su_H
                self.T_su = self.su_C
            else:
                print("Fin_Side shall either be 'C' or 'H'")
            
            "------- 2) Initialize the input values ------------------------"
            
            # Bundle Fluid
            self.T_matrix[0].fill(self.B_su.T)
            self.P_matrix[0].fill(self.B_su.p)
            self.H_matrix[0].fill(self.B_su.h)
            self.Q_dot_matrix[0].fill(self.B_su.h)
            
            # Tube Fluid
            for i in range(len(self.T_matrix[:,0])):
                if np.mod(i+1,2) == 0:
                    self.T_matrix[i,0] = self.T_su.T
                    self.P_matrix[i,0] = self.T_su.p
                    self.H_matrix[i,0] = self.T_su.h
                    
            "------- 3) Compute Heat transfer over rows ------------------------"
            
            i = 0
            j = 0
            
            for i in range(len(self.T_matrix[:,0])-1):
                if np.mod(i,2) == 0:
                    for j in range(len(self.T_matrix[0])):
                        
                        if j < len(self.T_matrix[0]) - 1:
                            self.H_matrix[i+2,j], self.H_matrix[i+1,j+1], self.P_matrix[i+2,j], self.P_matrix[i+1,j+1], self.T_matrix[i+2,j], self.T_matrix[i+1,j+1], self.Q_dot_matrix[i][j] = self.compute_cell(self.T_matrix[i,j], self.P_matrix[i,j], self.H_matrix[i,j], self.B_su.m_dot, self.T_matrix[i+1,j], self.P_matrix[i+1,j], self.H_matrix[i+1,j], self.T_su.m_dot/self.params['n_tubes'],j)
                        else:
                            self.H_matrix[i+2,j], _ , self.P_matrix[i+2,j], _ , self.T_matrix[i+2,j], _ , self.Q_dot_matrix[i][j] = self.compute_cell(self.T_matrix[i,j], self.P_matrix[i,j], self.H_matrix[i,j], self.B_su.m_dot, self.T_matrix[i+1,j], self.P_matrix[i+1,j], self.H_matrix[i+1,j], self.T_su.m_dot/self.params['n_tubes'],j)

            "------- 4) Compute Mass ------------------------"
            
            for i in range(len(self.T_matrix[:,0])):
                if np.mod(i,2) == 0: # Bank Fluid
                    self.D_matrix[i,:], self.M_matrix[i,:], self.x_matrix[i,:] = self.compute_mass(self.H_matrix[i,:], self.P_matrix[i,:], self.B_su.fluid,i)
                    self.M_bank = self.M_bank + sum(self.M_matrix[i,:])
                else: # Internal Fluid
                    self.D_matrix[i,:], self.M_matrix[i,:], self.x_matrix[i,:] = self.compute_mass(self.H_matrix[i,:], self.P_matrix[i,:], self.T_su.fluid,i)
                    self.M_tube = self.M_tube + sum(self.M_matrix[i,:])
                
            "------- 5) Compute Outputs ------------------------"
            
            "4.1) Heat Rate and Mass"
            
            self.Q_dot = self.Q_dot_matrix.sum()

            "4.2) Outlet Conditions - Bundle Side"
            
            h_out_mean_bundle = self.H_matrix[-1,:].mean()
            p_out_mean_bundle = self.P_matrix[-1,:].mean()
            T_out_mean_bundle = self.T_matrix[-1,:].mean()
            
            if self.params['Fin_Side'] == 'C':
                self.B_ex = self.ex_C
                self.T_ex = self.ex_H
            elif self.params['Fin_Side'] == 'H':
                self.B_ex = self.ex_H
                self.T_ex = self.ex_C
            else:
                print("Fin_Side shall either be 'C' or 'H'")

            self.B_ex.set_fluid(self.B_su.fluid)
            self.B_ex.set_m_dot(self.B_su.m_dot)
            
            self.B_ex.set_h(h_out_mean_bundle)
            self.B_ex.set_p(p_out_mean_bundle)
           
            
            "4.3) Outlet Conditions - Tube Side"

            p_out_tube = self.P_matrix[1,-1]
            #This assumes the pressure at the outlet of the first tube row represents the overall tube-side outlet pressure
            
            h_out_mean_tube = 0
            
            for i in range(len(self.H_matrix[:,0])-1):
                if np.mod(i,2) == 1:
                    h_out_mean_tube = h_out_mean_tube + self.H_matrix[i,-1]
                    
            h_out_mean_tube = h_out_mean_tube / self.params['n_rows']
            #By assuming an homogeneous repartition of the total flow rate over the tubes
            
                      
            self.T_ex.set_fluid(self.T_su.fluid)
            self.T_ex.set_m_dot(self.T_su.m_dot)
            
            self.T_ex.set_h(h_out_mean_tube)
            self.T_ex.set_p(p_out_tube)
            T_out_mean_tube=self.T_ex.T
                       
            self.solved = True
            
            
        return
