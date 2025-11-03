# -*- coding: utf-8 -*-

"""
Created on Thu Mar 20 11:17:07 2025

@authors: Samuel Gendebien & Basile Chaudoir
"""

"""
Source for the model : Manufacturing cost model for heat exchangers optimization

Antonio C. Caputo, Pacifico M. Pelagagge, Paolo Salini
"""

import numpy as np
import matplotlib.pyplot as plt

from toolbox.economics.cpi_data import actualize_price

class EconomicParameters:
    """Economic parameters for cost calculations"""
    def __init__(self):
        # Energy cost (€/kWh)
        self.Ce = 0.12
        
        # Interest rate 
        self.r = 0.05
        
        # Hours per years (h/year)
        self.h=2000
        
        self.operation_per_part = {
            'Rolled Shell' : 
                ['Cutting', 'Beveling', 'Welding', 'Rolling', 'Welding Check'],
            'Standard Tube Shell' : 
                ['Cutting', 'Welding', 'Welding Check'],
            'Baffle' : 
                ['Cutting', 'Beveling', 'Drilling'],
            'Flange' : 
                ['Cutting', 'Beveling', 'Drilling'], 
            'Tubesheet' : 
                ['Cutting', 'Beveling', 'Drilling'], 
            'Rolled Channel' : 
                ['Cutting', 'Beveling', 'Welding', 'Rolling', 'Welding Check'],
            'Standard Channel' : 
                ['Cutting', 'Welding', 'Welding Check'],
            'Removable Cover' :
                ['Cutting', 'Beveling', 'Drilling'],
            'Tube' :
                # ['Assembly', 'Cutting', 'Expansion', 'Welding', 'Welding Check'],
                ['Assembly', 'Expansion'],
            'Bolts' : 
                ['Assembly'],
            'Spacer' : 
                ['Assembly'],
            'Tie rods' : 
                ['Assembly']
            }
                        
        # Auxiliary costs [€/su]
        C_aux_b  = 0  # Beveling
        C_aux_c  = 0  # Cutting
        C_aux_cw = 5  # Welding check
        C_aux_d  = 1  # Drilling
        C_aux_e  = 0  # Tube Expansion
        C_aux_r  = 0  # Rolling
        C_aux_w  = 10 # Welding 
        
        self.C_aux = {
            'Beveling': C_aux_b,
            'Cutting': C_aux_c,
            'Welding Check' : C_aux_cw,
            'Drilling' : C_aux_d,
            'Rolling': C_aux_r,
            'Welding': C_aux_w,
            'Expansion': C_aux_e,
            'Assembly' : 0}

        # Material costs [€/kg]
        C_mat_B = 2    # Baffle [€/kg]
        C_mat_bt = 2.5 # Bolt [€/kg]
        C_mat_FL = 2   # Flange [€/kg]
        C_mat_S = 2    # Shell [€/kg]
        C_mat_SP = 2   # Spacer [€/kg]
        C_mat_T = 3    # Tube [€/kg]
        C_mat_TR = 2   # Tie rod [€/kg]
        C_mat_TS = 3.5 # Tube sheet [€/kg]

        self.C_mat = {
            'Baffle'              : C_mat_B,
            'Bolts'               : C_mat_bt,
            'Flange'              : C_mat_FL,
            'Removable Cover'     : C_mat_B,
            'Rolled Channel'      : C_mat_B,
            'Standard Channel'    : C_mat_B,
            'Rolled Shell'        : C_mat_S,
            'Standard Tube Shell' : C_mat_S,
            'Spacer'              : C_mat_SP,
            'Tube'                : C_mat_T,
            'Tie rods'            : C_mat_TR,
            'Tubesheet'           : C_mat_TS
            }
        
        # Load/Unload times [s]
        T_LU_b   = 30  # Beveling
        T_LU_c   = 40  # Cutting
        T_LU_cw  = 180 # Weilding check
        T_LU_db  = 240 # Baffle drilling
        T_LU_dTS = 120 # Tubesheet drilling
        T_LU_r   = 120 # Rolling
        T_LU_w   = 300 # Welding
        
        self.T_LU = {
            'Rolled Shell' : 
                {'Cutting': T_LU_c, 'Beveling': T_LU_b, 'Welding': T_LU_w, 'Rolling': T_LU_r, 'Welding Check' : T_LU_cw},
            'Standard Tube Shell' : 
                {'Cutting' : T_LU_c, 'Welding': T_LU_w, 'Welding Check' : T_LU_cw},
            'Baffle' : 
                {'Cutting': T_LU_c, 'Beveling': T_LU_b, 'Drilling' : T_LU_db},
            'Flange' : 
                {'Cutting': T_LU_c, 'Beveling': T_LU_b, 'Drilling' : T_LU_dTS}, 
            'Tubesheet' : 
                {'Cutting': T_LU_c, 'Beveling': T_LU_b, 'Drilling' : T_LU_dTS}, 
            'Rolled Channel' : 
                {'Cutting': T_LU_c, 'Beveling': T_LU_b, 'Welding': T_LU_w, 'Rolling': T_LU_r, 'Welding Check' : T_LU_cw},
            'Standard Channel' : 
                {'Cutting' : T_LU_c, 'Welding' : T_LU_w, 'Welding Check' : T_LU_cw},
            'Removable Cover' :
                {'Cutting': T_LU_c, 'Beveling': T_LU_b, 'Drilling' : T_LU_dTS},
            'Tube' :
                {'Cutting' : T_LU_c, 'Welding' : T_LU_w, 'Welding Check' : T_LU_cw}
            }
        
        # Setup time [min]
        T_SU_b   = 10  # Beveling
        T_SU_c   = 15  # Cutting
        T_SU_cw  = 5   # Weilding check
        T_SU_db  = 1   # Baffle drilling
        T_SU_dTS = 1   # Tubesheet drilling
        T_SU_r   = 10  # Rolling
        T_SU_w   = 25  # Welding
        
        self.T_SU = {
            'Rolled Shell' : 
                {'Cutting': T_SU_c, 'Beveling': T_SU_b, 'Welding': T_SU_w, 'Rolling': T_SU_r, 'Welding Check' : T_SU_cw},
            'Standard Tube Shell' : 
                {'Cutting' : T_SU_c, 'Welding': T_SU_w, 'Welding Check' : T_SU_cw},
            'Baffle' : 
                {'Cutting': T_SU_c, 'Beveling': T_SU_b, 'Drilling' : T_SU_db},
            'Flange' : 
                {'Cutting': T_SU_c, 'Beveling': T_SU_b, 'Drilling' : T_SU_dTS}, 
            'Tubesheet' : 
                {'Cutting': T_SU_c, 'Beveling': T_SU_b, 'Drilling' : T_SU_dTS}, 
            'Rolled Channel' : 
                {'Cutting': T_SU_c, 'Beveling': T_SU_b, 'Welding': T_SU_w, 'Rolling': T_SU_r, 'Welding Check' : T_SU_cw},
            'Standard Channel' : 
                {'Cutting' : T_SU_c, 'Welding' : T_SU_w, 'Welding Check' : T_SU_cw},
            'Removable Cover' :
                {'Cutting': T_SU_c, 'Beveling': T_SU_b, 'Drilling' : T_SU_dTS},
            'Tube' :
                {'Cutting' : T_SU_c, 'Welding' : T_SU_w,'Welding Check' : T_SU_cw}
            }        
        
        # Specific times [s/item]
        ST_e = 15      # Tube expansion time [s/item]
        ST_in_bt = 30  # Bolt insertion [s/item]
        ST_in_sp = 15  # Spacer insertion [s/item]
        ST_in_T = 3    # Tube insertion [s/item]
        ST_in_td = 3   # Tie rod insertion [s/item]

        self.ST_in = {
            'Tube' : {'Assembly' : ST_in_T, 'Expansion' : ST_e},
            'Bolts' : {'Assembly' : ST_in_bt},
            'Spacer' : {'Assembly' : ST_in_sp},
            'Tie rods' : {'Assembly' : ST_in_td},            
            }
        
        # Worker(s) per operation [-]
        m_b = 1    # Beveling machine [-]
        m_c = 1    # Cutting machine [-]
        m_cw = 1   # Welding Check Machine [-]
        m_d = 1    # Drilling machine [-]]
        m_r = 1    # Rolling machine [-]
        m_w = 1    # Welding machine [-]
        m_e = 1    # Tube Expansion machine [-]
        
        self.m = {
            'Beveling' : m_b,
            'Cutting'  : m_c,
            'Welding Check' : m_cw,
            'Drilling' : m_d,
            'Rolling' : m_r,
            'Welding' : m_w,
            'Expansion' : m_e,
            'Assembly' : 1
            }
        
        # Equipment investment costs [€]
        I_b = 20000    # Beveling machine [€]
        I_c = 120000   # Cutting machine [€]
        I_cw = 30000   # Channel welding [€]
        I_d = 50000    # Drilling machine [€]
        I_r = 100000   # Rolling machine [€]
        I_w = 75000    # Welding machine [€]
        I_e = 5000     # Tube Expansion machine [€]
        
        # Amortization time spans [yr]
        NYr_b = 5      # Beveling [yr]
        NYr_c = 5      # Channel [yr]
        NYr_cw = 10    # Welding Check [yr]
        NYr_d = 5      # Drilling [yr]
        NYr_r = 5      # Rolling [yr]
        NYr_w = 5      # Welding [yr]
        NYr_e = 5      # Equipment [yr]
        
        self.equ = {
            'Beveling' : {'Invest' : I_b, 'Years' : NYr_b},
            'Cutting'  : {'Invest' : I_c, 'Years' : NYr_c},
            'Welding Check' : {'Invest' : I_cw, 'Years' : NYr_cw},
            'Drilling' : {'Invest' : I_d, 'Years' : NYr_d},
            'Rolling' : {'Invest' : I_r, 'Years' : NYr_r},
            'Welding' : {'Invest' : I_w, 'Years' : NYr_w},
            'Expansion' : {'Invest' : I_e, 'Years' : NYr_e},
            'Assembly' : {'Invest' : 0, 'Years' : 0}
            }
        
        # Power consumption [kW]
        P_b = 15        # Beveling [kW]
        P_c = 100       # Cutting [kW]
        P_cw = 10       # Welding Check [kW]
        P_d = 10        # Drilling [kW]
        P_r = 100       # Rolling [kW]
        P_w = 100       # Welding [kW]
        P_e = 5         # Equipment [kW]
        
        self.Power = {
            'Beveling' : P_b,
            'Cutting'  : P_c,
            'Welding Check' : P_cw,
            'Drilling' : P_d,
            'Rolling' : P_r,
            'Welding' : P_w,
            'Expansion' : P_e,
            'Assembly' : 0
            }
        
        # Processing speeds [m/min] : Results from an optimization to fit operational costs from the source article 
        v_b = 0.3085   # Beveling [m/min]
        v_c = 0.1      # Cutting [m/min] assumed to be equal to 1 (to be better assessed)
        v_cw= 0.09     # Welding Check [m/min]
        v_d = 0.0562   # Drilling [m/min]
        v_r = 0.6189   # Rolling [m/min]
        v_w = 0.6145   # Welding [m/min]


        self.v = {
            'Beveling' : v_b,
            'Cutting'  : v_c,
            'Welding Check' : v_cw,
            'Drilling' : v_d,
            'Rolling' : v_r,
            'Welding' : v_w,
            'Expansion' : 0,
            'Assembly' : 0 
            }
        
        # Surface treatment cost [€/m²]
        self.C_gr = 2       # Grinding [€/m²]

        # Surface preparation costs [€/m²]
        self.C_pa = 4       # Painting [€/m²]
        self.C_pk = 5       # Pickling [€/m²]
        self.C_sb = 3       # Sandblasting [€/m²]

        # Electrode & welding parameters
        self.EC = 1.2       # Electrode cost [€/kg]
        self.EMY = 0.97     # Deposition efficiency ratio [-]
        self.EWL = 0.0053   # Electrode weight per unit length [kg/m]

        # Gas cost & consumption
        self.GC = 15        # Gas cost [€/m³]
        self.GFR = 1.4      # Gas flow rate [m³/h]

        # Dimensional parameters (length or spacing)
        self.L_l = 3        # Length of the lead [mm]
        self.L_ot = 5       # Overtravel distance [mm]
        self.L_Pst = 6      # Plate standard length [m]
        self.L_pt = 5       # Pretravel [mm]
        
        # Labor rate [€/h]
        self.LR = 22 #[€/h]

        # Operating factor [-]
        self.OF = 0.5       # Operating factor [-]
        
        # Electrical parameters
        self.WC = 150       # Welding current [A]
        self.WE = 0.9       # Welder electric efficiency [-]
        self.WFR = 2.5      # Wire feed rate [m/min]
        self.WP_Pst = 1.5   # Welding pipe [m]
        self.WP_pk = 5      # Welding plate thickness [cm]
        self.WV = 20        # Welding voltage [V]
        
        # Hourly operational costs [€/h]
        self.C_HO_b = 0.3   # Baffle [€/h] 
        self.C_HO_c = 10    # Channel [€/h]
        self.C_HO_cw = 5    # Channel welding [€/h]
        self.C_HO_d = 1     # Drilling [€/h]
        self.C_HO_e = 0     # Equipment [€/h]
        self.C_HO_r = 0     # Rolling [€/h]
        self.C_HO_w = 0     # Welding [€/h]
        
class HeatExchangerCost:
    def __init__(self, D_S_i, t_S ,r, L_B, t_B, B_c, N_TS, D_r, L_T, D_T_o, D_T_i, N_T, pitch_r, L_CH, N_CH, t_TS, N_FL, t_FL, t_RC ,D_SP_e, D_SP_i, N_TR, D_TR, L_TR, N_Bt):
        
        """
        Initialize with heat exchanger parameters
        D_S_i: Inlet shell diameter (m)
        t_S: Shell thickness (m)
        L_T: Tube length (m)
        N_T: Number of tubes (-)
        N_TS: Number of tube sheets (-)
        B_c: Baffle cut ratio (-)
        D_r: Diameter ratio for flanges and tube sheets (-)
        r: Number of shell sections (-)
        D_t: Tube external diameter (m)
        D_t_i: Tube internal diameter (m)
        L_B: Baffle spacing (m)
        t_B: Baffle thickness (m)
        pitch: Pitch ratio (-)
        L_CH: channel length (m)
        N_CH: number of channel (m)
        t_TS: thickness of (m)
        t_RC: thickness of removable cover (m)
        D_sp_e: External spacers diameter (m)
        D_sp_i: Internal spacers diameter (m)
        N_TR: Number of tie rods (-) 
        D_TR: Diameter of tie rods (m)
        L_TR: Length of tie rods (m)
        """

        self.params = EconomicParameters()  # Linking to economic parameters

        # Geometry parameters
        
        #Shell characteristics
        self.D_S_i = D_S_i               # Inlet shell diameter (m)
        self.t_S = t_S                   # Shell thickness (m)
        self.D_S_o=self.D_S_i+ self.t_S  # Inlet shell diameter (m)
        self.r = r                       # Number of shell sections (-)

        #Baffle characteristics
        self.L_B = L_B          # Baffle spacing (m)
        self.t_B = t_B          # Baffle thickness (m)
        self.B_c = B_c          # Baffle cut ratio (-)

        #Tube sheets
        self.N_TS = N_TS        # Number of tube sheets (-)
        self.D_r = D_r          # Diameter ratio for flanges and tube sheets (-)

        #Tube characteristics
        self.L_T = L_T          # Tube length (m)
        self.D_T_o = D_T_o      # Tube external diameter (m)
        self.D_T_i = D_T_i      # Tube internal diameter (m)
        self.N_T = N_T          # Number of tubes (-)
        self.pitch_r = pitch_r  # Pitch ratio (-)

        #Channel characteristics
        self.L_CH = L_CH        # Channel length (m)
        self.N_CH = N_CH        # Number of channel (m)
        self.t_TS = t_TS        # Channel thickness (m)
        
        #Removable cover characteristics
        self.t_RC = t_RC        # Thickness of removable cover (m)
        #Note: the number of removable cover is assumed to be equal to the number of channels (see Appendix B)
        
        #Flange characteristics
        self.N_FL= N_FL         # Number of flanges
        self.t_FL = t_FL        # Flange thickness (m)
        
        #Spacers characteristics
        self.D_SP_e = D_SP_e    # External spacers diameter (m)
        self.D_SP_i = D_SP_i    # Internal spacers diameter (m)
        
        #Tie rods characteristics
        self.N_TR = N_TR        # Number of tie rods (-) 
        self.D_TR = D_TR        # Diameter of tie rods (m)
        self.L_TR = L_TR        # Length of tie rods (m)
        
        #Bolts
        self.N_Bt=N_Bt          #Number of bolts (-)
        
        self.mass = {
            'Total' : 0
            }
        
        """
        Geometry parameters deduction from the defined inputs
        """
        # Calculate number of baffles
        self.N_B = np.floor(self.L_T / self.L_B) - 1
        
        # Calculate k4 = cos^-1((0.5-Bc)/0.5)
        self.k4 = np.arccos((0.5 - self.B_c) / 0.5)
        
        # Manufacturing consideration regarding the shell
        self.shell_type = "Rolled Shell" if self.D_S_i > 0.6 else "Standard Tube Shell"  # Shell is rolled if Ds > 600mm
        
        # Manufacturing consideration regarding the shell
        self.channel_type = "Rolled Channel" if self.D_S_i > 0.6 else "Standard Channel"  # Shell is rolled if Ds > 600mm
        
        # Tube pitch 
        self.pitch_T = pitch_r * self.D_T_o  # Tube pitch (m)
   
    
        # Operations
        self.operation_vect = ['Assembly', 'Beveling', 'Cutting', 'Drilling', 'Expansion', 'Rolling', 'Welding', 'Welding Check']

        self.operations = {}
        
        for operation in self.operation_vect:
            current_dict = {
                'investment': self.params.equ[operation]['Invest'],
                'life': self.params.equ[operation]['Years'],
                'power': self.params.Power[operation],
                'workers': self.params.m[operation],
                'aux_cost': self.params.C_aux[operation],
                'velocity': self.params.v[operation]}
            
            self.operations[operation] = current_dict
   
    def calculate_volume(self, component):
        """Calculate volume for each heat exchanger subassembly based on Appendix B formulas"""
        if component == 'Rolled Shell' or component == 'Standard Tube Shell' :
            return np.pi * ((self.D_S_i + 2*self.t_S)**2 - self.D_S_i**2) / 4 * self.L_T  # Shell volume
        elif component == 'Baffle':
            return self.D_S_i**2 * (np.pi / 4 * (1 - 1 / (np.pi) * self.k4) + 0.5 * np.sin(self.k4) * (0.5 - self.B_c)) * self.t_B * self.N_B  # Baffle volume
        elif component == 'Tubesheet':
            return np.pi * ((self.D_S_i + 2*self.D_r)**2 - self.D_S_i**2) / 4 * self.t_TS * self.N_TS  # Tube-sheet volume
        elif component == 'Tube':
            return np.pi * (self.D_T_o**2 - self.D_T_i**2) / 4 * self.L_T * self.N_T  # Tube volume
        elif component == 'Rolled Channel' or component == 'Standard Channel':
            return np.pi * ((self.D_S_i + 2*self.t_S)**2 - self.D_S_i**2) / 4 * self.L_CH * self.N_CH  # Channel volume
        elif component == 'Removable Cover':
            return np.pi * ((self.D_S_i + 2*self.D_r)**2) / 4 * self.t_RC * self.N_CH  # Removable cover volume
        elif component == 'Flange':
            return np.pi * ((self.D_S_i + 2*self.D_r)**2 - self.D_S_i**2) / 4 * self.t_FL * self.N_FL  # Flange volume (Equation 51)
        elif component == 'Tie rods':
            return self.N_TR * np.pi * (self.D_TR**2 / 4) * self.L_TR  # Tie rod volume (Equation 51 bis)
        elif component == 'Spacer':
            return self.N_B * np.pi * ((self.D_SP_e**2 - self.D_SP_i**2) / 4) * self.L_B * self.N_TR  # Spacer volume (Equations 5 ter)
        elif component == 'Bolts':
            return self.N_Bt * 3E-6   # Bolts volume (Equations 5 ter)
        else: 
            raise ValueError(f"Unknown component '{component}'")
                

    def calculate_material_cost(self, component):
        """Calculate material cost based on volume and density"""
        volume = self.calculate_volume(component)
        density = 7850  # kg/m³ for steel
                
        if component in self.params.C_mat:
            cost = volume * density * self.params.C_mat[component]
            self.mass[component] = volume*density
            self.mass['Total'] += volume*density
        else:
            cost = 0
            
        return cost

    def compute_processing_length(self):
        """"Assumptions about the shell rolled plate"""
        self.N_RP=3
        self.W_RP=self.L_T/self.N_RP
        
        """"Assumptions about the standard tube shell"""
        self.N_T_tr=3 #the number of tube trunks
        
        """"Assumptions about the circonferential bolt distance"""
        self.bd= (np.pi*self.D_S_i*(self.D_r + 1))/self.N_Bt #The circonferential bolt distance is assumed to be equal to 5 cm
        
        """"Assumptions about the channel widthee"""
        self.W_CH=self.D_S_i # 0.15 #the channel width is assumed to be equal to 0.15[m]
        
        """"Assumptions about the number of removable covers"""
        self.N_RC=self.N_CH #the number of removable cover is equal to the number of channel 
        
        # Additional drilling lengths
        L_D_add = (self.params.L_l + self.params.L_ot + self.params.L_pt)*1e-3 # m
                
        self.processing_lengths = {
            "Rolled Shell": {
                "Beveling": 2*(np.pi * self.D_S_i + self.W_RP)*self.N_RP,
                "Cutting": 2*(np.pi * self.D_S_i + self.W_RP)*self.N_RP,
                "Welding": self.W_RP*self.N_RP + np.pi * self.D_S_i * (1 + 1/self.N_RP)*self.N_RP, #Sum of longitudinal and circumferntial welding
                "Rolling": np.pi * self.D_S_i*self.N_RP},
            
            "Standard Tube Shell": {
                "Cutting": np.pi * self.D_S_o*self.N_T_tr,
                "Welding": np.pi * self.D_S_i * (1 - 1/self.N_T_tr)*self.N_T_tr},
            
            "Baffle": {
                "Beveling": self.D_S_i * ((np.pi - self.k4) + np.sin(self.k4))*self.N_B,
                "Cutting": self.D_S_i * ((np.pi - self.k4) + np.sin(self.k4))*self.N_B,
                "Drilling": ((1 - self.k4/np.pi) + (2/np.pi) * np.sin(self.k4) * (1/2 - self.B_c)) * (self.t_B+L_D_add) * self.N_T*self.N_B},
            
            "Flange": {
                "Beveling": 2 * np.pi * self.D_S_i * (1 + self.D_r)*self.N_FL,
                "Cutting": 2 * np.pi * self.D_S_i * (1 + self.D_r)*self.N_FL,
                "Drilling": (np.pi * self.D_S_i * (1 + self.D_r) / self.bd) * (self.t_FL+L_D_add) *self.N_FL},
            
            "Tubesheet": {
                "Beveling": np.pi * self.D_S_i * (1 + 2 * self.D_r)*self.N_TS,
                "Cutting": np.pi * self.D_S_i * (1 + 2 * self.D_r)*self.N_TS,
                "Drilling": (self.N_T + self.D_S_i * (1 + self.D_r) / self.bd) * (self.t_TS*self.N_TS)}, #+L_D_add)*self.N_TS},
            
            "Rolled Channel": {
                "Beveling": 2 * (np.pi * self.D_S_i + self.W_CH)*self.N_CH,
                "Cutting": 2 * (np.pi * self.D_S_i + self.W_CH)*self.N_CH,
                "Welding": self.W_CH*self.N_CH + 2 * np.pi * self.D_S_o*self.N_CH,
                "Rolling": np.pi * self.D_S_i*self.N_CH},
            
            "Standard Channel": {
                "Cutting": np.pi * self.D_S_o * self.N_CH,
                "Welding": 2 * np.pi * self.D_S_o* self.N_CH},
            
            "Removable Cover": {
                "Beveling": np.pi * self.D_S_i * (1 + 2 * self.D_r)* self.N_RC,
                "Cutting": np.pi * self.D_S_i * (1 + 2 * self.D_r)* self.N_RC,
                "Drilling": (self.D_S_i * (1 + self.D_r) / self.bd) * (self.t_RC+L_D_add) * self.N_RC},
            
            "Tube": {
                "Cutting": np.pi * self.D_T_o * self.N_T,
                "Welding": np.pi * self.D_T_o * self.N_T_tr * self.N_T}
        }
        
        # self.processing_lengths = {
        #     "Rolled Shell": {
        #         "Beveling": 2*(np.pi * self.D_S_i + self.W_RP),
        #         "Cutting": 2*(np.pi * self.D_S_i + self.W_RP),
        #         "Welding": self.W_RP + np.pi * self.D_S_i * (1 + 1/self.N_RP), #Sum of longitudinal and circumferntial welding
        #         "Rolling": np.pi * self.D_S_i},
            
        #     "Standard Tube Shell": {
        #         "Cutting": np.pi * self.D_S_o,
        #         "Welding": np.pi * self.D_S_i * (1 - 1/self.N_T_tr)},
            
        #     "Baffle": {
        #         "Beveling": self.D_S_i * ((np.pi - self.k4) + np.sin(self.k4)),
        #         "Cutting": self.D_S_i * ((np.pi - self.k4) + np.sin(self.k4)),
        #         "Drilling": ((1 - self.k4/np.pi) + (2/np.pi) * np.sin(self.k4) * (1/2 - self.B_c)) * (self.t_B+L_D_add) * self.N_T},
            
        #     "Flange": {
        #         "Beveling": 2 * np.pi * self.D_S_i * (1 + self.D_r),
        #         "Cutting": 2 * np.pi * self.D_S_i * (1 + self.D_r),
        #         "Drilling": (np.pi * self.D_S_i * (1 + self.D_r) / self.bd) * (self.t_FL+L_D_add)},
            
        #     "Tubesheet": {
        #         "Beveling": np.pi * self.D_S_i * (1 + 2 * self.D_r),
        #         "Cutting": np.pi * self.D_S_i * (1 + 2 * self.D_r),
        #         "Drilling": (self.N_T + self.D_S_i * (1 + self.D_r) / self.bd) * (self.t_TS+L_D_add)},
            
        #     "Rolled Channel": {
        #         "Beveling": 2 * (np.pi * self.D_S_i + self.W_CH),
        #         "Cutting": 2 * (np.pi * self.D_S_i + self.W_CH),
        #         "Welding": self.W_CH + 2 * np.pi * self.D_S_o,
        #         "Rolling": np.pi * self.D_S_i},
            
        #     "Standard Channel": {
        #         "Cutting": np.pi * self.D_S_o,
        #         "Welding": 2 * np.pi * self.D_S_o},
            
        #     "Removable Cover": {
        #         "Beveling": np.pi * self.D_S_i * (1 + 2 * self.D_r),
        #         "Cutting": np.pi * self.D_S_i * (1 + 2 * self.D_r),
        #         "Drilling": (self.D_S_i * (1 + self.D_r) / self.bd) * (self.t_RC+L_D_add)},
            
        #     "Tube": {
        #         "Cutting": np.pi * self.D_T_o,
        #         "Welding": np.pi * self.D_T_o * self.N_T_tr}
        # }
        
        self.batches = {
            'Rolled Shell' : self.r*self.N_RP,
            'Standard Tube Shell' : self.r*self.N_T_tr,
            'Baffle' : self.N_B,
            'Tubesheet' : self.N_TS,
            'Tube' : self.N_T,
            'Rolled Channel' : self.N_CH,
            'Standard Channel' : self.N_CH,
            'Removable Cover' : self.N_RC,
            'Flange' : self.N_FL,
            'Tie rods' : self.N_TR,
            'Bolts' : self.N_Bt,
            'Spacer' : self.N_B*2
            }
        
        return self.processing_lengths, self.batches


    def hourly_depreciation_cost(self, equipment_investment, equipment_life):
        """
        Computes the hourly depreciation cost of processing equipment.
        
        Parameters:
        r : Interest (discount) rate (e.g., 0.05 for 5%)
        equipment_life: Equipment lifetime (years)
        h: Number of yearly working hours  
        Returns:
        Hourly depreciation cost (€/hour)
        """
        # Compute capital recovery factor (τ_k)
        tau_k = (self.params.r * (1 + self.params.r) ** equipment_life) / ((1 + self.params.r) ** equipment_life - 1)
        
        # Compute hourly depreciation cost (CH_E,k)
        return equipment_investment * tau_k / self.params.h


    def calculate_energy_cost(self, power):
        """Calculate energy cost
        CEk = Pk * Cew * t
        where:
        - Pk is power consumption (kW)
        - Cew is energy cost (€/kWh)
        """
        return power * self.params.Ce  # Power in kW, Ce in €/kWh

    def calculate_hourly_operation_cost(self):
        """Calculate cost for a manufacturing operation"""
        
        self.C_HO = {} 
        
        for operation in self.operation_vect:
            # Hourly Labour costs
            if operation in self.params.m:
                C_H_L = self.params.m[operation]*self.params.LR
            else:
                C_H_L = 0
            
            # Electricity Costs
            if operation in self.params.Power:
                C_H_E = self.params.Power[operation]*self.params.Ce
            else:
                C_H_E = 0
                
            # Auxiliary Costs
            if operation in self.params.C_aux:
                C_H_A = self.params.C_aux[operation]
            else:
                C_H_A = 0
            
            # Depreciation costs
            if operation in self.params.equ:
                if self.params.equ[operation]['Invest'] != 0:
                    C_H_D = self.hourly_depreciation_cost(self.params.equ[operation]['Invest'], self.params.equ[operation]['Years'])
                else:
                    C_H_D = 0
            else:
                C_H_D = 0
                        
            self.C_HO[operation] = {
                'Labor' : C_H_L,
                'Energy' : C_H_E,
                'Depreciation' : C_H_D,
                'Hourly' : C_H_L + C_H_E + C_H_D,
                'Auxiliaries' : C_H_A 
                }
            
        return self.C_HO

    def calculate_operation_cost(self):
        """Calculate cost for a manufacturing operation"""
        
        self.C_comp = {}
        C_op_tot = 0
        
        for component in [self.shell_type, 'Baffle', 'Tubesheet', 'Tube', self.channel_type, 'Removable Cover', 'Flange', 'Tie rods', 'Spacer', 'Bolts']: # self.params.operation_per_part.keys():
            # print(f"---------------------------------------------------")

            # print(f"Component : {component}")
            self.C_op = {}
            C_tot = 0 
            for operation in self.params.operation_per_part[component]:
                # print(f"----------------")

                # print(f"Operation : {operation}")

                C_HO = self.C_HO[operation]['Hourly']
                C_L = self.C_HO[operation]['Labor']
                C_D = self.C_HO[operation]['Depreciation']
                C_AUX = self.C_HO[operation]['Auxiliaries']
                
                batch_size = self.batches[component] # 1 HX built at a time is considered
                
                if operation != "Assembly" and operation != "Expansion":
                    # print(f"Operation 2 : {operation}")
                    if component in self.params.T_SU:
                        T_SU = self.params.T_SU[component][operation]/3600
                    else: 
                        T_SU = 0
    
                    if component in self.params.T_LU:
                        T_LU = self.params.T_LU[component][operation]/3600
                    else: 
                        T_LU = 0
                        
                    if operation in self.processing_lengths[component]:
                        # print(f"Operation 3 : {operation}")

                        process_len = self.processing_lengths[component][operation]
                        process_speed = self.params.v[operation]
                        self.C_op[operation] = ((process_len/(process_speed*60))*C_HO + (C_L + C_D)*(T_LU + T_SU/batch_size) + C_AUX/batch_size)/self.params.OF
                        C_tot += self.C_op[operation]
                    else:
                        self.C_op[operation] = ((C_L + C_D)*(T_LU + T_SU/batch_size) + C_AUX/batch_size)/self.params.OF
                        C_tot += self.C_op[operation]
                else:
                    if operation == "Assembly":
                        if component == "Tube":
                            T_in = self.params.ST_in[component][operation]*self.N_T*(self.N_TS + self.B_c*self.N_B)
                            self.C_op[operation] = ((C_L + C_D)*(T_in/3600))/self.params.OF
                            C_tot += self.C_op[operation]
                        else:
                            T_in = self.params.ST_in[component][operation]*batch_size
                            self.C_op[operation] = ((C_L + C_D)*(T_in/3600))/self.params.OF                            
                            C_tot += self.C_op[operation]
                    else:
                        T_e = self.params.ST_in[component][operation]*self.N_T*self.N_TS
                        self.C_op[operation] = ((C_L + C_D)*(T_e/3600))/self.params.OF          
                        C_tot += self.C_op[operation]
                
                # print(f"Operation cost : {self.C_op[operation]}")
                
            self.C_comp[component] = self.C_op.copy()
            self.C_comp[component]['Total'] = C_tot
            C_op_tot += C_tot
            
        return C_op_tot

    def calculate_total_cost(self):
        """Calculate total cost """
        total_material_cost = sum(self.calculate_material_cost(component) for component in [self.shell_type, 'Baffle', 'Tubesheet', 'Tube', self.channel_type, 'Removable Cover', 'Flange','Tie rods', 'Spacer', 'Bolts'])
        hourly_operation_costs = self.calculate_hourly_operation_cost()
        processing_lengths, batches = self.compute_processing_length()
        total_operation_cost = self.calculate_operation_cost()
        
        # print(f"Total costs [€] : {total_material_cost + total_operation_cost}")
        # print(f"Total material costs [€] : {total_material_cost}")
        # print(f"Total operation costs [€] : {total_operation_cost}")
        
        total_costs = total_material_cost + total_operation_cost
        
        total_costs_actual = actualize_price(total_costs, 2015, "EUR")
        self.total_material_cost = actualize_price(total_material_cost, 2015, "EUR")
        self.total_operation_cost = actualize_price(total_operation_cost, 2015, "EUR")
        
        return total_costs_actual

    def print_cost_breakdown(self):
        """Print cost breakdown by component and operation with labeled plots"""
        component_costs = {}
        total_material_cost = 0
        total_operation_cost = 0
        components = []
    
        # Iterate over each component
        for component in [self.shell_type, 'Baffle', 'Tubesheet', 'Tube', self.channel_type, 'Removable Cover', 'Flange', 'Tie rods', 'Spacer', 'Bolts']:
            material_cost = self.calculate_material_cost(component)
            component_costs[component] = {'material': material_cost, 'operations': {}}
            total_material_cost += material_cost

            # Iterate over each operation
            for operation in ['Beveling', 'Cutting', 'Welding', 'Drilling', 'Rolling']:
                length = self.compute_processing_length(operation, component)
                cost = self.calculate_operation_cost(operation, length, component)
                component_costs[component]['operations'][operation] = cost['total']
                total_operation_cost += cost['total']
    
            components.append(component)  # Add component to list for plotting
    
        grand_total = total_material_cost + total_operation_cost
    
        # Calculate data for plots
        material_costs = [details['material'] for details in component_costs.values()]
        material_percentages = [(cost / total_material_cost * 100) if total_material_cost > 0 else 0 for cost in material_costs]
        processing_percentages = [(sum(details['operations'].values()) / total_operation_cost * 100) if total_operation_cost > 0 else 0 for details in component_costs.values()]
        total_component_costs = [details['material'] + sum(details['operations'].values()) for details in component_costs.values()]
        total_cost_percentages = [(cost / grand_total * 100) if grand_total > 0 else 0 for cost in total_component_costs]
    
        print("\nCost Breakdown by Component:")
        print("---------------------------")
        for comp, details in component_costs.items():
            material_percentage = (details['material'] / total_material_cost * 100) if total_material_cost > 0 else 0
            operation_total = sum(details['operations'].values())
            operation_percentage = (operation_total / total_operation_cost * 100) if total_operation_cost > 0 else 0
            total_component_cost = details['material'] + operation_total
            component_percentage = (total_component_cost / grand_total * 100) if grand_total > 0 else 0
    
            print(f"{comp.capitalize()}:")
            print(f"  Material Cost: {details['material']:.2f} € ({material_percentage:.2f}%)")
            for op, op_cost in details['operations'].items():
                print(f"  {op.capitalize()}: {op_cost:.2f} €")
            print("  ---------------------")
            print(f"  Total Operations Cost: {operation_total:.2f} € ({operation_percentage:.2f}%)")
            print(f"  **Total Component Cost: {total_component_cost:.2f} € ({component_percentage:.2f}%)**")
            print("")
    
        print("\nTotal Costs:")
        print("---------------------------")
        print(f"Total Material Cost: {total_material_cost:.2f} € (100%)")
        print(f"Total Operations Cost: {total_operation_cost:.2f} € (100%)")
        print(f"Grand Total: {grand_total:.2f} € (100%)")
    
        # Plot 1: Material Cost Percentage per Component (in %, with labels)
        plt.figure(figsize=(10, 6))
        bars = plt.bar(components, material_percentages, color='green')
        plt.xlabel("Components")
        plt.ylabel("Percentage of Total Material Cost (%)")
        plt.title("Material Cost Percentage per Component")
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.ylim(0, 100)
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, height, f'{height:.1f}%', ha='center', va='bottom')
        plt.tight_layout()
        plt.show()
    
        # Plot 2: Percentage of Total Processing Cost per Component (with labels)
        plt.figure(figsize=(10, 6))
        bars = plt.bar(components, processing_percentages, color='purple')
        plt.xlabel("Components")
        plt.ylabel("Percentage of Total Processing Cost (%)")
        plt.title("Processing Cost Percentage Breakdown per Component")
        plt.ylim(0, 100)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, height, f'{height:.1f}%', ha='center', va='bottom')
        plt.tight_layout()
        plt.show()
    
        # Plot 3: Percentage of Total Cost (Material + Processing) per Component (with labels)
        plt.figure(figsize=(10, 6))
        bars = plt.bar(components, total_cost_percentages, color='teal')
        plt.xlabel("Components")
        plt.ylabel("Percentage of Total Cost (%)")
        plt.title("Total Cost (Material + Processing) Percentage per Component")
        plt.ylim(0, 100)
        plt.xticks(rotation=45, ha='right')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, height, f'{height:.1f}%', ha='center', va='bottom')
        plt.tight_layout()
        plt.show()
    
        # Plot 4: Stacked Bar Chart (Total Manufacturing Cost in €, styled like the image)
        plt.figure(figsize=(6, 6))
        bars_material = plt.bar(["S&T"], [total_material_cost], color='lightblue', hatch='//', label="Material Cost")
        bars_operation = plt.bar(["S&T"], [total_operation_cost], bottom=[total_material_cost], color='skyblue', hatch='--', label="Processing Cost")
        plt.ylabel("Total Manufacturing Cost (€)")
        plt.title("Breakdown of Total Manufacturing Costs")
        plt.legend(loc='upper right')
        plt.ylim(0, grand_total * 1.2)
        for bar in bars_material:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, height / 2, f'{height:,.1f}', ha='center', va='center', color='black')
        for bar in bars_operation:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, total_material_cost + height / 2, f'{height:,.1f}', ha='center', va='center', color='black')
        plt.figtext(0.5, -0.05, "Breakdown of total manufacturing costs.", ha="center", fontsize=10)
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    import numpy as np
    
    # Shell
    D_s = 0.762 # m
    t_s = 0.011 # m
    n_series = 1 
    
    # Baffles
    B = 0.7 # m
    t_B = 0.032 # m : # !!! This seems high, it is weird
    BC = 40/100 
    
    # Tubesheets
    sigma = 100 # MPa
    P_s = 2 # MPa
    t_TS = 0.5*D_s*np.sqrt(P_s/sigma)
    n_TS= 2 
    D_r = 2*0.05/D_s
    
    # Tubes
    n_t = 546
    pitch_ratio = 26/20 
    L_t = 7.2 
    D_t_o = 0.02
    D_t_i = 0.8*D_t_o
    
    # Channels
    n_CH = n_series*2
    L_CH = 0.3
    
    # Removable covers
    t_RC = t_s*2
    
    # Flanges
    t_FL = 0.015
    n_FL = n_series*2
    
    # Tie rods
    D_tr = 0.006
    L_TR = L_t
    n_TR = 10
    
    # Spacers
    D_sp_i = D_tr
    D_sp_o = 2*D_sp_i
    
    # Bolts 
    N_bt = 100
    
    # Create calculator with dimensions from Table 5
    calculator = HeatExchangerCost(
        D_S_i=D_s,  
        t_S=t_s , 
        r=n_series, 
        L_B=B, 
        t_B=t_B, 
        B_c = BC, 
        N_TS=n_TS, 
        D_r=D_r, 
        L_T=L_t, 
        D_T_o=D_t_o, 
        D_T_i=D_t_i, 
        N_T=n_t,
        pitch_r=pitch_ratio, 
        L_CH=L_CH, 
        N_CH=n_CH, 
        t_TS=t_TS, 
        N_FL=n_FL, 
        t_FL=t_FL, 
        t_RC=t_RC,
        D_SP_e=D_sp_o, 
        D_SP_i=D_sp_i, 
        N_TR=n_TR, 
        D_TR=D_tr, 
        L_TR=L_TR, 
        N_Bt=N_bt
    )
    
    Costs = calculator.calculate_total_cost()
    
    # # Cost_Breakdown=calculator.print_cost_breakdown()
    # # Test=calculator.compute_processing_length("Drilling","Baffle")
    
    # print(" Material Costs Breakdown (€):")
    # # First, calculate the total material cost without calculating percentages
    # total_material = 0
    # component_costs = {}
    
    # # Loop through each component to calculate its cost
    # for component in [calculator.shell_type, 'Baffle', 'Tubesheet', 'Tube', calculator.channel_type, 'Removable Cover', 'Flange','Tie rods', 'Spacer', 'Bolts']:
    #     cost = calculator.calculate_material_cost(component)
    #     component_costs[component] = cost
    #     total_material += cost
    
    
    # # Now calculate percentages and print the result
    # print("Material Costs Breakdown (€):")
    
    # for component, cost in component_costs.items():
    #     percentage = (cost / total_material) * 100
    #     print(f"{component:16} {cost:8.2f} € ({percentage:5.1f}%)")
    
    # print (f"Total material     {total_material} €")
    
    
    # # Print total material cost
    # print("-" * 70)

    # from scipy.optimize import minimize

    # def objective(speeds_vector):
    #     total_op_costs = 0.95*5715.2
        
    #     target_cost_rolled_shell = total_op_costs*(0.11)
    #     target_cost_baffle = total_op_costs*(0.53)
    #     target_cost_tube = total_op_costs*(0.023)
    #     target_cost_TS = total_op_costs*(0.164)
    #     target_cost_channel = total_op_costs*(0.062)
    #     target_cost_flanges = total_op_costs*(0.047)
        
    #     calculator.params.v['Beveling'] = speeds_vector[0]
    #     calculator.params.v['Cutting'] = speeds_vector[1]
    #     calculator.params.v['Welding'] = speeds_vector[2]
    #     calculator.params.v['Rolling'] = speeds_vector[3]
    #     calculator.params.v['Drilling'] = speeds_vector[4]
    #     calculator.params.v['Welding Check'] = speeds_vector[5]
        
    #     total = calculator.calculate_total_cost()
        
    #     predicted_costs_shell = calculator.C_comp['Rolled Shell']['Total']  # example
    #     predicted_costs_baffle = calculator.C_comp['Baffle']['Total']  # example
    #     predicted_costs_tube = calculator.C_comp['Tube']['Total']  # example
    #     predicted_costs_TS = calculator.C_comp['Tubesheet']['Total']  # example
    #     predicted_costs_channel = calculator.C_comp['Rolled Channel']['Total']  # example
    #     predicted_costs_flanges = calculator.C_comp['Flange']['Total']  # example

        
    #     res_rolled_shell = (predicted_costs_shell - target_cost_rolled_shell)**2  # squared error
    #     res_baffle = (predicted_costs_baffle - target_cost_baffle)**2  # squared error
    #     res_tube = (predicted_costs_tube - target_cost_tube)**2  # squared error
    #     res_std_TS = (predicted_costs_TS - target_cost_TS)**2  # squared error        
    #     res_channel = (predicted_costs_channel - target_cost_channel)**2  # squared error        
    #     res_flanges = (predicted_costs_flanges - target_cost_flanges)**2  # squared error        
    #     # res_op_cost = 
        
    #     res_tot = res_rolled_shell + res_baffle + res_tube + res_std_TS + res_channel + res_flanges

    #     return res_tot
    
    # x0=[1.5, 1.0, 0.2, 0.2, 0.3, 0.09]
    
    # # Define bounds: each tuple (lower, upper), set lower bound to a small positive number like 1e-3
    # bounds = [(1e-4, 5)] * len(x0)  # This applies positivity to all speed variables
    
    # result = minimize(objective, x0=x0, bounds=bounds, method='L-BFGS-B')

    # print(f"Beveling speed     : {result.x[0]:.4f}")
    # print(f"Cutting speed      : {result.x[1]:.4f}")
    # print(f"Welding speed      : {result.x[2]:.4f}")
    # print(f"Rolling speed      : {result.x[3]:.4f}")
    # print(f"Drilling speed     : {result.x[4]:.4f}")
    # print(f"Welding Check speed: {result.x[5]:.4f}")

