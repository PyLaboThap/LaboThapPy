# -*- coding: utf-8 -*-

"""
Created on Thu Mar 20 11:17:07 2025

@author: Samuel Gendebien
"""
import math
import numpy as math
import matplotlib.pyplot as plt

class EconomicParameters:
    """Economic parameters for cost calculations"""
    def __init__(self):
        # Energy cost (€/kWh)
        self.Ce = 0.12
        
        # Interest rate (%)
        self.r = 0.05
        
        # Hours per years (h/year)
        self.h=2000
        
       # Auxiliary costs [€/su]
        self.C_AUX_b = 0    # Baffle [€/su]
        self.C_AUX_c = 0    # Channel [€/su]
        self.C_AUX_cw = 5   # Channel welding [€/su]
        self.C_AUX_d = 1    # Drilling [€/su]
        self.C_AUX_r = 0    # Rolling [€/su]
        self.C_AUX_w = 10    # Welding [€/su]

        # Surface treatment cost [€/m²]
        self.C_gr = 2       # Grinding [€/m²]

        # Hourly labor costs [€/h]
        self.C_HO_b = 0.3   # Baffle [€/h] 
        self.C_HO_c = 10    # Channel [€/h]
        self.C_HO_cw = 5    # Channel welding [€/h]
        self.C_HO_d = 1     # Drilling [€/h]
        self.C_HO_e = 0     # Equipment [€/h]
        self.C_HO_r = 0     # Rolling [€/h]
        self.C_HO_w = 0     # Welding [€/h]

        # Material costs [€/kg]
        self.C_mat_B = 2    # Baffle [€/kg]
        self.C_mat_bt = 2.5 # Bolt [€/kg]
        self.C_mat_FL = 2   # Flange [€/kg]
        self.C_mat_S = 2    # Shell [€/kg]
        self.C_mat_SP = 2   # Spacer [€/kg]
        self.C_mat_T = 3    # Tube [€/kg]
        self.C_mat_TR = 2   # Tie rod [€/kg]
        self.C_mat_TS = 3.5 # Tube sheet [€/kg]

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

        # Equipment investment costs [€]
        self.I_b = 20000    # Beveling machine [€]
        self.I_c = 120000   # Cutting machine [€]
        self.I_cw = 30000   # Channel welding [€]
        self.I_d = 50000    # Drilling machine [€]
        self.I_e = 5000     # Equipment [€]
        self.I_r = 100000   # Rolling machine [€]
        self.I_w = 75000    # Welding machine [€]

        # Dimensional parameters (length or spacing)
        self.L_l = 3        # Length of the lead [mm]
        self.L_ot = 5       # Overtravel distance [mm]
        self.L_Pst = 6      # Plate standard length [m]
        self.L_pt = 5       # Pretravel [mm]
        
        # Worker(s) per operation [-]
        self.m_b = 1    # Beveling machine [-]
        self.m_c = 1    # Cutting machine [-]
        self.m_cw = 1   # Channel welding [-]
        self.m_d = 1    # Drilling machine [-]]
        self.m_e = 1    # Equipment [-]
        self.m_r = 1    # Rolling machine [-]
        self.m_w = 1    # Welding machine [-]
        
        # Labor rate [€/h]
        self.LR = 22 #[€/h]
        
        # Amortization time spans [yr]
        self.NYr_b = 5      # Beveling [yr]
        self.NYr_c = 5      # Channel [yr]
        self.NYr_cw = 10    # Channel welding [yr]
        self.NYr_d = 5      # Drilling [yr]
        self.NYr_e = 5      # Equipment [yr]
        self.NYr_r = 5      # Rolling [yr]
        self.NYr_w = 5      # Welding [yr]

        # Operating factor [-]
        self.OF = 0.5       # Operating factor [-]

        # Power consumption [kW]
        self.P_b = 15        # Beveling [kW]
        self.P_c = 100       # Cutting [kW]
        self.P_cw = 10       # Channel welding [kW]
        self.P_w = 100       # Welding [kW]
        self.P_d = 10        # Drilling [kW]
        self.P_e = 5         # Equipment [kW]
        self.P_r = 100       # Rolling [kW]

        # Setup time [s/item]
        self.ST_e = 15      # Tube [s/item]
        self.ST_in_bt = 30  # Bolt [s/item]
        self.ST_in_sp = 15  # Spacer [s/item]
        self.ST_in_T = 3    # Tube [s/item]
        self.ST_in_td = 3   # Tie rod [s/item]

        self.ST_in = {
            'Tube' : {'Assembly' : 3},
            'Bolts' : {'Assembly' : 30},
            'Spacer' : {'Assembly' : 15},
            'Tie rod' : {'Assembly' : 3},            
            }

        # Load/unload time [s]
        self.T_LU = {
            'Baffle' : {'Cutting' : 30, 'Drilling' : 240},
            'Channel' : {'Cutting' : 40, 'Welding' : 180},
            'Tubesheet' : {'Drilling' : 120},
            'Rolled shell' : {'Rolling' : 120, 'Welding' : 300},
            'Standard Tube Shell' : {'Welding' : 300},
            }

        # Setup time [min]
        self.T_SU = {
            'Baffle' : {'Cutting' : 10, 'Drilling' : 1},
            'Channel' : {'Cutting' : 15, 'Welding' : 5},
            'Tubesheet' : {'Drilling' : 1},
            'Rolled shell' : {'Rolling' : 10, 'Welding' : 25},
            'Standard Tube Shell' : {'Welding' : 25},
            }

        # Processing speeds [m/min]
        self.v_b = 3                   # Beveling [m/min]
        self.v_c = 1                   # Cutting [m/min] assumed to be equal to 1 (to be better assessed)
        self.v_cw_mm_s = 1.5             # Channel welding [mm/s]
        self.v_cw=self.v_cw_mm_s*0.06    # Channel welding in [m/min]
        self.v_d = 0.3  # self.v_d = 3   # Drilling [m/min] assumed to be equal to 3 (to be better assessed)
        self.v_r = 0.2                  # Rolling [m/min]
        self.v_w = 0.2                     # Welding [m/min]
        
        # Electrical parameters
        self.WC = 150       # Welding current [A]
        self.WE = 0.9       # Welder electric efficiency [-]
        self.WFR = 2.5      # Wire feed rate [m/min]
        self.WP_Pst = 1.5   # Welding pipe [m]
        self.WP_pk = 5      # Welding plate thickness [cm]
        self.WV = 20        # Welding voltage [V]
        
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
        
        """
        Geometry parameters deduction from the defined inputs
        """
        # Calculate number of baffles
        self.N_B = math.floor(self.L_T / self.L_B) - 1
        
        # Calculate k4 = cos^-1((0.5-Bc)/0.5)
        self.k4 = math.arccos((0.5 - self.B_c) / 0.5)
        
        # Manufacturing consideration regarding the shell
        self.shell_type = "Rolled Shell" if self.D_S_i > 0.6 else "Standard Tube Shell"  # Shell is rolled if Ds > 600mm
        
        # Manufacturing consideration regarding the shell
        self.channel_type = "Rolled Channel" if self.D_S_i > 0.6 else "Standard Channel"  # Shell is rolled if Ds > 600mm
        
        # Tube pitch 
        self.pitch_T = pitch_r * self.D_T_o  # Tube pitch (m)
   
        self.operations = {
            'Beveling':{
                'investment': self.params.I_b,
                'life': self.params.NYr_b,
                'power': self.params.P_b,
                'workers': self.params.m_b,
                'aux_cost': self.params.C_AUX_b,
                'velocity': self.params.v_b,
                
            },
            'Cutting': {
                'investment': self.params.I_c,
                'life': self.params.NYr_c,
                'power': self.params.P_c,
                'workers': self.params.m_c,
                'aux_cost': self.params.C_AUX_c,
                'velocity': self.params.v_c,
                
            },
            'Welding': {
                'investment': self.params.I_w,
                'life': self.params.NYr_w,
                'power': self.params.P_w,
                'workers': self.params.m_w,
                'aux_cost': self.params.C_AUX_w,
                'velocity': self.params.v_w,
                
            },
            'Drilling': {
                'investment': self.params.I_d,
                'life': self.params.NYr_d,
                'power': self.params.P_d,
                'workers': self.params.m_d,
                'aux_cost': self.params.C_AUX_d,
                'velocity': self.params.v_d,
                
            },
            'Rolling':{
                'investment': self.params.I_r,
                'life': self.params.NYr_r,
                'power': self.params.P_r,
                'workers': self.params.m_r,
                'aux_cost': self.params.C_AUX_r,
                'velocity': self.params.v_r,
                
            }
        } 
   
    def calculate_volume(self, component):
        """Calculate volume for each heat exchanger subassembly based on Appendix B formulas"""
        if component == 'Rolled Shell' or component == 'Standard Tube Shell' :
            return math.pi * ((self.D_S_i + 2*self.t_S)**2 - self.D_S_i**2) / 4 * self.L_T  # Shell volume
        elif component == 'Baffle':
            return self.D_S_i**2 * (math.pi / 4 * (1 - 1 / (math.pi) * self.k4) + 0.5 * math.sin(self.k4) * (0.5 - self.B_c)) * self.t_B * self.N_B  # Baffle volume
        elif component == 'Tubesheet':
            return math.pi * ((self.D_S_i + 2*self.D_r)**2 - self.D_S_i**2) / 4 * self.t_TS * self.N_TS  # Tube-sheet volume
        elif component == 'Tube':
            return math.pi * (self.D_T_o**2 - self.D_T_i**2) / 4 * self.L_T * self.N_T  # Tube volume
        elif component == 'Rolled Channel' or component == 'Standard Channel':
            return math.pi * ((self.D_S_i + 2*self.t_S)**2 - self.D_S_i**2) / 4 * self.L_CH * self.N_CH  # Channel volume
        elif component == 'Removable Cover':
            return math.pi * ((self.D_S_i + 2*self.D_r)**2) / 4 * self.t_RC * self.N_CH  # Removable cover volume
        elif component == 'Flange':
            return math.pi * ((self.D_S_i + 2*self.D_r)**2 - self.D_S_i**2) / 4 * self.t_FL * self.N_FL  # Flange volume (Equation 51)
        elif component == 'Tie rods':
            return self.N_TR * math.pi * (self.D_TR**2 / 4) * self.L_TR  # Tie rod volume (Equation 51 bis)
        elif component == 'Spacer':
            return self.N_B * math.pi * ((self.D_SP_e**2 - self.D_SP_i**2) / 4) * self.L_B * self.N_TR  # Spacer volume (Equations 5 ter)
        elif component == 'Bolts':
            return self.N_Bt * 3E-6   # Bolts volume (Equations 5 ter)
        else: 
            raise ValueError(f"Unknown component '{component}'")
                

    def calculate_material_cost(self, component):
        """Calculate material cost based on volume and density"""
        volume = self.calculate_volume(component)
        density = 7850  # kg/m³ for steel
        material_costs = {
            'Rolled Shell': self.params.C_mat_S,
            'Standard Tube Shell': self.params.C_mat_S,
            'Baffle': self.params.C_mat_B,
            'Tubesheet': self.params.C_mat_TS,
            'Tube': self.params.C_mat_T,
            'Rolled Channel': self.params.C_mat_B, #Assumed to be equal to baffles
            'Standard Channel': self.params.C_mat_B, #Assumed to be equal to baffles
            'Flange': self.params.C_mat_FL,
            'Removable Cover': self.params.C_mat_B, #Assumed to be equal to baffles
            'Tie rods': self.params.C_mat_TR,
            'Spacer': self.params.C_mat_SP,
            'Bolts': self.params.C_mat_bt,
        }
        return volume * density * material_costs[component]

    def compute_processing_length(self, operation, component):
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
        
        self.processing_lengths = {
            "Rolled Shell": {
                "Beveling": 2*(math.pi * self.D_S_i + self.W_RP)*self.N_RP,
                "Cutting": 2*(math.pi * self.D_S_i + self.W_RP)*self.N_RP,
                "Welding": self.W_RP*self.N_RP+math.pi * self.D_S_i * (1 + 1/self.N_RP)*self.N_RP, #Sum of longitudinal and circumferntial welding
                "Rolling": math.pi * self.D_S_i*self.N_RP,
                
            },
            "Standard Tube Shell": {
                "Cutting": math.pi * self.D_S_o*self.N_T_tr,
                "Welding": math.pi * self.D_S_i * (1 - 1/self.N_T_tr)*self.N_T_tr,
                
            },
            "Baffle": {
                "Beveling": self.D_S_i * ((math.pi - self.k4) + math.sin(self.k4))*self.N_B,
                "Cutting": self.D_S_i * ((math.pi - self.k4) + math.sin(self.k4))*self.N_B,
                "Drilling": ((1 - self.k4/math.pi) + (2/math.pi) * math.sin(self.k4) * (1/2 - self.B_c)) * self.t_B * self.N_T*self.N_B,
                
            },
            "Flange": {
                "Beveling": 2 * math.pi * self.D_S_i * (1 + self.D_r)*self.N_FL,
                "Cutting": 2 * math.pi * self.D_S_i * (1 + self.D_r)*self.N_FL,
                "Drilling": (math.pi * self.D_S_i * (1 + self.D_r) / self.bd) * self.t_FL*self.N_FL,
                
            },
            "Tubesheet": {
                "Beveling": math.pi * self.D_S_i * (1 + 2 * self.D_r)*self.N_TS,
                "Cutting": math.pi * self.D_S_i * (1 + 2 * self.D_r)*self.N_TS,
                "Drilling": (self.N_T + self.D_S_i * (1 + self.D_r) / self.bd) * self.t_TS*self.N_TS,
               
            },
            "Rolled Channel": {
                "Beveling": 2 * (math.pi * self.D_S_i + self.W_CH)*self.N_CH,
                "Cutting": 2 * (math.pi * self.D_S_i + self.W_CH)*self.N_CH,
                "Welding": self.W_CH*self.N_CH + 2 * math.pi * self.D_S_o*self.N_CH,
                "Rolling": math.pi * self.D_S_i*self.N_CH,
                
            },
            "Standard Channel": {
                "Beveling": math.pi * self.D_S_o * self.N_CH,
                "Cutting": math.pi * self.D_S_o * self.N_CH,
                "Welding": 2 * math.pi * self.D_S_o* self.N_CH, #Only circumferential
            },
            "Removable Cover": {
                "Beveling": math.pi * self.D_S_i * (1 + 2 * self.D_r)* self.N_RC,
                "Cutting": math.pi * self.D_S_i * (1 + 2 * self.D_r)* self.N_RC,
                "Drilling": (self.D_S_i * (1 + self.D_r) / self.bd) * self.t_RC * self.N_RC,
               
            },
            "Tube": {
                "Cutting": math.pi * self.D_T_o* self.N_T,
                "Welding": math.pi * self.D_T_o * self.N_T_tr* self.N_T,#Only circumferential
               
            }
        }
        
        self.batches = {
            'Rolled Shell' : self.r,
            'Baffle' : self.N_B,
            'Tubesheet' : self.N_TS,
            'Tube' : self.N_T,
            'Rolled channel' : self.N_CH,
            'Removable Cover' : self.N_RC,
            'Flange' : self.N_FL,
            'Tie rods' : self.N_TR,
            }
        
        result = self.processing_lengths.get(component, {}).get(operation, 0)

        return result


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

    def calculate_operation_cost(self, operation, length, component):
        """Calculate cost for a manufacturing operation"""
        
        if operation not in self.operations:
            return {'total': 0, 'depreciation': 0, 'labor': 0, 'energy': 0, 'auxiliary': 0}
    
        op = self.operations[operation]
        process_time = length / (op['velocity'] * 60)  # In hours
        
        if component in self.params.T_SU:
            if operation in self.params.T_SU[component]:
                process_setup_time = process_time + self.params.T_SU[component][operation]/60
            else:
                process_setup_time = process_time
        else:
            process_setup_time = process_time

        if component in self.params.T_LU:
            if operation in self.params.T_LU[component]:
                process_total_time = process_setup_time + self.params.T_LU[component][operation]/3600 * self.batches[component]
            else:
                process_total_time = process_setup_time
        else:
            process_total_time = process_setup_time
        
        print(f"{operation} : {process_total_time}")
        
        # Calculate costs
        depreciation = self.hourly_depreciation_cost(op['investment'], op['life']) * process_total_time
        labor = self.params.LR * op['workers'] * process_total_time
        energy = self.calculate_energy_cost(op['power']) * process_time
        auxiliary = op['aux_cost'] * process_time
        
        total = depreciation + labor + energy + auxiliary

        oper_costs =  {
            'total': total,
            'depreciation': depreciation,
            'labor': labor,
            'energy': energy,
            'auxiliary': auxiliary, 
        }       

        print(operation)        
        print(oper_costs)
        return oper_costs


    def calculate_total_cost(self):
        """Calculate total cost """
        total_material_cost = sum(self.calculate_material_cost(component) for component in [self.shell_type, 'Baffle', 'Tubesheet', 'Tube', self.channel_type, 'Removable Cover', 'Flange','Tie rods', 'Spacer', 'Bolts'])
        total_manufacturing_cost = sum(self.calculate_operation_cost(op, self.compute_processing_length(op, component), component)['total']
                                       for component in [self.shell_type, 'Baffle', 'Tubesheet', 'Tube', self.channel_type, 'Removable Cover', 'Flange','Tie rods', 'Spacer', 'Bolts']
                                       for op in ['Beveling','Cutting', 'Welding', 'Drilling', 'Rolling'])
        total_cost = (total_material_cost + total_manufacturing_cost)
        return total_cost


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

            print("---------------")
            
            print(component)
            
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
    t_B = 0.032 # m
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

    Cost_Breakdown=calculator.print_cost_breakdown()
    #Test=calculator.compute_processing_length("Drilling","Baffle")
    
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


