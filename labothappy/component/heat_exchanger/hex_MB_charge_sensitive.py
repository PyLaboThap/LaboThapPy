
"""
Modification w/r to previous version:
    - Putting some order in the Objective Function "for" loops. Sparing some
    lines of code.
    - x_di_c correct calculation.
    - Implementation of various correlations for htc and DP
    - Implementation of total DP linearly interpolated on all discretizations
    - Implemenation of the code in the LaboThap Python Library
"""


"EXTERNAL IMPORTS"

import sys
import os
    
import correlations

import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
import copy

"INTERNAL IMPORTS"

from correlations.heat_exchanger.f_lmtd2 import f_lmtd2, F_shell_and_tube
from correlations.heat_exchanger.find_2P_boundaries import find_2P_boundaries

# HTC Correlations
from correlations.convection.plate_htc import han_BPHEX_DP, water_plate_HTC, martin_BPHEX_HTC, muley_manglik_BPHEX_HTC, han_boiling_BPHEX_HTC, han_cond_BPHEX_HTC, thonon_plate_HTC, kumar_plate_HTC, martin_holger_plate_HTC, amalfi_plate_HTC, shah_condensation_plate_HTC
from correlations.convection.pipe_htc import gnielinski_pipe_htc, boiling_curve, horizontal_tube_internal_condensation, horizontal_flow_boiling, flow_boiling_gungor_winterton, Liu_sCO2, Cheng_sCO2
from correlations.convection.shell_and_tube_htc import shell_bell_delaware_htc, shell_htc_kern
from correlations.convection.tube_bank_htc import ext_tube_film_condens
from correlations.convection.fins_htc import htc_tube_and_fins
from correlations.convection.printed_circuit_htc import PCHE_Lee, PCHE_conv

# DP Correlations 
from correlations.pressure_drop.shell_and_tube_DP import shell_DP_kern, shell_DP_donohue, shell_bell_delaware_DP
from correlations.pressure_drop.pipe_DP import gnielinski_pipe_DP , Churchill_DP, Choi_DP, Muller_Steinhagen_Heck_DP, Cheng_CO2_DP, Darcy_Weisbach

# Fluid Correlations
from correlations.properties.thermal_conductivity import conducticity_R1233zd

# Phase related correlations
from correlations.properties.void_fraction import void_fraction
from correlations.heat_exchanger.kim_dry_out_incipience import kim_dry_out_incipience

# Connectors
from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

# Component base frame
from component.base_component import BaseComponent

#%%
# Set to True to enable some debugging output to screen
debug = False

HTC_correlations = {
    "han_BPHEX_DP": han_BPHEX_DP,
    "water_plate_HTC": water_plate_HTC,
    "martin_BPHEX_HTC": martin_BPHEX_HTC,
    "muley_manglik_BPHEX_HTC": muley_manglik_BPHEX_HTC,
    "han_boiling_BPHEX_HTC": han_boiling_BPHEX_HTC,
    "han_cond_BPHEX_HTC": han_cond_BPHEX_HTC,
#    "thonon_plate_HTC": thonon_plate_HTC,
    "gnielinski_pipe_htc": gnielinski_pipe_htc,
    "shell_bell_delaware_htc": shell_bell_delaware_htc,
    "ext_tube_film_condens": ext_tube_film_condens,
    "htc_tube_and_fins": htc_tube_and_fins
}

def propsfluid_AS(T_mean, P_mean, T_wall, fluid, incompr_flag, AS):
    
    AS.update(CP.PT_INPUTS, P_mean, T_mean)    
    mu = AS.viscosity()
    cp = AS.cpmass()
    
    if fluid == 'R1233zd(E)':
        k = conducticity_R1233zd(T_mean, P_mean)
        Pr = mu * cp / k
    else:
        k = AS.conductivity()    
        Pr = AS.Prandtl()

    AS.update(CP.PT_INPUTS, P_mean, T_wall)    

    mu_wall = AS.viscosity()
    mu_rat = mu/mu_wall
    cp_wall = AS.cpmass()

    if fluid == 'R1233zd(E)':
        k_wall = conducticity_R1233zd(T_mean, P_mean)
        Pr_wall = mu_wall * cp_wall / k_wall
    else:
        k_wall = AS.conductivity()    
        Pr_wall = AS.Prandtl()
    
    return mu, Pr, k, mu_wall, mu_rat, Pr_wall, 0

class HeatExchangerMB(BaseComponent):
    """
        Component: Heat Exchanger

        Model: The model is based on the the work of Ian Bell "A generalized moving-boundary algorithm to predict the heat transfer
                rate of counterflow heat exchangers for any phase configuration". 

        **Descritpion**:

                It is a moving-boundary model allowing for the precise estimation of its performance, pressure drops and fluid charges inside it (using a void fraction correlation for two-phase flow conditions). 
                For now, the method has been adapted for many geometries : 
                - Brazed-Plate heat exchangers
                - Shell-and-Tube heat exchangers
                - Cross-CounterFlow tube and fin heat exchangers

                The precision of the model and the geometry dependent charge estimation requires to know the geometry. 
                The required parameters are adapted to the used heat transfer and pressure drop correlations that can be chosen from the ones implemented.

        **Assumptions**:

            - Steady-state operation.
            - Pressure drop proportional to the heat transfer area for each discretization.
            - Volume used for the charge computation are also assumed proportional to the heat transfer area for each discretization.
            - No loss to the ambient is considered.

        **Connectors**:

            su_H (MassConnector): Mass connector for the hot suction side.
            su_C (MassConnector): Mass connector for the cold suction side.

            ex_H (MassConnector): Mass connector for the hot exhaust side.
            ex_C (MassConnector): Mass connector for the cold exhaust side.

            Q_dot (HeatConnector): Heat connector for the heat transfer between the fluids

        **Parameters**:

            Parameters can be differentiated between geometrical parameters and simulation parameters. 

            Simulation required parameters are: 

            Flow_Type : Flow configuration of the fluid ('CounterFlow', 'CrossFlow', 'Shell&Tube', 'ParallelFlow') [-]

            htc_type  : Heat Transfer coefficient type ('User-Defined' or 'Correlation') [-]. If the 'User-Defined' option is chosen,
            values for the different heat transfer coefficients can be manually passed by the user for different state conditions: 
            'Liquid', 'Vapor', 'Two-Phase', 'Vapor-wet', 'Dryout', 'Transcritical'. If the 'Correlation' option is chosen, then the names
            of implemented correlations shall be passed. For 1 phase heat transfer, 2 phases heat transfer and for pressure drops. 

            H_DP_ON   : Hot side pressure drop considered or not [-]
            
            C_DP_ON   : Cold side pressure drop considered or not [-]
            
            n_disc    : number of heat exchanger discretizations [-]
        
            Geometry required parameters depend on the type of heat exchanger and on especially the used correlations. 
            More information can be found in the documentation of the correlations themselves. 

        **Inputs**:

            P_su_H: Hot suction side pressure. [Pa]

            h_su_H: Hot suction side enthalpy. [J/kg]

            fluid_H: Hot suction side fluid. [-]

            m_dot_H: Hot suction side mass flowrate. [kg/s]

            P_su_C: Cold suction side pressure. [Pa]

            h_su_C: Cold suction side enthalpy. [J/kg]

            fluid_C: Cold suction side fluid. [-]

            m_dot_C: Cold suction side mass flowrate. [kg/s]

        **Ouputs**:

            h_ex_H: Hot exhaust side specific enthalpy. [J/kg]

            P_ex_H: Hot exhaust side pressure. [Pa]

            h_ex_C: Cold exhaust side specific enthalpy. [J/kg]

            P_ex_C: Cold exhaust side pressure. [Pa]

            Q: Heat transfer rate [W]

            M_H : Hot fluid charge [kg]

            M_C : Cold fluid charge [kg]
    """
    "Creation of oclass hot and cold otherwise, I had some problems with shared references with creating several heat exchanger objects"
    class H:
        def __init__(self):
            self.Correlation_1phase = None
            self.Correlation_2phase = None
            self.HeatExchange_Correlation = None  
            self.PressureDrop_Correlation = None
            self.h_liq = None
            self.h_vap = None
            self.h_twophase = None
            self.h_vapwet = None
            self.h_tpdryout = None
            self.h_transcrit = None
            self.f_dp = None

    class C:
        def __init__(self):
            self.Correlation_1phase = None
            self.Correlation_2phase = None
            self.HeatExchange_Correlation = None
            self.PressureDrop_Correlation = None
            self.h_liq = None
            self.h_vap = None
            self.h_twophase = None
            self.h_vapwet = None
            self.h_tpdryout = None
            self.h_transcrit = None
            self.f_dp = None
    
    def __init__(self, HTX_Type):
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
                
        self.Q_HX = HeatConnector()
        self.F_fun = None
        self.w = [None]
        
        self.AS_C = None
        self.AS_H = None
        
        if HTX_Type == 'Plate' or HTX_Type == 'Shell&Tube' or HTX_Type == 'Tube&Fins' or HTX_Type == 'PCHE':
            self.HTX_Type = HTX_Type
        else:
            print("Heat exchanger types implemented for this model are : 'Plate', 'Shell&Tube', 'Tube&Fins'.")

        self.H = self.H()
        self.C = self.C()
        
    #%% INPUTS AND PARAMETERS RELATED METHODS
    
    def get_required_inputs(self): # Used in check_calculable to see if all of the required inputs are set

        self.sync_inputs()
        # Return a list of required inputs
        return['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']

    def get_required_parameters(self):
        """
        General Parameters : 
            
            - Flow_Type : Flow configuration of the fluid ('CounterFlow', 'CrossFlow', 'Shell&Tube', 'ParallelFlow')
            - htc_type  : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - H_DP_ON   : Hot side pressure drop considered or not
            - C_DP_ON   : Cold side pressure drop considered or not
            - n_disc    : number of discretizations
        
        Geometry Parameters depend on specific geometry python files.
            
        """
        general_parameters = ['Flow_Type', 'htc_type', 'H_DP_ON', 'C_DP_ON','n_disc']
        geometry_parameters = []

        if self.HTX_Type == 'Plate':     
            geometry_parameters = ['A_c', 'A_h', 'h', 'l', 'l_v',
                                   'C_CS', 'C_Dh', 'C_V_tot', 'C_canal_t', 'C_n_canals',
                                   'H_CS', 'H_Dh', 'H_V_tot', 'H_canal_t', 'H_n_canals',
                                   'casing_t', 'chevron_angle', 'fooling', 
                                   'n_plates', 'plate_cond', 'plate_pitch_co', 't_plates', 'w']

        elif self.HTX_Type == 'Shell&Tube':
                
            if self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC" or self.C.Correlation_1phase == "Shell_Bell_Delaware_HTC":

                geometry_parameters = ['A_eff', 'Baffle_cut', 'D_OTL', 'N_strips', 'S_V_tot',
                                    'Shell_ID', 'T_V_tot', 'Tube_L', 'Tube_OD', 'Tube_pass',
                                    'Tube_t', 'Tubesheet_t', 'central_spacing', 'clear_BS', 'clear_TB',
                                    'cross_passes', 'foul_s', 'foul_t', 'inlet_spacing', 'n_series',
                                    'n_tubes', 'outlet_spacing', 'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']

            if self.H.Correlation_1phase == "Shell_Kern_HTC" or self.C.Correlation_1phase == "Shell_Kern_HTC":

                geometry_parameters = ['A_eff', 'Baffle_cut', 'S_V_tot', 'Shell_ID', 'T_V_tot',
                                    'Tube_L', 'Tube_OD', 'Tube_pass','Tube_t', 'central_spacing',
                                    'cross_passes', 'foul_s', 'foul_t', 'n_tubes', 'pitch_ratio', 
                                    'tube_cond', 'tube_layout', 'Shell_Side']
        
        elif self.HTX_Type == 'Tube&Fins':
            
            geometry_parameters = ['A_finned', 'A_flow', 'A_in_tot', 'A_out_tot', 'A_unfinned',
                                   'B_V_tot', 'Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type',
                                   'Finned_tube_flag', 'L', 'T_V_tot', 'Tube_L', 'Tube_OD',
                                   'Tube_cond', 'Tube_t', 'fouling', 'h', 'k_fin',
                                   'Tube_pass', 'n_rows', 'n_tubes', 'pitch', 'pitch_ratio', 'tube_arrang',
                                   'w','Fin_Side']
                
        elif self.HTX_Type == 'PCHE':
            
            geometry_parameters = ['alpha', 'D_c', 'H_V_tot', 'C_V_tot', 'k_cond', 'L_c', 'N_c', 'N_p', 'R_p', 't_2', 't_3']
            
        return general_parameters + geometry_parameters
    
    #%% PRINT RELATED METHODS

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

    #%% CORRELATION RELATED METHODS

    def set_htc(self, htc_type = "Correlation", Corr_H = None, Corr_C = None, UD_H_HTC = None, UD_C_HTC = None):
        """
        General Parameters : 
            
            - htc_type : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - Corr_H   : Correlations for hot side
            - Corr_C   : Correlations for cold side
            - UD_H_HTC : User-Defined HTC for hot side
            - UD_C_HTC : User-Defined HTC for cold side
            
        """
        
        self.params['htc_type'] = htc_type
        # self.check_calculable()
        
        if htc_type == "User-Defined":
            
            self.H.HeatExchange_Correlation = "User-Defined"
            self.C.HeatExchange_Correlation = "User-Defined"
            
            # User-Defined Heat Transfer Coefficients (hot):
            self.H.h_liq = UD_H_HTC['Liquid']
            self.H.h_vap = UD_H_HTC['Vapor']
            self.H.h_twophase = UD_H_HTC['Two-Phase']
            self.H.h_vapwet = UD_H_HTC['Vapor-wet']
            self.H.h_tpdryout = UD_H_HTC['Dryout']
            self.H.h_transcrit = UD_H_HTC['Transcritical']
        
            # User-Defined Heat Transfer Coefficients (cold):
            self.C.h_liq = UD_C_HTC['Liquid']
            self.C.h_vap = UD_C_HTC['Vapor']
            self.C.h_twophase = UD_C_HTC['Two-Phase']
            self.C.h_vapwet = UD_C_HTC['Vapor-wet']
            self.C.h_tpdryout = UD_C_HTC['Dryout']
            self.C.h_transcrit = UD_C_HTC['Transcritical']
        
        else: 
            # Type 
            self.H.HeatExchange_Correlation = "Correlation"
            self.C.HeatExchange_Correlation = "Correlation"
            
            if self.HTX_Type == 'Plate' or self.HTX_Type == 'Shell&Tube' or self.HTX_Type == 'Tube&Fins' or self.HTX_Type == 'PCHE':
                
                self.H.Correlation_1phase = Corr_H["1P"]
                if "2P" in Corr_H:
                    self.H.Correlation_2phase = Corr_H["2P"]
                else:
                    self.H.Correlation_2phase = None
                    
                if "SC" in Corr_H:
                    self.H.Correlation_TC = Corr_H["SC"]
                else:
                    self.H.Correlation_TC = None
                
                self.C.Correlation_1phase = Corr_C["1P"]
                if "2P" in Corr_C:
                    self.C.Correlation_2phase = Corr_C["2P"]
                else:
                    self.C.Correlation_2phase = None

                if "SC" in Corr_C:
                    self.C.Correlation_TC = Corr_C["SC"]
                else:
                    self.C.Correlation_TC = None

            if self.C.Correlation_2phase == "Boiling_curve": # Compute the fluid boiling curve beforehand
                try:
                    self.AS_C = CP.AbstractState("BICUBIC&HEOS", self.su_C.fluid)   
                    self.AS_C.update(CP.PQ_INPUTS, self.su_C.p, 0)
                    
                    T_sat = self.AS_C.T()
                    (h_boil, DT_vect) = boiling_curve(self.params['Tube_OD'], self.su_C.fluid, T_sat, self.su_C.p)
                    self.C_f_boiling = interp1d(DT_vect,h_boil)
                except:
                    self.C_f_boiling = interp1d([0,10000],[20000,20000])

    def set_DP(self, DP_type = None, Corr_H = None, Corr_C = None, UD_H_DP = None, UD_C_DP = None):
        """
        General Parameters : 
            
            - htc_type : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - Corr_H   : Correlations for hot side
            - Corr_C   : Correlations for cold side
            - UD_H_HTC : User-Defined HTC for hot side
            - UD_C_HTC : User-Defined HTC for cold side
            
        """

        self.params['DP_type'] = DP_type

        # self.check_calculable()
        
        if self.params['H_DP_ON']:
            if self.params['DP_type'] == "User-Defined":
                # !!! User-Defined Pressure Drop Value (hot side):
                self.H.PressureDrop_Correlation = "User-Defined"            
                self.H.DP_val = UD_H_DP
            elif self.params['DP_type'] == "Correlation":
                if self.HTX_Type == 'Plate' or self.HTX_Type == 'Shell&Tube' or self.HTX_Type == 'Tube&Fins' or self.HTX_Type == 'PCHE':
                    self.H.PressureDrop_Correlation = Corr_H
                else:
                    print("No Pressure drop correlation for that HTX Type.")
            else:
                self.H.f_dp = {"K": 14.14, "B": 1.892}
        else: 
            self.H.f_dp = {"K": 0, "B": 0}        

        if self.params['C_DP_ON']:
            if self.params['DP_type'] == "User-Defined":
                # !!! User-Defined Pressure Drop Value (cold side):
                self.C.PressureDrop_Correlation = "User-Defined"
                self.C.DP_val = UD_C_DP
            elif self.params['DP_type'] == "Correlation":
                if self.HTX_Type == 'Plate' or self.HTX_Type == 'Shell&Tube' or self.HTX_Type == 'Tube&Fins' or self.HTX_Type == 'PCHE':
                    self.C.PressureDrop_Correlation = Corr_C
                else:
                    print("No Pressure drop correlation for that HTX Type.")
            else:
                self.C.f_dp = {"K": 14.14, "B": 1.892}
        else: 
            self.C.f_dp = {"K": 0, "B": 0}        

    #%% PINCH AND DISCRETIZATION RELATED METHODS
            
    def external_pinching(self):
        "Determine the maximum heat transfer rate based on the external pinching analysis"
        
        "1) Set temperature bound values" # !!! Find out why      
        T_hmin = 218 
        T_cmax = 273.15+260 # 481 # 
        
        "2) Hot fluid side pinch"
        
        # Computation of outlet lowest possible enthalpy of the hot fluid using either the entrance cold fluid inlet temperature, or the arbitrary minimal
        self.AS_H.update(CP.PT_INPUTS, self.p_hi, max(self.T_ci, T_hmin))
        h_ho_ideal = self.AS_H.hmass() # Equation 5 (Bell et al. 2015)
        
        # Q_max computation
        Q_dot_maxh = self.mdot_h*(self.h_hi-h_ho_ideal) # Equation 4 (Bell et al. 2015)
        
        if debug:
            print("Inlet Hot Enthalpy:", self.h_hi, "\n")
            print("Minimum Outlet Hot Enthalpy:", h_ho_ideal,"\n")
            print("mdot_h = ", self.mdot_h, "\n")
            print("Qmaxh =", Q_dot_maxh, "\n")

        "3) Cold fluid side pinch"
                
        # Computation of highest possible outlet enthalpy of the hot fluid using either the entrance hot fluid inlet temperature, or the arbitrary maximal
        self.AS_C.update(CP.PT_INPUTS, self.p_ci, min(self.T_hi, T_cmax))
        h_co_ideal = self.AS_C.hmass() # Equation  (Bell et al. 2015)
        
        # Q_max computation
        Q_dot_maxc = self.mdot_c*(h_co_ideal-self.h_ci) # Equation 6 (Bell et al. 2015)

        if debug:
            print("Inlet Cold Enthalpy:", self.h_ci, "\n")
            print("Maximum Outlet Cold Enthalpy:", h_co_ideal,"\n")
            print("mdot_c = ", self.mdot_c, "\n")
            print("Qmaxc =", Q_dot_maxc, "\n")
            
        "4) Final pinch and cell boundaries computation"
        
        Q_dot_max = min(Q_dot_maxh, Q_dot_maxc)
        
        if debug:
            print('Qmax (external pinching) is', Q_dot_max)

        self.calculate_cell_boundaries(Q_dot_max) # call calculate_cell_boundaries procedure

        return Q_dot_max

    def calculate_cell_boundaries(self, Q):
        """ Calculate the cell boundaries for each fluid """
        
        "1) Re-calculate the outlet enthalpies of each stream"
                
        self.h_co = self.h_ci + Q/self.mdot_c
        self.h_ho = self.h_hi - Q/self.mdot_h

        "2) Calculate the dew and bubble pressures and enthalpies by accounting for pressure drops"
        
        "2.1) For the cold fluid"
        
        if not self.Transcritical_c: # if not in transcritical
            if (not (self.C.f_dp == {"K": 0, "B": 0}) or self.DP_c < 1e-2):
                # If no pressure drop is calculated or if the pressure drop is neglectable, assign saturation conditions to the ideal case:
                self.p_cdew = self.p_ci
                self.p_cbubble = self.p_co
                self.T_cbubble = self.T_cbubble_ideal
                self.T_cdew = self.T_cdew_ideal
                self.h_cbubble = self.h_cbubble_ideal
                self.h_cdew = self.h_cdew_ideal
            else:
                h_cbubble, h_cdew, p_cbubble, p_cdew, _, _, = find_2P_boundaries(self.C_su.fluid, self.h_ci, self.h_co, self.p_ci, self.p_co)
                self.p_cdew = p_cdew
                self.p_cbubble = p_cbubble
                self.h_cdew = h_cdew
                self.h_cbubble = h_cbubble

                self.AS_C.update(CP.PQ_INPUTS, self.p_cdew, 1)
                self.T_cdew = self.AS_C.T()

                self.AS_C.update(CP.PQ_INPUTS, self.p_cbubble, 0)
                self.T_cbubble = self.AS_C.T()
                
        "2.2) For the hot fluid"
        
        if not self.Transcritical_h:
            if (not (self.H.f_dp == {"K": 0, "B": 0}) or self.DP_h < 1e-2) and not self.h_incomp_flag:
                #If no pressure drop is calculated or if the pressure drop is neglectable, assign saturation conditions to the ideal case:
                self.p_hdew = self.p_hi
                self.p_hbubble = self.p_ho    
                self.T_hbubble = self.T_hbubble_ideal
                self.T_hdew = self.T_hdew_ideal
                self.h_hbubble = self.h_hbubble_ideal
                self.h_hdew = self.h_hdew_ideal
            elif not self.h_incomp_flag:
                h_hbubble, h_hdew, p_hdew, p_hbubble, _, _, = find_2P_boundaries(self.H_su.fluid, self.h_hi, self.h_ho, self.p_hi, self.p_ho)
                self.p_hdew = p_hdew
                self.p_hbubble = p_hbubble
                self.h_hdew = h_hdew
                self.h_hbubble = h_hbubble
                
                self.AS_H.update(CP.PQ_INPUTS, self.p_hdew, 1)
                self.T_hdew = self.AS_H.T()

                self.AS_H.update(CP.PQ_INPUTS, self.p_hbubble, 0)
                self.T_hbubble = self.AS_H.T()

        "3) Discretization of the heat exchanger, user defined"
        # The minimum number of cells desired is n_disc
        # n_disc is increased by 1. Now it is the minimum number of cells boundaries.
        n_disc = self.params['n_disc'] + 1
        
        # Force the minimum number of discretizations needed
        n_disc = max(n_disc, 2)
        
        # The original Bell2015 code uses: "if h_cdew is not None". I don't see when that could happen without an error before.
        "3.1) Create a vector of enthalpies : initiates the cell boundaries. If phase change, dew and bubble point are added"
        
        # Cold side
        if not self.Transcritical_c:
            self.hvec_c = np.append(np.linspace(self.h_ci, self.h_co, n_disc), [self.h_cdew, self.h_cbubble])
        elif self.Transcritical_c:
            #In transcritical, just dont append the unexisting dew and bubble enthalpies
            self.hvec_c = np.linspace(self.h_ci, self.h_co, n_disc)

        self.hvec_c = np.sort(np.unique(self.hvec_c))
        # Ensure that the enthalpy boundaries are all enclosed by the inlet and outlet conditions
        self.hvec_c = [h for h in self.hvec_c if (h >= self.h_ci and h <= self.h_co)]
        
        # Hot side
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.hvec_h = np.append(np.linspace(self.h_hi, self.h_ho, n_disc), [self.h_hdew, self.h_hbubble])
        elif self.Transcritical_h or self.h_incomp_flag:
            #In transcritical, just dont append the unexisting dew and bubble enthalpies
            self.hvec_h = np.linspace(self.h_hi, self.h_ho, n_disc)
            
        self.hvec_h = np.sort(np.unique(self.hvec_h))
        # Ensure that the enthalpy boundaries are all enclosed by the inlet and outlet conditions
        self.hvec_h = [h for h in self.hvec_h if (h >= self.h_ho and h <= self.h_hi)]

        "4) Filling of the complementary cell boundaries"
        
        # Complementary cell boundaries serve to keep heat balances in check
        # Start at the first element in the vector
        k = 0
        
        EPS_REL = 1e-8     # relative tol
        EPS_ABS = 1e-3     # absolute tol in Watts (tune if needed)
        
        while k < len(self.hvec_c)-1 or k < len(self.hvec_h)-1:
            if len(self.hvec_c) == 2 and len(self.hvec_h) == 2: # If there is only one cell
                break

            # Determine which stream is the limiting next cell boundary
            Qcell_hk = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])
            Qcell_ck = self.mdot_c*(self.hvec_c[k+1]-self.hvec_c[k])

            if debug:
              print("k+1 Hot Enthalpy:", self.hvec_h[k+1], "\n")
              print("k Hot Enthalpy:", self.hvec_h[k],"\n")
              print("mdot_h = ", self.mdot_h, "\n")
              print("Qcell_hk =", Qcell_hk, "\n-------\n")
              print("k+1 Cold Enthalpy:", self.hvec_c[k+1], "\n")
              print("k Cold Enthalpy:", self.hvec_c[k],"\n")
              print("mdot_c = ", self.mdot_c, "\n")
              print("Qcell_hk =", Qcell_ck, "\n-------\n")
            
            diff = Qcell_hk - Qcell_ck
            tol  = max(EPS_ABS, EPS_REL * max(abs(Qcell_hk), abs(Qcell_ck)))
            
            if diff > tol:
                # hot has more heat -> add boundary on hot
                self.hvec_h.insert(k+1, self.hvec_h[k] + Qcell_ck / self.mdot_h)
            
            elif diff < -tol:
                # cold has more capacity -> add boundary on cold
                self.hvec_c.insert(k+1, self.hvec_c[k] + Qcell_hk / self.mdot_c)

            if debug:
                print(k,len(self.hvec_c),len(self.hvec_h),Qcell_hk, Qcell_ck)

            # Increment index
            k += 1

        if debug:
             print("Modified Length of hvec_c is ", len(self.hvec_c))
             print("Modified Length of hvec_h is ", len(self.hvec_h))

        assert(len(self.hvec_h) == len(self.hvec_c)) # Ensures that the two fluid streams have the same number of cell boundaries

        # Calculate the vectors of exchanged heat in each cell.
        self.Qvec_h = np.array([self.mdot_h*(self.hvec_h[i+1]-self.hvec_h[i]) for i in range(len(self.hvec_h)-1)])
        self.Qvec_c = np.array([self.mdot_c*(self.hvec_c[i+1]-self.hvec_c[i]) for i in range(len(self.hvec_c)-1)])

        if debug:
            if np.max(np.abs(self.Qvec_c/self.Qvec_h))<1e-5:
                print(self.Qvec_h, self.Qvec_c)
            
        "5) Computation of the thermal states at the cell boundaries"
        
        # Build the normalized enthalpy vectors (normalized by the total exchanged enthalpy)
        self.hnorm_h = (np.array(self.hvec_h)-self.hvec_h[0])/(self.hvec_h[-1]-self.hvec_h[0])
        self.hnorm_c = (np.array(self.hvec_c)-self.hvec_c[0])/(self.hvec_c[-1]-self.hvec_c[0])

        # Pressure distribution accross the HX. Method by RDickes (not validated).
        # Discretization of the pressure as a linear interpolation on enthalpy between in and out conditions # !!! Is this good ? 
        self.pvec_c = (1-self.hnorm_c)*self.p_ci + self.hnorm_c*self.p_co
        self.pvec_h = (1-self.hnorm_h)*self.p_ho + self.hnorm_h*self.p_hi

        self.Tvec_c = np.zeros(np.size(self.pvec_c))
        self.Tvec_h = np.zeros(np.size(self.pvec_h))
        
        self.svec_c = np.zeros(np.size(self.pvec_c))
        self.svec_h = np.zeros(np.size(self.pvec_h))
        
        self.LMTD = np.zeros(np.size(self.pvec_h)-1)
        self.F = np.zeros(np.size(self.pvec_h)-1)
        
        # Calculate the temperature and entropy at each cell boundary
        for i in range(len(self.hvec_c)):
            
            # try and except because some CoolProp convergence problems could occur when evaluating the propertie
            # Error Code : unable to solve 1phase PY flash with Tmin=179.699, Tmax=481.825 due to error: 
            # HSU_P_flash_singlephase_Brent could not find a solution because Hmolar [27458.3 J/mol] is above 
            # the maximum value of 27404.148193 J/mol : PropsSI("T","H",391244,"P",3005929,"Cyclopentane")
                        
            try:
                self.AS_C.update(CP.HmassP_INPUTS, self.hvec_c[i], self.pvec_c[i])
            except:
                self.AS_C.update(CP.HmassP_INPUTS, (self.hvec_c[i-1] + self.hvec_c[i])/2, (self.pvec_c[i-1] + self.pvec_c[i])/2)

            self.Tvec_c[i] = self.AS_C.T()
            self.svec_c[i] = self.AS_C.smass()

            try:
                self.AS_H.update(CP.HmassP_INPUTS, self.hvec_h[i], self.pvec_h[i])
            except:
                self.AS_H.update(CP.HmassP_INPUTS, (self.hvec_h[i-1] + self.hvec_h[i])/2, self.pvec_h[i])

            self.Tvec_h[i] = self.AS_H.T()
            self.svec_h[i] = self.AS_H.smass()

        "6) Vapour quality calculation if in the two-phase zone. If not, force, the value to -2 or 2"
        # Take Transcritical into account. in that Case, force x = 3
        
        "6.1) Hot Fluid"
        self.x_vec_h = []
        for i, h_h in enumerate(self.hvec_h):
            if not self.Transcritical_h and not self.h_incomp_flag:
                if h_h < self.h_hbubble:
                    self.x_vec_h.append(-2)        
                elif h_h > self.h_hdew:
                    self.x_vec_h.append(2)
                else:
                    self.x_vec_h.append(min(1, max(0, (h_h - self.h_hbubble)/(self.h_hdew - self.h_hbubble))))
            else:
                self.x_vec_h.append(3)
                
        "6.2) Cold Fluid"
        self.x_vec_c = []
        for i, h_c in enumerate(self.hvec_c):
            if not self.Transcritical_c:
                if h_c < self.h_cbubble:
                    self.x_vec_c.append(-2)        
                elif h_c > self.h_cdew:
                    self.x_vec_c.append(2)
                else:
                    self.x_vec_c.append(min(1, max(0, (h_c - self.h_cbubble)/(self.h_cdew - self.h_cbubble))))
            else:
                self.x_vec_c.append(3)

        "7) Vector of saturation temperatures at quality = 0.5, considering the pressure drop at each cell"
        
        self.Tvec_sat_pure_c = np.zeros(np.size(self.pvec_c))
        self.Tvec_sat_pure_h = np.zeros(np.size(self.pvec_h))

        for i in range(len(self.pvec_c)):
            if not self.Transcritical_c:
                self.AS_C.update(CP.PQ_INPUTS, self.pvec_c[i], 0.5)
                self.Tvec_sat_pure_c[i] = self.AS_C.T()

        for i in range(len(self.pvec_h)):    
            if not self.Transcritical_h and not self.h_incomp_flag:
                self.AS_H.update(CP.PQ_INPUTS, self.pvec_h[i], 0.5)
                self.Tvec_sat_pure_h[i] = self.AS_H.T()
            
        "8) Calculate pinch and heat exchanger border temperature deltas"
        
        self.DT_pinch = np.min((np.array(self.Tvec_h) - np.array(self.Tvec_c)))
        self.DT_ho_ci = self.Tvec_h[0] - self.Tvec_c[0]
        self.DT_hi_co = self.Tvec_h[-1] - self.Tvec_c[-1]

    def internal_pinching(self, stream):
        """
        Determine the maximum heat transfer rate based on the internal pinching analysis
        NB : external pinching analysis has already been done
        """      
        # Try to find the dew point enthalpy as one of the cell boundaries
        # that is not the inlet or outlet
        # Do this only if the fluid is sub-critical

        "1) Check for the hot stream"
        
        if stream == 'hot':
            if not self.Transcritical_h and not self.h_incomp_flag: # If the hot side is not transcritical
                for i in range(1,len(self.hvec_h)-1):
                    # Check if enthalpy is equal to the dewpoint enthalpy of hot stream and hot stream is colder than cold stream (impossible)
                    if (abs(self.hvec_h[i] - self.h_hdew) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]):
                        # Enthalpy of the cold stream at the pinch temperature
                        self.AS_C.update(CP.PT_INPUTS, self.pvec_c[i], self.T_hdew)
                        h_c_pinch = self.AS_C.hmass() # Equation 10 (Bell et al. 2015)

                        # Heat transfer in the cell
                        Qright = self.mdot_h*(self.h_hi-self.h_hdew) # Equation 9 (Bell et al. 2015)
                        
                        # New value for the limiting heat transfer rate
                        Qmax = self.mdot_c*(h_c_pinch-self.h_ci) + Qright # Equation 12

                        # Recalculate the cell boundaries
                        self.calculate_cell_boundaries(Qmax)
                
                        return Qmax
                    
            elif self.Transcritical_h:
                #If the hot side is transcritical, do nothing
                pass
            
            "2) Check for the cold stream"
            
        elif stream == 'cold':
            if not self.Transcritical_c:
                # Check for the cold stream pinch point
                for i in range(1,len(self.hvec_c)-1):
                    
                    # Check if enthalpy is equal to the bubblepoint enthalpy of cold stream and hot stream is colder than cold stream (impossible)
                    if (abs(self.hvec_c[i] - self.h_cbubble) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]):
                        # Enthalpy of the hot stream at the pinch temperature
                        self.AS_H.update(CP.PT_INPUTS, self.pvec_h[i], self.T_cbubble)
                        h_h_pinch = self.AS_H.hmass() # Equation 10 (Bell et al. 2015)

                        # Heat transfer in the cell
                        Qleft = self.mdot_c*(self.h_cbubble-self.h_ci) # Equation 13 (Bell et al. 2015)

                        # New value for the limiting heat transfer rate
                        Qmax = Qleft + self.mdot_h*(self.h_hi-h_h_pinch) # Equation 16 (Bell et al. 2015)

                        # Recalculate the cell boundaries
                        self.calculate_cell_boundaries(Qmax)

                        return Qmax
                    
            elif self.Transcritical_c:
                #If the hot side is transcritical, do nothing
                pass
        else:
            raise ValueError

    #%% HTC CORRELATION CHOOSING RELATED METHODS
    
    def compute_H_1P_HTC(self, k, Th_mean, p_h_mean, T_wall_h, G_h, havg_h):
        try:
            mu_h, Pr_h, k_h, mu_h_w, mu_rat, Pr_h_w, _ = propsfluid_AS(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, False, self.AS_H)       
        except (ValueError):
            if self.phases_h[k] == "liquid":
                mu_h, Pr_h, k_h, mu_h_w, mu_rat, Pr_h_w, _ = propsfluid_AS(Th_mean, p_h_mean, T_wall_h-1, self.H_su.fluid, False, self.AS_H)
            elif self.phases_h[k] == "vapor":
                mu_h, Pr_h, k_h, mu_h_w, mu_rat, Pr_h_w, _ = propsfluid_AS(Th_mean, p_h_mean, T_wall_h+1, self.H_su.fluid, False, self.AS_H)                
        if self.H.Correlation_1phase == "Gnielinski":
            if self.HTX_Type == 'Plate':
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_h, mu_h_w, Pr_h, k_h, G_h, self.geom.H_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_h, Pr_h, k_h, G_h, self.geom.H_Dh) # 
            elif self.HTX_Type == 'Shell&Tube':
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'Tube&Fins':
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'PCHE':
                Dh = np.pi*self.params['D_c']/(2+np.pi)
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, Dh, self.params['L_c']) 
        elif self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC":
            alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean, T_wall_h, p_h_mean, self.H_su.fluid, self.params)
        elif self.H.Correlation_1phase == 'Shell_Kern_HTC':
            alpha_h, self.Re_h[k], self.Pr_h[k] = shell_htc_kern(self.mdot_h, T_wall_h, Th_mean, p_h_mean, self.AS_H, self.params)      
        elif self.H.Correlation_1phase == 'Tube_And_Fins':
            alpha_h = htc_tube_and_fins(self.H_su.fluid, self.params, p_h_mean, havg_h, self.mdot_h, self.params['Fin_type'])[0]
        elif self.H.Correlation_1phase == 'water_plate_HTC':
            alpha_h = water_plate_HTC(mu_h, Pr_h, k_h, G_h, self.params['H_Dh'])
        elif self.H.Correlation_1phase == 'martin_holger_plate_HTC':
            alpha_h = martin_holger_plate_HTC(mu_h, Pr_h, k_h, self.mdot_h, self.params['H_n_canals'], Th_mean, p_h_mean, self.H_su.fluid, self.params['H_Dh'], self.params['l'], self.params['w'], self.params['amplitude'], self.params['chevron_angle'])
    
        return alpha_h

    def compute_C_1P_HTC(self, k, Tc_mean, p_c_mean, T_wall_c, G_c, havg_c):
        
        try:
            mu_c, Pr_c, k_c, mu_c_w, mu_rat, Pr_c_w, _ = propsfluid_AS(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False, self.AS_C)
        except (ValueError):
            if self.phases_c[k] == "liquid":
                mu_c, Pr_c, k_c, mu_c_w, mu_rat, Pr_c_w, _ = propsfluid_AS(Tc_mean-0.1, p_c_mean, T_wall_c, self.C_su.fluid, False, self.AS_C)
            elif self.phases_c[k] == "vapor":
                mu_c, Pr_c, k_c, mu_c_w, mu_rat, Pr_c_w, _ = propsfluid_AS(Tc_mean+0.1, p_c_mean, T_wall_c, self.C_su.fluid, False, self.AS_C)

        if self.HTX_Type == 'Plate' and (self.C_su.fluid == 'water' or self.C_su.fluid == 'Water'):
            alpha_c = water_plate_HTC(mu_c, Pr_c, k_c, G_c, self.params['C_Dh'])
        
        elif self.C.Correlation_1phase  == "Gnielinski":
            if self.HTX_Type == 'Plate':
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'Shell&Tube':
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']**self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'Tube&Fins':
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'PCHE':
                Dh = np.pi*self.params['D_c']/(2+np.pi)
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, Dh, self.params['L_c']) 
        elif self.C.Correlation_1phase == 'Shell_Bell_Delaware_HTC':
            alpha_c = shell_bell_delaware_htc(self.mdot_c, Tc_mean, T_wall_c, p_c_mean, self.C_su.fluid, self.params)
        elif self.C.Correlation_1phase == 'Shell_Kern_HTC':
            alpha_c, self.Re_c[k], self.Pr_c[k] = shell_htc_kern(self.mdot_c, T_wall_c, Tc_mean, p_c_mean, self.AS_C, self.params)
        elif self.C.Correlation_1phase == 'Tube_And_Fins':
            alpha_c = htc_tube_and_fins(self.C_su.fluid, self.params, p_c_mean, havg_c, self.mdot_c, self.params['Fin_type'])[0]
        elif self.C.Correlation_1phase == 'thonon_plate_HTC':
            alpha_c = thonon_plate_HTC(mu_c, Pr_c, k_c, G_c, self.params['C_Dh'], self.params['chevron_angle'])
        elif self.C.Correlation_1phase == 'kumar_plate_HTC':
            alpha_c = kumar_plate_HTC(mu_c, Pr_c, k_c, G_c, self.params['C_Dh'], mu_c_w, self.params['chevron_angle'])
        elif self.C.Correlation_1phase == 'martin_holger_plate_HTC':
            alpha_c = martin_holger_plate_HTC(mu_c, Pr_c, k_c, self.mdot_c, self.params['C_n_canals'], Tc_mean, p_c_mean, self.C_su.fluid, self.params['C_Dh'], self.params['l'], self.params['w'], self.params['amplitude'], self.params['chevron_angle'])

        return alpha_c

    # --------------------------------------------------------------    

    def compute_H_TC_HTC(self, k, Th_mean, p_h_mean, T_wall_h, G_h, havg_h):
        try:
            mu_h, Pr_h, k_h, mu_h_w, mu_rat, Pr_h_w, _ = propsfluid_AS(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, False, self.AS_H)       
        except (ValueError):
            if self.phases_h[k] == "liquid":
                mu_h, Pr_h, k_h, mu_h_w, mu_rat, Pr_h_w, _ = propsfluid_AS(Th_mean, p_h_mean, T_wall_h-1, self.H_su.fluid, False, self.AS_H)
            elif self.phases_h[k] == "vapor":
                mu_h, Pr_h, k_h, mu_h_w, mu_rat, Pr_h_w, _ = propsfluid_AS(Th_mean, p_h_mean, T_wall_h+1, self.H_su.fluid, False, self.AS_H)  
                
        if self.H.Correlation_TC == "Gnielinski":
            if self.HTX_Type == 'Plate':
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_h, mu_h_w, Pr_h, k_h, G_h, self.geom.H_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_h, Pr_h, k_h, G_h, self.geom.H_Dh) # 
            elif self.HTX_Type == 'Shell&Tube':
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'Tube&Fins':
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'PCHE':
                Dh = np.pi*self.params['D_c']/(2+np.pi)
                alpha_h, self.Re_h[k], self.Pr_h[k] = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, Dh, self.params['L_c']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 

        elif self.H.Correlation_TC == 'Liu_sCO2':
            self.AS_H.update(CP.PT_INPUTS, p_h_mean, Th_mean)
            rho_h = self.AS_H.rhomass()
            cp_h = self.AS_H.cpmass()
            
            alpha_h = Liu_sCO2(G_h, p_h_mean, T_wall_h, k_h, rho_h, mu_h, cp_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.H_su.fluid)
        
        elif self.H.Correlation_TC == 'Cheng_sCO2':
            q = self.Qvec_h[k]/(self.params['A_eff']*self.w[k])
            
            # print(f"q : {q}")
            # print(f"G : {G_h}")
            
            if k+1 > len(self.hvec_h):
                h_next = self.hvec_h[k+1]
            else:
                h_next = self.hvec_h[k]                
                
            alpha_h = Cheng_sCO2(G_h, q, T_wall_h, p_h_mean, self.hvec_h[k], h_next, mu_h, k_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.H_su.fluid)
            
            # print(f"-----------------")
            
        elif self.H.Correlation_TC == 'Lee':
            self.AS_H.update(CP.PT_INPUTS, p_h_mean, Th_mean)
            rho_h = self.AS_H.rhomass()
            alpha_h = PCHE_Lee(self.params['alpha'], self.params['D_c'], G_h, k_h, self.params['L_c'], mu_h, Pr_h, rho_h)

        elif self.H.Correlation_TC == "PCHE_Lee":
            alpha_h = PCHE_conv(self.params['alpha'], self.params['D_c'], G_h, k_h, self.params['L_c'], mu_h, mu_h_w, Pr_h, Th_mean, self.params['type_channel'])

        return alpha_h

    def compute_C_TC_HTC(self, k, Tc_mean, p_c_mean, T_wall_c, G_c, havg_c):
        
        try:
            mu_c, Pr_c, k_c, mu_c_w, mu_rat, Pr_c_w, _ = propsfluid_AS(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False, self.AS_C)
        except (ValueError):
            if self.phases_c[k] == "liquid":
                mu_c, Pr_c, k_c, mu_c_w, mu_rat, Pr_c_w, _ = propsfluid_AS(Tc_mean-0.1, p_c_mean, T_wall_c, self.C_su.fluid, False, self.AS_C)
            elif self.phases_c[k] == "vapor":
                mu_c, Pr_c, k_c, mu_c_w, mu_rat, Pr_c_w, _ = propsfluid_AS(Tc_mean+0.1, p_c_mean, T_wall_c, self.C_su.fluid, False, self.AS_C)
        
        if self.C.Correlation_TC  == "Gnielinski":
            if self.HTX_Type == 'Plate':
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'Shell&Tube':
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']**self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'Tube&Fins':
                alpha_c, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
            elif self.HTX_Type == 'PCHE':
                Dh = np.pi*self.params['D_c']/(2+np.pi)
                alpha_c, self.Re_c[k], self.Pr_c[k] = gnielinski_pipe_htc(mu_c, Pr_c, mu_c_w, k_c, G_c, Dh, self.params['L_c']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
 
        elif self.C.Correlation_TC == 'Liu_sCO2':
            self.AS_C.update(CP.PT_INPUTS, p_c_mean, Tc_mean)
            rho_c = self.AS_C.rhomass()
            cp_c = self.AS_C.cpmass()
            
            alpha_c = Liu_sCO2(G_c, p_c_mean, T_wall_c, k_c, rho_c, mu_c, cp_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.C_su.fluid)
        
        elif self.C.Correlation_TC == 'Cheng_sCO2':
            q = self.Qvec_c[k]/(self.params['A_eff']*self.w[k])
            
            if k+1 > len(self.hvec_c):
                h_next = self.hvec_c[k+1]
            else:
                h_next = self.hvec_c[k]                
                
            alpha_c = Cheng_sCO2(G_c, q, T_wall_c, p_c_mean, self.hvec_c[k], h_next, mu_c, k_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.C_su.fluid)
        
        elif self.C.Correlation_TC == 'Lee':
            self.AS_C.update(CP.PT_INPUTS, p_c_mean, Tc_mean)
            rho_c = self.AS_C.rhomass()
            alpha_c = PCHE_Lee(self.params['alpha'], self.params['D_c'], G_c, k_c, self.params['L_c'], mu_c, Pr_c, rho_c)
        
        elif self.H.Correlation_TC == "PCHE_Lee":
            alpha_c = PCHE_conv(self.params['alpha'], self.params['D_c'], G_c, k_c, self.params['L_c'], mu_c, mu_c_w, Pr_c, Tc_mean, self.params['type_channel'])
        
        return alpha_c

    # --------------------------------------------------------------    

    def compute_H_2P_HTC(self, k, Th_mean, p_h_mean, T_wall_h, G_h, havg_h, Th_sat_mean):
        
        if self.phases_h[k] == "two-phase":
            x_h = min(1, max(0, 0.5*(self.x_vec_h[k] + self.x_vec_h[k])))
        elif self.phases_h[k] == "vapor-wet":
            x_h = 1
                    
        # Thermodynamical variables
        self.AS_H.update(CP.PQ_INPUTS, p_h_mean, 0)
        
        try:
            mu_h_l = self.AS_H.viscosity()
        except:
            mu_h_l = CP.PropsSI('V', 'P', p_h_mean, 'Q', 0, self.su_H.fluid)
            
        rho_h_l = self.AS_H.rhomass()

        if self.H_su.fluid != 'R1233zd(E)':
            k_h_l = self.AS_H.conductivity()
            Pr_h_l = self.AS_H.Prandtl()
        else:
            k_h_l = conducticity_R1233zd(Th_sat_mean, p_h_mean)
            cp_h_l = self.AS_H.cpmass()
            Pr_h_l = mu_h_l * cp_h_l / k_h_l

        self.AS_H.update(CP.PQ_INPUTS, p_h_mean, 1)
        rho_h_v = self.AS_H.rhomass()
                    
        if self.H.Correlation_2phase == "Han_cond_BPHEX":
            alpha_h_2phase, _, DP_H = han_cond_BPHEX_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, rho_h_l, rho_h_v, G_h, self.params['H_Dh'], self.params['plate_pitch_co'], self.params['chevron_angle'], self.params['l_v'], self.params['H_n_canals'], self.H_su.m_dot, self.params['H_canal_t'])
        if self.H.Correlation_2phase == 'ext_tube_film_condens':
            self.AS_H.update(CP.HmassP_INPUTS, havg_h, p_h_mean)
            V_flow = G_h/self.AS_H.rhomass()
            if self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC":
                try: 
                    alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean, T_wall_h, p_h_mean, self.H_su.fluid, self.params)
                except:
                    alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean-0.1, T_wall_h, p_h_mean, self.H_su.fluid, self.params)
            
            if self.H.Correlation_1phase == "Shell_Kern_HTC":
                try: 
                    alpha_h, self.Re_h[k], self.Pr_h[k] = shell_htc_kern(self.mdot_h, T_wall_h, Th_mean, p_h_mean, self.AS_H, self.params) 
                except:
                    alpha_h, self.Re_h[k], self.Pr_h[k] = shell_htc_kern(self.mdot_h, T_wall_h, Th_mean-0.1, p_h_mean, self.AS_H, self.params)       
            
            elif self.H.Correlation_1phase == 'Tube_And_Fins':
                alpha_h = htc_tube_and_fins(self.H_su.fluid, self.params, p_h_mean, havg_h, self.mdot_h, self.params['Fin_type'])[0]
            alpha_h_2phase = ext_tube_film_condens(self.params['Tube_OD'], self.H_su.fluid, Th_mean, T_wall_h, V_flow)
        if self.H.Correlation_2phase == 'Horizontal_Tube_Internal_Condensation':
            self.AS_H.update(CP.HmassP_INPUTS, havg_h, p_h_mean)
            mu_h = self.AS_H.viscosity()
            
            if self.H_su.fluid != 'R1233zd(E)':
                k_h = self.AS_H.conductivity()
                Pr_h = self.AS_H.Prandtl()
            else:
                k_h = conducticity_R1233zd(Th_mean, p_h_mean)
                cp_h = self.AS_H.cpmass()
                Pr_h = mu_h * cp_h / k_h_l

            self.AS_H.update(CP.PT_INPUTS, p_h_mean, T_wall_h)
            Pr_h_w = self.AS_H.Prandtl()
            mu_h_w = self.AS_H.viscosity()   
                        
            alpha_h, self.Re_h[k], self.Pr_h[k]  = gnielinski_pipe_htc(mu_h, Pr_h, mu_h_w, k_h, G_h, self.params['Tube_OD'] - 2*self.params['Tube_t'], self.params['Tube_L']*self.params["Tube_pass"])
            alpha_h_2phase = horizontal_tube_internal_condensation(self.H_su.fluid , G_h, p_h_mean, self.x_vec_h[k], T_wall_h, self.params['Tube_OD'] - 2*self.params['Tube_t'])
        if self.H.Correlation_2phase == 'shah_condensation_plate_HTC':
            alpha_h_2phase = shah_condensation_plate_HTC(self.params['H_Dh'], self.params['l_v'], self.params['w_v'], self.params['amplitude'], self.params['phi'], self.mdot_h, p_h_mean, self.params['H_n_canals'], self.H_su.fluid)
        
        if self.H.Correlation_2phase == 'Tube_And_Fins':
            alpha_h_2phase = htc_tube_and_fins(self.H_su.fluid, self.params, p_h_mean, havg_h, self.mdot_h, self.params['Fin_type'])[0]
            
        if self.phases_h[k] == "two-phase" or self.phases_h[k] == "vapor-wet":
            alpha_h = alpha_h_2phase

        # elif self.phases_h[k] == "vapor-wet":
        #     w_vap_wet = (Th_mean - Th_sat_mean)/(Th_mean - T_wall_h)
        #     # The line below re-calculates alpha_h in case of having a vapor-wet condition
        #     alpha_h = alpha_h_2phase - w_vap_wet*(alpha_h_2phase - alpha_h) # the last alpha_h in this equation is the 1 Phase calculation

        return alpha_h

    def compute_C_2P_HTC(self, k, Tc_mean, p_c_mean, T_wall_c, G_c, Tc_sat_mean, alpha_h, LMTD, havg_c):

        try:
            a = self.x_c_calc 
        except:
            self.x_c_calc = []
            
        if self.phases_c[k] == "two-phase": 
            if len(self.phases_c[k]) == 1: 
                x_c = min(1, max(0, (self.x_vec_c[k]+1)/2))
            else:
                try: 
                    if self.phases_c[k+1] != "two-phase":
                        x_c = min(1, max(0, (self.x_vec_c[k]+1)/2))
                    else:
                        x_c = min(1, max(0, (self.x_vec_c[k]+self.x_vec_c[k+1])/2))
                except:
                    x_c = min(1, max(0, (self.x_vec_c[k]+1)/2))
                    
        elif self.phases_c[k] == "vapor-wet":
            x_c = 1
                
        # Thermodynamical variables
        self.AS_C.update(CP.PQ_INPUTS, p_c_mean, 0)
        
        try:
            mu_c_l = self.AS_C.viscosity()
        except:
            mu_c_l = CP.PropsSI('V', 'P', p_c_mean, 'Q', 0, self.su_C.fluid)
            
        rho_c_l = self.AS_C.rhomass()
        h_sat_c_l = self.AS_C.hmass()   
          
        if self.C_su.fluid == 'R1233zd(E)':
            k_c_l = conducticity_R1233zd(Tc_mean, p_c_mean)
            cp_c_l = self.AS_C.cpmass()
            Pr_c_l = mu_c_l * cp_c_l / k_c_l
        else:
            k_c_l = self.AS_C.conductivity()
            Pr_c_l = self.AS_C.Prandtl()
            
        self.AS_C.update(CP.PQ_INPUTS, p_c_mean, 1)
        rho_c_v = self.AS_C.rhomass()
        h_sat_c_v = self.AS_C.hmass()
        i_fg_c = h_sat_c_v - h_sat_c_l 
        
        # This boolean serves to spare to recalculate this properties in the dry-out analysis further ahead.
        self.ColdSide_Props_Calculated = True
        # !!! Include different types of Correlation HERE
        if self.C.Correlation_2phase == "Han_cond_BPHEX":
            alpha_c_2phase, _ = han_cond_BPHEX_HTC(x_c, mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v, G_c, self.params['C_Dh'], self.params['plate_pitch_co'], self.params['chevron_angle'])
        elif self.C.Correlation_2phase == "Han_Boiling_BPHEX_HTC":
            alpha_c_2phase, _ = han_boiling_BPHEX_HTC(min(x_c, self.x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, LMTD*self.F[k], self.Qvec_c[k], alpha_h, self.params['C_Dh'], self.params['chevron_angle'], self.params['plate_pitch_co'])
        elif self.C.Correlation_2phase == "Boiling_curve":              
            alpha_c_2phase = self.C_f_boiling(abs(T_wall_c - Tc_mean))
        elif self.C.Correlation_2phase == "gnielinski_pipe_htc":
            alpha_c_2phase, self.Re_c[k], self.Pr_c[k]  = gnielinski_pipe_htc(mu_c_l, Pr_c_l, mu_c_l, k_c_l, G_c, self.params['C_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) #
        elif self.C.Correlation_2phase == "amalfi_plate_HTC":
            alpha_c_2phase = amalfi_plate_HTC(self.params['C_Dh'], self.params['l'], self.params['w'], self.params['amplitude'], self.params['chevron_angle'], self.params['C_n_canals'], self.params['A_c'], self.mdot_c, p_c_mean, self.C_su.fluid)
        elif self.C.Correlation_2phase == "Flow_boiling":
            q = self.Qvec_c[k]/(self.params['A_eff']*self.w[k])
            D_in = self.params['Tube_OD']-2*self.params['Tube_t']
            
            alpha_c_2phase = horizontal_flow_boiling(self.AS_C, G_c, p_c_mean, x_c, D_in, q)
        
        elif self.C.Correlation_2phase == "Flow_boiling_gungor_winterton":
            q = self.Qvec_c[k]/(self.params['A_eff']*self.w[k])
            D_in = self.params['Tube_OD']-2*self.params['Tube_t']
            
            alpha_c_2phase = flow_boiling_gungor_winterton(self.su_C.fluid, G_c, p_c_mean, x_c, D_in, q, mu_c_l, Pr_c_l, k_c_l)
        
        else:
            raise ValueError("Correlation not found for Cold Side 2-Phase")
            
        if self.phases_c[k] == "two-phase":
            alpha_c = alpha_c_2phase
        elif self.phases_c[k] == "vapor-wet":
            w_vap_wet = (Tc_mean - Tc_sat_mean)/(Tc_mean - T_wall_c)
            #The line below re-calculates alpha_h in case of having a vapor-wet condition
            alpha_c = alpha_c_2phase - w_vap_wet*(alpha_c_2phase - alpha_c) #the last alpha_c in this equation is the 1 Phase calculation

        return alpha_c

    #%% PRESSURE DROP CORRELATION CHOOSING RELATED METHODS

    def compute_H_DP(self):
        
        if self.H.PressureDrop_Correlation == "Shell_Bell_Delaware_DP":
            DP_H = shell_bell_delaware_DP(self.su_H.m_dot, self.su_H.h, self.su_H.p, self.su_H.fluid, self.params)*self.params["n_series"]
        elif self.H.PressureDrop_Correlation == "Shell_Kern_DP":
            DP_H = shell_DP_kern(self.su_H.m_dot, (self.su_H.T + self.su_C.T)/2, self.su_H.h, self.su_H.p, self.su_H.fluid, self.params)*self.params["n_series"]
        elif self.H.PressureDrop_Correlation == "Gnielinski_DP":
            mu_h_in = CP.PropsSI('V', 'H', self.su_H.h, 'P', self.su_H.p, self.su_H.fluid)
            G_c, G_h = self.G_h_c_computation()

            if self.HTX_Type == 'PCHE':
                Dh = np.pi*self.params['D_c']/(2+np.pi)
                DP_H = gnielinski_pipe_DP(mu_h_in, self.su_H.D, G_h, Dh, self.params["L_c"], type_HX= 'PCHE')  
            else:
                DP_H = gnielinski_pipe_DP(mu_h_in, self.su_H.D, G_h, self.params["Tube_OD"]-2*self.params["Tube_t"], self.params["Tube_L"]*self.params["Tube_pass"])  

        elif self.H.PressureDrop_Correlation == 'Darcy_Weisbach':
            mu_h_in = CP.PropsSI('V', 'H', self.su_H.h, 'P', self.su_H.p, self.su_H.fluid)
            G_c, G_h = self.G_h_c_computation()
        
            Dh = np.pi*self.params['D_c']/(2+np.pi)
            DP_H = Darcy_Weisbach(mu_h_in, self.su_H.D, G_h, Dh, self.params["L_c"])  


        elif self.H.PressureDrop_Correlation == "Cheng_CO2_DP":
            G_c, G_h = self.G_h_c_computation()
            mu_h_in = CP.PropsSI('V', 'H', self.su_H.h, 'P', self.su_H.p, self.su_H.fluid)
            
            DP_H = Cheng_CO2_DP(G_h, self.params["Tube_OD"]-2*self.params["Tube_t"], self.params["Tube_L"]*self.params["Tube_pass"], self.su_H.p, self.su_H.h, mu_h_in, self.su_H.fluid)
        else:
            DP_H = 0
            
        return DP_H

    def compute_cell_H_DP(self):
        if self.H.PressureDrop_Correlation == "Shell_Bell_Delaware_DP":
            DP_H = shell_bell_delaware_DP(self.su_H.m_dot, self.su_H.h, self.su_H.p, self.su_H.fluid, self.params)*self.params["n_series"]/self.params['n_disc']
        elif self.H.PressureDrop_Correlation == "Shell_Kern_DP":
            DP_H = shell_DP_kern(self.su_H.m_dot, (self.su_H.T + self.su_C.T)/2, self.su_H.h, self.su_H.p, self.su_H.fluid, self.params)*self.params["n_series"]/self.params['n_disc']
        else:
            DP_H = 0
            
        return DP_H

    def compute_C_DP(self):

        if self.C.PressureDrop_Correlation == "Shell_Bell_Delaware_DP":
            DP_C = shell_bell_delaware_DP(self.su_C.m_dot, self.su_C.h, self.su_C.p, self.su_C.fluid, self.params)*self.params["n_series"]
        elif self.C.PressureDrop_Correlation == "Shell_Kern_DP":
            DP_C = shell_DP_kern(self.su_C.m_dot, (self.su_H.T + self.su_C.T)/2, self.su_C.h, self.su_C.p, self.su_C.fluid, self.params)*self.params["n_series"]
        elif self.C.PressureDrop_Correlation == "Gnielinski_DP":
            
            mu_c_in = CP.PropsSI('V', 'H', self.su_C.h, 'P', self.su_C.p, self.su_C.fluid)
            G_c, G_h = self.G_h_c_computation()
            
            # print(f"mu_c_in : {mu_c_in}")
            # print(f"rho_c_in : {self.su_C.D}")
            # print(f"G_c : {G_c}")
            # print(f"D_in : {self.params['Tube_OD']-2*self.params['Tube_t']}")
            # print(f"L : {self.params['Tube_L']*self.params['Tube_pass']}")
            # print(f"Fluid : {self.su_C.fluid}")
            
            if self.HTX_Type == 'PCHE':
                Dh = np.pi*self.params['D_c']/(2+np.pi)
                DP_C = gnielinski_pipe_DP(mu_c_in, self.su_C.D, G_c, Dh, self.params["L_c"], type_HX= 'PCHE')   
            else:
                DP_C = gnielinski_pipe_DP(mu_c_in, self.su_C.D, G_c, self.params["Tube_OD"]-2*self.params["Tube_t"], self.params["Tube_L"]*self.params["Tube_pass"])   
        
        elif self.H.PressureDrop_Correlation == 'Darcy_Weisbach':
            mu_c_in = CP.PropsSI('V', 'H', self.su_C.h, 'P', self.su_C.p, self.su_C.fluid)
            G_c, G_h = self.G_h_c_computation()
        
            Dh = np.pi*self.params['D_c']/(2+np.pi)
            DP_C = Darcy_Weisbach(mu_c_in, self.su_C.D, G_c, Dh, self.params["L_c"])  

        
        elif self.C.PressureDrop_Correlation == "Choi_DP":

            G_c, G_h = self.G_h_c_computation()
            rho_out = CP.PropsSI('D', 'P', self.su_C.p, 'Q', 1, self.su_C.fluid)
            
            DP_C = Choi_DP(self.su_C.fluid, G_c, rho_out, self.su_C.D, self.su_C.p, 1, self.su_C.x, self.params["Tube_L"]*self.params["Tube_pass"], self.params["Tube_OD"]-2*self.params["Tube_t"])
        
        elif self.C.PressureDrop_Correlation == "Muller_Steinhagen_Heck_DP":
            G_c, G_h = self.G_h_c_computation()
            rho_out = CP.PropsSI('D', 'P', self.su_C.p, 'Q', 1, self.su_C.fluid)

            DP_C = Muller_Steinhagen_Heck_DP(self.su_C.fluid, G_c, self.su_C.p, 0.37, self.params["Tube_L"]*self.params["Tube_pass"], self.params["Tube_OD"]-2*self.params["Tube_t"])
        
        elif self.C.PressureDrop_Correlation == "Cheng_CO2_DP":
            G_c, G_h = self.G_h_c_computation()
            mu_c_in = CP.PropsSI('V', 'H', self.su_C.h, 'P', self.su_C.p, self.su_C.fluid)
            
            DP_C = Cheng_CO2_DP(G_c, self.params["Tube_OD"]-2*self.params["Tube_t"], self.params["Tube_L"]*self.params["Tube_pass"], self.su_C.p, self.su_C.h, mu_c_in, self.su_C.fluid)
        
        else:
            DP_C = 0
        
        return DP_C

    def compute_cell_C_DP(self):

        if self.C.PressureDrop_Correlation == "Gnielinski_DP":
            mu_c_in = CP.PropsSI('V', 'H', self.su_C.h, 'P', self.su_C.p, self.su_C.fluid)
            G_c, G_h = self.G_h_c_computation()
            
            DP_C = gnielinski_pipe_DP(mu_c_in, self.su_C.D, G_c, self.params["Tube_OD"]-2*self.params["Tube_t"], self.params["Tube_L"]*self.params["Tube_pass"]/self.params['n_disc'])    
        else:
            DP_C = 0
        
        return DP_C

    #%% SOLVE RELATED METHODS

    def G_h_c_computation(self):
        if self.HTX_Type == 'Plate':          
            G_c = (self.mdot_c/self.params['C_n_canals'])/self.params['C_CS']
            G_h = (self.mdot_h/self.params['H_n_canals'])/self.params['H_CS']
        elif self.HTX_Type == 'Shell&Tube':
            A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            G_h = (self.params["Tube_pass"]/self.params["n_series"])*self.mdot_h/(A_in_one_tube*self.params['n_tubes'])
            G_c = (self.params["Tube_pass"]/self.params["n_series"])*self.mdot_c/(A_in_one_tube*self.params['n_tubes'])
        elif self.HTX_Type == 'Tube&Fins':
            A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            G_h = (self.params["Tube_pass"]/self.params["n_series"])*(self.mdot_h/self.params['n_tubes'])/A_in_one_tube
            G_c = (self.params["Tube_pass"]/self.params["n_series"])*(self.mdot_c/self.params['n_tubes'])/A_in_one_tube
        elif self.HTX_Type == 'PCHE': 
            A_in_one_channel = np.pi*(self.params['D_c']**2)/8
            G_h = self.mdot_h/(A_in_one_channel*self.params['N_c']*self.params['N_p']*(self.params['R_p']/(1+self.params['R_p'])))
            G_c = self.mdot_c/(A_in_one_channel*self.params['N_c']*self.params['N_p']*(1/(1+self.params['R_p'])))
            
        return G_c, G_h

    def solve(self, only_external = False, and_solve = True):
        # OK
        """
        Parameters
        ----------
            - mdot_h : Hot fluid flowrate [kg/s]
            - p_hi : Hot fluid pressure [Pa]
            - h_hi : Hot fluid specific enthalpy [J/kg]
            - mdot_c : Cold fluid flowrate [kg/s]
            - p_ci : Cold fluid pressure [Pa]
            - h_ci : Cold fluid specific enthalpy [J/kg]
            - only_external (optional) : calls only exxternal_pinching function to determine Q_dot max if set to True (if there is no phase change)
            - and_solve (optional) : if set to true, solves the heat exchanger

        Returns
        -------
            - Q : Exchanged heat rate [W]
        """
        
        self.check_calculable()
        self.check_parametrized()
                
        if not self.calculable:
            print("Component not calculable, check input")
            
        if not self.parametrized:
            print("Component not parametrized, check parameters")            
            
        "1) Main Input variables"
        
        self.H_su = self.su_H
        self.C_su = self.su_C
            
        self.H_ex = self.ex_H
        self.C_ex = self.ex_C
        
        # Incompressible fluid tag setting
        self.h_incomp_flag = (self.H_su.fluid.find("INCOMP") != -1)
        self.c_incomp_flag = (self.C_su.fluid.find("INCOMP") != -1)

        # Hot fluid
        self.mdot_h = self.H_su.m_dot
        self.h_hi = self.H_su.h
        self.p_hi = self.H_su.p
        
        # Cold fluid 
        self.mdot_c = self.C_su.m_dot
        self.h_ci = self.C_su.h
        self.p_ci = self.C_su.p
        
        # Instantiate the abstractstates
        if self.h_incomp_flag:
            rest, fluid_name = self.H_su.fluid.split("::")
            self.AS_H = CP.AbstractState("INCOMP", fluid_name)
        else:
            self.AS_H = CP.AbstractState("BICUBIC&HEOS", self.H_su.fluid)  
            
        if self.c_incomp_flag:
            rest, fluid_name = self.C_su.fluid.split("::")
            self.AS_C = CP.AbstractState("INCOMP", fluid_name)
        else:
            self.AS_C = CP.AbstractState("BICUBIC&HEOS", self.C_su.fluid)    

        self.AS_C.update(CP.HmassP_INPUTS, self.h_ci, self.p_ci)
        self.AS_H.update(CP.HmassP_INPUTS, self.h_hi, self.p_hi)

        self.T_ci = self.AS_C.T()
        self.T_hi = self.AS_H.T()
        
        # Determine the inlet temperatures from the pressure/enthalpy pairs
        x_ci = self.AS_C.Q()
        if x_ci > 0.999 or x_ci < 0.001: # Verify that the quality is well between 0 and 1
            self.T_ci = self.AS_C.T()
        else:
            try:
                self.AS_C.update(CP.HmassQ_INPUTS, self.h_ci, 0.5)
                self.T_ci = self.AS_C.T()
            except:
                self.AS_C.update(CP.PQ_INPUTS, self.p_ci, 0.5)
                self.T_ci = self.AS_C.T()
                
        if not self.h_incomp_flag:
        
            x_hi = self.AS_H.Q()

            if x_hi > 0.999 or x_hi < 0.001:
                self.AS_H.T()
            else:
                try:    
                    self.AS_H.update(CP.HQ_INPUTS, 0.5, self.h_hi)
                    self.T_hi = self.AS_H.T()
                except:
                    self.AS_H.update(CP.PQ_INPUTS, self.p_hi, 0.5)
                    self.T_hi = self.AS_H.T() # If the fluid is R1233zd(E) the T = f(H) is not implemented in CoolProp yet
        
        else: 
            self.AS_H.T()
                
        "2) Determine if the streams come in at a higher pressure than the transcritical pressure"
        
        # Initialization of the flags:    
        self.Transcritical_c = False
        self.Transcritical_h = False
        
        if not self.h_incomp_flag and (self.p_hi - self.AS_H.p_critical()) >= 1e-06:
            self.Transcritical_h = True

        if not self.c_incomp_flag and (self.p_ci - self.AS_C.p_critical()) >= 1e-06:
            self.Transcritical_c = True

        "3) Calculate the ideal bubble and dew temperatures/enthalpies for each stream IF the fluid is not transcritical"
        
        if not self.Transcritical_c:
            self.AS_C.update(CP.PQ_INPUTS, self.p_ci, 0)
            self.T_cbubble_ideal = self.AS_C.T()
            self.h_cbubble_ideal = self.AS_C.hmass()

            self.AS_C.update(CP.PQ_INPUTS, self.p_ci, 1)
            self.T_cdew_ideal    = self.AS_C.T()
            self.h_cdew_ideal    = self.AS_C.hmass()
            
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.AS_H.update(CP.PQ_INPUTS, self.p_hi, 0)
            self.T_hbubble_ideal = self.AS_H.T()
            self.h_hbubble_ideal = self.AS_H.hmass()

            self.AS_H.update(CP.PQ_INPUTS, self.p_hi, 1)
            self.T_hdew_ideal    = self.AS_H.T()
            self.h_hdew_ideal    = self.AS_H.hmass()
            
        "4) Calculate pressure drops"

        if self.params['H_DP_ON'] == True: # if the pressure drop are not neglected
            if self.H.PressureDrop_Correlation == "User-Defined":
                self.p_ho = self.p_hi - self.H.DP_val
            else:
                self.DP_h = self.compute_H_DP()
                self.p_ho = self.p_hi - self.DP_h   
        else:
            self.DP_h = 0
            self.p_ho = self.p_hi - self.DP_h  
        
            # if self.H_su.fluid == 'Water':
            #     self.DP_h = self.H.f_dp["K"] * (self.mdot_h/self.H_su.D)**2
            #     self.p_ho = self.p_hi - self.DP_h
            # else:
            #     self.DP_h = 0 # self.H.f_dp["K"] * self.mdot_h**(self.H.f_dp["B"]) # Empirical correlations : DP = K*m_dot**B # Han_BPHEX_DP(mu_h, G_h, self.geom.H_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co, rho_v, rho_l, self.geom.l_v, self.geom.H_n_canals, self.mdot_h, self.geom.H_canal_t) # 
            #     self.p_ho = self.p_hi - self.DP_h
            
        if self.params['C_DP_ON'] == True: # if the pressure drop are not neglected
            if self.C.PressureDrop_Correlation == "User-Defined":
                self.p_co = self.p_ci - self.C.DP_val
            else:
                self.DP_c = self.compute_C_DP()
                self.p_co = self.p_ci - self.DP_c  
            
        else:
            self.DP_c = 0
            self.p_co = self.p_ci - self.DP_c  
            
        # if self.C.f_dp != {"K": 0, "B": 0}:
        #     if self.C_su.fluid == 'Water':
        #         self.DP_c = self.C.f_dp["K"] * (self.mdot_c/self.C_su.D)**2
        #         self.p_co = self.p_ci - self.DP_c
        #     else:
        #         self.DP_c = 0 # self.C.f_dp["K"] * self.mdot_c**(self.C.f_dp["B"]) # Han_BPHEX_DP(mu_c, G_c, self.geom.H_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co, rho_v, rho_l, self.geom.l_v, self.geom.C_n_canals, self.mdot_c, self.geom.C_canal_t) # 
        #         self.p_co = self.p_ci - self.DP_c
        
        
        "5) Calculate maximum and actual heat rates"
                
        if (self.T_hi - self.T_ci) > 1e-2  and self.mdot_h  > 0 and self.mdot_c > 0: # Check that the operating conditions allow for heat transfer
            "5.1) Compute the external pinching & update cell boundaries"
            Qmax_ext = self.external_pinching() # Call to external-pinching procedure
            self.Qmax_ext = Qmax_ext
            Qmax = Qmax_ext
            
            if debug:
                print("External pinching calculation done. \n")
            
            "5.2) Compute the internal pinching & update cell boundaries"
            if not only_external: # If phase change is expected : Check the internal pinching
                for stream in ['hot','cold']:
                    Qmax_int = self.internal_pinching(stream) # Call to internal-pinching procedure
                    if Qmax_int is not None:
                        self.Qmax_int = Qmax_int
                        Qmax = Qmax_int
            
            # Maximum heat transfer rate determined by external or internal pinching
            self.Qmax = Qmax
            
            "5.3) Solve the heat exchanger to find the actual heat rate"
            if and_solve and not only_external:
                Q = self.solve_hx()
            self.epsilon_th = self.Q/self.Qmax # HTX efficiency
            self.residual = 1 - sum(self.w) # HTX residual # !!! (what is "w" ?)
            
            "5.4) Effective density computation for each cell taking into account void fraction"
            
            self.Dvec_h = np.empty(len(self.hvec_h))
            self.Dvec_c = np.empty(len(self.hvec_c))
            
            # Cross section void fraction : considered constant along each discretization in the void_fraction function
            self.eps_void_h = np.empty(len(self.hvec_c)) 
            self.eps_void_c = np.empty(len(self.hvec_c))
            
            for i in range(len(self.hvec_h)):
                if self.x_vec_h[i] <= 0 or self.x_vec_h[i] >= 1: 
                    self.AS_H.update(CP.HmassP_INPUTS, self.hvec_h[i], self.pvec_h[i])
                    self.Dvec_h[i] = self.AS_H.rhomass()
                    self.eps_void_h[i] = -1
                else:
                    self.AS_H.update(CP.PQ_INPUTS, self.pvec_h[i], 1)
                    rho_g = self.AS_H.rhomass()

                    self.AS_H.update(CP.PQ_INPUTS, self.pvec_h[i], 0)
                    rho_l = self.AS_H.rhomass()

                    self.eps_void_h[i], self.Dvec_h[i] = void_fraction(self.x_vec_h[i], rho_g, rho_l)
                    
            for i in range(len(self.hvec_c)-1):
                if self.x_vec_c[i] <= 0 or self.x_vec_c[i] >= 1: 
                    try:    
                        self.AS_C.update(CP.HmassP_INPUTS, self.hvec_c[i], self.pvec_c[i])
                    except:
                        self.AS_C.update(CP.HmassP_INPUTS, (self.hvec_c[i+1] + self.hvec_c[i-1])/2, self.pvec_c[i])
                        
                    self.Dvec_c[i] = self.AS_C.rhomass()
                    self.eps_void_c[i] = -1
                else:
                    self.AS_C.update(CP.PQ_INPUTS, self.pvec_c[i], 1)
                    rho_g = self.AS_C.rhomass()

                    self.AS_C.update(CP.PQ_INPUTS, self.pvec_c[i], 0)
                    rho_l = self.AS_C.rhomass()

                    self.eps_void_c[i], self.Dvec_c[i] = void_fraction(self.x_vec_c[i], rho_g, rho_l)
            
            "5.5) Computation of the fluid mass inside the HTX"
            
            if self.HTX_Type == 'Plate':
                # OK
                self.Vvec_h = self.params['H_V_tot']*np.array(self.w) # !!! Attention to this assumption
                self.Vvec_c = self.params['C_V_tot']*np.array(self.w)

            elif self.HTX_Type == 'Shell&Tube':
                if self.params['Shell_Side'] == 'H': # Shell Side is the hot side
                    self.Vvec_h = self.params['S_V_tot']*np.array(self.w) # !!! Attention to this assumption
                    self.Vvec_c = self.params['T_V_tot']*np.array(self.w)
                else:
                    self.Vvec_h = self.params['T_V_tot']*np.array(self.w) # !!! Attention to this assumption
                    self.Vvec_c = self.params['S_V_tot']*np.array(self.w)                    

            elif self.HTX_Type == 'Tube&Fins':
                if self.params['Fin_Side'] == 'H': # Shell Side is the hot side
                    self.Vvec_h = self.params['B_V_tot']*np.array(self.w) # !!! Attention to this assumption
                    self.Vvec_c = self.params['T_V_tot']*np.array(self.w)
                else:
                    self.Vvec_h = self.params['T_V_tot']*np.array(self.w) # !!! Attention to this assumption
                    self.Vvec_c = self.params['B_V_tot']*np.array(self.w) 
            
            elif self.HTX_Type == 'PCHE':
                self.Vvec_h = self.params['H_V_tot']*np.array(self.w) # !!! Attention to this assumption
                self.Vvec_c = self.params['C_V_tot']*np.array(self.w)
                
            # Initiates the mass vectors # !!! (for each cell ?)
            self.Mvec_h = np.empty(len(self.hvec_h)-1)
            self.Mvec_c = np.empty(len(self.hvec_c)-1)
            
            # Mass computation # Volume*Mean of the densities at the cell boundaries
            for i in range(len(self.hvec_h)-1):
                self.Mvec_h[i] = self.Vvec_h[i]*(self.Dvec_h[i] + self.Dvec_h[i+1])/2
                self.Mvec_c[i] = self.Vvec_c[i]*(self.Dvec_c[i] + self.Dvec_c[i+1])/2

            if and_solve:
                
                # Hot side
                self.ex_H.set_fluid(self.H_su.fluid)
                self.ex_H.set_p(self.pvec_h[0])
                self.ex_H.set_h(self.hvec_h[0])
                self.ex_H.set_T(self.Tvec_h[0])
                self.ex_H.set_m_dot(self.H_su.m_dot)
                self.H_ex = self.ex_H
                    
                # Cold side
                self.ex_C.set_fluid(self.C_su.fluid)
                self.ex_C.set_p(self.pvec_c[-1]) # Example temperature [K]
                self.ex_C.set_h(self.hvec_c[-1]) # Example Pressure [Pa]
                self.ex_C.set_T(self.Tvec_c[-1])
                self.ex_C.set_m_dot(self.C_su.m_dot) 
                self.C_ex = self.ex_C
                                            
                self.solved = True
                return Q
            
        else: # Just a flag if the heat exchanger is not solved
            self.Q = 1
            raise Exception("Hot and cold temperatures seem to be reversed or a flow rate is negative.")
 
    def objective_function(self, Q):
        "1) Initialize cell boundaries and results vectors"
        # OK
        # print('Objective Function')
        self.calculate_cell_boundaries(Q)
        # print("Cell boundaries calculated")
        # Initialize dry-out incipience identification
        self.x_di_c = 1
        self.dry_out_c = False
        self.ColdSide_Props_Calculated = False
        
        # Initialise results arrays
        # w = [] 
        
        if self.w[0] == None or abs(sum(self.w)-1) > 0.1:
            self.w = np.ones(len(self.hvec_h)-1)/(len(self.hvec_h)-1)
            
        self.Re_h = np.zeros(len(self.hvec_h)-1)
        self.Re_c = np.zeros(len(self.hvec_h)-1)
        self.Pr_h = np.zeros(len(self.hvec_h)-1)
        self.Pr_c = np.zeros(len(self.hvec_h)-1)

        self.Avec_h = []
        self.alpha_c = []
        self.alpha_h = []
        self.phases_h = []
        self.phases_c = []
        
        "2) Iteration over all cells"
        
        # The length of the hvec vectors is the same. hvec_h is arbitrarily taken below
        for k in range(len(self.hvec_h)-1): # iterate over each cell
            "2.1) Diverse temperature parameter definition"
            
            Thi = self.Tvec_h[k+1] 
            Tci = self.Tvec_c[k]
            Tho = self.Tvec_h[k]
            Tco = self.Tvec_c[k+1]
            
            DTA = max(Thi - Tco, 1e-02)
            DTB = max(Tho - Tci, 1e-02)
            
            if DTA == DTB:
                self.LMTD[k] = DTA
            else:
                try:
                    self.LMTD[k] = (DTA-DTB)/np.log(abs(DTA/DTB))
                except ValueError as VE:
                    print(Q, DTA, DTB)
                    raise
            
            # Wall Temperature:
            T_wall = (Thi + Tho + Tci + Tco)/4

            Tc_mean = 0.5*(Tci + Tco) # mean temperature over the cell
            Th_mean = 0.5*(Thi + Tho) # mean temperature over the cell
            
            if not self.Transcritical_c:
                Tc_sat_mean = 0.5*(self.Tvec_sat_pure_c[k] + self.Tvec_sat_pure_c[k+1]) # mean saturation temperature over the cell
            if not self.Transcritical_h and not self.h_incomp_flag:
                Th_sat_mean = 0.5*(self.Tvec_sat_pure_h[k] + self.Tvec_sat_pure_h[k+1]) # mean saturation temperature over the cell
                        
            "2.2) Hot side phase identification"            
            # If not transcritical
            if not self.Transcritical_h and not self.h_incomp_flag:
                havg_h = (self.hvec_h[k] + self.hvec_h[k+1])/2.0 # Average cell enthalpy over the cell
                
                if havg_h < self.h_hbubble: # below evaporation/bubble point
                    self.phases_h.append('liquid')
                    T_wall_h = T_wall
                    
                elif havg_h > self.h_hdew: # higher than condensation/dew point
                    T_wall_h = max(T_wall, Th_sat_mean)                    
                    if T_wall > Th_sat_mean:
                        self.phases_h.append('vapor')
                    else:
                        self.phases_h.append('vapor-wet')
                        
                else: # Between evaporation/bubble and condensation/dew point
                    self.phases_h.append('two-phase')
                    T_wall_h = T_wall
                    
            elif self.h_incomp_flag:
                self.phases_h.append('liquid')
                T_wall_h = T_wall
                havg_h = (self.hvec_h[k] + self.hvec_h[k+1])/2.0 # Average cell enthalpy over the cell
                
            elif self.Transcritical_h:
                self.phases_h.append("transcritical")
                T_wall_h = T_wall
                havg_h = (self.hvec_h[k] + self.hvec_h[k+1])/2.0 # Average cell enthalpy over the cell
                
            "2.3) Cold side phase identification"
            # If not transcritical
            if not self.Transcritical_c:
                havg_c = (self.hvec_c[k] + self.hvec_c[k+1])/2.0 # Average cell enthalpy over the cell
                
                if havg_c < self.h_cbubble: # below evaporation/bubble point
                    self.phases_c.append('liquid')
                    T_wall_c = T_wall
                    
                elif havg_c > self.h_cdew: # higher than condensation/dew point
                    self.phases_c.append("vapor")
                    T_wall_c = T_wall
                    
                else: # Between evaporation/bubble and condensation/dew point
                    T_wall_c = T_wall
                    if (0.5*self.x_vec_c[k] + 0.5*self.x_vec_c[k+1]) <= self.x_di_c:
                        self.phases_c.append('two-phase')  
                    else:
                        self.phases_c.append("two-phase-dryout")
                        self.dry_out_c = True
            
            elif self.c_incomp_flag:
                self.phases_c.append('liquid')
                T_wall_c = T_wall
                havg_c = (self.hvec_c[k] + self.hvec_c[k+1])/2.0 # Average cell enthalpy over the cell
                
            elif self.Transcritical_c:
                self.phases_c.append('transcritical')
                T_wall_c = T_wall
                havg_c = (self.hvec_c[k] + self.hvec_c[k+1])/2.0 # Average cell enthalpy over the cell
                            
            p_c_mean = 0.5*(self.pvec_c[k] + self.pvec_c[k+1]) # mean pressure over the cell
            p_h_mean = 0.5*(self.pvec_h[k] + self.pvec_h[k+1]) # mean pressure over the cell

            G_c, G_h = self.G_h_c_computation()
            
            # F correction factor for LMTD method:
            if self.params['Flow_Type'] != "CounterFlow":
                try:
                    self.AS_C.update(CP.PT_INPUTS, p_c_mean, Tc_mean)
                    C_c = self.AS_C.cpmass()
                except:
                    C_c = 20000
                try:
                    self.AS_H.update(CP.PT_INPUTS, p_h_mean, Th_mean)
                    C_h = self.AS_H.cpmass()
                except:
                    C_h = 20000
                
                C_min = min(C_c,C_h)
                C_max = max(C_c,C_h)
                C_r = C_min/C_max
                
                if self.params['n_series'] > 1:
                    self.R = (Thi - Tho)/(Tco - Tci)
                    self.P = (Tco - Tci)/(Thi - Tci)
                else:
                    
                    self.R = (self.Tvec_h[-1] - self.Tvec_h[0])/(self.Tvec_c[-1] - self.Tvec_c[0])
                    self.P = (self.Tvec_c[-1] - self.Tvec_c[0])/(self.Tvec_h[-1] - self.Tvec_c[0])                    
                                    
                if self.params['Flow_Type'] == 'Shell&Tube':
                    
                    if self.P < 0.1 or self.R > 10 or self.R < 0.4:
                        self.F[k] = 1
                    elif self.P > 0.99:
                        self.F[k] = 0
                    else:
                        if np.mod(self.params['Tube_pass'],2) == 0:
                            self.F[k] = F_shell_and_tube(self.R,self.P,self.params['n_series'])
                        else:
                            self.F[k] = f_lmtd2(self.R, self.P, self.params, C_r)
                else:    
                    self.F[k] = f_lmtd2(self.R, self.P, self.params, C_r)
                    
            else:
                self.F[k] = 1
            
            self.F[k] = max(self.F[k],0)    
                        
            #UA_req including F correction factor in case of cross flow
            UA_req = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])/(self.F[k]*self.LMTD[k])          

            "3) Cell heat transfer coefficients"            
            "3.1) Hot side - User defined"
            if self.H.HeatExchange_Correlation == "User-Defined": # User Defined
                if self.phases_h[k] == "liquid":
                    alpha_h = self.H.h_liq
                elif self.phases_h[k] == "vapor":
                    alpha_h = self.H.h_vap
                elif self.phases_h[k] == "two-phase":
                    alpha_h = self.H.h_twophase
                elif self.phases_h[k] == "vapor-wet":
                    alpha_h = self.H.h_vapwet
                elif self.phases_h[k] == "two-phase-dryout":
                    alpha_h = self.H.h_tpdryout
                elif self.phases_h[k] == "transcritical":
                    alpha_h = self.H.h_transcrit   
                    
            elif self.H.HeatExchange_Correlation == "Correlation": # Heat transfer coefficient calculated from Correlations.
                # 1 PHASE CORRELATION
                if self.phases_h[k] == "liquid" or self.phases_h[k] == "vapor":
                    alpha_h = self.compute_H_1P_HTC(k, Th_mean, p_h_mean, T_wall_h, G_h, havg_h)
                elif self.phases_h[k] == "transcritical":
                    alpha_h = self.compute_H_TC_HTC(k, Th_mean, p_h_mean, T_wall_h, G_h, havg_h)
                # 2 PHASE CORRELATION
                elif self.phases_h[k] == "two-phase" or self.phases_h[k] == "vapor-wet":
                    alpha_h = self.compute_H_2P_HTC(k, Th_mean, p_h_mean, T_wall_h, G_h, havg_h, Th_sat_mean)
                    
            "3.2) Cold side - User defined"
            if self.C.HeatExchange_Correlation == "User-Defined": # User Defined
                if self.phases_c[k] == "liquid":
                    alpha_c = self.C.h_liq
                elif self.phases_c[k] == "vapor":
                    alpha_c = self.C.h_vap
                elif self.phases_c[k] == "two-phase":
                    alpha_c = self.C.h_twophase
                elif self.phases_c[k] == "vapor-wet":
                    alpha_c = self.C.h_vapwet
                elif self.phases_c[k] == "two-phase-dryout":
                    alpha_c = self.C.h_tpdryout
                elif self.phases_c[k] == "transcritical":
                    alpha_c = self.C.h_transcrit
                    
            elif self.C.HeatExchange_Correlation == "Correlation": # Heat transfer coefficient calculated from Correlations:
                # 1 PHASE CORRELATION
                if self.phases_c[k] == "liquid" or self.phases_c[k] == "vapor":
                    alpha_c = self.compute_C_1P_HTC(k, Tc_mean, p_c_mean, T_wall_c, G_c, havg_c)
                # 2 PHASE CORRELATION
                elif self.phases_c[k] == "transcritical":
                    alpha_c = self.compute_C_TC_HTC(k, Tc_mean, p_c_mean, T_wall_c, G_c, havg_c)
                elif self.phases_c[k] == "two-phase" or self.phases_c[k] == "vapor-wet":
                    alpha_c = self.compute_C_2P_HTC(k, Tc_mean, p_c_mean, T_wall_c, G_c, Tc_sat_mean, alpha_h, self.LMTD[k], havg_c)
                          
            "3.3) Store the heat transfer coefficients"
            self.alpha_c.append(alpha_c)
            self.alpha_h.append(alpha_h)

            "4) Recover some parameters and conductivity calculation and compute UA_avail"
            
            if self.HTX_Type == 'Plate':     
                # In the equation below, thickness resistance is given with respect to A_h arbitrarely
                UA_avail = 1/((1+self.params['fooling'])/(alpha_h*self.params['A_h']) + 1/(alpha_c*self.params['A_c'])) #+ self.params['t_plates']/(self.params['plate_cond'])) 
            
            elif self.HTX_Type == 'Shell&Tube':       
                fact_cond_1 = np.log(self.params['Tube_OD']/(self.params['Tube_OD'] - 2*self.params['Tube_t']))
                fact_cond_2 = 2*np.pi*self.params['tube_cond']*self.params['Tube_L']*self.params['n_series']*self.params['n_tubes']
                R_cond = fact_cond_1/fact_cond_2                      

                self.A_in_tubes = self.params['n_series']*self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t']))
                self.A_out_tubes = self.params['n_series']*self.params['Tube_L']*self.params['n_tubes']*np.pi*(self.params['Tube_OD'])
            
                if self.params['foul_s'] != None:
                    R_fouling_s = self.params['foul_s'] / self.A_out_tubes
                else: 
                    R_fouling_s = 0
    
                if self.params['foul_t'] != None:
                    R_fouling_t = self.params['foul_t'] / self.A_in_tubes
                else: 
                    R_fouling_t = 0                
                    
                try: 
                    self.params["Overdesign"]
                except:
                    self.params["Overdesign"] = 0
                    
                    
                UA_avail = (1-self.params["Overdesign"])/(1/(alpha_c*self.A_in_tubes) + 1/(alpha_h*self.A_out_tubes) + R_fouling_s + R_fouling_t) # + R_fouling_s + R_fouling_t + R_cond) 
                
                # print("---------------")
                # print(UA_avail)
                # print("---------------")
                # print(1/(alpha_c*self.A_in_tubes))
                # print(1/(alpha_h*self.A_out_tubes))
                # print(R_fouling_s)
                # print(R_fouling_t)
                
            elif self.HTX_Type == 'Tube&Fins':
                
                fact_cond_1 = np.log(self.params['Tube_OD']/(self.params['Tube_OD'] - 2*self.params['Tube_t']))
                fact_cond_2 = 2*np.pi*self.params['Tube_cond']*self.params['Tube_L']*self.params['n_tubes']*self.params['n_series']
                R_cond = fact_cond_1/fact_cond_2                      
                
                self.A_in_tubes =  self.params['A_finned'] # self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t'])) # *self.geom.Tube_pass 
                self.A_out_tubes = self.params['A_finned']
                
                R_fouling = 0 # self.geom.fouling / self.A_out_tubes             
                        
                # In the equation below, thickness resistance is given with respect to A_h arbitrarely
            
                if self.params['Fin_Side'] == 'H': # Fin side is the hot side 
                    UA_avail =  1/(1/(alpha_c*self.A_in_tubes) + 1/(alpha_h*self.A_out_tubes) + R_fouling + R_cond) # 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.tube_cond)) 
                else: 
                    UA_avail =  1/(1/(alpha_h*self.A_in_tubes) + 1/(alpha_c*self.A_out_tubes) + R_fouling + R_cond) # 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.tube_cond)) 

            elif self.HTX_Type == 'PCHE':
                
                if 'foul' in self.params:
                    R_fouling = self.params['foul'] / self.A_in_tubes
                else: 
                    R_fouling = 0
                
                self.A_c = 1/(1+self.params['R_p'])*self.params['N_c']*self.params['N_p']*(np.pi/2)*self.params['D_c']*self.params['L_c']
                self.A_h = self.params['R_p']/(1+self.params['R_p'])*self.params['N_c']*self.params['N_p']*(np.pi/2)*self.params['D_c']*self.params['L_c']
                
                # self.t_e = ((self.params['D_c'] + self.params['t_3'])*(self.params['D_c']/2 + self.params['t_2']) - (1/8*np.pi*self.params['D_c']**2))/(self.params['D_c'] + self.params['t_3'])
                self.t_e = self.params['t_3'] - np.pi*self.params['D_c']/8
                
                R_cond = self.t_e/self.params['k_cond']
                
                UA_avail =  1/(1/(alpha_h*self.A_h) + 1/(alpha_c*self.A_c) + R_fouling + R_cond)

                
            "5) Calculate w, the main objective of this function. This variable serves for residual minimization in the solver"      

            if k > len(self.w)- 1:
                self.w = np.append(self.w, UA_req/UA_avail) 
            else:
                self.w[k] = UA_req/UA_avail                
                
            # w.append(UA_req/UA_avail)
            # self.w = w
            
            w = self.w
            
            if self.HTX_Type == 'Plate':
                self.Avec_h.append(w[k]*self.params['A_h']) 
                
            if self.HTX_Type == 'Shell&Tube':
                self.Avec_h.append(w[k]*self.params['A_eff'])

            if self.HTX_Type == 'Tube&Fins':
                self.Avec_h.append(w[k]*self.params['A_finned'])
     
            if self.HTX_Type == 'PCHE':
                self.Avec_h.append(w[k]*self.A_h)       
     
            "5.1) Calculate dryout incipience only if subcritical"
            if not self.dry_out_c and self.phases_h[k] == "two-phase" and not self.Transcritical_c:
                    if not self.ColdSide_Props_Calculated:
                        #If these properties where not calculated before, do it:
                        self.AS_C.update(CP.PQ_INPUTS, p_c_mean, 0)
                        mu_c_l = self.AS_C.viscosity()
                        rho_c_l = self.AS_C.rhomass()
                        h_sat_c_l = self.AS_C.hmass()

                        self.AS_C.update(CP.PQ_INPUTS, p_c_mean, 1)
                        rho_c_v = self.AS_C.rhomass()
                        h_sat_c_v = self.AS_C.hmass()
                        i_fg_c = h_sat_c_v - h_sat_c_l
                    
                    P_star_c = p_c_mean/self.AS_C.p_critical()
                    q_c = self.Qvec_c[k]/self.Avec_h[k]
                    try:
                        self.AS_C.update(CP.PQ_INPUTS, p_c_mean, 0)
                        sigma_c_l = self.AS_C.surface_tension()
                        
                        if self.HTX_Type == 'Plate':
                            self.x_di_c = kim_dry_out_incipience(G_c, q_c, self.params['C_Dh'], P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
                            
                        elif self.HTX_Type == 'Shell&Tube':
                            self.x_di_c = kim_dry_out_incipience(G_c, q_c,self.params['Tube_OD']-2*self.params['Tube_t'], P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
                    except:
                        self.x_di_c = -1
        
        if k < len(self.w) - 1:
            self.w = self.w[:-1]
                
                        
        if debug:
            print(Q, 1-sum(w))
                
        return 1-sum(w)

#%% 
    def solve_hx(self):
        """ 
        Solve the objective function using Brent's method and the maximum heat transfer 
        rate calculated from the pinching analysis
        """
        self.Q = scipy.optimize.brentq(self.objective_function, 1e-5, self.Qmax-1e-10, rtol = 1e-14, xtol = 1e-10)
        # print('OUT of evap', self.Q)
        return self.Q
    
#%% 
    def plot_objective_function(self, N = 100):
        """ Plot the objective function """
        Q = np.linspace(1e-5,self.Qmax,N)
        r = np.array([self.objective_function(_Q) for _Q in Q])
        fig, ax = plt.subplots()
        ax.plot(Q, r)
        ax.grid(True)
        # plt.show()
        
#%%
    def plot_ph_pair(self):
        import warnings
        warnings.simplefilter("ignore")
        """ Plot p-h plots for the pair of working fluids """
        Ph_diagram_h = PropertyPlot(self.H_su.fluid, "PH", unit_system = "EUR")
        plt.plot(0.001*np.array(self.hvec_h), self.pvec_h*1e-05,'s-')
        Ph_diagram_h.calc_isolines()
        Ph_diagram_h.title("Hot side P-h diagram. Fluid: " + str(self.H_su.fluid))
        Ph_diagram_h.show()
        #----------------------------------------------------------------------
        Ph_diagram_c = PropertyPlot(self.C_su.fluid, "PH", unit_system = "EUR")
        plt.plot(0.001*np.array(self.hvec_c), self.pvec_c*1e-05,'s-')
        Ph_diagram_c.calc_isolines()
        Ph_diagram_c.title("Cold side P-h diagram. Fluid: " + str(self.C_su.fluid))
        Ph_diagram_c.show()
        warnings.simplefilter("default")

#%%
    def plot_Ts_pair(self):
        import warnings
        warnings.simplefilter("ignore")
        """ Plot a T-s plot for the pair of working fluids """
        Ph_diagram_h = PropertyPlot(self.H_su.fluid, "TS", unit_system = "EUR")
        plt.plot(0.001*self.svec_h,self.Tvec_h - 273.15,'s-')
        Ph_diagram_h.calc_isolines()
        Ph_diagram_h.title("Hot side T-s diagram. Fluid: " + str(self.H_su.fluid))
        Ph_diagram_h.show()
        # #----------------------------------------------------------------------
        Ph_diagram_c = PropertyPlot(self.C_su.fluid, "TS", unit_system = "EUR")
        plt.plot(0.001*self.svec_c,self.Tvec_c - 273.15,'s-')
        Ph_diagram_c.calc_isolines()
        Ph_diagram_c.title("Cold side T-s diagram. Fluid: " + str(self.C_su.fluid))
        Ph_diagram_c.show()
        warnings.simplefilter("default")

#%%
    def plot_cells(self, fName = '', dpi = 400):
        """ Plot the cells of the heat exchanger """
        plt.figure(figsize = (4,3))
        plt.plot(self.hnorm_h, self.Tvec_h, 'rs-')
        plt.plot(self.hnorm_c, self.Tvec_c, 'bs-')
        plt.xlim(0,1)
        plt.ylabel('T [K]') 
        plt.xlabel(r'$\hat h$ [-]')
        plt.grid(True)
        plt.tight_layout(pad = 0.2)
        plt.show()
        if fName != '':
            plt.savefig(fName, dpi = dpi)

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")
        print("======================")

