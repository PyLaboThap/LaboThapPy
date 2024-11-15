"""
Supplemental code for paper:
I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of 
Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014
"""

"""
Modification w/r to previous version:
    - Putting some order in the Objective Function "for" loops. Sparing some
    lines of code.
    - x_di_c correct calculation.
    - Implementation of various correlations for htc and DP
    - Implementation of total DP linearly interpolated on all discretizations
    - Implemenation of the code in the LaboThap Python Library
"""

# External Toolbox 
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d


# Internal Toolbox 
import __init__
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.f_lmtd2 import f_lmtd2
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.propsfluid import propsfluid
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.find_2P_boundaries import find_2P_boundaries

# HTC and DP Correlations
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.plate_htc import han_BPHEX_DP, water_plate_HTC, martin_BPHEX_HTC, muley_manglik_BPHEX_HTC, han_boiling_BPHEX_HTC, han_cond_BPHEX_HTC
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.pipe_htc import gnielinski_pipe_htc, boiling_curve, horizontal_tube_internal_condensation
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.shell_and_tube_htc import shell_bell_delaware_htc
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.tube_bank_htc import ext_tube_film_condens
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.fins import htc_tube_and_fins

# Phase related correlations
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.void_fraction import void_fraction
from component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.modules.kim_dry_out_incipience import kim_dry_out_incipience

# Connectors
from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

# Component base frame
from component.base_component import BaseComponent

#%%
# Set to True to enable some debugging output to screen
debug = False

#%%

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

            su_H_p: Hot suction side pressure. [Pa]

            su_H_h: Hot suction side enthalpy. [J/kg]

            su_H_fluid: Hot suction side fluid. [-]

            su_H_m_dot: Hot suction side mass flowrate. [kg/s]

            su_C_p: Cold suction side pressure. [Pa]

            su_C_h: Cold suction side enthalpy. [J/kg]

            su_C_fluid: Cold suction side fluid. [-]

            su_C_m_dot: Cold suction side mass flowrate. [kg/s]

        **Ouputs**:

            ex_H_h: Hot exhaust side specific enthalpy. [J/kg]

            ex_H_p: Hot exhaust side pressure. [Pa]

            ex_C_h: Cold exhaust side specific enthalpy. [J/kg]

            ex_C_p: Cold exhaust side pressure. [Pa]

            Q_dot: Heat transfer rate [W]

            M_H : Hot fluid charge [kg]

            M_C : Cold fluid charge [kg]
    """

    class H():
        pass
    class C():
        pass  
    
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
        
        self.Q_dot = HeatConnector()
        
        if HTX_Type == 'Plate' or HTX_Type == 'Shell&Tube' or HTX_Type == 'Tube&Fins':
            self.HTX_Type = HTX_Type
        else:
            print("Heat exchanger types implemented for this model are : 'Plate', 'Shell&Tube', 'Tube&Fins'.")

    #%%    
    
    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        """
        Hot side required inputs : 
            
            - Hsu_T or Hsu_h : Hot supply temperature or enthalpy
            - Hsu_p          : Hot supply pressure
            - Hsu_fluid      : Hot supply fluid
            - Hsu_m_dot      : Hot supply flow rate
            
        Cold side required inputs : 
            
            - Csu_T or Hsu_h : Cold supply temperature or enthalpy
            - Csu_p          : Cold supply pressure
            - Csu_fluid      : Cold supply fluid
            - Csu_m_dot      : Cold supply flow rate
        """
        self.sync_inputs()
        # Return a list of required inputs
        return['Hsu_p', 'Hsu_T', 'Hsu_m_dot', 'Hsu_fluid', 'Csu_p', 'Csu_T', 'Csu_m_dot', 'Csu_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        # Hot Fluid
        if self.su_H.T is not None:
            self.inputs['Hsu_T'] = self.su_H.T
        elif self.su_H.h is not None:
            self.inputs['Hsu_h'] = self.su_H.h
        if self.su_H.p is not None:
            self.inputs['Hsu_p'] = self.su_H.p
        if self.su_H.fluid is not None:
            self.inputs['Hsu_fluid'] = self.su_H.fluid
        if self.su_H.m_dot is not None:
            self.inputs['Hsu_m_dot'] = self.su_H.m_dot
            
        # Cold Fluid                
        if self.su_C.T is not None:
            self.inputs['Csu_T'] = self.su_C.T
        elif self.su_C.h is not None:
            self.inputs['Csu_h'] = self.su_C.h
        if self.su_C.p is not None:
            self.inputs['Csu_p'] = self.su_C.p
        if self.su_C.fluid is not None:
            self.inputs['Csu_fluid'] = self.su_C.fluid
        if self.su_C.m_dot is not None:
            self.inputs['Csu_m_dot'] = self.su_C.m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        # Hot Fluid
        self.su_H.set_fluid(self.inputs['Hsu_fluid'])
        if 'Hsu_T' in self.inputs:
            self.su_H.set_T(self.inputs['Hsu_T'])
        elif 'Hsu_h' in self.inputs:
            self.su_H.set_h(self.inputs['Hsu_h'])
        if 'Hsu_p' in self.inputs:
            self.su_H.set_p(self.inputs['Hsu_p'])
        if 'Hsu_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['Hsu_m_dot'])

        # Cold Fluid
        self.su_C.set_fluid(self.inputs['Csu_fluid'])
        if 'Csu_T' in self.inputs:
            self.su_C.set_T(self.inputs['Csu_T'])
        elif 'Csu_h' in self.inputs:
            self.su_C.set_h(self.inputs['Csu_h'])
        if 'Csu_p' in self.inputs:
            self.su_C.set_p(self.inputs['Csu_p'])
        if 'Csu_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['Csu_m_dot'])

        return['fluid_wf', 'su_wf_h', 'su_wf_m_dot', 'fluid_sf', 'su_sf_T', 'su_sf_cp', 'su_sf_m_dot']

#%%

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
        general_parameters = ['Flow_Type','htc_type', 'H_DP_ON', 'C_DP_ON','n_disc']
        
        if self.HTX_Type == 'Plate':     
            geometry_parameters = ['A_c', 'A_h', 'h', 'l', 'l_v',
                                   'C_CS', 'C_Dh', 'C_V_tot', 'C_canal_t', 'C_n_canals',
                                   'H_CS', 'H_Dh', 'H_V_tot', 'H_canal_t', 'H_n_canals',
                                   'casing_t', 'chevron_angle', 'fooling', 
                                   'n_plates', 'plate_cond', 'plate_pitch_co', 't_plates', 'w']

        elif self.HTX_Type == 'Shell&Tube':
                
            geometry_parameters = ['A_eff', 'Baffle_cut', 'D_OTL', 'N_strips', 'S_V_tot',
                                   'Shell_ID', 'T_V_tot', 'Tube_L', 'Tube_OD', 'Tube_pass',
                                   'Tube_t', 'Tubesheet_t', 'central_spacing', 'clear_BS', 'clear_TB',
                                   'cross_passes', 'foul_s', 'foul_t', 'inlet_spacing', 'n_series',
                                   'n_tubes', 'outlet_spacing', 'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']
        
        elif self.HTX_Type == 'Tube&Fins':
            
            geometry_parameters = ['A_finned', 'A_flow', 'A_in_tot', 'A_out_tot', 'A_unfinned',
                                   'B_V_tot', 'Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type',
                                   'Finned_tube_flag', 'L', 'T_V_tot', 'Tube_L', 'Tube_OD',
                                   'Tube_cond', 'Tube_t', 'fouling', 'h', 'k_fin',
                                   'n_passes', 'n_rows', 'n_tubes', 'pitch', 'pitch_ratio', 'tube_arrang',
                                   'w','Fin_Side']
        
        return general_parameters + geometry_parameters
    
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
        
        self.check_calculable()
        
        if htc_type == "User-Defined":
            
            self.H.HeatExchange_Correlation = "User-Defined"
            self.C.HeatExchange_Correlation = "User-Defined"
            
            # !!! User-Defined Heat Transfer Coefficients (hot):
            self.H.h_liq = UD_H_HTC['Liquid']
            self.H.h_vap = UD_H_HTC['Vapor']
            self.H.h_twophase = UD_H_HTC['Two-Phase']
            self.H.h_vapwet = UD_H_HTC['Vapor-wet']
            self.H.h_tpdryout = UD_H_HTC['Dryout']
            self.H.h_transcrit = UD_H_HTC['Transcritical']
        
            # !!! User-Defined Heat Transfer Coefficients (cold):
            self.C.h_liq = UD_C_HTC['Liquid']
            self.C.h_vap = UD_C_HTC['Vapor']
            self.C.h_twophase = UD_C_HTC['Two-Phase']
            self.C.h_vapwet = UD_C_HTC['Vapor-wet']
            self.C.h_tpdryout = UD_C_HTC['Dryout']
            self.C.h_transcrit = UD_C_HTC['Transcritical']
            
            if self.params['H_DP_ON'] == True: # !!! arbitrarily set
                self.H.f_dp = {"K": 14.14, "B": 1.892}
            else:
                self.H.f_dp = {"K": 0, "B": 0}        
                    
            if self.params['C_DP_ON'] == True:
                self.C.f_dp = {"K": 14.14, "B": 1.892}
            else:
                self.C.f_dp = {"K": 0, "B": 0}
            
        else:
            # Type 
            self.H.HeatExchange_Correlation = "Correlation"
            self.C.HeatExchange_Correlation = "Correlation"
            
            if self.HTX_Type == 'Plate' or self.HTX_Type == 'Shell&Tube' or self.HTX_Type == 'Tube&Fins':
                
                self.H.Correlation_1phase = Corr_H["1P"]
                self.H.Correlation_2phase = Corr_H["2P"]
                
                self.C.Correlation_1phase = Corr_C["1P"]
                self.C.Correlation_2phase = Corr_C["2P"]

                if self.params['H_DP_ON'] == True: # !!! arbitrarily set
                    self.H.f_dp = {"K": 14.14, "B": 1.892}
                else:
                    self.H.f_dp = {"K": 0, "B": 0}        
                    
                if self.params['C_DP_ON'] == True:
                    self.C.f_dp = {"K": 14.14, "B": 1.892}
                else:
                    self.C.f_dp = {"K": 0, "B": 0}
                    
            if self.C.Correlation_2phase == "Boiling_curve": # Compute the fluid boiling curve beforehand
                try:
                    T_sat = CP.PropsSI('T','P', self.su_C.p,'Q',0,self.su_C.fluid)
                    (h_boil, DT_vect) = boiling_curve(self.params['Tube_OD'], self.su_C.fluid, T_sat, self.su_C.p)
                    self.C_f_boiling = interp1d(DT_vect,h_boil)
                except:
                    self.C_f_boiling = interp1d([0,10000],[20000,20000])
            
#%%
    def external_pinching(self):
        "Determine the maximum heat transfer rate based on the external pinching analysis"
        
        "1) Set of arbitrary bound values" # !!! Find out why      
        T_hmin = 218 # above CO2 freezing point (194.7 K)
            
        T_cmax = 273.15+250
        
        "2) Hot fluid side pinch"
        
        # Computation of outlet lowest possible enthalpy of the hot fluid using either the entrance cold fluid inlet temperature, or the arbitrary minimal
        self.h_ho = CP.PropsSI('H','T', max(self.T_ci, T_hmin),'P',self.p_hi,self.H_su.fluid) # Equation 5 (Bell et al. 2015)
        
        # Q_max computation
        Qmaxh = self.mdot_h*(self.h_hi-self.h_ho) # Equation 4 (Bell et al. 2015)
        
        if debug:
            print("Inlet Hot Enthalpy:", self.h_hi, "\n")
            print("Outlet Hot Enthalpy:", self.h_ho,"\n")
            print("mdot_h = ", self.mdot_h, "\n")
            print("Qmaxh =", Qmaxh, "\n")

        "3) Cold fluid side pinch"
                
        # Computation of highest possible outlet enthalpy of the hot fluid using either the entrance hot fluid inlet temperature, or the arbitrary maximal
        self.h_co = CP.PropsSI('H','T',min(self.T_hi, T_cmax),'P',self.p_ci,self.C_su.fluid) # Equation  (Bell et al. 2015)
        
        # Q_max computation
        Qmaxc = self.mdot_c*(self.h_co-self.h_ci) # Equation 6 (Bell et al. 2015)

        if debug:
            print("Inlet Cold Enthalpy:", self.h_ci, "\n")
            print("Outlet Cold Enthalpy:", self.h_co,"\n")
            print("mdot_c = ", self.mdot_c, "\n")
            print("Qmaxh =", Qmaxc, "\n")
            
        "4) Final pinch and cell boundaries computation"
        
        Qmax = min(Qmaxh, Qmaxc)
        
        if debug:
            print('Qmax (external pinching) is', Qmax)

        self.calculate_cell_boundaries(Qmax) # call calculate_cell_boundaries procedure

        return Qmax
    
#%%
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
                self.T_cdew = CP.PropsSI("T", "P", self.p_cdew, "Q", 1, self.C_su.fluid)
                self.T_cbubble = CP.PropsSI("T", "P", self.p_cbubble, "Q", 0, self.C_su.fluid)
                
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
                self.T_hdew = CP.PropsSI("T", "P", self.p_hdew, "Q", 1, self.H_su.fluid)
                self.T_hbubble = CP.PropsSI("T", "P", self.p_hbubble, "Q", 0, self.H_su.fluid)

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

            if round(Qcell_hk, 4) > round(Qcell_ck, 4):
                # Hot stream needs a complementary cell boundary as the available heat is higher than the cold cell heat absorption
                self.hvec_h.insert(k+1, self.hvec_h[k] + Qcell_ck/self.mdot_h)
                
            elif round(Qcell_hk, 4) < round(Qcell_ck, 4):
                # Cold stream needs a complementary cell boundary as the available heat absoprtion is higher than the hot cell heat 
                self.hvec_c.insert(k+1, self.hvec_c[k] + Qcell_hk/self.mdot_c)

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
        self.pvec_h = (1-self.hnorm_h)*self.p_hi + self.hnorm_h*self.p_ho

        # Calculate the temperature and entropy at each cell boundary
        self.Tvec_c = CP.PropsSI('T','H',self.hvec_c,'P',self.pvec_c,self.C_su.fluid)
        
        self.Tvec_h = CP.PropsSI('T','H',self.hvec_h,'P',self.pvec_h,self.H_su.fluid)
        self.svec_h = CP.PropsSI('S','H',self.hvec_h,'P',self.pvec_h,self.H_su.fluid)
        
        self.svec_c = CP.PropsSI('S','H',self.hvec_c,'P',self.pvec_c,self.C_su.fluid)


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
        
        if not self.Transcritical_c:
            self.Tvec_sat_pure_c = CP.PropsSI("T", "Q", 0.5, "P", self.pvec_c, self.C_su.fluid)
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.Tvec_sat_pure_h = CP.PropsSI("T", "Q", 0.5, "P", self.pvec_h, self.H_su.fluid)
            
        "8) Calculate pinch and heat exchanger border temperature deltas"
        
        self.DT_pinch = np.min((np.array(self.Tvec_h) - np.array(self.Tvec_c)))
        self.DT_ho_ci = self.Tvec_h[0] - self.Tvec_c[0]
        self.DT_hi_co = self.Tvec_h[-1] - self.Tvec_c[-1]

#%%
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
                        h_c_pinch = CP.PropsSI('H','T',self.T_hdew,'P',self.pvec_c[i], self.C_su.fluid) # Equation 10 (Bell et al. 2015)

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
                        
                        # Enthalpy of the cold stream at the pinch temperatur
                        h_h_pinch = CP.PropsSI('H','T',self.T_cbubble,'P',self.pvec_h[i], self.H_su.fluid) # Equation 14 (Bell et al. 2015)

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
    
    #%%
    def solve(self, only_external = False, and_solve = True):
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
            
        self.h_incomp_flag = (self.H_su.fluid.find("INCOMP") != -1)
                
        # Hot fluid
        self.mdot_h = self.H_su.m_dot
        self.h_hi = self.H_su.h
        self.p_hi = self.H_su.p
        
        # Cold fluid 
        self.mdot_c = self.C_su.m_dot
        self.h_ci = self.C_su.h
        self.p_ci = self.C_su.p
                
        # Determine the inlet temperatures from the pressure/enthalpy pairs
        self.T_ci = CP.PropsSI('T', 'P', self.p_ci, 'H', self.h_ci, self.C_su.fluid)
        self.T_hi = CP.PropsSI('T', 'P', self.p_hi, 'H', self.h_hi, self.H_su.fluid)

        "2) Determine if the streams come in at a higher pressure than the transcritical pressure"
        
        # Initialization of the flags:    
        self.Transcritical_c = False
        self.Transcritical_h = False
        
        if not self.h_incomp_flag and (self.p_hi - CP.PropsSI("PCRIT", self.H_su.fluid)) >= 1e-06:
            self.Transcritical_h = True

        "3) Calculate the ideal bubble and dew temperatures/enthalpies for each stream IF the fluid is not transcritical"
        
        
        if not self.Transcritical_c:
            self.T_cbubble_ideal = CP.PropsSI('T', 'P', self.p_ci, 'Q', 0, self.C_su.fluid)
            self.T_cdew_ideal    = CP.PropsSI('T', 'P', self.p_ci, 'Q', 1, self.C_su.fluid)
            self.h_cbubble_ideal = CP.PropsSI('H', 'T', self.T_cbubble_ideal, 'Q', 0, self.C_su.fluid)
            self.h_cdew_ideal    = CP.PropsSI('H', 'T', self.T_cdew_ideal, 'Q', 1, self.C_su.fluid)  
            
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.T_hbubble_ideal = CP.PropsSI('T', 'P', self.p_hi, 'Q', 0, self.H_su.fluid)
            self.T_hdew_ideal    = CP.PropsSI('T', 'P', self.p_hi, 'Q', 1, self.H_su.fluid)
            self.h_hbubble_ideal = CP.PropsSI('H', 'T', self.T_hbubble_ideal, 'Q', 0, self.H_su.fluid)
            self.h_hdew_ideal    = CP.PropsSI('H', 'T', self.T_hdew_ideal, 'Q', 1, self.H_su.fluid)
            
        "4) Calculate pressure drops"
        
        if self.params['H_DP_ON'] == True: # if the pressure drop are not neglected
            if self.H_su.fluid == 'Water':
                self.DP_h = self.H.f_dp["K"] * (self.mdot_h/self.H_su.D)**2
                self.p_ho = self.p_hi - self.DP_h
            else:
                # rho_v = CP.PropsSI('D','Q',1,'P',self.p_hi,self.H_su.fluid)
                # rho_l = CP.PropsSI('D','Q',0,'P',self.p_hi,self.H_su.fluid)
                # mu_h = CP.PropsSI('V','H',self.h_hi,'P',self.p_hi,self.H_su.fluid)
                
                # G_h = (self.mdot_h/self.geom.H_n_canals)/self.geom.H_CS
                self.DP_h = 0 # self.H.f_dp["K"] * self.mdot_h**(self.H.f_dp["B"]) # Empirical correlations : DP = K*m_dot**B # Han_BPHEX_DP(mu_h, G_h, self.geom.H_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co, rho_v, rho_l, self.geom.l_v, self.geom.H_n_canals, self.mdot_h, self.geom.H_canal_t) # 
                self.p_ho = self.p_hi - self.DP_h
            
        if self.C.f_dp != {"K": 0, "B": 0}:
            if self.C_su.fluid == 'Water':
                self.DP_c = self.C.f_dp["K"] * (self.mdot_c/self.C_su.D)**2
                self.p_co = self.p_ci - self.DP_c
            else:
                # rho_v = CP.PropsSI('D','Q',1,'P',self.p_ci,self.C_su.fluid)
                # rho_l = CP.PropsSI('D','Q',0,'P',self.p_ci,self.C_su.fluid)
                # mu_c = CP.PropsSI('V','H',self.h_ci,'P',self.p_ci,self.C_su.fluid)               
                
                # G_c = (self.mdot_c/self.geom.C_n_canals)/self.geom.C_CS
                self.DP_c = 0 # self.C.f_dp["K"] * self.mdot_c**(self.C.f_dp["B"]) # Han_BPHEX_DP(mu_c, G_c, self.geom.H_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co, rho_v, rho_l, self.geom.l_v, self.geom.C_n_canals, self.mdot_c, self.geom.C_canal_t) # 
                self.p_co = self.p_ci - self.DP_c
            
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
                    self.Dvec_h[i] = CP.PropsSI('D', "P", self.pvec_h[i], "H", self.hvec_h[i], self.H_su.fluid)
                    self.eps_void_h[i] = -1
                else:
                    rho_g = CP.PropsSI('D', 'P', self.pvec_h[i], "Q", 1, self.H_su.fluid)
                    rho_l = CP.PropsSI('D', 'P', self.pvec_h[i], "Q", 0, self.H_su.fluid)
                    self.eps_void_h[i], self.Dvec_h[i] = void_fraction(self.x_vec_h[i], rho_g, rho_l)
                    
            for i in range(len(self.hvec_c)-1):
                if self.x_vec_c[i] <= 0 or self.x_vec_c[i] >= 1: 
                    self.Dvec_c[i] = CP.PropsSI('D', "P", self.pvec_c[i], "H", self.hvec_c[i], self.C_su.fluid)
                    self.eps_void_c[i] = -1
                else:
                    rho_g = CP.PropsSI('D', 'P', self.pvec_c[i], "Q", 1, self.C_su.fluid)
                    rho_l = CP.PropsSI('D', 'P', self.pvec_c[i], "Q", 0, self.C_su.fluid)
                    self.eps_void_c[i], self.Dvec_c[i] = void_fraction(self.x_vec_c[i], rho_g, rho_l)
            
            "5.5) Computation of the fluid mass inside the HTX"
            
            if self.HTX_Type == 'Plate':
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
                self.ex_H.set_T(self.Tvec_h[0])
                self.ex_H.set_p(self.pvec_h[0])
                self.ex_H.set_m_dot(self.H_su.m_dot)
                
                self.H_ex = self.ex_H
                    
                # Cold side
                self.ex_C.set_fluid(self.C_su.fluid)
                self.ex_C.set_T(self.Tvec_c[-1]) # Example temperature [K]
                self.ex_C.set_p(self.pvec_c[-1]) # Example Pressure [Pa]
                self.ex_C.set_m_dot(self.C_su.m_dot) 
        
                self.C_ex = self.ex_C
                                            
                self.defined = True
                return Q
            
        else: # Just a flag if the heat exchanger is not solved
            self.Q = 1
            
#%%    
    def objective_function(self, Q):
        
        "1) Initialize cell boundaries and results vectors"
        
        # Cell boundaries
        self.calculate_cell_boundaries(Q)

        # Initialize dry-out incipience identification
        self.x_di_c = 1
        self.dry_out_c = False
        self.ColdSide_Props_Calculated = False
        
        # Initialise results arrays
        w = [] # The wet-bulb potential temperature
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
                LMTD = DTA
            else:
                try:
                    LMTD = (DTA-DTB)/np.log(abs(DTA/DTB))
                except ValueError as VE:
                    print(Q, DTA, DTB)
                    raise
            
            # Wall Temperature:
            T_wall = (Thi + Tho + Tci + Tco)/4
            
            # Cold side
            Tc_mean = 0.5*(Tci + Tco) # mean temperature over the cell
            
            if not self.Transcritical_c:
                Tc_sat_mean = 0.5*(self.Tvec_sat_pure_c[k] + self.Tvec_sat_pure_c[k+1]) # mean saturation temperature over the cell
                
            p_c_mean = 0.5*(self.pvec_c[k] + self.pvec_c[k+1]) # mean pressure over the cell
            
            # Hot side
            Th_mean = 0.5*(Thi + Tho) # mean temperature over the cell
            
            if not self.Transcritical_h and not self.h_incomp_flag:
                Th_sat_mean = 0.5*(self.Tvec_sat_pure_h[k] + self.Tvec_sat_pure_h[k+1]) # mean saturation temperature over the cell
            
            p_h_mean = 0.5*(self.pvec_h[k] + self.pvec_h[k+1]) # mean pressure over the cell

            "2.2) Hot side phase identification"
            
            # If not transcritical
            if not self.Transcritical_h and not self.h_incomp_flag:
                
                havg_h = (self.hvec_h[k] + self.hvec_h[k+1])/2.0 # Average cell enthalpy over the cell
                
                if havg_h < self.h_hbubble: # below evaporation/bubble point
                    self.phases_h.append('liquid')
                    T_wall_h = T_wall
                    
                elif havg_h > self.h_hdew: # higher than condensation/dew point
                    #T_wall_h not used so far but it could appear in future Correlations implemented.
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
                
            elif self.Transcritical_h:
                self.phases_h.append("transcritical")
                T_wall_h = T_wall
            
            "2.3) Cold side phase identification"

            # If not transcritical
            if not self.Transcritical_c:
                
                havg_c = (self.hvec_c[k] + self.hvec_c[k+1])/2.0 # Average cell enthalpy over the cell
                
                if havg_c < self.h_cbubble: # below evaporation/bubble point
                    self.phases_c.append('liquid')
                    #The line comented below is the original line on the code. It is causing trouble.
                    #T_wall_c = min(T_wall, Tc_sat_mean)
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
                        
            elif self.Transcritical_c:
                self.phase_c.append('transcritical')
            
            if self.HTX_Type == 'Plate':                
                G_c = (self.mdot_c/self.params['C_n_canals'])/self.params['C_CS']
                G_h = (self.mdot_h/self.params['H_n_canals'])/self.params['H_CS']
                
            elif self.HTX_Type == 'Shell&Tube':
                A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            
                # Mass flow in tubes if the fluid is in the tubes
                G_h = (self.mdot_h/self.params['n_tubes'])/A_in_one_tube
                G_c = (self.mdot_c/self.params['n_tubes'])/A_in_one_tube

            elif self.HTX_Type == 'Tube&Fins':
                A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            
                # Mass flow in tubes if the fluid is in the tubes
                G_h = (self.mdot_h/self.params['n_tubes'])/A_in_one_tube
                G_c = (self.mdot_c/self.params['n_tubes'])/A_in_one_tube
                
            # F correction factor for LMTD method:
            if self.params['Flow_Type'] != "CounterFlow":
                
                try:
                    C_c = CP.PropsSI('C', 'T', Tc_mean,'P',p_c_mean,self.C_su.fluid)
                except:
                    C_c = 20000
                    
                try:
                    C_h = CP.PropsSI('C', 'T', Th_mean,'P',p_h_mean,self.H_su.fluid)
                except:
                    C_h = 20000
                
                C_min = min(C_c,C_h)
                C_max = max(C_c,C_h)
                
                C_r = C_min/C_max
                
                self.R = (Thi - Tho)/(Tco - Tci)
                self.P = (Tco - Tci)/(Thi - Tci)
                self.F = f_lmtd2(self.R, self.P, self.params, C_r)
            else:
                self.F = 1
            
            #UA_req including F correction factor in case of cross flow
            UA_req = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])/(self.F*LMTD)          
            
            "3) Cell heat transfer coefficients"
            
            "3.1) Hot side - User defined"
            
            if self.H.HeatExchange_Correlation == "User-Defined":
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
                    
            elif self.H.HeatExchange_Correlation == "Correlation": 
                
                # Heat transfer coefficient calculated from Correlations.
                # Evaluate phases
                
                # 1 phase case
                if self.phases_h[k] == "liquid" or self.phases_h[k] == "vapor" or self.phases_h[k] == "transcritical":
                                        
                    if self.H.Correlation_1phase == "Gnielinski":
                        try:
                            mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = propsfluid(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, False)
                            Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', Th_mean,'P',p_h_mean,self.H_su.fluid)
                        except (ValueError):
                            if self.phases_h[k] == "liquid":
                                mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = propsfluid(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, False)
                                Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', Th_mean-1,'P',p_h_mean,self.H_su.fluid)
                            elif self.phases_h[k] == "vapor":
                                mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = propsfluid(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, False)
                                Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', Th_mean+1,'P',p_h_mean,self.H_su.fluid)

                        if self.HTX_Type == 'Plate':
                            alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_h, mu_h_w, Pr_h, k_h, G_h, self.geom.H_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_h, Pr_h, k_h, G_h, self.geom.H_Dh) # 
                        
                        elif self.HTX_Type == 'Shell&Tube':
                            alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                        
                        elif self.HTX_Type == 'Tube&Fins':
                            alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['n_passes']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                    
                    elif self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC":
                        alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean, T_wall, p_h_mean, self.H_su.fluid, self.params)
                        
                    elif self.H.Correlation_1phase == 'Tube_And_Fins':
                        alpha_h = htc_tube_and_fins(self.H_su.fluid, self.params, p_h_mean, Th_mean, self.mdot_h, self.params['Fin_type'])[0]
                        
                # 2 phases case
                elif self.phases_h[k] == "two-phase" or self.phases_h[k] == "vapor-wet":
                    if self.phases_h[k] == "two-phase":
                        x_h = min(1, max(0, 0.5*(self.x_vec_h[k] + self.x_vec_h[k])))
                    elif self.phases_h[k] == "vapor-wet":
                        x_h = 1
                    
                    # Thermodynamical variables
                    mu_h_l = CP.PropsSI("V", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                    k_h_l = CP.PropsSI("L", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                    Pr_h_l = CP.PropsSI("Prandtl", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                    rho_h_l = CP.PropsSI("D", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                    rho_h_v = CP.PropsSI("D", "Q", 1, "P", p_h_mean, self.H_su.fluid)
                    # i_fg_h = CP.PropsSI('H', 'Q', 1, 'P', p_h_mean, self.H_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_h_mean, self.H_su.fluid)
                    
                    # !!! Include different types of Correlation HERE
                    if self.H.Correlation_2phase == "Han_cond_BPHEX":
                        alpha_h_2phase, _, DP_H = han_cond_BPHEX_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, rho_h_l, rho_h_v, G_h, self.params['H_Dh'], self.params['plate_pitch_co'], self.params['chevron_angle'], self.params['l_v'], self.params['H_n_canals'], self.H_su.m_dot, self.params['H_canal_t'])
                    
                    if self.H.Correlation_2phase == 'ext_tube_film_condens':

                        V_flow = G_h/CP.PropsSI('D','H',havg_h,'P',p_h_mean,self.H_su.fluid)
                        
                        if self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC":
                            try: 
                                alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean, T_wall, p_h_mean, self.H_su.fluid, self.params)
                            except:
                                alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean-0.1, T_wall, p_h_mean, self.H_su.fluid, self.params)
                            
                        elif self.H.Correlation_1phase == 'Tube_And_Fins':
                            alpha_h = htc_tube_and_fins(self.C_su.fluid, self.params, p_c_mean, Tc_mean, self.mdot_c, self.params['Fin_type'])[0]
                        
                        alpha_h_2phase = ext_tube_film_condens(self.params['Tube_OD'], self.H_su.fluid, Th_mean, T_wall, V_flow)
                        
                    if self.H.Correlation_2phase == 'Horizontal_Tube_Internal_Condensation':
                        mu_h, Pr_h, k_h = CP.PropsSI(('V','PRANDTL','L'), 'H', havg_h,'P',p_h_mean,self.H_su.fluid)
                        Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', T_wall,'P',p_h_mean,self.H_su.fluid)
                                
                        alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['Tube_OD'] - 2*self.params['Tube_t'], self.params['Tube_L'])
                        alpha_h_2phase = horizontal_tube_internal_condensation(self.H_su.fluid , G_c, p_h_mean, self.x_vec_h[k], T_wall, self.params['Tube_OD'] - 2*self.params['Tube_t'])
                                        
                    if self.phases_h[k] == "two-phase":
                        alpha_h = alpha_h_2phase
                        
                    elif self.phases_h[k] == "vapor-wet":
                        w_vap_wet = (Th_mean - Th_sat_mean)/(Th_mean - T_wall)
                        # The line below re-calculates alpha_h in case of having a vapor-wet condition
                        alpha_h = alpha_h_2phase - w_vap_wet*(alpha_h_2phase - alpha_h) # the last alpha_h in this equation is the 1 Phase calculation

            "3.2) Cold side - User defined"
            
            # User Defined
            if self.C.HeatExchange_Correlation == "User-Defined":
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
                    
            elif self.C.HeatExchange_Correlation == "Correlation":

                # Heat transfer coefficient calculated from Correlations:
                # Evaluate phases
                
                # 1 phase case
                if self.phases_c[k] == "liquid" or self.phases_c[k] == "vapor" or self.phases_c[k] == "transcritical":
                    mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False)
            
                    if self.HTX_Type == 'Plate' and self.C_su.fluid == 'water':
                        alpha_c = water_plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.H_Dh)

                    elif self.C.Correlation_1phase  == "Gnielinski":
                        try:
                            mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False)
                            Pr_c_w, mu_c_w = CP.PropsSI(('PRANDTL','V'), 'T', Tc_mean,'P',p_c_mean,self.C_su.fluid)
                        except (ValueError):
                            if self.phases_c[k] == "liquid":
                                mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean-0.1, p_c_mean, T_wall_c, self.C_su.fluid, False)
                                Pr_c_w, mu_c_w = CP.PropsSI(('PRANDTL','V'), 'T', Tc_mean-1,'P',p_c_mean,self.C_su.fluid)
                            elif self.phases_c[k] == "vapor":
                                mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean+0.1, p_c_mean, T_wall_c, self.C_su.fluid, False)
                                Pr_c_w, mu_c_w = CP.PropsSI(('PRANDTL','V'), 'T', Tc_mean+1,'P',p_c_mean,self.C_su.fluid)
                        
                        if self.HTX_Type == 'Plate':
                            alpha_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c_w, k_c, G_c, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                        
                        elif self.HTX_Type == 'Shell&Tube':
                            alpha_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                        
                        elif self.HTX_Type == 'Tube&Fins':
                            alpha_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['n_passes']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                            
                    elif self.C.Correlation_1phase == 'Shell_Bell_Delaware_HTC':
                        alpha_c = shell_bell_delaware_htc(self.mdot_c, Tc_mean, T_wall, p_c_mean, self.C_su.fluid, self.params)
                         
                    elif self.C.Correlation_1phase == 'Tube_And_Fins':
                        alpha_c = htc_tube_and_fins(self.C_su.fluid, self.params, p_c_mean, Tc_mean, self.mdot_c, self.params['Fin_type'])[0]
                    
                    # if self.C.Correlation_1phase  == "Gnielinski":
                    #     alpha_c, _ = Gnielinski_Pipe_HTC(mu_c, Pr_c, k_c, G_c, self.geom.H_Dh, self.geom.l)

                # 2 phases case
                elif self.phases_c[k] == "two-phase" or self.phases_c[k] == "vapor-wet":
                    if self.phases_c[k] == "two-phase":
                        x_c = min(1, max(0, 0.5*(self.x_vec_c[k] + self.x_vec_c[k])))
                    elif self.phases_c[k] == "vapor-wet":
                        x_c = 1
                        
                    # Thermodynamical variables
                    mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                    k_c_l = CP.PropsSI("L", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                    Pr_c_l = CP.PropsSI("Prandtl", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                    rho_c_l = CP.PropsSI("D", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                    rho_c_v = CP.PropsSI("D", "Q", 1, "P", p_c_mean, self.C_su.fluid)
                    i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                    
                    # This boolean serves to spare to recalculate this properties in the dry-out analysis further ahead.
                    self.ColdSide_Props_Calculated = True
                    
                    # !!! Include different types of Correlation HERE
                    if self.C.Correlation_2phase == "Han_cond_BPHEX":
                        alpha_c_2phase, _ = han_cond_BPHEX_HTC(x_c, mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v, G_c, self.params['C_Dh'], self.params['plate_pitch_co'], self.params['chevron_angle'])
                    
                    elif self.C.Correlation_2phase == "Han_Boiling_BPHEX_HTC":
                        alpha_c_2phase, _ = han_boiling_BPHEX_HTC(min(x_c, self.x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, LMTD*self.F, self.Qvec_c[k], alpha_h, self.params['C_Dh'], self.params['chevron_angle'], self.params['plate_pitch_co'])
                    
                    elif self.C.Correlation_2phase == "Boiling_curve":                  
                        alpha_c_2phase = self.C_f_boiling(abs(T_wall - Tc_mean))
                    
                    else:
                        raise ValueError("Correlation not found for Cold Side 2-Phase")
                            
                    if self.phases_c[k] == "two-phase":
                        alpha_c = alpha_c_2phase
                    elif self.phases_c[k] == "vapor-wet":
                        w_vap_wet = (Tc_mean - Tc_sat_mean)/(Tc_mean - T_wall)
                        #The line below re-calculates alpha_h in case of having a vapor-wet condition
                        alpha_c = alpha_c_2phase - w_vap_wet*(alpha_c_2phase - alpha_c) #the last alpha_c in this equation is the 1 Phase calculation
            
            "3.3) Store the heat transfer coefficients"
            
            self.alpha_c.append(alpha_c)
            self.alpha_h.append(alpha_h)
            
            "4) Recover some parameters and conductivity calculation and compute UA_avail"
            
            if self.HTX_Type == 'Plate':     
                # In the equation below, thickness resistance is given with respect to A_h arbitrarely
                UA_avail = 1/((1+self.params['fooling'])/(alpha_h*self.params['A_h']) + 1/(alpha_c*self.params['A_c']) + self.params['t_plates']/(self.params['plate_cond'])) 
                
            elif self.HTX_Type == 'Shell&Tube':       
                
                fact_cond_1 = np.log(self.params['Tube_OD']/(self.params['Tube_OD'] - 2*self.params['Tube_t']))
                fact_cond_2 = 2*np.pi*self.params['tube_cond']*self.params['Tube_L']*self.params['n_series']*self.params['n_tubes']
                R_cond = fact_cond_1/fact_cond_2                      

                self.A_in_tubes = self.params['Tube_pass']*self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t']))
                self.A_out_tubes = self.params['Tube_pass']*self.params['Tube_L']*self.params['n_tubes']*np.pi*(self.params['Tube_OD'])
            
                if self.params['foul_s'] != None:
                    R_fouling_s = self.params['foul_s'] / self.A_out_tubes
                else: 
                    R_fouling_s = 0
    
                if self.params['foul_t'] != None:
                    R_fouling_t = self.params['foul_t'] / self.A_in_tubes
                else: 
                    R_fouling_t = 0                
            
                UA_avail = 1/(1/(alpha_c*self.A_in_tubes) + 1/(alpha_h*self.A_out_tubes) + R_fouling_s + R_fouling_t + R_cond) 
            
            elif self.HTX_Type == 'Tube&Fins':
                
                fact_cond_1 = np.log(self.params['Tube_OD']/(self.params['Tube_OD'] - 2*self.params['Tube_t']))
                fact_cond_2 = 2*np.pi*self.params['Tube_cond']*self.params['Tube_L']*self.params['n_tubes']*self.params['n_passes']
                R_cond = fact_cond_1/fact_cond_2                      
                
                self.A_in_tubes = self.params['n_passes']*self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t'])) # *self.geom.n_passes 
                self.A_out_tubes = self.params['A_finned']
                            
                #☺ self.A_out_tubes = 8705 #self.geom.n_passes*self.geom.Tube_L*self.geom.n_tubes*np.pi*(self.geom.Tube_OD)
            
                R_fouling = 0 # self.geom.fouling / self.A_out_tubes             
                        
                # In the equation below, thickness resistance is given with respect to A_h arbitrarely
            
                if self.params['Fin_Side'] == 'H': # Fin side is the hot side 
                    UA_avail =  1/(1/(alpha_c*self.A_in_tubes) + 1/(alpha_h*self.A_out_tubes) + R_fouling + R_cond) # 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.tube_cond)) 
                else: 
                    UA_avail =  1/(1/(alpha_h*self.A_in_tubes) + 1/(alpha_c*self.A_out_tubes) + R_fouling + R_cond) # 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.tube_cond)) 
      
            # R_fooling_h = self.params.H.R_fooling 
            # R_fooling_c = self.params.C.R_fooling
            
            "5) Calculate w, the main objective of this function. This variable serves for residual minimization in the solver"
            
            w.append(UA_req/UA_avail)
            
            self.w = w
            
            if self.HTX_Type == 'Plate':
                self.Avec_h.append(w[k]*self.params['A_h']) 
                
            if self.HTX_Type == 'Shell&Tube':
                self.Avec_h.append(w[k]*self.params['A_eff'])

            if self.HTX_Type == 'Tube&Fins':
                self.Avec_h.append(w[k]*self.params['A_finned'])
            
            "5.1) Calculate dryout incipience only if subcritical"
            if not self.dry_out_c and self.phases_h[k] == "two-phase" and not self.Transcritical_c:
                    if not self.ColdSide_Props_Calculated:
                        #If these properties where not calculated before, do it:
                        mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        rho_c_l = CP.PropsSI('D', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                        rho_c_v = CP.PropsSI('D', 'Q', 1, 'P', p_c_mean, self.H_su.fluid)
                        i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                    
                    P_star_c = p_c_mean/CP.PropsSI("Pcrit", "Q", 1, "P", p_c_mean, self.C_su.fluid)
                    q_c = self.Qvec_c[k]/self.Avec_h[k]
                    try:
                        sigma_c_l = CP.PropsSI("I", 'Q', 0, "P", p_c_mean, self.C_su.fluid)
                        
                        if self.HTX_Type == 'Plate':
                            self.x_di_c = kim_dry_out_incipience(G_c, q_c, self.params['C_Dh'], P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
                            
                        elif self.HTX_Type == 'Shell&Tube':
                            self.x_di_c = kim_dry_out_incipience(G_c, q_c,self.params['Tube_OD']-2*self.params['Tube_t'], P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
                    except:
                        self.x_di_c = -1
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
        return self.Q
    
#%%
    def plot_objective_function(self, N = 100):
        """ Plot the objective function """
        Q = np.linspace(1e-5,self.Qmax,N)
        r = np.array([self.objective_function(_Q) for _Q in Q])
        fig, ax = plt.subplots()
        ax.plot(Q, r)
        ax.grid(True)
        # plt.show()²
        
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

# """
# Supplemental code for paper:
# I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of 
# Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014
# """

# """
# Modification w/r to previous version:
#     - Putting some order in the Objective Function "for" loops. Sparing some
#     lines of code.
#     - x_di_c correct calculation.
# """

# # from __future__ import division, print_function

# # import sys
# # sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')

# """
# Supplemental code for paper:
# I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of 
# Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014
# """

# """
# Modification w/r to previous version:
#     - Putting some order in the Objective Function "for" loops. Sparing some
#     lines of code.
#     - x_di_c correct calculation.
#     - Implementation of various correlations for htc and DP
#     - Implementation of total DP linearly interpolated on all discretizations
#     - Implemenation of the code in the LaboThap Python Library
# """

# # External Toolbox 
# import CoolProp.CoolProp as CP
# from CoolProp.Plots import PropertyPlot
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.optimize
# from scipy.interpolate import interp1d

# import os, sys
# sys.path.append(os.path.abspath(".."))

# # Internal Toolbox 
# from modules.f_lmtd2 import f_lmtd2
# from modules.propsfluid import propsfluid
# from modules.find_2P_boundaries import find_2P_boundaries

# # HTC and DP Correlations
# from modules.plate_htc import han_BPHEX_DP, water_plate_HTC, martin_BPHEX_HTC, muley_manglik_BPHEX_HTC, han_boiling_BPHEX_HTC, han_cond_BPHEX_HTC
# from modules.pipe_htc import gnielinski_pipe_htc, boiling_curve, horizontal_tube_internal_condensation
# from modules.shell_and_tube_htc import shell_bell_delaware_htc
# from modules.tube_bank_htc import ext_tube_film_condens
# from modules.fins import htc_tube_and_fins

# # Phase related correlations
# from modules.void_fraction import void_fraction
# from modules.kim_dry_out_incipience import kim_dry_out_incipience

# # Connectors
# from connector.mass_connector import MassConnector
# from connector.heat_connector import HeatConnector

# # Component base frame
# from component.base_component import BaseComponent

# #%%
# # Set to True to enable some debugging output to screen
# debug = False

# #%%

# class MB_HeatExchanger(BaseComponent):
       
#     class H():
#         pass
    
#     class C():
#         pass
    
#     def __init__(self, HTX_Type):
#         super().__init__()
        
#         self.su = {'H' : MassConnector(),
#                    'C' : MassConnector()}
        
#         self.ex = {'H' : MassConnector(),
#                    'C' : MassConnector()} # Mass_connector
        
#         self.Q_dot = HeatConnector()
        
#         if HTX_Type == 'Plate' or HTX_Type == 'Shell&Tube' or HTX_Type == 'Tube&Fins':
#             self.HTX_Type = HTX_Type
#         else:
#             print("Heat exchanger types implemented for this model are : 'Plate', 'Shell&Tube', 'Tube&Fins'.")

#     def get_required_inputs(self):
        
#         if self.inputs == {}:
#             # Hot Fluid
#             if self.su['H'].T is not None:
#                 self.inputs['Hsu_T'] = self.su['H'].T
#             elif self.su['H'].h is not None:
#                 self.inputs['Hsu_h'] = self.su['H'].h
#             if self.su['H'].p is not None:
#                 self.inputs['Hsu_p'] = self.su['H'].p
#             if self.su['H'].fluid is not None:
#                 self.inputs['Hsu_fluid'] = self.su['H'].fluid
#             if self.su['H'].m_dot is not None:
#                 self.inputs['Hsu_m_dot'] = self.su['H'].m_dot
                
#             # Cold Fluid                
#             if self.su['C'].T is not None:
#                 self.inputs['Csu_T'] = self.su['C'].T
#             elif self.su['C'].h is not None:
#                 self.inputs['Csu_h'] = self.su['C'].h
#             if self.su['C'].p is not None:
#                 self.inputs['Csu_p'] = self.su['C'].p
#             if self.su['C'].fluid is not None:
#                 self.inputs['Csu_fluid'] = self.su['C'].fluid
#             if self.su['C'].m_dot is not None:
#                 self.inputs['Csu_m_dot'] = self.su['C'].m_dot
                
#         if self.inputs != {}:
#             # Hot Fluid
#             self.su['H'].set_fluid(self.inputs['Hsu_fluid'])
#             if 'Hsu_T' in self.inputs:
#                 self.su['H'].set_T(self.inputs['Hsu_T'])
#             elif 'Hsu_h' in self.inputs:
#                 self.su['H'].set_h(self.inputs['Hsu_h'])
#             if 'Hsu_p' in self.inputs:
#                 self.su['H'].set_p(self.inputs['Hsu_p'])
#             if 'Hsu_m_dot' in self.inputs:
#                 self.su['H'].set_m_dot(self.inputs['Hsu_m_dot'])

#             # Cold Fluid
#             self.su['C'].set_fluid(self.inputs['Csu_fluid'])
#             if 'Csu_T' in self.inputs:
#                 self.su['C'].set_T(self.inputs['Csu_T'])
#             elif 'Csu_h' in self.inputs:
#                 self.su['C'].set_h(self.inputs['Csu_h'])
#             if 'Csu_p' in self.inputs:
#                 self.su['C'].set_p(self.inputs['Csu_p'])
#             if 'Csu_m_dot' in self.inputs:
#                 self.su['C'].set_m_dot(self.inputs['Csu_m_dot'])

#         return ['Hsu_p', 'Hsu_T', 'Hsu_m_dot', 'Hsu_fluid', 'Csu_p', 'Csu_T', 'Csu_m_dot', 'Csu_fluid']

#     def get_required_parameters(self):
        
#         general_parameters = ['Flow_Type','htc_type', 'H_DP_ON', 'C_DP_ON','n_disc']
        
#         if self.HTX_Type == 'Plate':     
#             geometry_parameters = ['A_c', 'A_h', 'h', 'l', 'l_v',
#                                    'C_CS', 'C_Dh', 'C_V_tot', 'C_canal_t', 'C_n_canals',
#                                    'H_CS', 'H_Dh', 'H_V_tot', 'H_canal_t', 'H_n_canals',
#                                    'casing_t', 'chevron_angle', 'fooling', 
#                                    'n_plates', 'plate_cond', 'plate_pitch_co', 't_plates', 'w']

#         elif self.HTX_Type == 'Shell&Tube':
                
#             geometry_parameters = ['A_eff', 'Baffle_cut', 'D_OTL', 'N_strips', 'S_V_tot',
#                                    'Shell_ID', 'T_V_tot', 'Tube_L', 'Tube_OD', 'Tube_pass',
#                                    'Tube_t', 'Tubesheet_t', 'central_spacing', 'clear_BS', 'clear_TB',
#                                    'cross_passes', 'foul_s', 'foul_t', 'inlet_spacing', 'n_series',
#                                    'n_tubes', 'outlet_spacing', 'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']
        
#         elif self.HTX_Type == 'Tube&Fins':
            
#             geometry_parameters = ['A_eff', 'Baffle_cut', 'D_OTL', 'N_strips', 'S_V_tot',
#                                    'Shell_ID', 'T_V_tot', 'Tube_L', 'Tube_OD', 'Tube_pass',
#                                    'Tube_t', 'Tubesheet_t', 'central_spacing', 'clear_BS', 'clear_TB',
#                                    'cross_passes', 'foul_s', 'foul_t', 'inlet_spacing', 'n_series',
#                                    'n_tubes', 'outlet_spacing', 'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']
        
#         return general_parameters + geometry_parameters
    
#     def print_setup(self):
#         print("=== Heat Exchanger Setup ===")
#         print("Connectors:")
#         print(f"  - H_su: fluid={self.su['H'].fluid}, T={self.su['H'].T}, p={self.su['H'].p}, m_dot={self.su['H'].m_dot}")
#         print(f"  - C_su: fluid={self.su['C'].fluid}, T={self.su['C'].T}, p={self.su['C'].p}, m_dot={self.su['C'].m_dot}")

#         print("\nInputs:")
#         for input in self.get_required_inputs():
#             if input in self.inputs:
#                 print(f"  - {input}: {self.inputs[input]}")
#             else:
#                 print(f"  - {input}: Not set")


#         print("\nParameters:")
#         for param in self.get_required_parameters():
#             if param in self.params:
#                 print(f"  - {param}: {self.params[param]}")
#             else:
#                 print(f"  - {param}: Not set")

#         print("======================")
           
# #%%
#     def set_HTC(self, htc_type = "Correlation", Corr_H = None, Corr_C = None, UD_H_HTC = None, UD_C_HTC = None):
        
#         self.params['htc_type'] = htc_type
        
#         self.check_calculable()
        
#         if htc_type == "User-Defined":
            
#             self.H.HeatExchange_Correlation = "User-Defined"
#             self.C.HeatExchange_Correlation = "User-Defined"
            
#             # !!! User-Defined Heat Transfer Coefficients (hot):
#             self.H.h_liq = UD_H_HTC['Liquid']
#             self.H.h_vap = UD_H_HTC['Vapor']
#             self.H.h_twophase = UD_H_HTC['Two-Phase']
#             self.H.h_vapwet = UD_H_HTC['Vapor-wet']
#             self.H.h_tpdryout = UD_H_HTC['Dryout']
#             self.H.h_transcrit = UD_H_HTC['Transcritical']
        
#             # !!! User-Defined Heat Transfer Coefficients (cold):
#             self.C.h_liq = UD_C_HTC['Liquid']
#             self.C.h_vap = UD_C_HTC['Vapor']
#             self.C.h_twophase = UD_C_HTC['Two-Phase']
#             self.C.h_vapwet = UD_C_HTC['Vapor-wet']
#             self.C.h_tpdryout = UD_C_HTC['Dryout']
#             self.C.h_transcrit = UD_C_HTC['Transcritical']
            
#             if self.params['H_DP_ON'] == True: # !!! arbitrarily set
#                 self.H.f_dp = {"K": 14.14, "B": 1.892}
#             else:
#                 self.H.f_dp = {"K": 0, "B": 0}
                    
#             if self.params['C_DP_ON'] == True:
#                 self.C.f_dp = {"K": 14.14, "B": 1.892}
#             else:
#                 self.C.f_dp = {"K": 0, "B": 0}
            
#         else:
#             # Type 
#             self.H.HeatExchange_Correlation = "Correlation"
#             self.C.HeatExchange_Correlation = "Correlation"
            
#             if self.HTX_Type == 'Plate' or self.HTX_Type == 'Shell&Tube' or self.HTX_Type == 'Tube&Fins':
                
#                 self.H.Correlation_1phase = Corr_H["1P"]
#                 self.H.Correlation_2phase = Corr_H["2P"]
                
#                 self.C.Correlation_1phase = Corr_C["1P"]
#                 self.C.Correlation_2phase = Corr_C["2P"]

#                 if self.params['H_DP_ON'] == True: # !!! arbitrarily set
#                     self.H.f_dp = {"K": 14.14, "B": 1.892}
#                 else:
#                     self.H.f_dp = {"K": 0, "B": 0}        
                    
#                 if self.params['C_DP_ON'] == True:
#                     self.C.f_dp = {"K": 14.14, "B": 1.892}
#                 else:
#                     self.C.f_dp = {"K": 0, "B": 0}
                    
#             if self.C.Correlation_2phase == "Boiling_curve":
#                 try:
#                     T_sat = CP.PropsSI('T','P', self.su['C'].p,'Q',0,self.su['C'].fluid)
#                     (h_boil, DT_vect) = boiling_curve(self.params['Tube_OD'], self.su['C'].fluid, T_sat, self.su['C'].p)
#                     self.C_f_boiling = interp1d(DT_vect,h_boil)
#                 except:
#                     self.C_f_boiling = interp1d([0,10000],[20000,20000])
            
# #%%
#     def external_pinching(self):
#         "Determine the maximum heat transfer rate based on the external pinching analysis"
        
#         "1) Set of arbitrary bound values" # !!! Find out why      
#         T_hmin = 218 # above CO2 freezing point (194.7 K)
            
#         T_cmax = 273.15+250
        
#         "2) Hot fluid side pinch"
        
#         # Computation of outlet lowest possible enthalpy of the hot fluid using either the entrance cold fluid inlet temperature, or the arbitrary minimal
#         self.h_ho = CP.PropsSI('H','T', max(self.T_ci, T_hmin),'P',self.p_hi,self.H_su.fluid) # Equation 5 (Bell et al. 2015)
        
#         # Q_max computation
#         Qmaxh = self.mdot_h*(self.h_hi-self.h_ho) # Equation 4 (Bell et al. 2015)
        
#         if debug:
#             print("Inlet Hot Enthalpy:", self.h_hi, "\n")
#             print("Outlet Hot Enthalpy:", self.h_ho,"\n")
#             print("mdot_h = ", self.mdot_h, "\n")
#             print("Qmaxh =", Qmaxh, "\n")

#         "3) Cold fluid side pinch"
                
#         # Computation of highest possible outlet enthalpy of the hot fluid using either the entrance hot fluid inlet temperature, or the arbitrary maximal
#         self.h_co = CP.PropsSI('H','T',min(self.T_hi, T_cmax),'P',self.p_ci,self.C_su.fluid) # Equation  (Bell et al. 2015)
        
#         # Q_max computation
#         Qmaxc = self.mdot_c*(self.h_co-self.h_ci) # Equation 6 (Bell et al. 2015)

#         if debug:
#             print("Inlet Cold Enthalpy:", self.h_ci, "\n")
#             print("Outlet Cold Enthalpy:", self.h_co,"\n")
#             print("mdot_c = ", self.mdot_c, "\n")
#             print("Qmaxh =", Qmaxc, "\n")
            
#         "4) Final pinch and cell boundaries computation"
        
#         Qmax = min(Qmaxh, Qmaxc)
        
#         if debug:
#             print('Qmax (external pinching) is', Qmax)

#         self.calculate_cell_boundaries(Qmax) # call calculate_cell_boundaries procedure

#         return Qmax
    
# #%%
#     def calculate_cell_boundaries(self, Q):
#         """ Calculate the cell boundaries for each fluid """
        
#         "1) Re-calculate the outlet enthalpies of each stream"
                
#         self.h_co = self.h_ci + Q/self.mdot_c
#         self.h_ho = self.h_hi - Q/self.mdot_h

#         "2) Calculate the dew and bubble pressures and enthalpies by accounting for pressure drops"
        
#         "2.1) For the cold fluid"
        
#         if not self.Transcritical_c: # if not in transcritical
#             if (not (self.C.f_dp == {"K": 0, "B": 0}) or self.DP_c < 1e-2):
#                 # If no pressure drop is calculated or if the pressure drop is neglectable, assign saturation conditions to the ideal case:
#                 self.p_cdew = self.p_ci
#                 self.p_cbubble = self.p_co
#                 self.T_cbubble = self.T_cbubble_ideal
#                 self.T_cdew = self.T_cdew_ideal
#                 self.h_cbubble = self.h_cbubble_ideal
#                 self.h_cdew = self.h_cdew_ideal
#             else:
#                 h_cbubble, h_cdew, p_cbubble, p_cdew, _, _, = find_2P_boundaries(self.C_su.fluid, self.h_ci, self.h_co, self.p_ci, self.p_co)
#                 self.p_cdew = p_cdew
#                 self.p_cbubble = p_cbubble
#                 self.h_cdew = h_cdew
#                 self.h_cbubble = h_cbubble
#                 self.T_cdew = CP.PropsSI("T", "P", self.p_cdew, "Q", 1, self.C_su.fluid)
#                 self.T_cbubble = CP.PropsSI("T", "P", self.p_cbubble, "Q", 0, self.C_su.fluid)
                
#         "2.2) For the hot fluid"
        
#         if not self.Transcritical_h:
#             if (not (self.H.f_dp == {"K": 0, "B": 0}) or self.DP_h < 1e-2) and not self.h_incomp_flag:
#                 #If no pressure drop is calculated or if the pressure drop is neglectable, assign saturation conditions to the ideal case:
#                 self.p_hdew = self.p_hi
#                 self.p_hbubble = self.p_ho    
#                 self.T_hbubble = self.T_hbubble_ideal
#                 self.T_hdew = self.T_hdew_ideal
#                 self.h_hbubble = self.h_hbubble_ideal
#                 self.h_hdew = self.h_hdew_ideal
#             elif not self.h_incomp_flag:
#                 h_hbubble, h_hdew, p_hdew, p_hbubble, _, _, = find_2P_boundaries(self.H_su.fluid, self.h_hi, self.h_ho, self.p_hi, self.p_ho)
#                 self.p_hdew = p_hdew
#                 self.p_hbubble = p_hbubble
#                 self.h_hdew = h_hdew
#                 self.h_hbubble = h_hbubble
#                 self.T_hdew = CP.PropsSI("T", "P", self.p_hdew, "Q", 1, self.H_su.fluid)
#                 self.T_hbubble = CP.PropsSI("T", "P", self.p_hbubble, "Q", 0, self.H_su.fluid)

#         "3) Discretization of the heat exchanger, user defined"
#         # The minimum number of cells desired is n_disc
#         # n_disc is increased by 1. Now it is the minimum number of cells boundaries.
#         n_disc = self.params['n_disc'] + 1
        
#         # Force the minimum number of discretizations needed
#         n_disc = max(n_disc, 2)
        
#         # The original Bell2015 code uses: "if h_cdew is not None". I don't see when that could happen without an error before.
#         "3.1) Create a vector of enthalpies : initiates the cell boundaries. If phase change, dew and bubble point are added"
        
#         # Cold side
#         if not self.Transcritical_c:
#             self.hvec_c = np.append(np.linspace(self.h_ci, self.h_co, n_disc), [self.h_cdew, self.h_cbubble])
#         elif self.Transcritical_c:
#             #In transcritical, just dont append the unexisting dew and bubble enthalpies
#             self.hvec_c = np.linspace(self.h_ci, self.h_co, n_disc)
            
#         self.hvec_c = np.sort(np.unique(self.hvec_c))
#         # Ensure that the enthalpy boundaries are all enclosed by the inlet and outlet conditions
#         self.hvec_c = [h for h in self.hvec_c if (h >= self.h_ci and h <= self.h_co)]
        
#         # Hot side
#         if not self.Transcritical_h and not self.h_incomp_flag:
#             self.hvec_h = np.append(np.linspace(self.h_hi, self.h_ho, n_disc), [self.h_hdew, self.h_hbubble])
#         elif self.Transcritical_h or self.h_incomp_flag:
#             #In transcritical, just dont append the unexisting dew and bubble enthalpies
#             self.hvec_h = np.linspace(self.h_hi, self.h_ho, n_disc)
            
#         self.hvec_h = np.sort(np.unique(self.hvec_h))
#         # Ensure that the enthalpy boundaries are all enclosed by the inlet and outlet conditions
#         self.hvec_h = [h for h in self.hvec_h if (h >= self.h_ho and h <= self.h_hi)]

#         "4) Filling of the complementary cell boundaries"
        
#         # Complementary cell boundaries serve to keep heat balances in check
#         # Start at the first element in the vector
#         k = 0
#         while k < len(self.hvec_c)-1 or k < len(self.hvec_h)-1:
#             if len(self.hvec_c) == 2 and len(self.hvec_h) == 2: # If there is only one cell
#                 break

#             # Determine which stream is the limiting next cell boundary
#             Qcell_hk = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])
#             Qcell_ck = self.mdot_c*(self.hvec_c[k+1]-self.hvec_c[k])

#             if debug:
#               print("k+1 Hot Enthalpy:", self.hvec_h[k+1], "\n")
#               print("k Hot Enthalpy:", self.hvec_h[k],"\n")
#               print("mdot_h = ", self.mdot_h, "\n")
#               print("Qcell_hk =", Qcell_hk, "\n-------\n")
#               print("k+1 Cold Enthalpy:", self.hvec_c[k+1], "\n")
#               print("k Cold Enthalpy:", self.hvec_c[k],"\n")
#               print("mdot_c = ", self.mdot_c, "\n")
#               print("Qcell_hk =", Qcell_ck, "\n-------\n")

#             if round(Qcell_hk, 4) > round(Qcell_ck, 4):
#                 # Hot stream needs a complementary cell boundary as the available heat is higher than the cold cell heat absorption
#                 self.hvec_h.insert(k+1, self.hvec_h[k] + Qcell_ck/self.mdot_h)
                
#             elif round(Qcell_hk, 4) < round(Qcell_ck, 4):
#                 # Cold stream needs a complementary cell boundary as the available heat absoprtion is higher than the hot cell heat 
#                 self.hvec_c.insert(k+1, self.hvec_c[k] + Qcell_hk/self.mdot_c)

#             if debug:
#                 print(k,len(self.hvec_c),len(self.hvec_h),Qcell_hk, Qcell_ck)

#             # Increment index
#             k += 1

#         if debug:
#              print("Modified Length of hvec_c is ", len(self.hvec_c))
#              print("Modified Length of hvec_h is ", len(self.hvec_h))

#         assert(len(self.hvec_h) == len(self.hvec_c)) # Ensures that the two fluid streams have the same number of cell boundaries

#         # Calculate the vectors of exchanged heat in each cell.
#         self.Qvec_h = np.array([self.mdot_h*(self.hvec_h[i+1]-self.hvec_h[i]) for i in range(len(self.hvec_h)-1)])
#         self.Qvec_c = np.array([self.mdot_c*(self.hvec_c[i+1]-self.hvec_c[i]) for i in range(len(self.hvec_c)-1)])

#         if debug:
#             if np.max(np.abs(self.Qvec_c/self.Qvec_h))<1e-5:
#                 print(self.Qvec_h, self.Qvec_c)
            
#         "5) Computation of the thermal states at the cell boundaries"
        
#         # Build the normalized enthalpy vectors (normalized by the total exchanged enthalpy)
#         self.hnorm_h = (np.array(self.hvec_h)-self.hvec_h[0])/(self.hvec_h[-1]-self.hvec_h[0])
#         self.hnorm_c = (np.array(self.hvec_c)-self.hvec_c[0])/(self.hvec_c[-1]-self.hvec_c[0])

#         # Pressure distribution accross the HX. Method by RDickes (not validated).
#         # Discretization of the pressure as a linear interpolation on enthalpy between in and out conditions # !!! Is this good ? 
#         self.pvec_c = (1-self.hnorm_c)*self.p_ci + self.hnorm_c*self.p_co
#         self.pvec_h = (1-self.hnorm_h)*self.p_hi + self.hnorm_h*self.p_ho

#         # Calculate the temperature and entropy at each cell boundary
#         self.Tvec_c = CP.PropsSI('T','H',self.hvec_c,'P',self.pvec_c,self.C_su.fluid)
        
#         self.Tvec_h = CP.PropsSI('T','H',self.hvec_h,'P',self.pvec_h,self.H_su.fluid)
#         self.svec_h = CP.PropsSI('S','H',self.hvec_h,'P',self.pvec_h,self.H_su.fluid)
        
#         self.svec_c = CP.PropsSI('S','H',self.hvec_c,'P',self.pvec_c,self.C_su.fluid)


#         "6) Vapour quality calculation if in the two-phase zone. If not, force, the value to -2 or 2"
#         # Take Transcritical into account. in that Case, force x = 3
        
#         "6.1) Hot Fluid"
#         self.x_vec_h = []
#         for i, h_h in enumerate(self.hvec_h):
#             if not self.Transcritical_h and not self.h_incomp_flag:
#                 if h_h < self.h_hbubble:
#                     self.x_vec_h.append(-2)        
#                 elif h_h > self.h_hdew:
#                     self.x_vec_h.append(2)
#                 else:
#                     self.x_vec_h.append(min(1, max(0, (h_h - self.h_hbubble)/(self.h_hdew - self.h_hbubble))))
#             else:
#                 self.x_vec_h.append(3)
                
#         "6.2) Cold Fluid"
#         self.x_vec_c = []
#         for i, h_c in enumerate(self.hvec_c):
#             if not self.Transcritical_c:
#                 if h_c < self.h_cbubble:
#                     self.x_vec_c.append(-2)        
#                 elif h_c > self.h_cdew:
#                     self.x_vec_c.append(2)
#                 else:
#                     self.x_vec_c.append(min(1, max(0, (h_c - self.h_cbubble)/(self.h_cdew - self.h_cbubble))))
#             else:
#                 self.x_vec_c.append(3)

#         "7) Vector of saturation temperatures at quality = 0.5, considering the pressure drop at each cell"
        
#         if not self.Transcritical_c:
#             self.Tvec_sat_pure_c = CP.PropsSI("T", "Q", 0.5, "P", self.pvec_c, self.C_su.fluid)
#         if not self.Transcritical_h and not self.h_incomp_flag:
#             self.Tvec_sat_pure_h = CP.PropsSI("T", "Q", 0.5, "P", self.pvec_h, self.H_su.fluid)
            
#         "8) Calculate pinch and heat exchanger border temperature deltas"
        
#         self.DT_pinch = np.min((np.array(self.Tvec_h) - np.array(self.Tvec_c)))
#         self.DT_ho_ci = self.Tvec_h[0] - self.Tvec_c[0]
#         self.DT_hi_co = self.Tvec_h[-1] - self.Tvec_c[-1]

# #%%
#     def internal_pinching(self, stream):
#         """
#         Determine the maximum heat transfer rate based on the internal pinching analysis
#         NB : external pinching analysis has already been done
#         """
        
#         # Try to find the dew point enthalpy as one of the cell boundaries
#         # that is not the inlet or outlet
#         # Do this only if the fluid is sub-critical

#         "1) Check for the hot stream"
        
#         if stream == 'hot':
#             if not self.Transcritical_h and not self.h_incomp_flag: # If the hot side is not transcritical
#                 for i in range(1,len(self.hvec_h)-1):
                    
#                     # Check if enthalpy is equal to the dewpoint enthalpy of hot stream and hot stream is colder than cold stream (impossible)
#                     if (abs(self.hvec_h[i] - self.h_hdew) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]):
                        
#                         # Enthalpy of the cold stream at the pinch temperature
#                         h_c_pinch = CP.PropsSI('H','T',self.T_hdew,'P',self.pvec_c[i], self.C_su.fluid) # Equation 10 (Bell et al. 2015)

#                         # Heat transfer in the cell
#                         Qright = self.mdot_h*(self.h_hi-self.h_hdew) # Equation 9 (Bell et al. 2015)
                        
#                         # New value for the limiting heat transfer rate
#                         Qmax = self.mdot_c*(h_c_pinch-self.h_ci) + Qright # Equation 12

#                         # Recalculate the cell boundaries
#                         self.calculate_cell_boundaries(Qmax)
                
#                         return Qmax
                    
#             elif self.Transcritical_h:
#                 #If the hot side is transcritical, do nothing
#                 pass
            
#             "2) Check for the cold stream"
            
#         elif stream == 'cold':
#             if not self.Transcritical_c:
#                 # Check for the cold stream pinch point
#                 for i in range(1,len(self.hvec_c)-1):
                    
#                     # Check if enthalpy is equal to the bubblepoint enthalpy of cold stream and hot stream is colder than cold stream (impossible)
#                     if (abs(self.hvec_c[i] - self.h_cbubble) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]):
                        
#                         # Enthalpy of the cold stream at the pinch temperatur
#                         h_h_pinch = CP.PropsSI('H','T',self.T_cbubble,'P',self.pvec_h[i], self.H_su.fluid) # Equation 14 (Bell et al. 2015)

#                         # Heat transfer in the cell
#                         Qleft = self.mdot_c*(self.h_cbubble-self.h_ci) # Equation 13 (Bell et al. 2015)

#                         # New value for the limiting heat transfer rate
#                         Qmax = Qleft + self.mdot_h*(self.h_hi-h_h_pinch) # Equation 16 (Bell et al. 2015)

#                         # Recalculate the cell boundaries
#                         self.calculate_cell_boundaries(Qmax)

#                         return Qmax
                    
#             elif self.Transcritical_c:
#                 #If the hot side is transcritical, do nothing
#                 pass
#         else:
#             raise ValueError
    
#     #%%
#     def solve(self, only_external = False, and_solve = True):
#         """
#         Parameters
#         ----------
#             - mdot_h : Hot fluid flowrate [kg/s]
#             - p_hi : Hot fluid pressure [Pa]
#             - h_hi : Hot fluid specific enthalpy [J/kg]
#             - mdot_c : Cold fluid flowrate [kg/s]
#             - p_ci : Cold fluid pressure [Pa]
#             - h_ci : Cold fluid specific enthalpy [J/kg]
#             - only_external (optional) : calls only exxternal_pinching function to determine Q_dot max if set to True (if there is no phase change)
#             - and_solve (optional) : if set to true, solves the heat exchanger

#         Returns
#         -------
#             - Q : Exchanged heat rate [W]

#         """
        
#         self.check_calculable()
#         self.check_parametrized()
        
#         if not self.calculable:
#             print("Component not calculable, check input")
            
#         if not self.parametrized:
#             print("Component not parametrized, check parameters")            
            
#         "1) Main Input variables"
        
#         self.H_su = self.su['H']
#         self.C_su = self.su['C']
            
#         self.H_ex = self.ex['H']
#         self.C_ex = self.ex['C']   
            
#         self.h_incomp_flag = (self.H_su.fluid.find("INCOMP") != -1)
                
#         # Hot fluid
#         self.mdot_h = self.H_su.m_dot
#         self.h_hi = self.H_su.h
#         self.p_hi = self.H_su.p
        
#         # Cold fluid 
#         self.mdot_c = self.C_su.m_dot
#         self.h_ci = self.C_su.h
#         self.p_ci = self.C_su.p
                
#         # Determine the inlet temperatures from the pressure/enthalpy pairs
#         self.T_ci = CP.PropsSI('T', 'P', self.p_ci, 'H', self.h_ci, self.C_su.fluid)
#         self.T_hi = CP.PropsSI('T', 'P', self.p_hi, 'H', self.h_hi, self.H_su.fluid)

#         "2) Determine if the streams come in at a higher pressure than the transcritical pressure"
        
#         # Initialization of the flags:    
#         self.Transcritical_c = False
#         self.Transcritical_h = False
        
#         if not self.h_incomp_flag and (self.p_hi - CP.PropsSI("PCRIT", self.H_su.fluid)) >= 1e-06:
#             self.Transcritical_h = True

#         "3) Calculate the ideal bubble and dew temperatures/enthalpies for each stream IF the fluid is not transcritical"
        
        
#         if not self.Transcritical_c:
#             self.T_cbubble_ideal = CP.PropsSI('T', 'P', self.p_ci, 'Q', 0, self.C_su.fluid)
#             self.T_cdew_ideal    = CP.PropsSI('T', 'P', self.p_ci, 'Q', 1, self.C_su.fluid)
#             self.h_cbubble_ideal = CP.PropsSI('H', 'T', self.T_cbubble_ideal, 'Q', 0, self.C_su.fluid)
#             self.h_cdew_ideal    = CP.PropsSI('H', 'T', self.T_cdew_ideal, 'Q', 1, self.C_su.fluid)  
            
#         if not self.Transcritical_h and not self.h_incomp_flag:
#             self.T_hbubble_ideal = CP.PropsSI('T', 'P', self.p_hi, 'Q', 0, self.H_su.fluid)
#             self.T_hdew_ideal    = CP.PropsSI('T', 'P', self.p_hi, 'Q', 1, self.H_su.fluid)
#             self.h_hbubble_ideal = CP.PropsSI('H', 'T', self.T_hbubble_ideal, 'Q', 0, self.H_su.fluid)
#             self.h_hdew_ideal    = CP.PropsSI('H', 'T', self.T_hdew_ideal, 'Q', 1, self.H_su.fluid)
            
#         "4) Calculate pressure drops"
        
#         if self.params['H_DP_ON'] == True: # if the pressure drop are not neglected
#             if self.H_su.fluid == 'Water':
#                 self.DP_h = self.H.f_dp["K"] * (self.mdot_h/self.H_su.D)**2
#                 self.p_ho = self.p_hi - self.DP_h
#             else:
#                 # rho_v = CP.PropsSI('D','Q',1,'P',self.p_hi,self.H_su.fluid)
#                 # rho_l = CP.PropsSI('D','Q',0,'P',self.p_hi,self.H_su.fluid)
#                 # mu_h = CP.PropsSI('V','H',self.h_hi,'P',self.p_hi,self.H_su.fluid)
                
#                 # G_h = (self.mdot_h/self.geom.H_n_canals)/self.geom.H_CS
#                 self.DP_h = 0 # self.H.f_dp["K"] * self.mdot_h**(self.H.f_dp["B"]) # Empirical correlations : DP = K*m_dot**B # Han_BPHEX_DP(mu_h, G_h, self.geom.H_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co, rho_v, rho_l, self.geom.l_v, self.geom.H_n_canals, self.mdot_h, self.geom.H_canal_t) # 
#                 self.p_ho = self.p_hi - self.DP_h
            
#         if self.C.f_dp != {"K": 0, "B": 0}:
#             if self.C_su.fluid == 'Water':
#                 self.DP_c = self.C.f_dp["K"] * (self.mdot_c/self.C_su.D)**2
#                 self.p_co = self.p_ci - self.DP_c
#             else:
#                 # rho_v = CP.PropsSI('D','Q',1,'P',self.p_ci,self.C_su.fluid)
#                 # rho_l = CP.PropsSI('D','Q',0,'P',self.p_ci,self.C_su.fluid)
#                 # mu_c = CP.PropsSI('V','H',self.h_ci,'P',self.p_ci,self.C_su.fluid)               
                
#                 # G_c = (self.mdot_c/self.geom.C_n_canals)/self.geom.C_CS
#                 self.DP_c = 0 # self.C.f_dp["K"] * self.mdot_c**(self.C.f_dp["B"]) # Han_BPHEX_DP(mu_c, G_c, self.geom.H_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co, rho_v, rho_l, self.geom.l_v, self.geom.C_n_canals, self.mdot_c, self.geom.C_canal_t) # 
#                 self.p_co = self.p_ci - self.DP_c
            
#         "5) Calculate maximum and actual heat rates"
        
#         if (self.T_hi - self.T_ci) > 1e-2  and self.mdot_h  > 0 and self.mdot_c > 0: # Check that the operating conditions allow for heat transfer
            
#             "5.1) Compute the external pinching & update cell boundaries"
#             Qmax_ext = self.external_pinching() # Call to external-pinching procedure
#             self.Qmax_ext = Qmax_ext
#             Qmax = Qmax_ext
            
#             if debug:
#                 print("External pinching calculation done. \n")
            
#             "5.2) Compute the internal pinching & update cell boundaries"
#             if not only_external: # If phase change is expected : Check the internal pinching
#                 for stream in ['hot','cold']:
#                     Qmax_int = self.internal_pinching(stream) # Call to internal-pinching procedure
#                     if Qmax_int is not None:
#                         self.Qmax_int = Qmax_int
#                         Qmax = Qmax_int
            
#             # Maximum heat transfer rate determined by external or internal pinching
#             self.Qmax = Qmax
            
#             "5.3) Solve the heat exchanger to find the actual heat rate"
#             if and_solve and not only_external:
#                 Q = self.solve_HX()

#             self.epsilon_th = self.Q/self.Qmax # HTX efficiency
#             self.residual = 1 - sum(self.w) # HTX residual # !!! (what is "w" ?)
            
#             "5.4) Effective density computation for each cell taking into account void fraction"
            
#             self.Dvec_h = np.empty(len(self.hvec_h))
#             self.Dvec_c = np.empty(len(self.hvec_c))
            
#             # Cross section void fraction : considered constant along each discretization in the void_fraction function
#             self.eps_void_h = np.empty(len(self.hvec_c)) 
#             self.eps_void_c = np.empty(len(self.hvec_c))
            
#             for i in range(len(self.hvec_h)):
#                 if self.x_vec_h[i] <= 0 or self.x_vec_h[i] >= 1: 
#                     self.Dvec_h[i] = CP.PropsSI('D', "P", self.pvec_h[i], "H", self.hvec_h[i], self.H_su.fluid)
#                     self.eps_void_h[i] = -1
#                 else:
#                     rho_g = CP.PropsSI('D', 'P', self.pvec_h[i], "Q", 1, self.H_su.fluid)
#                     rho_l = CP.PropsSI('D', 'P', self.pvec_h[i], "Q", 0, self.H_su.fluid)
#                     self.eps_void_h[i], self.Dvec_h[i] = void_fraction(self.x_vec_h[i], rho_g, rho_l)
                    
#             for i in range(len(self.hvec_c)-1):
#                 if self.x_vec_c[i] <= 0 or self.x_vec_c[i] >= 1: 
#                     self.Dvec_c[i] = CP.PropsSI('D', "P", self.pvec_c[i], "H", self.hvec_c[i], self.C_su.fluid)
#                     self.eps_void_c[i] = -1
#                 else:
#                     rho_g = CP.PropsSI('D', 'P', self.pvec_c[i], "Q", 1, self.C_su.fluid)
#                     rho_l = CP.PropsSI('D', 'P', self.pvec_c[i], "Q", 0, self.C_su.fluid)
#                     self.eps_void_c[i], self.Dvec_c[i] = void_fraction(self.x_vec_c[i], rho_g, rho_l)
            
#             "5.5) Computation of the fluid mass inside the HTX"
            
#             if self.HTX_Type == 'Plate':
#                 self.Vvec_h = self.params['H_V_tot']*np.array(self.w) # !!! Attention to this assumption
#                 self.Vvec_c = self.params['C_V_tot']*np.array(self.w)

#             elif self.HTX_Type == 'Shell&Tube':
#                 if self.params['Shell_Side'] == 'H': # Shell Side is the hot side
#                     self.Vvec_h = self.params['S_V_tot']*np.array(self.w) # !!! Attention to this assumption
#                     self.Vvec_c = self.params['T_V_tot']*np.array(self.w)
#                 else:
#                     self.Vvec_h = self.params['T_V_tot']*np.array(self.w) # !!! Attention to this assumption
#                     self.Vvec_c = self.params['S_V_tot']*np.array(self.w)                    

#             elif self.HTX_Type == 'Tube&Fins':
#                 if self.params['Fin_Side'] == 'H': # Shell Side is the hot side
#                     self.Vvec_h = self.params['B_V_tot']*np.array(self.w) # !!! Attention to this assumption
#                     self.Vvec_c = self.params['T_V_tot']*np.array(self.w)
#                 else:
#                     self.Vvec_h = self.params['T_V_tot']*np.array(self.w) # !!! Attention to this assumption
#                     self.Vvec_c = self.params['B_V_tot']*np.array(self.w) 
                
#             # Initiates the mass vectors # !!! (for each cell ?)
#             self.Mvec_h = np.empty(len(self.hvec_h)-1)
#             self.Mvec_c = np.empty(len(self.hvec_c)-1)
            
#             # Mass computation # Volume*Mean of the densities at the cell boundaries
#             for i in range(len(self.hvec_h)-1):
#                 self.Mvec_h[i] = self.Vvec_h[i]*(self.Dvec_h[i] + self.Dvec_h[i+1])/2
#                 self.Mvec_c[i] = self.Vvec_c[i]*(self.Dvec_c[i] + self.Dvec_c[i+1])/2

#             if and_solve:
                
#                 # Hot side
#                 self.ex['H'].set_fluid(self.H_su.fluid)
#                 self.ex['H'].set_T(self.Tvec_h[0])
#                 self.ex['H'].set_p(self.pvec_h[0])
#                 self.ex['H'].set_m_dot(self.H_su.m_dot)
                
#                 self.H_ex = self.ex['H']
                    
#                 # Cold side
#                 self.ex['C'].set_fluid(self.C_su.fluid)
#                 self.ex['C'].set_T(self.Tvec_c[-1]) # Example temperature [K]
#                 self.ex['C'].set_p(self.pvec_c[-1]) # Example Pressure [Pa]
#                 self.ex['C'].set_m_dot(self.C_su.m_dot) 
        
#                 self.C_ex = self.ex['C']
                                            
#                 self.defined = True
#                 return Q
            
#         else: # Just a flag if the heat exchanger is not solved
#             self.Q = 1
            
# #%%    
#     def objective_function(self, Q):
        
#         "1) Initialize cell boundaries and results vectors"
        
#         # Cell boundaries
#         self.calculate_cell_boundaries(Q)

#         # Initialize dry-out incipience identification
#         self.x_di_c = 1
#         self.dry_out_c = False
#         self.ColdSide_Props_Calculated = False
        
#         # Initialise results arrays
#         w = [] # The wet-bulb potential temperature
#         self.Avec_h = []
#         self.alpha_c = []
#         self.alpha_h = []
#         self.phases_h = []
#         self.phases_c = []

#         "2) Iteration over all cells"
        
#         # The length of the hvec vectors is the same. hvec_h is arbitrarily taken below
#         for k in range(len(self.hvec_h)-1): # iterate over each cell
            
#             "2.1) Diverse temperature parameter definition"
            
#             Thi = self.Tvec_h[k+1] 
#             Tci = self.Tvec_c[k]
#             Tho = self.Tvec_h[k]
#             Tco = self.Tvec_c[k+1]
            
#             DTA = max(Thi - Tco, 1e-02)
#             DTB = max(Tho - Tci, 1e-02)
            
#             if DTA == DTB:
#                 LMTD = DTA
#             else:
#                 try:
#                     LMTD = (DTA-DTB)/np.log(abs(DTA/DTB))
#                 except ValueError as VE:
#                     print(Q, DTA, DTB)
#                     raise
            
#             # Wall Temperature:
#             T_wall = (Thi + Tho + Tci + Tco)/4
            
#             # Cold side
#             Tc_mean = 0.5*(Tci + Tco) # mean temperature over the cell
            
#             if not self.Transcritical_c:
#                 Tc_sat_mean = 0.5*(self.Tvec_sat_pure_c[k] + self.Tvec_sat_pure_c[k+1]) # mean saturation temperature over the cell
                
#             p_c_mean = 0.5*(self.pvec_c[k] + self.pvec_c[k+1]) # mean pressure over the cell
            
#             # Hot side
#             Th_mean = 0.5*(Thi + Tho) # mean temperature over the cell
            
#             if not self.Transcritical_h and not self.h_incomp_flag:
#                 Th_sat_mean = 0.5*(self.Tvec_sat_pure_h[k] + self.Tvec_sat_pure_h[k+1]) # mean saturation temperature over the cell
            
#             p_h_mean = 0.5*(self.pvec_h[k] + self.pvec_h[k+1]) # mean pressure over the cell

#             "2.2) Hot side phase identification"
            
#             # If not transcritical
#             if not self.Transcritical_h and not self.h_incomp_flag:
                
#                 havg_h = (self.hvec_h[k] + self.hvec_h[k+1])/2.0 # Average cell enthalpy over the cell
                
#                 if havg_h < self.h_hbubble: # below evaporation/bubble point
#                     self.phases_h.append('liquid')
#                     T_wall_h = T_wall
                    
#                 elif havg_h > self.h_hdew: # higher than condensation/dew point
#                     #T_wall_h not used so far but it could appear in future Correlations implemented.
#                     T_wall_h = max(T_wall, Th_sat_mean)
                    
#                     if T_wall > Th_sat_mean:
#                         self.phases_h.append('vapor')
#                     else:
#                         self.phases_h.append('vapor-wet')
                        
#                 else: # Between evaporation/bubble and condensation/dew point
#                     self.phases_h.append('two-phase')
#                     T_wall_h = T_wall
                    
#             elif self.h_incomp_flag:
#                 self.phases_h.append('liquid')
#                 T_wall_h = T_wall
                
#             elif self.Transcritical_h:
#                 self.phases_h.append("transcritical")
#                 T_wall_h = T_wall
            
#             "2.3) Cold side phase identification"

#             # If not transcritical
#             if not self.Transcritical_c:
                
#                 havg_c = (self.hvec_c[k] + self.hvec_c[k+1])/2.0 # Average cell enthalpy over the cell
                
#                 if havg_c < self.h_cbubble: # below evaporation/bubble point
#                     self.phases_c.append('liquid')
#                     #The line comented below is the original line on the code. It is causing trouble.
#                     #T_wall_c = min(T_wall, Tc_sat_mean)
#                     T_wall_c = T_wall
                    
#                 elif havg_c > self.h_cdew: # higher than condensation/dew point
#                     self.phases_c.append("vapor")
#                     T_wall_c = T_wall
                    
#                 else: # Between evaporation/bubble and condensation/dew point
#                     T_wall_c = T_wall
                    
#                     if (0.5*self.x_vec_c[k] + 0.5*self.x_vec_c[k+1]) <= self.x_di_c:
#                         self.phases_c.append('two-phase')  
                        
#                     else:
#                         self.phases_c.append("two-phase-dryout")
#                         self.dry_out_c = True
                        
#             elif self.Transcritical_c:
#                 self.phase_c.append('transcritical')
            
#             if self.HTX_Type == 'Plate':                
#                 G_c = (self.mdot_c/self.params['C_n_canals'])/self.params['C_CS']
#                 G_h = (self.mdot_h/self.params['H_n_canals'])/self.params['H_CS']
                
#             elif self.HTX_Type == 'Shell&Tube':
#                 A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            
#                 # Mass flow in tubes if the fluid is in the tubes
#                 G_h = (self.mdot_h/self.params['n_tubes'])/A_in_one_tube
#                 G_c = (self.mdot_c/self.params['n_tubes'])/A_in_one_tube

#             elif self.HTX_Type == 'Tube&Fins':
#                 A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            
#                 # Mass flow in tubes if the fluid is in the tubes
#                 G_h = (self.mdot_h/self.params['n_tubes'])/A_in_one_tube
#                 G_c = (self.mdot_c/self.params['n_tubes'])/A_in_one_tube
                
#             # F correction factor for LMTD method:
#             if self.params['Flow_Type'] != "CounterFlow":
                
#                 try:
#                     C_c = CP.PropsSI('C', 'T', Tc_mean,'P',p_c_mean,self.C_su.fluid)
#                 except:
#                     C_c = 20000
                    
#                 try:
#                     C_h = CP.PropsSI('C', 'T', Th_mean,'P',p_h_mean,self.H_su.fluid)
#                 except:
#                     C_h = 20000
                
#                 C_min = min(C_c,C_h)
#                 C_max = max(C_c,C_h)
                
#                 C_r = C_min/C_max
                
#                 self.R = (Thi - Tho)/(Tco - Tci)
#                 self.P = (Tco - Tci)/(Thi - Tci)
#                 self.F = f_lmtd2(self.R, self.P, self.params, C_r)
#             else:
#                 self.F = 1
            
#             #UA_req including F correction factor in case of cross flow
#             UA_req = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])/(self.F*LMTD)          
            
#             "3) Cell heat transfer coefficients"
            
#             "3.1) Hot side - User defined"
            
#             if self.H.HeatExchange_Correlation == "User-Defined":
#                 if self.phases_h[k] == "liquid":
#                     alpha_h = self.H.h_liq
#                 elif self.phases_h[k] == "vapor":
#                     alpha_h = self.H.h_vap
#                 elif self.phases_h[k] == "two-phase":
#                     alpha_h = self.H.h_twophase
#                 elif self.phases_h[k] == "vapor-wet":
#                     alpha_h = self.H.h_vapwet
#                 elif self.phases_h[k] == "two-phase-dryout":
#                     alpha_h = self.H.h_tpdryout
#                 elif self.phases_h[k] == "transcritical":
#                     alpha_h = self.H.h_transcrit
                    
#             elif self.H.HeatExchange_Correlation == "Correlation": 
                
#                 # Heat transfer coefficient calculated from Correlations.
#                 # Evaluate phases
                
#                 # 1 phase case
#                 if self.phases_h[k] == "liquid" or self.phases_h[k] == "vapor" or self.phases_h[k] == "transcritical":
                                        
#                     if self.H.Correlation_1phase == "Gnielinski":
#                         try:
#                             mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = propsfluid(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, False)
#                             Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', Th_mean,'P',p_h_mean,self.H_su.fluid)
#                         except (ValueError):
#                             if self.phases_h[k] == "liquid":
#                                 mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = propsfluid(Th_mean-0.1, p_h_mean, T_wall_h, self.H_su.fluid, False)
#                                 Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', Th_mean-1,'P',p_h_mean,self.H_su.fluid)
#                             elif self.phases_h[k] == "vapor":
#                                 mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = propsfluid(Th_mean+0.1, p_h_mean, T_wall_h, self.H_su.fluid, False)
#                                 Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', Th_mean+1,'P',p_h_mean,self.H_su.fluid)

#                         if self.HTX_Type == 'Plate':
#                             alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_h, mu_h_w, Pr_h, k_h, G_h, self.geom.H_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_h, Pr_h, k_h, G_h, self.geom.H_Dh) # 
                        
#                         elif self.HTX_Type == 'Shell&Tube':
#                             alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                        
#                         elif self.HTX_Type == 'Tube&Fins':
#                             alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['n_passes']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                    
#                     elif self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC":
#                         alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean, T_wall, p_h_mean, self.H_su.fluid, self.params)
                        
#                     elif self.H.Correlation_1phase == 'Tube_And_Fins':
#                         alpha_h = htc_tube_and_fins(self.H_su.fluid, self.params, p_h_mean, Th_mean, self.mdot_h, self.params['Fin_type'])[0]
                        
#                 # 2 phases case
#                 elif self.phases_h[k] == "two-phase" or self.phases_h[k] == "vapor-wet":
#                     if self.phases_h[k] == "two-phase":
#                         x_h = min(1, max(0, 0.5*(self.x_vec_h[k] + self.x_vec_h[k])))
#                     elif self.phases_h[k] == "vapor-wet":
#                         x_h = 1
                    
#                     # Thermodynamical variables
#                     mu_h_l = CP.PropsSI("V", "Q", 0, "P", p_h_mean, self.H_su.fluid)
#                     k_h_l = CP.PropsSI("L", "Q", 0, "P", p_h_mean, self.H_su.fluid)
#                     Pr_h_l = CP.PropsSI("Prandtl", "Q", 0, "P", p_h_mean, self.H_su.fluid)
#                     rho_h_l = CP.PropsSI("D", "Q", 0, "P", p_h_mean, self.H_su.fluid)
#                     rho_h_v = CP.PropsSI("D", "Q", 1, "P", p_h_mean, self.H_su.fluid)
#                     # i_fg_h = CP.PropsSI('H', 'Q', 1, 'P', p_h_mean, self.H_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_h_mean, self.H_su.fluid)
                    
#                     # !!! Include different types of Correlation HERE
#                     if self.H.Correlation_2phase == "Han_cond_BPHEX":
#                         alpha_h_2phase, _, DP_H = han_cond_BPHEX_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, rho_h_l, rho_h_v, G_h, self.params['H_Dh'], self.params['plate_pitch_co'], self.params['chevron_angle'], self.params['l_v'], self.params['H_n_canals'], self.H_su.m_dot, self.params['H_canal_t'])
                    
#                     if self.H.Correlation_2phase == 'ext_tube_film_condens':
#                         V_flow = G_h/CP.PropsSI('D','H',havg_h,'P',p_h_mean,self.H_su.fluid)
                        
#                         if self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC":
#                             alpha_h = shell_bell_delaware_htc(self.mdot_h, Th_mean, T_wall, p_h_mean, self.H_su.fluid, self.params)
                            
#                         elif self.H.Correlation_1phase == 'Tube_And_Fins':
#                             alpha_h = htc_tube_and_fins(self.C_su.fluid, self.params, p_c_mean, Tc_mean, self.mdot_c, self.params['Fin_type'])[0]
                            
#                         alpha_h_2phase = ext_tube_film_condens(self.params['Tube_OD'], self.H_su.fluid, Th_mean, T_wall, V_flow)
                        
#                     if self.H.Correlation_2phase == 'Horizontal_Tube_Internal_Condensation':
#                         mu_h, Pr_h, k_h = CP.PropsSI(('V','PRANDTL','L'), 'H', havg_h,'P',p_h_mean,self.H_su.fluid)
#                         Pr_h_w, mu_h_w = CP.PropsSI(('PRANDTL','V'), 'T', T_wall,'P',p_h_mean,self.H_su.fluid)
                                
#                         alpha_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h_w, k_h, G_h, self.params['Tube_OD'] - 2*self.params['Tube_t'], self.params['Tube_L'])
#                         alpha_h_2phase = horizontal_tube_internal_condensation(self.H_su.fluid , G_c, p_h_mean, self.x_vec_h[k], T_wall, self.params['Tube_OD'] - 2*self.params['Tube_t'])
                                        
#                     if self.phases_h[k] == "two-phase":
#                         alpha_h = alpha_h_2phase
                        
#                     elif self.phases_h[k] == "vapor-wet":
#                         w_vap_wet = (Th_mean - Th_sat_mean)/(Th_mean - T_wall)
#                         # The line below re-calculates alpha_h in case of having a vapor-wet condition
#                         alpha_h = alpha_h_2phase - w_vap_wet*(alpha_h_2phase - alpha_h) # the last alpha_h in this equation is the 1 Phase calculation

#             "3.2) Cold side - User defined"
            
#             # User Defined
#             if self.C.HeatExchange_Correlation == "User-Defined":
#                 if self.phases_c[k] == "liquid":
#                     alpha_c = self.C.h_liq
#                 elif self.phases_c[k] == "vapor":
#                     alpha_c = self.C.h_vap
#                 elif self.phases_c[k] == "two-phase":
#                     alpha_c = self.C.h_twophase
#                 elif self.phases_c[k] == "vapor-wet":
#                     alpha_c = self.C.h_vapwet
#                 elif self.phases_c[k] == "two-phase-dryout":
#                     alpha_c = self.C.h_tpdryout
#                 elif self.phases_c[k] == "transcritical":
#                     alpha_c = self.C.h_transcrit
                    
#             elif self.C.HeatExchange_Correlation == "Correlation":

#                 # Heat transfer coefficient calculated from Correlations:
#                 # Evaluate phases
                
#                 # 1 phase case
#                 if self.phases_c[k] == "liquid" or self.phases_c[k] == "vapor" or self.phases_c[k] == "transcritical":
#                     mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False)
            
#                     if self.HTX_Type == 'Plate' and self.C_su.fluid == 'water':
#                         alpha_c = water_plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.H_Dh)

#                     elif self.C.Correlation_1phase  == "Gnielinski":
#                         try:
#                             mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False)
#                             Pr_c_w, mu_c_w = CP.PropsSI(('PRANDTL','V'), 'T', Tc_mean,'P',p_c_mean,self.C_su.fluid)
#                         except (ValueError):
#                             if self.phases_c[k] == "liquid":
#                                 mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean-0.1, p_c_mean, T_wall_c, self.C_su.fluid, False)
#                                 Pr_c_w, mu_c_w = CP.PropsSI(('PRANDTL','V'), 'T', Tc_mean-1,'P',p_c_mean,self.C_su.fluid)
#                             elif self.phases_c[k] == "vapor":
#                                 mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = propsfluid(Tc_mean+0.1, p_c_mean, T_wall_c, self.C_su.fluid, False)
#                                 Pr_c_w, mu_c_w = CP.PropsSI(('PRANDTL','V'), 'T', Tc_mean+1,'P',p_c_mean,self.C_su.fluid)
                        
#                         if self.HTX_Type == 'Plate':
#                             alpha_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c_w, k_c, G_c, self.params['H_Dh'], self.params['l']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                        
#                         elif self.HTX_Type == 'Shell&Tube':
#                             alpha_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                        
#                         elif self.HTX_Type == 'Tube&Fins':
#                             alpha_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c_w, k_c, G_c, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['n_passes']) # Muley_Manglik_BPHEX_HTC(mu_c, mu_c_w, Pr_c, k_c, G_c, self.geom.C_Dh, self.geom.chevron_angle) # Simple_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.C_Dh) # 
                            
#                     elif self.C.Correlation_1phase == 'Shell_Bell_Delaware_HTC':
#                         alpha_c = shell_bell_delaware_htc(self.mdot_c, Tc_mean, T_wall, p_c_mean, self.C_su.fluid, self.params)
                         
#                     elif self.C.Correlation_1phase == 'Tube_And_Fins':
#                         alpha_c = htc_tube_and_fins(self.C_su.fluid, self.params, p_c_mean, Tc_mean, self.mdot_c, self.params['Fin_type'])[0]
                    
#                     # if self.C.Correlation_1phase  == "Gnielinski":
#                     #     alpha_c, _ = Gnielinski_Pipe_HTC(mu_c, Pr_c, k_c, G_c, self.geom.H_Dh, self.geom.l)

#                 # 2 phases case
#                 elif self.phases_c[k] == "two-phase" or self.phases_c[k] == "vapor-wet":
#                     if self.phases_c[k] == "two-phase":
#                         x_c = min(1, max(0, 0.5*(self.x_vec_c[k] + self.x_vec_c[k])))
#                     elif self.phases_c[k] == "vapor-wet":
#                         x_c = 1
                        
#                     # Thermodynamical variables
#                     mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
#                     k_c_l = CP.PropsSI("L", "Q", 0, "P", p_c_mean, self.C_su.fluid)
#                     Pr_c_l = CP.PropsSI("Prandtl", "Q", 0, "P", p_c_mean, self.C_su.fluid)
#                     rho_c_l = CP.PropsSI("D", "Q", 0, "P", p_c_mean, self.C_su.fluid)
#                     rho_c_v = CP.PropsSI("D", "Q", 1, "P", p_c_mean, self.C_su.fluid)
#                     i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                    
#                     # This boolean serves to spare to recalculate this properties in the dry-out analysis further ahead.
#                     self.ColdSide_Props_Calculated = True
                    
#                     # !!! Include different types of Correlation HERE
#                     if self.C.Correlation_2phase == "Han_cond_BPHEX":
#                         alpha_c_2phase, _ = han_cond_BPHEX_HTC(x_c, mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v, G_c, self.params['C_Dh'], self.params['plate_pitch_co'], self.params['chevron_angle'])
                    
#                     elif self.C.Correlation_2phase == "Han_Boiling_BPHEX_HTC":
#                         alpha_c_2phase, _ = han_boiling_BPHEX_HTC(min(x_c, self.x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, LMTD*self.F, self.Qvec_c[k], alpha_h, self.params['C_Dh'], self.params['chevron_angle'], self.params['plate_pitch_co'])
                    
#                     elif self.C.Correlation_2phase == "Boiling_curve":                  
#                         alpha_c_2phase = self.C_f_boiling(abs(T_wall - Tc_mean))
                    
#                     else:
#                         raise ValueError("Correlation not found for Cold Side 2-Phase")
                            
#                     if self.phases_c[k] == "two-phase":
#                         alpha_c = alpha_c_2phase
#                     elif self.phases_c[k] == "vapor-wet":
#                         w_vap_wet = (Tc_mean - Tc_sat_mean)/(Tc_mean - T_wall)
#                         #The line below re-calculates alpha_h in case of having a vapor-wet condition
#                         alpha_c = alpha_c_2phase - w_vap_wet*(alpha_c_2phase - alpha_c) #the last alpha_c in this equation is the 1 Phase calculation
            
#             "3.3) Store the heat transfer coefficients"
            
#             self.alpha_c.append(alpha_c)
#             self.alpha_h.append(alpha_h)
            
#             "4) Recover some parameters and conductivity calculation and compute UA_avail"
            
#             if self.HTX_Type == 'Plate':     
#                 # In the equation below, thickness resistance is given with respect to A_h arbitrarely
#                 UA_avail = 1/((1+self.params['fooling'])/(alpha_h*self.params['A_h']) + 1/(alpha_c*self.params['A_c']) + self.params['t_plates']/(self.params['plate_cond'])) 
                
#             elif self.HTX_Type == 'Shell&Tube':       
                
#                 fact_cond_1 = np.log(self.params['Tube_OD']/(self.params['Tube_OD'] - 2*self.params['Tube_t']))
#                 fact_cond_2 = 2*np.pi*self.params['tube_cond']*self.params['Tube_L']*self.params['n_series']*self.params['n_tubes']
#                 R_cond = fact_cond_1/fact_cond_2                      

#                 self.A_in_tubes = self.params['Tube_pass']*self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t']))
#                 self.A_out_tubes = self.params['Tube_pass']*self.params['Tube_L']*self.params['n_tubes']*np.pi*(self.params['Tube_OD'])
            
#                 if self.params['foul_s'] != None:
#                     R_fouling_s = self.params['foul_s'] / self.A_out_tubes
#                 else: 
#                     R_fouling_s = 0
    
#                 if self.params['foul_t'] != None:
#                     R_fouling_t = self.params['foul_t'] / self.A_in_tubes
#                 else: 
#                     R_fouling_t = 0                
            
#                 UA_avail = 1/(1/(alpha_c*self.A_in_tubes) + 1/(alpha_h*self.A_out_tubes) + R_fouling_s + R_fouling_t + R_cond) 
            
#             elif self.HTX_Type == 'Tube&Fins':
                
#                 fact_cond_1 = np.log(self.params['Tube_OD']/(self.params['Tube_OD'] - 2*self.params['Tube_t']))
#                 fact_cond_2 = 2*np.pi*self.params['Tube_cond']*self.params['Tube_L']*self.params['n_tubes']*self.params['n_passes']
#                 R_cond = fact_cond_1/fact_cond_2                      
                
#                 self.A_in_tubes = self.params['n_passes']*self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t'])) # *self.geom.n_passes 
#                 self.A_out_tubes = self.params['A_finned']
                            
#                 #☺ self.A_out_tubes = 8705 #self.geom.n_passes*self.geom.Tube_L*self.geom.n_tubes*np.pi*(self.geom.Tube_OD)
            
#                 R_fouling = 0 # self.geom.fouling / self.A_out_tubes             
                        
#                 # In the equation below, thickness resistance is given with respect to A_h arbitrarely
            
#                 if self.params['Fin_Side'] == 'H': # Fin side is the hot side 
#                     UA_avail =  1/(1/(alpha_c*self.A_in_tubes) + 1/(alpha_h*self.A_out_tubes) + R_fouling + R_cond) # 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.tube_cond)) 
#                 else: 
#                     UA_avail =  1/(1/(alpha_h*self.A_in_tubes) + 1/(alpha_c*self.A_out_tubes) + R_fouling + R_cond) # 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.tube_cond)) 
      
#             # R_fooling_h = self.params.H.R_fooling 
#             # R_fooling_c = self.params.C.R_fooling
            
#             "5) Calculate w, the main objective of this function. This variable serves for residual minimization in the solver"
            
#             w.append(UA_req/UA_avail)
            
#             self.w = w
            
#             if self.HTX_Type == 'Plate':
#                 self.Avec_h.append(w[k]*self.params['A_h']) 
                
#             if self.HTX_Type == 'Shell&Tube':
#                 self.Avec_h.append(w[k]*self.params['A_eff'])

#             if self.HTX_Type == 'Tube&Fins':
#                 self.Avec_h.append(w[k]*self.params['A_finned'])
            
#             "5.1) Calculate dryout incipience only if subcritical"
#             if not self.dry_out_c and self.phases_h[k] == "two-phase" and not self.Transcritical_c:
#                     if not self.ColdSide_Props_Calculated:
#                         #If these properties where not calculated before, do it:
#                         mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
#                         rho_c_l = CP.PropsSI('D', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
#                         rho_c_v = CP.PropsSI('D', 'Q', 1, 'P', p_c_mean, self.H_su.fluid)
#                         i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                    
#                     P_star_c = p_c_mean/CP.PropsSI("Pcrit", "Q", 1, "P", p_c_mean, self.C_su.fluid)
#                     q_c = self.Qvec_c[k]/self.Avec_h[k]
#                     try:
#                         sigma_c_l = CP.PropsSI("I", 'Q', 0, "P", p_c_mean, self.C_su.fluid)
                        
#                         if self.HTX_Type == 'Plate':
#                             self.x_di_c = kim_dry_out_incipience(G_c, q_c, self.params['C_Dh'], P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
                            
#                         elif self.HTX_Type == 'Shell&Tube':
#                             self.x_di_c = kim_dry_out_incipience(G_c, q_c,self.params['Tube_OD']-2*self.params['Tube_t'], P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
#                     except:
#                         self.x_di_c = -1
#         if debug:
#             print(Q, 1-sum(w))
            
#         return 1-sum(w)

# #%%
#     def solve_HX(self):
#         """ 
#         Solve the objective function using Brent's method and the maximum heat transfer 
#         rate calculated from the pinching analysis
#         """
#         self.Q = scipy.optimize.brentq(self.objective_function, 1e-5, self.Qmax-1e-10, rtol = 1e-14, xtol = 1e-10)
#         return self.Q
    
# #%%
#     def plot_objective_function(self, N = 100):
#         """ Plot the objective function """
#         Q = np.linspace(1e-5,self.Qmax,N)
#         r = np.array([self.objective_function(_Q) for _Q in Q])
#         fig, ax = plt.subplots()
#         ax.plot(Q, r)
#         ax.grid(True)
#         # plt.show()²
        
# #%%
#     def plot_ph_pair(self):
#         import warnings
#         warnings.simplefilter("ignore")
#         """ Plot p-h plots for the pair of working fluids """
#         Ph_diagram_h = PropertyPlot(self.H_su.fluid, "PH", unit_system = "EUR")
#         plt.plot(0.001*np.array(self.hvec_h), self.pvec_h*1e-05,'s-')
#         Ph_diagram_h.calc_isolines()
#         Ph_diagram_h.title("Hot side P-h diagram. Fluid: " + str(self.H_su.fluid))
#         Ph_diagram_h.show()
#         #----------------------------------------------------------------------
#         Ph_diagram_c = PropertyPlot(self.C_su.fluid, "PH", unit_system = "EUR")
#         plt.plot(0.001*np.array(self.hvec_c), self.pvec_c*1e-05,'s-')
#         Ph_diagram_c.calc_isolines()
#         Ph_diagram_c.title("Cold side P-h diagram. Fluid: " + str(self.C_su.fluid))
#         Ph_diagram_c.show()
#         warnings.simplefilter("default")

# #%%
#     def plot_Ts_pair(self):
#         import warnings
#         warnings.simplefilter("ignore")
#         """ Plot a T-s plot for the pair of working fluids """
#         Ph_diagram_h = PropertyPlot(self.H_su.fluid, "TS", unit_system = "EUR")
#         plt.plot(0.001*self.svec_h,self.Tvec_h - 273.15,'s-')
#         Ph_diagram_h.calc_isolines()
#         Ph_diagram_h.title("Hot side T-s diagram. Fluid: " + str(self.H_su.fluid))
#         Ph_diagram_h.show()
#         # #----------------------------------------------------------------------
#         Ph_diagram_c = PropertyPlot(self.C_su.fluid, "TS", unit_system = "EUR")
#         plt.plot(0.001*self.svec_c,self.Tvec_c - 273.15,'s-')
#         Ph_diagram_c.calc_isolines()
#         Ph_diagram_c.title("Cold side T-s diagram. Fluid: " + str(self.C_su.fluid))
#         Ph_diagram_c.show()
#         warnings.simplefilter("default")

# #%%
#     def plot_cells(self, fName = '', dpi = 400):
#         """ Plot the cells of the heat exchanger """
#         plt.figure(figsize = (4,3))
#         plt.plot(self.hnorm_h, self.Tvec_h, 'rs-')
#         plt.plot(self.hnorm_c, self.Tvec_c, 'bs-')
#         plt.xlim(0,1)
#         plt.ylabel('T [K]') 
#         plt.xlabel(r'$\hat h$ [-]')
#         plt.grid(True)
#         plt.tight_layout(pad = 0.2) 
#         plt.show()
#         if fName != '':
#             plt.savefig(fName, dpi = dpi)
