import __init__

"""
import sys
import os

# Get the absolute path of the directory that contains the script (simulation_model.py)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Determine the project root directory (which contains both 'connector' and 'component')
project_root = os.path.abspath(os.path.join(current_dir, '..', '..')) 

# Add the project root to sys.path if it's not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)
"""

#%%

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

from CoolProp.CoolProp import AbstractState
import CoolProp.CoolProp as CoolProp
from scipy.optimize import fsolve, root, minimize
import numpy as np
import math

class HexCstPinch(BaseComponent):
    """
    Component: Heat Exchanger with constant pinch point.

    **Description**:
    
        Simulates a Heat Exchanger with a constant pinch point.
        The Pinch Point is understood as being the location on the heat exchanger length where
        the temperature difference between the hot and the cold fluid is minimal.

    **Assumptions**:

        - Steady-state operation
        - No pressure drops considered
        - No loss to the ambient considered.

    **Connectors**:

        su_H (MassConnector): Mass connector for the hot suction side.
        su_C (MassConnector): Mass connector for the cold suction side.

        ex_H (MassConnector): Mass connector for the hot exhaust side.
        ex_C (MassConnector): Mass connector for the cold exhaust side.

        Q (HeatConnector): Heat connector for the heat transfer between the fluids

    **Parameters**:

        Pinch: Pinch point temperature difference [K] or [°C]
        
        Delta_T_sh_sc: Superheating or subcooling, depending if the HEX is an evaporator (superheating) or a condenser (subcooling)
            
        HX_type: HX type, i.e. evaporator or condenser

    **Inputs**:
        
        fluid_C: Cold suction side fluid. [-]

        h_su_C: Cold suction side enthalpy. [J/kg]

        P_su_C: Cold suction side pressure. [Pa]

        m_dot_C: Cold suction side mass flow rate. [kg/s]

        fluid_H: Hot suction side fluid. [-]

        h_su_H: Hot suction side enthalpy. [J/kg]

        P_su_H: Hot suction side pressure. [Pa]

        m_dot_H: Hot suction side mass flow rate. [kg/s]

    **Outputs**:

        h_ex_C: Cold exhaust side enthalpy. [J/kg]

        P_ex_C: Cold exhaust side pressure. [Pa]

        h_ex_H: Hot exhaust side enthalpy. [J/kg]

        P_ex_H: Hot exhaust side pressure. [Pa]

        Q_dot: Heat Exchanger's heat duty. [W]
    """

    def __init__(self):
        super().__init__()
        self.su_C = MassConnector() # Working fluid supply
        self.su_H = MassConnector() # Secondary fluid supply
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.DP_c = 0
        self.DP_h = 0

        self.Q = HeatConnector()
        self.guesses = {}

    def get_required_inputs(self):
        # Return a list of required inputs
        return ['fluid_C', 'h_su_C', 'P_su_C', 'm_dot_C', 'fluid_H', 'h_su_H', 'P_su_H', 'm_dot_H']
    
    def get_required_parameters(self):
        return [
            'Pinch', # pinch point
            'Delta_T_sh_sc', # Superheating or subcooling
            'HX_type' # Evaporator or condenser
        ]
    
    def get_required_guesses(self):
        return [ 'P_sat' ]

    def equivalent_effectiveness(self):
        
        # External Pinching (your current Q_dot_max)
        if self.params['HX_type'] == "evaporator":
            self.AS_H.update(CoolProp.PT_INPUTS, self.su_H.p, self.su_C.T)
            H_h_id = self.AS_H.hmass()
            
            self.AS_C.update(CoolProp.PT_INPUTS, self.P_sat, self.su_H.T)
            H_c_id = self.AS_C.hmass()
        
        elif self.params['HX_type'] == "condenser":
            self.AS_H.update(CoolProp.PT_INPUTS, self.P_sat, self.su_C.T)
            H_h_id = self.AS_H.hmass()
            
            self.AS_C.update(CoolProp.PT_INPUTS, self.su_C.p, self.su_H.T)
            H_c_id = self.AS_C.hmass()
        
        Q_dot_maxh = self.su_H.m_dot*(self.su_H.h- H_h_id)
        Q_dot_maxc = self.su_C.m_dot*(H_c_id-self.su_C.h)
        self.Q_dot_max_ext = np.min([Q_dot_maxh,Q_dot_maxc])   # rename to *_ext

        # ----------------------------------------------------------------------
        # Internal Pinching: find Q_dot_max_internal by setting pinch ~ 0 K
        # ----------------------------------------------------------------------
        Q_dot_save = self.Q  # current duty at design pinch

        # Save current state
        Pinch_save = self.params['Pinch']
        P_sat_save = self.P_sat

        su_C_p_save = self.su_C.p
        su_H_p_save = self.su_H.p
        su_C_T_save = self.su_C.T
        su_H_T_save = self.su_H.T
        su_C_h_save = self.su_C.h
        su_H_h_save = self.su_H.h

        # Use a very small pinch as "zero" to avoid numerical issues
        self.params['Pinch'] = 1e-2  # K

        # Determine the type of heat exchanger and set the initial guess for pressure
        if self.params['HX_type'] == 'evaporator':
            guess_T_sat = self.su_H.T - self.params['Pinch'] - self.params['Delta_T_sh_sc']
            
            # print(f"guess_T_sat: {guess_T_sat}")
            self.AS_C.update(CoolProp.QT_INPUTS,0.5,guess_T_sat)
            P_ev_guess = self.AS_C.p() # Guess the saturation pressure, first checks if P_sat is in the guesses dictionary, if not it calculates it
            x = [P_ev_guess]
        
            try:
                """EVAPORATOR MODEL"""
                root(self.system_evap, x, method = 'lm', tol=1e-7)
                    
            except Exception as e:
                # Handle any errors that occur during solving
                print(f"Convergence problem in pinch analysis of evaporator model: {e}")
        
        elif self.params['HX_type'] == 'condenser':
            guess_T_sat = self.su_C.T + self.params['Pinch'] + self.params['Delta_T_sh_sc']
            
            self.AS_H.update(CoolProp.QT_INPUTS,0.5,guess_T_sat)
            P_cd_guess = self.AS_H.p() # Guess the saturation pressure, first checks if P_sat is in the guesses dictionary, if not it calculates it
            x = [P_cd_guess]
        
            try:
                """CONDENSER MODEL"""
                fsolve(self.system_cond, x)
                
            except Exception as e:
                # Handle any errors that occur during solving
                print(f"Convergence problem in pinch analysis of evaporator model: {e}")
            
        self.Q_dot_max_int = self.Q

        # Restore original pinch and state by re-solving at design pinch
        self.params['Pinch'] = Pinch_save

        # Restore connector primitive states
        self.su_C.set_p(su_C_p_save)
        self.su_H.set_p(su_H_p_save)
        self.su_C.set_T(su_C_T_save)
        self.su_H.set_T(su_H_T_save)
        self.su_C.h = su_C_h_save
        self.su_H.h = su_H_h_save
        self.P_sat = P_sat_save
        
        self.Q_dot_max = np.min([self.Q_dot_max_ext , self.Q_dot_max_int])

        self.Q = Q_dot_save
        
        if np.isfinite(self.Q_dot_max) and self.Q_dot_max > 0:
            self.epsilon = self.Q / self.Q_dot_max
        else:
            self.epsilon = np.nan

    def system_evap(self, x):
        P_ev = x[0]
        
        PP_list = []
        
        # Ensure the pressure is non-negative
        P_triple = self.AS_C.trivial_keyed_output(CoolProp.iP_triple)
        if P_ev < P_triple:
            P_ev = P_triple + 1
        
        P_crit = self.AS_C.trivial_keyed_output(CoolProp.iP_critical)
        if P_ev > P_crit:
            P_ev = P_crit - 1000
        # Get the temperature of the evaporator based on the pressure and quality
        # Vapor zone
        
        self.AS_C.update(CoolProp.PQ_INPUTS, P_ev, 0.5)
        T_sat_ev = self.AS_C.T()
        self.T_sat_ev = T_sat_ev
        self.su_C.p = P_ev
        
        "Refrigerant side calculations"
        # Liquid zone
        try:
            self.AS_C.update(CoolProp.PT_INPUTS,P_ev,self.su_C.T)
        except:
            self.AS_C.update(CoolProp.PQ_INPUTS,P_ev,0)
            
        h_C_su = self.AS_C.hmass()
        
        self.AS_C.update(CoolProp.PQ_INPUTS,P_ev,0)
        h_C_x0 = self.AS_C.hmass()
        
        self.Q_dot_sc = self.su_C.m_dot * (h_C_x0 - h_C_su)

        # Two-phase zone
        self.AS_C.update(CoolProp.PQ_INPUTS,P_ev,1)
        h_C_x1 = self.AS_C.hmass()

        if self.Q_dot_sc > 0:
            self.Q_dot_tp = self.su_C.m_dot * (h_C_x1 - h_C_x0)
        else:
            self.Q_dot_sc = 0
            self.Q_dot_tp = self.su_C.m_dot * (h_C_x1 - self.su_C.h)
            
        # Vapor zone
        self.T_C_ex = T_sat_ev + self.params['Delta_T_sh_sc']
        self.AS_C.update(CoolProp.PT_INPUTS,P_ev,self.T_C_ex)
        h_C_ex = self.AS_C.hmass()

        if self.Q_dot_tp > 0:
            self.Q_dot_sh = self.su_C.m_dot * (h_C_ex - h_C_x1)
        else:
            self.Q_dot_tp = 0
            self.Q_dot_sh = self.su_C.m_dot * (h_C_ex - self.su_C.h)

        # Total heat transfer
        Q_dot_ev = self.Q_dot_sc + self.Q_dot_tp + self.Q_dot_sh
        
        "Secondary fluid side calculations"
        # First zone - Superheating
        self.h_H_x1 = self.su_H.h - self.Q_dot_sh/self.su_H.m_dot
        self.AS_H.update(CoolProp.HmassP_INPUTS,self.h_H_x1,self.su_H.p)
        self.T_H_x1 = self.AS_H.T()
        
        # Second zone - Phase change
        self.h_H_x0 = self.h_H_x1 - self.Q_dot_tp/self.su_H.m_dot
        self.AS_H.update(CoolProp.HmassP_INPUTS,self.h_H_x0,self.su_H.p)
        self.T_H_x0 = self.AS_H.T()

        # Third zone - Subcooling
        self.h_H_ex = self.h_H_x0 - self.Q_dot_sc/self.su_H.m_dot
        self.AS_H.update(CoolProp.HmassP_INPUTS,self.h_H_ex,self.su_H.p)
        self.T_H_ex = self.AS_H.T()    
        
        PP_list.append(self.T_H_ex - self.su_C.T)
        
        if self.Q_dot_sc > 0:
            PP_list.append(self.T_H_x0 - T_sat_ev)            

        if self.Q_dot_sh > 0:
            PP_list.append(self.su_H.T - self.T_C_ex)
                
        # Calculate pinch point and residual$
        self.PP_array = np.array(PP_list)
        
        self.PPTD = min(self.PP_array)

        # PPTD = min(self.T_H_ex - self.su_C.T, self.T_H_x0 - T_sat_ev, self.T_H_x1 - T_sat_ev, self.su_H.T - self.T_C_ex)

        self.res = abs(self.PPTD - self.params['Pinch'])
        
        # Update the state of the working fluid
        self.Q = Q_dot_ev
        self.P_sat = P_ev
        
        return self.res
    
    def system_cond(self, x):
        P_cd = x[0]
                
        PP_list = []
        
        # Ensure the pressure is non-negative
        P_triple = self.AS_H.trivial_keyed_output(CoolProp.iP_triple)
        if P_cd < P_triple:
            P_cd = P_triple*2
        
        P_crit = self.AS_H.trivial_keyed_output(CoolProp.iP_critical)
        if P_cd > P_crit:
            P_cd = P_crit - 1000
                
        # Get the temperature of the condenser based on pressure and quality
        self.AS_H.update(CoolProp.PQ_INPUTS,P_cd,0.5)
        T_sat_cd = self.AS_H.T()            
        self.T_sat_cd = T_sat_cd
        self.su_H.p = P_cd

        "Refrigerant side calculations"
        # Vapor zone
        try:
            self.AS_H.update(CoolProp.PT_INPUTS,P_cd,self.su_H.T)
        except:
            self.AS_H.update(CoolProp.PQ_INPUTS,P_cd,1)
        
        h_H_su = self.AS_H.hmass()
        self.su_H.h = h_H_su
        
        self.AS_H.update(CoolProp.PQ_INPUTS,P_cd,1)
        h_H_x1 = self.AS_H.hmass()
        
        if T_sat_cd < self.su_H.T:
            self.Q_dot_sh = self.su_H.m_dot * (h_H_su - h_H_x1)
        else:
            self.Q_dot_sh = 0

        # Two-phase zone
        self.AS_H.update(CoolProp.PQ_INPUTS,P_cd,0)
        h_H_x0 = self.AS_H.hmass()
        
        if self.Q_dot_sh > 0:
            self.Q_dot_tp = self.su_H.m_dot * (h_H_x1 - h_H_x0)
        else:
            self.Q_dot_sh = 0
            self.Q_dot_tp = self.su_H.m_dot * (self.su_H.h - h_H_x0)
        
        # Liquid zone
        self.T_H_ex = T_sat_cd - self.params['Delta_T_sh_sc']
        self.AS_H.update(CoolProp.PT_INPUTS,P_cd,self.T_H_ex)
        h_H_ex = self.AS_H.hmass()
        
        if self.Q_dot_tp > 0:
            self.Q_dot_sc = self.su_H.m_dot * (h_H_x0 - h_H_ex)
        else:
            self.Q_dot_tp = 0
            self.Q_dot_sc = self.su_H.m_dot * (self.su_H.h - h_H_ex)

        # Total heat transfer
        Q_dot_cd = self.Q_dot_sh + self.Q_dot_tp + self.Q_dot_sc
        
        "Secondary fluid side calculations"
        # First zone
        self.h_C_x0 = self.su_C.h + self.Q_dot_sc/self.su_C.m_dot
        self.AS_C.update(CoolProp.HmassP_INPUTS,self.h_C_x0,self.su_C.p)
        self.T_C_x0 = self.AS_C.T()
        
        # Second zone
        self.h_C_x1 = self.h_C_x0 + self.Q_dot_tp/self.su_C.m_dot
        self.AS_C.update(CoolProp.HmassP_INPUTS,self.h_C_x1,self.su_C.p)
        self.T_C_x1 = self.AS_C.T()

        # Third zone
        self.h_C_ex = self.h_C_x1 + self.Q_dot_sh/self.su_C.m_dot
        self.AS_C.update(CoolProp.HmassP_INPUTS,self.h_C_ex,self.su_C.p)
        self.T_C_ex = self.AS_C.T()
        
        PP_list.append(self.su_H.T - self.T_C_ex)
        
        if self.Q_dot_sh > 0:
            PP_list.append(T_sat_cd - self.T_C_x1)            

        if self.Q_dot_sc > 0:
            PP_list.append(self.T_H_ex - self.su_C.T)
                
        # Calculate pinch point and residual$
        self.PP_array = np.array(PP_list)
        
        self.PPTD = min(self.PP_array)
        
        # Calculate residual
        self.res = abs(self.PPTD  - self.params['Pinch'])
        
        # Update the state of the working fluid
        self.Q = Q_dot_cd
        self.P_sat = P_cd
        
        return self.res

    def solve(self):
        # Ensure all required checks are performed

        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            print("HTX IS NOT CALCULABLE")
            return
        
        fluid_C = self.su_C.fluid  # Extract cold fluid name
        self.AS_C = AbstractState("BICUBIC&HEOS", fluid_C)  # Create a reusable state object
        
        fluid_H = self.su_H.fluid  # Extract hot fluid name
        self.AS_H = AbstractState("BICUBIC&HEOS", fluid_H)  # Create a reusable state object

        if 'DP_h' in self.params:
            self.DP_h = self.params['DP_h']

        if 'DP_c' in self.params:
            self.DP_c = self.params['DP_c']

        # Determine the type of heat exchanger and set the initial guess for pressure
        if self.params['HX_type'] == 'evaporator':
            guess_T_sat = self.su_H.T - self.params['Pinch'] - self.params['Delta_T_sh_sc']
            
            # print(f"guess_T_sat: {guess_T_sat}")
            self.AS_C.update(CoolProp.QT_INPUTS,0.5,guess_T_sat)
            P_ev_guess = self.AS_C.p() # Guess the saturation pressure, first checks if P_sat is in the guesses dictionary, if not it calculates it
            x = [P_ev_guess]

            try:
                """EVAPORATOR MODEL"""
                root(self.system_evap, x, method = 'lm', tol=1e-7)
                
                """Update connectors after the calculations"""
                self.update_connectors()

                # Mark the model as solved if successful
                if self.res < 1e-2:
                    self.solved = True
                else:
                    print("System not solved according to specified tolerance in Evaporator")
                    self.solved = False
                    
            except Exception as e:
                # Handle any errors that occur during solving
                self.solved = False
                print(f"Convergence problem in evaporator model: {e}")

        elif self.params['HX_type'] == 'condenser':
            guess_T_sat = self.su_C.T + self.params['Pinch'] + self.params['Delta_T_sh_sc']
            
            self.AS_H.update(CoolProp.QT_INPUTS,0.5,guess_T_sat)
            P_cd_guess = self.AS_H.p() # Guess the saturation pressure, first checks if P_sat is in the guesses dictionary, if not it calculates it
            x = [P_cd_guess]

            try:
                """CONDENSER MODEL"""
                fsolve(self.system_cond, x)

                """Update connectors after the calculations"""
                self.update_connectors()

                # Mark the model as solved if successful
                self.solved = True
            except Exception as e:
                # Handle any errors that occur during solving
                self.solved = False
                print(f"Convergence problem in condenser model: {e}")

        """Compute HX Equivalent Efficiency"""
        # External Pinching
        

        # Internal Pinching
        


    def update_connectors(self):
        
        "Mass Connectors"

        if self.params['HX_type'] == 'evaporator':

            self.su_C.set_p(self.P_sat + self.DP_c)

            self.ex_C.set_fluid(self.su_C.fluid)
            self.ex_C.set_T(self.T_C_ex)
            self.ex_C.set_p(self.P_sat)
            self.ex_C.set_m_dot(self.su_C.m_dot)

            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_m_dot(self.su_H.m_dot)
            self.ex_H.set_p(self.su_H.p - self.DP_h)
            self.ex_H.set_T(self.T_H_ex)
            
            "Heat conector"
            self.Q.set_Q_dot(self.Q)

        else: 

            self.su_H.set_p(self.P_sat + self.DP_h)

            self.ex_H.set_fluid(self.su_H.fluid)
            self.ex_H.set_T(self.T_H_ex)
            self.ex_H.set_p(self.P_sat)
            self.ex_H.set_m_dot(self.su_H.m_dot)

            self.ex_C.set_fluid(self.su_C.fluid)
            self.ex_C.set_m_dot(self.su_C.m_dot)
            self.ex_C.set_p(self.su_C.p - self.DP_c)
            self.ex_C.set_T(self.T_C_ex)
            
            "Heat conector"
            self.Q.set_Q_dot(self.Q)
            
    def print_results(self):
        print("=== Heat Exchanger Results ===")
        print(f"Q: {self.Q.Q_dot}")

        if self.params['HX_type'] == 'evaporator':
            print(f"P_sat: {self.su_C.p}")
        else:
            print(f"P_sat: {self.su_H.p}")

        print("======================")

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, P={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, P={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, P={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, P={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q.Q_dot}")
        print("======================")

    def plot_disc(self):
        import matplotlib.pyplot as plt
        
        if self.params['HX_type'] == 'condenser':
            plt.figure()
            
            plt.plot([0, self.Q_dot_sh]                          , [self.su_H.T, self.T_sat_cd]  , 'r', label='H')
            plt.plot([self.Q_dot_sh, self.Q_dot_sh+self.Q_dot_tp], [self.T_sat_cd, self.T_sat_cd], 'r')
            plt.plot([self.Q_dot_sh+self.Q_dot_tp, self.Q]       , [self.T_sat_cd, self.ex_H.T], 'r')

            plt.plot([0, self.Q_dot_sh]                          , [self.ex_C.T, self.T_C_x1], 'b', label='C')
            plt.plot([self.Q_dot_sh, self.Q_dot_sh+self.Q_dot_tp], [self.T_C_x1, self.T_C_x0], 'b')
            plt.plot([self.Q_dot_sh+self.Q_dot_tp, self.Q]       , [self.T_C_x0, self.su_C.T], 'b')

            plt.grid()
            plt.legend()
            plt.show()

        if self.params['HX_type'] == 'evaporator':
            plt.figure()
            
            plt.plot([0, self.Q_dot_sc]                          , [self.su_C.T, self.T_sat_ev]  , 'b', label='C')
            plt.plot([self.Q_dot_sc, self.Q_dot_sc+self.Q_dot_tp], [self.T_sat_ev, self.T_sat_ev], 'b')
            plt.plot([self.Q_dot_sc+self.Q_dot_tp, self.Q]       , [self.T_sat_ev, self.ex_C.T]  , 'b')

            plt.plot([0, self.Q_dot_sc]                          , [self.ex_H.T, self.T_H_x0], 'r', label='H')
            plt.plot([self.Q_dot_sc, self.Q_dot_sc+self.Q_dot_tp], [self.T_H_x0, self.T_H_x1], 'r')
            plt.plot([self.Q_dot_sc+self.Q_dot_tp, self.Q]       , [self.T_H_x1, self.su_H.T], 'r')

            plt.grid()
            plt.legend()
            plt.show()