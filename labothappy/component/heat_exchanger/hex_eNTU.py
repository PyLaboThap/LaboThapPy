# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:32:39 2024

@author: Basile
"""

from component.base_component import BaseComponent

from correlations.convection.pipe_htc import gnielinski_pipe_htc 
from toolbox.geometries.heat_exchanger.e_NTU import e_NTU

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI


class HXeNTU(BaseComponent):
    
    """
    Component: Heat Exchanger

    Model: ε-NTU (Effectiveness - Number of Transfer Units) method.

    **Description**:

        This component models a heat exchanger using the ε-NTU method, a widely used approach for estimating heat transfer performance 
        in steady-state conditions when outlet temperatures are not known a priori. It calculates heat transfer based on fluid properties, 
        flow configuration, and geometry using thermal resistances and heat transfer correlations.

        The model is applicable to various geometries (e.g., pipe-type, plate-type) and requires fluid and geometric properties. 
        The thermal effectiveness is computed via an external ε-NTU correlation, which supports multiple flow configurations 
        (e.g., CounterFlow, ParallelFlow, CrossFlow).

    **Assumptions**:

        - Steady-state operation.
        - No heat loss to the environment.
        - No pressure drop considered (isenthalpic mixing assumed).
        - Thermophysical properties are evaluated at average temperatures.
        - No phase change within the exchanger.

    **Connectors**:

        su_hot (MassConnector): Hot fluid inlet connector.
        su_cold (MassConnector): Cold fluid inlet connector.
        ex_hot (MassConnector): Hot fluid outlet connector.
        ex_cold (MassConnector): Cold fluid outlet connector.
        Q_dot (HeatConnector): Connector for total heat transfer rate.

    **Parameters**:
            
        Flow_Type : Flow configuration of the fluid ('CounterFlow', 'CrossFlow', 'Shell&Tube', 'ParallelFlow') [-]
        A_htx: Total heat exchange area [m²]
        L_HTX: Length of the heat exchanger [m]
        V_HTX: Volume of the heat exchanger [m³]
        A_canal_H: Cross-sectional area of hot fluid channels [m²]
        A_canal_C: Cross-sectional area of cold fluid channels [m²]
        D_h: Hydraulic diameter [m]
        k_plate: Thermal conductivity of the separating plate [W/m.K]
        t_plate: Thickness of the separating plate [m]
        n_plates: Number of plates [-]
        co_pitch: Plate corrugation pitch [m]
        chevron_angle: Plate chevron angle [degrees]
        fouling: Fouling resistance [m².K/W]

    **Inputs**:

        T_su_H: Hot fluid inlet temperature [K]
        p_su_H: Hot fluid inlet pressure [Pa]
        h_su_H: Hot fluid inlet enthalpy [J/kg]
        fluid_su_H: Hot fluid identifier [-]
        m_dot_su_H: Hot fluid mass flow rate [kg/s]

        T_su_C: Cold fluid inlet temperature [K]
        p_su_C: Cold fluid inlet pressure [Pa]
        h_su_C: Cold fluid inlet enthalpy [J/kg]
        fluid_su_C: Cold fluid identifier [-]
        m_dot_su_C: Cold fluid mass flow rate [kg/s]

    **Outputs**:

        h_ex_H: Hot fluid outlet enthalpy [J/kg]
        p_ex_H: Hot fluid outlet pressure [Pa]
        h_ex_C: Cold fluid outlet enthalpy [J/kg]
        p_ex_C: Cold fluid outlet pressure [Pa]
        Q_dot: Heat transfer rate [W]
"""
    
    
    def __init__(self):
        super().__init__()
        self.su_hot = MassConnector()
        self.su_cold = MassConnector()

        self.ex_hot = MassConnector()
        self.ex_cold = MassConnector() # Mass_connector

        self.Q_dot = HeatConnector()

    def get_required_inputs(self):
        
        if self.inputs == {}:
            # Hot Fluid
            if self.su_hot.T is not None:
                self.inputs['Hsu_T'] = self.su_hot.T
            elif self.su_hot.h is not None:
                self.inputs['Hsu_h'] = self.su_hot.h
            if self.su_hot.p is not None:
                self.inputs['Hsu_p'] = self.su_hot.p
            if self.su_hot.fluid is not None:
                self.inputs['Hsu_fluid'] = self.su_hot.fluid
            if self.su_hot.m_dot is not None:
                self.inputs['Hsu_m_dot'] = self.su_hot.m_dot
                
            # Cold Fluid                
            if self.su_cold.T is not None:
                self.inputs['Csu_T'] = self.su_cold.T
            elif self.su_cold.h is not None:
                self.inputs['Csu_h'] = self.su_cold.h
            if self.su_cold.p is not None:
                self.inputs['Csu_p'] = self.su_cold.p
            if self.su_cold.fluid is not None:
                self.inputs['Csu_fluid'] = self.su_cold.fluid
            if self.su_cold.m_dot is not None:
                self.inputs['Csu_m_dot'] = self.su_cold.m_dot
                
        if self.inputs != {}:
            # Hot Fluid
            self.su_hot.set_fluid(self.inputs['Hsu_fluid'])
            if 'Hsu_T' in self.inputs:
                self.su_hot.set_T(self.inputs['Hsu_T'])
            elif 'Hsu_h' in self.inputs:
                self.su_hot.set_h(self.inputs['Hsu_h'])
            if 'Hsu_p' in self.inputs:
                self.su_hot.set_p(self.inputs['Hsu_p'])
            if 'Hsu_m_dot' in self.inputs:
                self.su_hot.set_m_dot(self.inputs['Hsu_m_dot'])

            # Cold Fluid
            self.su_cold.set_fluid(self.inputs['Csu_fluid'])
            if 'Csu_T' in self.inputs:
                self.su_cold.set_T(self.inputs['Csu_T'])
            elif 'Csu_h' in self.inputs:
                self.su_cold.set_h(self.inputs['Csu_h'])
            if 'Csu_p' in self.inputs:
                self.su_cold.set_p(self.inputs['Csu_p'])
            if 'Csu_m_dot' in self.inputs:
                self.su_cold.set_m_dot(self.inputs['Csu_m_dot'])

        return ['Hsu_p', 'Hsu_T', 'Hsu_m_dot', 'Hsu_fluid', 'Csu_p', 'Csu_T', 'Csu_m_dot', 'Csu_fluid']

    def get_required_parameters(self):
        """ Returns the list of required parameters to describe the geometry and physical configuration """
        return ['A_htx', 'L_HTX', 'V_HTX', 'Flow_Type',
                'A_canal_h', 'A_canal_c', 'D_h',
                'k_plate', 't_plate', 'n_plates',
                'co_pitch', 'chevron_angle', 'fouling']
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - H_su: fluid={self.su_hot.fluid}, T={self.su_hot.T}, p={self.su_hot.p}, m_dot={self.su_hot.m_dot}")
        print(f"  - C_su: fluid={self.su_cold.fluid}, T={self.su_cold.T}, p={self.su_cold.p}, m_dot={self.su_cold.m_dot}")

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

    def solve(self):
        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:
            
            # Detect Phase change
            # self.detect_phase_change()
            
            # Calcul de C_r
            cp_h = PropsSI('C', 'H', self.su_hot.h, 'P', self.su_hot.p, self.su_hot.fluid)
            cp_c = PropsSI('C', 'H', self.su_cold.h, 'P', self.su_cold.p, self.su_cold.fluid)
            
            C_h = cp_h*self.su_hot.m_dot #Heat capacity rate
            C_c = cp_c*self.su_cold.m_dot
            
            C_min = min(C_h, C_c)
            C_max = max(C_h, C_c)
            C_r = C_min/C_max # Heat capacity ratio
                        
            # Calcul de NTU
            T_w = (self.su_hot.T + self.su_cold.T)/2
            
            # --- Heat transfer coefficient estimation using Gnielinski correlation ---
            mu_h, Pr_h, k_h = PropsSI(('V','PRANDTL','L'), 'H', self.su_hot.h, 'P', self.su_hot.p, self.su_hot.fluid)
            mu_c, Pr_c, k_c = PropsSI(('V','PRANDTL','L'), 'H', self.su_cold.h, 'P', self.su_cold.p, self.su_cold.fluid)
            
            G_h = self.su_hot.m_dot/self.params['A_canal_h']
            G_c = self.su_cold.m_dot/self.params['A_canal_c']
            
            
            h_h = gnielinski_pipe_htc(mu_h, Pr_h, Pr_h, k_h, G_h, self.params['D_h'], self.params['L_HTX'])[0]
            h_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c, k_c, G_c, self.params['D_h'], self.params['L_HTX'])[0]

            # --- Global heat transfer coefficient (AU)  ---
            AU = (1/(self.params['A_htx']*h_h) + 1/(self.params['A_htx']*h_c) + self.params['t_plate']/(self.params['k_plate']*self.params['A_htx']) + self.params['fouling']/self.params['A_htx'])**(-1)         

            NTU = AU/C_min
                        
            # --- Calculate effectiveness from NTU correlation ---
            eps = e_NTU(NTU, C_r, self.params)

                        
            # --- Estimate maximum heat transfer Q(ideal case with infinite area) ---
            h_c_Th = PropsSI('H','T',self.su_hot.T,'P',self.su_cold.p,self.su_cold.fluid)
            h_h_Tc = PropsSI('H','T',self.su_cold.T,'P',self.su_hot.p,self.su_hot.fluid)
            
            DH_pc_c = PropsSI('H','Q',1,'P',self.su_cold.p,self.su_cold.fluid) - PropsSI('H','Q',0,'P',self.su_cold.p,self.su_cold.fluid)

            # Special case for incompressibles
            if "INCOMP" not in self.su_hot.fluid:
                DH_pc_h = PropsSI('H','Q',1,'P',self.su_hot.p,self.su_hot.fluid) - PropsSI('H','Q',0,'P',self.su_hot.p,self.su_hot.fluid)
            else:
                DH_pc_h = 0
            
            Qmax_c = self.su_cold.m_dot*((h_c_Th - self.su_cold.h))
            Qmax_h = self.su_hot.m_dot*((self.su_hot.h - h_h_Tc))
                        
            Qmax = min(Qmax_c, Qmax_h)
            
            Q = eps*Qmax  # Actual heat exchanged

            # --- Set exhaust states (new enthalpies) and link to connectors ---
            self.ex_hot.set_properties(H = self.su_hot.h - Q/self.su_hot.m_dot, fluid = self.su_hot.fluid, m_dot = self.su_hot.m_dot, P = self.su_hot.p)
            self.ex_cold.set_properties(H = self.su_cold.h + Q/self.su_cold.m_dot, fluid = self.su_cold.fluid, m_dot = self.su_cold.m_dot, P = self.su_cold.p)

            self.Q_dot.set_Q_dot(Q)

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
            print(f"  - H_ex: fluid={self.ex_hot.fluid}, T={self.ex_hot.T}, p={self.ex_hot.p}, m_dot={self.ex_hot.m_dot}")
            print(f"  - C_ex: fluid={self.ex_cold.fluid}, T={self.ex_cold.T}, p={self.ex_cold.p}, m_dot={self.ex_cold.m_dot}")
            print(f"  - Q_dot: {self.Q_dot.Q_dot}")

        else:
            print("Heat Exchanger component is not defined. Ensure it is solved first.")

