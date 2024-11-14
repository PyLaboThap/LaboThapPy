# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:09:18 2024

@author: Elise Neven
@email: elise.neven@uliege.be

"""

from CoolProp.CoolProp import PropsSI

class MassConnector:
    """
    A class to handle and calculate fluid properties based on given state variables.

    **Attributes**:
    
        fluid : str, optional
            The name of the fluid.
        m_dot : float, optional
            Mass flow rate in kg/s.
        V_dot : float, optional
            Volume flow rate in m^3/s.
        T : float, optional
            Temperature in Kelvin.
        p : float, optional
            Pressure in Pascal.
        h : float, optional
            Specific enthalpy in J/kg.
        s : float, optional
            Specific entropy in J/kg/K.
        D : float, optional
            Mass density in kg/m^3.
        x : float, optional
            Quality in kg/kg.
        cp : float, optional
            Specific heat capacity in J/kg/K.
        completely_known : bool
            True if all properties and mass flow rate are known, otherwise False.
        state_known : bool
            True if all state properties are known, otherwise False.
        variables_input : list of lists
            A list of the state variables used to define the fluid state. Each entry is a list of [variable_name, value].

    **Methods**:

        __init__(self, fluid=None):
            Initializes the MassConnector object with optional fluid.
        
        reset(self):
            Resets all properties to None and flags to False.

        check_completely_known(self):
            Checks if all necessary properties and mass flow rate are known to determine if the state is completely known.

        calculate_properties(self):
            Calculates fluid properties using CoolProp based on known variables.

        is_two_phase_PT(self):
            Determines if the current state is in the two-phase region based on temperature and pressure.

        set_properties(self, kwargs):
            Sets various properties of the fluid. Accepts named arguments for 'fluid', 'm_dot', 'V_dot', 'T', 'P', 'H', 'S', 'D', 'x', and 'cp'.

        set_fluid(self, value):
            Sets the fluid type and checks if it's valid.

        set_m_dot(self, value):
            Sets the mass flow rate and ensures volume flow rate is not set simultaneously.

        set_V_dot(self, value):
            Sets the volume flow rate and ensures mass flow rate is not set simultaneously.

        set_T(self, value):
            Sets the temperature and updates the list of known variables.

        set_p(self, value):
            Sets the pressure and updates the list of known variables.

        set_h(self, value):
            Sets the specific enthalpy and updates the list of known variables.

        set_s(self, value):
            Sets the specific entropy and updates the list of known variables.

        set_D(self, value):
            Sets the mass density and updates the list of known variables.

        set_x(self, value):
            Sets the quality and updates the list of known variables.

        set_cp(self, value):
            Sets the specific heat capacity without calculating other properties.

        print_resume(self, unit_T='K', unit_p='Pa'):
            Prints a summary of the current properties in specified units.

    **Exceptions**:
    
    ValueError
        Raised if an unsupported fluid name is used or if required variables for property calculation are missing.
    """

    def __init__(self, fluid=None):

        self.completely_known = False # True if all the properties and the mass flow rate are known
        self.state_known = False      # True if all the properties are known
        self.variables_input = []     # List of the variables used to define the state of the fluid
        
        self.fluid = fluid           # Fluid name
        self.m_dot = None           # Mass flow rate [kg/s]
        self.V_dot = None           # Volume flow rate [m^3/s]
        self.T = None               # Temperature [K]
        self.p = None               # Pressure [Pa]
        self.h = None               # Specific enthalpy [J/kg]
        self.s = None               # Specific entropy [J/kg/K]
        self.D = None               # Mass density [kg/m^3]
        self.x = None               # Quality [kg/kg]
        self.cp = None              # Specific heat capacity [J/kg/K]

    def reset(self):
        self.completely_known = False
        self.state_known = False
        self.variables_input = []
        self.m_dot = None
        self.V_dot = None
        self.T = None
        self.p = None
        self.h = None
        self.s = None
        self.D = None
        self.x = None
        self.cp = None
    
    def check_completely_known(self):
        if self.fluid != None:
            if len(self.variables_input)>2:
                print("Error: Too many state variables")
            elif len(self.variables_input)<2:
                pass
            elif len(self.variables_input)==2:
                self.calculate_properties()

            if (self.m_dot != None or self.V_dot !=None) and self.state_known:
                if self.m_dot != None:
                    self.V_dot = self.m_dot / self.D * 3600
                elif self.V_dot != None:
                    self.m_dot = self.V_dot * self.D / 3600
                self.completely_known = True
        else:
            pass

    def calculate_properties(self):
        if 'INCOMP' in self.fluid: # If the fluid is incompressible, the quality makes no sense
            try:
                self.T = PropsSI('T', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.p = PropsSI('P', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.h = PropsSI('H', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.s = PropsSI('S', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.D = PropsSI('D', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.cp = PropsSI('CPMASS', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.state_known = True
            except:
                print("Error: This pair of inputs is not yet supported.")

        else:
            try:
                self.x = PropsSI('Q', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.T = PropsSI('T', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.p = PropsSI('P', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.h = PropsSI('H', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.s = PropsSI('S', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.D = PropsSI('D', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.cp = PropsSI('CPMASS', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                self.state_known = True
            except:
                try:
                    if self.is_two_phase_PT():
                        print('Cannot use the temperature and the pressure to calculate the properties in two-phase.')
                    else:
                        self.x = PropsSI('Q', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                        self.h = PropsSI('H', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                        self.p = PropsSI('P', 'Q', self.x, 'H', self.h, self.fluid)
                        self.T = PropsSI('T', 'Q', self.x, 'H', self.h, self.fluid)
                        self.s = PropsSI('S', 'Q', self.x, 'H', self.h, self.fluid)
                        self.D = PropsSI('D', 'Q', self.x, 'H', self.h, self.fluid)
                        self.cp = PropsSI('CPMASS', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                        self.state_known = True
                except:
                    print("Error: This pair of inputs is not yet supported.")
    
    def is_two_phase_PT(self):
        if 'T' in [var[0] for var in self.variables_input] and 'P' in [var[0] for var in self.variables_input]:
            T = self.T if self.T else [var[1] for var in self.variables_input if var[0] == 'T'][0]
            p = self.p if self.p else [var[1] for var in self.variables_input if var[0] == 'P'][0]
            
            T_sat = PropsSI('T', 'P', p, 'Q', 0.5, self.fluid)
            return T <= T_sat
        return False

    def set_properties(self, **kwargs):
        for key, value in kwargs.items():
            if key == 'fluid':
                self.set_fluid(value)
            elif key == 'm_dot':
                self.set_m_dot(value)
            elif key == 'V_dot':
                self.set_V_dot(value)
            elif key == 'T':
                self.set_T(value)
            elif key == 'P':
                self.set_p(value)
            elif key == 'H':
                self.set_h(value)
            elif key == 'S':
                self.set_s(value)
            elif key == 'D':
                self.set_D(value)
            elif key == 'x':
                self.set_x(value)
            elif key == 'cp':
                self.set_cp(value)
            else:
                print(f"Error: Invalid property '{key}'")

    def set_fluid(self, value):
        if self.fluid != None:
            pass
        elif self.fluid == None:
            try:
                if 'INCOMP' in value:
                    # Check if the fluid name is correct for incompressible fluids
                    PropsSI('D', 'T', 300.0, 'P', 101325, value)
                else:
                    # Check if the fluid name is correct for compressible fluids
                    PropsSI('M', value)
                self.fluid = value
            except:
                print("Error: Incorrect fluid name:", value)
        self.check_completely_known()


    def set_m_dot(self, value):
        self.m_dot = value
        self.V_dot = None #makes sure that the volume flow rate and the mass flow rate are not both known
        self.check_completely_known()

    def set_V_dot(self, value):
        self.V_dot = value
        self.m_dot = None #makes sure that the volume flow rate and the mass flow rate are not both known
        self.check_completely_known()
        
    def set_T(self, value):
        # print('set_T', value)
        if self.T != None: # If the temperature is already known, update the value and the corresponding variable in the list
            self.T = value
            for i, var in enumerate(self.variables_input):
                if var[0] == 'T':
                    self.variables_input[i][1] = value
                    break
        else:              # If the temperature is not known, set the value and add the variable to the list
            self.T = value
            self.variables_input = self.variables_input+[['T',value]]
        self.check_completely_known()
        
    def set_p(self, value):
        # print('set_p', value)
        if self.p != None: # If the pressure is already known, update the value and the corresponding variable in the list
            self.p = value
            for i, var in enumerate(self.variables_input):
                if var[0] == 'P':
                    self.variables_input[i][1] = value
                    break
        else:             # If the pressure is not known, set the value and add the variable to the list
            self.p = value
            self.variables_input = self.variables_input+[['P',value]]
        self.check_completely_known()
        
    def set_h(self, value):
        if self.h != None: # If the specific enthalpy is already known, update the value and the corresponding variable in the list
            self.h = value
            for i, var in enumerate(self.variables_input):
                if var[0] == 'H':
                    self.variables_input[i][1] = value
                    break
        else:            # If the specific enthalpy is not known, set the value and add the variable to the list
            self.h = value
            self.variables_input = self.variables_input+[['H',value]]
        self.check_completely_known()
        
    def set_s(self, value):
        if self.s != None: # If the specific entropy is already known, update the value and the corresponding variable in the list
            self.s = value
            for i, var in enumerate(self.variables_input):
                if var[0] == 'S':
                    self.variables_input[i][1] = value
                    break
        else:           # If the specific entropy is not known, set the value and add the variable to the list
            self.s = value
            self.variables_input = self.variables_input+[['S',value]]
        self.check_completely_known()
        
    def set_D(self, value):
        if self.D != None: # If the mass density is already known, update the value and the corresponding variable in the list
            self.D = value
            for i, var in enumerate(self.variables_input):
                if var[0] == 'D':
                    self.variables_input[i][1] = value
                    break
        else:           # If the mass density is not known, set the value and add the variable to the list
            self.D = value
            self.variables_input = self.variables_input+[['D',value]]
        self.check_completely_known()
            
    def set_x(self, value):
        if self.x != None: # If the quality is already known, update the value and the corresponding variable in the list
            self.x = value
            for i, var in enumerate(self.variables_input):
                if var[0] == 'Q':
                    self.variables_input[i][1] = value
                    break
        else:          # If the quality is not known, set the value and add the variable to the list
            self.x = value
            self.variables_input = self.variables_input+[['Q',value]]
        self.check_completely_known()

    def set_cp(self, value):
        if self.cp != None: # If the cp is already known, update the value and the corresponding variable in the list
            self.cp = value
        else:          # If the quality is not known, set the value and add the variable to the list
            self.cp = value # We don't want to calculate in this case
            
            
    def print_resume(self, unit_T='K', unit_p='Pa'):
        """
        Parameters
        ----------
        unit_T = Temperature unit: 'K' or 'C'
        unit_p = Temperature unit: 'Pa' or 'bar'
        
        """
        
        print("Fluid: " + self.fluid + "")
        print("Mass flow rate: " + str(self.m_dot) + "[kg/s]")
        print("Volume flow rate: " + str(self.V_dot) + "[m^3/h]")
        
        if unit_T == 'K':
            print("Temperature: " + str(self.T) + "[K]")
        elif unit_T == 'C':
            print("Temperature: " + str(self.T-273.15) + "[°C]")
        else:
            print("Error: Wrong argument unit_T in the method print_resume")

        if unit_p == 'Pa':
            print("Pressure: " + str(self.p) + "[Pa]")
        elif unit_p == 'bar':
            print("Pressure: " + str(self.p/1e5) + "[bar]")
        else:
            print("Error: Wrong argument unit_p in the method print_resume")
            
        
        print("Spec. enthalpy: " + str(self.h) + "[J/kg]")
        
        print("Spec. entropy: " + str(self.s) + "[J/kg/K]")
        
        print("Mass density: " + str(self.D) + "[kg/m^3]")
        print("Quality: " + str(self.x) + "[-]")

