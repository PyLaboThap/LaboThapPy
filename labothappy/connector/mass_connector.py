# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:09:18 2024

@author: Elise Neven
@email: elise.neven@uliege.be

"""

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import warnings

class MassConnector:
    """
    MassConnector is a class that handle and calculate fluid properties based on given state variables.

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
        SC : float, optional
            Subcooling in K.
        SH : float, optional
            Superheating in K.
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

        set_SC(self, value):
            Sets the subcooling and updates the list of known variables.
        
        set_SH(self, value):
            Sets the superheating and updates the list of known variables.

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
        
        if fluid != None:
            try:
                try: 
                    self.AS = CP.AbstractState("BICUBIC&HEOS", fluid) 
                    self.fluid = fluid
                except:
                    if "INCOMP" in fluid:
                        fluid = fluid.split(":")[-1]
                    self.AS = CP.AbstractState("INCOMP", fluid)
                    self.fluid = fluid
            except:
                warnings.warn("Error: Incorrect fluid name:", fluid)
        else:
            self.fluid = None       # Fluid name
            
        self.m_dot = None           # Mass flow rate [kg/s]
        self.V_dot = None           # Volume flow rate [m^3/s]
        self.T = None               # Temperature [K]
        self.p = None               # Pressure [Pa]
        self.h = None               # Specific enthalpy [J/kg]
        self.s = None               # Specific entropy [J/kg/K]
        self.D = None               # Mass density [kg/m^3]
        self.x = None               # Quality [kg/kg]
        self.cp = None              # Specific heat capacity [J/kg/K]
        self.SC = None              # Subcooling [K]
        self.SH = None              # Superheating [K]
         
        self.CP_map = {
            # --- Mass-based inputs ---
            ('D', 'H'): CP.DmassHmass_INPUTS, ('D', 'P'): CP.DmassP_INPUTS, ('D', 'Q'): CP.DmassQ_INPUTS, ('D', 'S'): CP.DmassSmass_INPUTS, ('D', 'T'): CP.DmassT_INPUTS, ('D', 'U'): CP.DmassUmass_INPUTS,
        
            ('H', 'P'): CP.HmassP_INPUTS, ('H', 'Q'): CP.HmassQ_INPUTS, ('H', 'S'): CP.HmassSmass_INPUTS, ('H', 'T'): CP.HmassT_INPUTS,
        
            ('P', 'Q'): CP.PQ_INPUTS, ('P', 'S'): CP.PSmass_INPUTS, ('P', 'T'): CP.PT_INPUTS, ('P', 'U'): CP.PUmass_INPUTS,
        
            ('Q', 'S'): CP.QSmass_INPUTS, ('Q', 'T'): CP.QT_INPUTS, ('S', 'T'): CP.SmassT_INPUTS, ('S', 'U'): CP.SmassUmass_INPUTS, ('T', 'U'): CP.TUmass_INPUTS,
        
            # --- Molar-based inputs ---
            ('Dmolar', 'Hmolar'): CP.DmolarHmolar_INPUTS, ('Dmolar', 'P'): CP.DmolarP_INPUTS, ('Dmolar', 'Q'): CP.DmolarQ_INPUTS, ('Dmolar', 'Smolar'): CP.DmolarSmolar_INPUTS, ('Dmolar', 'T'): CP.DmolarT_INPUTS, ('Dmolar', 'Umolar'): CP.DmolarUmolar_INPUTS,
        
            ('Hmolar', 'P'): CP.HmolarP_INPUTS, ('Hmolar', 'Q'): CP.HmolarQ_INPUTS, ('Hmolar', 'Smolar'): CP.HmolarSmolar_INPUTS, ('Hmolar', 'T'): CP.HmolarT_INPUTS,
        
            ('P', 'Smolar'): CP.PSmolar_INPUTS, ('P', 'Umolar'): CP.PUmolar_INPUTS, ('Q', 'Smolar'): CP.QSmolar_INPUTS, ('Smolar', 'T'): CP.SmolarT_INPUTS, ('Smolar', 'Umolar'): CP.SmolarUmolar_INPUTS, ('T', 'Umolar'): CP.TUmolar_INPUTS,
            }

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
                warnings.warn("Error: Too many state variables")
            elif (len(self.variables_input) == 2 or (len(self.variables_input) == 1 and (self.SC is not None or self.SH is not None))):
                self.calculate_properties()
            elif len(self.variables_input)<2:
                pass


            if (self.m_dot != None or self.V_dot !=None) and self.state_known:
                if self.m_dot != None:
                    self.V_dot = self.m_dot / self.D * 3600
                elif self.V_dot != None:
                    self.m_dot = self.V_dot * self.D / 3600
                self.completely_known = True
        else:
            pass

    def get_AS_inputs(self, key1, key2):
        """
        Return:
            - The appropriate CoolProp AbstractState INPUT constant.
            - A flag (0 if original order kept, 1 if reversed).
        
        Example:
            self.get_AS_inputs('P', 'T')  ->  (CP.PT_INPUTS, 0)
            self.get_AS_inputs('T', 'P')  ->  (CP.PT_INPUTS, 1)
        """
        # Determine if we reversed key order during sorting
        pair_sorted = tuple(sorted([key1, key2]))
        reversed_flag = int(pair_sorted != (key1, key2))  # 1 if reversed
    
        if pair_sorted not in self.CP_map:
            raise ValueError(f"Unsupported input pair: {pair_sorted}. "
                             f"Valid keys: {list(self.CP_map.keys())}")
    
        return self.CP_map[pair_sorted], reversed_flag
    
    def _replace_state_variable(self, key, value):
        """
        Replace or insert a state variable in variables_input.
        Ensures only valid CoolProp inputs are used.
        """
        for i, var in enumerate(self.variables_input):
            if var[0] == key:
                self.variables_input[i][1] = value
                return
        self.variables_input.append([key, value])


    def calculate_properties(self):
        # 1) Resolve superheating/subcooling into T or P
        try:
            if self.SH is not None and self.SC is not None:
                raise ValueError("Cannot specify both SH and SC simultaneously.")

            # ---- SUPERHEATING ----
            if self.SH is not None:
                if self.p is not None:
                    # p known -> compute T
                    self.AS.update(CP.PQ_INPUTS, self.p, 1) # Saturation temp at given p
                    T_sat = self.AS.T()
                    self.T = T_sat + self.SH
                elif self.T is not None:
                    # T known -> compute p
                    self.AS.update(CP.QT_INPUTS, 1, self.T - self.SH) # Saturation p at given T
                    self.p = self.AS.p()
                else:
                    pass  # Cannot resolve SH without T or p
            # ---- SUBCOOLING ----
            elif self.SC is not None:
                if self.p is not None:
                    # p known -> compute T
                    self.AS.update(CP.PQ_INPUTS, self.p, 0) # Saturation temp at given p
                    T_sat = self.AS.T()
                    self.T = T_sat - self.SC
                elif self.T is not None:
                    # T known -> compute p
                    self.AS.update(CP.QT_INPUTS, 0, self.T + self.SC) # Saturation p at given T
                    self.p = self.AS.p()
                else:
                    pass  # Cannot resolve SC without T or p

        except Exception as e:
            warnings.warn(f"SH/SC resolution failed: {e}")
            return
        
        # --- After resolving SH / SC ---
        if self.SH is not None or self.SC is not None:

            # Ensure variables_input contains T or P, not SC/SH
            if self.T is not None:
                self._replace_state_variable('T', self.T)

            if self.p is not None:
                self._replace_state_variable('P', self.p)

            # Remove any illegal keys just in case
            self.variables_input = [
                v for v in self.variables_input if v[0] in ['T', 'P', 'H', 'S', 'D', 'Q']
            ]

        
        # 2) Calculate properties based on two known variables
        AS_inputs, reversed_flag = self.get_AS_inputs(self.variables_input[0][0], self.variables_input[1][0])
             
        try:
            if not reversed_flag:
                self.AS.update(AS_inputs, self.variables_input[0][1], self.variables_input[1][1])
            else:
                self.AS.update(AS_inputs, self.variables_input[1][1], self.variables_input[0][1])
                            
            self.T = self.AS.T()
            self.p = self.AS.p()
            self.h = self.AS.hmass()
            self.s = self.AS.smass()
            self.D = self.AS.rhomass()
            self.cp = self.AS.cpmass()
            self.x = self.AS.Q()
            self.state_known = True
            
        except:
            warnings.warn("Error: This pair of inputs is not yet supported.")
            
    def set_properties(self, **kwargs):
        
        if 'fluid' in kwargs:
            self.set_fluid(kwargs['fluid'])
        
        for key, value in kwargs.items():
            if key.lower() == 'fluid':
                self.set_fluid(value)
            elif key == 'm_dot':
                self.set_m_dot(value)
            elif key == 'V_dot':
                self.set_V_dot(value)
            elif key.upper() == 'T':
                self.set_T(value)
            elif key.upper() == 'P':
                self.set_p(value)
            elif key.upper() == 'H':
                self.set_h(value)
            elif key.upper() == 'S':
                self.set_s(value)
            elif key.upper() == 'D':
                self.set_D(value)
            elif key.lower() == 'x':
                self.set_x(value)
            elif key.lower() == 'cp':
                self.set_cp(value)
            elif key.upper() == 'SC':
                self.set_SC(value)
            elif key.upper() == 'SH':
                self.set_SH(value)
            else:
                warnings.warn(f"Error: Invalid property '{key}'")
        
    def set_fluid(self, value):
        if self.fluid != None:
            pass
        elif self.fluid == None:
            try:
                try: 
                    self.AS = CP.AbstractState("BICUBIC&HEOS", value) 
                    self.fluid = value
                except:
                    if "INCOMP" in value:
                        value = value.split(":")[-1]
                    self.AS = CP.AbstractState("INCOMP", value)
                    self.fluid = value
            except:
                warnings.warn("Error: Incorrect fluid name:", value)
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
                    
                    if i == 1:
                        self.variables_input[0], self.variables_input[1] = self.variables_input[1], self.variables_input[0]

                    self.check_completely_known()
                    return
                
            self.variables_input.pop()
            self.variables_input.insert(0, ['T', value])
            self.check_completely_known()
            return
                
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
                    
                    if i == 1:
                        self.variables_input[0], self.variables_input[1] = self.variables_input[1], self.variables_input[0]
                    
                    self.check_completely_known()                        
                    return
            
            self.variables_input.pop()
            self.variables_input.insert(0, ['P', value])
            self.check_completely_known()       
            return
            
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
                    
                    if i == 1:
                        self.variables_input[0], self.variables_input[1] = self.variables_input[1], self.variables_input[0]

                    self.check_completely_known()                        
                    return
            
            self.variables_input.pop()
            self.variables_input.insert(0, ['H', value])
            self.check_completely_known()

            return
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
                    
                    if i == 1:
                        self.variables_input[0], self.variables_input[1] = self.variables_input[1], self.variables_input[0]
                        
                    self.check_completely_known()
                    return

            self.variables_input.pop()
            self.variables_input.insert(0, ['S', value])
            self.check_completely_known()
            return                
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
                    
                    if i == 1:
                        self.variables_input[0], self.variables_input[1] = self.variables_input[1], self.variables_input[0]
                    
                    self.check_completely_known()
                    return

            self.variables_input.pop()
            self.variables_input.insert(0, ['D', value])
            self.check_completely_known()
            
            return
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
                    
                    if i == 1:
                        self.variables_input[0], self.variables_input[1] = self.variables_input[1], self.variables_input[0]
                    
                    self.check_completely_known()
                    return

            self.variables_input.pop()
            self.variables_input.insert(0, ['Q', value])
            self.check_completely_known()

            return

        else:          # If the quality is not known, set the value and add the variable to the list
            self.x = value
            self.variables_input = self.variables_input+[['Q',value]]
        self.check_completely_known()

    def set_cp(self, value):
        if self.cp != None: # If the cp is already known, update the value and the corresponding variable in the list
            self.cp = value
        else:          # If the quality is not known, set the value and add the variable to the list
            self.cp = value # We don't want to calculate in this case

    def set_SC(self, value):
        self.SC = value
        self.SH = None  # enforce exclusivity
        self.check_completely_known()

    def set_SH(self, value):
        self.SH = value
        self.SC = None  # enforce exclusivity
        self.check_completely_known()

            
    def print_resume(self, unit_T='K', unit_p='Pa'):
        """
        Parameters
        ----------
        unit_T = Temperature unit: 'K' or 'C'
        unit_p = Temperature unit: 'Pa' or 'bar'
        
        """
        if self.fluid is not None:
            print("Fluid: " + self.fluid + "")
        else:
            print("Fluid: None")            
            
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

