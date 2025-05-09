# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:12:17 2025

@author: Alanis Zeoli
"""

from CoolProp.CoolProp import HAPropsSI
import pandas as pd
import warnings

class HAConnector:
    """
    A class to handle and calculate fluid properties based on given state variables.

    **Attributes**:
    
        m_dot : float, optional
            Mass flow rate [kg/s].
        V_dot : float, optional
            Volume flow rate [m^3/s].
        T : float, optional
            Temperature [K].
        P : float, optional
            Pressure [Pa].
        h : float, optional
            Specific enthalpy [J/kg].
        w : float, optional
            Specific humidity [kg/kg].
        RH: float, optional
            Relative humidity [-].
        v : float, optional
            Specific volume [m^3/kg].
        Tdp : float, optional
            Dew point temperature [K].
        Twb : float, optional
            Wet bulb temperature [K].
        cp : float, optional
            Specific heat capacity in J/kg/K.

    **Methods**:

        __init__(self):
            Initializes the MassConnector object with fluid = humid air.
        
        reset(self):
            Resets all properties to None and flags to False.

        set_properties(self,**kwargs):
            Sets fluid properties based on given inputs.

        get_properties(self,properties):
            Computes asked properties using CoolProp based on known variables.
            'properties' is a list of the properties that should be computed.

        print_resume(self, unit_T='K', unit_p='Pa'):
            Prints a summary of the current properties in specified units.

    **Exceptions**:
    
    ValueError
        Raised if required variables for property calculation are missing.
    """

    def __init__(self,**kwargs):
        # Definition of attributes
        # Add attribute table with name of variable, def of variable and units
        attr_data = {'m_dot':['Mass flow rate','kg/s'],
                     'V_dot':['Volumetric flow rate','m^3/s'],
                     'P':['Pressure','Pa'],
                     'T':['Temperature','K'],
                     'h':['Spec. enthalpy','J/K'],
                     'w':['Humidity ratio','kg/kg'],
                     'RH':['Relative humidity',''],
                     'Tdp':['Dew point','K'],
                     'Twb':['Wet bulb','K'],
                     'v':['Spec. volume','m^3/kg'],
                     'cp':['Spec. heat','J/kg.K']}
        self.attributes = pd.DataFrame(data=attr_data,index=['Name','Units'])
        
        # Initialisation of all attributes to None
        for attr in self.attributes.columns:
            setattr(self,attr,None)
        
        # Definition of properties usable in HAPropsSI
        self.HA_attributes = ['T','h','w','RH','Tdp','Twb','v','cp']
        self.HA_prop = ['T','H','W','RH','D','B','V','CP']
        
        self.unknown_prop = self.HA_prop.copy()
        self.known_prop = []
        
        if len(kwargs):
            self.set_properties(**kwargs)

    # Resets all variables
    def reset(self):
        for prop in self.attributes:
            setattr(self,prop,None)
            
        self.unknown_prop = self.HA_prop.copy()
        self.known_prop = []

    def set_properties(self,**kwargs):
        for key, value in kwargs.items():
           if key in self.attributes:
               attr = key
               if key in self.HA_attributes:
                   prop = self.HA_prop[self.HA_attributes.index(attr)]
               else: 
                   prop = 0
               
           elif key in self.HA_prop:
               prop = key
               attr = self.HA_attributes[self.HA_prop.index(prop)]
           else:
               print(f"Error:{key} is not a valid property")
               return
           
           setattr(self,attr,value)
           
           if prop in self.HA_prop:
               if prop in self.known_prop:
                   # Place the property at the end of the list to use it as a reference to compute the fluid state
                   self.known_prop.remove(prop)
                   self.known_prop.append(prop)
               else:
                   self.unknown_prop.remove(prop)
                   self.known_prop.append(prop)
               
        if getattr(self,'P') is None:
            self.P = 101325         # Default: atmospheric value

    def get_properties(self,properties):
        
        # Properties can be computed only if 3 properties (including presssure) are known
        if len(self.known_prop) > 1:
            # Property 1
            arg1 = self.known_prop[-1]
            attr1 = self.HA_attributes[self.HA_prop.index(arg1)]
            val1 = getattr(self,attr1)
            
            # Property 2
            arg2 = self.known_prop[-2]
            attr2 = self.HA_attributes[self.HA_prop.index(arg2)]
            val2 = getattr(self,attr2)
            
            # Property 3 is always the pressure
            P = self.P
            
            # 1. Computation of non CoolProp properties
            if 'm_dot' in properties:
                if self.V_dot is float:
                    if self.v is None:
                        volume = HAPropsSI('V',arg1,val1,arg2,val2,'P',P)
                        self.set_properties(v=volume)
                        
                    m_dot = self.V_dot/self.v
                    self.set_properties(m_dot=m_dot)     
                    
                else:
                    print("Impossible to compute mass flow rate")
                    
                properties.remove('m_dot')
                    
            elif 'V_dot' in properties:       
                if self.m_dot is not None:
                    if self.v is None:
                        volume = HAPropsSI('V',arg1,val1,arg2,val2,'P',P)
                        self.set_properties(v=volume)
                        
                    V_dot = self.m_dot*self.v
                    self.set_properties(V_dot=V_dot)
                else:
                    print("Impossible to compute volumetric flow rate")
            
                properties.remove('V_dot')
            
            # Set other properties
            props = []
            props_value = []
            
            for arg in properties:
                if arg in self.HA_attributes:
                    prop = self.HA_prop[self.HA_attributes.index(arg)]
                    
                elif arg in self.HA_prop:
                    prop = arg
                    
                else:
                    print(f"Error:{arg} is not a valid property")
                    return
            
                try:    
                    props_value.append(HAPropsSI(prop,arg1,val1,arg2,val2,'P',P))
                except:
                    warnings.warn("Error: Check consistency of input variables.")
                    return
                
                props.append(prop)
                
            props_dict = dict(zip(props,props_value))
            self.set_properties(**props_dict)
             
        else:
            print("Not enough properties known.")
        
            
    def print_resume(self, unit_T='K', unit_P='Pa', unit_w='kg/kg', unit_RH=''):
        """
        Parameters
        ----------
        unit_T = Temperature unit: 'K' or 'C'
        unit_P = Temperature unit: 'Pa' or 'bar'
        unit_w = Temperature unit: 'kg/kg' or 'g/kg'
        unit_RH = Temperature unit: '' or '%' 
        """
        
        for attr in self.attributes.columns:
            name = self.attributes[attr]['Name']
            value = getattr(self,attr)
            if value is not None:
                units = self.attributes[attr]['Units']
                
                match units:
                    case 'Pa':
                        if unit_P == 'bar':
                            units = unit_P
                            value = value/1e5
                            
                    case 'K':
                        if unit_T == 'C':
                            units = unit_T
                            value = value-273.15
                            
                    case 'kg/kg':
                        if unit_w == 'g/kg':
                            units = unit_w
                            value = value*1e3
                            
                    case '':
                        if unit_RH == '%':
                            units = unit_RH
                            value = value*1e2
                
                print(name + " : " + str(format(value,'.2f')) + " " + units)
