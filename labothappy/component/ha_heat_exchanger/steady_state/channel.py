# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:57:09 2025

@author: Alanis Zeoli
"""

from connector.humid_air_connector import HAConnector
from connector.mass_connector import MassConnector

from CoolProp.CoolProp import HAPropsSI
from CoolProp.CoolProp import PropsSI

class Channel:
    def __init__(self,name,**inputs):
        self.name = name
        self.su = HAConnector()
        self.ex = HAConnector()
        
        if 'channel_type' in inputs:
            self.channel_type = inputs.pop('channel_type')
        else:
            self.channel_type = 'wet'
            
        self.su.set_properties(**inputs)
        
        if 'T' not in inputs:
            self.su.get_properties('T')
        if 'w' not in inputs:
            self.su.get_properties('W')        
        
    def add_output(self,T_ex,w_ex,P_ex=101325):
        self.ex.set_properties(T=T_ex,w=w_ex,P=P_ex)
        
    def update_input(self,T_su,w_su,P_su=101325):
        self.su.set_properties(T=T_su,w=w_su,P=P_su)
        
    def __repr__(self):
        T_in = self.su.T-273.15
        T_out = self.ex.T-273.15
        return f"Inlet: T = {T_in}°C, w = {self.su.w}kg/kg \nOutlet: T = {T_out}°C, w = {self.ex.w}kg/kg"