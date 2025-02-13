# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:31:35 2025

@author: Alanis Zeoli
"""

# /!\ Also implement an easy way to update values (iteration)

from component.base_component import BaseComponent
from channel import Channel

from CoolProp.CoolProp import HAPropsSI
    
class DEC_KPI(BaseComponent):
    
    def __init__(self):
        super().__init__() 
        self.channels = []
        
    def add_output(self,channel_number,T_out,w_out,P_out=101325):
        self.channels[channel_number].add_output(T_out,w_out,P_out)
                
    def display_channels(self):
        for channel in self.channels:
            print(f"\nChannel {channel.name}:")
            print(channel)
    
    # ------------------------------------------------------------------------
    def add_input(self,name,**inputs):
        """Create a new channel and update connector properties."""          
        if 'channel_type' not in inputs:
            inputs['channel_type'] = 'wet'
        new_channel = Channel(name,**inputs)
        self.channels.append(new_channel)   

    def get_required_parameters(self):
        return ['epsilon_wb']
        
    def solve(self):
        # Definition of inputs
        channel_number = 0 # There is only one channel layer in a DEC
        T_in = self.channels[channel_number].su.T
        w_in = self.channels[channel_number].su.w
        P_in = self.channels[channel_number].su.P
        
        # Definition of parameters
        epsilon_wb = self.params['epsilon_wb']
        
        # Output temperature computation
        T_lim = HAPropsSI('B','T',T_in,'W',w_in,'P',P_in)
        T_out = T_in - epsilon_wb*(T_in-T_lim)
        
        # Output humidity computation
        w_max = HAPropsSI('W','T',T_lim,'RH',1,'P',P_in)
        w_out = w_in + epsilon_wb*(w_max-w_in)
        
        # No head losses
        P_out = P_in
        
        # Add output
        self.add_output(channel_number,T_out,w_out,P_out)
        
        
DEC = DEC_KPI()
DEC.add_input("p", T=22+273.15, w=0.012)
DEC.set_parameters(epsilon_wb=0.85)
DEC.solve()

DEC.display_channels()