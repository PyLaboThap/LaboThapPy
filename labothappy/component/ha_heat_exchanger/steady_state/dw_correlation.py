# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:01:38 2025

@author: Alanis Zeoli

Solving a DW based on correlations developed by Panaras et al. (2010)
"""

from component.base_component import BaseComponent
from channel import Channel

from CoolProp.CoolProp import HAPropsSI
import numpy as np
from scipy.optimize import fsolve

class DW_corr(BaseComponent):
    
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
        return ['eta_F1','eta_F2']
        
    def solve(self):
        # Definition of inputs
        channel_pro = min(range(len(self.channels)), key=lambda i: self.channels[i].su.T)
        T_p_in = self.channels[channel_pro].su.T
        w_p_in = self.channels[channel_pro].su.w
        P_p_in = self.channels[channel_pro].su.P
        
        channel_reg = max(range(len(self.channels)), key=lambda i: self.channels[i].su.T)
        T_r_in = self.channels[channel_reg].su.T
        w_r_in = self.channels[channel_reg].su.w
        # P_r_in = self.channels[channel_reg].su.P
        
        
        # Definition of parameters
        eta_F1 = self.params['eta_F1']
        eta_F2 = self.params['eta_F2']
        
        # Output temperature computation
        F1_p = -2865/(T_p_in**1.49)+4.344*w_p_in**0.8624;
        F2_p = (T_p_in**1.49)/6360-1.127*w_p_in**0.07969;

        F1_r = -2865/(T_r_in**1.49)+4.344*w_r_in**0.8624;
        F2_r = (T_r_in**1.49)/6360-1.127*w_r_in**0.07969;

        F1_2 = F1_p + eta_F1*(F1_r - F1_p);
        F2_2 = F2_p + eta_F2*(F2_r - F2_p);
        
        
        def eqn(w_p_out,F1_2,F2_2):            
            return F1_2*F2_2 + F1_2*1.127*w_p_out**0.07969 - F2_2*4.344*w_p_out**0.8624 - 4.896*w_p_out**0.9421 + 0.45
        
        # Initial guess
        w_p_out_guess = w_p_in-0.003
        w_p_out = fsolve(eqn,w_p_out_guess,args=(F1_2,F2_2));
    
        F2_2 = F2_p + eta_F2*(F2_r - F2_p);
        T_p_out = (6360*(F2_2+1.127*w_p_out**0.07969))**(1/1.49);
        
        # No head losses
        P_p_out = P_p_in
        # P_r_out = P_r_in
        
        # Add output
        self.add_output(channel_pro,T_p_out,w_p_out,P_p_out)
        self.add_output(channel_reg,T_p_out,w_p_out,P_p_out)
        
        
DW = DW_corr()
DW.add_input("pro", T=32+273.15, w=0.014)
DW.add_input("reg", T=70+273.15, w=0.01)
DW.set_parameters(eta_F1=0.08,eta_F2=0.6)
DW.solve()

DW.display_channels()