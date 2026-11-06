#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: Basile Chaudoir
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def _wavelength_from_phi(phi, a):
    """Corrugation wavelength Lambda [m] of a sinusoidal profile of amplitude a [m] whose area
    enlargement factor is phi (Martin: X = 2 pi a / Lambda,
    phi = 1/6 (1 + sqrt(1 + X^2) + 4 sqrt(1 + X^2/2)))."""
    if not 1.0 < phi < 30.0:
        raise ValueError(f"Enlargement factor phi = {phi:.4g} outside ]1, 30[: check the heat transfer "
                         "area, the number of plates and the projected plate dimensions.")
    f = lambda X: (1 + np.sqrt(1 + X**2) + 4*np.sqrt(1 + X**2 / 2)) / 6 - phi
    X = brentq(f, 1e-9, 50.0)
    return 2*np.pi*a / X

class PlateGeomSWEP(object):
    
    # Plate data needed to compute the derived quantities
    _RAW_KEYS = ('A_h', 'h', 'l_v', 'w', 'n_plates', 't_plate', 'n_channel_c', 'n_channel_h')
    
    def __init__(self, **kwargs):

        "Heat transfer area"
        self.A_h = None # Hot side
        self.A_c = None # Cold side
        
        "HX Volume"
        self.H_V_tot = None # Hot side volume
        self.C_V_tot = None # Cold side volume
        
        "HX Dimensions"
        self.l = None # length
        self.w = None # width
        self.h = None # height
        self.l_v = None # length between ports
        self.w_v = None # width between ports
        
        "Plate related parameters"
        self.n_plates = None # number of plates
        self.t_plate = None # plate thickness
        self.plate_cond = None # Plate conduction
        self.fooling = None # Fooling Factor
        self.alpha_ch = None # chevron angle
        
        "Number of channels"
        self.n_channel_c = None # cold side
        self.n_channel_h = None # hot side
        self.n_channels = None # total
        
        "Corrugation (computed from the heat transfer area)"
        self.b = None # channel gap
        self.a = None # corrugation amplitude
        self.amplitude = None # same as a, kept for the scripts that still pass 'amplitude'
        self.phi = None # enlargement factor
        self.Lambda = None # corrugated pitch
        self.X = None # dimensionless corrugation parameter
        
        "Hydraulic diameter and flow area of each side"
        self.C_Dh = None
        self.H_Dh = None
        self.C_A_flow = None
        self.H_A_flow = None

    def _derived_geometry(self):
        """Quantities derived from the plate data: channel gap, corrugation amplitude, enlargement
        factor and corrugated pitch (from the heat transfer area), hydraulic diameter, flow areas."""
        self.n_channels = self.n_channel_h + self.n_channel_c
        
        # Channel gap and corrugation amplitude
        self.b = (self.h - self.n_plates * self.t_plate) / self.n_channels # channel gap
        self.a = self.b/2 # corrugation amplitude
        self.amplitude = self.a
        
        # Enlargement factor and corrugation pitch, computed from the heat transfer area
        n_active = self.n_plates - 2 # plates exchanging heat (the two end plates do not)
        A_projected = n_active * self.l_v * self.w # TODO: projected area as defined by the manufacturer
        self.phi = self.A_h / A_projected # enlargement factor
        self.Lambda = _wavelength_from_phi(self.phi, self.a) # corrugated pitch (sinusoidal profile)
        self.X = 2*np.pi*self.a/self.Lambda # dimensionless corrugation parameter
        
        # Dh : hydraulic diameter
        self.C_Dh = 2*self.b
        self.H_Dh = 2*self.b
        
        # Canal Surface
        self.C_A_flow = self.b*self.w*self.n_channel_c
        self.H_A_flow = self.b*self.w*self.n_channel_h

    def set_parameters(self, name, **kwargs):
        
        if name == "B35TM0x10/1P":
            
            # Heat transfer area  
            self.A_h = 0.752 # m^2
            self.A_c = 0.752 # m^2
            
            # HTX dimensions
            self.l = 0.393 # [m] : length
            self.w = 0.243 # [m] : width
            self.h = 0.0446 # [m] : height
            self.l_v = 0.324 # [m] : length between ports
            
            # Number and thickness of plates
            self.n_plates = 10
            self.t_plate = 0.0008 # [m] # !!! arbitrary
            
            # Fooling factor
            self.fooling = 98.51/100
        
            # Number of canals
            self.n_channel_c = 5
            self.n_channel_h = 4
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction
            self.alpha_ch = 30*np.pi/180 # !!! arbitrary
        
            # Total volume of each part
            self.H_V_tot = 0.9*1e-3 # [m^3]
            self.C_V_tot = 0.9*1e-3 # [m^3]
        
        elif name == "B20Hx24/1P":
            
            # Heat transfer area  
            self.A_h = 0.616 # m^2
            self.A_c = 0.616 # m^2
            
            # HTX dimensions
            self.l = 0.324 # [m] : length
            self.w = 0.094 # [m] : width
            self.h = 0.041 # [m] : height
            self.l_v = 0.2682 # [m] : length between ports
            # self.w_v = 0.094 # [m] : width between ports
            
            # Number and thickness of plates
            self.n_plates = 24
            self.t_plate = 0.0005 # [m] # !!! arbitrary
            
            # Fooling factor
            self.fooling = 0.705/100
        
            # Number of canals
            self.n_channel_c = 12
            self.n_channel_h = 11
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction
            self.alpha_ch = 20*np.pi/180 # !!! arbitrary
        
            # Total volume of each part
            self.H_V_tot = 0.34*1e-3 # [m^3]
            self.C_V_tot = 0.34*1e-3 # [m^3]        

        elif name == "P200THx140/1P_Evaporator":
            
            # POUR EVAPORATEUR
            # Heat transfer area  
            self.A_h = 17.8 # m^2 area hot side
            self.A_c = 17.8 # m^2 area cold side
            
            # HTX dimensions
            self.l = 0.525 # [m] :  length
            self.w = 0.243 # [m] :  width
            self.h = 0.3426 # [m] : height F??
            self.l_v = 0.450 # [m] : length between ports H??
            self.w_v = 0.171 # [m] : width between ports H??
            
            # Number and thickness of plates
            self.n_plates = 140 # OK
            self.t_plate = 0.0003 #0.0019 # [m] # !!! arbitrary épaisseur plaque (certainement plus grand)

            # Fooling factor
            self.fooling = 0.111/100 # [m^2, °C/kW] OK
        
            # Number of canals per pass
            self.n_channel_c = 69 # Cold side OK
            self.n_channel_h = 70 # Hot side OK
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction acier inoxidable 
            self.alpha_ch = 30*np.pi/180 # !!! arbitrary 20, 30 ou 45°
        
            # Total volume of each part
            self.H_V_tot = 16.63*1e-3 # [m^3] OK hot volume
            self.C_V_tot = 16.87*1e-3 # [m^3] OK cold volume  
        
        elif name == "P200THx140/1P_Condenser":
            
            # Heat transfer area  
            self.A_h = 17.8 # m^2 hot side OK
            self.A_c = 17.8 # m^2 cold side OK
            
            # HTX dimensions
            self.l = 0.525 # [m] :  OK
            self.w = 0.243 # [m] :  OK
            self.h = 0.346 # [m] : height F??
            self.l_v = 0.450 # [m] : length between ports H??
            self.w_v = 0.171 # [m] : width between ports

            # Number and thickness of plates
            self.n_plates = 140 # OK
            # ORC MODE:
            self.t_plate = 0.0001 # [m] # !!! arbitrary épaisseur plaque (certainement plus grand)
            # HP MODE:
            # self.t_plate = 0.001 # [m] # !!! arbitrary épaisseur plaque (certainement plus grand)

            # Fooling factor
            self.fooling = 0.258/100 # [m^2, °C/kW] OK
        
            # Number of canals per pass
            self.n_channel_c = 70 # Cold side OK
            self.n_channel_h = 69 # Hot side OK
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction acier inoxidable OK
            # ORC MODE:
            self.alpha_ch = 60*np.pi/180 # !!! arbitrary 20, 30 ou 45°
            # HP MODE:
            # self.alpha_ch = 30*np.pi/180 # !!! arbitrary 20, 30 ou 45°

            # Total volume of each part
            self.H_V_tot = 16.63*1e-3 # [m^3] OK hot volume
            self.C_V_tot = 16.63*1e-3 # [m^3] OK cold volume

        # Manual setting, or modification of the values above
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")

        # Derived quantities (corrugation, hydraulic diameter, flow areas), once the plate data is complete
        if all(getattr(self, key) is not None for key in self._RAW_KEYS):
            self._derived_geometry()