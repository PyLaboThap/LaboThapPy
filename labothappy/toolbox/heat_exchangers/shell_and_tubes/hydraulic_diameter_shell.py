# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 18:05:51 2026

@author: titouanjanod
"""

import math

"""From 
A. M. Flynn, Kern’s Process Heat Transfer, 2nd ed. Somerset: John Wiley & Sons, Incorporated, 2019.
Available: https://udghoshna.wordpress.com/wp-content/uploads/2013/06/35306996-process-heat-transfer-donald-q-kern1.pdf
Page 139
"""

def hydraulic_diameter(pitch, OD, tubes_pattern):
    """
    Computes the shell "side hydraulic diameter" in a shell and tube heat exchanger, 
    also depicted as "equivalent diameter" in literature
    
    
    Inputs:
        pitch: distance between the center of each tube of the bundle [m]
        OD: Outside diameter of the tubes of the bundle [m]
        
    Parameter:
        pattern: organisation of the tubes, either 'square' or 'triangle'
        
    Output:
        D_h_shell: Shell-side hydraulic diameter
    """
    if tubes_pattern == "square":
        D_h_shell = 4*(pitch**2 - math.pi * OD**2/4)/ (math.pi * OD)
    
    
    elif tubes_pattern == "triangle":
        D_h_shell = 4*(0.86*pitch**2 - math.pi * OD**2/4)/ (math.pi * OD)
    
    return D_h_shell