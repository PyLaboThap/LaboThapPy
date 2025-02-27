# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 11:21:49 2021

@author: jvega
"""


def kim_dry_out_incipience(G, q, Dh, P_star, rho_l, rho_v, mu_l, sigma, i_fg):
    """
    Inputs
    ------
    ?
    
    Outputs
    -------
    ?
    
    Reference
    ---------
    ?
    
    """
    
    Bo = q/G/i_fg
    We_lo = (Dh*G**2)/(rho_l*sigma)
    Ca = (mu_l*G)/(rho_l*sigma)
    x_di = 1.4*(We_lo**0.03)*(P_star**0.08) - 15*(Bo**0.15)*(Ca**0.35)*(rho_v/rho_l)**0.06
    return x_di