# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 11:24:10 2021

Description:
    This function is an implementation of the find_2P_boundaries written by Rdickes.
    All the functions involved are grouped in the "master" find_2P_boundaries function.
    
    Objective:
        -Find saturation liquid and vapour pressures and enthalpies when pressure drop occurs in a HX.
    
    Important notice:
        1) Only pure fluids are considered
        2) find_saturated_liquid_pressure and find_saturated_vapour_pressure are condensed
        into a single find_saturated_pressure function that uses a selector.

@author: jvega
"""

from scipy.optimize import fminbound
import CoolProp.CoolProp as CP

#%%
def find_2P_boundaries(fluid, h_su, h_ex, P_su, P_ex):
    #-------------------------------------------------------------------------
    def find_saturated_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, sat):
        """This function has a selector to switch between saturated liquid or vapor"""
        P_guess = x*ub#[a*ub for a in x]
        #---------------------------------------------------------------------
        if sat == "liq":
            h_sat = CP.PropsSI("H", "P", P_guess, "Q", 0, fluid)
        elif sat == "vap":
            h_sat = CP.PropsSI("H", "P", P_guess, "Q", 1, fluid)
        else:
            raise ValueError("Typo when assigning the saturated condition")
        #---------------------------------------------------------------------   
        ratio_h =(h_su - h_sat)/max(0.01,(h_su - h_ex))
        P_bis = (1-ratio_h)*P_su + ratio_h*P_ex
        res = abs(P_guess - P_bis)/P_bis
        return res
    #-------------------------------------------------------------------------
    def res_find_saturated_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, sat):
        res = find_saturated_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, sat)
        return res
    #-------------------------------------------------------------------------
    
    if abs(P_su-P_ex)<1e1:
        P_l = 0.5*P_su + 0.5*P_ex
        P_v = 0.5*P_su + 0.5*P_ex
        flag_l_bd = 1
        flag_v_bd = 1
    else:
        lb = 0.999*min(P_ex, P_su)
        ub = 1.001*max(P_ex, P_su)
        #----------------------------------------------------------------------
        f_l = lambda x: res_find_saturated_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, "liq")
        out_fminbnd = fminbound(f_l,lb/ub ,ub/ub, full_output=1)
        x_P_l, _, flag_l_bd = out_fminbnd[0],  out_fminbnd[1], abs(out_fminbnd[2]-1)
        P_l = ub*x_P_l #[a*ub for a in x_P_l]
        #----------------------------------------------------------------------
        f_v = lambda x: res_find_saturated_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, "vap")
        out_fminbnd =  fminbound(f_v, lb/ub, ub/ub, full_output=1);
        x_P_v, _, flag_v_bd = out_fminbnd[0],  out_fminbnd[1], abs(out_fminbnd[2]-1)
        P_v = ub*x_P_v #[a*ub for a in x_P_v]
        #----------------------------------------------------------------------
    h_l = CP.PropsSI("H", "P", P_l, "Q", 0, fluid)
    h_v = CP.PropsSI("H", "P", P_v, "Q", 1, fluid)
    #----------------------------------------------------------------------
    return h_l, h_v, P_l, P_v, flag_l_bd, flag_v_bd