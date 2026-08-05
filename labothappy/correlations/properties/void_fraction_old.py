# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 14:22:15 2024

@author: Basile
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

"""
Source : 

John R. Thome and Andrea Cioncolini, 2015. Encyclopedia of Two-Phase Heat Transfer and Flow, Set 1: 
Fundamentals and Methods, Volume 3: Flow Boiling in Macro and Microchannels, World Scientific 
Publishing - Chapter 4 : Void Fraction
    
"""

#%%

def generate_boolean_vector(length, ranges):
    boolean_vector = [np.nan] * length  # Initialize the vector with 0s

    for start, end in ranges:
        boolean_vector[start:end+1] = [1] * (end - start + 1)

    return boolean_vector

def rms_between_arrays(arr1, arr2):
    # Check if the arrays have the same length
    if len(arr1) != len(arr2):
        raise ValueError("Arrays must have the same length")

    # Calculate the squared differences
    squared_diff = (arr1 - arr2) ** 2

    # Calculate the mean of squared differences
    mean_squared_diff = np.mean(squared_diff)

    # Calculate the RMS value
    rms_value = np.sqrt(mean_squared_diff)

    return rms_value

def void_fraction(Q, rho_g, rho_l):
    """
    Inputs
    ----------
    Q : Vapor Quality [-]
    
    rho_g : saturated gas density [kg/m^3]
        
    rho_l : saturated liquid density [kg/m^3]

    Outputs
    -------
    eps : Void fraction [-] : Actual volumetric fraction of gas phase and void in a 2 phase flow

    """
    
    
    "Test of Cioncolini and Thome (2012)"
    # Annular 2P flow 
    # Better for 0 < x < 1, 1e-3 < rho_g/rho_l < 1, 0.7 < eps < 1
    # Here : rho_g/rho_l = 0.3 at 36 bars
    
    h_CT = -2.129 + 3.129*(rho_g/rho_l)**(-0.2186)
    n_CT = 0.3487 + 0.6513*(rho_g/rho_l)**(0.515)
    
    eps = (h_CT*Q**n_CT)/(1+(h_CT-1)*Q**n_CT)    
    
    # if eps_CT < 0.7:
    #     G = m_dot/A_channel
        
    #     J_l = (1-Q)*G/rho_l # Liquid phase superficial velocities 
    #     J_g = Q*G/rho_g # Gas phase superficial velocities 
        
    #     "Rouhani and Axelsson (1970)"
    #     C_0 = 1 + 0.2*(1-Q)
    #     V_drift = 1.18*((g*sigma*(rho_l-rho_g))/(rho_l**2))**(0.25)
        
    #     eps_RA = J_g/(C_0*(J_l + J_g) + V_drift)
        
    #     eps = eps_RA
    # else: 
    #     eps = eps_CT
    
    rho_2 = rho_g*eps + rho_l*(1-eps)
    
    return eps, rho_2

test_rho = 0

if test_rho == 1:
    Q = np.linspace(0,1,100) # Vap quality
    rho = np.zeros(len(Q))
    rho_corr = np.zeros(len(Q))
    eps = np.zeros(len(Q))
    P_sat = 0.8*1e5
    
    rho_g = PropsSI('D','P',P_sat,'Q',1,'Cyclopentane')
    rho_l = PropsSI('D','P',P_sat,'Q',0,'Cyclopentane')
    
    for i in range(len(Q)):
        rho[i] = PropsSI('D','P',0.8*1e5,'Q',Q[i],'Cyclopentane')
        eps[i],rho_corr[i] = void_fraction(Q[i],rho_g,rho_l)

    plt.plot(Q, rho, label = 'Density')
    plt.plot(Q, rho_corr, label = 'Corrected Density')
    plt.legend()
    plt.grid()
    
    plt.show()

    Q = np.linspace(0,1,100) # Vap quality
    rho = np.zeros(len(Q))
    rho_corr = np.zeros(len(Q))
    eps = np.zeros(len(Q))
    P_sat = 35*1e5
    
    rho_g = PropsSI('D','P',P_sat,'Q',1,'Cyclopentane')
    rho_l = PropsSI('D','P',P_sat,'Q',0,'Cyclopentane')
    
    for i in range(len(Q)):
        rho[i] = PropsSI('D','P',P_sat,'Q',Q[i],'Cyclopentane')
        eps[i],rho_corr[i] = void_fraction(Q[i],rho_g,rho_l)

    plt.plot(Q, rho, label = 'Density')
    plt.plot(Q, rho_corr, label = 'Corrected Density')
    plt.legend()
    plt.grid()
    
    plt.show()

#%%

test = 0
if test == 1:

    wf = "Cyclopentane"
    P_vect = [34*1e5,0.8*1e5] # Sat pressure
    
    rms_SRH = 0
    rms_Fauske = 0
    rms_Moody = 0
    rms_AT = 0
    rms_Bankoff = 0
    rms_RA = 0
    rms_DIX = 0
    rms_DIX_WG = 0
    rms_LM = 0
    rms_CT = 0
    
    for P in P_vect:
        Q = np.linspace(0,1,100) # Vap quality
        m_dot = 0.014 # kg/s
        g = 9.81 # m/s^2 
        
        sigma = PropsSI('I','P',P,'Q',Q,wf) # Surface Tension 
        rho = PropsSI('D','P',P,'Q',Q,wf) # Surface Tension 
        T = PropsSI('T','P',P,'Q',Q,wf) # Sat T
        
        (rho_l,mu_l) = PropsSI(('D','L'),'P',P,'Q',0,wf) # Density and viscosity
        
        (rho_g,mu_g) = PropsSI(('D','L'),'P',P,'Q',1,wf) # Density and viscosity
        
        #%%
        "1) Slip ratio Methods"
        
        v_g = 0.2 # gas phase velocity : m/s # !!! : arbitrary
        v_l = 0.1 # liquid phase velocity : m/s # !!! : arbitrary
        
        S_gl = v_g/v_l
        
        eps_SR = (1 + S_gl*((1-Q)/Q)*(rho_g/rho_l))**(-1)
        eps_SR_plot = eps_SR
        
        "1.1) Homogeneous model"
        # v_g = v_l 
        # Good when flow is 1 phase
        #x_valid_SRH = 
        
        ranges_SRH = [(0, 20), (80, 99)]
        x_valid_SRH = generate_boolean_vector(len(Q), ranges_SRH)
        
        eps_SRH = (1 + ((1-Q)/Q)*(rho_g/rho_l))**(-1)
        eps_SRH_plot = eps_SRH*x_valid_SRH
        
        "1.2) Fauske Method (1962)"
        eps_SR_Fauske = (1 + ((1-Q)/Q)*(rho_g/rho_l)**(1/2))**(-1)
        eps_SR_Fauske_plot = eps_SR_Fauske*x_valid_SRH
        
        "1.3) Moody Method (1966)"
        eps_SR_Moody = (1 + ((1-Q)/Q)*(rho_g/rho_l)**(2/3))**(-1)
        eps_SR_Moody_plot = eps_SR_Moody*x_valid_SRH
        
        #%%
        "2) K eps_h methods"
        
        ranges_K_eps_h = [(20,80)]
        x_valid_K_eps_h = generate_boolean_vector(len(Q), ranges_K_eps_h)
        
        "2.1) Armand and Treschev (1947)"
        K_AT = 0.833 + 0.167*Q
        eps_AT = K_AT*eps_SRH
        eps_AT_plot = eps_AT*x_valid_K_eps_h
        
        "2.2) Bankoff (1960)"
        K_Bankoff = 0.71 + 0.0145*(P/1e6)
        eps_Bankoff = K_Bankoff*eps_SRH
        eps_Bankoff_plot = eps_Bankoff*x_valid_K_eps_h
        
        #%%
        "3) Drift_flux Correlations"
        
        D = 0.001 # Channel Diameter
        W = 0.1 # Channel Width
        A_channel = D*W
        
        G = m_dot/A_channel
        
        J_l = (1-Q)*G/rho_l # Liquid phase superficial velocities 
        J_g = Q*G/rho_g # Gas phase superficial velocities 
        
        ranges_DF = [(0,20)]
        x_valid_DF = generate_boolean_vector(len(Q), ranges_DF)
        
        "3.1) Rouhani and Axelsson (1970)"
        C_0_RA = 1 + 0.2*(1-Q)
        V_drift_RA = 1.18*((g*sigma*(rho_l-rho_g))/(rho_l**2))**(0.25)
        
        eps_DF_RA = J_g/(C_0_RA*(J_l + J_g) + V_drift_RA)
        eps_DF_RA_plot = eps_DF_RA*x_valid_DF
        
        "3.2) DIX model, Chexal et al. (1989)"
        theta = np.pi/2 # Vertical flow inclination
        
        n_DIX = (rho_g/rho_l)**0.1
        C_0_DIX = J_g/(J_l + J_g)*(1+(J_l/J_g)**n_DIX) 
        V_drift_DIX = 2.9*((g*sigma*(rho_l-rho_g))/(rho_l**2))**(0.25)
        
        eps_DF_DIX = J_g/(C_0_DIX*(J_l + J_g) + V_drift_DIX)
        eps_DF_DIX_plot = eps_DF_DIX*x_valid_DF
        
        "3.3) Woldesemayat and Ghajar (2007)"
        
        V_drift_DIX_WG = 2.9*((g*sigma*D*(1 + np.cos(theta))*(rho_l-rho_g))/(rho_l**2))**(0.25)*(1.22 + 1.22*np.sin(theta))**(101325/P)
        eps_DF_DIX_WG = J_g/(C_0_DIX*(J_l + J_g) + V_drift_DIX_WG)
        
        eps_DF_DIX_WG_plot = eps_DF_DIX_WG*x_valid_DF
        
        #%%
        "4) Various Empirical Correlations"
        
        "4.1) Lockhart and Martinelli"
        
        L = 0.01
        
        u_g = 1
        u_l = 0.2
        
        # Reynolds numbers
        
        Re_g = (rho_g*u_g*L)/mu_g #
        Re_l = (rho_l*u_l*L)/mu_l #
        
        # Friction coef (Hagen-Poiseuille and Mac Adams)
        
        if Re_g < 1500:
            
            f_g = 16/Re_g
            
        else:
            
            f_g = 0.046/Re_g**0.2
         
        if Re_l < 1500:
            
            f_l = 16/Re_l
            
        else:
            
            f_l = 0.046/Re_l**0.2
        
        # Lockhart and Martinelli function (1949)
        X = (f_l/f_g)**(0.5) * ((1 - Q)/Q) * (rho_l/rho_g)**(0.5)
        
        # Void fraction
        eps_LM = 1/(1 + 0.28*X**(0.71))
        eps_LM_plot = eps_LM
        
        #%%
        
        "4.2) Cioncolini and Thome (2012)"
        # Annular 2P flow 
        # Better for 0 < x < 1, 1e-3 < rho_g/rho_l < 1, 0.7 < eps < 1
        # Here : rho_g/rho_l = 0.3 at 36 bars
        
        h_CT = -2.129 + 3.129*(rho_g/rho_l)**(-0.2186)
        n_CT = 0.3487 + 0.6513*(rho_g/rho_l)**(0.515)
        
        eps_CT = (h_CT*Q**n_CT)/(1+(h_CT-1)*Q**n_CT)
        eps_CT_plot = eps_CT
        
        #%%
        
        "5) Plot everything"
        
        plt.plot(Q,eps_SR_plot,color = 'darkgreen')
        plt.plot(Q,eps_SRH_plot,color = 'limegreen')
        plt.plot(Q,eps_SR_Fauske_plot,color = 'olivedrab')
        plt.plot(Q,eps_SR_Moody_plot,color = 'lightseagreen')
        
        plt.plot(Q,eps_AT_plot,color = 'red')
        # plt.plot(Q,eps_Bankoff_plot,color = 'orangered')
        
        plt.plot(Q,eps_DF_RA_plot,color = 'blue')
        plt.plot(Q,eps_DF_DIX_plot,color = 'deepskyblue')
        plt.plot(Q,eps_DF_DIX_WG_plot,color = 'navy')
        
        # plt.plot(Q,eps_LM_plot,color = 'purple')
        plt.plot(Q,eps_CT_plot,color = 'Orchid')
        
        plt.grid()
        plt.show()
        
        #%%
        
        "6) Make a final curve"
        
        eps_mean = np.zeros(len(Q))
        
        for i in range(len(Q)):
            n_var = 0
        
            if not (np.isnan(eps_SRH_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_SRH_plot[i] 
                n_var = n_var + 1 
        
            if not (np.isnan(eps_SR_Fauske_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_SR_Fauske_plot[i] 
                n_var = n_var + 1 
        
            if not (np.isnan(eps_SR_Moody_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_SR_Moody_plot[i] 
                n_var = n_var + 1 
        
            if not (np.isnan(eps_AT_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_AT_plot[i] 
                n_var = n_var + 1 
        
            # if not (np.isnan(eps_Bankoff_plot[i])):
            #     eps_mean[i] = eps_mean[i] + eps_Bankoff_plot[i] 
            #     n_var = n_var + 1 
        
            if not (np.isnan(eps_DF_RA_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_DF_RA_plot[i] 
                n_var = n_var + 1 
        
            if not (np.isnan(eps_DF_DIX_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_DF_DIX_plot[i] 
                n_var = n_var + 1 
                
            if not (np.isnan(eps_DF_DIX_WG_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_DF_DIX_WG_plot[i] 
                n_var = n_var + 1 
                
            # if not (np.isnan(eps_LM_plot[i])):
            #     eps_mean[i] = eps_mean[i] + eps_LM_plot[i] 
            #     n_var = n_var + 1 
        
            if not (np.isnan(eps_CT_plot[i])):
                eps_mean[i] = eps_mean[i] + eps_CT_plot[i] 
                n_var = n_var + 1 
                
            eps_mean[i] = eps_mean[i]/n_var
        
        plt.plot(Q,eps_mean,'red')
        
        #%%
        
        "7) Find best method"
        
        rms_SRH = rms_SRH + rms_between_arrays(eps_mean, eps_SRH)
        rms_Fauske = rms_Fauske + rms_between_arrays(eps_mean, eps_SR_Fauske)
        rms_Moody = rms_Moody + rms_between_arrays(eps_mean, eps_SR_Moody)
        rms_AT = rms_AT + rms_between_arrays(eps_mean, eps_AT)
        rms_Bankoff = rms_Bankoff + rms_between_arrays(eps_mean, eps_Bankoff)
        rms_RA = rms_RA + rms_between_arrays(eps_mean, eps_DF_RA)
        rms_DIX = rms_DIX + rms_between_arrays(eps_mean[1:-1], eps_DF_DIX[1:-1])
        rms_DIX_WG = rms_DIX_WG + rms_between_arrays(eps_mean[1:-1], eps_DF_DIX_WG[1:-1])
        rms_LM = rms_LM + rms_between_arrays(eps_mean, eps_LM)
        rms_CT = rms_CT + rms_between_arrays(eps_mean, eps_CT)
    
        plt.plot(Q,Q,label = 'Bissectriss')
        plt.plot(Q,eps_mean, label = 'Mean of all methods')
        plt.plot(Q,eps_DF_RA, label = 'Rouhani and Axelsson')
        plt.plot(Q,eps_CT, label = 'Cioncolini and Thome')
        
        plt.legend()
        plt.show()
    
    print("Choice of Rouhani and Axelsson method for drift flux. For annular flows, choose Cioncolini and Thome")