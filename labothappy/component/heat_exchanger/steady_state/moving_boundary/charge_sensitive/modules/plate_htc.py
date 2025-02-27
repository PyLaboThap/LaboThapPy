# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:47:52 2023

@author: Basile
"""

import numpy as np
from scipy.optimize import fsolve

from CoolProp.CoolProp import PropsSI

# For one phase water heat transfer:
def water_plate_HTC(mu, Pr, k, G, Dh):
    """
    Calibrated heat transfer coefficient correlation for water side
    
    Inputs
    ----------
    mu : Viscosity [kg/(m*s)]
    
    Pr : Prandtl Number [/]
    
    k : thermal conductivity [W/(m*K)]
        
    G : Mass flux [kg/(m^2 * s)]
    
    Dh : Spacing between plates [m]

    Outputs
    -------
    h_conv : HTC in convection
    
    Reference
    -------
    Refrigerant R134a vaporisation heat transfer and pressure drop inside a small brazed plate heat exchanger
    G.A. Longo, A. Gasparella

    """
    
    # Bounds on Re (arbitrary) # !!! shall be checked
    Re_max = 1e6
    Re_min = 5
    # Reynolds number
    Re = G*Dh/mu
    if Re <= Re_max and Re >= Re_min:
        # Nusselt number
        Nu = 0.277*Re**(0.766)*Pr**(0.333)
    else: 
        print("Reynolds number for water out of bounds to use the water plate correlation.")
        return 0
        
    # HTC
    h_conv = (k/Dh)*Nu

    return h_conv

# For the R1233zd(E) one phase:
def martin_holger_plate_HTC(mu, Pr, k, m_dot, nb_channels, T_mean, P_mean, fluid, D_h, length, width, amplitude, chevron_angle):
    """
    Reference:
    -----------
    H. M. Hofmann, M. Kind, and H. Martin, ‘Measurements on steady state heat transfer and flow structure and new 
    correlations for heat and mass transfer in submerged impinging jets’, International Journal of Heat and Mass Transfer, 
    
    Correlation from the VDI Atlas
    """

    " Thermodynamics properties:"
    # mu =  PropsSI('V','P', P_mean,'T',T_mean, refrigerant)          # Viscosity, Pa s
    rho = PropsSI('D','P', P_mean,'T',T_mean, fluid)          # Density, kg m^-3
    # cp = PropsSI('C','P', P_mean,'T',T_mean, refrigerant)           # Specific Heat, J/kg-K
    
    

    " Mass flow rate per plate"
    m_dot_ch = m_dot/nb_channels                                      # Mass flow rate in the chanels, kg s^-1
    
    " Mass velocity for chanel CHEK!!!!"
    w_ch = m_dot_ch/(amplitude*width*rho)                                   # Mass velocity in chanels, kg m^2 s^-1
    Re = rho*w_ch*D_h/mu                                          # Reynolds Number, -

    " Factor for correlations: provided by Focke et al."
    if Re >= 2000:                                                  # Regimen: Turbulent
        xhi_0   = (1.8*np.log(Re)-1.5)**-2
        xhi_1_0 = 39/Re**0.289
    elif Re <2000:                                                  # Regime: Laminar
        xhi_0   = 64/Re
        xhi_1_0 = 597/Re +3.85
    
    " Constant given by Martin"
    a = 3.8
    b = 0.18 
    c = 0.36    
    
    " Factor xhi" 
    xhi_1 = a*xhi_1_0

    " Beta angle from degree to radians"
    chevron_angle = 90*np.pi/180 - chevron_angle
    beta_r = chevron_angle     #(90*np.pi/180)-chevron_angle                                   # Chevron angle, in Radians   
     
    " Friction factor equation 18"    
    f = (np.cos(beta_r)/np.sqrt(b*np.tan(beta_r) + c*np.sin(beta_r) + xhi_0/np.cos(beta_r)) +(1 - np.cos(beta_r))/np.sqrt(xhi_1))**(-2)
 
    " Hagen number "
    Hg = f*Re**2/2

    " Pressure Drop:"    
    DeltaP = Hg * (mu**2*length)/(rho*D_h**3)
    
    " Extracted from the comparison with Heavear et al. [10]"
    c_q = 0.122
    q = 0.374
    
    # " Wall temperature "
    # T_wall = T_mean-5                                                # Wall Temperature, K
    # if T_wall<273.15:
    #     T_wall = 274.15
        
    # mu_w =  PropsSI('V','P', P_mean,'T',T_wall, refrigerant)          # Viscosity at wall T, Pa s
    
    " Nusslet number: "
    # Nu = c_q*Pr**(1/3)*(mu/mu_w)**(1/6)*(2*Hg*np.sin(2*beta_r))**q
    Nu = c_q*Pr**(1/3)*(2*Hg*np.sin(2*beta_r))**q
    h_conv = Nu*k/D_h
    # print(f' Nu: {Nu}')
    # print(f' mu_w: {mu_w}')
    # print(f' T_wall : {T_wall}')
    # print(f' Result: {(mu/mu_w)**(1/6)}')
    
    return h_conv

# For the R1233zd(E) two phase (evaporation):
def amalfi_plate_HTC(D_h, length, width, amplitude, chevron_angle, nb_channels, A_tot, m_dot, P_mean, fluid):

    def PHX_EV_Amalfi(x, D_h, length, width, amplitude, chevron_angle, A_tot, m_dot, P_mean, fluid):
        " Gas properties "
        rho_g = PropsSI('D','P', P_mean, 'Q', 1, fluid)
        mu_g = PropsSI('V','P', P_mean, 'Q', 1, fluid)
        h_g = PropsSI('H','P',P_mean,'Q',1, fluid)

        " Liquid properties "    
        rho_l = PropsSI('D','P', P_mean, 'Q', 0, fluid)
        mu_l = PropsSI('V','P', P_mean, 'Q', 0, fluid)
        DeltaH_lg =  PropsSI('H','P',P_mean,'Q',1, fluid)-PropsSI('H','P',P_mean,'Q',0, fluid)

        " Instant quality-properties "
        rho_x = PropsSI('D','P', P_mean, 'Q', x, fluid)  
        if fluid == 'R1233zd(E)':
            T_mean = PropsSI('T','P',P_mean,'Q',1,fluid)
            k_x = (0.09513*T_mean - 17.963)/1000
        else:
            k_x  = PropsSI('L','P', P_mean, 'Q', x, fluid)       
        sigma_x = PropsSI('surface_tension', 'P', P_mean, 'Q', x, fluid)
        h_r_x = PropsSI('H','P',P_mean,'Q', x, fluid)   

        " Average conditions "
        x_m = 0.5*(0 + x)
        rho_m = (x_m/rho_g + (1-x_m)/rho_l)**-1
    
        " Mass velocity"
        G = m_dot/(2*amplitude*width)

        " Heat Flux:"
        Q_dot = m_dot*(h_g - h_r_x)
        q_dot = Q_dot/A_tot 

        " Boiling Number"
        Bo = q_dot/(G*DeltaH_lg)
      
        " Weber number"
        We = G**2*D_h/(rho_m*sigma_x)

        " Reynolds as fully Liquid"
        Re_LT = G*D_h/mu_l
    
        " Reynolds with only the gas fraction "                    
        Re_GS = G*x*D_h/mu_g             
    
        " Bond Number "
        g = 9.81
        Bd = (rho_l-rho_x)*g*D_h**2/sigma_x   
         
        " Normalization of the density"
        ratio_rho  = rho_x/rho_l 
        
        " Normalization of the Beta angles, Maximum 70"
        chevron_angle = chevron_angle
        ratio_beta = chevron_angle/70  

        " Boiling condition for compute the Nusslet Number:"
        if Bd<4:
            Nu_TP = 982*ratio_beta**(1.101)*We**(0.315)*Bo**(0.320)*ratio_rho**(-0.224)
            
        elif Bd>=4:
            Nu_TP = 18.495**ratio_beta**0.248 * Re_GS**0.135 * Re_LT**0.351 * Bd**0.235 * Bo**0.198 * ratio_rho**-0.223      
        hcv_r = Nu_TP*k_x/D_h    

        " Pressure"
        v_g = 1/rho_g                                                 # Liquid Volume Q = 1, kg m^-3
        v_l = 1/rho_l                                                 # Liquid Volume Q = 0, kg m^-3

        " Mass velocity per total"  
        G = m_dot/nb_channels/(width*amplitude)

        " Pressure Drop: Extracted from Longo et al, and based in his ref 16.(Collier et Thome)"
        " Decelaration pressure:"
        DeltaP_a = G**2*(v_g-v_l)*(1-0)
        
        " Gravity pressure rise (elevation):"
        g = 9.81
        DeltaP_g = g*rho_m*length
        
        " Manifold and ports pressure drops (Shah and Focke correlation 1988):"
        DeltaP_c = 1.5*G**2/(2*rho_m)
        
        " Friction pressure drops Longo et al."
        KEV = G**2/(2*rho_m)
        DeltaP_f = 1.95*KEV*1000   # kPa, This values depends on the refrigerant in Longo et al. he uses 2.0
        
        " Total pressure Drop:"
        DeltaP_t = DeltaP_a + DeltaP_g + DeltaP_c + DeltaP_f
        
        " Total Pressure drop:"
        DeltaP_t = DeltaP_a + DeltaP_g + DeltaP_c + DeltaP_f    
        # print(f' hcv_r: {hcv_r}')
        return hcv_r, DeltaP_t
    
    " Quality subdivision "
    N_subdiv = 100
    x_vec = np.linspace(0.1, 1, N_subdiv)

    " h_cv vector "
    hcv = np.zeros(N_subdiv)
    
    " Heat transffer calculation:"
    for i in range (N_subdiv):
        hcv[i], DeltaP_t = PHX_EV_Amalfi(x_vec[i], D_h, length, width, amplitude, chevron_angle, A_tot, m_dot, P_mean, fluid)
    " Average calculation:"
    hcv_mean = np.mean(hcv)

    return hcv_mean

# For the R1233zd(E) two phase (condensation):
def shah_condensation_plate_HTC(D_h, L_p, B_p, b, phi, M_dot_r, P_mean, N_ch, refrigerant):
    """ Shah correlation for condensation in Plate heat Exchanger 2021, for Beta valid 35 and 70°
    
    This is correlation is based in the Longo et al. 2015 correlation. This works presents an 
    imporvement in the prediction.
 
        Shah (2021) Corr. for Condensation in corrugated PHX
        Range: Water, ammmonia, R1234ZE ... (18 fluids)
        Corrugation height (b): [1.2 - 5.0] mm
        Beta : [30-75]°
        Corr. Pitch (lambda): [4.9-12.7]mm
        G : [2.3 - 165] kg/(m2 s)
        x : [0-1] 
        µ

    Inputs:
    -------
    D_h : Hydraulic Diameter, m
    L_p : Port-port centerline distance (Length), m
    B_p : Port-port centerline distance (Width), m
    b : Corrugation height or amplitude, m
    phi : Enlargement factor, -
    M_dot_r : Refrigerant mass flow rate, kg s^-1
    P_mean : Condensation pressure, Pa,
    N_ch : Number of plates of the fluid, -
    refrigerant : Fluid in condensation
"""
    def PHX_CD_Shah_x(x = 0.5, D_hyd = D_h, L_p = L_p, B_p = B_p, b = b, phi = phi, M_dot_r = M_dot_r, P_mean = P_mean, N_ch = N_ch, refrigerant = refrigerant):
        
        " Saturation conditions "
        rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)           # Vapor density Q =1, kg m^-3
        rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3
        v_g = 1/rho_g                                                   # Liquid Volume Q = 1, kg m^-3
        v_l = 1/rho_l                                                   # Liquid Volume Q = 0, kg m^-3
        mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
        cp_l = PropsSI('C','P', P_mean, 'Q', 0, refrigerant)            # Liquid Specific Heat Q = 0, J/kg-K
        T_sat = PropsSI('T','P', P_mean, 'Q', 0, refrigerant)           # Liquid Temperature Q = 0, K
        
        if refrigerant == 'R1233ZD(E)' or refrigerant == 'R1233zd(E)':
            k_l = (-0.2614*T_sat + 159.19)/1000
            Pr_l = cp_l*mu_l/k_l

        else:
            k_l  = PropsSI('L','P', P_mean, 'Q', 0, refrigerant)        # Liquid Conductivity Q = 0, W m^-2 K^-1
            Pr_l = PropsSI('Prandtl', 'P', P_mean,'Q', 0,refrigerant)   # Liquid Prandt number Q =0, -

        " Refrigerant mass flow rate per chanell"
        M_dot_r_ch = M_dot_r/N_ch                                      # Mass flow rate, kg s^-1
        
        " Mass velocity per total"  
        G = M_dot_r_ch/(B_p*b)
        
        " Reynolds number  in liquid phase flow:"
        Re_LS = G*(1-x)*D_hyd/mu_l

        " Equivalent mass velocity: "
        G_eq = G*((1-x)+x*(rho_l/rho_g)**0.5)
        
        " Equivalent Reynolds number:"
        Re_eq = G_eq*D_hyd/mu_l

        " Heat transffer coefficient in forced convection regime (Longo et al. 2015):"
        hcv_fc = 1.875*(k_l/D_hyd)*phi*Re_eq**0.445*Pr_l**(1/3)
        
        " Heat transffer coefficient in gravity controlled regime (Shah et al. 2021):"
        g = 9.81
        hcv_grav = 1.32*phi*Re_LS**(-1/3) *((rho_l*(rho_l-rho_g)*g*k_l**3)/mu_l**2)**(1/3)

        " Shah did not find a consistent trend for the following control:"
        if Re_eq >= 1600:
            hcv_r = hcv_fc
        else:
            hcv_r = max(hcv_grav, hcv_fc)
            
        " Pressure Drop: Extracted from Longo et al, and based in his ref 16.(Collier et Thome)"
        " Decelaration pressure:"
        DeltaP_a = G**2*(v_g-v_l)*(1-0)
        
        " Gravity pressure rise (elevation):"
        x_m = 0.5
        rho_m = (x_m/rho_g + (1-x_m)/rho_l)**(-1)
        DeltaP_g = g*rho_m*L_p
        
        " Manifold and ports pressure drops (Shah and Focke correlation 1988):"
        DeltaP_c = 1.5*G**2/(2*rho_m)
        
        " Friction pressure drops Longo et al."
        KEV = G**2/(2*rho_m)
        DeltaP_f = 1.95*KEV*1000   # kPa, This values depends on the refrigerant in Longo et al. he uses 2.0
        
        " Total pressure Drop:"
        DeltaP_t = DeltaP_a + DeltaP_g + DeltaP_c + DeltaP_f

        return hcv_r, DeltaP_t
    
    " Quality  calculation: "
    
    " Number of subdivision "
    N_subdiv = 100
    
    " Vector of quality, x=1, problem in Re_L"
    x_vec = np.linspace(0, 0.99, N_subdiv) 
    hcv = np.zeros(N_subdiv)
    
    " Calculation along the differents qualities"
    for i in range (N_subdiv):
        " Instant heat transfer calculation "
        hcv[i], DeltaP_t = PHX_CD_Shah_x(x_vec[i],
                                D_h,
                                L_p,
                                B_p,
                                b, 
                                phi, 
                                M_dot_r, 
                                P_mean,
                                N_ch, 
                                refrigerant)
 
    hcv_mean = np.mean(hcv)
    return hcv_mean



# POUR LE 1 PHASE
def thonon_plate_HTC(mu, Pr, k, G, Dh, chevron_angle):
    """
    Reference
    ----------
    Thonon, B., Design Method for Plate Evaporators and Condensers, 1st International Conference on Process Intensification
    for the Chemical Industry, BHR Group Conference Series Publication, no. 18, pp. 37–47, 1995.
    """
    Re = G*Dh/mu
    print('Re_c', Re)
    print("G_c", G)
    print("Dh_c", Dh)
    print("mu_c", mu)

    if Re < 50 or Re > 15000:
        raise ValueError(f"Reynolds number ({Re:.2f}) is out of bounds for the Thonon correlation (50 <= Re <= 15,000).")

    if chevron_angle == 75*np.pi/180:
        C1 = 0.1
        m = 0.687
        if Re <= 1000:
            C2 = 28.21
            p = 0.9
        else:
            C2 = 0.872
            p = 0.392
    if chevron_angle == 60*np.pi/180:
        C1 = 0.2267
        m = 0.631
        if Re <= 550:
            C2 = 26.34
            p = 0.830
        else:
            C2 = 0.572
            p = 0.217
    if chevron_angle == 45*np.pi/180:
        C1 = 0.2998
        m = 0.645
        if Re <= 200:
            C2 = 18.19
            p = 0.682
        else:
            C2 = 0.6857
            p = 0.172
    if chevron_angle == 30*np.pi/180:
        C1 = 0.2946
        m = 0.700
        if Re <= 160:
            C2 = 45.57
            p = 0.670
        else:
            C2 = 0.370
            p = 0.172

    Nu = C1*Re**m*Pr**1/3
    h_conv = (k/Dh)*Nu
    return h_conv

def kumar_plate_HTC(mu, Pr, k, G, Dh, mu_wall, chevron_angle):
    """
    Reference
    ----------
    Kumar, H., The Plate Heat Exchanger: Construction and Design, Institute of Chemical Engineering Symposium Series,
    no. 86, pp. 1275–1288, 1984.
    """

    Re = G*Dh/mu

    if chevron_angle <=30*np.pi/180:
        if Re<=10:
            C1 = 0.718
            m = 0.349
            C2 = 50
            p = 1
        if Re >10:
            C1 = 0.348
            m = 0.663
            if Re < 100:
                C2 = 19.4
                p = 0.589
            if Re >= 100:
                C2 = 2.99
                p = 0.183
    if chevron_angle == 45*np.pi/180:
        if Re<10:
            C1 = 0.718
            m = 0.349
        if Re>10 and Re<100:
            C1 = 0.4
            m = 0.598
        if Re>100:
            C1 = 0.3
            m = 0.663
        if Re<15:
            C2 = 47
            p = 1
        if Re>15 and Re<300:
            C2 = 18.29
            p = 0.652
        if Re > 300:
            C2 = 1.441
            p = 0.206
    if chevron_angle == 50*np.pi/180:
        if Re<20:
            C1 = 0.63
            m = 0.333
            C2 = 34
            p = 1
        if Re>20 and Re<300:
            C1 = 0.291
            m = 0.591
            C2 = 11.25
            p = 0.631
        if Re>300:
            C1 = 0.13
            m = 0.732
            C2 = 0.772
            p = 0.161
    if chevron_angle == 60*np.pi/180:
        if Re<20:
            C1 = 0.562
            m = 0.326
        if Re>20 and Re<400:
            C1 = 0.306
            m = 0.529
        if Re>400:
            C1 = 0.108
            m = 0.703
        if Re<40:
            C2 = 24
            p = 1
        if Re>40 and Re<400:
            C2 = 3.24
            p = 0.457
        if Re>400:
            C2 = 0.76
            p = 0.215
    if chevron_angle>=65*np.pi/180:
        if Re<20:
            C1 = 0.562
            m = 0.326
        if Re>20 and Re<500:
            C1 = 0.331
            m = 0.503
        if Re>500:
            C1 = 0.087
            m = 0.718
        if Re<50:
            C2 = 24
            p = 1
        if Re>50 and Re<500:
            C2 = 2.8
            p = 0.451
        if Re>500:
            C2 = 0.639
            p = 0.213

    Nu = C1*(Re**m)*(Pr**0.33)*((mu/mu_wall)**0.17) # Water, herringbone plates, φ = 1.17
    f = C2/((Re)**p)

    h_conv = (k/Dh)*Nu
    return h_conv
        

def simple_plate_HTC(mu, Pr, k, G, Dh):
    # Reynolds number
    Re = G*Dh/mu
    
    if Re < 5*1e5:
        Nu = 0.3387*Re**(1/2)*Pr**(1/3)/(1+(0.0468/Pr)**(2/3))**(1/4)
    else:
        Nu = 0.0296*Re**(4/5)*Pr**(1/3)
        
    h_conv = (k/Dh)*Nu
    
    return h_conv


def muley_manglik_BPHEX_HTC(mu, mu_w, Pr, k, G, Dh, chevron_angle):
    # Reynolds number
    Re = G*Dh/mu
    
    beta = 180*chevron_angle/np.pi
    
    C = 0.2668 - 0.006967*beta + 7.244*1e-5*beta**2
    C_2 = Re**(0.728 + 0.0543*np.sin((2*np.pi*beta/90) + 3.7))
    
    Nu = C * C_2 * Pr**(1/3) * (mu/mu_w)**(0.14)
        
    h_conv = (k/Dh)*Nu
    
    return h_conv

def martin_BPHEX_HTC(mu, mu_w, Pr, k, G, Dh, chevron_angle):
    "Martin Holger: Correlation from the VDI Atlas p.1515"
 
    beta = chevron_angle
    Re = G*Dh/mu
    
    "Factor for correlations: provided by Focke et al."
    if Re >= 2000: # Regime : Turbulent
        xhi_0   = (1.8*np.log(Re)-1.5)**-2
        xhi_1_0 = 39/Re**0.289
    elif Re < 2000: # Regime: Laminar
        xhi_0   = 64/Re
        xhi_1_0 = 597/Re +3.85
    
    "Constant given by Martin"
    a = 3.8
    b = 0.18 
    c = 0.36    
    
    "Factor xhi"
    xhi_1 = a*xhi_1_0
    
    "Friction factor"    
    f = (np.cos(beta)/np.sqrt(b*np.tan(beta) + c*np.sin(beta) + xhi_0/np.cos(beta)) +(1 - np.cos(beta))/np.sqrt(xhi_1))**(-2)
    
    "Hagen number"
    Hg = f*Re**2/2
    
    "Extracted from the comparison with Heavear et al. [10]"
    c_q = 0.122
    q = 0.374
    
    "Nusslet number:"
    Nu = c_q*Pr**(1/3)*(mu/mu_w)**(1/6)*(2*Hg*np.sin(2*beta))**q
    
    "Heat Transffer Coefficient [W m^-2]:"
    hcv = Nu*k/Dh
    
    return hcv 


#%%
def han_boiling_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, i_fg, G, DT_log, Qdot, hconv_h, Dh, theta, pitch_co):
    """
    Inputs
    ------
    x        : Vapor quality [-]
    mu_l     : Liquid viscosity [Pa*s]
    k_l      : Liquid thermal conductivity [W/(m*K)]
    Pr_l     : Liquid Prandtl Number [-]
    rho_l    : Liquid density [kg/m^3]
    rho_v    : Vapor density [kg/m^3]
    i_fg     : Vaporisation heat [J/(kg)]
    G        : Mass flux [kg/(m^2 * s)]
    DT_log   : Logarithmic temperature difference [K]    
    Qdot     : Heat Transfer rate [W]
    hconv_h  : Convection heat transfer coefficient [W/(m^2*K)]
    Dh       : Plate spacing [m]
    theta    : Chevron angle [°]
    pitch_co : Corrugated pitch [m]
    
    Outputs
    -------
    h_cond : Condensation HTC [W/(m^2*K)]
    Nu     : Nusselt Number [-]
    DP_tot : Total Pressure Drop [Pa]
    
    Reference
    ---------
    HFC-410A vaporisation inside a commercial brazed plate heat exchanger
    Giovanni A. Longo, Andrea Gasparella
    """
    
    def iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, hconv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        Bo_g = max(Bo_g, 1e-8)
        Nu = Ge1*Re_eq**Ge2*Bo_g**0.3*Pr_l**0.4
        h = Nu*k_l/Dh
        U = (1/h +  1/hconv_h)**-1
        A_tp = AU_tp/U
        q = Qdot/A_tp
        Bo = q/(G_eq*i_fg)
        res_Bo = (Bo-Bo_g)/Bo_g
        return res_Bo, Nu, h, U, A_tp, q, Bo
    
    def res_iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        res_Bo,_,_,_,_,_,_ = iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
        
        return res_Bo
    
    G_eq = G * ( (1 - x) + x * (rho_l/rho_v)**0.5)
    Re_eq = G_eq*Dh/mu_l
    AU_tp = Qdot/DT_log
    Ge1 = 2.81*(pitch_co/Dh)**(-0.041)*(theta)**(-2.83)
    Ge2 = 0.746*(pitch_co/Dh)**(-0.082)*(theta)**(0.61)
    Bo_0 = 0.5
    f_Bo = lambda xx: res_iter_Han_boiling(xx, Ge1, Ge2,Re_eq, Pr_l, k_l, hconv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    sol = fsolve(f_Bo, Bo_0,)
    Bo_sol = sol[0]
    _, Nu, h, _, _, _, _ = iter_Han_boiling(Bo_sol, Ge1, Ge2,Re_eq, Pr_l, k_l, hconv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    
    h_boiling = h
    
    return h_boiling, Nu

def han_cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, G, Dh, pitch_co, beta, L_v, N_cp, m_dot, D_p):
    """
    Inputs
    ------
    x        : Vapor quality [-]
    mu_l     : Liquid viscosity [Pa*s]
    k_l      : Liquid thermal conductivity [W/(m*K)]
    Pr_l     : Liquid Prandtl Number [-]
    rho_l    : Liquid density [kg/m^3]
    rho_v    : Vapor density [kg/m^3]
    G        : Mass flux [kg/(m^2 * s)]
    Dh       : Plate spacing [m]
    pitch_co : Corrugated pitch [m]
    beta     : Chevron angle [°]
    L_v      : Vertical length between fluid ports [m] 
    N_cp     : Number of canals [-]
    m_dot    : Flowrate [kg/s]
    
    Outputs
    -------
    h_cond : Condensation HTC [W/(m^2*K)]
    Nu     : Nusselt Number [-]
    DP_tot : Total Pressure Drop [Pa]
    
    Reference
    ---------
    "The caracteristics of condensation in brazed plate heat exchangers with different chevron angles", Journal of the Korean Physical Society, 2003
    Han et. al
    
    """
    
    # Preliminary calculations
    theta = np.pi/2 - beta
    g = 9.81 # gravity acceleration
    
    G_eq = G * ( (1-x) + x * (rho_l/rho_v)**0.5)
    Re_eq = G_eq*Dh/mu_l
    
    # Heat Transfer
    Ge1 = 11.22*(pitch_co/Dh)**(-2.83)*(theta)**(-4.5)
    Ge2 = 0.35*(pitch_co/Dh)**(0.23)*(theta)**(1.48)
    
    Nu = Ge1*Re_eq**Ge2*Pr_l**(1/3)
    h_cond = Nu*k_l/Dh
    
    # Pressure drop
    Ge3 = 3521.1*(pitch_co/Dh)**(4.17)*(theta)**(-7.75)
    Ge4 = -1.024*(pitch_co/Dh)**(0.0925)*(theta)**(-1.3)
    
    f = Ge3*Re_eq**Ge4
    
    # Two phase related pressure drop
    DP_tp = f*(L_v*N_cp/Dh)*G_eq**2*rho_l
    
    # Port pressure drop
    m_dot_eq = m_dot*(1 - x + x*(rho_l/rho_v)**0.5)
    G_p = 4*(m_dot_eq/(np.pi*D_p**2))
    rho_m = 1/( (x/rho_v) + (1 - x)/rho_l )

    DP_port = 1.4*G_p**2/(2*rho_m)

    # Static head loss
    DP_stat = -rho_m*g*L_v # negative because downward flow <-> Condenser

    # The acceleration pressure drop for condensation is expressed as : ??? 
    
    DP_tot = DP_tp + DP_port + DP_stat
    
    return h_cond, Nu, DP_tot

def han_BPHEX_DP(mu, G, Dh, beta, pitch_co, rho_v, rho_l, L_v, N_cp, m_dot, D_p): 
    """
    Inputs
    ------
    mu       : viscosity [Pa*s]
    rho_l    : Liquid density [kg/m^3]
    rho_v    : Vapor density [kg/m^3]
    G        : Mass flux [kg/(m^2 * s)]
    Dh       : Plate spacing [m]
    pitch_co : Corrugated pitch [m]
    beta     : Chevron angle [°]
    L_v      : Vertical length between fluid ports [m]
    N_cp     : Number of canals [-]
    m_dot    : Flowrate [kg/s]
    D_p      : Port diameter [m]
    
    Outputs
    -------
    h_cond : Condensation HTC [W/(m^2*K)]
    Nu     : Nusselt Number [-]
    DP_tot : Total Pressure Drop [Pa]
    
    Reference
    ---------
    "The caracteristics of condensation in brazed plate heat exchangers with different chevron angles", Journal of the Korean Physical Society, 2003
    Han et. al
    
    """
    
    # Preliminary calculations
    theta = np.pi/2 - beta
    g = 9.81 # gravity acceleration
    
    Re_eq = G*Dh/mu
    
    # Pressure drop
    Ge3 = 3521.1*(pitch_co/Dh)**(4.17)*(theta)**(-7.75)
    Ge4 = -1.024*(pitch_co/Dh)**(0.0925)*(theta)**(-1.3)
    
    f = Ge3*Re_eq**Ge4
    
    rho = rho_v + rho_l
    
    # Two phase related pressure drop
    DP_tp = f*(L_v*N_cp/Dh)*G**2*rho
    
    # Port pressure drop
    G_p = 4*(m_dot/(np.pi*D_p**2))

    DP_port = 1.4*G_p**2/(2*rho)

    # Static head loss
    DP_stat = -rho*g*L_v # negative because downward flow <-> Condenser

    # The acceleration pressure drop for condensation is expressed as : ??? 
    
    DP_tot = DP_tp + DP_port + DP_stat
    
    return DP_tot