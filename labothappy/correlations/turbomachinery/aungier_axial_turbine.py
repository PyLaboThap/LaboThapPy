# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 21:00:27 2025

@author: basil
"""

import numpy as np

def profile_losses(alpha1_pr, alpha2_pr, beta1, c, e, mu2, M1_pr, M2_pr, R_c, rho2, s, t_max, w2):#
    """
    alpha_1_pr : inlet flow angle with respect to tangent in the rotating coordinate system [°] 
    alpha_2_pr : outlet flow angle with respect to tangent in the rotating coordinate system [°] 
    beta_1     : inlet blade angle with respect to tangent [°] 
    c          : blade cord [m]
    e          : blade roughness [m]
    mu2        : outlet dynamic viscosity evaluated with static conditions [Pa*s]
    M1_pr      : inlet Mach number in the rotating coordinate system (using w2) [-]
    M2_pr      : outlet Mach number in the rotating coordinate system (using w2) [-]
    R_c        : suction surface radius of curvature [m]
    rho2       : outlet mass density evaluated with static conditions [kg/m^3]
    s          : blade pitch [m]
    t_max      : blade maximum thickness [m]
    w2         : outlet flow speed in the rotating coordinate system [m/s]
    """

    xi = (90 - beta1)/(90 - alpha2_pr)

    "1 Experience factor (Kacker-Okapuu)"
    Kmod = 0.67 # For modern optimized design (2006) - =1 for old designs
    
    "2) Off-design incidence correction"
    i = beta1 - alpha1_pr
    
    "2.1) Compute i_s"
    
    # 2.1.1) i_sr : reference is
    i_s0 = 20 - (xi + 1)/0.11
    
    if alpha2_pr <= 40:
        A = 61.8 - (1.6-alpha2_pr/165)*alpha2_pr
        B = 71.9 - 1.69*alpha2_pr
        C = 7.8 - (0.28-alpha2_pr/320)*alpha2_pr
        D = 14.2 - (0.16+alpha2_pr/160)*alpha2_pr # written alpha_2 in the book and not alpha_2_pr but it seems like the total expression depends only on alpha_2_pr
        
        i_sr = i_s0 + A - B*xi**2 + C*xi**3 + D*xi**4
    else:
        A_40 = 61.8 - (1.6-40/165)*40
        B_40 = 71.9 - 1.69*40
        C_40 = 7.8 - (0.28-40/320)*40
        D_40 = 14.2 - (0.16+40/160)*40
        
        i_sr_40 = i_s0 + A_40 - B_40*xi**2 + C_40*xi**3 + D_40*xi**4
        
        i_sr = i_s0 + abs(i_sr_40-i_s0) * abs(55-alpha2_pr) / 15
        
    # 2.1.2) Dis : Correction factor from reference is
    c = 1
    X = s/c - 0.75
    
    if s/c <= 0.8: 
        Dis = -38*X - 53.5*X**2 - 29 * X**3
    else:
        Dis = -2.0374 - (s/c - 0.8)*(69.58 - (alpha2_pr/14.48)**3.1)
    
    i_s = i_sr + Dis
    
    if i/i_s < -3:
        Kinc = -1.39214 - 1.90738*(i/i_s)
    elif i/i_s < 0:
        Kinc = 1 + 0.52*abs(i/i_s)**1.7
    elif i/i_s < 1.7:
        Kinc = 1 + (i/i_s)**(2.3+0.5*(i/i_s))
    else:
        Kinc = 6.23 + 9.8577*((i/i_s)-1.7)
        
    Kinc = min(Kinc,20)
    
    "3) Mach number correction"
    if M2_pr <= 0.6:
        Km = 1
    else:
        Km = 1 + (1.65*(M2_pr-0.6) + 240*(M2_pr-0.06)**2)*(s/R_c)**(3*M2_pr - 0.6)
    
    "4) Compressibility correction"
    M1_corr = (M1_pr + 0.566 - abs(0.566-M1_pr))/2
    M2_corr = (M2_pr + 1 - abs(M2_pr - 1))/2
    
    X = 2*M1_corr/(M1_corr+M2_corr+abs(M2_corr - M1_corr))
    
    K1 = 1-0.625*(M2_corr-0.2+abs(M2_corr-0.2))
    
    Kp = 1 - (1-K1)*X**2
    
    "5) Reynolds number correction"
    Re_c = rho2*w2*c/mu2
    
    if e == 0:
        Re_r = -1
    else:
        Re_r = 100*c/e

    coef_r = (np.log10(5*1e5)/np.log10(Re_r))**2.58

    if Re_c < 1e5: # laminar regime
        Kre = np.sqrt(1e5/Re_c)
    else:
        if Re_c < 5*1e5: # Transition regime    
            Kre = 1
        else: # Turbulent regime
            if Re_c > Re_r > 0:
                Kre = coef_r
            else:
                Kre = (np.log10(5*1e5)/np.log10(Re_c))**2.58
                
        
    "6) Profile loss coefficient for nozzle blades (beta_1 = 90°)"

    if alpha2_pr <= 27:
        A = 0.025 + (27-alpha2_pr)/530
    else:
        A = 0.025 + (27-alpha2_pr)/3085
    
    B = 0.1583 - alpha2_pr/1640
    
    s_c = s/c
    
    if alpha2_pr <= 30:
        s_c_min = 0.46 + alpha2_pr/77
        X = s_c - s_c_min
        
        C = 0.08*((alpha2_pr/30)**2 - 1)

        Yp1 = A + B*X**2 + C*X**3
    else:
        s_c_min = 0.614 + alpha2_pr/130
        X = s_c - s_c_min
        
        n = 1 + alpha2_pr/30
        
        Yp1 = A + B*abs(X)**n
    
    "7 Profile loss coefficient for impulse blades (alpha_2' = 90°)"
    s_c_min = 0.224 + 1.575*(alpha2_pr/90) - (alpha2_pr/90)**2
    X = s_c - s_c_min
    
    A = 0.242 - alpha2_pr/151 + (alpha2_pr/127)**2
    
    if alpha2_pr <= 30:
        B = 0.3 + (30-alpha2_pr)/50
    else:
        B = 0.3 + (30-alpha2_pr)/275
        
    C = 0.88 - alpha2_pr/42.4 + (alpha2_pr/72.8)**2

    Yp2 = A + B*X**2 - C*X**3
    
    "8) DY_TE Substract trailing edge loss for t = 0.02*s"
    o = s * np.sin(alpha2_pr) # throat opening
    beta_g = np.arcsin(o/s) # gauging angle
    
    t2 = 0.02*s
    DPt = 1/2 * rho2 * w2**2 * (s*np.sin(beta_g)/(s*np.sin(beta_g)-t2)-1)**2
    
    DY_TE = 2*DPt/(rho2*w2**2)
    
    "9) Total"
    fact_1 = Yp1 + xi**2*(Yp2 - Yp1)
    fact_2 = (5*t_max/c)**xi 
    
    Yp = Kmod * Kinc * Km * Kp * Kre * (fact_1 * fact_2 - DY_TE)    

    return Yp

def secondary_losses(alpha1_pr, alpha2_pr, beta1, b_z, c, e, h, mu2, M1_pr, M2_pr, rho2, s, w2):#
    """
    alpha_1_pr : inlet flow angle with respect to tangent in the rotating coordinate system [°]
    alpha_2_pr : outlet flow angle with respect to tangent in the rotating coordinate system [°]
    beta_1     : inlet blade angle with respect to tangent [°] 
    b_z        : axial projection of cord length [m]
    c          : blade cord [m]
    e          : blade roughness [m]
    h          : blade height [m]
    mu2        : outlet dynamic viscosity evaluated with static conditions [Pa*s]
    M1_pr      : inlet Mach number in the rotating coordinate system (using w2) [-]
    M2_pr      : outlet Mach number in the rotating coordinate system (using w2) [-]
    rho2       : outlet mass density evaluated with static conditions [kg/m^3]
    s          : blade pitch [m]
    w2         : outlet flow speed in the rotating coordinate system [m/s]
    """
    
    "1) Lift Coefficient"
    C_L = 2*(np.cot(alpha1_pr)+np.cot(alpha2_pr))*s/c

    "2) Ainley Loading Parameter"
    alpha_m_pr = 90 - np.arctan((np.cot(alpha1_pr)-np.cot(alpha2_pr))/2)

    Z = (C_L*c/s)**2 * np.sin(alpha2_pr)**2 / np.sin(alpha_m_pr)**3

    "3) Aspect ratio correction"
    if h/c >= 2:
        F_AR = c/h
    else:
        F_AR = 0.5*(2*c/h)**0.7

    "4) Preliminary Estimate of Ys"
    Ys_est = 0.0334*F_AR*Z*np.sin(alpha2_pr)/np.sin(beta1)
    
    "5) Reynolds number correction"

    Re_c = rho2*w2*c/mu2
        
    if e == 0:
        Re_r = -1
    else:
        Re_r = 100*c/e

    coef_r = (np.log10(5*1e5)/np.log10(Re_r))**2.58

    if Re_c < 1e5: # laminar regime
        Kre = np.sqrt(1e5/Re_c)
    else:
        if Re_c < 5*1e5: # Transition regime    
            Kre = 1
        else: # Turbulent regime
            if Re_c > Re_r > 0:
                Kre = coef_r
            else:
                Kre = (np.log10(5*1e5)/np.log10(Re_c))**2.58
    
    "6) Compressibility correction"
    M1_corr = (M1_pr + 0.566 - abs(0.566-M1_pr))/2
    M2_corr = (M2_pr + 1 - abs(M2_pr - 1))/2
    
    X = 2*M1_corr/(M1_corr+M2_corr+abs(M2_corr - M1_corr))
    K1 = 1-0.625*(M2_corr-0.2+abs(M2_corr-0.2))
    
    Kp = 1 - (1-K1)*X**2
    
    "7) Modified Compressibility Correction"
    Ks = 1 - (1-Kp) * (b_z/h)**2 / (1+(b_z/h)**2)
    
    Ys = Ks*Kre*np.sqrt(Ys_est**2 / (1+7.5*Ys_est**2))

    return Ys

def aungier_loss_model(alpha1_pr, alpha2_pr, beta1, c, delta, D_lw, e, h, i, mu2, M1_pr, M2_pr, N_lw, R_c, rho2, s, t_TE, t_max, vm2, w2):
    """
    alpha_1_pr : inlet flow angle with respect to tangent in the rotating coordinate system [rad] = beta_1 ? 
    alpha_2_pr : outlet flow angle with respect to tangent in the rotating coordinate system [rad] = beta_2 ? 
    beta_1     : inlet blade angle with respect to tangent [rad] = khi 2 ? 
    c          : blade cord [m]
    delta      : blade tip clearance [m]
    D_lw       : Lashing wire diameter [m]
    e          : blade roughness [m]
    h          : blade height [m]
    i          : incidence angle [rad]
    mu2        : outlet dynamic viscosity evaluated with static conditions [Pa*s]
    M1_pr      : inlet Mach number in the rotating coordinate system (using w2) [-]
    M2_pr      : outlet Mach number in the rotating coordinate system (using w2) [-]
    N_lw       : number of lashing wires [-]
    R_c        : suction surface radius of curvature [m]
    rho2       : outlet mass density evaluated with static conditions [kg/m^3]
    s          : blade pitch [m]
    t_TE       : trailing edge thickness [m]
    t_max      : blade maximum thickness [m]
    vm2        : outlet axial velocity [m/s]
    w2         : outlet flow speed in the rotating coordinate system [m/s]
    """
    
    alpha1_pr = (180/np.pi)*alpha1_pr # [rad] to [°]    
    alpha2_pr = (180/np.pi)*alpha2_pr # [rad] to [°]
    beta1 = (180/np.pi)*beta1 # [rad] to [°]
    
    #%% "1) Total Pressure loss"
    
    "1.1 Profile Losses"
    
    Yp = profile_losses(alpha1_pr, alpha2_pr, beta1, c, e, i, mu2, M1_pr, M2_pr, R_c, rho2, s, t_max, w2)
    
    "1.2 Secondary Losses"
    Ys = secondary_losses(alpha1_pr, alpha2_pr, beta1, c, e, h, i, mu2, M1_pr, M2_pr, R_c, rho2, s, t_max, w2)
    
    "1.3 Blade Clearance Losses"
    # 1) Lift Coefficient
    C_L = 2*(np.cot(alpha1_pr)+np.cot(alpha2_pr))*s/c

    # 2) Ainley Loading Parameter
    alpha_m_pr = 90 - np.arctan((np.cot(alpha1_pr)-np.cot(alpha2_pr))/2)
    Z = (C_L*c/s)**2 * np.sin(alpha2_pr)**2 / np.sin(alpha_m_pr)**3
    
    Ycl = 0.47*Z*(c/h)*(delta/c)**0.78
    
    "1.4 Trailing Edge Losses"
    o = s * np.sin(alpha2_pr) # throat opening
    beta_g = np.arcsin(o/s) # gauging angle
    
    DPt = 1/2 * rho2 * w2**2 * (s*np.sin(beta_g)/(s*np.sin(beta_g)-t_TE)-1)**2
    
    Yte = 2*DPt/(rho2*w2**2)    
    
    "1.5 Supersonic Expansion Losses"
    if M1_pr <= 1:
        Yex = 0
    else:
        Yex = ((M2_pr-1)/M2_pr)**2
        
    "1.6 Shock Loss"
    
    if M1_pr <= 0.4:
        X1 = 0 
    else:
        X1 = M1_pr - 0.4
    
    if M1_pr <= M2_pr:
        X2 = 0
    else:
        X2 = M1_pr/M2_pr - 1
    
    Ysh_est = 0.8*X1**2 + X2**2
    Ysh = np.sqrt(Ysh_est**2 / (1+Ysh_est**2))
    
    "1.7 Lashing Wire Losses"
    Re = rho2*vm2*D_lw/mu2
    
    if Re <= 5*1e5:
        C_D = 1
    else:
        C_D = 0.35
    
    Ylw = N_lw*C_D*D_lw*vm2**2 / (h*w2**2)
    
    "1.8 Total"
    # Y = (Pt1' - Pt2')/(Pt2' - P2)
    Y = Yp + Ys + Ycl + Yte + Yex + Ysh + Ylw
    
    #%% "2) Parasitic loss"
    
    "2.1 Partial admission work (rotor)"
    DHadm = 0

    "2.2 Disk friction work (diaphragm-disk rotor)"
    DHDF = 0
    
    "2.3 Leakage bypass loss (shrouded rotors and rotor balance holes)"
    DHlea = 0
    
    "2.4 Clearance gap windage loss (shrouded blades and nozzle of diaphragm-disk rotors)"
    DHgap = 0
    
    "2.5 Moisture work loss (rotors)"
    DHQ = 0
    
    "2.6 Total"
    DHpar = DHadm + DHDF + DHlea + DHgap + DHQ
    
    return




