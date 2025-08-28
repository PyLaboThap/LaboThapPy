# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 21:00:27 2025

@author: basil
"""

import numpy as np

def sind(x):
    return np.sin(np.deg2rad(x))

def cosd(x):
    return np.cos(np.deg2rad(x))

def tand(x):
    return np.tan(np.deg2rad(x))

def cotd(x):
    return 1/np.tan(np.deg2rad(x))

def profile_losses(alpha1_pr, alpha2_pr, beta_g, beta1, c, e, mu2, M1_pr, M2_pr, R_c, rho2, s, t_max, w2):#
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

    xi = (90-beta1)/(90-alpha2_pr) # (beta1)/(alpha2_pr)
        
    if abs(xi) > 1.3:
        raise ValueError("abs(xhi) value above 1.2")

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
    elif M2_pr <= 1:
        Km = 1 + (1.65*(M2_pr-0.6) + 240*(M2_pr-0.06)**2)*(s/R_c)**(3*M2_pr - 0.6)
    else:
        Km = 1
    
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
    t2 = 0.02*s
    denom = s*sind(beta_g) - t2
    denom = np.sign(denom) * max(abs(denom), 1e-12)
    DPt = 0.5*rho2*w2**2 * ( (s*sind(beta_g)/denom) - 1 )**2    
    
    DY_TE = 2*DPt/(rho2*w2**2)
    
    "9) Total"
    fact_1 = Yp1 + xi**2*(Yp2 - Yp1)
    tau = np.clip(5*t_max/c, 1e-3, 1e3)
    fact_2 = tau**xi 
    
    Yp = Kmod * Kinc * Km * Kp * Kre * (fact_1 * fact_2 - DY_TE)    

    if __name__ == "__main__":

        print(f"alpha1_pr: {alpha1_pr}")
        print(f"alpha2_pr: {alpha2_pr}")
        print(f"beta1: {beta1}")        
    
        print(f"Yp: {Yp}")
        print(f"Kmod: {Kmod}")
        print(f"Yp2: {Yp2}")
        print(f"Yp1: {Yp1}")
        print(f"xi: {xi}")
        print(f"isr: {i_sr}")
    
        print(f"Kinc: {Kinc}")
        print(f"Km: {Km}")
        print(f"Kp: {Kp}")
        print(f"Kre: {Kre}")
        print(f"fact_1*fact_2: {fact_1*fact_2}")
        print(f"fact_1: {fact_1}")
        print(f"fact_2: {fact_2}")
        print(f"DY_TE: {DY_TE}")
        print(f"---------------------")

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
    C_L = 2*(1/tand(alpha1_pr)+1/tand(alpha2_pr))*s/c

    "2) Ainley Loading Parameter"
    alpha_m_pr = 180*np.arctan((1/tand(alpha1_pr)+1/tand(alpha2_pr))/2)/np.pi # 90 -

    Z = (C_L*c/s)**2 * sind(alpha2_pr)**2 / sind(alpha_m_pr)**3

    # if __name__ == "__main__":
    #     print(f"sin_ratio : {sind(alpha2_pr)**2 / sind(alpha_m_pr)**3}")
    #     print(f"cos_ratio : {cosd(alpha2_pr)**2 / cosd(alpha_m_pr)**3}")
    #     print(f"Z : {Z}")
    #     print(f"C_L : {C_L}")
    #     print(f"c/s : {c/s}")
    #     print(f"---------------------")

    "3) Aspect ratio correction"
    if h/c >= 2:
        F_AR = c/h
    else:
        F_AR = 0.5*(2*c/h)**0.7

    "4) Preliminary Estimate of Ys"
    Ys_est = 0.0334*F_AR*Z*sind(alpha2_pr)/sind(beta1)
    
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

    # if __name__ == "__main__":
    #     print(f"Ks : {Ks}")
    #     print(f"F_AR : {F_AR}")
    #     print(f"Kre : {Kre}")
    #     print(f"fact : {np.sqrt(Ys_est**2 / (1+7.5*Ys_est**2))}")
    #     print(f"Ys_est : {Ys_est}")
    
    #     print(f"Ys : {Ys}")

    #     print(f"---------------------")


    return Ys

# a2, delta_seal, D_bh, D_lw, epsilon, h1, h2, i, mdot, N_bh, N_lw, N_seal, P1, P2, R, r_m, r_seal, s_seal, t, T2, u_m, vm2

def aungier_loss_model(alpha1_pr, alpha2_pr, beta_g, beta1, c, delta, D_lw, e, h, mu2, M1_pr, M2_pr, N_lw, R_c, rho2, s, t_max, t_TE, vm2, w2, rad_flag):
    """
    alpha_1_pr : inlet flow angle with respect to tangent in the rotating coordinate system [rad] = beta_1 ? 
    alpha_2_pr : outlet flow angle with respect to tangent in the rotating coordinate system [rad] = beta_2 ? 
    beta_1     : inlet blade angle with respect to tangent [rad] = khi 2 ? 
    c          : blade cord [m]
    delta      : blade tip clearance [m]
    delta_seal : seal clearance [m]
    D_bh       : balance holes diameter [m]
    D_lw       : Lashing wire diameter [m]
    e          : blade roughness [m]
    epsilon    : part of blocked flow through the rotor [-]
    h          : blade height [m]
    h1         : inlet enthalpy [J/kg]
    h2         : outlet enthalpy [J/kg]
    i          : incidence angle [rad]
    mdot       : mass flow rate [kg/s]
    mu2        : outlet dynamic viscosity evaluated with static conditions [Pa*s]
    M1_pr      : inlet Mach number in the rotating coordinate system (using w2) [-]
    M2_pr      : outlet Mach number in the rotating coordinate system (using w2) [-]
    N_bh       : number of balancing holes [-]
    N_lw       : number of lashing wires [-]
    N_seal     : number of sealing strips [-]
    P1         : inlet pressure [Pa]
    P2         : outlet pressure [Pa]
    R          : gas constant [J/(mol*K)]
    R_c        : suction surface radius of curvature [m]
    r_m        : rotor mean radius [m]
    r_seal     : seal radius [m]
    rho2       : outlet mass density evaluated with static conditions [kg/m^3]
    s          : blade pitch [m]
    s_seal     : seal pitch [m]
    t          : blade local thickness [m]
    T2         : outlet temperature [K]
    t_TE       : trailing edge thickness [m]
    t_max      : blade maximum thickness [m]
    u_m        : mean rotor blade velocity [m/s]
    vm2        : outlet axial velocity [m/s]
    w2         : outlet flow speed in the rotating coordinate system [m/s]
    """
    
    if rad_flag:
        alpha1_pr = (180/np.pi)*alpha1_pr # [rad] to [°]    
        alpha2_pr = (180/np.pi)*alpha2_pr # [rad] to [°]
        beta1 = (180/np.pi)*beta1 # [rad] to [°]
        
    alpha1_pr = 90 - alpha1_pr
    alpha2_pr = 90 - alpha2_pr # + 90
    beta1 = 90 - beta1
    
    # if __name__ == "__main__":
    #     print(f"alpha1_pr:{alpha1_pr}")
    #     print(f"alpha2_pr:{alpha2_pr}")
    #     print(f"beta1:{beta1}")
    
    #%% "1) Total Pressure loss"
    
    "1.1 Profile Losses"
    
    Yp = profile_losses(alpha1_pr, alpha2_pr, beta_g, beta1, c, e, mu2, M1_pr, M2_pr, R_c, rho2, s, t_max, w2)
    
    b_z = c*sind(beta1)
    
    "1.2 Secondary Losses"
    Ys = secondary_losses(alpha1_pr, alpha2_pr, beta1, b_z, c, e, h, mu2, M1_pr, M2_pr, rho2, s, w2)
    
    "1.3 Blade Clearance Losses"
    # 1) Lift Coefficient
    C_L = 2*(1/tand(alpha1_pr)+1/tand(alpha2_pr))*s/c

    # 2) Ainley Loading Parameter
    alpha_m_pr = 90 - 180*np.arctan((1/tand(alpha1_pr)+1/tand(alpha2_pr))/2)/np.pi
    Z = (C_L*c/s)**2 * sind(alpha2_pr)**2 / sind(alpha_m_pr)**3
    
    Ycl = 0.47*Z*(c/h)*(delta/c)**0.78
    
    # if __name__ == "__main__":
    #     print(f"sin_ratio : {sind(alpha2_pr)**2 / sind(alpha_m_pr)**3}")
    #     print(f"cos_ratio : {cosd(alpha2_pr)**2 / cosd(alpha_m_pr)**3}")
    #     print(f"Ycl : {Ycl}")
    #     print(f"Z : {Z}")
    #     print(f"C_L : {C_L}")
    #     print(f"c/s : {c/s}")
    #     print(f"---------------------")
    
    "1.4 Trailing Edge Losses"
    DPt = 1/2 * rho2 * w2**2 * (s*np.sin(beta_g)/(s*np.sin(beta_g)-t_TE)-1)**2
    
    Yte = 2*DPt/(rho2*w2**2)    
    
    "1.5 Supersonic Expansion Losses"
    if M2_pr <= 1:
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
    
    Y_dict = {'Y_tot': Y,
              'Yp'   : Yp,
              'Ys'   : Ys,
              'Ycl'  : Ycl,
              'Yte'  : Yte,
              'Yex'  : Yex,
              'Ysh'  : Ysh,
              'Ylw'  : Ylw
              }
        
    return Y_dict

if __name__ == "__main__":
    
    from CoolProp.CoolProp import PropsSI
    
    study_case = 'Test'
    
    if study_case == 'CO2':    
        
        # States
        
        fluid = 'CO2'
        p_in_stator = 15*1e5 # Pa
        T_in_stator = 773 # K
        
        p_out_stator = 15*1e5/1.55 # Pa
        T_out_stator = 378.5 # K
        
        rho_out_stat, mu_out_stat, a_out = PropsSI(('D','V','A'),'P', p_out_stator, 'T', T_out_stator, fluid)
        a_in = PropsSI(('A'),'P', p_in_stator, 'T', T_in_stator, fluid)
        
        # Speed triangles
        
        R = 0.5
        psi = 1
        phi = 0.5
        
        vu2OverU = (2*(1-R) + psi)/2
        wu2OverU  = vu2OverU - 1
        
        vu3OverU = (2*(1-R) - psi)/2
        wu3OverU  = vu3OverU - 1
        
        beta1 = beta3 = np.arctan(wu3OverU/phi)
        beta2 = np.arctan(wu2OverU/phi)
        
        alpha1 = alpha3 = np.arctan(vu3OverU/phi)
        alpha2 = np.arctan(vu2OverU/phi)
        
        # General params
        omega = 3000 # RPM
        r_m = 209.8 * 1e-3 # m
        
        omega_rad = omega*(2*np.pi)/60
        u = r_m*omega_rad
        
        vm = phi*u
        wu2 = wu2OverU*u
        wu3 = wu3OverU*u
        vu2 = vu2OverU*u
        vu3 = vu3OverU*u
        
        v2 = np.sqrt(vu2**2 + vm**2)
        w2 = np.sqrt(wu2**2 + vm**2)
    
        wu2 = wu2OverU*u
        vu2 = vu2OverU*u
        
        v3 = np.sqrt(vu3**2 + vm**2)
        w3 = np.sqrt(wu3**2 + vm**2)
        
        # Geometry
        
        h_rot = 29.8*1e-3 # m
        c_rot = 18.15*1e-3 # m
        xhi_rot = alpha2 # rad
        t_max_rot = 3.55*1e-3 # m
        
        M2_pr_rot = w2/a_out
        M3_pr_rot = w3/a_out
        
        s_rot = 16.5*1e-3 # m
        
        R_c = c_rot/(2*np.sin((alpha2 - alpha3)))
        
        # Params
        
        D_lw = 0
        N_lw = 0
        t_TE = 0.4*1e-3 # m
        e = 0.002*1e-3 # m : roughness
        delta_tip_rot = 0*1e-3 # m : rotor tip clearance 
            
        beta_g = 0.5*(alpha1 + alpha2)   # mid-passage metal angle
        beta_gr = 0.5*(beta2 + beta3)   # mid-passage metal angle
        
        # results_stator
        results_rotor = aungier_loss_model(beta2, beta3, beta_g, alpha2, c_rot, delta_tip_rot, D_lw, 
                                           e, h_rot, mu_out_stat, M2_pr_rot, M3_pr_rot, N_lw, R_c, 
                                           rho_out_stat, s_rot, t_max_rot, t_TE, vm, w3, 1)
        
        
        # stage.Y_vec_R = aungier_loss_model(self.Vel_Tri['beta2'], self.Vel_Tri['beta3'], self.Vel_Tri['beta2'], stage.chord_R, 
        #                        self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.static_states['V'][3], 
        #                        stage.M2_R, stage.M3_R, self.params['N_lw'], R_c, stage.static_states['D'][3], stage.pitch_R, stage.t_blade_R, self.params['t_TE'],
        #                        self.Vel_Tri['vm'], self.Vel_Tri['w3'],1)
    
    if study_case == 'ORC':    
        
        # States
        
        fluid = 'R245fa'
        p_in_stator = 12.6*1e5 # Pa
        T_in_stator = 373.7 # K
        
        mdot = 50
        
        Wdot = 1.52*1e6 # W
        
        PR = 6.4
        p_out_stator = p_in_stator/PR # Pa
        
        h_in_stator = PropsSI('H', 'T', T_in_stator, 'P', p_in_stator, fluid)        
        h_out_stator = h_in_stator
        
        T_out_stator = PropsSI('T', 'H', h_out_stator, 'P', p_out_stator, fluid)  
        
        rho_out_stat, mu_out_stat, a_out = PropsSI(('D','V','A'),'P', p_out_stator, 'T', T_out_stator, fluid)
        a_in = PropsSI(('A'),'P', p_in_stator, 'T', T_in_stator, fluid)
        
        h_in_rotor = h_out_stator
        h_out_rotor = h_in_rotor - Wdot/mdot
        
        p_out_rotor = p_out_stator

        T_out_rotor = PropsSI('T', 'H', h_out_rotor, 'P', p_out_rotor, fluid)  
        
        rho_out_rot, mu_out_rot, a_out = PropsSI(('D','V','A'),'P', p_out_rotor, 'T', T_out_rotor, fluid)

        
        "Speed triangles"
        
        # Stator
        alpha1 = 3
        alpha2 = 71 # [°]
        
        khi1 = alpha1

        beta2 = -24

        v2 = 202
        w2 = 71 # [m/s]
        vm = 65
        
        M1 = 0.53
        M2 = 1.47

        # Rotor
        
        alpha3 = 3 # [°]
        beta3 = 68 # [°]
        
        v3 = 65 # [m/s]
        w3 = 172 # [m/s]
        
        M2r = 0.52
        M3 = 1.23
        
        "General params"
        omega = 10000 # RPM
        
        cs = 44.34*1e-3 # m
        cr = 53.26*1e-3 # m
        
        e = 0.002*1e-3
        
        ZW = 0.8
                
        solidityStator = 2*cosd(alpha2)/cosd(alpha1) * sind(alpha2-alpha1)/ZW
        
        ss = (cs/solidityStator)/2
        
        hs = cs*0.43
        hr = cr*0.9
        
        solidityRotor = 2*cosd(beta3)/cosd(beta2) * sind(beta3-beta2)/ZW
        sr = (cr/solidityRotor)/2
        
        d = 0.31
        tip_cl = 0.0177*hr# max(0.001, 5*1e-4*d)
        
        "Assumptions"
        
        R_c = cs/(2*abs(np.sin((alpha1 - alpha2))))
        
        o_s = ss*sind(alpha2)
        o_r = ss*sind(beta3)
        
        t_TE_s = 0.0806*o_s
        t_TE_r = 0.0575*o_r
        
        t_max = 0.08*cs
        t_maxr = 0.3*cr
        
        beta_g = alpha2 #0.5*(alpha1 + alpha2)   
        beta_gr = beta3 # 0.5*(beta2 + beta3)   
                
        YdictS = aungier_loss_model(alpha1, alpha2, beta_g, khi1, cs, 0, 0, e, hs, mu_out_stat, M1, M2, 0, R_c, rho_out_stat, ss, t_max, t_TE_s, vm, v2, 0)

        YdictR = aungier_loss_model(beta2, beta3, beta_gr, beta2, cr, tip_cl, 0, e, hr, mu_out_rot, M2r, M3, 0, R_c, rho_out_rot, sr, t_maxr, t_TE_r, vm, w3, 0)
        
        ref_value_S = {
            'Y_tot' : 0.191,
            'Yp' : 0.013,
            'Ys' : 0.055,
            'Yte' : 0.008,
            'Ysh' : 0.013,
            'Yex' : 0.101,
            'Ycl' : 0
            }
        
        ref_value_R = {
            'Y_tot' : 0.253,
            'Yp' : 0.016,
            'Ys' : 0.05,
            'Yte' : 0.004,
            'Ysh' : 0.011,
            'Yex' : 0.034,
            'Ycl' : 0.139
            }

        def relative_errors(computed, reference):
            errors = {}
            for key, ref_val in reference.items():
                comp_val = computed.get(key, None)
                if comp_val is None:
                    errors[key] = None
                elif ref_val != 0:
                    errors[key] = (comp_val - ref_val) / ref_val
                else:
                    errors[key] = float('inf') if comp_val != 0 else 0.0
            return errors
        
        # Example usage:
        err_S = relative_errors(YdictS, ref_value_S)
        err_R = relative_errors(YdictR, ref_value_R)
        
        print("Relative errors - Stator:")
        for k, v in err_S.items():
            print(f"{k}: {v:.2%}" if v is not None else f"{k}: N/A")
        
        print("\nRelative errors - Rotor:")
        for k, v in err_R.items():
            print(f"{k}: {v:.2%}" if v is not None else f"{k}: N/A")
 
    if study_case == 'Test':    
    
        # States
        
        fluid = 'CO2'
        p_in_stator = 24679496 # Pa
        T_in_stator = PropsSI('T', 'H', 1159373, 'P', p_in_stator, fluid)  # K
        
        mdot = 50
        
        Wdot = 1.52*1e6 # W
        p_out_stator = 22409929 # Pa
        
        h_in_stator = PropsSI('H', 'T', T_in_stator, 'P', p_in_stator, fluid)        
        h_out_stator = h_in_stator
        
        T_out_stator = PropsSI('T', 'H', h_out_stator, 'P', p_out_stator, fluid)  
        
        rho_out_stat, mu_out_stat, a_out = PropsSI(('D','V','A'),'P', p_out_stator, 'T', T_out_stator, fluid)
        a_in = PropsSI(('A'),'P', p_in_stator, 'T', T_in_stator, fluid)
        
        h_in_rotor = h_out_stator
        h_out_rotor = h_in_rotor - Wdot/mdot
        
        p_out_rotor = p_out_stator
    
        T_out_rotor = PropsSI('T', 'H', h_out_rotor, 'P', p_out_rotor, fluid)  
        
        rho_out_rot, mu_out_rot, a_out = PropsSI(('D','V','A'),'P', p_out_rotor, 'T', T_out_rotor, fluid)
    
        
        "Speed triangles"
        
        # Stator
        alpha1 = -0.691573555766352
        alpha2 = 1.205642036706006
        
        khi1 = alpha1
    
        beta2 = 0.8130386970159275
    
        v2 = 192.77
        w2 = 100.16 # [m/s]
        vm = 68.84
        
        M1 = 0.18245508233270885
        M2 = 0.3991590667470393
    
        # Rotor
        
        alpha3 = -0.69157355576635 # [rad]
        beta3 = -1.1740541019416815 # [rad]
        
        v3 = 89.4 # [m/s]
        w3 = 178.14 # [m/s]
        
        M2r = 0.20738868745930328
        M3 = 0.3729579990457117
        
        "General params"
        omega = 3000 # RPM
        
        cs = 0.02654026363206464 # m
        cr = 0.02506353058148693 # m
        
        e = 0.002*1e-3
        
        ZW = 0.8
                
        solidityStator = 2*cosd(alpha2)/cosd(alpha1) * sind(alpha2-alpha1)/ZW
        
        ss = (cs/solidityStator)/2
        
        hs = 0.019044317415287667
        hr = 0.017786116071656515
        
        solidityRotor = 2*cosd(beta3)/cosd(beta2) * sind(beta3-beta2)/ZW
        sr = (cr/solidityRotor)/2
        
        tip_cl = 0.4*1e-3
        
        "Assumptions"
        
        R_c = cs/(2*abs(np.sin((alpha1 - alpha2))))
        
        o_s = ss*sind(alpha2)
        o_r = ss*sind(beta3)
        
        t_TE_s = 0.0005
        t_TE_r = 0.0005
        
        t_max = 0.0010331685787306332
        t_maxr = 0.0012531765290743466
        
        beta_g = alpha2 #0.5*(alpha1 + alpha2)   
        beta_gr = beta3 # 0.5*(beta2 + beta3)   
                
        YdictS = aungier_loss_model(alpha1, alpha2, beta_g, khi1, cs, 0, 0, e, hs, mu_out_stat, M1, M2, 0, R_c, rho_out_stat, ss, t_max, t_TE_s, vm, v2, 1)
    
        YdictR = aungier_loss_model(-beta2, -beta3, beta_gr, -beta2, cr, tip_cl, 0, e, hr, mu_out_rot, M2r, M3, 0, R_c, rho_out_rot, sr, t_maxr, t_TE_r, vm, w3, 1)
        