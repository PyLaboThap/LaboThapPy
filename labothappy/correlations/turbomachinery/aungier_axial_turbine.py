# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 21:00:27 2025

@author: basil
"""

import numpy as np

def profile_losses(alpha1_pr, alpha2_pr, beta1, c, M2_pr, R_c, s, t_max):#
    """
    alpha_1_pr : inlet flow angle with respect to tangent in the rotating coordinate system [rad] = beta_1 ? 
    alpha_2_pr : outlet flow angle with respect to tangent in the rotating coordinate system [rad] = beta_2 ? 
    beta_1     : inlet blade angle with respect to tangent [rad] = khi 2 ? 
    c          : blade cord [m]
    M2_pr      : outlet Mach number in the rotating coordinate system (using w2) [-]
    R_c        : suction surface radius of curvature [m]
    s          : blade pitch [m]
    t_max      : blade maximum thickness [m]
    """

    # alpha_1_pr = (180/np.pi)*alpha_1_pr # [rad] to [°]    
    # alpha_2_pr = (180/np.pi)*alpha_2_pr # [rad] to [°]
    # beta1 = (180/np.pi)*beta1 # [rad] to [°]

    xi = (90 - beta1)/(90 - alpha2_pr)

    # 1 Experience factor (Kacker-Okapuu)
    Kmod = 1
    
    "2) Off-design incidence correction"
    # i = beta1 - alpha1_pr
    
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
    
    # 3 Mach number correction
    if M2_pr <= 0.6:
        Km = 1
    else:
        Km = 1 + (1.65*(M2_pr-0.6) + 240*(M2_pr-0.06)**2)*(s/R_c)**(3*M2_pr - 0.6)
    
    # 4 Compressibility correction
    Kp = 1
    
    # 5 Reynolds number correction
    Kre = 1
    
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
    
    # 8 DY_TE
    DY_TE = 0
    
    "9) Total"
    fact_1 = Yp1 + xi**2*(Yp2 - Yp1)
    fact_2 = (5*t_max/c)**xi 
    
    Yp = Kmod * Kinc * Km * Kp * Kre * (fact_1 * fact_2 - DY_TE)    

    return Dis

s_vec = np.linspace(0.4, 1, 20)
a2_vec = np.array([30,40,50])

Dis = np.zeros([len(s_vec), len(a2_vec)])

for i in range(len(s_vec)):
    for j in range(len(a2_vec)):
        Dis[i][j]= profile_losses(a2_vec[j], s_vec[i])

import matplotlib.pyplot as plt

# Plotting
plt.figure(figsize=(8, 6))
for j, a2 in enumerate(a2_vec):
    plt.plot(s_vec, Dis[:, j], label=f'α₂ = {a2}°')

plt.xlabel('Input (xi or s/c)', fontsize=12)
plt.ylabel('Profile Loss Yₚ₁', fontsize=12)
plt.title('Profile Loss vs Input for Different α₂', fontsize=14)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend(title='α₂ Values')
plt.tight_layout()
plt.show()

def aungier_loss_model():
    
    #%% "1) Total Pressure loss"
    
    "1.1 Profile Losses"
    
    Yp = profile_losses()
    
    "1.2 Secondary Losses"
    Ys = 0
    
    "1.3 Blade Clearance Losses"
    Ycl = 0
    
    "1.4 Trailing Edge Losses"
    Yte = 0
    
    "1.5 Supersonic Expansion Losses"
    Yex = 0
    
    "1.6 Shock Loss"
    Ysh = 0
    
    "1.7 Lashing Wire Losses"
    Ylw = 0
    
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




