from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import curve_fit, fsolve
from math import log10, sqrt, inf
import warnings
import matplotlib.pyplot as plt

" -- Constants values: "
g = 9.81                         # m s^-2
g_lb = 32.2

# =============================================================================
" ----------------- Finned and tube Heat Exchanger Correlations-------------- "
# =============================================================================

# =============================================================================
" ----------------------- Single phase interflow ---------------------------- "
# =============================================================================
def Single_phase_internalflow(fluid,                # Refrigerant 
                              t_su = 273.15+57.43,  # Supply T
                              t_ex = 273.15+52.43,  # Exit T
                              P_su = 36.31e+05,     # Supply P
                              P_ex = 36.3e+05,      # Exit P
                              rho_su = 841.6,       # Supply density
                              rho_ex = 900,         # Exit supply
                              M_dot = 0.05,         # Refrigerant Mass flow
                              N_circuits = 2,       # Circuit number
                              L = 0.474,            # Total large of the tubes
                              D = 0.0065,           # Internal diameter
                              W = 0.4571):         
    
    """
    Convective heat transfer coeficient and pressure frop for a single phase
    inside a pipe. THe equation are separed for two regimen Laminar (<2300)
    and for turbulent.
    """
    
    " Average temperature and pressure of the air"
    t_mean = (t_su + t_ex)/2
    P_mean = (P_su + P_ex)/2
    
    " Fluid properties "
    # rho_su = PropsSI('D','P', P_su,'T',t_su, fluid)             # Density
    # rho_ex = PropsSI('D','P', P_ex,'T',t_ex, fluid)             # Density
    rho = PropsSI('D','P', P_mean,'T',t_mean, fluid)             # Density
    mu = PropsSI('V','P', P_mean,'T',t_mean, fluid)              # Visocity
    Pr = PropsSI('PRANDTL','P', P_mean,'T',t_mean, fluid)        # Prandtl number
    c_p = PropsSI('C','P', P_mean,'T',t_mean, fluid)             # Mass specific heat   
    K = PropsSI('L','P', P_mean,'T',t_mean, fluid)               # Conductivity     
    nu = mu/rho                                                  # Cinematic viscocity
    
    " Mass flow per tube "
    M_dot_tube = M_dot/N_circuits
    
    " Cross sectonal area "
    A = np.pi*D**2/4
    
    " Fluid velocity "
    C = M_dot_tube/(rho*A)

    " Mass velocity "
    G = M_dot_tube/(A)
    
    " Reynolds number "
    Re = C*D/nu
    
    " Graetz number "
    Gz = D*Re*Pr/L
    Gz_inv = 1/Gz
    
    " Nusselt number and friction factor "
    if Re<2300:
        f = 64/Re
        Nu_0 = 3.657/(np.tanh(2.264*Gz_inv**(1/3))+1.7*Gz_inv**(2/3)) + 0.0499*(1/Gz_inv)*np.tanh(Gz_inv)
        Nu = Nu_0/(np.tanh(2.432*Pr**(1/6)*Gz_inv**(1/6)))
    else:
        f = (0.79*np.log(Re)-1.64)**(-2)
        Nu = ((f/8)*(Re-1000)*Pr)/(1+12.7*(f/8)**(1/2)*(Pr**(2/3)-1))
    
    " Convective transfer coefficient "
    hcv = Nu*K/D
    
    " Friction "
    DELTAP_f = f*L/D*rho*C**2/2     
    
    " Momentum "                    # Density change
    DELTAP_m = G**2 *(1/rho_ex - 1/rho_su)
    
    " Singular "                    # Pressure loss due the singularities
    n_U = np.trunc(L/W)
    k_U = 1.23                      # Singular coefficient for elbows
    DELTAP_s =  k_U*n_U*rho*C**2/2
    
    DELTA_P = DELTAP_f +DELTAP_m + DELTAP_s
    
    " Refrigerant mass inside the section "
    M_bar_r = rho*A*L*N_circuits
    
    return hcv, DELTA_P, M_bar_r
#%%
# =============================================================================
" ------------- Two phase internal flow: Condensation  ---------------------- "
# =============================================================================
def Two_Phase_InternalFlow_CD(fluid,
                            M_dot = 0.05,
                            N_circuits = 2,
                            D = 0.0065,
                            L = 4.401,
                            W = 0.4571,
                            S_T = 0.032,
                            S_L = 0.028,
                            P_sat = 36.32e+05,
                            x_su = 1,
                            x_ex = 0):
    
    " This Escurrimient corresponds to a"
    
    " Average quality "
    x=(x_su+x_ex)/2   
    
    " Refrigerant properties "
    t_sat_xsu = PropsSI('T','P', P_sat,'Q',x_su, fluid)
    t_sat_xex = PropsSI('T','P', P_sat,'Q',x_ex, fluid)      #
    t_sat = (t_sat_xsu+t_sat_xex)/2                          # T
    rho_l =  PropsSI('D','P', P_sat,'Q',0, fluid)            # Density
    rho_g =  PropsSI('D','P', P_sat,'Q',1, fluid)            # Density
    mu_l = PropsSI('V','P', P_sat,'Q',0, fluid)              # Visocity
    mu_g = PropsSI('V','P', P_sat,'Q',1, fluid)              # Visocity
    Pr_l = PropsSI('PRANDTL','P', P_sat,'Q',0, fluid)        # Prandtl
    k_l = PropsSI('L','P', P_sat,'Q',0, fluid)               # Conductivity
    sigma = PropsSI('I','P', P_sat,'Q',0,fluid)              # Surface tension
    rho_bar = (rho_l +rho_g)/2
    
    " Homogeneous density "
    rho_tp = (x/rho_g+(1-x)/rho_l)**(-1)
    
    " Cross sectional area "
    A = np.pi/4*D**2
    
    " Mass velocity "
    G = M_dot/(N_circuits*A)
    
    " Dimensionless numbers "
    Re_lo = G*D/mu_l
    Re_go = G*D/mu_g
    Fr_tp = G**2/(g*D*rho_tp**2)
    We_tp = G**2*D/(rho_tp*sigma)
    
    " ------ Pressure drop ------ "
    " Friedel "
    f_lo_fd = 0.079/Re_lo**0.25
    f_go_fd = 0.079/Re_go**0.25
    A_1fd = (1-x)**2+x**2*(rho_l/rho_g)*(f_go_fd/f_lo_fd)
    A_2fd = x**0.78*(1-x)**0.24*(rho_l/rho_g)**0.91*(mu_g/mu_l)**0.19*(1-mu_g/mu_l)**0.7
    phi_lo_Friedel_2=(A_1fd+3.24*A_2fd/(Fr_tp**0.045*We_tp**0.035))
        
    C_lo = G/rho_l
    f_lo = f_lo_fd
    DELTAP_lo = f_lo*L/D*rho_l*C_lo**2/2
    DELTAP_f_tp = phi_lo_Friedel_2*DELTAP_lo 
    
    " Momentum "                    # Density change
    DELTAP_m_tp = G**2 *(1/rho_l - 1/rho_g)

    " Singular pressure drop "
    n_U = np.trunc(L/W)
    D_U = np.sqrt(S_T**2+S_L**2)
    L_U = np.pi/2*D_U
    lambda_U = 0.68
    C_2 = 1+20*D/L_U
    C_sh = (lambda_U+(C_2-lambda_U)*(1-rho_g/rho_l)**0.5)*((rho_l/rho_g)**0.5+(rho_g/rho_l)**0.5)
    X_tt = ((1-x)/x)**0.9*(rho_g/rho_l)**0.5*(mu_l/mu_g)**0.1
    phi_l_2 = 1+C_sh/X_tt+1/X_tt**2
    C_l = G*(1-x)/rho_l
    Re_l = C_l*D/(mu_l/rho_l)
    f_l = 0.079/Re_l**0.25  
    DELTAP_l = f_l*L_U/D*rho_l*C_l**2/2  
    DELTAP_s_tp = n_U*phi_l_2*DELTAP_l

    DELTAP_tp = DELTAP_s_tp + DELTAP_m_tp + DELTAP_f_tp
    
    " ------ Convective heat transfer coefficient ------ "
    " Thome et al. correlation "
    alpha_h = (1+((1-x)/x)*(rho_g/rho_l))**(-1)
    alpha_ra = x/rho_g*((1+0.12*(1-x))*(x/rho_g+(1-x)/rho_l)+1.18*(1-x)*(g*sigma*(rho_l-rho_g))**0.25/(G*rho_l**0.5))**(-1)
    alpha = (alpha_h-alpha_ra)/np.log(alpha_h/alpha_ra)
    A_l = (1-alpha)*A
    
    value = (D**2-A_l*8/(2*np.pi))
    
    if value>=0:
        delta = (D-np.sqrt((D**2-A_l*8/(2*np.pi))))/2
    else:
        delta = D/2

    C_l = G*(1-x)/(rho_l*(1-alpha))
    C_g = G*x/(rho_g*alpha)
    f_i = 1+(C_g/C_l)**(1/2)*((rho_l-rho_g)*g*delta**2/sigma)**(1/4)
    Re_l_Thome = 4*G*(1-x)*delta/((1-alpha)*mu_l)
    Nu_Thome = 0.003*Re_l_Thome**0.74*Pr_l**0.5*f_i
    h_Thome = Nu_Thome*k_l/delta

    " Refrigerant mass inside the section "
    M_bar_r = rho_bar*A*L*N_circuits
    
    return h_Thome, DELTAP_tp, M_bar_r
#%%
# =============================================================================
" --------------- Two phase internal flow: Boilling  ------------------------ "
# =============================================================================
def Two_Phase_InternalFlow_EV(fluid,
                            M_dot = 0.05,
                            q = 11.67,
                            N_circuits =3,
                            D = 0.0065,
                            L = 10.285714285714286,
                            W = 0.6857142857142857,
                            S_T = 0.032,
                            S_L = 0.028,
                            P_sat = 7.928e+05,
                            x_su = 0.36,
                            x_ex = 1):
    " Average quality "
    x = (x_su+x_ex)/2
    
    " Refrigerant properties "
    t_sat_xsu = PropsSI('T','P', P_sat,'Q',x_su, fluid)
    t_sat_xex = PropsSI('T','P', P_sat,'Q',x_ex, fluid)      #
    t_sat = (t_sat_xsu+t_sat_xex)/2                          # T
    rho_l =  PropsSI('D','P', P_sat,'Q',0, fluid)            # Density
    rho_g =  PropsSI('D','P', P_sat,'Q',1, fluid)            # Density
    mu_l = PropsSI('V','P', P_sat,'Q',0, fluid)              # Visocity
    mu_g = PropsSI('V','P', P_sat,'Q',1, fluid)              # Visocity
    h_l = PropsSI('H','P', P_sat,'Q',0, fluid)               # Enthalpy
    h_g = PropsSI('H','P', P_sat,'Q',1, fluid)               # Enthalpy
    Pr_l = PropsSI('PRANDTL','P', P_sat,'Q',0, fluid)        # Prandtl
    Pr_g = PropsSI('PRANDTL','P', P_sat,'Q',1, fluid)        # Prandtl
    k_l = PropsSI('L','P', P_sat,'Q',0, fluid)               # Conductivity
    k_g = PropsSI('L','P', P_sat,'Q',1, fluid)               # Conductivity
    sigma = PropsSI('I','P', P_sat,'Q',0,fluid)              # Surface tension
    
    " Vaporization enthalpy "
    h_lg = h_g-h_l
    
    " Homogeneous density "
    rho_tp = (x/rho_g+(1-x)/rho_l)**(-1)

    " Cross sectional area "
    A = np.pi/4*D**2
    
    " Mass velocity "
    G = M_dot/(N_circuits*A)

    " Dimensionless numbers "
    Re_lo = G*D/mu_l
    Re_go = G*D/mu_g
    Fr_tp = G**2/(g*D*rho_tp**2)
    We_tp = G**2*D/(rho_tp*sigma)    
    
    "  Pressure drop "
    " Friedel "
    f_lo_fd = 0.079/(Re_lo**0.25)
    f_go_fd = 0.079/(Re_go**0.25)
    A_1fd = (1-x)**2+x**2*(rho_l/rho_g)*(f_go_fd/f_lo_fd)
    A_2fd = x**0.78*(1-x)**0.24*(rho_l/rho_g)**0.91*(mu_g/mu_l)**0.19*(1-mu_g/mu_l)**0.7
    phi_lo_Friedel_2 = (A_1fd+3.24*A_2fd/(Fr_tp**0.045*We_tp**0.035))
    
    C_lo = G/rho_l
    f_lo = f_lo_fd
    DELTAP_lo = (f_lo*L/D)*rho_l*(C_lo**2)/2
    DELTAP_f_tp = phi_lo_Friedel_2*DELTAP_lo 
    
    " Momentum pressure drop "
    C_0_su = 1+0.2*(1-x_su)
    C_r_su = 1.18*(1-x_su)*(sigma*g*(rho_l-rho_g)/rho_l**2)**(0.25)
    alpha_rh_su = x_su/rho_g*((C_0_su*(x_su/rho_g+(1-x_su)/rho_l))+C_r_su/G)**(-1)
    
    DELTAP_m_tp = G**2*((1/rho_g)-((1-x_su)**2/(rho_l*(1-alpha_rh_su))+x_su**2/(rho_g*alpha_rh_su)))   
    
    " Singular pressure drop "
    n_U = np.trunc(L/W)
    D_U = np.sqrt(S_T**2+S_L**2)
    L_U = np.pi/2*D_U
    lambda_p = 0.68
    C_2 = 1+20*D/L_U
    C_sh = (lambda_p+(C_2-lambda_p)*(1-rho_g/rho_l)**0.5)*((rho_l/rho_g)**0.5+(rho_g/rho_l)**0.5)
    X_tt = ((1-x)/x)**0.9*(rho_g/rho_l)**0.5*(mu_l/mu_g)**0.1
    phi_l_2 = 1+C_sh/X_tt+1/X_tt**2
    C_l = G*(1-x)/rho_l
    Re_l = C_l*D/(mu_l/rho_l)
    f_l = 0.079/Re_l**0.25
    DELTAP_l = f_l*L_U/D*rho_l*C_l**2/2
    DELTAP_s_tp = n_U*phi_l_2*DELTAP_l
    
    " Total pressure drop "
    DELTAP_tp = DELTAP_f_tp + DELTAP_m_tp + DELTAP_s_tp
    " Shah correlation "
    Bo = q/(G*h_lg/1000)
    Fr_l = G**2/(rho_l**2*g*D)   
    if Bo>=11e-04:
        F = 14.7
    else:
        F = 15.43        
    # hcv_lo = 0.023*Re_lo**0.8*Pr_l**0.4*k_l/D
    # hcv_go = 0.023*Re_go**0.8*Pr_g**0.4*k_g/D      
    C_l = G*(1-x)/rho_l
    Re_l = C_l*D/(mu_l/rho_l)
    hcv_l = 0.023*Re_l**0.8*Pr_l**0.4*k_l/D
    Co = (1/x-1)**0.8*(rho_g/rho_l)**0.5

    if Fr_l>=0.04:
        N = Co
    else:
        N = 0.38*Fr_l**(-0.3)*Co
    psi_cb = 1.8/N**0.8
    
    if N>1.0:
        if Bo>0.3e-4:
            psi_nb = 230*Bo**0.5
        else:
            psi_nb = 1 +46*Bo**0.5   
        psi = max(psi_nb, psi_cb)
    if N>0. and N<=1:
        psi_bs = F*Bo**0.5*np.exp(2.74*N**(-0.1))
        psi = max(psi_bs,psi_cb) 
    
    hcv_Shah = psi*hcv_l
    
    " Void Fraction "
    C_0 = 1 + 0.2*(1-x)
    C_r = 1.18*(1-x)*(sigma*g*(rho_l-rho_g)/(rho_l)**2)**0.25
    alpha_rh = x/rho_g*((C_0*(x/rho_g+(1-x)/rho_l))+ C_r/G)**(-1)
    rho_rh = alpha_rh*rho_g+(1-alpha_rh)*rho_l

    " Refrigerant mass inside this section "
    M_bar_r = rho_rh*A*L*N_circuits    
    
    return DELTAP_tp, hcv_Shah, M_bar_r

if __name__ == "__main__":
    cc = Single_phase_internalflow('R410A')  
    d = Two_Phase_InternalFlow_CD('R410A')
    e = Two_Phase_InternalFlow_EV('R410A')
    
#%%
# =============================================================================
"  ----------------  Air side heat transffer coefficient -------------------- "
# =============================================================================
def Air_side_HCV(t_su,                   # Supplt T
                  t_ex = 273.15+32.31,    # Exit T
                  P_su = 101325,          # Supply P
                  P_ex = 101000,          # Exit P
                  r_h = 0.0009746,        # Hydraulic ratio
                  G = 1.914):             # Mass velocity
    
    """
    Convective heat transffer coeficient for the air Side. The correlation are
    extracted from Heat exchanger (Kays and London) [1]. That equations allows
    determinated the heat transfer un the Chilton Colburn factor. The equation
    for j-factor corresponds to an acuarracy curve from the paper
    """
    
    " Average temperature and pressure of the air"
    t_mean = (t_su + t_ex)/2
    P_mean = (P_su + P_ex)/2
    
    " Fluid properties "
    mu = PropsSI('V','P', P_mean,'T',t_mean, 'Air')       # Visocity
    c_p = PropsSI('C','P', P_mean,'T',t_mean, 'Air')      # Mass specific heat
    Pr = PropsSI('PRANDTL','P', P_mean,'T',t_mean, 'Air') # Prandtl number
    
    " Reynolds number "
    Re = 4*r_h*G/mu
    
    " Chilton colburn Factor dependence on Re [correspond an interpolation of graphic given by Kays and London] "
    j = 0.0962*Re**(-0.351)   
    
    " Staton Number "
    St = j/(Pr**(2/3))
    
    " Convective heat transfer coefficient "
    hcv = St*G*c_p    
    
    return hcv

#%%
# =============================================================================
" --------------------------- Air pressure drop  --------------------------- "
# =============================================================================
def Air_drop_pressure(t_su,                       # Supply T
                      t_ex = 273.15+32.31,        # Exit T
                      P_su = 101325,              # Supply P
                      P_ex = 101000,              # Exit P
                      r_h = 0.0009746,            # Hydraulic ratio
                      A = 18.01,                  # Area
                      A_free = 0.2089,            # Free area
                      sigma = 0.6348,             # Sigma of condenser
                      G =1.914):                  # Mass velocity
    """
    Air pressure drop, the equation are extracted from the Kays and London, 
    however the equation are presented by McQuiston in Heating Ventilind and
    Air conditioning
    
    """
    
    " Average temperature and pressure of the air"
    t_mean = (t_su + t_ex)/2
    P_mean = (P_su + P_ex)/2
    
    " Fluid properties "
    rho_su = PropsSI('D','P', P_su,'T',t_su, 'Air')          # Density
    rho_ex = PropsSI('D','P', P_ex,'T',t_ex, 'Air')          # Density
    rho_mean = (rho_su + rho_ex)/2                               # Density
    mu = PropsSI('V','P', P_mean,'T',t_mean, 'Air')              # Visocity
    c_p = PropsSI('C','P', P_mean,'T',t_mean, 'Air')             # Mass specific heat
    Pr = PropsSI('PRANDTL','P', P_mean,'T',t_mean, 'Air')        # Prandtl number
    
    " Reynolds number "
    Re = 4*r_h*G/mu       
    
    " Friction Factor "
    f = 0.0993*Re**(-0.231)
    
    " Concentration and expansion factors"
    K_c = -0.4108*sigma**2+0.0141*sigma+1.1804
    K_e = 1.003*sigma**2-2.7677*sigma+0.9978
    
    " Contraction "
    Delta_Pc = (G**2/(2*rho_su))*(1-sigma**2+K_c)
    
    " Friction"
    Delta_Pf = f*A/A_free *G**2/2*rho_mean
    
    " Momentum"
    Delta_Pm = G**2*(1/rho_ex - 1/rho_su)
    
    " Expansion"
    Delta_Pe = G**2/(2*rho_ex)*(1-sigma**2-K_e)
    
    " Air side pressure drop "
    DeltaP = Delta_Pc +Delta_Pf + Delta_Pm +Delta_Pe
    return DeltaP

#%%
# =============================================================================
" -------------- Air heat transffer coeficient McQuiston -------------------- " 
# =============================================================================

def Air_McQuiston(D, t_in, p_in, velocity, sigma, x_a, x_b, D_h):
    rho_a = PropsSI('D','P',p_in,'T', t_in,'Air' )           # kg m^-3
    V_a = velocity                                           # m s^-1
    G_a_fr = rho_a*V_a                                       # kg s^-1 m^-2
    G_a_c = G_a_fr/sigma                                     # kg s^-1 m^-2
    mu = PropsSI('V', 'T', t_in, 'P', p_in, 'Air')           # Pa s
    Re_o = G_a_c*D/mu                                        # Reynolds number in the D
    R_a = 4*x_b*x_a*sigma/(np.pi*D_h*D)                      # Ratio Area
    JP = Re_o**(-0.4) * (R_a)**(-0.15) 
    c_1 = 1.264955
    j  = (0.2675*JP +1.325e-06)*c_1                          # j_4, -0.001549676030974099*Re_o/1000 + 0.012953563726170852
    c_p = PropsSI('C','P', p_in,'T',t_in,'Air')              # J kg^-1 K
    Pr = PropsSI('PRANDTL', 'T', t_in, 'P', p_in, 'Air')     # Number of Prandt
    h_a = (G_a_c*j*c_p)/(Pr**(2/3))
    return h_a, G_a_c

#%%
# =============================================================================
" ---------- Air heat transffer coeficient Kays and London ------------------- "
# =============================================================================
def KaysandLondon(x,b0,b1):
    return b0*x + b1

j = np.array([0.015, 0.0125, 0.011, 0.009, 0.008, 0.006, 0.005])
Re_KyL = np.array([0.4, 0.7, 1, 1.5, 2, 4, 6])
res, cov = curve_fit(KaysandLondon, Re_KyL, j)
a0 = res[0]
a1 = res[1]

def OutletFin_KaysLond(t_in, p_in, velocity, sigma, D_h):
    rho_a = PropsSI('D','P',p_in,'T', t_in,'Air' )           # kg m^-3
    V_a = velocity                                           # m s^-1
    G_a_fr = rho_a*V_a                                       # kg s^-1 m^-2
    G_a_c = G_a_fr/sigma                                     # kg s^-1 m^-2
    mu = PropsSI('V', 'T', t_in, 'P', p_in, 'Air')           # Pa s
    Re_o = G_a_c*D_h/mu                                      # Reynolds number in the D
    j  =  float(a0)*Re_o/1000 + float(a1)                    
    c_p = PropsSI('C','P', p_in,'T',t_in,'Air')              # J kg^-1 K
    Pr = PropsSI('PRANDTL', 'T', t_in, 'P', p_in, 'Air')     # Number of Prandt
    h_a = (G_a_c*j*c_p)/(Pr**(2/3))

    return h_a
#%%
# =============================================================================
" ---------------- Overall Efificienci in fin surfaces --------------------- "
# =============================================================================
def overall_surface_efficiency(A_f,             # Fin area
                                A_T = 7.898,     # Total area
                                S_T = 0.032,     # Transversal distance between fins
                                S_L = 0.028,     # Longuitud distance between fins
                                h_cv = 47.05,  # heat transfer coefficient kW/m
                                e_fin = 0.0003,   # Fin thickness
                                D_i = 0.0085):    # Internal diameter of pipes
    " Connections "
    a = S_T
    b = S_L
    
    " Fin internal radius "
    r_i = D_i/2
    
    " Fin conductivity (Aluminium)"
    k_fin = 236.6/1000       

    " Fin efficiency "
    M = min(a/2, b)
    L = 0.5*np.sqrt(a**2/4 +b**2)    
    m_f = (2*h_cv/1000/(k_fin*e_fin))**0.5
    beta = L/M
    PSI = M/r_i
    r_e = r_i*(1.27*PSI*np.sqrt(beta-0.3))
    PHI = ((r_e/r_i)-1)*(1+0.35*np.log(r_e/r_i))
    eta_fin = np.tanh(m_f*r_i*PHI)/(m_f*r_i*PHI)

    " Overall surface efficiency "
    eta_0 = 1-A_f/A_T*(1-eta_fin)
    
    return eta_fin, eta_0

