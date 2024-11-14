# =============================================================================
" ---------------- HEAT TRANSFER COEFFICIENT CORRELATIONS ------------------- "
# =============================================================================
" ------------------------------ Import Modules ----------------------------- "
# =============================================================================
from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import curve_fit, fsolve
from math import log10, sqrt, inf
import warnings
import matplotlib.pyplot as plt
" -- Constants values: "
g = 9.81                         # m s^-2
g_lb = 32.2

''' Module where all the heat transffer coefficient correlations are depicts,
    this module take account two type of heat exchanger:
        1.- Plate Heat Exchanger (PTH)
        2.- Finned Tube Heat Exchanger (FTHX)
    
    Each correlation are based on the literature, and into each one is possible
    find the corresponding reference.
    
    Nomenclature:
        
         * Type of functions: HXType_Transformation_Author_year
             - HXType: PHX or FTHX
             - Transformation:   1PH -> Single Phase
                                 CD -> Condensation
                                 EV -> Evaporation or boiling
            - Author: Corresponds to the name of who created the corelation
            - Year: Year of the literature

'''
#%%
# =============================================================================
" -------------------------------- Historial ---------------------------------"
# =============================================================================    
'''
    * 2024-05-16: Martin Holger correlation: Checked
                : Is added the correlation for the computaion of k and PR
                for the R1233ZD(E).
    * 2024-05-17: The Shah correlation for condensation was modified:
                G: corresponds to the TOTAL mass flow for then should be higher
                than 2.3 and lower than 165.
                    - An update is carried out based on the definition of Logo
                    et al. 2015.
                
                : Specific volume was wrong written.
                
                * Doubt according to the heat transffer coeficcient.
                
                Amalfi correlation:
                    * The ratio_rho was modified as the rho_x/rho_l
    *2024-05-21: Comparation with the correlation proposed by Basil.
    
    *2024-07-09: Integration of a new correaltion for R1233ZD(E) for the
                 thermal conductivity
                          
'''   
#%%
# =============================================================================
" ---------------------- Plate Heat Exchanger Correlations------------------- "
# =============================================================================

# =============================================================================
" -- PHX Single Phase  ------------------------------------------------------ "
# =============================================================================
" Martin Holger: Correlation from the VDI Atlas "
def PHX_1PH_Martin_VDI (D_hyd = 0.00321,
                        L_p = 470e-3  , 
                        B_p = 119e-3 ,
                        b = 0.00184,
                        Beta = 86, 
                        M_dot_r = 0.5,
                        P_mean = 8e+05,
                        T_mean = 273.15 + 70,
                        N_ch = 86,
                        refrigerant = 'R1234ZE',
                        state = 'liq'):

    " Thermodynamics properties:"
    mu =  PropsSI('V','P', P_mean,'T',T_mean, refrigerant)          # Viscosity, Pa s
    rho = PropsSI('D','P', P_mean,'T',T_mean, refrigerant)          # Density, kg m^-3
    cp = PropsSI('C','P', P_mean,'T',T_mean, refrigerant)           # Specific Heat, J/kg-K
    
    if refrigerant == 'R1233ZD(E)':
        """ Richard A. Perkins and Marcia L. Huber
    
            "Measurement and Correlation of the Thermal Conductivity of trans-1-Chloro-3,3,3-trifluoropropene (R1233zd(E))"
            J. Chem. Eng. Data 2017, 62, 2659-2665
            Note that the critical enhancement contribution is not implemented.
    
            Same source as EES
           """
        # Dilute gas thermal conductivity
        T_c = 382.52 #[K]
        A_0 = -0.0103589
        A_1 = 0.0308929
        A_2 = 0.000230348
    
        k_0 = A_0 + A_1*(T_mean/T_c) + A_2*(T_mean/T_c)**2
    
        # Residual gas thermal conductivity
        try:
            Rho = PropsSI('D', 'T', T_mean, 'P', P_mean, 'R1233zd(E)')
        except:
            Rho = PropsSI('D', 'T', T_mean, 'Q', 0, 'R1233zd(E)')
        Rho_c = 489.24 #[kg/m^3]
    
        B_11 = -0.0428296
        B_12 = 0.0434288
        B_21 = 0.0927099
        B_22 = -0.0605844
        B_31 = -0.0702107
        B_32 = 0.0440187
        B_41 = 0.0249708
        B_42 = -0.0155082
        B_51 = -0.00301838
        B_52 = 0.0021019
    
        Delta_kr = (B_11 + B_12*(T_mean/T_c))*Rho/Rho_c + (B_21 + B_22*(T_mean/T_c))*(Rho/Rho_c)**2 + (B_31 + B_32*(T_mean/T_c))*(Rho/Rho_c)**3 + (B_41 + B_42*(T_mean/T_c))*(Rho/Rho_c)**4 + (B_51 + B_52*(T_mean/T_c))*(Rho/Rho_c)**5
    
        k = k_0 + Delta_kr
        Pr = cp*mu/k
            
    else:
        k  =  PropsSI('L','P', P_mean,'T',T_mean, refrigerant)          # Thermal conductivity, W m^-1 K^-1
        Pr = PropsSI('PRANDTL', 'P', P_mean,'T',T_mean, refrigerant)    # Number of Prandtl, -

    " Mass flow rate per plate"
    M_dot_r_ch = M_dot_r/N_ch                                       # Mass flow rate in the chanels, kg s^-1
    
    " Mass velocity for chanel CHEK!!!!"
    w_ch = M_dot_r_ch/(b*B_p*rho)                                   # Mass velocity in chanels, kg m^2 s^-1
    Re = rho*w_ch*D_hyd/mu                                          # Reynolds Number, -

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
    beta_r = Beta*np.pi/180                                         # Chevron angle, in Radians   
     
    " Friction factor [13]"    
    f = (np.cos(beta_r)/np.sqrt(b*np.tan(beta_r) + c*np.sin(beta_r) + xhi_0/np.cos(beta_r)) +(1 - np.cos(beta_r))/np.sqrt(xhi_1))**(-2)
 
    " Hagen number "
    Hg = f*Re**2/2
    
    " Pressure Drop:"    
    DeltaP = Hg * (mu**2*L_p)/(rho*D_hyd**3)
    
    " Extracted from the comparison with Heavear et al. [10]"
    c_q = 0.122
    q = 0.374
    
    " Wall temperature "
    T_wall = T_mean-10                                                # Wall Temperature, K
    if T_wall<273.15:
        T_wall = 274.15
        
    mu_w =  PropsSI('V','P', P_mean,'T',T_wall, refrigerant)          # Viscosity at wall T, Pa s
    
    " Nusslet number: "
    Nu = c_q*Pr**(1/3)*(mu/mu_w)**(1/6)*(2*Hg*np.sin(2*beta_r))**q

    " Heat Transffer Coefficient [W m^-2]:"
    hcv = Nu*k/D_hyd 
    
    return hcv, DeltaP

#%%
# =============================================================================
" -- Gienlinski - VDI -page 696 --------------------------------------------- " 
# =============================================================================

def Gnielinski_Pipe_HTC(mu, Pr, k, G, Dh, L):
    #-------------------------------------------------------------------------
    def Gnielinski_Laminar(Re, Pr, Dh, L):
        Nu_1 = 4.364
        Nu_2 = 1.953*(Re*Pr*Dh/L)**(1/3)
        Nu = (Nu_1**3 + 0.6**3 + (Nu_2 - 0.6)**3)**(1/3)
        return Nu
    def Gnielinski_Turbulent(Re, Pr):
        f = (1.8*log10(Re) - 1.5)**(-2)
        Nu = ((f/8)*(Re-1000)*Pr)/(1+12.7*sqrt(f/8)*(Pr**(2/3)-1))
        return Nu
    #-------------------------------------------------------------------------
    Re_min = 0
    Re_max = 1e06
    Re = G*Dh/mu
    #-------------------------------------------------------------------------
    if Re > 1e4: #fully turbulent
        Pr_min = 0.1
        Pr_max = 1000
        Nu = Gnielinski_Turbulent(Re, Pr)
    elif Re < 2300: #fully laminar
        Pr_min = 0.6
        Pr_max = inf
        Nu = Gnielinski_Laminar(Re, Pr, Dh, L)
    else: #transition zone
        Pr_min = 0.1
        Pr_max = 1000
        gamma = (Re - 2300)/(1e4 - 2300)
        Nu_lam2300 = Gnielinski_Laminar(2300, Pr, Dh, L)
        Nu_turb10000 = Gnielinski_Turbulent(1e4, Pr)
        Nu = (1-gamma)*Nu_lam2300 + gamma*Nu_turb10000
    #-------------------------------------------------------------------------
    hConv = Nu*k/Dh
    #-------------------------------------------------------------------------
    if Re >= Re_max or Pr <=Re_min:
        warnings.warn('Gnielinski singe-phase: Reynolds Out of validity range !!!')
    if Pr >= Pr_max or Pr <= Pr_min:
        warnings.warn('Gnielinski singe-phase: Prandtl Out of validity range  !!!')
    #-------------------------------------------------------------------------
    return hConv, Nu

#%%
def water_Plate_HTC(mu, Pr, k, G, Dh):
    """
    Calibrated heat transfer coeffcient correlation for water side
    
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
        print("Reynolds number for water out of bounds.")
        return 0
        
    # HTC
    h_conv = (k/Dh)*Nu

    return h_conv

def Simple_Plate_HTC(mu, Pr, k, G, Dh):
    # Reynolds number
    Re = G*Dh/mu
    
    if Re < 5*1e5:
        Nu = 0.3387*Re**(1/2)*Pr**(1/3)/(1+(0.0468/Pr)**(2/3))**(1/4)
    else:
        Nu = 0.0296*Re**(4/5)*Pr**(1/3)
        
    h_conv = (k/Dh)*Nu
    
    return h_conv

def Muley_Manglik_BPHEX_HTC(mu, mu_w, Pr, k, G, Dh, chevron_angle):
    # Reynolds number
    Re = G*Dh/mu
    
    beta = 180*chevron_angle/np.pi
    
    C = 0.2668 - 0.006967*beta + 7.244*1e-5*beta**2
    C_2 = Re**(0.728 + 0.0543*np.sin((2*np.pi*beta/90) + 3.7))
    
    Nu = C * C_2 * Pr**(1/3) * (mu/mu_w)**(0.14)
        
    h_conv = (k/Dh)*Nu
    
    return h_conv

def Martin_BPHEX_HTC(mu, mu_w, Pr, k, G, Dh, chevron_angle):
    "Martin Holger: Correlation from the VDI Atlas"
 
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

def Bentao_single(b, W, N, D_hyd, t_mean_r, P_mean, t_mean_s, refrigerant, m_dot_r):
    
    T_wall=(t_mean_s+t_mean_r)/2
    mu_wall = PropsSI('V','P', P_mean, 'T', T_wall, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    v_mass_r = m_dot_r/(b*W)/(N-1)*2    
    
    " Thermodynamics properties:"
    mu =  PropsSI('V','P', P_mean,'T',T_mean, refrigerant)          # Viscosity, Pa s
    rho = PropsSI('D','P', P_mean,'T',T_mean, refrigerant)          # Density, kg m^-3
    cp = PropsSI('C','P', P_mean,'T',T_mean, refrigerant)           # Specific Heat, J/kg-K
    k  =  PropsSI('L','P', P_mean,'T',T_mean, refrigerant)          # Thermal conductivity, W m^-1 K^-1
    Pr = PropsSI('PRANDTL', 'P', P_mean,'T',T_mean, refrigerant)    # Number of Prandtl, -
        
    nu_r = mu / rho
    v_r = v_mass_r /rho
    Re = v_r* D_hyd /nu_r    
    h = 0.2092*(k/D_hyd)*Re**0.78*Pr**(1/3)*(mu/mu_wall)**0.14
    
    return h
    
#%%
if __name__ == "__main__":
    
    " Common inpunts:"
    refrigerant = 'R134a'
    D_hyd = 0.00321
    L_p = 470e-3
    B_p = 119e-3
    b = 0.00184
    M_dot_r = 0.35
    Beta = 60
    N_ch = 40
    state = 'liq'
    T_sc = 10
    P_mean = 18e5
    T_sat = PropsSI('T', 'Q',0, 'P', P_mean, refrigerant)
    T_mean = T_sat - T_sc
    print(' --------------------------------------------------')
    print(' |                SINGLE PHASE                    |')
    print(' --------------------------------------------------')
    print(f' Fluid: {refrigerant}')
    print(f' D_hyd: {D_hyd} [m]')    
    print(f' Port-port centerline distance (large) L_p: {L_p} [m]')   
    print(f' Port-port centerline distance (Width) B_p: {B_p} [m]')   
    print(f' Corrugation height or amplitude b : {b} [m]')   
    print(f' Mass flow rate: {M_dot_r} [kg/s]')   
    print(f' Chevron angle: {Beta} [°]')   
    print(f' Channel number N_ch: {N_ch} [-]')   
    print(f' Mean Pressure P_mean: {P_mean/1e5} [bar]')   
    print(f' Mean Temeprature T_mean N_ch: {T_mean-273.15} [°C]')       
    
    " Martin Holger correlation"
    hcv_mh, DeltaP_hg = PHX_1PH_Martin_VDI (D_hyd = D_hyd,
                            L_p = L_p , 
                            B_p = B_p,
                            b = b,
                            Beta = Beta, 
                            M_dot_r = M_dot_r,
                            P_mean = P_mean,
                            T_mean = T_mean,
                            N_ch = N_ch,
                            refrigerant = refrigerant,
                            state = state)
 
    " Genieliski"
    " Thermodynamics properties:"
    mu =  PropsSI('V','P', P_mean,'T',T_mean, refrigerant)          # Viscosity, Pa s
    rho = PropsSI('D','P', P_mean,'T',T_mean, refrigerant)          # Density, kg m^-3
    cp = PropsSI('C','P', P_mean,'T',T_mean, refrigerant)           # Specific Heat, J/kg-K
    k  =  PropsSI('L','P', P_mean,'T',T_mean, refrigerant)          # Thermal conductivity, W m^-1 K^-1
    Pr = PropsSI('PRANDTL', 'P', P_mean,'T',T_mean, refrigerant)    # Number of Prandtl, -

    " Mass velocity"
    G = M_dot_r/(b*B_p*N_ch)   # Which is the definition of G?
    print(f' Mass velocity G (per plate): {G } [kg/m^2-s]')    
    print(f' Reynold Re = G*d_hyd/mu: {G*D_hyd/mu}')

    hcv_gnielisnki = Gnielinski_Pipe_HTC(mu, Pr, k, G, D_hyd, L_p)

    " Simple HTC"
    hcv_s = Simple_Plate_HTC(mu, Pr, k, G, D_hyd)

    " Muley Maglik"
    hcv_m = Muley_Manglik_BPHEX_HTC(mu, mu, Pr, k, G, D_hyd, Beta*np.pi/180)

    " Holger Basil:"
    hcv_mh2 = Martin_BPHEX_HTC(mu, mu, Pr, k, G, D_hyd, Beta*np.pi/180)
    
    " Bentao_single:"
    hcv_bentap_s = Bentao_single(b, B_p, N_ch, D_hyd, T_mean, P_mean, T_mean-10, refrigerant, M_dot_r)
    
    print(' --------------------------------------------------')
    print(' | -- Results -------------------------------------')
    print(' --------------------------------------------------')
    print(f' Martin Holger correlation: {hcv_mh}')
    print(f' Gnieliski: {hcv_gnielisnki}')
    print(f' Simple HTC: {hcv_s}')
    print(f' Muley: {hcv_m}')
    print(f' Martin Holger II: {hcv_mh2}')
    print(f' Bentao: {hcv_bentap_s}')
    print(' --------------------------------------------------')
    print('\n')
    
#%%
# =============================================================================
" -- Condensation in PHX   -   SHAH 2021 -------------------------------------"                            
# =============================================================================
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
"""
def PHX_CD_Shah_2021 (D_hyd = 0.00221,                # Hydraulic Diameter, m
                      L_p = 0.55,                     # Port-port centerline distance (Width), m
                      B_p = 119e-3,                   # Port-port centerline distance (Width), m
                      b = 0.00184,                    # Corrugation height or amplitude, m
                      phi = 1.147,                     # Enlargement factor, -
                      M_dot_r = 0.5,                  # Refrigerant mass flow rate, kg s^-1
                      P_mean = 6.9e+05,                 # Condensation pressure, Pa,
                      N_ch = 86,                      # Number of plates of the fluid, -
                      refrigerant = 'R1234ZE'):        # Fluid in condensation  
    
    def PHX_CD_Shah_x(x = 0.5,                        # Quality, -
                          D_hyd = D_hyd,
                          L_p = L_p,
                          B_p = B_p,
                          b = b,
                          phi = phi,
                          M_dot_r = M_dot_r,
                          P_mean = P_mean,
                          N_ch = N_ch,
                          refrigerant = refrigerant):
        
        " Saturation conditions "
        rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)           # Vapor density Q =1, kg m^-3
        rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3
        v_g = 1/rho_g                                                   # Liquid Volume Q = 1, kg m^-3
        v_l = 1/rho_l                                                   # Liquid Volume Q = 0, kg m^-3
        mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
        cp_l = PropsSI('C','P', P_mean, 'Q', 0, refrigerant)            # Liquid Specific Heat Q = 0, J/kg-K
        T_sat = PropsSI('T','P', P_mean, 'Q', 0, refrigerant)           # Liquid Temperature Q = 0, K
        
        if refrigerant == 'R1233ZD(E)':
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
                                D_hyd,
                                L_p,
                                B_p,
                                b, 
                                phi, 
                                M_dot_r, 
                                P_mean,
                                N_ch, 
                                refrigerant)
 
    hcv_mean = np.mean(hcv)
    return hcv_mean, DeltaP_t, hcv


def Han_Cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, G, Dh, pitch_co, beta, L_v, N_cp, m_dot, D_p):
    """
    Inputs
    ------
    x : Vapor quality
    mu_l : Liquid viscosity
    k_l : Liquid thermal conductivity
    Pr_l : Liquid Prandtl Number
    rho_l : Liquid density
    rho_v : Vapor density
    G : Mass flux
    Dh : Plate spacing
    pitch_co : Corrugated pitch
    beta : Chevron angle
    L_v : Vertical length between fluid ports (268,2 mm) 
    N_cp : Number of canals
    m_dot : Flowrate
    
    Outputs
    -------
    h_cond : Condensation HTC [W/(m^2*K)]
    Nu : Nusselt Number [-]
    DP_tot : Total Pressure Drop [Pa]
    
    Reference
    ---------
    "The caracteristics of condensation in brazed plate heat exchangers with different chevron angles", Journal of the Korean Physical Society, 2003
    Han et. al
    
    """
    
    # Preliminary calculations
    theta = np.pi/2 - beta

    G_eq = G*( (1-x) + x * (rho_l/rho_v)**0.5)
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

    " Decelaration pressure:"
    v_g = 1/rho_g                                                   # Liquid Volume Q = 1, kg m^-3
    v_l = 1/rho_l                                                   # Liquid Volume Q = 0, kg m^-3
    DeltaP_a = G**2*(v_g-v_l)*(1-0)
    
    DP_tot = DP_tp + DP_port + DP_stat + DeltaP_a
    
    return h_cond, Nu#, DP_tot
#%%
def Bentao(x, b, W, L, N, D_hyd, t_mean_r, P_mean, t_mean_s, refrigerant, m_dot_r):

    T_wall=(t_mean_s+t_mean_r)/2
    mu_wall = PropsSI('V','P', P_mean, 'T', T_wall, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    v_mass_r = m_dot_r/(b*W)/(N-1)*2
    A = L*W *(N - 2)
    
    " Saturation conditions "
    rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)           # Vapor density Q =1, kg m^-3
    rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3
    v_g = 1/rho_g                                                   # Liquid Volume Q = 1, kg m^-3
    v_l = 1/rho_l                                                   # Liquid Volume Q = 0, kg m^-3
    mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    mu_g = PropsSI('V','P', P_mean, 'Q', 1, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    h_l = PropsSI('H','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3    
    h_g = PropsSI('H','P', P_mean, 'Q', 1, refrigerant)           # Liquid density Q =0, kg m^-3
    cp_l = PropsSI('C','P', P_mean, 'Q', 0, refrigerant)            # Liquid Specific Heat Q = 0, J/kg-K
    T_sat = PropsSI('T','P', P_mean, 'Q', 0, refrigerant)           # Liquid Temperature Q = 0, K
    k_l  = PropsSI('L','P', P_mean, 'Q', 0, refrigerant)            # Liquid Conductivity Q = 0, W m^-2 K^-1
    Pr_l = PropsSI('Prandtl', 'P', P_mean,'Q', 0,refrigerant)       # Liquid Prandt number Q =0, -
    
    "Heat flux:"
    Q_dot = m_dot_r*(h_g-h_l)
    q_dot = Q_dot/A    
    
    G_eq = v_mass_r * ((1 - x) + x * (rho_l / rho_g)**0.5)
    Re_eq = G_eq * D_hyd / mu_l
    Co = (rho_g / rho_l) * ((1 - x)/x)**0.8
    h_fg = h_g - h_l
    Bo = q_dot / (v_mass_r * h_fg)
    Re = v_mass_r*D_hyd/mu_l

    mu_m=mu_l*(1-x)+mu_g*x
    
    " Liquid Heat transffer"
    h_l = 0.2092*(k_l/D_hyd)*Re**0.78*Pr_l**(1/3)*(mu_m/mu_wall)**0.14
    
    " Correction Factor Quingfa Lin:"
    Nu = 4.118*Re_eq**0.4 * Pr_l**(1/3)
    h_tp = Nu *k_l/D_hyd
    
    return h_tp

#%%
if __name__ == "__main__":

    " Common inpunts:"
    refrigerant = 'R134a'
    D_hyd = 0.00321
    L_p = 470e-3
    B_p = 119e-3
    b = 0.00184
    M_dot_r = 0.35
    Beta = 60
    N_ch = 40
    D_p = 60.3e-3         # Port Diameter, m 
    p = 3e-3            # Plate pitch (b+t), m
    Lambda = 7.2e-3     # Plate corrugation wavelength, m, in kaka{c} the name corresponds to Pc       
    X = b*np.pi/Lambda                                # X number for compute the enlargement factor
    phi = 1/6*(1+np.sqrt(1+X**2)+4*np.sqrt(1+X**2/2)) # Enlargement factor (Effective area with corrug/Proyected Area Withoud Corrug) 
    
    P_mean = 18e5
    T_sat = PropsSI('T', 'Q',1, 'P', P_mean, refrigerant)
    T_mean = T_sat
    
    print(' --------------------------------------------------')
    print(' |                CONDENSATION                    |')
    print(' --------------------------------------------------')    
    print(f' Fluid: {refrigerant}')
    print(f' D_hyd: {D_hyd} [m]')    
    print(f' Port-port centerline distance (large) L_p: {L_p} [m]')   
    print(f' Port-port centerline distance (Width) B_p: {B_p} [m]')   
    print(f' Corrugation height or amplitude b : {b} [m]')   
    print(f' Plate Pitch p : {p} [m]')   
    print(f' Mass flow rate: {M_dot_r} [kg/s]')   
    print(f' Chevron angle: {Beta} [°]')   
    print(f' Channel number N_ch: {N_ch} [-]')   
    print(f' Mean Pressure P_mean: {P_mean/1e5} [bar]')   
    print(f' Mean Temeprature T_mean: {T_mean-273.15} [°C]')           
    
    " Shah Correlation:"
    hcv_shah, DeltaP, hcv_t = PHX_CD_Shah_2021 (D_hyd = D_hyd,              
                          L_p = L_p,                     
                          B_p = B_p,                  
                          b = b,                  
                          phi = phi,                   
                          M_dot_r = M_dot_r,                 
                          P_mean = P_mean,            
                          N_ch = N_ch,                    
                          refrigerant = refrigerant)
    
    
    " Hanh Correlation:"
    
    rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)           # Vapor density Q =1, kg m^-3
    rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3
    mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    cp_l = PropsSI('C','P', P_mean, 'Q', 0, refrigerant)            # Liquid Specific Heat Q = 0, J/kg-K
    T_sat = PropsSI('T','P', P_mean, 'Q', 0, refrigerant)           # Liquid Temperature Q = 0, K
    k_l  = PropsSI('L','P', P_mean, 'Q', 0, refrigerant)            # Liquid Conductivity Q = 0, W m^-2 K^-1
    Pr_l = PropsSI('Prandtl', 'P', P_mean,'Q', 0,refrigerant)       # Liquid Prandt number Q =0, -

    " Mass velocity"
    G = M_dot_r/(b*B_p*N_ch)   # Which is the definition of G?    
    print(f' Mass velocity G (per plate): {G } [kg/m^2-s]')    
    print(f' Reynold Re = G*d_hyd/mu: {G*D_hyd/mu}')
    
    " Number of subdivision "
    N_subdiv = 100
    
    " Vector of quality, x=1, problem in Re_L"
    x_vec = np.linspace(0, 0.99, N_subdiv) 
    hcv_hanh_i = np.zeros(N_subdiv)

    for i in range(len(x_vec)):
        hcv_hanh_i[i], Nu = Han_Cond_BPHEX_HTC(x_vec[i],
                                  mu_l,
                                  k_l,
                                  Pr_l,
                                  rho_l,
                                  rho_g,
                                  G,    # IMPORTANT DISCUSS WITH basil the difference of the G
                                  D_hyd,   
                                  p,
                                  Beta*np.pi/180, 
                                  L_p,
                                  N_ch,
                                  M_dot_r,                              
                                  b)
    hcv_hanh = np.mean(hcv_hanh_i)
   
    " Bentao Correlation "
    hcv_bentao_i = np.zeros(N_subdiv)
    for i in range(len(x_vec)):  
        hcv_bentao_i[i] = Bentao(x_vec[i],
                          b,
                          B_p,
                          L_p,
                          N_ch,
                          D_hyd,
                          T_mean, 
                          P_mean,
                          T_mean-10,
                          refrigerant,
                          M_dot_r)
    
    hcv_bentao = np.mean(hcv_bentao_i)
    print(' --------------------------------------------------')
    print(' | -- Results -------------------------------------')
    print(' --------------------------------------------------')
    print(f' Shah correlation: {hcv_shah}')
    print(f' Hanh: {hcv_hanh}')
    print(f' Bentao: {hcv_bentao}')
    print(' --------------------------------------------------')
    print('\n')
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$hcv_{tp}$ $\left[\frac{W}{m^{2} K}\right]$',fontsize=25)
    ax.set_xlabel(r'$x$ $[-]$',fontsize=25)
    ax.grid('on',alpha = 0.5)
    ax.tick_params(axis='both', which='major', labelsize= 25)
    ax.plot(x_vec, hcv_t, color='red', label=r'Shah et al.', marker = 'D')
    ax.plot(x_vec, hcv_hanh_i, color='blue', label=r'Hanh et al.', marker = 'D')
    ax.plot(x_vec, hcv_bentao_i, color='orange', label=r'W.S Kuo et al.', marker = 'D')
    plt.yticks(np.arange(0,25000,2500)) 
    plt.ylim(0, 22500)
    plt.xticks(np.arange(0, 1.1, 0.2))
    # plt.xlim(0, 24)
    plt.legend(prop={'size': 18}, loc='center right', ncol = 1)
    plt.tight_layout()
    # plt.savefig('Condensation.svg')
    
    
#%%
# =============================================================================
" -- Boiling in PHX   - Amalfi 2015. ---------------------------------------- "
# =============================================================================
'''    Saturated Boiling with Forced Flow - Amalfi (2016b) Correlation - Recommended by Shah
Fluids: Water; Halogenated Refrigerants (CFCs, HCFCs, and HFCs); Water-Ammonia
    beta: [27-70]°
    D_hyd : [1.7-8] mm
    Re_LT : [41-5360]
    Re_LS : [12-5320]
    rho_f/rho_g : [19-79]
    Bond N°: [1.9-79]
    Weber N°: [0.026-162]'''
    
" Amalfi correlation: "
def PHX_EV_Amalfi_2015 (D_hyd = 0.00321,                # Hydraulic Diameter, m  
                        L_p = 0.47,                     # Port-port centerline distance (Width), m
                        B_p = 119e-3,                   # Port-port centerline distance (Width), m
                        b = 0.00184,                    # Corrugation height or amplitude, m
                        Beta = 60,                      # Chevron angle, ° 
                        A_tot = 5.517,                  # Area of transffer, m^2
                        M_dot_r = 0.5,                 # Refrigerant mass flow rate, kg s^-1
                        P_mean = 4e+05,                 # Condensation pressure, Pa,
                        N_ch = 86,                      # Number of plates of the fluid, -
                        refrigerant = 'R1233ZD(E)'):        # Fluid in condensation    
  
    def PHX_EV_Amalfi(x, 
                      D_hyd = D_hyd,
                      L_p = L_p,
                      B_p = B_p,
                      b = b,
                      Beta = Beta,
                      A_tot = A_tot,
                      M_dot_r = M_dot_r,
                      P_mean = P_mean,
                      N_ch = N_ch, 
                      refrigerant = refrigerant):  
        
        " Gas properties "
        rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)
        mu_g = PropsSI('V','P', P_mean, 'Q', 1, refrigerant)
        h_g = PropsSI('H','P',P_mean,'Q',1,refrigerant)
        
        " Liquid properties "    
        rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)
        mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)
        DeltaH_lg =  PropsSI('H','P',P_mean,'Q',1,refrigerant)-PropsSI('H','P',P_mean,'Q',0,refrigerant)
        
        " Instant quality-properties "
        rho_x = PropsSI('D','P', P_mean, 'Q', x, refrigerant)  
        if refrigerant == 'R1233ZD(E)':
            T_mean = PropsSI('T','P',P_mean,'Q',1,refrigerant)
            k_x = (0.09513*T_mean - 17.963)/1000

        else:
            k_x  = PropsSI('L','P', P_mean, 'Q', x, refrigerant)       
        sigma_x = PropsSI('surface_tension', 'P', P_mean, 'Q', x, refrigerant)
        h_r_x = PropsSI('H','P',P_mean,'Q', x, refrigerant)    
        
        " Average conditions "
        x_m = 0.5*(0 + x)
        rho_m = (x_m/rho_g + (1-x_m)/rho_l)**-1
    
        " Mass velocity"
        G = M_dot_r/(2*b*B_p)

        " Heat Flux:"
        Q_dot = M_dot_r*(h_g - h_r_x)
        q_dot = Q_dot/A_tot
        
        " Boiling Number"
        Bo = q_dot/(G*DeltaH_lg)
      
        " Weber number"
        We = G**2*D_hyd/(rho_m*sigma_x)

        " Reynolds as fully Liquid"
        Re_LT = G*D_hyd/mu_l
    
        " Reynolds with only the gas fraction "                    
        Re_GS = G*x*D_hyd/mu_g             
    
        " Bond Number "
        Bd = (rho_l-rho_x)*g*D_hyd**2/sigma_x   
         
        " Normalization of the density"
        ratio_rho  = rho_x/rho_l 
        
        " Normalization of the Beta angles, Maximum 70"
        ratio_beta = Beta/70  
        
        " Boiling condition for compute the Nusslet Number:"
        if Bd<4:
            Nu_TP = 982*ratio_beta**(1.101)*We(**0.315)*Bo**(0.320)*ratio_rho**(-0.224)
        elif Bd>=4:
            Nu_TP = 18.495**ratio_beta**0.248 * Re_GS**0.135 * Re_LT**0.351 * Bd**0.235 * Bo**0.198 * ratio_rho**-0.223      
        hcv_r = Nu_TP*k_x/D_hyd    

        " Pressure"
        v_g = 1/rho_g                                                 # Liquid Volume Q = 1, kg m^-3
        v_l = 1/rho_l                                                 # Liquid Volume Q = 0, kg m^-3

        " Mass velocity per total"  
        G = M_dot_r/N_ch/(B_p*b)

        " Pressure Drop: Extracted from Longo et al, and based in his ref 16.(Collier et Thome)"
        " Decelaration pressure:"
        DeltaP_a = G**2*(v_g-v_l)*(1-0)
        
        " Gravity pressure rise (elevation):"
        DeltaP_g = g*rho_m*L_p
        
        " Manifold and ports pressure drops (Shah and Focke correlation 1988):"
        DeltaP_c = 1.5*G**2/(2*rho_m)
        
        " Friction pressure drops Longo et al."
        KEV = G**2/(2*rho_m)
        DeltaP_f = 1.95*KEV*1000   # kPa, This values depends on the refrigerant in Longo et al. he uses 2.0
        
        " Total pressure Drop:"
        DeltaP_t = DeltaP_a + DeltaP_g + DeltaP_c + DeltaP_f
        
        " Total Pressure drop:"
        DeltaP_t = DeltaP_a + DeltaP_g + DeltaP_c + DeltaP_f    
        return hcv_r, DeltaP_t
    
    " Quality subdivision "
    N_subdiv = 100
    x_vec = np.linspace(0.1, 1, N_subdiv)

    " h_cv vector "
    hcv = np.zeros(N_subdiv)
    
    " Heat transffer calculation:"
    for i in range (N_subdiv):
        hcv[i], DeltaP_t = PHX_EV_Amalfi(x_vec[i], 
                                          D_hyd,
                                          L_p,
                                          B_p,
                                          b,
                                          Beta,
                                          A_tot,
                                          M_dot_r,
                                          P_mean,
                                          N_ch, 
                                          refrigerant)
                    
    " Average calculation:"
    hcv_mean = np.mean(hcv)

    return hcv_mean, DeltaP_t, hcv
#%%
def Han_Boiling_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v,  i_fg, G, DT_log, Qdot, honv_h, Dh, theta, pitch_co):
    
    def iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        Bo_g = max(Bo_g, 1e-8)
        Nu = Ge1*Re_eq**Ge2*Bo_g**0.3*Pr_l**0.4
        h = Nu*k_l/Dh
        U = (1/h +  1/honv_h)**-1
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
    f_Bo = lambda xx: res_iter_Han_boiling(xx, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    sol = fsolve(f_Bo, Bo_0,)
    Bo_sol = sol[0]
    _, Nu, h, _, _, _, _ = iter_Han_boiling(Bo_sol, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    
    h_boiling = h
    
    return h_boiling, Nu

#%%

def Bentao_ev(b, W, L, N, D_hyd, t_mean_r, P_mean, t_mean_s, refrigerant, m_dot_r):

    T_wall=(t_mean_s+t_mean_r)/2
    mu_wall = PropsSI('V','P', P_mean, 'T', T_wall, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    v_mass_r = m_dot_r/(b*W)/(N-1)*2
    A = L*W *(N - 2)
    
    " Saturation conditions "
    rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)           # Vapor density Q =1, kg m^-3
    rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3
    v_g = 1/rho_g                                                   # Liquid Volume Q = 1, kg m^-3
    v_l = 1/rho_l                                                   # Liquid Volume Q = 0, kg m^-3
    mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    mu_g = PropsSI('V','P', P_mean, 'Q', 1, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    h_l = PropsSI('H','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3    
    h_g = PropsSI('H','P', P_mean, 'Q', 1, refrigerant)           # Liquid density Q =0, kg m^-3
    cp_l = PropsSI('C','P', P_mean, 'Q', 0, refrigerant)            # Liquid Specific Heat Q = 0, J/kg-K
    T_sat = PropsSI('T','P', P_mean, 'Q', 0, refrigerant)           # Liquid Temperature Q = 0, K
    k_l  = PropsSI('L','P', P_mean, 'Q', 0, refrigerant)            # Liquid Conductivity Q = 0, W m^-2 K^-1
    Pr_l = PropsSI('Prandtl', 'P', P_mean,'Q', 0,refrigerant)       # Liquid Prandt number Q =0, -
    
    "Heat flux:"
    Q_dot = m_dot_r*(h_g-h_l)
    q_dot = Q_dot/A    

    delta_X = Q_dot/m_dot_r/(h_g-h_l)
    X_ex = 1
    x_m = X_ex-delta_X/2        
    G_eq = v_mass_r * ((1 - x_m) + x_m * (rho_l / rho_g)**0.5)
    Re_eq = G_eq * D_hyd / mu_l
    
    h_fg = h_g - h_l
    Bo = q_dot / (v_mass_r * h_fg)
    Bo_eq = q_dot / (G_eq * h_fg)    
    Re = v_mass_r*D_hyd/mu_l

    mu_m=mu_l*(1-x_m)+mu_g*x_m
    
    " Liquid Heat transffer"
    h_l = 0.2092*(k_l/D_hyd)*Re**0.78*Pr_l**(1/3)*(mu_m/mu_wall)**0.14
    
    " Correction Factor Quingfa:"
    Nu = 2.17*Re_eq**0.495*Re**0.05*Pr_l**(1/3)
    h_tp = Nu *k_l/D_hyd  
    return h_tp

#%%

if __name__ == "__main__":
    
    " Common inpunts:"
    refrigerant = 'R134a'
    D_hyd = 0.00321
    L_p = 470e-3
    B_p = 119e-3
    b = 0.00184
    A_tot = 5.8
    M_dot_r = 0.35
    Beta = 60
    N_ch = 40
    state = 'liq'
    T_sc = 10
    P_mean = 4e5
    T_sat = PropsSI('T', 'Q',0, 'P', P_mean, refrigerant)
    T_mean = T_sat 
    
    print(' --------------------------------------------------')
    print(' |                  BOILING                       |')
    print(' --------------------------------------------------')      
    print(f' Fluid: {refrigerant}')
    print(f' D_hyd: {D_hyd} [m]')    
    print(f' Port-port centerline distance (large) L_p: {L_p} [m]')   
    print(f' Port-port centerline distance (Width) B_p: {B_p} [m]')   
    print(f' Corrugation height or amplitude b : {b} [m]')   
    print(f' Mass flow rate: {M_dot_r} [kg/s]')   
    print(f' Chevron angle: {Beta} [°]')   
    print(f' Channel number N_ch: {N_ch} [-]')   
    print(f' Mean Pressure P_mean: {P_mean/1e5} [bar]')   
    print(f' Mean Temeprature T_mean N_ch: {T_mean-273.15} [°C]')  
    
    hcv_amalfi, DeltaP_amalfi, hcv_a_t = PHX_EV_Amalfi_2015 (D_hyd = D_hyd,           # Hydraulic Diameter, m  
                            L_p = L_p,                        # Port-port centerline distance (Width), m
                            B_p = B_p,                        # Port-port centerline distance (Width), m
                            b =b,                             # Corrugation height or amplitude, m
                            Beta = Beta,                      # Chevron angle, ° 
                            A_tot = A_tot,                    # Area of transffer, m^2
                            M_dot_r = M_dot_r,                # Refrigerant mass flow rate, kg s^-1
                            P_mean = P_mean,                  # Condensation pressure, Pa,
                            N_ch = N_ch,                      # Number of plates of the fluid, -
                            refrigerant = refrigerant)        # Fluid in condensation    


    " Hanh Correlation:"
    rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)           # Vapor density Q =1, kg m^-3
    rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)           # Liquid density Q =0, kg m^-3
    mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)            # Liquid Viscosity Q = 0, Pa s^-1
    cp_l = PropsSI('C','P', P_mean, 'Q', 0, refrigerant)            # Liquid Specific Heat Q = 0, J/kg-K
    T_sat = PropsSI('T','P', P_mean, 'Q', 0, refrigerant)           # Liquid Temperature Q = 0, K
    k_l  = PropsSI('L','P', P_mean, 'Q', 0, refrigerant)            # Liquid Conductivity Q = 0, W m^-2 K^-1
    Pr_l = PropsSI('Prandtl', 'P', P_mean,'Q', 0,refrigerant)       # Liquid Prandt number Q =0, -

    " Mass velocity"
    G = M_dot_r/(b*B_p*N_ch)   # Which is the definition of G?    
    print(f' Mass velocity G (per plate): {G } [kg/m^2-s]')    
    print(f' Reynold Re = G*d_hyd/mu: {G*D_hyd/mu}')
 
    " Gas properties "
    rho_g = PropsSI('D','P', P_mean, 'Q', 1, refrigerant)
    mu_g = PropsSI('V','P', P_mean, 'Q', 1, refrigerant)
    h_g = PropsSI('H','P',P_mean,'Q',1,refrigerant)
    
    " Liquid properties "    
    rho_l = PropsSI('D','P', P_mean, 'Q', 0, refrigerant)
    mu_l = PropsSI('V','P', P_mean, 'Q', 0, refrigerant)
    h_l = PropsSI('H','P',P_mean,'Q',0,refrigerant)
    
    DeltaH_lg =  PropsSI('H','P',P_mean,'Q',1,refrigerant)-PropsSI('H','P',P_mean,'Q',0 ,refrigerant)   
    
    G = M_dot_r/(b*B_p*N_ch)   # Which is the definition of G?    
    print(f' Mass velocity G (per plate): {G } [kg/m^2-s]')    
    print(f' Reynold Re = G*d_hyd/mu: {G*D_hyd/mu}')
   
    " Number of subdivision "
    N_subdiv = 100
    
    " Vector of quality, x=1, problem in Re_L"
    x_vec = np.linspace(0, 0.99, N_subdiv) 
    hcv_hanh_b_i = np.zeros(N_subdiv)
    
    for i in range(len(x_vec)):
        
        x_i = x_vec[i]        
        h_r_x = PropsSI('H','P',P_mean,'Q', x_i, refrigerant)  
        Q_dot = M_dot_r*(h_g - h_r_x)    
         
        hcv_hanh_b_i[i], Dp_hanh = Han_Boiling_BPHEX_HTC(x_i ,
                                                mu_l,
                                                k_l,
                                                Pr_l,
                                                rho_l,
                                                rho_g, 
                                                DeltaH_lg,
                                                G,
                                                5,
                                                Q_dot,
                                                2377,
                                                D_hyd,
                                                Beta*np.pi/180,
                                                p)



    hcv_hanh_b = np.mean(hcv_hanh_b_i)
    
    
    " Bentao correlation"
    h_bentao_b = Bentao_ev(b, B_p, L_p, N_ch, D_hyd, T_mean, P_mean, T_mean-10, refrigerant, M_dot_r)
     
    print(' --------------------------------------------------')
    print(' | -- Results -------------------------------------')
    print(' --------------------------------------------------')
    print(f' Amalfi correlation: {hcv_amalfi}')
    print(f' Hanh: {hcv_hanh_b}')
    print(f' Bentao: {h_bentao_b}')
    print(' --------------------------------------------------')
    print('\n')
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$hcv_{tp}$ $\left[\frac{W}{m^{2} K}\right]$',fontsize=25)
    ax.set_xlabel(r'$x$ $[-]$',fontsize=25)
    ax.grid('on',alpha = 0.5)
    ax.tick_params(axis='both', which='major', labelsize= 25)
    ax.plot(x_vec, hcv_a_t, color='red', label=r'Amalfi et al.', marker = 'D')
    ax.plot(x_vec, hcv_hanh_b_i, color='blue', label=r'Hanh et al.', marker = 'D')
    # plt.yticks(np.arange(0,25000,2500)) 
    # plt.ylim(0, 22500)
    plt.xticks(np.arange(0, 1.1, 0.2))
    # plt.xlim(0, 24)
    plt.legend(prop={'size': 18}, loc='center right', ncol = 1)
    plt.tight_layout()
    # plt.savefig('Boiling.svg')
    
#%%
# =============================================================================
" ----------------------------- Water side -----------------------------------" 
# =============================================================================
""" Correlation for water side Magnekick et al.[]
    I should complety this part with more information
"""
def PHX_Water(B_p = 119e-3 ,                     # Port-port centerline distance (Width), m
                Beta = 86,                         # beta angles chevron   
                b = 0.00184,                       # Corrugation height or amplitude             
                M_dot_w = 5,                       # Water mass flow rate
                P_mean = 1.5e+05,                  # Water Pressure
                T_mean = 273.15 + 40,              # Water Temperature
                N_ch = 86):                        # Chanel number
    
    fluid = 'Water' 

    " Beta angle from degree to radians"
    beta_r = Beta*np.pi/180                                         # Chevron angle, in Radians
    
    " Hydraulic diameter definition"
    D_hyd = B_p
    
    " Thermodynamics properties:"
    mu =  PropsSI('V','P', P_mean,'T',T_mean, fluid)          # Viscosity, Pa s
    rho = PropsSI('D','P', P_mean,'T',T_mean, fluid)          # Density, kg m^-3
    k  =  PropsSI('L','P', P_mean,'T',T_mean, fluid)          # Thermal conductivity, W m^-1 K^-1
    Pr = PropsSI('PRANDTL', 'P', P_mean,'T',T_mean, fluid)    # Number of Prandtl, -

    " Mass velocity for chanel"
    G_ch = M_dot_w/(b*B_p*rho*N_ch)                                 # Mass velocity in chanels, kg m^2 s^-1
    Re = rho*G_ch*D_hyd/mu                                          # Reynolds Number, -
    
    " Factor for correlations: Magnekick et al."
    if Re >= 1000:                                                  # Regimen: Turbulent
        Nu = (0.2668-0.006967*beta_r+7.244e-05*beta_r**2)*Re**(0.728+0.0543*np.sin(np.pi*beta_r/45+3.71))*Pr**(1/3)*(1)**0.14
        hcv_w = Nu*k/D_hyd           
        
    elif Re <1000:                                                  # Regime: Laminar
        hcv_w = 0.277*(k/D_hyd)*Re**(2/3)*Pr**(1/3)
    
    return hcv_w

if __name__ == "__main__":
    cp = PHX_Water()
   
#%%    
#%%
#%%    
#%%
#%%    
#%%
#%%
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

