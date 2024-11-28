# ============================================================================
" ----------------------------- Import Libraries ---------------------------- "
# =============================================================================

from scipy.optimize import root
import numpy as np
import sys

" Common Modules from BaseClases"
base_clases_path = r'../'
sys.path.append(base_clases_path)
from BaseClases import HTCorrelation as HTC
from BaseClases import Classes as Point
from BaseClases import VoidFraction as VF

#%%
# =============================================================================
" -------------------------------- Historial ---------------------------------"
# =============================================================================    
'''
    * 2024-05-20: Geometry verification.
                  1.- Hydraulic diameter: Ok
                  2.- Geometry correlations it is okay, is important remark that
                  here usually is used 2 a as a b
                  
                  # HTC: For the water the Martin Holger is employed
    * 2024-06-14: * Geometry verification.
                  * Numerical verification.    
                  
    * 2024-06-16: * Code verification
    
    * 2024-06-16: * Validation with the ELISE data
    
    * 2024-07-09: * Integration of the foulling coefficient
                        f_fouling = f_thickness/k_p
                
'''  
#%%
# =============================================================================
" ------------ Geometry of the Brazed Plate Heat Exchanger -------------------"
# =============================================================================     
Beta = 60             # Chevron angle, °
L_p = 593e-3          # Port-port centerline distance (Large), m
B_p = 300e-3          # Port-port centerline distance (Width), m
D_p = 150e-3          # Port Diameter, m
N_pas = 1             # Number of passes per channels
N_p = 160              # N° of Plates, - (only uses par number if impar -1)
t = 0.25e-3           # Plate thickness, m
L_c = (61.38 + 2.39*N_p)*1e-3        # Depth of the heat exchanger, m 
Lambda = 2.6e-3       # Plate corrugation wavelength, m, in kaka{c} the name corresponds to Pc
k_p = 385             # Conductivity of the plate (copper), W m^-1K^-1
V_t = 73.47/1000      # Total volume
f_fouling = 0.4578/1000  # Fouling factor
refrigerant = 'R1233ZD(E)'

# =============================================================================
" ----------------------- Geometry pre calculations ------------------------- " 
# =============================================================================
""" This predefinition are extracted from:
    - VDI Atlas: N6 Holger Martin: Pressure Drop and Heat Transfer in Plate 
      Heat exchanger
    - Kakac, 2020 " Heat Exchangers Selection, Rating, and Thermal Design, Fourth
      Edition "
    - Ian Bell, 2017, ACHP Documentation, Release 1.5 """
    
L_pv = L_p - D_p                                  # Partially vertical large, m
L_p_w = B_p + D_p                                 # Partially width large, m
p = L_c/N_p                                       # Plate pitch (b+t), m 
b = (p-t)                                         # Corrugation height or amplitude, m, Bell uses 2a, which is equivalent to b
X = b*np.pi/Lambda                                # Wavenumber
phi = 1/6*(1+np.sqrt(1+X**2)+4*np.sqrt(1+X**2/2)) # Enlargement factor
D_hyd = 2*b/phi                                   # Hydraulic diameter   
N_p_active = N_p-2                                # Total number of active plates
A_1p = B_p*L_p*phi                                # Project area of 1 plate, m^2
N_ch = (N_p)/(2)                                  # Number of channels, -
A_h = 2*N_ch*A_1p                                 # Area of the refrigerant side, m^2
A_c = 2*N_ch*A_1p                                 # Area of the Water side, m^2
A_tot = max(A_h, A_c)                             # Heat Transfer Area of the enlarged surface, m^2
A_cd = A_tot

#%%
# =============================================================================
" -------------------------- Condenser Model ------------------------------- "
# =============================================================================
" Declaration of the point: "
r_su_sh_cd = Point.Refrigerant()
r_su_tp_cd = Point.Refrigerant()
r_su_sc_cd = Point.Refrigerant()
r_ex_sh_cd = Point.Refrigerant()
r_ex_tp_cd = Point.Refrigerant()
r_ex_sc_cd = Point.Refrigerant()
r_sat_cd = Point.Refrigerant()

" Set the refrigerant in the differents poins of the cycle"
r_sat_cd.type = r_su_sh_cd.type = r_ex_sh_cd.type \
    = r_su_tp_cd.type = r_ex_tp_cd.type = r_su_sc_cd.type\
        = r_ex_sc_cd.type = refrigerant

w_su_sh_cd = Point.Water()
w_su_tp_cd = Point.Water()
w_su_sc_cd = Point.Water()
w_ex_sh_cd = Point.Water()
w_ex_tp_cd = Point.Water()
w_ex_sc_cd = Point.Water()

cd_sc = Point.Evaporator()
cd_tp = Point.Evaporator()
cd_sh = Point.Evaporator()

def ThreePhasePHX(P_r_cd_i,
                  DELTAP_sc,
                  DELTAP_tp,
                  DELTAP_sh,
                  M_dot_r_cd = 2.37,             # Refrigerant mass flow rate, kg s^-1
                  T_r_su_cd = 324.911616527858,  # Refrigerant supply Temperature, K
                  P_w_su_cd = 1e+05,             # Water pressure supply, Pa
                  T_w_su_cd = 10+273.15,         # Water Temperature supply, K
                  M_dot_w_cd = 22.5,             # Water Mass flow rate, kg s^-1
                  T_glide = False,               # Water glide, K
                  DELTAT_sc = 3                  # Subcooling  
                  ):

    # ================================================================
    " ------------------- SF-INITIAL CONDITION --------------------- "
    # ================================================================    
    " Refrigerant Properties "
    " Supply"
    r_su_sh_cd.T = T_r_su_cd
    r_su_sh_cd.P = P_r_cd_i + DELTAP_sh
    r_su_sh_cd.M_dot = M_dot_r_cd
    r_ex_sc_cd.M_dot = r_su_sc_cd.M_dot = r_ex_tp_cd.M_dot = r_su_tp_cd.M_dot = r_ex_sh_cd.M_dot = r_su_sh_cd.M_dot
    r_su_sh_cd.set_PropsPT()
    
    "Exhaust"
    r_ex_tp_cd.P = P_r_cd_i - DELTAP_tp
    r_ex_tp_cd.set_PhasePQ_0()   
    
    r_ex_sc_cd.T = r_ex_tp_cd.T_sat-DELTAT_sc
    r_ex_sc_cd.P = r_su_sh_cd.P - DELTAP_sh - DELTAP_tp- DELTAP_sc
    r_ex_sc_cd.set_PropsPT()

    " Energy Balance"    
    Q_dot_r_cd = r_su_sh_cd.M_dot*(r_su_sh_cd.H-r_ex_sc_cd.H)

    " Water connection: "
    w_su_sc_cd.T = T_w_su_cd
    w_su_sc_cd.P = P_w_su_cd
    w_su_sc_cd.set_PropsPT()
    
    " Water exhaust: "
    if T_glide == False:
        w_ex_sh_cd.P = w_su_sc_cd.P
        w_ex_sc_cd.M_dot = w_su_tp_cd.M_dot = w_ex_tp_cd.M_dot = w_su_sh_cd.M_dot = w_ex_sh_cd.M_dot = w_su_sc_cd.M_dot = M_dot_w_cd
        w_ex_sh_cd.H = Q_dot_r_cd/w_su_sc_cd.M_dot + w_su_sc_cd.H
        w_ex_sh_cd.set_PropsPH()

    else:
         w_ex_sh_cd.T = w_su_sc_cd.T + T_glide
         w_ex_sh_cd.P = w_su_sc_cd.P
         w_ex_sh_cd.set_PropsPT()
         w_su_sc_cd.M_dot = Q_dot_r_cd/(w_ex_sh_cd.H-w_su_sc_cd.H)
         w_ex_sc_cd.M_dot = w_su_tp_cd.M_dot = w_ex_tp_cd.M_dot = w_su_sh_cd.M_dot = w_ex_sh_cd.M_dot = w_su_sc_cd.M_dot

    # ================================================================
    " ----------------------- SUBCOOLED ZONE ----------------------- "
    # ================================================================
    " Refrigerant properties "
    r_su_sc_cd.P = r_ex_sc_cd.P+DELTAP_sc
    r_su_sc_cd.set_PhasePQ_0()
    cp_r_bar_sc_cd = (r_su_sc_cd.cp + r_ex_sc_cd.cp)/2
    
    " Water properties "
    w_ex_sc_cd.P = w_su_sc_cd.P
    
    " Energy balance "
    cd_sc.Q_dot = r_su_sc_cd.M_dot*(r_su_sc_cd.H - r_ex_sc_cd.H)
    w_ex_sc_cd.H = cd_sc.Q_dot/w_su_sc_cd.M_dot + w_su_sc_cd.H
    w_ex_sc_cd.set_PropsPH()
    cp_w_bar_sc_cd = (w_su_sc_cd.cp + w_ex_sc_cd.cp)/2
    
    " Capacitive flows "
    C_dot_r_sc_cd = cp_r_bar_sc_cd*r_su_sc_cd.M_dot
    C_dot_w_sc_cd = cp_w_bar_sc_cd*w_su_sc_cd.M_dot
    C_min_sc_cd = min(C_dot_r_sc_cd, C_dot_w_sc_cd)                        
    C_max_sc_cd = max(C_dot_r_sc_cd, C_dot_w_sc_cd)                          
    C_r_sc_cd = C_min_sc_cd/C_max_sc_cd                                       
    
    " Effectivness- NTU "
    cd_sc.epsilon = min(1, cd_sc.Q_dot/(C_min_sc_cd*(r_su_sc_cd.T-w_su_sc_cd.T)))
    cd_sc.NTU = 1/(C_r_sc_cd-1)*np.log((cd_sc.epsilon-1)/(cd_sc.epsilon*C_r_sc_cd-1))   
    
    " -- Heat transfer coefficients -------------------------------- "
    hcv_r_sc_cd, DELTAP_sc_i = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                    L_p = L_p , 
                    B_p = B_p,
                    b = b,
                    Beta = Beta, 
                    M_dot_r = M_dot_r_cd,
                    P_mean = (r_su_sc_cd.P + r_ex_sc_cd.P)/2,
                    T_mean = (r_ex_sc_cd.T + r_su_sc_cd.T)/2,
                    N_ch = N_ch,
                    refrigerant = refrigerant)
       
    hcv_w_sc_cd, DELTAP_sc_w_cd = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                L_p = L_p , 
                B_p = B_p,
                b = b,
                Beta = Beta, 
                M_dot_r = w_su_sc_cd.M_dot,  
                P_mean = (w_su_sc_cd.P + w_ex_sc_cd.P)/2,
                T_mean = (w_su_sc_cd.T + w_ex_sc_cd.T)/2,
                N_ch = N_ch,
                refrigerant = 'Water')
   
    " -- Global Heat transfer coefficient  "
    R_sc_cd = 1/(hcv_w_sc_cd) + t/k_p + f_fouling + 1/(hcv_r_sc_cd)
    U_sc_cd = 1/R_sc_cd
    cd_sc.A = cd_sc.NTU*C_min_sc_cd/(U_sc_cd)
    cd_sc.alpha = cd_sc.A/A_cd
    
    # ================================================================
    " ----------------------- TWO PHASE ZONE ----------------------- "
    # ================================================================
    " Connections : "
    " Refrigerant:"
    r_ex_tp_cd.P = r_su_sc_cd.P
    r_ex_tp_cd.set_PhasePQ_0()
    
    " Water:"
    w_su_tp_cd.T = w_ex_sc_cd.T
    w_su_tp_cd.P = w_ex_sc_cd.P
    w_su_tp_cd.set_PropsPT()
    w_ex_tp_cd.P = w_su_tp_cd.P
    
    " Supply Ref."
    r_su_tp_cd.P = r_ex_tp_cd.P + DELTAP_tp
    r_su_tp_cd.set_PhasePQ_1()
    cp_r_bar_tp_cd = (r_su_tp_cd.cp + r_ex_tp_cd.cp)/2
    
    " Enery balance"
    cd_tp.Q_dot = r_ex_tp_cd.M_dot*(r_su_tp_cd.H-r_ex_tp_cd.H) 
    w_ex_tp_cd.H = cd_tp.Q_dot/w_su_tp_cd.M_dot + w_su_tp_cd.H
    w_ex_tp_cd.set_PropsPH()
    cp_w_bar_tp_cd = (w_su_tp_cd.cp + w_ex_tp_cd.cp)/2
    
    " Capacitive flows "
    C_dot_r_tp_cd = np.nan            
    C_dot_w_tp_cd = cp_w_bar_tp_cd*w_su_tp_cd.M_dot
    C_min_tp_cd = C_dot_w_tp_cd  
                    
    " Effectivness- NTU "
    cd_tp.epsilon = min(1, cd_tp.Q_dot/(C_min_tp_cd*(r_su_tp_cd.T-w_su_tp_cd.T)))   
    cd_tp.NTU = -np.log(1-cd_tp.epsilon) 
    
    " -- Heat transfer coefficients -------------------------------- "
    hcv_r_tp_cd, DELTAP_tp_i, hcv_r_tp_cd_x = HTC.PHX_CD_Shah_2021(D_hyd = D_hyd,               
                  L_p = L_p,                  
                  B_p = B_p,                  
                  b = b,                   
                  phi = phi,               
                  M_dot_r = M_dot_r_cd,    
                  P_mean = (r_su_tp_cd.P + r_ex_tp_cd.P)/2,         
                  N_ch = N_ch,           
                  refrigerant = refrigerant)     

    hcv_w_tp_cd, DELTAP_tp_w_cd = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                L_p = L_p , 
                B_p = B_p,
                b = b,
                Beta = Beta, 
                M_dot_r = w_su_tp_cd.M_dot,  
                P_mean = (w_su_tp_cd.P + w_ex_tp_cd.P)/2,
                T_mean = (w_su_tp_cd.T + w_ex_tp_cd.T)/2,
                N_ch = N_ch,
                refrigerant = 'Water')
    
    " -- Global Heat transfer coefficient "
    R_tp_cd = 1/hcv_w_tp_cd + t/k_p + f_fouling + 1/hcv_r_tp_cd
    U_tp_cd = 1/R_tp_cd
    cd_tp.A = cd_tp.NTU*C_min_tp_cd/(U_tp_cd)
    cd_tp.alpha = cd_tp.A/A_cd
    
    # ================================================================
    " ---------------------- SUPERHEATING ZONE --------------------- "
    # ================================================================
    " Connections : "
    " Refrigerant."
    r_ex_sh_cd.P = r_su_tp_cd.P
    r_ex_sh_cd.set_PhasePQ_1()
    
    cp_r_bar_sh_cd = (r_ex_sh_cd.cp + r_su_sh_cd.cp)/2
    
    " Water"
    w_su_sh_cd.T = w_ex_tp_cd.T
    w_su_sh_cd.P = w_ex_tp_cd.P
    w_su_sh_cd.set_PropsPT()   
    w_ex_sh_cd.P = w_su_sh_cd.P
    
    " Enery balance"
    cd_sh.Q_dot = r_su_sh_cd.M_dot*(r_su_sh_cd.H-r_ex_sh_cd.H)  
    w_ex_sh_cd.H = cd_sh.Q_dot/w_su_sh_cd.M_dot + w_su_sh_cd.H
    w_ex_sh_cd.set_PropsPH()
    cp_w_bar_sh_cd = (w_ex_sh_cd.cp + w_su_sh_cd.cp)/2
    
    " Capacitive flows "
    C_dot_r_sh_cd = cp_r_bar_sh_cd*r_su_sh_cd.M_dot
    C_dot_w_sh_cd = cp_w_bar_sh_cd*w_su_sh_cd.M_dot
    C_min_sh_cd = min(C_dot_r_sh_cd, C_dot_w_sh_cd)                         
    C_max_sh_cd = max(C_dot_r_sh_cd, C_dot_w_sh_cd)                          
    C_r_sh_cd = C_min_sh_cd/C_max_sh_cd                              
    
    " Effectivness- NTU "
    cd_sh.epsilon = min(1, cd_sh.Q_dot/(C_min_sh_cd*(r_su_sh_cd.T - w_su_sh_cd.T)))
    cd_sh.NTU = 1/(C_r_sh_cd-1)*np.log((cd_sh.epsilon-1)/(cd_sh.epsilon*C_r_sh_cd-1))   
    
    " -- Heat transfer coefficients -------------------------------- "
    hcv_r_sh_cd, DELTAP_sh_i = HTC.PHX_1PH_Martin_VDI(D_hyd =D_hyd,
                    L_p = L_p , 
                    B_p = B_p,
                    b = b,
                    Beta = Beta, 
                    M_dot_r = M_dot_r_cd,
                    P_mean = (r_ex_sh_cd.P + r_su_sh_cd.P)/2,
                    T_mean = (r_su_sh_cd.T + r_ex_sh_cd.T)/2,
                    N_ch = N_ch,
                    refrigerant = refrigerant)
  
    hcv_w_sh_cd, DELTAP_sh_w_cd = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                L_p = L_p , 
                B_p = B_p,
                b = b,
                Beta = Beta, 
                M_dot_r = w_su_sh_cd.M_dot,  
                P_mean = (w_su_sh_cd.P + w_ex_sh_cd.P)/2,
                T_mean = (w_su_sh_cd.T + w_ex_sh_cd.T)/2,
                N_ch = N_ch,
                refrigerant = 'Water')
    
    " -- Global Heat transfer coefficient "
    R_sh_cd = 1/(hcv_w_sh_cd) + t/k_p + f_fouling + 1/(hcv_r_sh_cd)
    U_sh_cd = 1/R_sh_cd
    cd_sh.A = cd_sh.NTU*C_min_sh_cd/(U_sh_cd)
    cd_sh.alpha = cd_sh.A/A_cd

    # ================================================================
    " ------------------ MASS CHARGUE COMPUTATION ------------------ "
    # ================================================================  
    " -- Void Fraction calculation:"
    V_sc_cd = cd_sc.alpha*V_t
    V_tp_cd = cd_tp.alpha*V_t
    V_sh_cd = cd_sh.alpha*V_t
    
    discretization_tp = 100
    rho_tp_vd_i = np.zeros(discretization_tp)
    M_bar_tp_i = np.zeros(discretization_tp)
    
    Q_ite = np.linspace(0,1,discretization_tp)
    P_ite = np.linspace(r_su_tp_cd.P, r_ex_tp_cd.P, discretization_tp)
    V_ite_tp = V_tp_cd/discretization_tp
    
    for i in range(0, discretization_tp):
        
        rho_tp_vd_i[i], _ = VF.Biphasic(Q = Q_ite[i],
                                P = P_ite[i],
                                refrigerant = r_su_tp_cd.type)
        
        M_bar_tp_i[i] = V_ite_tp*rho_tp_vd_i[i]  
    
    
    M_bar_sc = V_sc_cd*(r_su_sc_cd.rho + r_ex_sc_cd.rho)/2
    M_bar_tp = np.sum(M_bar_tp_i)    
    M_bar_sh = V_sh_cd*(r_su_sh_cd.rho + r_ex_sh_cd.rho)/2      
    
    M_bar_t = M_bar_sc + M_bar_tp + M_bar_sh
    
    # ================================================================
    " --------------------- OVERALL CONDITIONS --------------------- "
    # ================================================================ 
    "-- Area values:"    
    A_cd_i = cd_sc.A + cd_tp.A + cd_sh.A
    R_A_cd = A_cd-A_cd_i

    "-- Pressure Drop values:"
    R_dp_sc = DELTAP_sc - DELTAP_sc_i
    R_dp_tp = DELTAP_tp - DELTAP_tp_i/10
    R_dp_sh = DELTAP_sh - DELTAP_sh_i/50
    # R_dp_sc = 0
    # R_dp_tp = 0
    # R_dp_sh = 0

    cd_sc.DELTAP_r = DELTAP_sc
    cd_tp.DELTAP_r = DELTAP_tp
    cd_sh.DELTAP_r = DELTAP_sh
   
    Q_dot_cd = cd_sc.Q_dot + cd_tp.Q_dot + cd_sh.Q_dot    
    DELTAP_t = cd_sc.DELTAP_r + cd_tp.DELTAP_r + cd_sh.DELTAP_r 

    return R_A_cd, R_dp_sc, R_dp_tp, R_dp_sh, Q_dot_cd, DELTAP_t, r_ex_sc_cd.T,\
        (r_ex_tp_cd.P+ r_su_tp_cd.P)/2,r_ex_sc_cd.P, w_ex_sh_cd.T, w_ex_sh_cd.P,\
            w_ex_sh_cd.M_dot, cd_sc.alpha, cd_tp.alpha, cd_sh.alpha, M_bar_t, \
                M_bar_sc, M_bar_tp, M_bar_sh


if __name__ == "__main__":

    def Opt(x):
        # print('\n')
        print(f'iteration: {x}')
        a = ThreePhasePHX(x[0], x[1], x[2], x[3],   
                        M_dot_r_cd = 2.345,             # Refrigerant mass flow rate, kg s^-1
                        T_r_su_cd = 49.96 + 273.15,     # Refrigerant supply Temperature, K
                        P_w_su_cd = 1e+05,             # Water pressure supply, Pa
                        T_w_su_cd = 15+273.15,          # Water Temperature supply, K
                        M_dot_w_cd = False,            # Water glide, K
                        T_glide = 7,                   # Water glide, K
                        DELTAT_sc = 7.18               # Subcooling  
                        )
        
        print(f'Res:{a[0], a[1], a[2], a[3]}')
        print('\n')
        # print('\n')
        return a[0], a[1], a[2], a[3]
    
    x_0 = [128435.76958363,    578.24976693,  10761.63573897,    977.51164651]
    
    import time
    start = time.time()
    sol = root(Opt, x_0, method = 'hybr')
    stop = time.time()
    a = ThreePhasePHX(sol.x[0],
                      sol.x[1],
                      sol.x[2],
                      sol.x[3],
                        M_dot_r_cd = 2.345,             # Refrigerant mass flow rate, kg s^-1
                        T_r_su_cd = 49.96 + 273.15,     # Refrigerant supply Temperature, K
                        P_w_su_cd = 1e+05,             # Water pressure supply, Pa
                        T_w_su_cd = 15+273.15,          # Water Temperature supply, K
                        M_dot_w_cd = False,            # Water glide, K
                        T_glide = 7,                   # Water glide, K
                        DELTAT_sc = 7.18               # Subcooling  
                        )
    
    
    print('\n')
    print('|-----------------------------------------|')
    print('|-- Optimization Function ----------------|')
    print('|-----------------------------------------|')
    print(f' Optimization Time: {stop-start} [s]')
    print(f' Residual: {a[0]} [-]')
    print(f' Q_dot_cd: {a[4]/1000} [kW]')
    print(f' DeltaP_r: {a[5]/1e+05} [bar]')
    print(f' T_r_ex_cd: {a[6]-273.15} [°C]')
    print(f' P_r_cd: {a[7]/1e+05} [bar]')
    print(f' P_r_ex_cd: {a[8]/1e+05} [bar]')
    print(f' T_r_sat: {(r_ex_tp_cd.T-273.15 + r_su_tp_cd.T-273.15)/2} [°C]')
    print(f' T_w_ex_cd: {a[9]-273.15} [°C]')
    print(f' P_w_ex_cd: {a[10]/1e+05} [bar]')
    print(f' M_dot_w_cd: {a[11]} [kg/s]')
    print(f' alpha_sc: {a[12]} [kg/s]')
    print(f' alpha_tp: {a[13]} [kg/s]')
    print(f' alpha_sh: {a[14]} [kg/s]')
    print(f' Refrigerant Charge: {a[15]} [kg]')        
    print(f' Area condenser: {A_cd}')
    print(f' Pinch Point: {r_ex_sh_cd.T - w_ex_tp_cd.T}')
    print('|-----------------------------------------|')    
    print('|-----------------------------------------|')  
    
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8, 8))  
    ax = fig.add_subplot(111)
    ax.grid('on', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'$\alpha_{cd}$ [-]', fontsize=20, labelpad=15)
    ax.set_ylabel(r'$T_{}$ [°C]', fontsize=20, labelpad=15)
    plt.plot([ 1-a[12], 1], [r_su_sc_cd.T-273.15, r_ex_sc_cd.T-273.15], label='Subcooling', linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([1-(a[12] + a[13]),1- a[12]], [r_su_tp_cd.T-273.15, r_ex_tp_cd.T-273.15], label='Two-phase', linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([1-(a[12] + a[13] + a[14]), 1- (a[12] + a[13])], [r_su_sh_cd.T-273.15, r_ex_sh_cd.T-273.15], label='Superheating', linewidth=4)  # Cambia el ancho de la línea a 4
    # plt.plot([1, 0], [w_su_sc_cd.T-273.15, w_ex_sh_cd.T-273.15], label='Water', linewidth=4)  # Cambia el ancho de la línea a 4
    
    plt.plot([ 1-a[12], 1], [w_ex_sc_cd.T-273.15, w_su_sc_cd.T-273.15], linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([ 1-(a[12] + a[13]),1- a[12]], [w_ex_tp_cd.T-273.15, w_su_tp_cd.T-273.15], linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([1-(a[12] + a[13] + a[14]), 1- (a[12] + a[13])], [w_ex_sh_cd.T-273.15, w_su_sh_cd.T-273.15], linewidth=4)  # Cambia el ancho de la línea a 4
    

    ax.text(0, 45, f'Presión: {round(a[7]/1e+05,4)} [bar]', fontsize=15, color='green')
    ax.text(0, 42, f'Area: {round(A_cd,4)} [m^2]', fontsize=15, color='green')
    ax.text(0, 39, f'Pinch Point: {round(r_ex_sh_cd.T - w_ex_tp_cd.T,4)} [K]', fontsize=15, color='green')
    ax.text(0, 36, f'T_r_ex_cd : {round(a[6]-273.15,4)} [°C]', fontsize=15, color='green')
    ax.text(0, 33, f'T_w_ex_cd : {round(a[9]-273.15,4)} [°C]', fontsize=15, color='green')
    ax.text(0, 30, f'T_r_sat: {round(r_ex_tp_cd.T-273.15,4)} [°C]', fontsize=15, color='green')
    ax.text(0, 27, f'Q_dot_cd: {round(a[4]/1000,4)} [°C]', fontsize=15, color='green')    
    ax.text(0, 24, f'NTU_sc: {round(cd_sc.NTU,4)}, NTU_tp: {round(cd_tp.NTU,4)}, NTU_sh: {round(cd_sh.NTU,4)} [-]', fontsize=15, color='green') 
    ax.text(0, 20, f'N°Plates: {round(N_p,4)}', fontsize=15, color='red') 
    
    ax.legend( loc='upper right')
    plt.tight_layout()
    plt.show()
    plt.savefig('CD-ORC-B633x250.svg')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
