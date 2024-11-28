# =============================================================================
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
L_p = 854e-3          # Port-port centerline distance (Large), m
B_p = 179e-3          # Port-port centerline distance (Width), m
D_p = 100e-3          # Port Diameter, m
N_pas = 1             # Number of passes per channels
N_p = 220              # N° of Plates, - (only uses par number if impar -1)
t = 0.25e-3           # Plate thickness, m
L_c = (12 + 2.29*N_p)*1e-3        # Depth of the heat exchanger, m 
Lambda = 2.65e-3       # Plate corrugation wavelength, m, in kaka{c} the name corresponds to Pc
k_p = 385             # Conductivity of the plate (copper), W m^-1K^-1
V_t = 62.247/1000      # Total volume
f_fouling = 0.363/1000  # Fouling factor
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
A_ev = A_tot                              # 6.35

# %%
# =============================================================================
" -------------------------- Evaporator Model ------------------------------- "
# =============================================================================
" Declaration of the point: "
r_su_sh_ev = Point.Refrigerant()
r_su_tp_ev = Point.Refrigerant()
r_su_sc_ev = Point.Refrigerant()
r_ex_sh_ev = Point.Refrigerant()
r_ex_tp_ev = Point.Refrigerant()
r_ex_sc_ev = Point.Refrigerant()
r_sat_ev = Point.Refrigerant()

" Set the refrigerant in the differents poins of the cycle"
r_sat_ev.type = r_su_sh_ev.type = r_ex_sh_ev.type \
    = r_su_tp_ev.type = r_ex_tp_ev.type = r_su_sc_ev.type\
        = r_ex_sc_ev.type = refrigerant

w_su_sh_ev = Point.Water()
w_su_tp_ev = Point.Water()
w_su_sc_ev = Point.Water()
w_ex_sh_ev = Point.Water()
w_ex_tp_ev = Point.Water()
w_ex_sc_ev = Point.Water()

ev_sc = Point.Condenser()
ev_tp = Point.Condenser()
ev_sh = Point.Condenser()

def ThreePhasePHX(P_r_ev_i,
                  DELTAP_sc,
                  DELTAP_tp,
                  DELTAP_sh,
                M_dot_r_ev = 2.38,             # Refrigerant mass flow rate, kg s^-1
                T_r_su_ev = 11.39 + 273.15,    # Refrigerant supply Temperature, K
                P_w_su_ev = 4e+05,             # Water pressure supply, Pa
                T_w_su_ev = 120+273.15,        # Water Temperature supply, K
                M_dot_w_ev = 7.092,            # Water Mass flow rate, kg s^-1
                T_glide = False,               # Water glide, K
                DELTAT_sh = 7.16):             # Superheating  
    
    # ================================================================
    " ------------------- SF-INITIAL CONDITION --------------------- "
    # ================================================================
    " Supply refrigerant temperature "
    r_su_sc_ev.T = T_r_su_ev
    r_su_sc_ev.P = P_r_ev_i + DELTAP_sc + DELTAP_tp
    r_su_sc_ev.M_dot = M_dot_r_ev
    r_ex_sc_ev.M_dot = r_ex_tp_ev.M_dot = r_su_tp_ev.M_dot = r_ex_sh_ev.M_dot\
                          = r_su_sh_ev.M_dot = r_su_sc_ev.M_dot
    r_su_sc_ev.set_PropsPT()

    " Exhaust refrigerant conditions:"
    r_ex_tp_ev.P = P_r_ev_i
    r_ex_tp_ev.set_PhasePQ_1()

    r_ex_sh_ev.T = r_ex_tp_ev.T_sat + DELTAT_sh
    r_ex_sh_ev.P = P_r_ev_i - DELTAP_sh
    r_ex_sh_ev.set_PropsPT()
    
    " Energy Balance"
    Q_dot_r_ev = r_su_sh_ev.M_dot*(r_ex_sh_ev.H-r_su_sc_ev.H)

    " Water connection: "
    w_su_sh_ev.T = T_w_su_ev
    w_su_sh_ev.P = P_w_su_ev
    w_su_sh_ev.set_PropsPT()
    
    " Water exhaust: "
    if T_glide == False:
        w_ex_sc_ev.P = w_su_sh_ev.P
        w_su_sh_ev.M_dot = w_ex_sc_ev.M_dot = M_dot_w_ev
        w_ex_sc_ev.H = w_su_sh_ev.H - Q_dot_r_ev/w_su_sh_ev.M_dot
        w_ex_sc_ev.set_PropsPH()
        w_ex_sc_ev.M_dot = w_su_tp_ev.M_dot = w_ex_tp_ev.M_dot = w_ex_sh_ev.M_dot = w_su_sc_ev.M_dot= w_su_sh_ev.M_dot
    else:
        w_ex_sc_ev.T = w_su_sh_ev.T - T_glide
        w_ex_sc_ev.P = w_su_sh_ev.P
        w_ex_sc_ev.set_PropsPT()
        w_su_sh_ev.M_dot = Q_dot_r_ev/(w_su_sh_ev.H-w_ex_sc_ev.H)
        w_ex_sc_ev.M_dot = w_su_tp_ev.M_dot = w_ex_tp_ev.M_dot = w_ex_sh_ev.M_dot = w_su_sc_ev.M_dot= w_su_sh_ev.M_dot
   
    # ================================================================
    " ----------------------- SUBCOOLED ZONE ----------------------- "
    # ================================================================
    " Refrigerant properties "
    r_ex_sc_ev.P = P_r_ev_i + DELTAP_tp
    r_ex_sc_ev.set_PhasePQ_0()
    cp_r_bar_sc_ev = (r_su_sc_ev.cp + r_ex_sc_ev.cp)/2
    " Water properties "
    w_su_sc_ev.P = w_su_sh_ev.P
    
    " Energy balance "
    ev_sc.Q_dot = r_su_sc_ev.M_dot*(r_ex_sc_ev.H - r_su_sc_ev.H) 
    w_su_sc_ev.H = ev_sc.Q_dot/w_su_sc_ev.M_dot + w_ex_sc_ev.H
    w_su_sc_ev.set_PropsPH()
    cp_w_bar_sc_ev = (w_su_sc_ev.cp + w_ex_sc_ev.cp)/2

    " Capacitive flows "
    C_dot_r_sc_ev = cp_r_bar_sc_ev*r_su_sc_ev.M_dot
    C_dot_w_sc_ev = cp_w_bar_sc_ev*w_su_sc_ev.M_dot
    C_min_sc_ev = min(C_dot_r_sc_ev, C_dot_w_sc_ev)                           
    C_max_sc_ev = max(C_dot_r_sc_ev, C_dot_w_sc_ev)                           
    C_r_sc_ev = C_min_sc_ev/C_max_sc_ev                                       

    " Effectivness- NTU "
    ev_sc.epsilon = min(1, ev_sc.Q_dot/(C_min_sc_ev*(w_su_sc_ev.T-r_su_sc_ev.T)))
    ev_sc.NTU = 1/(C_r_sc_ev-1)*np.log((ev_sc.epsilon-1)/(ev_sc.epsilon*C_r_sc_ev-1))   

    " -- Heat transfer coefficients -------------------------------- "
    hcv_r_sc_ev, DELTAP_sc_i = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                    L_p = L_p , 
                    B_p = B_p,
                    b = b,
                    Beta = Beta, 
                    M_dot_r = M_dot_r_ev,
                    P_mean = (r_su_sc_ev.P + r_ex_sc_ev.P)/2,
                    T_mean = (r_ex_sc_ev.T + r_su_sc_ev.T)/2,
                    N_ch = N_ch,
                    refrigerant = refrigerant)

    hcv_w_sc_ev, DELTAP_sc_w_ev = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                L_p = L_p , 
                B_p = B_p,
                b = b,
                Beta = Beta, 
                M_dot_r = w_su_sc_ev.M_dot,  
                P_mean = (w_su_sc_ev.P + w_ex_sc_ev.P)/2,
                T_mean = (w_su_sc_ev.T + w_ex_sc_ev.T)/2,
                N_ch = N_ch,
                refrigerant = 'Water')

    " -- Global Heat transfer coefficient  "
    R_sc_ev = 1/(hcv_w_sc_ev) + t/k_p + f_fouling + 1/(hcv_r_sc_ev)
    ev_sc.U = 1/R_sc_ev
    ev_sc.A = ev_sc.NTU*C_min_sc_ev/(ev_sc.U)
    ev_sc.alpha = ev_sc.A/A_ev
    
    # ================================================================
    " ----------------------- TWO PHASE ZONE ----------------------- "
    # ================================================================
    " Connections : "
    " Refrigerant:"
    r_su_tp_ev.P = r_ex_sc_ev.P
    r_su_tp_ev.set_PhasePQ_0()
    
    " Water:"
    w_ex_tp_ev.T = w_su_sc_ev.T
    w_ex_tp_ev.P = w_su_sc_ev.P
    w_ex_tp_ev.set_PropsPT()
    w_su_tp_ev.P = w_su_sc_ev.P
    
    " Supply Refrifrerant"
    cp_r_bar_tp_ev = (r_su_tp_ev.cp + r_ex_tp_ev.cp)/2
    
    " Enery balance"
    ev_tp.Q_dot = r_ex_tp_ev.M_dot*(r_ex_tp_ev.H-r_su_tp_ev.H) 
    w_su_tp_ev.H = ev_tp.Q_dot/w_su_tp_ev.M_dot + w_ex_tp_ev.H
    w_su_tp_ev.set_PropsPH()
    cp_w_bar_tp_ev = (w_su_tp_ev.cp + w_ex_tp_ev.cp)/2

    " Capacitive flows "
    C_dot_r_tp_ev = np.nan            
    C_dot_w_tp_ev = cp_w_bar_tp_ev*w_su_tp_ev.M_dot
    C_min_tp_ev = C_dot_w_tp_ev  
                    
    " Effectivness- NTU "
    ev_tp.epsilon = min(1, ev_tp.Q_dot/(C_min_tp_ev*(w_su_tp_ev.T-r_su_tp_ev.T)))    
    ev_tp.NTU = -np.log(1-ev_tp.epsilon) 

    " -- Heat transfer coefficients -------------------------------- "
    hcv_r_tp_ev, DELTAP_tp_i, hcv_r_tp_ev_i = HTC.PHX_EV_Amalfi_2015(D_hyd = D_hyd,               
                  L_p = L_p,                  
                  B_p = B_p,                  
                  b = b, 
                  Beta = Beta,
                  A_tot = A_tot,
                  M_dot_r = M_dot_r_ev,    
                  P_mean = (r_su_tp_ev.P + r_ex_tp_ev.P)/2,         
                  N_ch = N_ch,           
                  refrigerant = refrigerant)       
    
    hcv_w_tp_ev, DELTAP_tp_w_ev = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                L_p = L_p , 
                B_p = B_p,
                b = b,
                Beta = Beta, 
                M_dot_r = w_su_tp_ev.M_dot,  
                P_mean = (w_su_tp_ev.P + w_ex_tp_ev.P)/2,
                T_mean = (w_su_tp_ev.T + w_ex_tp_ev.T)/2,
                N_ch = N_ch,
                refrigerant = 'Water')                    

    " -- Global Heat transfer coefficient "
    R_tp_ev = 1/hcv_w_tp_ev + t/k_p + f_fouling + 1/(hcv_r_tp_ev)
    ev_tp.U = 1/R_tp_ev
    ev_tp.A = ev_tp.NTU*C_min_tp_ev/(ev_tp.U)
    ev_tp.alpha = ev_tp.A/A_ev

    # ================================================================
    " ---------------------- SUPERHEATING ZONE --------------------- "
    # ================================================================
    " Connections : "
    " Refrigerant."
    r_su_sh_ev.P = r_ex_tp_ev.P 
    r_su_sh_ev.set_PhasePQ_1()
    cp_r_bar_sh_ev = (r_ex_sh_ev.cp + r_su_sh_ev.cp)/2
    
    " Water"
    w_ex_sh_ev.T = w_su_tp_ev.T
    w_ex_sh_ev.P = w_su_tp_ev.P
    w_ex_sh_ev.set_PropsPT()   
    w_su_sh_ev.P = w_ex_sh_ev.P
    
    " Enery balance"
    ev_sh.Q_dot = r_su_sh_ev.M_dot*(r_ex_sh_ev.H-r_su_sh_ev.H)  
    w_su_sh_ev.H = ev_sh.Q_dot/w_su_sh_ev.M_dot + w_ex_sh_ev.H
    w_su_sh_ev.set_PropsPH()
    cp_w_bar_sh_ev = (w_ex_sh_ev.cp + w_su_sh_ev.cp)/2

    " Capacitive flows "
    C_dot_r_sh_ev = cp_r_bar_sh_ev*r_su_sh_ev.M_dot
    C_dot_w_sh_ev = cp_w_bar_sh_ev*w_su_sh_ev.M_dot
    C_min_sh_ev = min(C_dot_r_sh_ev, C_dot_w_sh_ev)                      
    C_max_sh_ev = max(C_dot_r_sh_ev, C_dot_w_sh_ev)                         
    C_r_sh_ev = C_min_sh_ev/C_max_sh_ev                                       
    
    " Effectivness- NTU "
    ev_sh.epsilon = ev_sh.Q_dot/(C_min_sh_ev*(w_su_sh_ev.T - r_su_sh_ev.T))
    ev_sh.NTU = 1/(C_r_sh_ev-1)*np.log((ev_sh.epsilon-1)/(ev_sh.epsilon*C_r_sh_ev-1))   

    " -- Heat transfer coefficients -------------------------------- "
    hcv_r_sh_ev, DELTAP_sh_i = HTC.PHX_1PH_Martin_VDI(D_hyd =D_hyd,
                    L_p = L_p , 
                    B_p = B_p,
                    b = b,
                    Beta = Beta, 
                    M_dot_r = M_dot_r_ev,
                    P_mean = (r_ex_sh_ev.P + r_su_sh_ev.P)/2,
                    T_mean = (r_su_sh_ev.T + r_ex_sh_ev.T)/2,
                    N_ch = N_ch,
                    refrigerant = refrigerant)
    
    hcv_w_sh_ev, DELTAP_sh_w_ev = HTC.PHX_1PH_Martin_VDI(D_hyd = D_hyd,
                L_p = L_p , 
                B_p = B_p,
                b = b,
                Beta = Beta, 
                M_dot_r = w_su_sh_ev.M_dot,  
                P_mean = (w_su_sh_ev.P + w_ex_sh_ev.P)/2,
                T_mean = (w_su_sh_ev.T + w_ex_sh_ev.T)/2,
                N_ch = N_ch,
                refrigerant = 'Water')        
         
    " -- Global Heat transfer coefficient "
    R_sh_ev = 1/(hcv_w_sh_ev) + t/k_p + f_fouling + 1/(hcv_r_sh_ev)
    ev_sh.U = 1/R_sh_ev
    ev_sh.A = ev_sh.NTU*C_min_sh_ev/(ev_sh.U)
    ev_sh.alpha = ev_sh.A/A_ev
    
    # ================================================================
    " ------------------ MASS CHARGUE COMPUTATION ------------------ "
    # ================================================================    
    " -- Void Fraction calculation:"
    V_sc_ev = ev_sc.alpha*V_t
    V_tp_ev = ev_tp.alpha*V_t
    V_sh_ev = ev_sh.alpha*V_t
    
    discretization_tp = 100
    rho_tp_vd_i = np.zeros(discretization_tp)
    M_bar_tp_i = np.zeros(discretization_tp)
    
    Q_ite = np.linspace(0,1,discretization_tp)
    P_ite = np.linspace(r_su_tp_ev.P, r_ex_tp_ev.P, discretization_tp)
    V_ite_tp = V_tp_ev/discretization_tp
    
    for i in range(0, discretization_tp):
        
        rho_tp_vd_i[i], _ = VF.Biphasic(Q = Q_ite[i],
                                P = P_ite[i],
                                refrigerant = r_su_tp_ev.type)
        
        M_bar_tp_i[i] = V_ite_tp*rho_tp_vd_i[i]  
        

    M_bar_sc = V_sc_ev*(r_su_sc_ev.rho + r_ex_sc_ev.rho)/2
    M_bar_tp = np.sum(M_bar_tp_i)   
    M_bar_sh = V_sh_ev*(r_su_sh_ev.rho + r_ex_sh_ev.rho)/2  
    
    M_bar_t = M_bar_sc + M_bar_tp + M_bar_sh

    # ================================================================
    " --------------------- OVERALL CONDITIONS --------------------- "
    # ================================================================      
    "-- Area values:"
    A_ev_i = ev_sc.A + ev_tp.A + ev_sh.A
    R_A_ev = A_ev-A_ev_i
    
    "-- Pressure Drop values:"
    R_dp_sc = DELTAP_sc - DELTAP_sc_i
    R_dp_tp = DELTAP_tp - DELTAP_tp_i/10
    R_dp_sh = DELTAP_sh - DELTAP_sh_i/10
    
    ev_sc.DELTAP_r = DELTAP_sc
    ev_tp.DELTAP_r = DELTAP_tp
    ev_sh.DELTAP_r = DELTAP_sh
   
    Q_dot_ev = ev_sc.Q_dot + ev_tp.Q_dot + ev_sh.Q_dot    
    DELTAP_t = ev_sc.DELTAP_r + ev_tp.DELTAP_r + ev_sh.DELTAP_r 

    return R_A_ev, R_dp_sc, R_dp_tp, R_dp_sh , Q_dot_ev, DELTAP_t, r_ex_sh_ev.T,\
        (r_ex_tp_ev.P+ r_su_tp_ev.P)/2,r_ex_sh_ev.P, w_ex_sc_ev.T, w_ex_sc_ev.P, \
            w_ex_sc_ev.M_dot, ev_sc.alpha, ev_tp.alpha, ev_sh.alpha, M_bar_t,\
                M_bar_sc, M_bar_tp, M_bar_sh


if __name__ == "__main__":
    
    def Opt(x):
        print(f'iteration: {x}')
        a = ThreePhasePHX(x[0], x[1], x[2], x[3],
                          M_dot_r_ev = 2.355,         # Refrigerant mass flow rate, kg s^-1
                          T_r_su_ev = 17.9 + 273.15,    # Refrigerant supply Temperature, K
                          P_w_su_ev = 4e+05,      # Water pressure supply, Pa
                          T_w_su_ev = 90+273.15,    # Water Temperature supply, K
                          M_dot_w_ev = 13.63,        # Water Mass flow rate, kg s^-1
                          T_glide = False,          # Water glide, K
                          DELTAT_sh = 6.54)             # Superheating  

       
        print(f'Res:{a[0], a[1], a[2], a[3]}')
        print(f'\n')
        return a[0], a[1], a[2], a[3]

    x_0 = [699741.26324852,    739.94429034,   2234.35695919,   1360.79406851]
           
    import time
    start = time.time()
    sol = root(Opt, x_0, method = 'hybr')
    stop = time.time()
    print('\n')
    print(sol.x)
    a = ThreePhasePHX(sol.x[0],
                      sol.x[1],
                      sol.x[2],
                      sol.x[3],
                          M_dot_r_ev = 2.355,         # Refrigerant mass flow rate, kg s^-1
                          T_r_su_ev = 17.9 + 273.15,    # Refrigerant supply Temperature, K
                          P_w_su_ev = 4e+05,      # Water pressure supply, Pa
                          T_w_su_ev = 90+273.15,    # Water Temperature supply, K
                          M_dot_w_ev = 13.63,        # Water Mass flow rate, kg s^-1
                          T_glide = False,          # Water glide, K
                          DELTAT_sh = 6.54)             # Superheating  
  
    print('\n')
    print('|-----------------------------------------|')
    print('|-- Optimization Function ----------------|')
    print('|-----------------------------------------|')
    print(f' Optimization Time: {stop-start} [s]')
    print(f' Residual: {a[0]} [-]')
    print(f' Q_dot_ev: {a[4]/1000} [kW]')
    print(f' DeltaP_r: {a[5]/1e+05} [bar]')
    print(f' T_r_ex_ev: {a[6]-273.15} [°C]')
    print(f' P_r_ev: {a[7]/1e+05} [bar]')
    print(f' T_r_sat: {(r_ex_tp_ev.T-273.15 + r_su_tp_ev.T-273.15)/2} [°C]')

    print(f' P_r_ex_ev: {a[8]/1e+05} [bar]')
    print(f' T_w_ex_ev: {a[9]-273.15} [°C]')
    print(f' P_w_ex_ev: {a[10]/1e+05} [bar]')
    print(f' M_dot_w_ev: {a[11]} [kg/s]')
    print(f' alpha_sc: {a[12]} [kg/s]')
    print(f' alpha_tp: {a[13]} [kg/s]')
    print(f' alpha_sh: {a[14]} [kg/s]')    
    print(f' Refrigerant Charge: {a[15]} [kg]')    
    print(f' Area condenser: {A_ev}')
    print(f' Pinch Point: {w_ex_tp_ev.T-r_ex_sc_ev.T}')
    print('|-----------------------------------------|')    
    print('|-----------------------------------------|')  
    
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8, 8))  
    ax = fig.add_subplot(111)
    ax.grid('on', alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'$\alpha_{ev}$ [-]', fontsize=20, labelpad=15)
    ax.set_ylabel(r'$T_{}$ [°C]', fontsize=20, labelpad=15)
    plt.plot([0, a[12]], [r_su_sc_ev.T-273.15, r_ex_sc_ev.T-273.15], label='Subcooling', linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([a[12], a[12] + a[13]], [r_su_tp_ev.T-273.15, r_ex_tp_ev.T-273.15], label='Two-phase', linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([a[12] + a[13], a[12] + a[13] + a[14]], [r_su_sh_ev.T-273.15, r_ex_sh_ev.T-273.15], label='Superheating', linewidth=4)  # Cambia el ancho de la línea a 4
    # plt.plot([0, 1], [w_ex_sc_ev.T-273.15, w_su_sh_ev.T-273.15], label='Water', linewidth=4)  # Cambia el ancho de la línea a 4

    plt.plot([0, a[12]], [w_ex_sc_ev.T-273.15, w_su_sc_ev.T-273.15], linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([a[12], a[12] + a[13]], [w_ex_tp_ev.T-273.15, w_su_tp_ev.T-273.15],  linewidth=4)  # Cambia el ancho de la línea a 4
    plt.plot([a[12] + a[13], a[12] + a[13] + a[14]], [w_ex_sh_ev.T-273.15, w_su_sh_ev.T-273.15], linewidth=4)  # Cambia el ancho de la línea a 4
    
    ax.text(0, 75, f'Presión: {round(a[7]/1e+05,4)} [bar]', fontsize=15, color='green')
    ax.text(0, 70, f'Area: {round(A_ev,4)} [m^2]', fontsize=15, color='green')
    ax.text(0, 65, f'Pinch Point: {round(w_ex_tp_ev.T-r_ex_sc_ev.T,4)} [K]', fontsize=15, color='green')   
    ax.text(0, 60, f'T_r_ex_ev : {round(a[6]-273.15,4)} [°C]', fontsize=15, color='green')
    ax.text(0, 55, f'T_w_ex_ev : {round(a[9]-273.15,4)} [°C]', fontsize=15, color='green')
    ax.text(0, 50, f'T_r_sat: {(r_ex_tp_ev.T-273.15 + r_su_tp_ev.T-273.15)/2} [°C]', fontsize=15, color='green')


    ax.text(0, 45, f'Q_dot_ev: {round(a[4]/1000,4)} [°C]', fontsize=15, color='green')    
    ax.text(0, 40, f'NTU_sc: {round(ev_sc.NTU,4)}, NTU_tp: {round(ev_tp.NTU,4)}, NTU_sh: {round(ev_sh.NTU,4)} [-]', fontsize=15, color='green') 
    ax.text(0, 35, f'N°Plates: {round(N_p,4)}', fontsize=15, color='red')   
    ax.legend( loc='lower right')
    
    plt.tight_layout()
    plt.show()
    plt.savefig('EV-B439X110_2.svg')
    
    
