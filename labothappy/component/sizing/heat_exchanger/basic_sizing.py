import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def f_lmtd2(R,P,params):
    """
    Input
    -----
    
    R = (T_h_i - T_h_ex)/(T_c_ex - T_c_i) --> must be positive [0; infinte [
    P = (T_c_ex - T_c_i)/(T_h_i - T_c_i) --> must be within 0 and 1
    
    params : heat exchanger paramater dictionnary
    C_r = C_min/C_max
    
    Output
    ------
    
    F : Correction factor (0 to 1) for the LTMD method for non counter flow configurations
    
    Reference
    ---------
    EES code for LMTD method
    """
    
    def res_NTU_shell_and_tube(NTU, epsilon, C_r, params):
                
        NTU = max(1e-4, NTU)
        fact_1 = 1 + C_r + np.sqrt(1+C_r**2)
        fact_2 = 1+np.exp(-NTU*np.sqrt(1+C_r**2))
        fact_3 = 1-np.exp(-NTU*np.sqrt(1+C_r**2))
        
        if params['n_series'] == 1:
            epsilon_g = 2/(fact_1*(fact_2)/(fact_3))
        else:
            eps_1 = 2/(fact_1*(fact_2)/(fact_3))
            
            fact_4 = ((1-eps_1*C_r)/(1-eps_1))
            epsilon_g = (fact_4**(params['n_series']) - 1)/(fact_4**(params['n_series']) - C_r)
            
        res = (epsilon-epsilon_g)

        return res
    
    def res_NTU_crossflow(Ntu, epsilon, C_r):
        Ntu = max(1e-4, Ntu)
        epsilon_g = 1-np.exp((1/C_r)*Ntu**0.22*(np.exp(-C_r*Ntu**0.78)-1)) # Both fluids unmixed
        res = (epsilon-epsilon_g)
        return res
    
    def res_NTU_parallelflow(Ntu, epsilon, C_r):
        Ntu = max(1e-4, Ntu)
        epsilon_g = (1-np.exp(-Ntu*(1+C_r)))/(1+C_r)
        res = (epsilon-epsilon_g)
        return res    

    #R = (T_h_i - T_h_ex)/(T_c_ex - T_c_i) --> must be positive [0; infinte [
    #P = (T_c_ex - T_c_i)/(T_h_i - T_c_i) --> must be within 0 and 1
    
    #-----------------------------------------------------------------------------
    #R and P must be positive
    R = abs(R)
    P = abs(P)
    #-----------------------------------------------------------------------------
    R_min = 1e-4
    R_max = 100
    
    R = max(R_min,min(R_max,R))
    if R < 0.41:
        P_max_for_given_R = 0.99999;
    else:
        P_max_for_given_R = 0.995*((0.966475661367996)*R**3 + (-1.431274296912407)*R**2 + (0.247230065033875)*R + (0.309118607270897)) / (R**4 + (-1.766309686745371)*R**3 + (1.287179055148762)*R**2 + (-0.902512766020191)*R + (0.484880205333508))
    P = max(0,  min(P_max_for_given_R,P))
    
    if R <= 1: # Cdot_min = Cdot_h
        epsilon = P
        Pbis = P
        Rbis = R
        C_r = R
    else: # Cdot_max = Cdot_c
        epsilon = P*R
        Pbis = P*R
        Rbis = 1/R
        C_r = 1/R
        
    if params['Flow_Type'] == 'Shell&Tube':
        f = lambda x: res_NTU_shell_and_tube(x, epsilon, C_r,params)
    
    elif params['Flow_Type'] == 'CrossFlow':
        f = lambda x: res_NTU_crossflow(x, epsilon, C_r)        
    
    elif params['Flow_Type'] == 'ParallelFlow':
        f = lambda x: res_NTU_parallelflow(x, epsilon, C_r)   
    
    else:
        print("No e-NTU Correlation implemented for other Flow_Type than : 'Shell&Tube', 'CrossFlow', 'ParallelFlow' ")        
    
    out_fsolve = fsolve(f, 1, full_output= 1)#, args=(), fprime=None, full_output=0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)    
    NTU, res_NTU,flag_ntu = float(out_fsolve[0]), float(out_fsolve[1].get('fvec')), out_fsolve[2]
    
    if abs(res_NTU)< 1e-3  and flag_ntu > 0:
        if  R==1:
            F=P/NTU/(1-P)
        elif P <= 1e-04:
            F = 1
        else:
            F = epsilon*np.log((Pbis-1)/(Rbis*Pbis -1))/NTU/Pbis/(Rbis-1);
        flag = flag_ntu
    else:

        F = 1
        flag = flag_ntu
        
    return F

def find_UA(Q_dot, T_c_i, T_c_o, T_h_i, T_h_o, flow = 'CounterFlow', params = None):
    
    if flow == "ParallelFlow":
        if T_c_o == T_c_i or T_h_o == T_h_i:
            DTA = T_h_i - T_c_o
            DTB = T_h_o - T_c_i
        else:
            DTA = T_h_i - T_c_i
            DTB = T_h_o - T_c_o

        LMTD = (DTA - DTB)/np.log(DTA/DTB)
        UA = Q_dot/LMTD
        return UA    

    DTA = T_h_i - T_c_o
    DTB = T_h_o - T_c_i

    LMTD = (DTA - DTB)/np.log(DTA/DTB)
    
    if T_c_o == T_c_i or T_h_o == T_h_i or flow == "CounterFlow":
        F = 1

    else: 
        R = (T_h_i - T_h_o)/(T_c_o - T_c_i) 
        P = (T_c_o - T_c_i)/(T_h_i - T_c_i) 

        if flow == "CrossFlow": # unmixed
            F = f_lmtd2(R,P,params)
        elif flow == "Shell&Tube": 
            F = f_lmtd2(R,P,params)
        elif flow == "ParallelFlow":
            F = f_lmtd2(R,P,params)
        else:
            print("That flow configuration is not implemented.")

    F = min(1,F)

    UA = Q_dot/(LMTD*F)

    return UA

