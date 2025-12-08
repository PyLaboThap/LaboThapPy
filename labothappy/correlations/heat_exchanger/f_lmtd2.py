# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:50:11 2021

@author: jvega
"""

import numpy as np
from scipy.optimize import brentq

#%%

def F_shell_and_tube(R1, P1, N):
    # Limit P1 based on R1
    if R1 < 0.41:
        P_max = 0.99999
    else:
        num = 0.995 * (0.966475661367996*R1**3 - 1.431274296912407*R1**2 +
                       0.247230065033875*R1 + 0.309118607270897)
        den = (R1**4 - 1.766309686745371*R1**3 + 1.287179055148762*R1**2 -
               0.902512766020191*R1 + 0.484880205333508)
        P_max = num / den

    # Ensure P1 stays within [0, P_max]
    P1 = min(max(P1, 0), P_max)

    if R1 == 1:
        W_prime = (N - N * P1) / (N - N * P1 + P1)

        F = 2**0.5 * ((1 - W_prime) / W_prime) / np.log(
            (W_prime / (1 - W_prime) + 1 / 2**0.5) / (W_prime / (1 - W_prime) - 1 / 2**0.5)
        )
        
    else:
        S = np.sqrt(R1**2 + 1) / (R1 - 1)
        W = ((1 - P1 * R1) / (1 - P1))**(1 / N)
        F = S * np.log(W) / np.log((1 + W - S + S * W) / (1 + W + S - S * W))
    
    if np.isnan(F):
        F = 0
    
    return F

def f_lmtd2(R, P, params, C_r):
    R = abs(R)
    P = abs(P)

    R = min(max(R, 1e-4), 100)

    # Precompute frequently accessed values
    flow_type = params['Flow_Type']
    n_series = params.get('n_series', 1)
    tube_pass = params.get('Tube_pass', 1)

    # Limit P based on R
    if R < 0.41:
        P_max = 0.99999
    else:
        num = 0.995 * (0.966475661367996*R**3 - 1.431274296912407*R**2 +
                       0.247230065033875*R + 0.309118607270897)
        den = (R**4 - 1.766309686745371*R**3 + 1.287179055148762*R**2 -
               0.902512766020191*R + 0.484880205333508)
        P_max = num / den

    P = min(max(P, 0), P_max)

    # Compute epsilon
    if R <= 1:
        epsilon = P
        Pbis, Rbis = P, R
    else:
        epsilon = P * R
        Pbis, Rbis = P * R, 1 / R

    # Select correlation
    if flow_type == 'Shell&Tube':

        Cr = C_r
        sqrt_term = np.sqrt(1 + Cr*Cr)
        fact1 = 1 + Cr + sqrt_term

        def eps_g(NTU):
            NTU = max(NTU, 1e-4)
            fact2 = 1 + np.exp(-NTU*sqrt_term)
            fact3 = 1 - np.exp(-NTU*sqrt_term)

            eps1 = 2.0 / (fact1 * (fact2/fact3))
            if n_series == 1:
                return eps1
            f4 = (1 - eps1*Cr) / (1 - eps1)
            return (f4**n_series - 1) / (f4**n_series - Cr)

    elif flow_type == 'CrossFlow':
        def eps_g(NTU):
            NTU = max(NTU, 1e-4)
            return 1 - np.exp((1/C_r)*NTU**0.22*(np.exp(-C_r*NTU**0.78)-1))

    elif flow_type == 'ParallelFlow':
        def eps_g(NTU):
            NTU = max(NTU, 1e-4)
            return (1 - np.exp(-NTU*(1+C_r))) / (1+C_r)

    # Solve epsilon = eps_g(NTU)
    def root_fn(NTU):
        return epsilon - eps_g(NTU)

    try:
        NTU = brentq(root_fn, 1e-4, 50)   # NTU max = 50 is safe; can adjust
    except ValueError:
        NTU = 1.0  # fallback
        return 1.0

    # Compute F
    if R == 1:
        F = P/NTU/(1-P)
    elif P <= 1e-4:
        F = 1.0
    else:
        F = epsilon * np.log((Pbis-1)/(Rbis*Pbis - 1)) / (NTU*Pbis*(Rbis-1))

    # Special Shell&Tube case
    if flow_type == 'Shell&Tube' and tube_pass > 1 and P < 0.1:
        F = 1.0

    return F

if __name__ == "__main__":
    
    params = {
        'Flow_Type' : 'Shell&Tube',
        'Tube_pass' : 2,
        'n_series' : 1
        }
    
    F = F_shell_and_tube(0.2,10,1)
    F2 = f_lmtd2(0.6490414360002104,0.8000269409831248,params,1)

    print(F)

    print(F2)

    R_vec = np.arange(0.1, 10, 0.1)
    P_vec = np.arange(0.1, 1, 0.01)

    F_mat = np.zeros([len(R_vec), len(P_vec)])

    for i in range(len(R_vec)):
        for j in range(len(P_vec)):
            F_mat[i][j] = F = F_shell_and_tube(R_vec[i],P_vec[j],1)
            
