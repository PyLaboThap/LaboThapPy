
import numpy as np
from labothappy.component.compressor.compressor_radial_mean_line import CompressorRadialMeanLine

#%%
  
Comp = CompressorRadialMeanLine()

Comp.set_inputs(
    fluid = 'CO2',
    m_dot  = 2.15,
    T_su = 305.97,
    P_su = 76.9*1e5,
    N_rot  = 50000,
)

Comp.set_parameters(
    alpha1_des = 1.063  *180/np.pi,
    n_blade_R = 15, 
    t_b = 0.762*1e-3,
    
    b1 = 0.00684,
    b2 = 0.00216,
    b3 = 0.00216,        
    b5 = 0.00216,
    
    CP           = 0.44,

    eps_imp      = 0.254*1e-3,
    eps_bf_imp   = 0.254*1e-3,
    k_imp        = 0.01*1e-3,
            
    L_z = 0.1137,
    
    r1 = 0.0077,
    r1s = 0.011,
    r1h = 0.0042,
    r2 = 0.0217, 
    r3 = 0.0396,
    r5 = 0.0494,
    
    xhi1 = 46.89,
    xhi2 = 43.44,
    xhi4 = 1.433    *180/np.pi,
    xhi5 = 0.00036  *180/np.pi,
    )

Comp.solve()

Comp.print_results()
