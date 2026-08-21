




# ============================================================================
# ONE PHASE CORRELATIONS
# ============================================================================

# TO IMPLEMENT: LONGO??

# For the R1233zd(E) one phase:
# def martin_holger_plate_HTC(mu, Pr, k, m_dot, nb_channels, T_mean, P_mean, fluid, D_h, length, width, amplitude, chevron_angle):
#     """
#     Reference:
#     -----------
#     H. M. Hofmann, M. Kind, and H. Martin, ‘Measurements on steady state heat transfer and flow structure and new 
#     correlations for heat and mass transfer in submerged impinging jets’, International Journal of Heat and Mass Transfer, 
    
#     Correlation from the VDI Atlas
#     """

#     " Thermodynamics properties:"
#     # mu =  PropsSI('V','P', P_mean,'T',T_mean, refrigerant)          # Viscosity, Pa s
#     rho = PropsSI('D','P', P_mean,'T',T_mean, fluid)          # Density, kg m^-3
#     # cp = PropsSI('C','P', P_mean,'T',T_mean, refrigerant)           # Specific Heat, J/kg-K
    
    

#     " Mass flow rate per plate"
#     m_dot_ch = m_dot/nb_channels                                      # Mass flow rate in the chanels, kg s^-1
    
#     " Mass velocity for chanel CHEK!!!!"
#     w_ch = m_dot_ch/(amplitude*width*rho)                                   # Mass velocity in chanels, kg m^2 s^-1
#     Re = rho*w_ch*D_h/mu                                          # Reynolds Number, -

#     " Factor for correlations: provided by Focke et al."
#     if Re >= 2000:                                                  # Regimen: Turbulent
#         xhi_0   = (1.8*np.log(Re)-1.5)**-2
#         xhi_1_0 = 39/Re**0.289
#     elif Re <2000:                                                  # Regime: Laminar
#         xhi_0   = 64/Re
#         xhi_1_0 = 597/Re +3.85
    
#     " Constant given by Martin"
#     a = 3.8
#     b = 0.18 
#     c = 0.36    
    
#     " Factor xhi" 
#     xhi_1 = a*xhi_1_0

#     " Beta angle from degree to radians"
#     chevron_angle = 90*np.pi/180 - chevron_angle
#     beta_r = chevron_angle     #(90*np.pi/180)-chevron_angle                                   # Chevron angle, in Radians   
     
#     " Friction factor equation 18"    
#     f = (np.cos(beta_r)/np.sqrt(b*np.tan(beta_r) + c*np.sin(beta_r) + xhi_0/np.cos(beta_r)) +(1 - np.cos(beta_r))/np.sqrt(xhi_1))**(-2)
 
#     " Hagen number "
#     Hg = f*Re**2/2

#     " Pressure Drop:"    
#     DeltaP = Hg * (mu**2*length)/(rho*D_h**3)
    
#     " Extracted from the comparison with Heavear et al. [10]"
#     c_q = 0.122
#     q = 0.374
    
#     # " Wall temperature "
#     # T_wall = T_mean-5                                                # Wall Temperature, K
#     # if T_wall<273.15:
#     #     T_wall = 274.15
        
#     # mu_w =  PropsSI('V','P', P_mean,'T',T_wall, refrigerant)          # Viscosity at wall T, Pa s
    
#     " Nusslet number: "
#     # Nu = c_q*Pr**(1/3)*(mu/mu_w)**(1/6)*(2*Hg*np.sin(2*beta_r))**q
#     Nu = c_q*Pr**(1/3)*(2*Hg*np.sin(2*beta_r))**q
#     h_conv = Nu*k/D_h
#     # print(f' Nu: {Nu}')
#     # print(f' mu_w: {mu_w}')
#     # print(f' T_wall : {T_wall}')
#     # print(f' Result: {(mu/mu_w)**(1/6)}')
    
#     return h_conv




