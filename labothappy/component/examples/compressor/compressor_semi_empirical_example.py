from labothappy.component.compressor.compressor_semi_empirical import CompressorSE

"Example of the semi-empirical model compressor component"
# Inputs: N_rot, T_amb, P_su, h_su, P_ex, fluid 
#         OR m_dot, T_amb, P_su, h_su, P_ex, fluid
# Parameters: AU_amb, AU_su_n, AU_ex_n, d_su1, m_dot_n, A_leak, W_dot_loss_0, alpha, C_loss, 
# rv_in, V_s, mode

# Create an instance of the expander component
compressor = CompressorSE()

#-------------------------------------------------------------------------------------------#
"First case: m_dot known and N_rot unknown"
#-------------------------------------------------------------------------------------------#
# "Specify m_dot as an input"
compressor.set_parameters(
    AU_amb=9.96513290e+00, AU_su_n=1.02359773e+01, AU_ex_n=2.24133147e+00, d_ex=1.82304791e-02, m_dot_n=0.1, 
    A_leak=3.66336680e-07, W_dot_loss_0=9.05482168e-01, alpha=3.22395090e-03, C_loss=1.11169710e-061, rv_in=1.7,
    V_s=1.17889079e-04, mode = 'm_dot'
)

# "1. Inputs set through connectors"
# # Set properties for su connector
# compressor.su.set_properties(P=319296.5575177148, T=331.033964665788, m_dot=0.2, fluid='R1233zd(E)')
# compressor.ex.set_properties(P=606240.1433176235)

# # Set ambient temperature
# compressor.Q_amb.set_T_cold(293)


# "2. Inputs set directly"
# compressor.set_inputs(
#     # N_rot=6000,
#     m_dot=0.2,
#     T_amb=293,
#     P_su=319296.5575177148,
#     T_su=331.033964665788,
#     P_ex=606240.1433176235,
#     fluid='R1233ZDE'  # Make sure to include fluid information
# )

#-------------------------------------------------------------------------------------------#
"Second case: N_rot known and m_dot unknown"
#-------------------------------------------------------------------------------------------#
"Specify N_rot as an input"
# compressor.set_parameters(
#     AU_amb=9.96513290e+00, AU_su_n=1.02359773e+01, AU_ex_n=2.24133147e+00, d_ex=1.82304791e-02, m_dot_n=0.1, 
#     A_leak=3.66336680e-07, W_dot_loss_0=9.05482168e-01, alpha=3.22395090e-03, C_loss=1.11169710e-061, rv_in=1.7,
#     V_s=1.17889079e-04, mode = 'N_rot'
# )

"1. Inputs set through connectors"
# Set properties for su connector
compressor.su.set_properties(P=319296.5575177148, T=331.033964665788, m_dot=0.2, fluid='R1233zd(E)')
compressor.ex.set_properties(P=606240.1433176235)

# Set ambient temperature
compressor.Q_amb.set_T_cold(293)

# Set rotational speed
compressor.W.set_N_rot(6000)



#---------------------------------------------------------------------------------------------#
"Solve compressor component"
#---------------------------------------------------------------------------------------------#
# compressor.print_setup()
# Solve the compressor component
compressor.solve()
compressor.print_results()