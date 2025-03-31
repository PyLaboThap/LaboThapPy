from labothappy.component.volumetric_machine.compressor.steady_state.semi_empirical.calibration.calibration import calibrate_compressor_model
import pandas as pd
# from labothappy.component.volumetric_machine.compressor.steady_state.semi_empirical.calibration.calibration import plot_parity
data_path = "labothappy/component/volumetric_machine/compressor/steady_state/semi_empirical/calibration/Exp_data_cp.csv"

# Load data
data = pd.read_csv(data_path)
    
# Extract columns based on mapping
P_su = data['P_su [Pa]']
P_ex = data['P_ex [Pa]']
T_ex = data['T_ex [°C]']+273.15
N_cp_RPM = data['RPM']
SH = data['SH [K]']
W_dot_sh_meas = data['W_dot_shaft [W]']
m_dot_meas = data['m_dot_wf [kg/s]']

data = {
    "P_su": P_su,
    "P_ex": P_ex,
    "SH": SH,
    "N_cp_RPM": N_cp_RPM,
    "SH": SH,
    "W_dot_sh_meas": W_dot_sh_meas,
    "m_dot_meas": m_dot_meas,
    "T_ex_meas": T_ex
}

initial_params = [10, 5, 5, 6.5e-3, 2e-7, 1, 0.05, 0.5] # AU_amb, AU_su_n, AU_ex_n, d_ex, A_leak, W_dot_loss_0, alpha, C_loss = params_opt
bounds = ([0.1, 0.1, 0.1, 0.1e-3, 1e-10, 0.001, 0.001, 0.001], [100, 100, 100, 20e-3, 1e-5, 1000, 0.999, 50])
fitting_params = {
    "rv_in": 1.7,
    "V_s": 0.000121,
    "T_amb": 293,
    "m_dot_n": 0.1
}

result = calibrate_compressor_model(data, initial_params, bounds, fitting_params, fluid='R134a')

print("Optimal Parameters:", result["optimal_params"])

# # Example usage of the plot_parity function
# AU_amb, AU_su_n, AU_ex_n, d_su, A_leak, W_dot_loss_0, alpha, C_loss = result["optimal_params"]
