from labothappy.component.volumetric_machine.expander.steady_state.semi_empirical.calibration_fct.calibration import calibrate_expander_model
from labothappy.component.volumetric_machine.expander.steady_state.semi_empirical.calibration_fct.calibration import plot_parity
data_path = "labothappy/component/volumetric_machine/expander/steady_state/semi_empirical/calibration_fct/data_example.csv"

col_names = {
    'P_su': 'P_su',
    'rp': 'Press.Ratio',
    'N_exp_RPM': 'RPM',
    'T_su': 'T_su',
    'W_dot_sh_meas': 'W_shaft',
    'm_dot_meas': 'm_dot_wf',
}

initial_params = [10, 5, 5, 6.5e-3, 2e-7, 1, 0.05, 0.5] # AU_amb, AU_su_n, AU_ex_n, d_su, A_leak, W_dot_loss_0, alpha, C_loss = params_opt
bounds = ([0.1, 0.1, 0.1, 0.1e-3, 1e-10, 0.001, 0.001, 0.001], [100, 100, 100, 20e-3, 1e-5, 1000, 0.999, 50])
fitting_params = {
    "rv_in": 1.7,
    "V_s": 0.0000712,
    "T_amb": 293,
    "m_dot_n": 0.1
}

result = calibrate_expander_model(data_path, col_names, initial_params, bounds, fitting_params, fluid='R134a')

print("Optimal Parameters:", result["optimal_params"])

# Example usage of the plot_parity function
AU_amb, AU_su_n, AU_ex_n, d_su, A_leak, W_dot_loss_0, alpha, C_loss = result["optimal_params"]
