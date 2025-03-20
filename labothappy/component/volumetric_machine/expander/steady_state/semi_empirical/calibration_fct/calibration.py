import pandas as pd
from scipy.optimize import least_squares

from labothappy.component.volumetric_machine.expander.steady_state.semi_empirical.exp_semi_empirical import ExpanderSE
from toolbox.parity_plot import plot_parity
# from labothappy.connector. mass_connector import MassConnector

import pandas as pd
from scipy.optimize import least_squares
from labothappy.component.volumetric_machine.expander.steady_state.semi_empirical.exp_semi_empirical import ExpanderSE
from toolbox.parity_plot import plot_parity

def calibrate_expander_model(data_path, col_names, initial_params, bounds, fitting_params, fluid):
    """
    Generic function to calibrate the expander model using experimental data.
    
    Parameters:
    - data_path (str): Path to the CSV file containing the experimental data.
    - col_names (dict): Dictionary mapping parameter names to column names in the CSV.
                        Example: {'P_su': 'P_su', 'rp': 'Press.Ratio', ...}
    - initial_params (list): Initial guesses for the parameters to optimize.
    - bounds (tuple): Bounds for the optimization in the format (lower_bounds, upper_bounds).
    - fitting_params (dict): Dictionary of known/fixed parameters for the model.
    
    Returns:
    - dict: Contains optimized parameters and fitting results.
    """
    # Load data
    data = pd.read_csv(data_path)
    
    # Extract columns based on mapping
    P_su = data[col_names['P_su']]
    rp = data[col_names['rp']]
    P_ex = P_su / rp
    N_exp_RPM = data[col_names['N_exp_RPM']]
    T_su = data[col_names['T_su']]
    W_dot_sh_meas = data[col_names['W_dot_sh_meas']]
    m_dot_meas = data[col_names['m_dot_meas']]
    
    # Shared storage for calculated values
    results = {
        "m_dot_calc": [],
        "W_dot_sh_calc": []
    }
    
    # Objective function for optimization
    def objective(params_opt):
        # Unpack parameters to be optimized
        AU_amb, AU_su_n, AU_ex_n, d_su, A_leak, W_dot_loss_0, alpha, C_loss = params_opt
        
        # Clear previous calculations
        results["m_dot_calc"].clear()
        results["W_dot_sh_calc"].clear()
        
        diff = []  # Residuals
        for i in range(len(P_su)):
            # Initialize Expander object
            expander = ExpanderSE()
            expander.set_parameters(
                AU_amb=AU_amb, AU_su_n=AU_su_n, AU_ex_n=AU_ex_n, d_su1=d_su, m_dot_n=0.1, 
                A_leak=A_leak, W_dot_loss_0=W_dot_loss_0, alpha=alpha, C_loss=C_loss, rv_in=1.7, V_s=0.0000712
            )
            
            # Set fluid and connectors
            expander.su.set_fluid(fluid)
            expander.ex.set_fluid(fluid)
            expander.su.set_p(P_su[i])
            expander.su.set_T(T_su[i])
            expander.ex.set_p(P_ex[i])
            expander.W_exp.set_N(N_exp_RPM[i])
            expander.Q_amb.set_T_cold(293)
            
            # Solve
            expander.solve()
            
            # Store results
            results["m_dot_calc"].append(expander.m_dot)
            results["W_dot_sh_calc"].append(expander.W_dot_exp)
            
            # Residuals for least squares
            diff.append(((expander.W_dot_exp - W_dot_sh_meas[i]) / W_dot_sh_meas[i]) *
                        ((expander.m_dot - m_dot_meas[i]) / m_dot_meas[i]))
        
        return diff
    
    # Perform least-squares optimization
    result = least_squares(objective, initial_params, bounds=bounds)
    
    # Plot parity for mass flow rate
    plot_parity(
        measured=m_dot_meas,
        calculated=results["m_dot_calc"],
        label_measured='Measured Mass Flow',
        label_calculated='Calculated Mass Flow',
        title='Parity Plot: Mass Flow Rate',
        x_label=r'$\dot{m}_{meas} [kg/s]$',
        y_label=r'$\dot{m}_{calc} [kg/s]$',
        file_name='Calibration_m_dot'
    )
    
    # Plot parity for shaft work
    plot_parity(
        measured=W_dot_sh_meas,
        calculated=results["W_dot_sh_calc"],
        label_measured='Measured Shaft Work',
        label_calculated='Calculated Shaft Work',
        title='Parity Plot: Shaft Work',
        x_label=r'$\mathrm{\dot{W}_{sh, meas}} [W]$',
        y_label=r'$\mathrm{\dot{W}_{sh, calc}} [W]$',
        file_name='Calibration_W_dot'
    )
    
    return {
        "optimal_params": result.x,
        "optimization_result": result,
        "calculated_values": results  # Return calculated values for further analysis
    }
