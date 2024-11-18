import matplotlib.pyplot as plt

def plot_parity(measured, calculated, label_measured, label_calculated, title, x_label, y_label, file_name):
    """
    Generate a parity plot comparing measured and calculated data.

    Parameters:
    - measured: List or array of measured values.
    - calculated: List or array of calculated values.
    - label_measured: String label for measured data points.
    - label_calculated: String label for calculated data points.
    - title: Title of the plot.
    - x_label: X-axis label.
    - y_label: Y-axis label.
    - file_name: File name for saving the plot.
    """
    fig, ax = plt.subplots(figsize=(5.5, 4.3), constrained_layout=True)

    # 45° line
    plt.plot([min(measured), max(measured)], [min(measured), max(measured)], 'k', linewidth=1.7, label=r'45° line')

    # 10% error lines
    plt.plot([min(measured), max(measured)], [min(measured) * 1.1, max(measured) * 1.1], '--', linewidth=1.4, color='gray', label=r'10% error')
    plt.plot([min(measured), max(measured)], [min(measured) * 0.9, max(measured) * 0.9], '--', linewidth=1.4, color='gray')

    # Data points
    plt.scatter(measured, calculated, marker='*', label=label_calculated)

    # Labels and legend
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.title(title, fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)

    # Save plot
    plt.savefig(file_name + '.eps', format='eps', bbox_inches='tight')
    plt.savefig(file_name + '.svg', format='svg', bbox_inches='tight')
    plt.show()

# Example usage of the plot_parity function
# plot_parity(
#     measured=m_dot_meas,
#     calculated=m_dot_calc,
#     label_measured='Measured Mass Flow',
#     label_calculated='Calculated Mass Flow',
#     title='Parity Plot: Mass Flow Rate',
#     x_label=r'$\dot{m}_{meas} [kg/s]$',
#     y_label=r'$\dot{m}_{calc} [kg/s]$',
#     file_name='Calibration_m_dot'
# )

# plot_parity(
#     measured=W_dot_sh_meas,
#     calculated=W_dot_sh_calc,
#     label_measured='Measured Shaft Work',
#     label_calculated='Calculated Shaft Work',
#     title='Parity Plot: Shaft Work',
#     x_label=r'$\mathrm{\dot{W}_{sh, meas}}  [W]$',
#     y_label=r'$\mathrm{\dot{W}_{sh, calc}} [W]$',
#     file_name='Calibration_W_dot'
# )
