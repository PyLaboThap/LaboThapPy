import matplotlib.pyplot as plt
import numpy as np

from labothappy.component.expander.turbine_polyneff import TurbinePolynEff

p_in = np.array([10, 15, 20])  # 0.8
T_in = np.array([173, 195.6, 213.1])
m_dot_tune = np.array([4.74, 7.02, 9.29])

m_dot_model = np.zeros(len(T_in))

# %%
"--------- 1) DECAGONE Turbine ------------------------------------------------------------------------------------------"

"Create turb object"
Turb = TurbinePolynEff()

"Set params"
Turb.set_parameters(
    D_inlet=0.5,
    N_turb_rated=1500,
    turb_voltage=400,
    turb_phases=3,
    eta_max_motor=95,
    W_dot_el_rated=2e6,
    eta_m=50,
    eta_is_coefs=[0.1] * 10,
    eta_is_coefs_red=[0.1] * 10,
    A_th=0.05,
)

for i in range(len(p_in)):  # Loop for model flowrate verification (verify the input vs the ouput flowrate)

    "Set inputs"
    Turb.set_inputs(
        fluid="Cyclopentane",
        T_su=273.15 + T_in[i],
        P_su=p_in[i] * 1e5,
        P_ex=0.8 * 1e5,
        N_rot=50,
    )

    Turb.solve()
    m_dot_model[i] = Turb.m_dot

plt.figure()
plt.plot(m_dot_tune, m_dot_model)
plt.plot([0, 20], [0, 20])
plt.grid()
plt.show()

Turb.W.print_resume()


"""

Outputs shall be : 

"""
