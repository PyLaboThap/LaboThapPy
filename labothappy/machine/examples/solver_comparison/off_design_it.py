import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.expander.expander_semi_empirical import ExpanderSE
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity
from labothappy.toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP

class TS_curve_generator:
    
    def __init__(self, Fluid):
        self.Fluid = Fluid
        self.TS_curve()
        
    def TS_curve(self):
        # Calculate critical pressure and minimum pressure
        P_crit = 0.99 * PropsSI('Pcrit', 'T', 0, 'P', 0, self.Fluid)
        P_min = PropsSI('P', 'T', 273.15 - 10, 'Q', 0, self.Fluid)

        # Create an array of pressures
        P_lin = np.linspace(P_min, P_crit, 100)
        
        # Initialize arrays to store data
        s_liq = np.zeros(100)
        T_liq = np.zeros(100)
        s_vap = np.zeros(100)
        T_vap = np.zeros(100)

        # Loop through the pressures and collect data
        for i in range(100):
            s_liq[i] = PropsSI('S', 'P', P_lin[i], 'Q', 0, self.Fluid)
            T_liq[i] = PropsSI('T', 'P', P_lin[i], 'Q', 0, self.Fluid)
            s_vap[i] = PropsSI('S', 'P', P_lin[i], 'Q', 1, self.Fluid)
            T_vap[i] = PropsSI('T', 'P', P_lin[i], 'Q', 1, self.Fluid)

        # Create the full TS curve by combining liquid and vapor data
        self.s_TS_curve = np.concatenate((s_liq, np.flip(s_vap)))
        self.T_TS_curve = np.concatenate((T_liq, np.flip(T_vap)))
        plt.plot(self.s_TS_curve, self.T_TS_curve, color='black', linestyle='--')
        plt.xlabel('Entropy (J/kg K)')
        plt.ylabel('Temperature (K)')
        # plt.show()

    def points(self, s_array, T_array):

        plt.plot(s_array, T_array, 'bo-', markersize=4)
        # plt.plot(self.s_TS_curve, self.T_TS_curve, color='blue')
        # plt.xlabel('Entropy (J/kg K)')
        # plt.ylabel('Temperature (K)')
        # plt.title('TS Curve')
        plt.show()


# -------- 1) Instanciate Circuit --------
fluid = 'R1233zd(E)'
orc = IterativeCircuit(fluid)

# -------- 2) Create components --------
Expander = ExpanderSE()
Condenser = HexMBChargeSensitive('Plate')
Pump = PumpCurveSimilarity()
Evaporator = HexMBChargeSensitive('Plate')

# -------- 3) Set component parameters --------
# Expander
Expander.set_parameters(AU_amb=9.3, AU_su_n=4.75, AU_ex_n=17.7, d_su1=6.48e-3, m_dot_n=0.1, 
            A_leak=9.99e-06, W_dot_loss_0=2.37e+1 , alpha= 1.16e-1, C_loss=1.13, rv_in=1.7, V_s=2*0.0000712, mode='P_M')
# I take into account the fact that there's two expander by multiplying by two the Vs

# Evaporator
evaporator_geom = PlateGeomSWEP()
evaporator_geom.set_parameters("P200THx140/1P_Evaporator")

Evaporator.set_parameters(
    # Set the geometry of the evaporator
    A_c=evaporator_geom.A_c, A_h=evaporator_geom.A_h, h=evaporator_geom.h, l=evaporator_geom.l, l_v=evaporator_geom.l_v, w_v=evaporator_geom.w_v,
    C_CS=evaporator_geom.C_CS, C_Dh=evaporator_geom.C_Dh, C_V_tot=evaporator_geom.C_V_tot, C_canal_t=evaporator_geom.C_canal_t, C_n_canals=evaporator_geom.C_n_canals, 
    H_CS=evaporator_geom.H_CS, H_Dh=evaporator_geom.H_Dh, H_V_tot=evaporator_geom.H_V_tot, H_canal_t=evaporator_geom.H_canal_t, H_n_canals=evaporator_geom.H_n_canals,
    casing_t=evaporator_geom.casing_t, chevron_angle=evaporator_geom.chevron_angle, fooling=evaporator_geom.fooling, 
    n_plates=evaporator_geom.n_plates, plate_cond=evaporator_geom.plate_cond, plate_pitch_co=evaporator_geom.plate_pitch_co, t_plates=evaporator_geom.t_plates, w=evaporator_geom.w, 
    amplitude=evaporator_geom.amplitude, phi=evaporator_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0,
)

Evaporator.set_htc(htc_type="Correlation", Corr_H={"1P": "water_plate_HTC"}, Corr_C={"1P": "martin_holger_plate_HTC", "2P": "amalfi_plate_HTC"}) # 'User-Defined' or 'Correlation' # 28
Evaporator.set_DP()


# Condenser
condenser_geom = PlateGeomSWEP()
condenser_geom.set_parameters("P200THx140/1P_Condenser")

Condenser.set_parameters(
    # Set the geometry of the condenser
    A_c=condenser_geom.A_c, A_h=condenser_geom.A_h, h=condenser_geom.h, l=condenser_geom.l, l_v=condenser_geom.l_v, w_v=condenser_geom.w_v,
    C_CS=condenser_geom.C_CS, C_Dh=condenser_geom.C_Dh, C_V_tot=condenser_geom.C_V_tot, C_canal_t=condenser_geom.C_canal_t, C_n_canals=condenser_geom.C_n_canals,
    H_CS=condenser_geom.H_CS, H_Dh=condenser_geom.H_Dh, H_V_tot=condenser_geom.H_V_tot, H_canal_t=condenser_geom.H_canal_t, H_n_canals=condenser_geom.H_n_canals,
    casing_t=condenser_geom.casing_t, chevron_angle=condenser_geom.chevron_angle, fooling=condenser_geom.fooling,
    n_plates=condenser_geom.n_plates, plate_cond=condenser_geom.plate_cond, plate_pitch_co=condenser_geom.plate_pitch_co, t_plates=condenser_geom.t_plates, w=condenser_geom.w,
    amplitude=condenser_geom.amplitude, phi=condenser_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=10,
)

Condenser.set_htc(htc_type="Correlation", Corr_H={"1P": "martin_holger_plate_HTC", "2P": "shah_condensation_plate_HTC"}, Corr_C={"1P": "water_plate_HTC"}) # 'User-Defined' or 'Correlation' # 28
Evaporator.set_DP()

# Pump
V_dot_curve = np.array([20, 30, 40, 50, 60, 70, 80])   # m3/h
Delta_H_curve = np.array([57, 55, 52, 49, 45, 42, 36])  # m (head falls with flow)
eta_is_curve = np.array([0.45, 0.59, 0.69, 0.75, 0.79, 0.79, 0.75])  # eff peaks near mid-flow
NPSH_r_curve = np.array([1.1, 1.1, 1.4, 1.8, 2, 3, 4.7])  # m, increases again near max flow
N_rated = 2900 # RPM
# Reference point: water at 20°C and 1 atm for a rated speed of 2900 RPM

Pump.set_parameters(
    V_dot_curve = V_dot_curve,
    Delta_H_curve = Delta_H_curve,
    eta_is_curve = eta_is_curve,
    NPSH_r_curve = NPSH_r_curve,
    N_rot_rated = N_rated,
    mode = "P_M",  # Mode can be "M_N", "P_M", or "P_N"
)

# -------- 4) Add components to circuit --------
orc.add_component(Expander, "Expander")
orc.add_component(Condenser, "Condenser")
orc.add_component(Pump, "Pump")
orc.add_component(Evaporator, "Evaporator")

# -------- 5) Link components with mass connectors --------
orc.link_components("Expander", "m-ex", "Condenser", "m-su_H")
orc.link_components("Condenser", "m-ex_H", "Pump", "m-su")
orc.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
orc.link_components("Evaporator", "m-ex_C", "Expander", "m-su")

# -------- 6) Add fluid sources --------
CD_source = MassConnector('Water')
T_su_w_cd = 5 +273.15
P_su_w_cd = 2e5
m_dot_w_cd = 3  # kg/s
EV_source = MassConnector('Water')
T_su_w_ev = 70+273.15
P_su_w_ev = 2e5
m_dot_w_ev = 2.5  # kg/s

orc.add_source("CD_Water", CD_source, orc.components["Condenser"], "m-su_C")
orc.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
orc.add_source("EV_Water", EV_source, orc.components["Evaporator"], "m-su_H")
orc.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="EV_Water")

# -------- 7) Set cycle inputs --------
m_dot_ref = 0.4 # kg/s
SC_cd = 5  # K
N_exp = 6000 # RPM
T_amb = 293 # K

orc.set_cycle_input(target="Pump:su", m_dot = m_dot_ref, SC=SC_cd)
orc.set_cycle_input(target="Expander:W", N_rot = N_exp)
orc.set_cycle_input(target="Expander:Q_amb", T_amb=T_amb)

# Expander.print_setup()

# -------- 8) Set iteration variables --------
# P_pp_su_lb = max(PropsSI('P', 'Q', 0, 'T', T_su_w_cd, fluid), PropsSI('P_min', 'Q', 0, 'T', 273.15, fluid))
# P_pp_ex_ub = PropsSI('P', 'Q', 0, 'T', min(PropsSI('Tcrit', 'Q', 0, 'T', 273.15, fluid)-2, T_su_w_ev-1), fluid)
# rp_max = P_pp_ex_ub/P_pp_su_lb
# rp_min = min(1.01, rp_max)   #min(1.01, rp_max)
# Possible that the guesses are ill conditioned
# P_pp_su_ub = P_pp_ex_ub/rp_min
# -> Si jamais besoin d'essayer plusieurs Guesses

P_LP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
P_HP_guess = PropsSI("P", "T", T_su_w_ev-10, "Q", 1, fluid)

orc.set_iteration_variable(
    target=["Pump:su", "Expander:ex"],
    variable="p",
    guess=P_LP_guess,
    tolerance=1e-6
)

orc.set_iteration_variable(
    target="Pump:ex",
    variable="p",
    guess=P_HP_guess,
    tolerance=1e-6
)

# -------- 9) Set residual variables --------
orc.set_residual_variable(
    pre_target="Pump:su",
    post_target="Condenser:ex_H",
    variable="h",
    tolerance=1e-3
)

orc.set_residual_variable(
    pre_target="Expander:W",
    post_target="Expander:W",
    variable="N_rot",
    tolerance=1e-3
)

# -------- 10) Solve circuit --------
orc.solve()
print(f"Converged at P_HP = {Expander.su.p}, P_LP = {Expander.ex.p}")
orc.print_states()

"Graphs"
# Create array with point of the cycle
T_cd = PropsSI('T', 'P', Condenser.su_H.p, 'Q', 0.5, fluid)
T_ev = PropsSI('T', 'P', Evaporator.su_C.p, 'Q', 0.5, fluid)
s_cd_1 = PropsSI('S', 'P', Condenser.su_H.p, 'Q', 1, fluid)
s_cd_0 = PropsSI('S', 'P', Condenser.su_H.p, 'Q', 0, fluid)
s_ev_1 = PropsSI('S', 'P', Evaporator.su_C.p, 'Q', 1, fluid)
s_ev_0 = PropsSI('S', 'P', Evaporator.su_C.p, 'Q', 0, fluid)

T_array = [Expander.su.T, Condenser.su_H.T, T_cd, T_cd, Condenser.ex_H.T, Evaporator.su_C.T, T_ev, T_ev, Evaporator.ex_C.T, Expander.su.T]
s_array = [Expander.su.s, Condenser.su_H.s, s_cd_1, s_cd_0, Condenser.ex_H.s, Evaporator.su_C.s, s_ev_0, s_ev_1, Evaporator.ex_C.s, Expander.su.s]

T_c = np.linspace(Condenser.ex_C.T, Condenser.su_C.T, 100)
s_c = np.linspace(Condenser.su_H.s, Condenser.ex_H.s, 100)
T_h = np.linspace(Evaporator.su_H.T, Evaporator.ex_H.T, 100)
s_h = np.linspace(Evaporator.ex_C.s, Evaporator.su_C.s, 100)

plt.plot(s_h, T_h, color='blue', linestyle='-', label='Heat sink')
plt.plot(s_c, T_c, color='red', linestyle='-', label='Heat source')

# Add numbers next to points
for i, (s, T) in enumerate(zip(s_array, T_array)):
    plt.text(s, T, str(i+1), fontsize=12, color='black', ha='right', va='bottom')


TS_curve = TS_curve_generator(fluid)
TS_curve.points(s_array, T_array)




