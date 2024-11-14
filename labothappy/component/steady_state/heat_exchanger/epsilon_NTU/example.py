from library.component.steadystate.heat_exchanger.epsilon_NTU.simulation_model import HXeNTU
import numpy as np

HX = HXeNTU()

HX.set_inputs(
    # First fluid
    Hsu_fluid = 'Cyclopentane',
    Hsu_T = 205 + 273.15, # K
    Hsu_p = 1*1e5, # Pa
    Hsu_m_dot = 0.014, # kg/s

    # Second fluid
    Csu_fluid = 'Water',
    Csu_T = 12 + 273.15, # K
    Csu_p = 4*1e5, # Pa
    Csu_m_dot = 0.08, # kg/s  # Make sure to include fluid information
)

 
A_htx = 0.752 # m^2

# HTX dimensions
L = 0.393 # [m] : length
w = 0.243 # [m] : width
h = 0.0446 # [m] : height
l_v = 0.324 # [m] : length between ports
casing_t = 0.005 # [m] : casing thickness # !!! arbitrary

# Number and thickness of plates
n_plates = 10
t_plates = 0.0008 # [m] # !!! arbitrary

# Fooling factor
fooling = 98.51/1000 # (m^2*K)/W

# Number of canals
C_n_canals = 5
H_n_canals = 4

# Plate values 
plate_cond = 45 # [W/(m*K)] : plate conduction
plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
chevron_angle = 20*np.pi/180 # !!! arbitrary

# Total volume of each part
V_tot = 0.9*1e-3 # [m^3]

# Canal thickness
C_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*C_n_canals)
H_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*H_n_canals)

# Total Canal Surface
C_CS = C_canal_t*(w-2*casing_t)*C_n_canals
H_CS = H_canal_t*(w-2*casing_t)*H_n_canals

# Dh : hydraulic diameter
C_Dh = (4*C_canal_t*w)/(2*C_canal_t+2*w)
H_Dh = (4*H_canal_t*w)/(2*H_canal_t+2*w)
    
HX.set_parameters(
    A_htx=0.752, L_HTX=0.393, V_HTX=0.9*1e-3, Flow_Type = 'CounterFlow',
    A_canal_h=H_CS, A_canal_c=C_CS, D_h=H_Dh, 
    k_plate=45, t_plate=0.0008, n_plates = 10,
    co_pitch=0.005, chevron_angle=20*np.pi/180, fouling=98.51/1000
)

# Solve the expander component
HX.solve()
HX.print_setup()
HX.print_results()
