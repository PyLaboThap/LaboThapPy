
from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline

class PT_Collector_Geom(object):
    """
    This class aims to store geometrical data of parabolic through solar concentrators
    """
     
    def __init__(self, **kwargs):

        "Heat transfer area"
        self.L = None # Hot side
        self.A_c = None # Cold side
        
        "HX Volume"
        self.H_V_tot = None # Hot side volume
        self.C_V_tot = None # Cold side volume
        
        "HX Dimensions"
        self.l = None # length
        self.w = None # width
        self.h = None # height
        self.l_v = None # length between ports
        self.casing_t = None # casing thickness
        
        "Plate related parameters"
        self.n_plates = None # number of plates
        self.t_plates = None # plate thickness
        self.plate_cond = None # Plate conduction
        self.fooling = None # Fooling Factor
        self.plate_pitch_co = None # corrugated pitch 
        self.chevron_angle = None # chevron angle

    def set_parameters(self, name, **kwargs):
        
        if name == "Soponova_MicroCSP":
            
            def efficiency():
                # Data
                DT_C = np.array([
                    100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500,
                ]) * 5 / 9  # T_htf - T_amb (F° to C°)

                Irrad = np.array([500, 600, 700, 800, 900, 1000])  # Irradiation (W/m^2)

                eff_matrix = np.array([
                    [0.628, 0.615, 0.606, 0.593, 0.580, 0.566, 0.551, 0.535, 0.519, 0.499, 0.482, 0.461, 0.440, 0.419, 0.398, 0.375, 0.351],
                    [0.633, 0.622, 0.613, 0.604, 0.594, 0.582, 0.570, 0.558, 0.544, 0.527, 0.513, 0.497, 0.482, 0.462, 0.444, 0.426, 0.406],
                    [0.637, 0.628, 0.620, 0.611, 0.603, 0.593, 0.584, 0.573, 0.562, 0.548, 0.536, 0.522, 0.509, 0.493, 0.478, 0.460, 0.444],
                    [0.642, 0.634, 0.627, 0.621, 0.612, 0.603, 0.594, 0.584, 0.575, 0.564, 0.554, 0.543, 0.531, 0.516, 0.504, 0.490, 0.473],
                    [0.646, 0.640, 0.634, 0.628, 0.621, 0.612, 0.603, 0.596, 0.589, 0.578, 0.568, 0.557, 0.547, 0.536, 0.524, 0.509, 0.496],
                    [0.652, 0.645, 0.639, 0.635, 0.627, 0.620, 0.611, 0.605, 0.597, 0.588, 0.577, 0.567, 0.557, 0.547, 0.536, 0.523, 0.512]
                ])  # Collection efficiency (-)

                # Create the 2D interpolator using RectBivariateSpline
                fun_eff = RectBivariateSpline(Irrad, DT_C, eff_matrix)

                return fun_eff

            self.coll_eff = efficiency()

            # Collector - General Data  
            self.L = 3.657 # m : length
            self.W = 1.524 # m : width
            self.A = 5.574 # [m] : Area
            self.A_r = 5.07 # [m] : Reflective area
            self.m = 68 # [kg] : Collector weight
            self.L_f = 0.3048 # [m] : Focal Length

            # Collector Tube
            self.Tube_OD = 0.0254 # [m] : Tube external diameter
            self.Tube_V = 1.288*1e-3 # [m^3] : Tube Volume

            # Optical Parameters
            self.alpha_r = 0.92 # [-] : Receiver absorptivity
            self.refl_m = 0.91 # [-] : Mirror reflectivity
            self.epsilon_r = 0.23 # [-] : Receiver emittance @ 400°C
            self.envel_tau = 0.95 # [-] : Enveloppe transmitivity
            self.eta_other = 0.89 # [-] : Overall Optical Efficiency

            # Operational limits
            self.Vdot_min = 22.7/(60*1000) # [m^3/s] : minimum recommended flowrate
            self.Vdot_max = 45.4/(60*1000) # [m^3/s] : maximum recommended flowrate

            self.T_f_min = 50 + 273.15 # [K] : minimum operating fluid temperature
            self.T_f_max = 260 + 273.15 # [K] : minimum operating fluid temperature

            self.V_w_max_st = 161*(1000/3600) # [m/s] : maximum wind speed (stowed)
            self.V_w_max_tr = 54*(1000/3600) # [m/s] : maximum wind speed (tracking)

            self.t_lifetime = 30 # [y] : Life expectancy   
               
            # Empirical parameters for correlation in 
            # Semi-empirical correlation to model heat losses 
            # along solar parabolic trough collectors 
            # Rémi Dickes, Vincent Lemort and Sylvain Quoilin 

            self.a = np.array([ 2.062*1e1,
                               -2.893*1e-1,
                                1.472*1e-3,
                                2.240*1e-8,
                                1.198*1e-3,
                                1.403*1e-6,
                                1.045*1e0,
                               -3.043*1e-2,
                               -8.481*1e0,
                                2.073*1e-1])

        else: # manual setting
        
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")

