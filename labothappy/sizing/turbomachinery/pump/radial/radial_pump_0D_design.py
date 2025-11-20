
import numpy as np
import CoolProp.CoolProp as CP
from correlations.turbomachinery.correlations_0D import cordier_line
from toolbox.economics.cpi_data import actualize_price

class RadialPumpODDesign():
    
    def __init__(self, fluid):
        
        self.params = {}
        self.inputs = {}

        self.CAPEX = {}

        self.AS = CP.AbstractState('HEOS', fluid)
        
    #%% DATA HANDLING
    
    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value
               
    #%% CORDIER LINE
    def cordier_line(self):
        # x = specific speed, y = specific diameter
        x = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                      2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30], dtype=float)
    
        y = np.array([9.66, 6.98, 5.65, 4.70, 4.05, 3.57, 3.25,
                      3.01, 1.89, 1.57, 1.40, 1.26, 1.19, 1.11,
                      1.04, 0.97, 0.93, 0.64, 0.48], dtype=float)
    
        def f(xq):
            # np.interp uses linear interpolation and *linear extrapolation*
            return np.interp(xq, x, y)
    
        return f
    
    #%% ESTIMATE EFFICIENCY
    def pump_efficiency_0D_estimation(self, omega_s, Q_pp, *, enforce_limits=False, allow_extrap=False):
        """
        Returns eta(omega_s, Q_pp) using 2D linear interpolation.
        - Vectorized over omega_s and Q_pp with NumPy broadcasting.
        - If enforce_limits=True, points with omega_s outside the interpolated
          [low, high] band (vs Q_pp) are set to NaN.
        - If allow_extrap=True, values outside the data grid are snapped to the edge;
          otherwise they return NaN.
        """
    
        # ---------- Static data ----------
        Q_grid = np.array([18, 36, 72, 180, 360, 720, 1800, 18000], dtype=float)
        omega_grid = np.linspace(10, 90, 17, dtype=float)
    
        eta_18    = np.array([55.62,66.26,71.34,74.12,75.85,76.71,77.09,77.28,76.99,76.71,76.04,75.46,74.70,74.12,73.16,72.40,71.25])
        eta_36    = np.array([59.27,70.77,76.23,78.91,80.64,81.41,81.60,81.60,81.21,80.83,80.06,79.30,78.43,77.67,76.52,75.65,74.50])
        eta_72    = np.array([62.04,73.45,78.91,82.08,83.80,84.66,84.95,84.86,84.19,83.51,82.94,82.27,81.41,80.54,79.68,78.91,78.24])
        eta_180   = np.array([63.96,75.56,81.60,84.86,86.68,87.73,88.02,88.02,87.83,87.06,86.58,85.91,85.14,84.76,83.90,83.13,82.46])
        eta_360   = np.array([65.59,76.71,82.84,86.10,87.92,89.17,89.65,89.65,89.55,89.07,88.59,87.64,86.87,85.81,84.86,84.09,83.23])
        eta_720   = np.array([66.45,77.57,83.51,86.96,89.07,90.03,90.80,90.70,90.70,90.13,89.94,89.17,88.50,87.35,86.68,85.72,84.47])
        eta_1800  = np.array([68.21,78.89,84.79,88.03,90.00,90.94,91.62,91.79,91.62,91.20,90.77,90.17,89.66,89.06,88.29,87.44,86.75])
        eta_18000 = np.array([70.09,80.17,85.98,89.49,91.28,92.56,93.16,93.50,93.50,93.50,93.33,92.99,92.65,92.22,91.62,91.11,90.60])
    
        ETA_GRID = np.vstack([eta_18, eta_36, eta_72, eta_180, eta_360, eta_720, eta_1800, eta_18000]).astype(float)
    
        # omega_s min/max limits vs Q
        LIMITS = np.array([
            [10, 90],   # Q=18
            [10, 90],   # Q=36
            [10, 50],   # Q=72
            [15, 55],   # Q=180
            [15, 70],   # Q=360
            [15, 90],   # Q=720
            [20, 80],   # Q=1800
            [30, 90],   # Q=18000
        ], dtype=float)

        def _interp_limits(q):
            """Interpolate [low, high] omega_s limits for any q (vectorized)."""
            
            low  = np.interp(q, Q_grid, LIMITS[:, 0])
            high = np.interp(q, Q_grid, LIMITS[:, 1])
            return low, high
        
        def _bilinear_interp(omega, q, *, allow_extrap=True):
            """
            Vectorized bilinear interpolation over (Q_grid, omega_grid) -> ETA_GRID.
            Returns NaN outside domain unless allow_extrap=True (then snaps to edge).
            """
            
            omega = np.asarray(omega, dtype=float)
            q = np.asarray(q, dtype=float)
            O, Q = np.broadcast_arrays(omega, q)
        
            # Indices to the "left"
            iq = np.searchsorted(Q_grid, Q, side='right') - 1
            io = np.searchsorted(omega_grid, O, side='right') - 1
        
            # Treat exact-right-edge as inside (snap to last cell)
            iq = np.where(Q == Q_grid[-1], len(Q_grid)-2, iq)
            io = np.where(O == omega_grid[-1], len(omega_grid)-2, io)
        
            # Out-of-bounds mask
            oob = (iq < 0) | (iq >= len(Q_grid)-1) | (io < 0) | (io >= len(omega_grid)-1)
        
            # Clamp for safe indexing (will overwrite oob later)
            iqc = np.clip(iq, 0, len(Q_grid)-2)
            ioc = np.clip(io, 0, len(omega_grid)-2)
        
            q0 = Q_grid[iqc];   q1 = Q_grid[iqc+1]
            o0 = omega_grid[ioc]; o1 = omega_grid[ioc+1]
        
            eta00 = ETA_GRID[iqc,   ioc  ]
            eta01 = ETA_GRID[iqc,   ioc+1]
            eta10 = ETA_GRID[iqc+1, ioc  ]
            eta11 = ETA_GRID[iqc+1, ioc+1]
        
            tq = np.where(q1 > q0, (Q - q0) / (q1 - q0), 0.0)
            to = np.where(o1 > o0, (O - o0) / (o1 - o0), 0.0)
        
            eta0 = eta00*(1 - to) + eta01*to
            eta1 = eta10*(1 - to) + eta11*to
            out = eta0*(1 - tq) + eta1*tq
        
            if allow_extrap:
                # snap out-of-bounds to nearest edge rather than NaN
                # (already clamped indices/tq/to produce an edge value)
                return out
        
            # Otherwise, NaN outside domain
            out = np.where(oob, np.nan, out)
            return out
    
        eta = _bilinear_interp(omega_s, Q_pp, allow_extrap=allow_extrap)
        
        if enforce_limits:
            low, high = _interp_limits(np.asarray(Q_pp, dtype=float))
            # Broadcast to eta's shape
            O, L = np.broadcast_arrays(np.asarray(omega_s, float), low)
            _, H = np.broadcast_arrays(np.asarray(omega_s, float), high)
            in_band = (O >= L) & (O <= H)
            eta = np.where(in_band, eta, np.nan)
    
        return eta

    def cost_estimation(self):
        
        # CO2 1500⋅X0.8 + 50000 $ 2009 X – shaft power (kW)
        
        self.CAPEX['Total'] = actualize_price(self.n_parallel*(50000 + 1500 *(self.W_dot_pp/1000)), 2009, currency="USD") # dollars 2009
        
        #self.n_parallel * 124427*(self.Q_pp)**0.3895 # dollars 2019
                
        return self.CAPEX['Total']

    def Omega_system(self):
        
        import warnings
        warnings.filterwarnings("ignore")
        
        g = 9.81 # m/s
        cordier = cordier_line()

        self.omega_Hz = self.Omega/60 
        self.omega_rads = self.omega_Hz*(2*np.pi)
        
        self.AS.update(CP.PT_INPUTS, self.inputs['P_su'], self.inputs['T_su'])
        self.rho_in = self.AS.rhomass()

        # Height difference
        self.DH = (self.inputs['P_ex'] - self.inputs['P_su'])/(self.rho_in*g) + (self.inputs['H2'] - self.inputs['H1']) + (self.inputs['v2']**2 - self.inputs['v1']**2)/(2*g)
        self.Q = self.inputs['m_dot']/self.rho_in 
        
        self.Q_pp = self.Q/self.n_parallel
        
        # Specific speed
        self.Omega_s = (self.omega_rads*np.sqrt(self.Q_pp))/((g*self.DH)**(3/4))

        self.Omega_s_imp = 2733*self.Omega_s

        if self.Omega_s_imp < 1500:
            self.geom = "Radial-Vane Area"
        elif self.Omega_s_imp < 4000:
            self.geom = "Francis-Vane Area"
        elif self.Omega_s_imp < 7500:
            self.geom = "Mixed-Flow Area"
        else:
            self.geom = "Axial-Flow Area"
        
        # Specific Diameter using Cordier Line
        D_s = cordier(self.Omega_s)

        # Real diameter 
        self.D = D_s*((self.Q_pp)**(1/2))/(g*self.DH)**(1/4) # [m]
        self.eta_is = self.pump_efficiency_0D_estimation(self.Omega_s_imp, self.Q_pp * 3600)/100
        
        h_in = self.AS.hmass()
        
        self.AS.update(CP.PSmass_INPUTS, self.inputs['P_ex'], self.AS.smass())
        h_is = self.AS.hmass()
        self.h_ex = h_in + (h_is-h_in)/self.eta_is
        
        self.W_dot_pp = self.inputs['m_dot']*(self.h_ex - h_is)/self.n_parallel
        self.W_dot = self.inputs['m_dot']*(self.h_ex - h_is)
    
    def pick_npp_by_threshold(self, eta, pp_threshold=2.0):
        """
        Step 1: Identify the column (speed) with the highest efficiency (ignoring NaNs).
        Step 2: In that column only, choose the smallest n_pp such that
                adding another pump yields < pp_threshold %-points gain.
    
        Returns:
            best_col : index of selected column (speed)
            best_row : index of selected row (n_pp)
        """
        eta = np.asarray(eta)
    
        # --- Step 1: Find column with highest efficiency (ignore NaNs)
        max_per_col = np.nanmax(eta, axis=0)       # max in each column
        best_col = int(np.nanargmax(max_per_col))  # column with highest max
    
        # --- Step 2: Apply the marginal-gain rule on that chosen column
        col = eta[:, best_col]
        valid_idx = np.where(~np.isnan(col))[0]
    
        if len(valid_idx) == 0:
            return best_col, None  # no feasible design
    
        best_row = valid_idx[-1]  # default: highest valid n_pp
    
        for i in range(len(valid_idx) - 1):
            r = valid_idx[i]
            r_next = valid_idx[i+1]
            gain_pp = col[r_next] - col[r]
    
            if np.isnan(gain_pp):
                continue
    
            if gain_pp < pp_threshold:
                best_row = r
                break
    
        return best_col, best_row
    
    def design(self):
        
        self.eta_matrix = np.zeros([len(self.params["n_parallel_choices"]), len(self.params['Omega_choices'])])
        self.cost_matrix = np.zeros([len(self.params["n_parallel_choices"]), len(self.params['Omega_choices'])])
        self.power_matrix = np.zeros([len(self.params["n_parallel_choices"]), len(self.params['Omega_choices'])])
        
        for i in range(len(self.params["n_parallel_choices"])):
            self.n_parallel = self.params["n_parallel_choices"][i]
            for j in range(len(self.params['Omega_choices'])):
                self.Omega = self.params['Omega_choices'][j]
                self.Omega_system()
                self.eta_matrix[i][j] = self.eta_is
        
        index_omega, index_npp = self.pick_npp_by_threshold(self.eta_matrix)

        self.Omega = self.params['Omega_choices'][index_omega]
        self.n_parallel = self.params["n_parallel_choices"][index_npp]

        self.Omega_system()
        self.cost_estimation()
        
        return 

if __name__ == "__main__":

    PP_des = RadialPumpODDesign('CO2')
    
    PP_des.set_inputs(
        P_su = 40e5, # Pa
        P_ex = 160*1e5, # Pa
        
        T_su = 5.15 + 273.15, # T
        
        H1 = 0, # m
        H2 = 0, # m
        
        v1 = 0, # m/s
        v2 = 0, # m/s
        
        m_dot = 500, # kg/s
        )
    
    PP_des.set_parameters(
        Omega_choices = np.array([750, 1000, 1500, 3000]),
        n_parallel_choices = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        )
    
    PP_des.design()
    
    # val = PP_des.pump_efficiency_0D_estimation(omega_s=55.0, Q_pp=500.0)
    
    # # Vectors (broadcasted):
    # omegas = np.array([20, 40, 60, 80])
    # flows  = np.array([30, 200, 1000, 5000])
    # vals = PP_des.pump_efficiency_0D_estimation(omegas, flows)
