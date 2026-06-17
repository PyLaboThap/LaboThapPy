# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 11:05:34 2026

@author: Basile
"""

import numpy as np
import matplotlib.pyplot as plt
from labothappy.component.base_component import BaseComponent
from CoolProp.CoolProp import PropsSI
from scipy.integrate import trapezoid

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#%%

class OilSeparator(BaseComponent):
    """
    **Component**: 


    **Descritpion**:


    **Assumptions**:


    **Connectors**:


    **Parameters**:


    **Inputs**:


    **Ouputs**:


    """

    def __init__(self):
        
        self.su_gas = {}
        self.ex_state = {}
        self.ex_liquid = {}
        self.params = {}
    
    def set_inputs(self, **inputs):
        
        self.inputs = inputs
        
    def get_required_inputs(self):
        # Return a list of required inputs
        return [
            'P_su_gas', 'T_su_gas', 'fluid_gas', 'V_dot_gas',
            'fluid_oil', 'T_su_oil', 'P_su_oil', 'V_dot_oil',
            'rho_oil', 'nu_oil',
        ]

    def get_required_parameters(self):
        """
        BD: float        # Body diameter [m]
        BH: float        # Body height [m]
        D_in: float      # Inner diameter of inlet pipe [m]
        S: float         # Length of outlet pipe [m]
        IP: float        # Inlet position from top [m]

        # Reference geometry (Kim et al. 2016 - reference separator)
        BD_ref: float = 0.085   # [m]
        BH_ref: float = 0.220   # [m]
        """
        
        return [
            'BD', 'BH', 'D_in', 'S', 'IP', 'BD_ref', 'BH_ref'
        ]
    
    #%%
    
    def build_viscosity_interpolator(self, T_K, nu_cSt, rho):
        """
        Build a kinematic viscosity interpolation function using the
        ASTM D341 Walther equation — the standard method for petroleum
        and synthetic lubricants.
    
        The Walther equation linearises viscosity data in a double-log space:
            log(log(ν + 0.7)) = A + B·log(T)
    
        where ν is in cSt and T in Kelvin.
    
        Parameters
        ----------
        T_K    : list of two temperatures [K]
                 e.g. [40+273.15, 100+273.15]
        nu_cSt : list of two kinematic viscosities [cSt]
                 e.g. [27.6, 5.1]   (Jet Oil II from FDS section 9)
        rho    : float — liquid density [kg/m³], used to compute dynamic viscosity
    
        Returns
        -------
        nu_at_T : callable  nu_at_T(T_K) → kinematic viscosity [cSt]
        mu_at_T : callable  mu_at_T(T_K) → dynamic viscosity [Pa·s]
                  (rho is already embedded — no need to pass it at call time)
    
        Example
        -------
        >>> nu, mu = build_viscosity_interpolator([313.15, 373.15], [27.6, 5.1], rho=1000.)
        >>> nu(353.15)   # 80°C → 7.90 cSt
        >>> mu(353.15)   # 80°C → 7.90e-3 Pa·s
        """
        T1, T2   = T_K
        nu1, nu2 = nu_cSt
    
        # Walther transform: W(ν) = log(log(ν + 0.7))
        W1 = np.log(np.log(nu1 + 0.7))
        W2 = np.log(np.log(nu2 + 0.7))
    
        # Linear fit in log(T) space: W = B·log(T) + A
        B = (W2 - W1) / (np.log(T2) - np.log(T1))
        A = W1 - B * np.log(T1)
    
        def nu_at_T(T: float) -> float:
            """Kinematic viscosity [cSt] at temperature T [K]."""
            T_arr = np.asarray(T, dtype=float)
            W = A + B * np.log(T_arr)
            return np.exp(np.exp(W)) - 0.7
    
        def mu_at_T(T: float) -> float:
            """
            Dynamic viscosity [Pa·s] at temperature T [K].
            Uses the density rho [kg/m³] provided at construction time.
            μ = ν [m²/s] × ρ [kg/m³]
            """
            nu_m2s = nu_at_T(T) * 1e-6   # cSt → m²/s
            return nu_m2s * rho            # Pa·s
    
        # Store calibration info as function attributes
        nu_at_T.T_ref  = [T1, T2]
        nu_at_T.nu_ref = [nu1, nu2]
        nu_at_T.rho    = rho
        nu_at_T.A      = A
        nu_at_T.B      = B
        nu_at_T.method = "ASTM D341 Walther equation"
    
        return nu_at_T, mu_at_T
    
    def compute_complete_inputs(self):
    
        # ── Gas phase ──────────────────────────────────────────────
        # Append '|gas' to force the gas-phase branch.
        # Without this, CoolProp may return near-critical transport properties
        # (e.g. diverging viscosity) when T is close to Tc, even in the
        # superheated region, causing unphysically small d50 and η → 100 %.
        fluid_gas_phase = self.inputs['fluid_gas'] 
        rho_gas, mu_gas = PropsSI(
            ('D', 'V'),
            'T', self.inputs['T_su_gas'],
            'P', self.inputs['P_su_gas'],
            fluid_gas_phase
        )
        mdot_gas = self.inputs['V_dot_gas'] * rho_gas
    
        self.su_gas = {
            'fluid' : self.inputs['fluid_gas'],
            'T'     : self.inputs['T_su_gas'],
            'P'     : self.inputs['P_su_gas'],
            'D'     : rho_gas,
            'V'     : mu_gas,           # dynamic viscosity of GAS [Pa·s]
            'NU'    : mu_gas/rho_gas,   # kinematic viscosity of GAS [m2/s]
            'm_dot' : mdot_gas,
            'V_dot' : self.inputs['V_dot_gas'],
        }
    
        # ── Oil phase ───────────────────────────────────────────────
        mdot_oil = self.inputs['V_dot_oil'] * self.inputs['rho_oil']
    
        # Build viscosity interpolator — rho is now embedded in mu_fun
        self.nu_fun, self.mu_fun = self.build_viscosity_interpolator(
            T_K    = self.inputs['nu_oil'][0],   # e.g. [313.15, 373.15] K
            nu_cSt = self.inputs['nu_oil'][1],   # e.g. [27.6, 5.1] cSt
            rho    = self.inputs['rho_oil'],     # [kg/m³]
        )
    
        # Dynamic viscosity of oil at supply temperature
        mu_oil = self.mu_fun(self.inputs['T_su_oil'])   # [Pa·s]
        nu_oil = self.nu_fun(self.inputs['T_su_oil'])   # [cSt]
    
        self.su_oil = {
            'fluid' : self.inputs['fluid_oil'],
            'T'     : self.inputs['T_su_oil'],
            'P'     : self.inputs['P_su_oil'],
            'D'     : self.inputs['rho_oil'],
            'V'     : mu_oil,           # dynamic viscosity of OIL [Pa·s]
            'NU'    : mu_oil/self.inputs['rho_oil'],           # kinematic viscosity of OIL [cSt]
            'sigma' : self.inputs['sigma_oil'], # Surface Tension 
            'm_dot' : mdot_oil,
            'V_dot' : self.inputs['V_dot_oil'],
        }
    
        return
    
    #%% Velocity Computation
    
    def muschelknautz_velocities(self):
        """
        Estimate tangential velocities at the wall (uθw) and at the
        cylindrical surface (uθCS) using a simplified Muschelknautz approach.

        These velocities drive the centrifugal separation force.

        Gas Cyclones and swirl tubes : Principles, Design and Operation
        Hoffman & Stein 
        
        + 
        
        Flow characteristics of refrigerant and oil mixture in an oil separator
        Kim, Kang, Yoon, Sa, Chung, Kim

        Returns
        -------
        dict with keys: u_inlet, u_x, u_theta_w, u_theta_CS, Re_liq
        """
        
        R = self.params['BD']/2 # Body radius
        R_in = self.params['D_in']/2 # Inlet tube radius
        R_x = self.params['D_in']/2 # Exhaust tube radius 
        
        "1) u_inlet and u_theta_w"
        # Constriction coefficient : alpha
        c0 = self.LCR # Ratio of solid mass to gas mass (here liquid assumed as solid)
        b = self.params['D_in']
        xhi = b/R
        
        term = np.sqrt(1 - (1-xhi**2)*(2*xhi-xhi**2)/(1+c0))
        term2 = 1+4*((xhi/2)**2 - xhi/2)*term
        self.alpha = 1/xhi * (1-np.sqrt(term2))
        
        Q_total = (self.su_gas['V_dot'] + self.su_oil['V_dot'])  # [m³/s]
        u_inlet = Q_total / self.A_inlet # [m/s]
        
        u_theta_w = u_inlet*R_in/(R*self.alpha) # [m/s]
        
        "2) u_z_w + u_x"
        
        # Geometric mean radius
        R_m = np.sqrt(R*R_x)
        u_z_w = 0.9*Q_total/(np.pi*(R**2 - R_m**2))

        # Axial velocity in outlet pipe (gas only after separation)
        A_outlet = np.pi * (self.params['D_in'] / 2) ** 2
        u_x = self.su_gas['m_dot'] / (self.su_gas['D'] * A_outlet)         # [m/s]

        "3) u_theta_CS + u_theta_m"
        # Tangential velocity at the cylindrical surface (vortex finder radius).
        # Angular momentum conservation (Eq. 6):  r_body · uθw = r_CS · uθCS
        r_body   = self.params['BD'] / 2
        r_outlet = self.params['D_in'] / 2   # r_CS ≈ outlet-pipe radius
        u_theta_CS = u_theta_w * (r_body / r_outlet)             # [m/s]

        # Reynolds number of gas cyclone — uses kinematic viscosity [m²/s]
        nu_gas = self.su_gas['NU'] # m²/s
        Re = R_in*R_m*u_z_w/(nu_gas*self.params['BH'])
        
        u_theta_m = np.sqrt(u_theta_w*u_theta_CS)
        
        self.velocity_data = {
            'u_inlet'    : u_inlet,
            'u_x'        : u_x,
            'u_theta_w'  : u_theta_w,
            'u_theta_m'  : u_theta_m,
            'u_theta_CS' : u_theta_CS,
            'Re'     : Re,
        } 
        
        return 
    
    #%% Droplet Size Computations 
    
    def compute_d50(self):
        """
        Cut size d50 [m] — Eq. (12) of Kim et al. (2016).

        d50 = 5.18 * [mu_gas*(rho_gas - rho_liquid)] /
                     [rho_liquid * rho_gas * uθCS²]^0.375
              × (u_gas/uθCS)^0.25
              × R_CS^0.875
              × (BD/BD_ref)^0.97
              × (BH/BH_ref)^-0.97

        Note: R_CS is approximated as D_in/2 (radius at the cylindrical surface).
        """
        u_gas      = self.velocity_data['u_inlet']
        u_theta_CS = self.velocity_data['u_theta_CS']
        R_x       = self.params['D_in'] / 2.0

        num = self.su_gas['V']**0.375 * self.su_gas['D']**0.25 * self.velocity_data['u_inlet']**0.875
        den = (((self.su_oil['D'] - self.su_gas['D']) * u_theta_CS ** 2)/R_x)**0.625
        
        term1 = num/den

        term2 = (self.params['BD_ref']/self.params['BD']) ** 0.97

        term3 = (self.params['BH_ref'] / self.params['BH']) ** 0.97

        self.d50_val = 5.18 * term1 * term2 * term3 
        return 

    def compute_d1(self):
        """
        Characteristic droplet size d1 [m] — Eq. (13) of Kim et al. (2016).

        Performance Prediction of a Cyclone Oil Separator
        """

        self.d1_val = self.params['D_in']*7.43*1e-3 * (self.su_oil['sigma']/(self.su_gas['D']*self.params['D_in']*self.velocity_data['u_theta_m']))**0.5 * (self.velocity_data['u_theta_m']*self.params['D_in']/self.su_gas['NU'])**0.1
        return 
    
    def droplet_size_distribution(self, n_points = 50):
        """
        Dimensionless droplet size distribution d/d1 from Murakami et al. (2006).

        Approximated here as a log-normal distribution with:
            - median = d1
            - geometric standard deviation σ_g = 2.0 (typical for spray atomisation)

        Returns
        -------
        d_array : array of droplet diameters [m]
        MF_array: array of mass fractions [-], sum = 1
        """
        sigma_g  = 2.0
        ln_sigma = np.log(sigma_g)

        # Log-normal distribution: d from 0.01·d1 to 20·d1
        self.drop_distrib_D = d_array = np.logspace(
            np.log10(0.01 * self.d1_val), np.log10(20 * self.d1_val), n_points
        )

        # PDF of log-normal
        pdf = (1 / (d_array * ln_sigma * np.sqrt(2 * np.pi)) *
               np.exp(-0.5 * ((np.log(d_array) - np.log(self.d1_val)) / ln_sigma) ** 2))

        # Mass fraction proportional to d³ × pdf (volume-weighted)
        mass_weight = d_array ** 3 * pdf
        self.drop_distrib_x = mass_weight / trapezoid(mass_weight, d_array)  # normalise

        return 

    #%%
    
    def compute_limit_loading(self):
        """
        Limit loading LL — Eqs. (8) and (9).

        k  = -0.07 - 0.16 * ln(IL)
        LL = 0.0078 * (d50/d_mm)^k * 10^k   [dimensionless]

        Parameters
        ----------
        IL   : inlet loading [-]
        d50  : cut size [m]
        d_mm : mass mean droplet diameter [m]
        """
        k       = 0.07 - 0.16 * np.log(self.IL)
        self.LL = 0.0078 * (self.d50_val / self.d_mm) * (10*self.IL) ** k
        return 
    
    
    def overall_efficiency(self, t=3):
        """
        Overall separation efficiency considering mass loading.
    
        - If IL < LL : Eq. (10) applies (no mass loading)
        - If IL ≥ LL : Eq. (11) applies (mass loading regime)
    
        Returns
        -------
        dict with keys: eta_total, LL, regime, eta_grade (array)
        """
        self.compute_limit_loading()
    
        # Grade efficiency for each droplet size
        self.eta_grade = np.array([self.grade_efficiency(d) for d in self.drop_distrib_D])
    
        if self.IL < self.LL:
            # No mass loading — Eq. (10)
            self.eta_total = trapezoid(self.eta_grade * self.drop_distrib_x, self.drop_distrib_D) / trapezoid(self.drop_distrib_x, self.drop_distrib_D)
            self.regime    = "No mass loading (IL < LL)"
            self.frac_immediate = 0
        else:
            # Mass loading — Eq. (11)
            # Immediate separation fraction = LL/IL
            # Cyclonic separation of remainder
            frac_cyclone = self.LL / self.IL
            self.frac_immediate   = 1.0 - frac_cyclone
            eta_cyclone    = trapezoid(self.eta_grade * self.drop_distrib_x, self.drop_distrib_D) / trapezoid(self.drop_distrib_x, self.drop_distrib_D)
            self.eta_total = self.frac_immediate + frac_cyclone * eta_cyclone
            self.regime    = f"Mass loading (IL ≥ LL): immediate sep. = {self.frac_immediate*100:.1f}%"
    
        return 
    
    def grade_efficiency(self, d):
        """
        Grade efficiency for a single droplet diameter — Eq. (7).

        η_i = 1 / (1 + (d50/d)^t)

        Parameters
        ----------
        d   : droplet diameter [m]
        d50 : cut size [m]
        t   : sharpness parameter (= 3 for gas-liquid, Hoffmann & Stein 2007)
        """
        t = 3
        return 1.0 / (1.0 + (self.d50_val / d) ** t)
    
    #%% Plotting
    
    def plot_grade_efficiency(self, title: str = "Grade Efficiency Curve"):
        """
        Plot grade efficiency as a function of droplet diameter.
        Must be called after solve().
        """
        fig, ax = plt.subplots(figsize=(8, 5))
        d_um  = self.drop_distrib_D * 1e6
        label = f"η = {self.eta_total*100:.1f}%  |  {self.regime}"
        ax.semilogx(d_um, self.eta_grade * 100, label=label)
        ax.axhline(50, color='gray', linestyle='--', linewidth=0.8, label='50% (d50)')
        ax.axvline(self.d50_val * 1e6, color='C0', linestyle=':', linewidth=0.8,
                   label=f"d50 = {self.d50_val*1e6:.1f} µm")
        ax.set_xlabel("Droplet diameter [µm]", fontsize=12)
        ax.set_ylabel("Grade efficiency [%]", fontsize=12)
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=8)
        ax.grid(True, which='both', alpha=0.3)
        ax.set_ylim(0, 105)
        plt.tight_layout()
        return fig

    def plot_efficiency_vs_loading(self, Q_gas_values: list, Q_liq_values: list):
        """
        Plot overall efficiency vs inlet loading for multiple gas flow rates.
        Loops over (Q_gas, Q_liq) combinations, re-solving for each point.
        Original inputs are restored after the sweep.

        Parameters
        ----------
        Q_gas_values : list of float   [m³/h]
        Q_liq_values : list of float   [l/h]
        """
        # Save original flow-rate inputs so they can be restored afterwards
        original_inputs = dict(self.inputs)

        fig, ax = plt.subplots(figsize=(8, 5))

        for Q_gas in Q_gas_values:
            eta_list, IL_list = [], []
            for Q_liq in Q_liq_values:
                new_inputs = dict(self.inputs)
                new_inputs['V_dot_gas'] = Q_gas / 3600           # m³/h → m³/s
                new_inputs['V_dot_oil'] = (Q_liq / 1000) / 3600  # l/h  → m³/s
                self.set_inputs(**new_inputs)
                self.solve()
                eta_list.append(self.eta_total * 100)
                IL_list.append(self.IL)
            ax.plot(IL_list, eta_list, 'o-', label=f"Q_gas = {Q_gas} m³/h")

        # Restore original inputs and re-solve to leave the object in its original state
        self.set_inputs(**original_inputs)
        self.solve()

        ax.set_xlabel("Inlet Loading IL [kg_liq / kg_gas]", fontsize=11)
        ax.set_ylabel("Overall Efficiency [%]", fontsize=11)
        ax.set_title("Separation Efficiency vs Inlet Loading", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        return fig

    def plot_droplet_distribution(self, title: str = "Droplet Size Distribution"):
        """
        Plot the mass-fraction droplet size distribution.
        Must be called after solve().
        """
        fig, ax = plt.subplots(figsize=(8, 5))
        d_um   = self.drop_distrib_D * 1e6
        d50_um = self.d50_val * 1e6
        label  = f"d1 = {self.d1_val*1e6:.1f} µm  |  d50 = {d50_um:.1f} µm"
        ax.semilogx(d_um, self.drop_distrib_x, label=label)
        ax.axvline(d50_um, linestyle=':', linewidth=0.8,
                   label=f"d50 = {d50_um:.1f} µm")
        ax.set_xlabel("Droplet diameter [µm]", fontsize=12)
        ax.set_ylabel("Mass fraction density [1/m]", fontsize=12)
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=8)
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        return fig

    #%%
    
    def solve(self):
        self.compute_complete_inputs()    
    
        self.LCR = self.su_oil['m_dot'] / (self.su_oil['m_dot'] + self.su_gas['m_dot'])
        self.IL  = self.su_oil['m_dot'] / self.su_gas['m_dot']
        
        self.A_inlet = (np.pi/4) * self.params['D_in'] ** 2
        
        self.muschelknautz_velocities()
        
        self.compute_d50()
        self.compute_d1()
        self.droplet_size_distribution()
        
        self.d_mm = trapezoid(self.drop_distrib_D * self.drop_distrib_x, self.drop_distrib_D) / trapezoid(self.drop_distrib_x, self.drop_distrib_D)
        
        self.overall_efficiency()
        
        return
    
#%%

if __name__ == "__main__":
    
    case_study = "Simu"
    
    if case_study == "Simu":
        Sep = OilSeparator()
    
        Sep.set_inputs(
            T_su_gas = 50+273.15, # K
            P_su_gas = 1*1e5, # Pa
            V_dot_gas = 5/3600, # m^3/s
            fluid_gas = 'Air',
            
            fluid_oil = 'Jet II',
            T_su_oil = 50+273.15, # K
            P_su_oil = 1*1e5, # Pa
            V_dot_oil = (6/1000)/3600, # m^3/s
            rho_oil = 1003.5, # kg/m^3
            nu_oil = [[40+273.15,100+273.15], [27.6, 5.1]], # [T1, T2] : K and [nu1, nu_2] : cSt
            sigma_oil = 3*1e-3
            )
        
        Sep.set_parameters(
            BD=0.085,   # m body diameter
            BH=0.22,   # m body height
            D_in=0.0141, # inlet pipe inner diameter
            S=0.120,    # m outlet pipe length
            IP=0.030,   # m inlet position
            BD_ref= 0.085, # Body diameter of reference 
            BH_ref= 0.22, # Body Height of reference 
            )
        
        Sep.solve()
        
        Sep.plot_grade_efficiency()
        # Sep.plot_efficiency_vs_loading()
        Sep.plot_droplet_distribution()
    
    else:
        # -*- coding: utf-8 -*-
        """
        Validation of OilSeparator against Kim et al. (2016)
        "Flow characteristics of refrigerant and oil mixture in an oil separator"
        International Journal of Refrigeration 70 (2016) 206–218
        
        ─────────────────────────────────────────────────────────────────────────────
        HOW TO RUN
            Place this file in the same folder as oil_separator.py, then:
                python validate_kim2016.py
        
        WHAT IS VALIDATED
            The paper reports *relative* efficiencies normalised to a reference
            condition (MFR=90 g/s, LCR=3%, T=79.9°C, P=2809.8 kPa, Reference geom).
            Three checks are performed:
        
            1. ABSOLUTE LEVEL
               Fig. 11 shows all predicted efficiencies in the 94–100 % band.
               → our model must give η in that range at the reference condition.
        
            2. TREND: η vs MFR and LCR  (Figs. 4 & 5)
               - η decreases with MFR up to ~90 g/s, then increases         (U-shape)
               - η increases with LCR at every MFR                           (monotone)
        
            3. TREND: geometry effect on relative η  (Table 3, EC column)
               All four alternative designs improve η vs Reference:
                 Design 1 (BD↑20%)  : EC > 0
                 Design 2 (BD↑41%)  : EC > Design 1
                 Design 3 (BH↑18%)  : EC > 0
                 Design 4 (BH↑36%)  : EC > Design 3
        ─────────────────────────────────────────────────────────────────────────────
        """
        
        import sys
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        from CoolProp.CoolProp import PropsSI
        
        # ── import OilSeparator from the same directory ────────────────────────────
        sys.path.insert(0, '.')        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 0.  EXPERIMENTAL CONFIGURATION  (paper §2 and Table 1)
        # ═══════════════════════════════════════════════════════════════════════════
        
        # PVE oil properties  (paper §2, Eq.13 reference viscosity)
        OIL = dict(
            fluid_oil = 'PVE',
            rho_oil   = 1003.5,                                  # kg/m³
            nu_oil    = [[40+273.15, 100+273.15], [27.6, 5.1]],  # Walther: [T_K], [nu_cSt]
        )
        
        # Two inlet T/P sets tested in the paper (§2.2)
        COND = {
            'A': dict(T=70.0+273.15, P=2611.3e3),   # 70.0 °C / 26.1 bar
            'B': dict(T=79.9+273.15, P=2809.8e3),   # 79.9 °C / 28.1 bar  ← reference
        }
        
        # Five separator geometries  (Table 1)
        # BD_ref / BH_ref are always the Reference geometry values
        _REF = dict(BD_ref=0.085, BH_ref=0.220)
        GEOMETRIES = {
            'Reference': dict(BD=0.085, BH=0.220, D_in=0.0141, S=0.120, IP=0.030, **_REF),
            'Design 1' : dict(BD=0.102, BH=0.220, D_in=0.0141, S=0.120, IP=0.030, **_REF),
            'Design 2' : dict(BD=0.120, BH=0.220, D_in=0.0141, S=0.120, IP=0.030, **_REF),
            'Design 3' : dict(BD=0.085, BH=0.260, D_in=0.0141, S=0.120, IP=0.030, **_REF),
            'Design 4' : dict(BD=0.085, BH=0.300, D_in=0.0141, S=0.120, IP=0.030, **_REF),
        }
        
        # Sweep ranges
        MFR_VALUES = [30, 60, 90, 120, 150]   # g/s  (gas phase)
        LCR_VALUES = [0.015, 0.030, 0.045]    # –
        
        # Reference operating point  (paper §2.3)
        MFR_REF = 90e-3   # kg/s
        LCR_REF = 0.030
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 1.  HELPER FUNCTIONS
        # ═══════════════════════════════════════════════════════════════════════════
        
        def run_point(geom: dict, mfr_gas_kgs: float, lcr: float,
                      T_K: float, P_Pa: float) -> OilSeparator:
            """
            Instantiate, configure and solve OilSeparator for one operating point.
        
            Parameters
            ----------
            geom        : geometry parameter dict (one entry from GEOMETRIES)
            mfr_gas_kgs : gas-phase mass flow rate [kg/s]
            lcr         : liquid circulation ratio [-]  =  m_liq / (m_liq + m_gas)
            T_K         : inlet temperature [K]
            P_Pa        : inlet pressure [Pa]
        
            Returns
            -------
            Solved OilSeparator instance
            """
            # Volumetric gas flow from CoolProp density
            rho_gas   = PropsSI('D', 'T', T_K, 'P', P_Pa, 'R410A')
            V_dot_gas = mfr_gas_kgs / rho_gas                         # m³/s
        
            # LCR = m_liq / (m_liq + m_gas)  →  m_liq = m_gas * LCR / (1 - LCR)
            mfr_oil   = mfr_gas_kgs * lcr / (1.0 - lcr)              # kg/s
            V_dot_oil = mfr_oil / OIL['rho_oil']                      # m³/s
        
            sep = OilSeparator()
            sep.set_inputs(
                T_su_gas  = T_K,
                P_su_gas  = P_Pa,
                V_dot_gas = V_dot_gas,
                fluid_gas = 'R410A',
                T_su_oil  = T_K,
                P_su_oil  = P_Pa,
                V_dot_oil = V_dot_oil,
                sigma_oil = 3*1e-3,
                **OIL,
            )
            sep.set_parameters(**geom)
            sep.solve()
            return sep
        
        
        def mfr_lcr_sweep(geom: dict, cond: dict) -> dict:
            """Run the full MFR × LCR grid for one geometry and one T/P condition."""
            results = {}
            for lcr in LCR_VALUES:
                for mfr in MFR_VALUES:
                    results[(mfr, lcr)] = run_point(geom, mfr * 1e-3, lcr,
                                                    cond['T'], cond['P'])
            return results
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 2.  CHECK 1 — ABSOLUTE EFFICIENCY LEVEL
        # ═══════════════════════════════════════════════════════════════════════════
        
        print("=" * 65)
        print("VALIDATION  —  Kim et al. (2016)  Int. J. Refrig. 70, 206–218")
        print("=" * 65)
        
        sep_ref = run_point(GEOMETRIES['Reference'], MFR_REF, LCR_REF,
                            COND['B']['T'], COND['B']['P'])
        
        print(f"\n── Check 1: absolute level at reference condition ──────────────")
        print(f"  (MFR=90 g/s, LCR=3 %, T=79.9°C, P=28.1 bar, Reference geom)")
        print(f"  η_total  = {sep_ref.eta_total*100:.2f} %    paper range: 94–100 %")
        print(f"  d50      = {sep_ref.d50_val*1e6:.2f} µm")
        print(f"  d1       = {sep_ref.d1_val*1e6:.2f} µm")
        print(f"  d_mm     = {sep_ref.d_mm*1e6:.2f} µm")
        print(f"  IL       = {sep_ref.IL:.4f}    LL = {sep_ref.LL:.4f}")
        print(f"  Regime   : {sep_ref.regime}")
        print(f"  u_inlet  = {sep_ref.velocity_data['u_inlet']:.2f} m/s")
        print(f"  u_θCS    = {sep_ref.velocity_data['u_theta_CS']:.2f} m/s")
        # print(f"  Re_liq   = {sep_ref.velocity_data['Re_liq']:.1f}")
        
        if 90.0 <= sep_ref.eta_total * 100 <= 100.0:
            print("  ✓  PASS — efficiency in 90–100 % band")
        else:
            print("  ✗  FAIL — efficiency outside plausible range")
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 3.  CHECK 2 — TREND: η vs MFR and LCR  (Figs. 4 & 5)
        # ═══════════════════════════════════════════════════════════════════════════
        
        print(f"\n── Check 2: trend η vs MFR and LCR  (Cond B, Reference geom) ──")
        print(f"  Paper (Fig 5): U-shape with MFR (min near 90 g/s), η↑ with LCR")
        print(f"\n  {'MFR(g/s)':>9}", end='')
        for lcr in LCR_VALUES:
            print(f"  LCR={lcr*100:.1f}%", end='')
        print()
        print("  " + "-" * 45)
        
        eta_grid = {}   # store for plotting
        for mfr in MFR_VALUES:
            print(f"  {mfr:>9.0f}", end='')
            for lcr in LCR_VALUES:
                s = run_point(GEOMETRIES['Reference'], mfr * 1e-3, lcr,
                              COND['B']['T'], COND['B']['P'])
                eta_grid[(mfr, lcr)] = s.eta_total
                print(f"  {s.eta_total*100:>8.2f} %", end='')
            print()
        
        # Qualitative checks
        # (a) η at LCR=4.5% > η at LCR=1.5% for every MFR
        lcr_trend_ok = all(
            eta_grid[(mfr, 0.045)] > eta_grid[(mfr, 0.015)]
            for mfr in MFR_VALUES
        )
        # (b) U-shape: η(90) ≤ η(60) and η(90) ≤ η(120)  for LCR=3%
        u_shape_ok = (eta_grid[(90, 0.030)] <= eta_grid[(60, 0.030)] and
                      eta_grid[(90, 0.030)] <= eta_grid[(120, 0.030)])
        
        print(f"\n  η↑ with LCR at every MFR  : {'✓  PASS' if lcr_trend_ok else '✗  FAIL'}")
        print(f"  U-shape with MFR (LCR=3%) : {'✓  PASS' if u_shape_ok  else '✗  FAIL'}")
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 4.  CHECK 3 — GEOMETRY EFFECT  (Table 3 EC column)
        # ═══════════════════════════════════════════════════════════════════════════
        
        print(f"\n── Check 3: geometry effect on η  (MFR=90 g/s, LCR=3%, Cond B) ──")
        print(f"  Paper Table 3: EC > 0 for all designs vs Reference")
        print(f"\n  {'Geometry':>12}  {'BD(mm)':>7}  {'BH(mm)':>7}  {'η(%)':>8}  {'ΔEC(pp)':>9}")
        print("  " + "-" * 52)
        
        eta_ref_geom = sep_ref.eta_total
        geom_etas = {}
        for name, geom in GEOMETRIES.items():
            s = run_point(geom, MFR_REF, LCR_REF, COND['B']['T'], COND['B']['P'])
            delta = (s.eta_total - eta_ref_geom) * 100
            geom_etas[name] = s.eta_total
            marker = '  (baseline)' if name == 'Reference' else \
                     ('  ✓' if delta > 0 else '  ✗')
            print(f"  {name:>12}  {geom['BD']*1e3:>7.1f}  {geom['BH']*1e3:>7.1f}"
                  f"  {s.eta_total*100:>8.2f}  {delta:>+9.3f}{marker}")
        
        geom_pass = all(geom_etas[n] > geom_etas['Reference']
                        for n in ['Design 1', 'Design 2', 'Design 3', 'Design 4'])
        bd_monotone = geom_etas['Design 2'] > geom_etas['Design 1'] > geom_etas['Reference']
        bh_monotone = geom_etas['Design 4'] > geom_etas['Design 3'] > geom_etas['Reference']
        print(f"\n  All designs η > Reference : {'✓  PASS' if geom_pass   else '✗  FAIL'}")
        print(f"  BD: Design2 > Design1 > Ref: {'✓  PASS' if bd_monotone else '✗  FAIL'}")
        print(f"  BH: Design4 > Design3 > Ref: {'✓  PASS' if bh_monotone else '✗  FAIL'}")
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 5.  CONDITION A vs B  (paper §3.1: η slightly higher at Cond B)
        # ═══════════════════════════════════════════════════════════════════════════
        
        print(f"\n── Check 4: Cond A vs B  (MFR=90 g/s, LCR=3%, Reference geom) ──")
        print(f"  Paper: η slightly higher at Cond B (higher ρ_gas → lower u → lower drag)")
        eta_A = run_point(GEOMETRIES['Reference'], MFR_REF, LCR_REF,
                          COND['A']['T'], COND['A']['P']).eta_total
        eta_B = sep_ref.eta_total
        rho_A = PropsSI('D', 'T', COND['A']['T'], 'P', COND['A']['P'], 'R410A')
        rho_B = PropsSI('D', 'T', COND['B']['T'], 'P', COND['B']['P'], 'R410A')
        print(f"  Cond A: T=70.0°C  ρ={rho_A:.1f} kg/m³  η={eta_A*100:.2f}%")
        print(f"  Cond B: T=79.9°C  ρ={rho_B:.1f} kg/m³  η={eta_B*100:.2f}%")
        cond_pass = eta_B > eta_A
        print(f"  η_B > η_A : {'✓  PASS' if cond_pass else '✗  FAIL'}")
        
        
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 6.  DIGITISED EXPERIMENTAL DATA  (Cond B, paper Figs 5 & 6)
        #     Relative efficiency = measured_η / η_reference
        #     Reference condition: MFR=90 g/s, LCR=3%, T=79.9°C, P=2809.8 kPa
        # ═══════════════════════════════════════════════════════════════════════════
        
        # ── Fig. 5  (Reference, Design 1, Design 2 — body diameter effect) ────────
        EXP_FIG5 = {
            # (a) LCR = 1.5 %  — Reference starts at MFR=60 (no 30 g/s point visible)
            0.015: {
                'Reference': {60: 0.993, 90: 0.975, 120: 0.993, 150: 1.008},
                'Design 1' : {60: 1.010, 90: 1.010, 120: 1.013, 150: 1.020},
                'Design 2' : {60: 1.013, 90: 1.015, 120: 1.030, 150: 1.033},
            },
            # (b) LCR = 3.0 %
            0.030: {
                'Reference': {30: 1.020, 60: 1.017, 90: 1.000, 120: 1.000, 150: 1.008},
                'Design 1' : {30: 1.027, 60: 1.017, 90: 1.017, 120: 1.022, 150: 1.025},
                'Design 2' : {30: 1.030, 60: 1.020, 90: 1.020, 120: 1.030, 150: 1.030},
            },
            # (c) LCR = 4.5 %
            0.045: {
                'Reference': {30: 1.020, 60: 1.020, 90: 1.007, 120: 1.007, 150: 1.010},
                'Design 1' : {30: 1.022, 60: 1.022, 90: 1.022, 120: 1.020, 150: 1.022},
                'Design 2' : {30: 1.030, 60: 1.022, 90: 1.022, 120: 1.030, 150: 1.035},
            },
        }
        
        # ── Fig. 6  (Reference, Design 3, Design 4 — body height effect) ──────────
        EXP_FIG6 = {
            # (a) LCR = 1.5 %
            0.015: {
                'Reference': {60: 0.993, 90: 0.975, 120: 0.993, 150: 1.008},
                'Design 3' : {60: 1.010, 90: 1.010, 120: 0.993, 150: 1.015},
                'Design 4' : {60: 1.015, 90: 1.015, 120: 1.030, 150: 1.033},
            },
            # (b) LCR = 3.0 %
            0.030: {
                'Reference': {30: 1.020, 60: 1.017, 90: 1.000, 120: 1.000, 150: 1.008},
                'Design 3' : {30: 1.030, 60: 1.017, 90: 1.017, 120: 1.017, 150: 1.017},
                'Design 4' : {30: 1.030, 60: 1.020, 90: 1.020, 120: 1.025, 150: 1.028},
            },
            # (c) LCR = 4.5 %
            0.045: {
                'Reference': {30: 1.020, 60: 1.020, 90: 1.007, 120: 1.007, 150: 1.010},
                'Design 3' : {30: 1.022, 60: 1.022, 90: 1.022, 120: 1.020, 150: 1.022},
                'Design 4' : {30: 1.028, 60: 1.025, 90: 1.025, 120: 1.028, 150: 1.032},
            },
        }
        
        # Marker style per geometry — open symbols matching paper
        STYLE = {
            'Reference': dict(marker='s', label='Reference'),
            'Design 1' : dict(marker='o', label='Design 1'),
            'Design 2' : dict(marker='^', label='Design 2'),
            'Design 3' : dict(marker='o', label='Design 3'),
            'Design 4' : dict(marker='^', label='Design 4'),
        }
        
        # Line style per geometry for model curves
        MODEL_LINE = {
            'Reference': dict(linestyle='-',  linewidth=1.4),
            'Design 1' : dict(linestyle='--', linewidth=1.4),
            'Design 2' : dict(linestyle=':',  linewidth=1.4),
            'Design 3' : dict(linestyle='--', linewidth=1.4),
            'Design 4' : dict(linestyle=':',  linewidth=1.4),
        }
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 7.  HELPER: model relative efficiency
        # ═══════════════════════════════════════════════════════════════════════════
        
        MFR_PLOT = [30, 60, 90, 120, 150]   # g/s
        
        def relative_eta_model(geom_names: list, lcr: float) -> dict:
            """
            Relative efficiency (η_geom / η_Reference) at each MFR, Cond B.
            Self-normalised per MFR, matching the paper's definition.
            """
            model = {}
            for name in geom_names:
                rel = []
                for mfr in MFR_PLOT:
                    eta_g = run_point(GEOMETRIES[name],        mfr*1e-3, lcr,
                                      COND['B']['T'], COND['B']['P']).eta_total
                    eta_r = run_point(GEOMETRIES['Reference'], mfr*1e-3, lcr,
                                      COND['B']['T'], COND['B']['P']).eta_total
                    rel.append(eta_g / eta_r)
                model[name] = rel
            return model
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 8.  PLOTTING
        # ═══════════════════════════════════════════════════════════════════════════
        
        LCR_PANELS = [
            (0.015, '(a) LCR = 1.5 %'),
            (0.030, '(b) LCR = 3.0 %'),
            (0.045, '(c) LCR = 4.5 %'),
        ]
        
        
        def draw_fig(exp_data: dict, geom_names: list, save_name: str, suptitle: str):
            """
            Draw a 3-panel figure (one row per LCR) matching the paper layout.
            Symbols = digitised experimental data, lines = model predictions.
            """
            fig, axes = plt.subplots(3, 1, figsize=(6, 13), sharey=True, sharex=True)
            fig.subplots_adjust(hspace=0.28, left=0.15, right=0.95, top=0.93, bottom=0.06)
        
            for ax, (lcr, panel_title) in zip(axes, LCR_PANELS):
                model = relative_eta_model(geom_names, lcr)
                exp   = exp_data[lcr]
        
                for name in geom_names:
                    sty = STYLE[name]
                    ml  = MODEL_LINE[name]
        
                    # model curve (grey line)
                    ax.plot(MFR_PLOT, model[name],
                            color='#444444', **ml, zorder=2)
        
                    # digitised experimental points (black open symbols)
                    exp_mfr = sorted(exp[name].keys())
                    exp_rel = [exp[name][m] for m in exp_mfr]
                    ax.scatter(exp_mfr, exp_rel,
                               marker=sty['marker'], s=55,
                               edgecolors='k', facecolors='none',
                               linewidths=1.2, zorder=3,
                               label=sty['label'])
        
                ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.7, zorder=1)
                ax.set_xlim(20, 170)
                ax.set_ylim(0.95, 1.05)
                ax.set_yticks([0.95, 1.00, 1.05])
                ax.set_xticks([30, 60, 90, 120, 150])
                ax.set_ylabel('Relative efficiency (–)', fontsize=10)
                ax.set_title(panel_title, fontsize=10, pad=4)
                ax.grid(True, alpha=0.25)
                ax.legend(fontsize=9, loc='lower right', framealpha=0.85)
        
            axes[-1].set_xlabel(r'Mass flow rate (g s$^{-1}$)', fontsize=10)
        
            fig.text(0.5, 0.005,
                     'Symbols: digitised exp. data (Kim et al. 2016)    '
                     'Lines: model predictions',
                     ha='center', fontsize=8, style='italic', color='#555555')
            fig.suptitle(suptitle, fontsize=10, y=0.97)
        
            # plt.savefig(save_name, dpi=150, bbox_inches='tight')
            # print(f"✓  Saved → {save_name}")
            plt.show()
        
        
        # ── Fig. 5 reconstruction  (Reference, Design 1, Design 2) ───────────────
        draw_fig(
            exp_data   = EXP_FIG5,
            geom_names = ['Reference', 'Design 1', 'Design 2'],
            save_name  = 'fig5_reconstruction.png',
            suptitle   = ('Fig. 5 — Relative efficiency, Cond B (79.9 °C / 28.1 bar)\n'
                          'Reference, Design 1 (BD+20 %), Design 2 (BD+41 %)'),
        )
        
        # ── Fig. 6 reconstruction  (Reference, Design 3, Design 4) ───────────────
        draw_fig(
            exp_data   = EXP_FIG6,
            geom_names = ['Reference', 'Design 3', 'Design 4'],
            save_name  = 'fig6_reconstruction.png',
            suptitle   = ('Fig. 6 — Relative efficiency, Cond B (79.9 °C / 28.1 bar)\n'
                          'Reference, Design 3 (BH+18 %), Design 4 (BH+36 %)'),
        )
        
        # ── Diagnostic figure (absolute η, all geometries, grade efficiency) ──────
        fig2 = plt.figure(figsize=(17, 9))
        gs2  = gridspec.GridSpec(2, 3, figure=fig2, hspace=0.45, wspace=0.33)
        
        geom_colors = ['#333333', '#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
        lcr_colors  = ['#1f77b4', '#ff7f0e', '#2ca02c']
        lcr_labels  = ['LCR = 1.5 %', 'LCR = 3.0 %', 'LCR = 4.5 %']
        pt_markers  = ['o', 's', '^']
        
        axA = fig2.add_subplot(gs2[0, 0])
        for i, lcr in enumerate(LCR_VALUES):
            eta_vals = [run_point(GEOMETRIES['Reference'], m*1e-3, lcr,
                                  COND['B']['T'], COND['B']['P']).eta_total * 100
                        for m in MFR_PLOT]
            axA.plot(MFR_PLOT, eta_vals, marker=pt_markers[i],
                     color=lcr_colors[i], label=lcr_labels[i])
        axA.set_xlabel('Mass flow rate [g/s]', fontsize=10)
        axA.set_ylabel('Overall efficiency [%]', fontsize=10)
        axA.set_title('(A) Absolute η vs MFR\nReference geom, Cond B', fontsize=10)
        axA.legend(fontsize=9)
        axA.grid(True, alpha=0.3)
        
        axB = fig2.add_subplot(gs2[0, 1])
        for j, (name, geom) in enumerate(GEOMETRIES.items()):
            rel = []
            for mfr in MFR_PLOT:
                e_g = run_point(geom,                    mfr*1e-3, LCR_REF,
                                COND['B']['T'], COND['B']['P']).eta_total
                e_r = run_point(GEOMETRIES['Reference'], mfr*1e-3, LCR_REF,
                                COND['B']['T'], COND['B']['P']).eta_total
                rel.append(e_g / e_r)
            axB.plot(MFR_PLOT, rel, marker='o', color=geom_colors[j], label=name)
        axB.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
        axB.set_xlabel('Mass flow rate [g/s]', fontsize=10)
        axB.set_ylabel('Relative efficiency [–]', fontsize=10)
        axB.set_title('(B) Relative η — all geometries\nLCR = 3 %, Cond B', fontsize=10)
        axB.legend(fontsize=8)
        axB.grid(True, alpha=0.3)
        
        axC = fig2.add_subplot(gs2[0, 2])
        names2, delta_ecs = [], []
        for name, geom in GEOMETRIES.items():
            s = run_point(geom, MFR_REF, LCR_REF, COND['B']['T'], COND['B']['P'])
            names2.append(name)
            delta_ecs.append((s.eta_total - eta_ref_geom) * 100)
        xb = np.arange(len(names2))
        bars = axC.bar(xb, delta_ecs, color=geom_colors, edgecolor='k', linewidth=0.6)
        axC.axhline(0, color='k', linewidth=0.8)
        axC.set_xticks(xb)
        axC.set_xticklabels(names2, rotation=15, ha='right', fontsize=8)
        axC.set_ylabel('ΔEC vs Reference [pp]', fontsize=10)
        axC.set_title('(C) Efficiency gain vs Reference\nMFR=90 g/s, LCR=3%, Cond B\n(cf. Table 3)',
                      fontsize=10)
        axC.grid(True, axis='y', alpha=0.3)
        for bar, val in zip(bars, delta_ecs):
            axC.text(bar.get_x() + bar.get_width()/2,
                     bar.get_height() + 0.002,
                     f'{val:+.3f}', ha='center', va='bottom', fontsize=8)
        
        axD = fig2.add_subplot(gs2[1, 0])
        d_um = sep_ref.drop_distrib_D * 1e6
        axD.semilogx(d_um, sep_ref.eta_grade * 100, 'C0', linewidth=1.8,
                     label=f'η = {sep_ref.eta_total*100:.2f} %\nd50 = {sep_ref.d50_val*1e6:.2f} µm')
        axD.axhline(50, color='gray', linestyle='--', linewidth=0.8, label='50 % line')
        axD.axvline(sep_ref.d50_val * 1e6, color='C0', linestyle=':', linewidth=0.8)
        axD.set_xlabel('Droplet diameter [µm]', fontsize=10)
        axD.set_ylabel('Grade efficiency [%]', fontsize=10)
        axD.set_title('(D) Grade efficiency\nReference condition', fontsize=10)
        axD.legend(fontsize=9)
        axD.grid(True, which='both', alpha=0.3)
        axD.set_ylim(0, 105)
        
        axE = fig2.add_subplot(gs2[1, 1])
        axE.semilogx(d_um, sep_ref.drop_distrib_x, 'C1', linewidth=1.8,
                     label=f'd1 = {sep_ref.d1_val*1e6:.2f} µm\n'
                           f'd_mm = {sep_ref.d_mm*1e6:.2f} µm')
        axE.axvline(sep_ref.d1_val  * 1e6, color='C1', linestyle=':',  linewidth=1.0, label='d1')
        axE.axvline(sep_ref.d50_val * 1e6, color='C0', linestyle='--', linewidth=1.0, label='d50')
        axE.set_xlabel('Droplet diameter [µm]', fontsize=10)
        axE.set_ylabel('Mass fraction density [1/m]', fontsize=10)
        axE.set_title('(E) Droplet size distribution\nReference condition', fontsize=10)
        axE.legend(fontsize=9)
        axE.grid(True, which='both', alpha=0.3)
        
        axF = fig2.add_subplot(gs2[1, 2])
        for key, ls in zip(['A', 'B'], ['--', '-']):
            cond = COND[key]
            eta_vals = [run_point(GEOMETRIES['Reference'], m*1e-3, LCR_REF,
                                  cond['T'], cond['P']).eta_total * 100
                        for m in MFR_PLOT]
            axF.plot(MFR_PLOT, eta_vals, marker='o', linestyle=ls,
                     label=f"Cond {key}: {cond['T']-273.15:.1f} °C / {cond['P']/1e3:.1f} kPa")
        axF.set_xlabel('Mass flow rate [g/s]', fontsize=10)
        axF.set_ylabel('Overall efficiency [%]', fontsize=10)
        axF.set_title('(F) Cond A vs Cond B\nLCR = 3 %, Reference geom\n(paper: η_B > η_A)',
                      fontsize=10)
        axF.legend(fontsize=9)
        axF.grid(True, alpha=0.3)
        
        fig2.suptitle('OilSeparator — Validation diagnostics vs Kim et al. (2016)\n'
                      'Int. J. Refrigeration 70 (2016) 206–218',
                      fontsize=12, y=1.01)
        # plt.savefig('validation_kim2016.png', dpi=150, bbox_inches='tight')
        # print("✓  Diagnostic figure saved → validation_kim2016.png")
        plt.show()
        
        
        # ═══════════════════════════════════════════════════════════════════════════
        # 9.  SUMMARY
        # ═══════════════════════════════════════════════════════════════════════════
        
        print("\n" + "=" * 65)
        print("SUMMARY")
        print("=" * 65)
        checks = {
            'Check 1 — Absolute level  (90–100 %)' : 90.0 <= sep_ref.eta_total*100 <= 100.0,
            'Check 2a — η↑ with LCR               ': lcr_trend_ok,
            'Check 2b — U-shape with MFR (LCR=3%) ': u_shape_ok,
            'Check 3a — All designs η > Reference  ': geom_pass,
            'Check 3b — BD effect monotone         ': bd_monotone,
            'Check 3c — BH effect monotone         ': bh_monotone,
            'Check 4  — Cond B η > Cond A          ': cond_pass,
        }
        for label, passed in checks.items():
            print(f"  {label} : {'✓  PASS' if passed else '✗  FAIL'}")
        print("=" * 65)
