# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 14:51:32 2026

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

class DrumSeparator(BaseComponent):
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

#%% 

    def compute_terminal_velocity(self, dp, max_iter=50, tol=1e-6):
        """
        Terminal settling velocity of a liquid droplet in gas.
        Uses iterative Schiller-Naumann / Newton drag correlation.
    
        Parameters
        ----------
        dp : droplet diameter [m]
    
        Returns
        -------
        Updates self.velocity_data with Vt [m/s], Re_p [-], Cd [-], regime [str]
    
        References
        ----------
        Schiller & Naumann (1933) — valid for Re_p < 800
        Newton (1687)             — valid for Re_p > 1000
        Clift, Grace & Weber, "Bubbles, Drops and Particles", 1978
        """
        g         = 9.81                        # m/s²
        rho_g     = self.su_gas['D']            # kg/m³
        mu_g      = self.su_gas['V']            # Pa·s
        rho_l     = self.su_oil['D']            # kg/m³
        delta_rho = rho_l - rho_g              # kg/m³
    
        # ── Archimède number: allows direct regime pre-check ─────────────────
        Ar = rho_g * delta_rho * g * dp**3 / mu_g**2
    
        # ── Initial guess: always Stokes ──────────────────────────────────────
        Vt = (g * dp**2 * delta_rho) / (18.0 * mu_g)
    
        # ── Iterative solve ───────────────────────────────────────────────────
        solved = False
        for _ in range(max_iter):
    
            Re_p = rho_g * Vt * dp / mu_g
    
            # Drag coefficient — unified Schiller-Naumann up to Re=1000, Newton above
            if Re_p < 1e-6:
                Cd = 24.0 / 1e-6            # avoid division by zero
            elif Re_p <= 1000.0:
                Cd = (24.0 / Re_p) * (1.0 + 0.15 * Re_p**0.687)
            else:
                Cd = 0.44
    
            # Terminal velocity from force balance
            Vt_new = np.sqrt((4.0 / 3.0) * g * dp * delta_rho / (Cd * rho_g))
    
            # Convergence check
            if abs(Vt_new - Vt) / (abs(Vt) + 1e-30) < tol:
                Vt = Vt_new
                solved = True
                break
    
            Vt = 0.5 * (Vt + Vt_new)       # damped update
    
        # ── Final Re and regime ───────────────────────────────────────────────
        Re_p = rho_g * Vt * dp / mu_g
    
        if Re_p < 1.0:
            regime = "Stokes"
        elif Re_p < 1000.0:
            regime = "Intermediate (Schiller-Naumann)"
        else:
            regime = "Newton"
    
        if not solved:
            import warnings
            warnings.warn(f"compute_terminal_velocity: no convergence for dp={dp*1e6:.1f} µm")
    
        self.velocity_data = {
            'Vt'    : Vt,
            'Re_p'  : Re_p,
            'Cd'    : Cd,
            'Ar'    : Ar,
            'regime': regime,
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

        self.d1_val = self.params['D_in']*7.43*1e-3 * (self.su_oil['sigma']/(self.su_gas['D']*self.params['D_in']*self.u_in))**0.5 * (self.u_in*self.params['D_in']/self.su_gas['NU'])**0.1
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

    def overall_efficiency_vertical(self):
        """
        Overall separation efficiency of a vertical knock-out drum,
        integrated over the droplet size distribution.
    
        η_total = ∫ η(dp) · f(dp) d(dp)
    
        Parameters
        ----------
        d_array  : droplet diameter array [m]
        MF_array : mass fraction density [m⁻¹]
        Vg       : superficial gas velocity [m/s]
    
        Returns
        -------
        dict with eta_total, dp_100, eta_grade array, regime array
        """
        self.eta_grade  = np.zeros(len(self.drop_distrib_D))
        self.V_t  = np.zeros(len(self.drop_distrib_D))
        self.Vel_data   = []
    
        for i, dp in enumerate(self.drop_distrib_D):
            self.compute_terminal_velocity(dp)
            
            if self.velocity_data['Vt'] > self.u_g: # Flow does not entrain the drop : captured
                self.eta_grade[i] = 1
            else:
                self.eta_grade[i] = 0 # Flow entrains the drop : not captured
            
            self.V_t[i] = self.velocity_data['Vt']
            self.Vel_data.append(self.velocity_data)
    
        # Overall efficiency (mass-weighted integral)
        self.eta_total = trapezoid(self.eta_grade * self.drop_distrib_x, self.drop_distrib_D) / trapezoid(self.drop_distrib_x, self.drop_distrib_D)
    
        # Cut diameter
        for i in range(len(self.eta_grade)):
            eta = self.eta_grade[i]
            
            if eta == 1:
                break
        
        self.dp_100 = self.drop_distrib_D[i]
        
        
        # dp_100 = cut_diameter_vertical(rho_L, rho_G, mu_G, Vg, g)
    
        return 


#%%

    def solve(self):
        
        "1) Prepare Inputs"
        
        self.compute_complete_inputs()
        
        "2) Compute Flow and Inlet speeds"
        
        self.A_vessel = np.pi * (self.params['BD'] / 2.0)**2   # [m²]
        Vdot_tot = self.su_gas['V_dot'] + self.su_oil['V_dot']
        self.u_g = Vdot_tot / self.A_vessel              # superficial gas velocity [m/s]
        
        self.A_in = np.pi * (self.params['D_in'] / 2.0)**2   # [m²]
        self.u_in = Vdot_tot / self.A_in              # superficial gas velocity [m/s]

        "3) Average droplet diameter" # Shall maybe compute median
        
        self.compute_d1()
        self.droplet_size_distribution()
        
        self.d_mm = trapezoid(self.drop_distrib_D * self.drop_distrib_x, self.drop_distrib_D) / trapezoid(self.drop_distrib_x, self.drop_distrib_D)

        "4) Compute separation efficiency"        
        self.overall_efficiency_vertical()
        
        self.mdot_ex_oil = self.su_oil['m_dot']*(1-self.eta_total)
        
        return

#%%

if __name__ == "__main__":
    
    case_study = "Simu"
    
    if case_study == "Simu":
        Drum = DrumSeparator()
    
        Drum.set_inputs(
            T_su_gas = 50+273.15, # K
            P_su_gas = 1*1e5, # Pa
            V_dot_gas = 52/3600, # m^3/s
            fluid_gas = 'Air',
            
            fluid_oil = 'Jet II',
            T_su_oil = 50+273.15, # K
            P_su_oil = 1*1e5, # Pa
            V_dot_oil = (650/1000)/3600, # m^3/s
            rho_oil = 1003.5, # kg/m^3
            nu_oil = [[40+273.15,100+273.15], [27.6, 5.1]], # [T1, T2] : K and [nu1, nu_2] : cSt
            sigma_oil = 3*1e-3
            )
        
        Drum.set_parameters(
            BD=0.5,   # m body diameter
            # BH=0.22,   # m body height
            D_in=0.0141, # inlet pipe inner diameter
            # S=0.120,    # m outlet pipe length
            # IP=0.030,   # m inlet position
            # BD_ref= 0.085, # Body diameter of reference 
            # BH_ref= 0.22, # Body Height of reference 
            )
        
        Drum.solve()
        
    
    