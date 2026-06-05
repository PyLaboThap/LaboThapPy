import numpy as np
from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
import CoolProp.CoolProp as CP
from scipy.optimize import brentq
from labothappy.correlations.properties.void_fraction import void_fraction

class EjectorCstEff(BaseComponent):
    """
    **Component**: Ejector

    **Model**: Constant isentropic efficiency model

    **Reference**: Analysis of a Two Phase Flow Ejector For Transcritical CO2 Cycle (2008)
                   Fang Liu, Eckhard A. Groll

    **Description**:

        This component models an ejector using a constant isentropic efficiency. 
        Given the suction conditions (pressures, enthalpies, fluid) and the exhaust pressure,
        it calculates the exhaust specific enthalpy and exhaust temperature.
        This simple model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.
        - Isentropic efficiency stays constant for all the conditions.
        - Negligible heat losses and mechanical losses except those accounted for by efficiency.

    **Connectors**:
        
        - su_1 : Motive nozzle supply
        - su_2 : Suction chamber supply
        - ex : Diffuser exhaust

    **Parameters**:


    **Inputs**:



    **Ouputs**:


    """

    def __init__(self):
        super().__init__()
        
        # Define mass flow connectors for suction and exhaust
        self.su_1 = MassConnector()
        self.su_2 = MassConnector()
        self.ex = MassConnector()  
        
        # Define work connector for mechanical work output
        self.print_flag = 1

    def get_required_inputs(self):  # Used in check_calculablle to see if all of the required inputs are set
        # Return a list of required inputs
        return []

    def get_required_parameters(self):
        # Return a list of model parameters required for solving
        return []

    def _kornhauser(self, h_mb, P_b):
        """
        HEM speed of sound — Kornhauser (1990).
        Partial derivatives computed analytically via CoolProp.
        Returns None if fluid is not two-phase or if den <= 0.
        """
    
        # ── Saturation properties ──────────────────────────────────────────
        self.AS.update(CP.PQ_INPUTS, P_b, 0)
        s_f = self.AS.smass()
        h_f = self.AS.hmass()
        v_f = 1.0 / self.AS.rhomass()
    
        self.AS.update(CP.PQ_INPUTS, P_b, 1)
        s_g = self.AS.smass()
        h_g = self.AS.hmass()
        v_g = 1.0 / self.AS.rhomass()
    
        x = (h_mb - h_f) / (h_g - h_f)
        x = np.clip(x, 0.0, 1.0)
    
        # ── Partial derivatives along saturation curve ─────────────────────
        # dX/dP|_sat via chain rule:
        #   dX/dP|_sat = dX/dP|_T + dX/dT|_P * dT/dP|_sat
        #
        # dT/dP|_sat = Clausius-Clapeyron = T*(v_g - v_f)/(h_g - h_f)
    
        def sat_derivs(quality):
            self.AS.update(CP.PQ_INPUTS, P_b, quality)
            T_sat = self.AS.T()
        
            # Clausius-Clapeyron
            dT_dP_sat = T_sat * (v_g - v_f) / (h_g - h_f)
        
            dh_dP_T = self.AS.first_partial_deriv(CP.iHmass, CP.iP, CP.iT)
            dh_dT_P = self.AS.first_partial_deriv(CP.iHmass, CP.iT, CP.iP)
            ds_dP_T = self.AS.first_partial_deriv(CP.iSmass, CP.iP, CP.iT)
            ds_dT_P = self.AS.first_partial_deriv(CP.iSmass, CP.iT, CP.iP)
        
            # v = 1/rho  →  dv/dP = -1/rho² * drho/dP
            rho     = self.AS.rhomass()
            drho_dP_T = self.AS.first_partial_deriv(CP.iDmass, CP.iP, CP.iT)
            drho_dT_P = self.AS.first_partial_deriv(CP.iDmass, CP.iT, CP.iP)
            dv_dP_T   = -drho_dP_T / rho**2
            dv_dT_P   = -drho_dT_P / rho**2
        
            dh_dP_sat = dh_dP_T + dh_dT_P * dT_dP_sat
            ds_dP_sat = ds_dP_T + ds_dT_P * dT_dP_sat
            dv_dP_sat = dv_dP_T + dv_dT_P * dT_dP_sat
        
            return dh_dP_sat, ds_dP_sat, dv_dP_sat
    
        dh_f_dP, ds_f_dP, dv_f_dP = sat_derivs(0)
        dh_g_dP, ds_g_dP, dv_g_dP = sat_derivs(1)
    
        # ── Isentropic derivatives of mixture (Kornhauser 1990) ────────────
        dx_dP_s = -(x * ds_g_dP + (1-x) * ds_f_dP) / (s_g - s_f)
    
        h_prime = x*dh_g_dP + (1-x)*dh_f_dP + (h_g - h_f)*dx_dP_s
        v_prime = x*dv_g_dP + (1-x)*dv_f_dP + (v_g - v_f)*dx_dP_s
    
        # ── Kornhauser formula ─────────────────────────────────────────────
        v_mix = 1.0 / self.rho_mb
        num   = v_mix**2 * (h_g - h_f)
        den   = (v_g - v_f) * h_prime - v_prime * (h_g - h_f)
    
        if den <= 0:
            return None
    
        return np.sqrt(num / den)

    def motive_nozzle(self, P_t):
        """
        Iterates on P_t to find exit nozzle conditions (considered sonic)
        """
        # Outlet enthalpy
        self.AS.update(CP.PSmass_INPUTS, P_t, self.su_1.s)        
        h_mb_is = self.AS.hmass()
        h_mb    = self.su_1.h - self.params['eta_m'] * (self.su_1.h - h_mb_is)
    
        # Outlet speed
        dh = self.su_1.h - h_mb
        if dh <= 0:
            return -1e6   # pas de détente possible à cette P_b
        v_mb = np.sqrt(2 * dh)
    
        # Phase and density
        self.AS.update(CP.HmassP_INPUTS, h_mb, P_t)
        self.rho_mb = self.AS.rhomass()
        phase  = self.AS.phase()
        
        # Speed of sound
        if phase != CP.iphase_twophase:
            self.a_mb = self.AS.speed_sound()
        else:
            self.a_mb = self._kornhauser(h_mb, P_t)
    
        return v_mb - self.a_mb  # = 0 à la condition sonique

    def suction_chamber(self, P_b):
        """
        Suction chamber — iterates on P_b.
        
        P_b must be slightly below P_su2 to drive secondary flow.
        P_b is NOT equal to P_mix — the secondary stream is further
        compressed/expanded between P_b and P_mix in the mixing section.
        
        In design mode with ω imposed and ṁ_p known :
            ṁ_s = ω * ṁ_p  is known
            A_sb is unknown → iterate on P_b such that
            the secondary stream reaches P_b with the correct ṁ_s
            
        Residual : ṁ_s_energy(P_b) - ṁ_s_target = 0
        where ṁ_s_energy is computed from energy eq. + A_sb constraint.
        
        But without A_sb fixed, the correct residual is :
            P_b must be consistent with P_mix from motive stream.
        
        Here : P_b = P_mix (common mixing pressure, Li & Groll assumption)
        with P_b slightly below P_su2 as physical check.
        """
        fluid = self.AS.name()
    
        # Physical check : P_b must be below P_su2
        if P_b >= self.su_2.p:
            return 1e6   # no flow possible
    
        # ── Isentropic then real enthalpy ─────────────────────────────────
        self.AS.update(CP.PSmass_INPUTS, P_b, self.su_2.s)
        h_sb_is = self.AS.hmass()
        h_sb    = self.su_2.h - self.params['eta_s'] * (self.su_2.h - h_sb_is)
    
        dh = self.su_2.h - h_sb
        if dh <= 0:
            return 1e6
        v_sb = np.sqrt(2 * dh)
    
        # ── Density at (h_sb, P_b) ────────────────────────────────────────
        rho_sb = CP.PropsSI('D', 'H', h_sb, 'P', P_b, fluid)
    
        # ── Area of secondary stream at P_b ───────────────────────────────
        # A_sb = A_mix - A_mb  where A_mix is from motive stream
        # In design mode : A_sb emerges, not fixed
        # Residual : momentum balance at mixing section
        # Both streams at P_b → momentum balance gives P_mix
        # P_mix must equal P_b (constant pressure mixing assumption)
    
        # Momentum balance residual :
        # (ṁ_p + ṁ_s) * V_mix = ṁ_p * V_mb + ṁ_s * V_sb
        w     = self.params['w']
        V_mix = (self.v_mb + w * v_sb) / (1.0 + w)
    
        # Energy balance → h_mix
        h0_mix = (self.su_1.h + w * self.su_2.h) / (1.0 + w)
        h_mix  = h0_mix - V_mix**2 / 2.0
    
        # Density of mixture at P_b
        rho_mix = CP.PropsSI('D', 'H', h_mix, 'P', P_b, fluid)
    
        # Area balance : A_mix from mixture, A_mb + A_sb from individual streams
        A_mix_momentum = (1.0 + w) * self.su_1.m_dot / (rho_mix * V_mix)
        A_mb           = self.su_1.m_dot / (self.rho_mb * self.v_mb)
        A_sb           = self.su_2.m_dot / (rho_sb * v_sb)
        A_mix_geometry = A_mb + A_sb
    
        # Store for next step
        self.h_sb   = h_sb
        self.v_sb   = v_sb
        self.rho_sb = rho_sb
        self.A_sb   = A_sb
        self.A_mb   = A_mb
    
        # Residual : area from momentum balance = area from geometry
        return A_mix_momentum - A_mix_geometry   # = 0 at solution

    def mixing_section(self, v_mix):
        
        # Energy conservation
        h_mix = self.h0_mix-v_mix**2/2
        
        # Mass conservation
        self.A_mix = self.A_sb + self.A_mb
        self.rho_mix = (self.su_1.m_dot + self.su_2.m_dot)/(self.A_mix*v_mix)
        
        # Momentum conservation
        eta_mix = self.params['eta_mix']
        m_dot_t = self.su_1.m_dot + self.su_2.m_dot
        
        mom_p   = self.P_b_sonic   * self.A_mb
        mom_s   = self.P_b_suction * self.A_sb
        flux_p  = eta_mix * self.su_1.m_dot * self.v_mb
        flux_s  = eta_mix * self.su_2.m_dot * self.v_sb
        flux_mix= eta_mix * m_dot_t * v_mix
        mom_in  = mom_p + flux_p + mom_s + flux_s
        
        self.P_mix = (self.P_b_sonic*self.A_mb + eta_mix*self.rho_mb*self.A_mb*self.v_mb**2 + self.P_b_suction*self.A_sb + eta_mix*self.rho_sb*self.A_sb*self.v_sb**2)/self.A_mix - self.rho_mix*v_mix**2 
        
        # self.AS.update(CP.DmassP_INPUTS, self.rho_mix, self.P_mix)
        # self.h_mix = self.AS.hmass()
        self.h_mix = CP.PropsSI('H', 'D', self.rho_mix, 'P', self.P_mix, self.su_1.fluid)

        return h_mix - self.h_mix


    def solve(self):
        # Perform checks to ensure the model can be calculated and has parameters
        self.check_calculable()
        self.check_parametrized()

        self.AS = CP.AbstractState('HEOS', self.su_1.fluid)

        if not (self.calculable and self.parametrized):
            self.solved = False
            print(
                "EjectorCstEff could not be solved. It is not calculable and/or not parametrized"
            )
            return
        
        "1) Motive nozzle"        
        # ── Recherche de P_b par Brent ──────────────────────────────────────────
        # Bornes : P_evap < P_b < P_gc
        # À P_b proche de P_gc : v_mb ≈ 0 << a_mb → résidu négatif
        # À P_b proche de P_evap : v_mb grand >> a_mb → résidu positif
        P_low  = self.su_2.p * 1.01          # légèrement au-dessus de P_evap
        P_high = self.su_1.p * 0.99          # légèrement en dessous de P_gc
        
        r_low  = self.motive_nozzle(P_low)
        r_high = self.motive_nozzle(P_high)

        if r_low * r_high > 0:
            raise ValueError("Could not solve motive nozzle, verify inlet pressures")
        
        self.P_b_sonic = brentq(self.motive_nozzle, P_low, P_high, xtol=1e2, rtol=1e-6)
        
        # ── État final à P_b sonique ────────────────────────────────────────────
        self.AS.update(CP.PSmass_INPUTS, self.P_b_sonic, self.su_1.s)
        h_mb_is       = self.AS.hmass()
        self.h_mb     = self.su_1.h - self.params['eta_m'] * (self.su_1.h - h_mb_is)
        self.v_mb     = np.sqrt(2 * (self.su_1.h - self.h_mb))
        
        self.AS.update(CP.HmassP_INPUTS, self.h_mb, self.P_b_sonic)
        self.rho_mb   = self.AS.rhomass()
        self.P_b      = self.P_b_sonic
            
        self.A_mb = 1 / (self.rho_mb*self.v_mb*self.su_1.m_dot)

        "2) Suction chamber"
        
        self.params['w'] = self.su_2.m_dot/self.su_1.m_dot # Entrainement ratio
        T_triple = self.AS.Ttriple()
        
        self.AS.update(CP.QT_INPUTS, 0, T_triple)
        P_Triple = self.AS.p()
        
        self.P_b_suction = brentq(self.suction_chamber, P_Triple, self.su_2.p*0.999, xtol=1e-2, rtol=1e-6)
        
        "3) Mixing section"
        
        h0_mb  = self.h_mb + self.v_mb**2 / 2.0
        h0_sb  = self.h_sb + self.v_sb**2 / 2.0
        self.h0_mix = (self.su_1.m_dot * h0_mb + self.su_2.m_dot * h0_sb) \
                / (self.su_1.m_dot + self.su_2.m_dot)
                
        converged = False
        self.v_min = 1
        
        while converged is False:
            if self.v_min == self.v_mb:
                raise ValueError("Mixing section could not find a solution")
                return
    
            try:
                self.v_mix = brentq(self.mixing_section, self.v_min, self.v_mb, xtol=1e-5, rtol=1e-6)
                converged = True
            except:
                self.v_min += 1
        
        "4) Diffuser flow"
        
        self.P_ex = self.params['C_t']*(0.5*self.rho_mix*self.v_mix**2) + self.P_mix # pressure recovery coefficient
        self.h_ex = self.h_mix
        self.m_dot_ex = self.su_1.m_dot + self.su_2.m_dot
        
        self.AS.update(CP.PQ_INPUTS, self.P_mix, 0)
        rho_l_mix = self.AS.rhomass()
        
        self.AS.update(CP.PQ_INPUTS, self.P_mix, 1)
        rho_g_mix = self.AS.rhomass()

        self.AS.update(CP.HmassP_INPUTS, self.h_mix, self.P_mix)
        self.x_mix = self.AS.Q()

        eps = void_fraction(self.x_mix, rho_g_mix, rho_l_mix)[0]
        
        def res_C_t(A_d):
            # Rouhani (1969) - Diffusing a homogenized two-phase flow - I. OWEN, A. ABDUL-GHANI and A. M. AMINI
            C_t_comp = 2*self.A_mix/A_d * self.rho_mix * (1-self.A_mix/A_d)*(1/(rho_l_mix*(1-eps)**2))
            return self.params['C_t'] - C_t_comp
        
        A_d_id = self.A_mix/np.sqrt(1 - 0.8)
        self.A_d = brentq(res_C_t, self.A_mix, A_d_id, xtol=1e-10, rtol=1e-10)
        
        self.params['A_d'] = self.A_d
        self.params['A_mb'] = self.A_mb
        self.params['A_sb'] = self.A_sb
        self.params['A_mix'] = self.A_mix
        
        self.update_connectors(self.P_ex, self.h_ex, self.m_dot_ex)
        
        return

    def update_connectors(self, P_ex, h_ex, m_dot_ex):
        """Update the connectors with the calculated values."""
        
        self.ex.reset()
        self.ex.set_fluid(self.su_1.fluid)
        self.ex.set_h(h_ex)
        self.ex.set_p(P_ex)
        self.ex.set_m_dot(m_dot_ex)

    def print_results(self):
        print("=== Expander Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - W_dot_exp: {self.W.W_dot} [W]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Expander Results ===")
        print("Mass connectors:")
        print(
            f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]"
        )
        print(
            f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]"
        )
        print("=========================")
        print("Work connector:")
        print(f"  - W_dot_exp: {self.W.W_dot} [W]")
        print("=========================")

if __name__ == "__main__":

    ej = EjectorCstEff()
    
    P_su_1 = 100
    P_su_2 = 39
    
    ej.set_inputs(
        T_su_1 = 40 + 273.15, # K
        P_su_1 = P_su_1*1e5, # Pa
        m_dot_su_1 = 0.18, # kg/s
        fluid_su_1 = 'CO2',
        
        T_su_2 = 15 + 273.15, # K
        P_su_2 = P_su_2*1e5, # Pa
        m_dot_su_2 = 0.07, # kg/s
        fluid_su_2 = 'CO2',
        )
    
    ej.set_parameters(
        eta_m = 0.9,
        eta_s = 0.9,
        eta_mix = 1,
        C_t = 0.8,
        )

    ej.solve()
    