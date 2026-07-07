import __init__

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import minimize, brentq
from labothappy.correlations.turbomachinery.radial_compressor_losses import radial_compressor_rotor_losses, radial_compressor_stator_losses

class CompressorRadialMeanLine(BaseComponent):
    """
    **Component**: Compressor

    **Model**: Constant isentropic efficiency

    **Descritpion**:

        This model determines the exhaust specific enthalpy and the exhaust temperature of a compressor. This model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.
        - Isentropic efficiency stays constant for all the conditions.

    **Connectors**:

        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

        W (WorkConnector): Work connector for the mechanical work.

    **Parameters**:

        eta_is: Isentropic efficiency. [-]

    **Inputs**:

        P_su: Suction side pressure. [Pa]

        T_su: Suction side temperature. [K]

        P_ex: Exhaust side pressure. [Pa]

        fluid: Working fluid. [-]

        m_dot: Mass flow rate of working fluid. [kg/s]

    **Ouputs**:

        h_ex: Exhaust side specific enthalpy. [J/kg] 

        T_ex: Exhaust side temperature. [K]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Mass_connector for the suction side
        self.ex = MassConnector() # Mass_connector for the exhaust side
        self.W = WorkConnector()

        self._STATE_KEYS = ('H', 'S', 'P', 'D', 'A', 'V')

        self.total_states  = {k: {i: np.nan for i in range(1, 6)} for k in self._STATE_KEYS}
        self.static_states  = {k: {i: np.nan for i in range(1, 6)} for k in self._STATE_KEYS}

        # Velocity Triangle Data
        self.Vel_Tri_R = {}
        self.Vel_Tri_S = {}
        

    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su', 'T_su', 'N_rot', 'fluid', 'm_dot']

    def get_required_parameters(self):
        return []

#%%

    def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.total_states['H'][position] = self.AS.hmass()            
        self.total_states['S'][position] = self.AS.smass()            
        self.total_states['P'][position] = self.AS.p()            
        self.total_states['D'][position] = self.AS.rhomass()            
    
        try:        
            self.total_states['A'][position] = self.AS.speed_sound()            
        except:
            self.total_states['A'][position] = -1  
            
        self.total_states['V'][position] = self.AS.viscosity()            
        
        return
    
    def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.static_states['H'][position] = self.AS.hmass()            
        self.static_states['S'][position] = self.AS.smass()            
        self.static_states['P'][position] = self.AS.p()            
        self.static_states['D'][position] = self.AS.rhomass()    
        
        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1            
            
        self.static_states['V'][position] = self.AS.viscosity()            
    
        return

#%%

    def solve_Rotor(self):
        
        # 1) Rotor inlet
        
        self.impeller_inlet_blockage = (
            1 - self.params['n_blade_R'] * self.params['t_b']
            / (np.pi * self.params['r1'] * 2 * np.sin(np.pi * abs(self.params['xhi1']) / 180))
        )
        self.A1_th = np.pi * (self.params['r1s']**2 - self.params['r1h']**2) * self.impeller_inlet_blockage
        
        self.params['pitch1'] = pitch1 = 2*np.pi*self.params['r1'] / self.params['n_blade_R']

        self.A1 = self.params['b1'] * pitch1 * self.params['n_blade_R']
        self.o1 = pitch1 * np.cos(np.pi/180 * self.params['xhi1']) - self.params['t_b']
        
        self.Vel_Tri_R['alpha1'] = alpha1 = self.params['alpha1_des']*np.pi/180
        self.Vel_Tri_R['u1'] = u1 = self.inputs['N_rot'] * 2 * np.pi/60 * self.params['r1']
        self.u1h = self.inputs['N_rot'] * 2 * np.pi/60 * self.params['r1h'] 
        self.u1s = self.inputs['N_rot'] * 2 * np.pi/60 * self.params['r1s']
        self.Vel_Tri_R['u2'] = self.inputs['N_rot'] * 2 * np.pi/60 * self.params['r2']

        def compute_h1_new(h1): 
            self.update_static_AS(CP.HmassSmass_INPUTS, h1, self.total_states['S'][1], 1)
            
            self.Vel_Tri_R['vm1'] = vm1 = self.inputs['m_dot'] / (self.static_states['D'][1] * self.A1_th)
            self.Vel_Tri_R['vu1'] = vu1 = vm1*np.tan(self.Vel_Tri_R['alpha1'])
            self.Vel_Tri_R['v1']  = v1  = np.sqrt(vm1**2 + vu1**2)

            self.Vel_Tri_R['wu1'] = wu1 = vu1 - self.Vel_Tri_R['u1']
            self.Vel_Tri_R['w1']  = w1  = np.sqrt(vm1**2 + wu1**2)
            self.Vel_Tri_R['beta1'] = np.arccos(vm1/w1)
            
            self.beta1h = np.arctan(self.u1h / vm1)
            self.beta1s = np.arctan(self.u1s / vm1)
            
            self.w1s = vm1 / np.cos(self.beta1s)
            self.w1h = vm1 / np.cos(self.beta1h)
                        
            h1_new = self.total_states['H'][1] - v1**2 / 2
                        
            return h1_new

        def residual_h1(h1):
            return h1 - compute_h1_new(h1)
        
        
        h1_min = self.total_states['H'][1] * 0.98
        h1_max = self.total_states['H'][1]
        
        h1_solution = brentq(residual_h1, h1_min, h1_max, xtol=1e-6)

        # 2) Rotor outlet
        
        self.params['pitch2'] = pitch2 = 2*np.pi*self.params['r2'] / self.params['n_blade_R']
        self.A2_th = self.params['b2'] * pitch2 * self.params['n_blade_R']
        
        self.A2_th = (
            (2*np.pi*self.params['r2']*self.params['b2'] - self.params['n_blade_R']*self.params['b2']*self.params['t_b'])
            * np.cos(np.pi * abs(self.params['xhi2']) / 180)
        )
        
        # Slip Factor
        self.sigma = 1 - np.sqrt(np.cos(self.params['xhi2']*np.pi/180)) / self.params['n_blade_R']**0.7
        angle_star_deg = 19 + 0.2*(90 - self.params['xhi2'])
        sigma_star = np.sin(np.pi/180 * angle_star_deg)
        
        r1_r2_lim = (self.sigma - sigma_star) / (1 - sigma_star)
        
        if self.params['r1']/self.params['r2'] > r1_r2_lim:
            num  = (self.params['r1']/self.params['r2']) - r1_r2_lim**(np.sqrt((90 - self.inputs['xhi2'])/10))
            fact = 1 - (num / (1 - sigma_star))
            self.sigma = self.sigma * fact
        
        
        def system_rotor(vm2):

            self.Vel_Tri_R['vm2'] = vm2
        
            rho2 = self.inputs['m_dot'] / (vm2 * self.A2_th)
            
            self.Vel_Tri_R['vu2'] = vu2 = (
                self.sigma * self.Vel_Tri_R['u2']
                - vm2 * np.tan(self.params['xhi2'] * np.pi/180)
            )
            
            self.dh0 = vu2*self.Vel_Tri_R['u2'] - self.Vel_Tri_R['vu1']*self.Vel_Tri_R['u1']
            
            h02 = self.total_states['H'][1] + self.dh0
            
            self.Vel_Tri_R['v2'] = v2 = np.sqrt(vu2**2 + vm2**2)
            h2 = h02 - v2**2 / 2
            
            try:
                self.AS.update(CP.DmassHmass_INPUTS, h2, rho2)
                s2 = self.AS.smass()
            except:   
                try:
                    s2 = PropsSI('S', 'H', h2, 'D', rho2, self.AS.fluid_names()[0]) 
                except:
                    return 1e6

            self.update_static_AS(CP.HmassSmass_INPUTS, h2, s2, 2)
            self.update_total_AS(CP.HmassSmass_INPUTS, h02, self.static_states['S'][2], 2)

            self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][2], self.total_states['S'][1])
            h_is = self.AS.hmass()
            
            self.Vel_Tri_R['alpha2'] = alpha2 = np.arccos(vm2 / v2)
            
            L_z = self.params['L_z']
                
            self.Vel_Tri_R['wu2'] = wu2 = self.Vel_Tri_R['vu2'] - self.Vel_Tri_R['u2']
            self.Vel_Tri_R['w2']  = w2  = np.sqrt(wu2**2 + vm2**2)
            self.Vel_Tri_R['beta2'] = beta2 = np.arccos(vm2 / w2)
            
            self.rotor_losses = radial_compressor_rotor_losses(
                A1=self.A1, A1_th=self.A1_th, alpha2=alpha2,
                beta1=self.Vel_Tri_R['beta1'], beta1h=self.beta1h, beta1s=self.beta1s,
                beta2=self.Vel_Tri_R['beta2'], b2=self.params['b2'], C_df=0.004, C_fi=0.004,
                Dh0=self.dh0, eps_a=self.params['eps_imp'], eps_b=self.params['eps_bf_imp'],
                eps_r=self.params['eps_imp'], k_roughness=self.params['k_imp'], L_z=L_z,
                mdot=self.inputs['m_dot'], mu1=self.static_states['V'][1],
                mu2=self.static_states['V'][2], n_bl_r=self.params['n_blade_R'],
                rho1=self.static_states['D'][1], rho2=self.static_states['D'][2],
                r1h=self.params['r1h'], r1s=self.params['r1s'], r2=self.params['r2'],
                u2=self.Vel_Tri_R['u2'], vu2=vu2, v1m=self.Vel_Tri_R['vm1'], v2=v2,
                w1=self.Vel_Tri_R['w1'], w1_th=self.Vel_Tri_R['w1'], w1h=self.w1h,
                w1s=self.w1s, w2=w2, xhi1=self.params['xhi1']*np.pi/180,
                xhi2=self.params['xhi2']*np.pi/180,
            )
            
            h02_new = self.rotor_losses['tot'] + h_is
            h2_new  = h02_new - v2**2 / 2
            
            self.update_static_AS(CP.HmassSmass_INPUTS, h2_new, self.static_states['S'][2], 2)
            
            try:
                s2 = self.static_states['S'][2]
                self.update_total_AS(CP.HmassSmass_INPUTS, h02, s2, 2)
            except Exception:
                return 1e6
            
            res = (h2 - h2_new) / h2_new
            self.res_rotor_ex = res
            
            return res
        
        rho2_max_phys = 3.0 * self.static_states['D'][1]
        vm2_min = self.inputs['m_dot'] / (rho2_max_phys * self.A2_th)
        vm2_max = self.inputs['m_dot'] / (self.static_states['D'][1] * self.A2_th)
        
        res_min = system_rotor(vm2_min)
        res_max = system_rotor(vm2_max)
        
        if res_min * res_max > 0:
            raise ValueError()
        
        self.sol_rotor_ex = brentq(system_rotor, vm2_min, vm2_max, xtol=1e-6)
        
        if self.rotor_losses['tot'] > self.dh0:
            raise ValueError()
        
        "Compute constraint terms"
        
        self.M1s_rel = self.w1s / self.static_states['A'][1]
        self.M1_rel  = self.Vel_Tri_R['w1'] / self.static_states['A'][1]
        self.W2_W1s  = self.Vel_Tri_R['w2'] / self.w1s
        
        vu1 = self.Vel_Tri_R['vu1']
        vu2 = self.Vel_Tri_R['vu2']
        u2  = self.Vel_Tri_R['u2']
        self.DR = 1.0 - (vu2 + vu1) / (2.0 * u2)
        
        return

#%%

    def solve_Stator(self):

        "S3) Vaneless Space Exhaust"
        self.static_states['P'][3] = p3 = (
            self.static_states['P'][2]
            + self.params["CP"] * (self.total_states['P'][2] - self.static_states['P'][2])
        )
        self.total_states['H'][3] = h03 = self.total_states['H'][2]
        
        def vaneless_system(x):
        
            s3 = x
        
            self.update_total_AS(CP.HmassSmass_INPUTS, self.total_states['H'][3], s3, 3)
            self.update_static_AS(CP.PSmass_INPUTS, self.static_states['P'][3], s3, 3)
        
            p03 = self.total_states['P'][3]
            K   = (self.total_states['P'][2] - p03) / (self.total_states['P'][2] - self.static_states['P'][2])
            CP_id = self.params['CP'] + K 
        
            self.Vel_Tri_S['v3']    = v3  = np.sqrt(2*(self.total_states['H'][3] - self.static_states['H'][3]))
            self.Vel_Tri_S['alpha3'] = self.Vel_Tri_R['alpha2']
            self.Vel_Tri_S['vm3']   = vm3 = v3 * np.cos(self.Vel_Tri_S['alpha3'])
            self.Vel_Tri_S['vu3']   = vu3 = v3 * np.sin(self.Vel_Tri_S['alpha3'])
            self.Vel_Tri_S['u3']    = u3  = self.Vel_Tri_R['u2'] * self.params['r3']/self.params['r2']
            self.Vel_Tri_S['wu3']   = wu3 = vu3 - u3
            self.Vel_Tri_S['wm3']   = wm3 = vm3
            self.Vel_Tri_S['w3']    = np.sqrt(wm3**2 + wu3**2)
            self.Vel_Tri_S['beta3'] = beta3 = np.arctan(wu3 / wm3)
            
            res = (self.inputs['m_dot'] - (vm3 * self.static_states['D'][3] * self.A2_th) * (1-CP_id)**(-0.5))**2
            
            return res
            
        self.AS.update(CP.HmassP_INPUTS, h03, p3)
        s3_max = self.AS.smass()
        s3_min = self.static_states['S'][2]
        
        sol = minimize(vaneless_system, s3_min, method='L-BFGS-B', bounds=[(s3_min, s3_max)],
                       options={'ftol': 1e-8, 'gtol': 1e-8})
    
        "S4) Vaned Diffuser Inlet"
        s4  = self.static_states['S'][3]
        h04 = self.total_states['H'][3]
        
        self.update_total_AS(CP.HmassSmass_INPUTS, h04, s4, 4)
        self.Vel_Tri_R['u3'] = self.Vel_Tri_R['u2'] * self.params['r3']/self.params['r2']

        self.params['pitch3'] = pitch3 = 2*np.pi*self.params['r3'] / self.params['n_blade_R']
        self.A3 = pitch3 * self.params['b3'] * self.params['n_blade_R']
        
        def vaned_diffuser_inlet_system(x):
            h4 = x
            self.update_static_AS(CP.HmassSmass_INPUTS, h4, s4, 4)

            self.Vel_Tri_S['v4']  = v4  = np.sqrt(2*(self.total_states['H'][4] - h4))
            self.Vel_Tri_S['vu4'] = vu4 = self.Vel_Tri_S['vu3']
            self.Vel_Tri_S['vm4'] = vm4 = np.sqrt(max(v4**2 - vu4**2, 0))
            self.Vel_Tri_S['alpha4'] = alpha4 = np.arctan(vu4 / vm4)
            
            self.Vel_Tri_S['u4']   = u4  = self.Vel_Tri_R['u3']
            self.Vel_Tri_S['wu4']  = wu4 = vu4 - u4
            self.Vel_Tri_S['wm4']  = wm4 = vm4
            self.Vel_Tri_S['w4']   = np.sqrt(wm4**2 + wu4**2)
            self.Vel_Tri_S['beta4']= beta4 = np.arctan(wu4 / wm4)
            
            res = vm4 - self.inputs['m_dot'] / (self.static_states['D'][4] * self.A3)
            self.res_diff_in = res
            return res
        
        h4_min = self.static_states['H'][3] * 0.9
        h4_max = self.total_states['H'][4]
        
        self.sol_vaned_diff = brentq(vaned_diffuser_inlet_system, h4_min, h4_max, xtol=1e-6)
        
        self.o4    = pitch3 * np.cos(self.Vel_Tri_S['alpha4']) - self.params['t_b']
        self.A3_th = self.o4 * self.params['b3'] * self.params['n_blade_R']
        
        "S5) Vaned Diffuser Outlet"
        
        self.A5 = self.A3_th * self.params['b5']/self.params['b3'] * self.params['r5']/self.params['r3']
        
        def system_stator(x):
            
            v5  = x[0]
            p5  = x[1]
        
            self.Vel_Tri_S['alpha5'] = alpha5 = self.params['xhi5'] * np.pi / 180
            vu5 = v5 * np.sin(alpha5)
            vm5 = v5 * np.cos(alpha5)
        
            self.Vel_Tri_S['vm5'] = vm5
            self.Vel_Tri_S['v5']  = v5
            self.Vel_Tri_S['vu5'] = vu5
        
            h05 = self.total_states['H'][4]
            h5  = h05 - v5**2 / 2
        
            self.stator_losses = radial_compressor_stator_losses(
                A4=self.A3, A4_th=self.A3_th, beta4=self.Vel_Tri_S['beta4'], C_f=0.004,
                r3=self.params['r3'], r4=self.params['r3'], r5=self.params['r5'],
                vm=vm5, w4=self.Vel_Tri_S['w4'], xhi3=self.Vel_Tri_S['alpha3'],
                xhi4=self.params['xhi4'], xhi5=self.params['xhi5'],
            )
        
            # Isentropic reference now uses the FREE p5, not p_ex
            self.AS.update(CP.PSmass_INPUTS, p5, self.static_states['S'][4])
            h5is = self.AS.hmass()
        
            h5_new = h5is + self.stator_losses['tot']
        
            self.update_static_AS(CP.HmassP_INPUTS, h5_new, p5, 5)
            rho5 = self.AS.rhomass()
        
            self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)
        
            # res1 : enthalpy consistency (same as before)
            res1 = (h5 - h5_new) / h05
        
            # res2 : continuity at diffuser exit — A5 and mdot are known geometry/operating point
            # mdot = rho5 * vm5 * A5
            res2 = (rho5 * vm5 * self.A5 - self.inputs['m_dot']) / self.inputs['m_dot']
        
            return [res1, res2]
        
        from scipy.optimize import least_squares

        v5_guess     = self.inputs['m_dot'] / (self.static_states['D'][4] * self.A5)
        v5_upper    = self.inputs['m_dot'] / (self.static_states['D'][4] * self.A5)

        p5_guess = 1.1*self.inputs['P_su']

        p5_lower    = 1.01*self.inputs['P_su']
        p5_upper    = 1.5*self.inputs['P_su']

        self.sol_sys_stator = least_squares(
            system_stator,
            x0     = [v5_guess, p5_guess],
            bounds = ([1.0, p5_lower], [v5_upper, p5_upper]),
            method = 'trf',
        )
        
        v5_sol, p5_sol = self.sol_sys_stator.x
        
        return

#%%

    def solve(self):
        
        # 0) Setup the solving
        
        self.check_calculable()
        self.check_parametrized()
        # self.print_setup()
        
        if not self.calculable:
            self.solved = False
            print("CompressorRadialMeanLine could not be solved. It is not calculable.")
            return
            
        if not self.parametrized:
            self.solved = False
            print("CompressorRadialMeanLine could not be solved. It is not parametrized.")
            return
        
        # 1) Setup the inlet

        self.AS = CP.AbstractState('HEOS', self.su.fluid)
        self.update_total_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
        
        # 2) Rotor 
        self.solve_Rotor()
        
        # 3) Stator
        self.solve_Stator()

        # Total-static efficiency
        self.AS.update(CP.PSmass_INPUTS, self.static_states['P'][5], self.total_states['S'][1])
        hout_is = self.AS.hmass()
        self.eta_is = (hout_is - self.total_states['H'][1]) / \
                      (self.static_states['H'][5] - self.total_states['H'][1])
        
        # Total-total efficiency
        self.AS.update(CP.PSmass_INPUTS, self.total_states['P'][5], self.total_states['S'][1])
        hout_is = self.AS.hmass()
        self.eta_is_tt = (hout_is - self.total_states['H'][1]) / \
                         (self.total_states['H'][5] - self.total_states['H'][1])
               
        return 
    
    def update_connectors(self, h_ex, w, W_dot):
        
        self.ex.reset()
        
        self.ex.set_h(h_ex)
        self.ex.set_p(self.inputs['P_ex'])
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        
        self.W.set_w(w)
        self.W.set_W_dot(W_dot)

    def print_results(self):
        print("=== Compressor Results ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, h={self.su.h}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, h={self.ex.h}")

        print("\nResults:")
        print(f"  - h_ex: {self.ex.h}")
        print(f"  - T_ex: {self.ex.T}")
        print("=========================")

    def print_states_connectors(self):
        print("=== Compressor Results ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")

if __name__ == "__main__":
    
    Comp = CompressorRadialMeanLine()

    Comp.set_inputs(
        fluid = 'CO2',
        m_dot  = 1.5*2.15,
        T_su = 305.97,
        P_su = 76.9*1e5,
        N_rot  = 50000,
    )

    Comp.set_parameters(
        alpha1_des = 1.063  *180/np.pi,
        n_blade_R = 15, 
        t_b = 0.762*1e-3,
        
        b1 = 0.00684,
        b2 = 0.00216,
        b3 = 0.00216,        
        b5 = 0.00216,
        
        CP           = 0.44,

        eps_imp      = 0.254*1e-3,
        eps_bf_imp   = 0.254*1e-3,
        k_imp        = 0.01*1e-3,
                
        L_z = 0.1137,
        
        r1 = 0.0077,
        r1s = 0.011,
        r1h = 0.0042,
        r2 = 0.0217, 
        r3 = 0.0396,
        r5 = 0.0494,
        
        xhi1 = 46.89,
        xhi2 = 43.44,
        xhi4 = 1.433    *180/np.pi,
        xhi5 = 0.00036  *180/np.pi,
        )

    Comp.solve()
