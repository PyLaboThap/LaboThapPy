
from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import minimize, fsolve
import matplotlib.pyplot as plt

class axial_turbine_mean_line_design(object):

    def __init__(self):

        self.fluid = None
        self.m_dot = 0

        self.n_stage = 0
        self.W_dot_stage = [0]

        self.W_dot_req = 0
        self.p_end = 0
        self.Rp_tot = 0
        self.Rp_stage = 0

        self.R = [0]
        self.psi = [0]
        self.phi = [0]
        self.eta_s_tt = [0]

        self.omega = 0
        self.r = [0]
        self.AR = [0]
        self.h = [0]
        self.Re = 0

        self.v_m = [0]
        self.u = [0]

        self.v = {  '1': [0],
                    '2': [0],
                    '3': [0]}
        
        self.w = {'1': [0],
                  '2': [0],
                  '3': [0]}

        self.alpha = {'1': [0],
                      '2': [0],
                      '3': [0]}
        
        self.beta = {'1': [0],
                     '2': [0],
                     '3': [0]}

        self.p_tot = {'1': [0],
                     '2': [0],
                     '3': [0]}
        
        self.T_tot = {'1': [0],
                     '2': [0],
                     '3': [0]}
        
        self.h_tot = {'1': [0],
                     '2': [0],
                     '3': [0]}
        
        self.s =    {'1': [0],
                     '2': [0],
                     '3': [0]}      

    def angle_velocity_triangles_axial(self, R, psi, phi, alpha_1, blade_type):
        # Axial in/out flow
        # alpha_1 in rad
        
        # Inlet relative angle
        beta_1 = alpha_1 + np.atan(1/phi)

        # Outlet relative angle
        beta_2 = np.atan((1-2*R)/phi) - alpha_1

        # Outlet absolute angle
        if blade_type == 'Stator':
            alpha_2 = alpha_1 - np.atan(psi/phi)
        elif blade_type == 'Rotor':
            alpha_2 = alpha_1 + np.atan(psi/phi)
        else:
            print("'blade_type' shall be either 'Stator' or 'Rotor'")

        return beta_1, beta_2, alpha_2

    def axial_turbine_first_stage_opt(self, local = False):
        # Axial in/out flow considered
        self.n_stage = 1

        h0 = PropsSI('H','T', self.T_tot['1'][0],'P', self.p_tot['1'][0], self.fluid) # J/kg
        self.s['1'][0] = PropsSI('S','T', self.T_tot['1'][0],'P', self.p_tot['1'][0], self.fluid) # J/K
        self.h_tot['1'][0] = h0 # J/kg

        # ---- 1 -> 2 Stator blade --------------------------
        self.h_tot['2'][0] = self.h_tot['1'][0]
        self.p_tot['2'][0] = self.p_tot['1'][0]
        self.T_tot['2'][0] = self.T_tot['1'][0]
        self.s['2'][0] = self.s['1'][0]

        # Speed of sound 
        a_2 = PropsSI('A','H', self.h_tot['2'][0],'S', self.s['2'][0], self.fluid) # J/kg

        # Rotor inlet absolute speed 
        self.v['2'][0] = self.M_secu*a_2

        # Speed triangles
        beta_1, beta_2, alpha_2 = self.angle_velocity_triangles_axial(self.R, self.psi, self.phi, self.alpha['1'][0], 'Stator') # abs & rel angles

        self.beta['1'][0] = beta_1
        self.beta['2'][0] = beta_2
        self.alpha['2'][0] = alpha_2

        alpha_2, beta_3, alpha_3 = self.angle_velocity_triangles_axial(self.R, self.psi, self.phi, self.alpha['2'][0], 'Rotor') # abs & rel angles

        self.beta['2'][0] = beta_2
        self.beta['3'][0] = beta_3
        self.alpha['3'][0] = alpha_3

        self.v_m[0] = self.v['2'][0]*np.cos(self.alpha['2'][0]) # meridional speed
        self.w['2'][0] = self.v_m[0]/np.cos(self.beta['2'][0]) # relative speed

        self.v['1'][0] = self.v_m[0]/np.cos(self.alpha['1'][0]) # axial inlet flow : meridional speed is conserved accros the turbine

        w_u2 = np.sqrt(self.w['2'][0]**2 - self.v_m[0]**2) # relative speed projection in rotor speed direction
        v_u2 = np.sqrt(self.v['2'][0]**2 - self.v_m[0]**2) # absolute speed projection in rotor speed direction

        self.u[0] = v_u2 - w_u2 # rotor speed

        # ---- 2 -> 3 Rotor blade --------------------------
        self.v['3'][0] = self.v_m[0]/np.cos(self.alpha['3'][0])
        self.w['3'][0] = self.v_m[0]/np.cos(self.beta['3'][0])

        # Work production
        Dh_tot = self.psi*self.u[0]**2

        self.W_dot_stage[0] = Dh_tot*self.m_dot

        self.h_tot['3'][0] = self.h_tot['2'][0] - Dh_tot

        # ---- Computation of rotor radius and rotation speed --------------------------
        h1 = self.h_tot['1'][0] - self.v['1'][0]**2 /2
        h3 = self.h_tot['3'][0] - self.v['3'][0]**2 /2

        h3s = h1 + (h3 - h1)/self.eta_s_tt[0]
        h3s_tot = self.h_tot['1'][0] + (self.h_tot['3'][0] - self.h_tot['1'][0])/self.eta_s_tt[0]

        self.p_tot['3'][0] = PropsSI('P','H', h3s_tot,'S', self.s['2'][0], self.fluid)

        self.s['3'][0] = PropsSI('S', 'H', self.h_tot['3'][0], 'P', self.p_tot['3'][0], self.fluid)
        self.T_tot['3'][0] = PropsSI('T', 'H', self.h_tot['3'][0], 'P', self.p_tot['3'][0], self.fluid)

        # Static density and viscosity
        rho3 = PropsSI('D', 'H', h3, 'S', self.s['3'], self.fluid)
        mu3 = PropsSI('V', 'H', h3, 'S', self.s['3'], self.fluid)

        # Rotor blade cord
        Q = self.m_dot/rho3
        A = Q/self.v_m[0]
        c = (self.Re*mu3)/(rho3*self.v_m[0])

        # Rotor blade height
        self.h[0] = c*self.AR[0]
        
        # Rotor radius
        self.r[0] = A/(2*np.pi*self.h[0])

        # Hub diameter
        self.r_h = self.r[0] - self.h[0]

        # Rotation speed
        self.omega = self.u[0]/self.r[0]
        self.Rp_stage = self.p_tot['1'][0]/self.p_tot['3'][0]

        return self.W_dot_stage[0], self.Rp_stage

    def add_axial_turbine_stage(self, stage):
        # Axial in/out flow considered

        self.h_tot['1'].append(self.h_tot['3'][stage-1])  # J/kg
        self.T_tot['1'].append(self.T_tot['3'][stage-1])  # J/kg
        self.p_tot['1'].append(self.p_tot['3'][stage-1])  # J/kg
        self.s['1'].append(self.s['3'][stage-1]) # J/K

        # ---- 1 -> 2 Stator blade --------------------------
        self.h_tot['2'].append(self.h_tot['1'][stage])
        self.s['2'].append(self.s['1'][stage])

        # Speed of sound 
        a_2 = PropsSI('A','H', self.h_tot['2'][stage],'S', self.s['2'][stage], self.fluid) # J/kg

        # Rotor inlet absolute speed 
        self.v['2'].append(self.M_secu*a_2)

        # Speed triangles
        self.alpha['1'].append(self.alpha['3'][stage-1])

        beta_1, beta_2, alpha_2 = self.angle_velocity_triangles_axial(self.R, self.psi, self.phi, self.alpha['1'][stage], 'Stator') # abs & rel angles

        self.beta['1'].append(beta_1)
        self.beta['2'].append(beta_2)
        self.alpha['2'].append(alpha_2)

        beta_2, beta_3, alpha_3 = self.angle_velocity_triangles_axial(self.R, self.psi, self.phi, alpha_2, 'Rotor') # abs & rel angles

        self.beta['2'].append(beta_2)
        self.beta['3'].append(beta_3)
        self.alpha['3'].append(alpha_3)

        self.v_m.append(self.v['2'][stage]*np.cos(alpha_2)) # meridional speed
        self.w['2'].append(self.v_m[stage]/np.cos(beta_2)) # relative speed

        self.v['1'].append(self.v_m[stage]/np.cos(self.alpha['1'][stage])) # axial inlet flow : meridional speed is conserved accros the turbine

        w_u2 = np.sqrt(self.w['2'][stage]**2 - self.v_m[stage]**2) # relative speed projection in rotor speed direction
        v_u2 = np.sqrt(self.v['2'][stage]**2 - self.v_m[stage]**2) # absolute speed projection in rotor speed direction

        self.u.append(v_u2 - w_u2) # rotor speed

        # ---- 2 -> 3 Rotor blade --------------------------
        self.v['3'].append(self.v_m[stage]/np.cos(alpha_3))
        self.w['3'].append(self.v_m[stage]/np.cos(beta_3))

        # Work production
        Dh_tot = self.psi*self.u[stage]**2
        self.W_dot_stage.append(Dh_tot*self.m_dot)

        self.h_tot['3'].append(self.h_tot['2'][stage] - Dh_tot)

        # ---- Computation of rotor radius and rotation speed --------------------------
        self.p_tot['3'].append(self.p_tot['1'][stage]/self.Rp_stage)
        self.T_tot['3'].append(PropsSI('T', 'H', self.h_tot['3'][stage], 'P', self.p_tot['3'][stage], self.fluid))
        self.s['3'].append(PropsSI('S', 'H', self.h_tot['3'][stage], 'P', self.p_tot['3'][stage], self.fluid))

        h3s_tot = PropsSI('H', 'S', self.s['2'][stage], 'P',  self.p_tot['3'][stage], self.fluid)
        self.eta_s_tt.append((self.h_tot['2'][stage] - self.h_tot['3'][stage])/(self.h_tot['2'][stage] - h3s_tot))

        self.r.append(self.u[stage]/self.omega)

        # ---- Blade height and aspect ratio --------------------------
        self.h.append(self.r[stage] - self.r_h)

        # Static density and viscosity
        h3 = self.h_tot['3'][stage] - self.v['3'][stage]**2
        rho3 = PropsSI('D', 'H', h3, 'S', self.s['3'][stage], self.fluid)
        mu3 = PropsSI('V', 'H', h3, 'S', self.s['3'][stage], self.fluid)

        # Rotor blade cord
        c = (self.Re*mu3)/(rho3*self.v_m[stage])

        self.AR.append(self.h[stage]/c)

        return 

    def design(self, inputs, params):

        # Load inputs and parameters

        self.fluid = inputs['fluid']
        self.m_dot = inputs['m_dot']
        self.W_dot_req = inputs['W_dot']

        self.v['1'][0] = inputs['v0']
        self.T_tot['1'][0] = inputs['T0']
        self.p_tot['1'][0] = inputs['p0']
        self.alpha['1'][0]  = params['alpha_1']

        self.AR[0] = params['AR']

        self.R = params['R']
        self.psi = params['psi']
        self.phi = params['phi']
        self.eta_s_tt[0] = params['eta_s_tt']

        self.Rp_tot = inputs['p0']/inputs['p_end']
        self.M_secu = params['M_secu']
        self.Re = params['Re']

        def W_dot_res(vars):
            eta_s_tt_1, phi = vars
            self.eta_s_tt[0] = eta_s_tt_1
            self.phi = phi

            W_dot_stage, Rp_stage = self.axial_turbine_first_stage_opt()
            self.n_stage = np.log(self.Rp_tot)/np.log(Rp_stage)

            n_stage_round = round(self.n_stage,0)
            return abs(n_stage_round-self.n_stage)

        bounds = [(0.3, 0.8),
                  (0.4, 0.6)]     # y bounds

        result = minimize(W_dot_res, x0=[self.eta_s_tt[0], self.phi], bounds=bounds, method='SLSQP')
        
        self.n_stage = int(round(self.n_stage,0))

        for i in range(self.n_stage-1):
            self.add_axial_turbine_stage(i+1)

        self.W_dot_tot = sum(self.W_dot_stage)
        
        h_tot_s_out = PropsSI('H','S', self.s['1'][0], 'P', self.p_tot['3'][-1],self.fluid)
        h_tot_in = self.h_tot['1'][0]
        h_tot_out = self.h_tot['3'][-1]
        
        self.eta_s_t = (self.h_tot['1'][0]-self.h_tot['3'][-1])/(self.h_tot['1'][0]-h_tot_s_out)

        return

first_stage_params = {  'M_secu'  :     0.8,
                        'R'       :     0.5, # Degree of reaction 
                        'psi'     :       1, # Work coefficient
                        'phi'     :     0.5, # Flow coefficient
                        'Re'      :   1*1e6, # Required Reynolds number at blade outlet for turbulent flow
                        'AR'      :       3, # Blade aspect ratio 
                        'eta_s_tt':     0.8, # Stage static-to-static isentropic efficiency
                        'alpha_1' :       0  # Stator inlet angle [rad]
                    }

# ---- Cycle inputs --------------------------

first_stage_inputs = {'fluid':'Cyclopentane',
          'm_dot':         34.51, # kg/s
          'W_dot':      2506*1e3, # W
          'v0'   :             0, # m/s^2
          'T0'   :    273.15+131, # K
          'p0'   :     767.8*1e3, # Pa
          'p_end'   :        82*1e3 # Pa
          }

Sehrene_turbine = axial_turbine_mean_line_design()
Sehrene_turbine.design(first_stage_inputs, first_stage_params)

print("Design Results:")
print("-------------------------")
print(f"Turbine imposed pressure ratio: {Sehrene_turbine.Rp_tot} [-]")
print(f"Turbine number of stages: {Sehrene_turbine.n_stage} [-]")
print(f"Turbine stage pressure ratio: {Sehrene_turbine.Rp_stage} [-]")
print("-------------------------")
print(f"Turbine required power: {Sehrene_turbine.W_dot_req*1e-3} [kW]")
print(f"Turbine total power: {round(Sehrene_turbine.W_dot_tot*1e-3,1)} [kW]")
print("-------------------------")
print(f"Turbine required isentropic efficiency: {0.8} [-]")
print(f"Turbine design isentropic efficiency: {round(Sehrene_turbine.eta_s_t,2)} [-]")
print("-------------------------")
print(f"Turbine rotating speed: {Sehrene_turbine.omega*60/(2*np.pi)} [RPM]")
print(f"Maximum Mach number imposed between stator and rotors: {Sehrene_turbine.M_secu} [-]")
print("-------------------------")

stages_vec = np.linspace(1,Sehrene_turbine.n_stage, Sehrene_turbine.n_stage)

fontsize = 16

plt.figure()
plt.plot(stages_vec, Sehrene_turbine.r)
plt.ylabel("Rotor Radius [m]", fontsize = fontsize)
plt.axis([0.8,4.2, 0.14,0.24])
plt.grid()
plt.show()

plt.figure()
plt.plot(stages_vec, np.array(Sehrene_turbine.W_dot_stage)*1e-3)
plt.ylabel("Stage Work [kW]", fontsize = fontsize)
plt.axis([0.8,4.2, 600, 750])
plt.grid()
plt.show()

plt.figure()
plt.plot(stages_vec, np.array(Sehrene_turbine.h)*1e3)
plt.ylabel("Blade height [mm]", fontsize = fontsize)
plt.axis([0.8,4.2, 20, 60])
plt.grid()
plt.show()

plt.figure()
plt.plot(stages_vec, np.array(Sehrene_turbine.AR))
plt.ylabel("Blade Aspect Ratio [-]", fontsize = fontsize)
plt.axis([0.8,4.2, 0, 4])
plt.grid()
plt.show()
