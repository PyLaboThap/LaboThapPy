# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2026

@author: titouanjanod

PCM-TES with an aluminium Kelvin-cell lattice (module of the CEA Grenoble), 1D lumped model.
Stepwise version of titouanjanod_perso/WP3/dynamic_modelling/lattice_TES_duCEA.py (1D 'lumped' part only).

Three classes, one per level of the dynamic architecture:
    TES_step       : the component. Computes ONE time step and keeps the current state h_PCM.
                     solve() computes the step without changing the state, advance() accepts it.
    TES_lattice    : the memory of a TES_step. Stores the history of all the steps (heat maps, Q(t)...),
                     saves / loads it and makes all the plots in one command (plot_all).
    TES_simulation : standalone simulation of the TES alone: constant HTF inlet with optional bypass (simulate,
                     i.e. the old solve(Dt_ch_dis, ...)) or HTF supply m_dot(t), T(t) imposed, e.g. measured
                     (simulate_measured), and the sizing of the TES.
                     It is NOT the generic integrator of a whole dynamic cycle.

Use (see component/examples/storage/TES_lattice_CEA_example.py):
    SIM = TES_simulation()
    SIM.set_parameters(T_PCM=..., cp_sol=..., ..., n_segment_box=..., n_box_serie=..., n_pass_box=..., n_box_parallel=...)
    history = SIM.simulate(Dt_ch_dis, tau, m_dot_tot, T_HTF_su, p_HTF, T_PCM_eq_init, fluid)
    history.plot_all()
"""

import __init__

from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

from component.base_component import BaseComponent

from correlations.convection.pipe_htc import gnielinski_pipe_htc, horizontal_tube_internal_condensation

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.lines import Line2D
from pathlib import Path
from datetime import datetime
from time import perf_counter

debug = False


# =============================================================================
# Plot style (Doumont-like: no box, direct labels, few ticks, title on top)
# =============================================================================
line_width = 3            # main curves
line_width_thin = 1.5     # reference / secondary lines
axis_line_width = 1.5     # axes spines and ticks
marker_size_exp = 40      # experimental points (scatter)
skip_exp = 60             # only one experimental point out of skip_exp is plotted (too many points)
fontsize_label = 18       # axis labels
fontsize_tick = 15        # tick labels
fontsize_annot = 16       # direct labels of curves
fontsize_title = 18       # figure titles
fig_size = (9, 6)
fig_size_heat_map = (10, 4.5)
dy_stack = 0.075          # vertical spacing of stacked labels [axes fraction]
white_mix = (0.55, 0.0)   # fraction of white for the lightest / darkest shade [-]

colors = {"Q_dot": "#1f77b4", "Q_TES": "#d62728",
          "HTF_in": "red", "HTF_out": "#ff7f0e", "HTF_mix": "#2ca02c", "m_dot": "#6BAED6",
          "PCM_0": "indigo", "PCM_50": "mediumorchid", "PCM_100": "pink",
          "cp": "red", "h": "blue", "T_PCM": "#E07A5F"}
color_key = "#3A3A3A"     # neutral text (key of the line / marker styles)


def new_axes():
    fig, ax = plt.subplots(figsize=fig_size)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(axis_line_width)
    ax.tick_params(labelsize=fontsize_tick, width=axis_line_width, length=6, direction="out")
    return fig, ax


def color_y1(ax, color):
    """y1 axis (spine, ticks, label) in the colour of its curve, when there is a y2 axis."""
    ax.spines["left"].set_color(color)
    ax.tick_params(axis="y", colors=color)
    ax.yaxis.label.set_color(color)


def add_y2(ax, color):
    """Secondary y axis on the right, spine / ticks / label in the colour of its curve."""
    ax2 = ax.twinx()
    for side in ("top", "left", "bottom"):
        ax2.spines[side].set_visible(False)
    ax2.spines["right"].set_linewidth(axis_line_width)
    ax2.spines["right"].set_color(color)
    ax2.tick_params(labelsize=fontsize_tick, width=axis_line_width, length=6, direction="out",
                    colors=color)
    ax2.yaxis.label.set_color(color)
    return ax2


def label_at_end(ax, x, y, text, color, dy_pt=0.0):
    """Direct label at the last finite point of a curve (instead of a legend)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.any():
        return
    ax.annotate(text, xy=(x[ok][-1], y[ok][-1]), xytext=(8, dy_pt), textcoords="offset points",
                color=color, fontsize=fontsize_annot, va="center", ha="left", annotation_clip=False)


def label_stack(ax, items, x=0.97, y_top=0.95, ha="right"):
    """
    Coloured labels stacked in axes coordinates, in the same vertical order as the curves
    (used when curve ends are too close for label_at_end). items: [(text, color), ...] top -> bottom.
    """
    for k, (text, color) in enumerate(items):
        ax.text(x, y_top - k * dy_stack, text, transform=ax.transAxes, color=color,
                fontsize=fontsize_annot, ha=ha, va="top")


def model_exp_key(ax, x, y_top, ha="right"):
    """
    Key of the two styles used on the figures, drawn instead of written: a line = model, a cross = exp.
    No frame: the key is in neutral dark colour (valid for every curve), the data are coloured.
    (x, y_top): top corner of the frame in axes coordinates, top right corner if ha = "right", top left if "left".
    On a figure with a y2 axis, give the top axis (twinx drawn last), otherwise the curves are drawn over the frame.
    """
    handles = [Line2D([], [], color=color_key, lw=line_width),
               Line2D([], [], color=color_key, marker="x", ls="none", markersize=np.sqrt(marker_size_exp),
                      markeredgewidth=1.5)]
    leg = ax.legend(handles, ["model", "exp."], loc="upper right" if ha == "right" else "upper left",
                    bbox_to_anchor=(x, y_top), bbox_transform=ax.transAxes, borderaxespad=0,
                    fontsize=fontsize_annot, labelcolor=color_key, handlelength=1.8, handletextpad=0.6,
                    frameon=False)
    leg.set_zorder(10)
    return leg


def shades(color, n):
    """n shades of one colour, light -> dark (mixed with white)."""
    c = np.array(mc.to_rgb(color))
    return [tuple(w + (1 - w) * c) for w in np.linspace(*white_mix, n)]


def create_result_folder(root_path, name_simu):

    folder_name = datetime.now().strftime(f"%Y_%m_%d-%H_%M_%S simulation - {name_simu}")

    result_folder = Path(root_path) / folder_name

    result_folder.mkdir(parents=True, exist_ok=True)

    return result_folder


# =============================================================================
# Level 1: one time step
# =============================================================================
class TES_step(BaseComponent):
    """
    **Component**: PCM-TES with an aluminium Kelvin-cell lattice cast around 16 tubes (module of the CEA Grenoble)

    **Model**: 1D lumped, enthalpy method, ONE time step (explicit Euler)

    **Description**:

        The TES is made of boxes (n_box_serie in series, n_box_parallel in parallel), each box has n_tubes_box tubes
        connected in n_pass_box passes. Only one tube branch is simulated (all branches are identical), discretised
        in Nx = n_segment_box * n_box_serie meshes along the HTF flow.
        In each mesh, the PCM and the aluminium lattice are lumped in one equivalent material (one temperature).
        The heat exchanged between the HTF and the equivalent material is Q_dot = UA * LMTD.

        The state of the TES (h_PCM_state, one value per mesh) is kept in the object.
        solve() computes the step t_n -> t_n+1 WITHOUT changing the state (it can be called several times,
        e.g. by an integrator that iterates on a cycle), advance() accepts the new state.

    **Assumptions**:

        - Equivalent material (PCM + lattice) per mesh, no radial gradient, conductivity lambda_eq
        - Latent heat taken into account with a gaussian apparent cp -> one h_PCM(T_PCM) curve
        - Charge or discharge decided in each mesh (HTF hotter or colder than the PCM)
        - HTF: no thermal inertia, no pressure drop, no heat loss to the ambient
        - Explicit Euler: the HTF profile of the step is computed with the PCM state at t_n

    **Connectors**:

        su (MassConnector): HTF supply.

        ex (MassConnector): HTF exhaust.

        Q (HeatConnector): heat given by the HTF to the whole TES [W] (> 0 in charge, < 0 in discharge).

    **Parameters**:

        T_PCM: Melting temperature of the PCM [K]

        cp_sol, cp_liq: Specific heat of the solid / liquid PCM [J/(kg.K)]

        Dh_latent: Latent heat of the PCM [J/kg]

        rho_PCM_sol, rho_PCM_liq: Density of the solid / liquid PCM [kg/m^3]

        lambda_eq: Conductivity of the equivalent material PCM + lattice [W/(m.K)]

        n_segment_box: Number of meshes along one box [-]

        n_box_serie, n_box_parallel: Number of boxes in series / in parallel [-]

        n_pass_box: Number of passes in one box [-]

        dt: Time step [s]

        (optional) Geometry of the box (De_tube, thk_tube, L_tube, V_TES, volume fractions...), see __init__,
        htc_correlation ('Gnielinski' or 'Hausen'),
        n_subsegment_cell (number of HTF sub-segments per cell in solve_h_HTF_ex, default 5)

    **Inputs**:

        fluid: HTF [-]

        h_su: HTF supply specific enthalpy [J/kg]

        P_su: HTF supply pressure [Pa]

        m_dot: HTF mass flow rate entering the whole TES [kg/s]

    **Outputs**:

        h_ex, P_ex: HTF exhaust specific enthalpy [J/kg] and pressure [Pa] (= P_su)

        Q_dot_tot: Heat given by the HTF to the TES [W]

        h_PCM_new, T_PCM_new: State of the equivalent material at t_n+1 [J/kg], [K] (accepted by advance())

        h_HTF_nodes, T_HTF_nodes: HTF at the inlet of each mesh and at the outlet of the last one (Nx+1 values)

    **Use**:

        set_parameters(...) -> set_initial_state(T_PCM_eq_init) -> at each step: set_inputs(...) -> solve() -> advance()
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector() # HTF supply
        self.ex = MassConnector() # HTF exhaust
        self.Q = HeatConnector()  # heat given by the HTF to the TES [W]

        self.wf = None
        self.AS = None

        self.AS_on = True

        self.sign_mode = None
        self.htc_correlation = 'Gnielinski' # 'Gnielinski' or 'Hausen'

        self.rho_alu = 2698 # kg/m^3
        self.cp_alu = 930 # J/kg/K
        self.lambda_alu = 160 # W/(m*K)
        self.V_foam_mesh = 2.89*1e-4
        self.thk_alu = 3e-3 # m

        self.lambda_PCM = None # W/(m*K)
        self.V_PCM_mesh = 1.609*1e-3

        self.rho_tube = 7800 # kg/m^3
        self.cp_tube = 500 # J/kg/K
        self.lambda_tube = 45 # W/(m*K) # 
        self.thk_tube = 2e-3 # m
        self.De_tube = 17.2e-3 # m
        self.L_tube = 1 #m
        # number of tubes in a standard box
        self.n_tubes_box = 16

        self.V_TES = 0.25 #m^3

        self.T_PCM = None # K
        self.cp_liq = None # J/kg/K
        self.Dh_latent = None # J/kg
        self.cp_sol = None # J/kg/K
        self.thk_PCM = 58.9e-3 # m
        self.rho_PCM_sol = None # kg/m^3
        self.rho_PCM_liq = None # kg/m^3

        self.V_frac_PCM = 78.7
        self.V_frac_shell = 4.8
        self.V_frac_alu = 14.9
        self.V_frac_tube = 0.7
        self.V_frac_HTF = 0.9

        self.lambda_eq = None # W/m·K

        # discretisation of the TES
        self.n_segment_box = None # number of meshes along one box
        self.n_box_serie = None
        self.n_pass_box = None
        self.n_box_parallel = None
        self.n_subsegment_cell = 5 # number of HTF sub-segments per cell (RK2 in solve_h_HTF_ex)

        self.dt = None # time step [s]

        # state of the TES (one value per mesh) and state at t_n+1 computed by solve()
        self.h_PCM_state = None # J/kg
        self.T_PCM_state = None # K
        self.h_PCM_new = None # J/kg
        self.T_PCM_new = None # K
        self.t = 0 # s


    def get_required_inputs(self):
        return ['fluid', 'h_su', 'P_su', 'm_dot']

    def get_required_parameters(self):
        return ['T_PCM', 'cp_sol', 'cp_liq', 'Dh_latent', 'rho_PCM_sol', 'rho_PCM_liq', 'lambda_eq',
                'n_segment_box', 'n_box_serie', 'n_pass_box', 'n_box_parallel', 'dt']


    def set_parameters(self, **kwargs):
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise ValueError(f"Unknown parameter '{key}'")
            setattr(self, key, value)
        super().set_parameters(**kwargs)

        # derived parameters, once the PCM and the discretisation are known (dt not needed)
        if all(self.params.get(param) is not None for param in self.get_required_parameters() if param != 'dt'):
            self.update_derived_parameters()


    def update_derived_parameters(self):

        self.Di_tube = self.De_tube - 2*self.thk_tube
        self.De_alu = self.De_tube + 2*self.thk_alu

        self.V_non_porosity = self.L_tube * self.n_tubes_box * self.De_alu**2/4 * np.pi
        self.rho_PCM_av = (self.rho_PCM_liq + self.rho_PCM_sol)/2

        self.V_frac_wo_shell = self.V_frac_HTF + self.V_frac_tube + self.V_frac_alu + self.V_frac_PCM

        self.V_PCM = self.V_TES / self.V_frac_wo_shell * self.V_frac_PCM
        self.V_alu = self.V_TES / self.V_frac_wo_shell * self.V_frac_alu
        self.V_tube = self.V_TES / self.V_frac_wo_shell * self.V_frac_tube
        self.V_HTF = self.V_TES / self.V_frac_wo_shell * self.V_frac_HTF

        self.V_alu_coating = self.L_tube * self.n_tubes_box * (self.De_alu**2/4 - self.De_tube**2/4) * np.pi
        self.V_foam = self.V_alu - self.V_alu_coating

        self.m_PCM = self.V_PCM * self.rho_PCM_av
        self.m_foam = self.V_foam * self.rho_alu

        self.A_eq = np.pi * ((self.De_alu/2+self.thk_PCM)**2 - (self.De_alu/2)**2)

        self.V_eq = self.V_PCM + self.V_foam
        self.m_eq = self.m_PCM + self.m_foam

        self.rho_eq = (self.rho_alu * self.V_foam + self.rho_PCM_av * self.V_PCM) / (self.V_foam + self.V_PCM)

        # discretisation (was in the old solve)
        self.Nx = math.floor(self.n_segment_box * self.n_box_serie)
        self.L_mesh = self.L_tube * self.n_pass_box / self.n_segment_box
        # number of tube branches in parallel in the whole TES (only one is simulated)
        self.n_tubes_parallel = self.n_tubes_box/self.n_pass_box * self.n_box_parallel
        m_PCM_branch = (self.m_PCM+self.m_foam)*(self.n_pass_box/self.n_tubes_box)*self.n_box_serie
        self.m_PCM_mesh = m_PCM_branch / self.Nx

        # h_PCM(T_PCM) curve of the equivalent material, computed once
        self.T_profile, self.cp_profile, self.h_profile = self.enthalpy_profile_creation()


    def U_eq(self, L_mesh, htf_htc):

        r_tube_i = self.Di_tube / 2
        r_tube_o = r_tube_i + self.thk_tube
        r_alu_o = r_tube_o + self.thk_alu
        r_pcm_o = r_alu_o + self.thk_PCM

        A_i = 2 * np.pi * r_tube_i * L_mesh

        R_conv = 1 / (htf_htc * A_i)

        R_tube = np.log(r_tube_o / r_tube_i) / (2 * np.pi * self.lambda_tube * L_mesh)

        R_alu = np.log(r_alu_o / r_tube_o) / (2 * np.pi * self.lambda_alu * L_mesh)

        R_pcm = np.log(r_pcm_o / r_alu_o) / (2 * np.pi * self.lambda_eq * L_mesh)

        R_tot = R_conv + R_tube + R_alu + R_pcm/2

        UA = 1 / R_tot

        Ueq = UA / A_i

        return Ueq, UA, R_tot, R_conv, R_tube, R_alu, R_pcm


    def enthalpy_profile_creation(self):

        # width of the latent T gaussian
        DT_gaus = 12

        # steps for the sensible and latent parts
        DT_sens = 0.1
        DT_lat = 0.01
        # Total width of the PCM profile
        DT_tot = 350


        T_solid = np.arange(self.T_PCM-(DT_tot/2),self.T_PCM-(DT_gaus/2), DT_sens)

        T_melting = np.arange(self.T_PCM-(DT_gaus/2), self.T_PCM+(DT_gaus/2), DT_lat)

        T_liquid = np.arange(self.T_PCM+(DT_gaus/2), self.T_PCM+(DT_tot/2), DT_sens)

        T_profile = np.concatenate((T_solid, T_melting, T_liquid))

        # cp of the PCM without the latent heat: cp_sol below the melting zone, cp_liq above,
        # and a linear transition from cp_sol to cp_liq over the melting zone (width DT_gaus)
        cp_background = np.piecewise(T_profile,
            [T_profile < self.T_PCM-(DT_gaus/2), (T_profile >= self.T_PCM-(DT_gaus/2))
                & (T_profile <= self.T_PCM+(DT_gaus/2)), T_profile > self.T_PCM+(DT_gaus/2)],
            [self.cp_sol,lambda T:self.cp_sol + (self.cp_liq-self.cp_sol) *(T-(self.T_PCM-(DT_gaus/2)))/DT_gaus,self.cp_liq])

        sigma = (DT_gaus/2)/3

        A = self.Dh_latent / (sigma*np.sqrt(2*np.pi))

        cp_latent = A*np.exp(-(T_profile-self.T_PCM)**2/(2*sigma**2))

        cp_PCM = cp_background + cp_latent

        # equivalent material: mass-weighted mean of the PCM and the aluminium lattice
        cp_profile = (self.m_PCM * cp_PCM + self.m_foam * self.cp_alu) / (self.m_PCM + self.m_foam)

        # h(T) = integral of cp dT (h = 0 at the lowest T of the profile)
        h_profile = cumulative_trapezoid(cp_profile, T_profile, initial=0)

        return T_profile, cp_profile, h_profile


    def interpolate_T_PCM(self, h_PCM_eq, h_profile, T_profile):
        return np.interp(h_PCM_eq, h_profile, T_profile)

    def interpolate_h_PCM(self, T_PCM_eq, h_profile, T_profile):
        return np.interp(T_PCM_eq, T_profile, h_profile)



    def HTC_calculation_Hausen(self, h_HTF_bulk, T_wall, p_HTF, m_dot_tube, L_mesh):


        self.AS.update(CP.HmassP_INPUTS, h_HTF_bulk, p_HTF)
        mu_f = self.AS.viscosity()
        cp_f = self.AS.cpmass()
        lambda_f = self.AS.conductivity()
        phase_bulk = self.AS.phase()

        self.AS.update(CP.PT_INPUTS, p_HTF, T_wall)
        mu_w = self.AS.viscosity()

        # wall on the other side of the saturation (subcooled boiling or wet-wall desuperheating, not modelled):
        # saturated viscosity on the bulk side
        if self.AS.phase() != phase_bulk:
            if phase_bulk == CP.iphase_liquid:
                self.AS.update(CP.PQ_INPUTS, p_HTF, 0)
                mu_w = self.AS.viscosity()
            elif phase_bulk == CP.iphase_gas:
                self.AS.update(CP.PQ_INPUTS, p_HTF, 1)
                mu_w = self.AS.viscosity()

        Re = 4*m_dot_tube/(np.pi*self.Di_tube*mu_f)
        Pr = cp_f*mu_f/lambda_f

        Nu = 0.116*(1+(self.Di_tube/L_mesh)**(2/3))*(Re**(2/3)-125)*Pr**(1/3)*(mu_f/mu_w)**0.14

        if debug:
            print(f"Re = {Re}")
            print(f"Nu={Nu}")

        htf_htc = Nu*lambda_f/self.Di_tube

        return htf_htc



    def HTC_calculation_Gnielinski(self, h_HTF_bulk, T_wall, p_HTF, m_dot_tube, L_mesh):


        self.AS.update(CP.HmassP_INPUTS, h_HTF_bulk, p_HTF)
        mu_f = self.AS.viscosity()
        cp_f = self.AS.cpmass()
        lambda_f = self.AS.conductivity()
        phase_bulk = self.AS.phase()

        self.AS.update(CP.PT_INPUTS, p_HTF, T_wall)
        mu_w = self.AS.viscosity()

        # wall on the other side of the saturation (subcooled boiling or wet-wall desuperheating, not modelled):
        # saturated viscosity on the bulk side
        if self.AS.phase() != phase_bulk:
            if phase_bulk == CP.iphase_liquid:
                self.AS.update(CP.PQ_INPUTS, p_HTF, 0)
                mu_w = self.AS.viscosity()
            elif phase_bulk == CP.iphase_gas:
                self.AS.update(CP.PQ_INPUTS, p_HTF, 1)
                mu_w = self.AS.viscosity()

        Pr = cp_f*mu_f/lambda_f

        S_i = math.pi * self.Di_tube**2/4
        G = m_dot_tube / S_i
        alpha = gnielinski_pipe_htc(mu_f, Pr, mu_w, lambda_f, G, self.Di_tube, L_mesh)[0]

        return alpha



    def solve_h_HTF_ex(self, h_HTF_in, T_PCM_eq, p_HTF, m_dot_tube, L_mesh, n_subsegment_cell=5):
        """ To calculate the outlet enthalpy of a cell, by marching on n_subsegment_cell sub-segments (midpoint method, RK2) """

        self.AS.update(CP.HmassP_INPUTS, h_HTF_in, p_HTF)
        T_HTF_in = self.AS.T()

        DT1 = self.sign_mode*(T_HTF_in - T_PCM_eq)

        if DT1 < 1e-1:
            return T_HTF_in, h_HTF_in

        L_subsegment = L_mesh / n_subsegment_cell
        A_subsegment = np.pi * self.Di_tube * L_subsegment

        # the HTF cannot go beyond T_PCM_eq
        self.AS.update(CP.PT_INPUTS, p_HTF, T_PCM_eq)
        h_HTF_PCM = self.AS.hmass()

        def Q_dot_subsegment(h_HTF, T_wall):

            self.AS.update(CP.HmassP_INPUTS, h_HTF, p_HTF)
            T_HTF = self.AS.T()
            is_two_phase = (self.AS.phase() == CP.iphase_twophase)
            if is_two_phase:
                x_HTF = self.AS.Qmass()

            for _ in range(30):

                # L_mesh for the entrance effect, as in solve_T_HTF_ex_LMTD
                if is_two_phase:
                    alpha = horizontal_tube_internal_condensation(self.wf, m_dot_tube, p_HTF, x_HTF, T_wall, self.Di_tube)
                elif self.htc_correlation == 'Hausen':
                    alpha = self.HTC_calculation_Hausen(h_HTF, T_wall, p_HTF, m_dot_tube, L_mesh)
                else:
                    alpha = self.HTC_calculation_Gnielinski(h_HTF, T_wall, p_HTF, m_dot_tube, L_mesh)

                UA_subsegment = self.U_eq(L_subsegment, alpha)[1]

                # short sub-segment: local DT ~ LMTD of the sub-segment
                Q_dot = UA_subsegment*(T_HTF - T_PCM_eq)

                T_wall_new = T_HTF - Q_dot/(alpha*A_subsegment)

                if abs(T_wall_new - T_wall) < 1e-3:
                    T_wall = T_wall_new
                    break

                T_wall = T_wall_new

            return Q_dot, T_wall

        h = h_HTF_in
        T_wall = T_PCM_eq

        for k in range(n_subsegment_cell):

            # prediction of the middle of the sub-segment with the inlet DT
            Q_dot_in, T_wall = Q_dot_subsegment(h, T_wall)
            h_mid = h - (Q_dot_in/2)/m_dot_tube

            # Q_dot with the DT in the middle of the sub-segment
            Q_dot_mid, T_wall = Q_dot_subsegment(h_mid, T_wall)
            h = h - Q_dot_mid/m_dot_tube

            if self.sign_mode*(h - h_HTF_PCM) <= 0:
                h = h_HTF_PCM
                break

        h_HTF_ex = h
        self.AS.update(CP.HmassP_INPUTS, h_HTF_ex, p_HTF)
        T_HTF_ex = self.AS.T()

        return T_HTF_ex, h_HTF_ex



    def set_initial_state(self, T_PCM_eq_init):
        """
        Initial temperature of the equivalent material (PCM + lattice) [K]:
        one value (uniform) or one value per mesh (Nx values).
        """
        self.T_PCM_state = np.ones(self.Nx) * T_PCM_eq_init
        self.h_PCM_state = self.interpolate_h_PCM(self.T_PCM_state, self.h_profile, self.T_profile)
        self.t = 0


    def solve(self):
        
        """Step t_n -> t_n+1. Does NOT change h_PCM_state: call advance() to accept the new state."""
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False
            print("TES_step could not be solved. It is not calculable and/or not parametrized")
            return

        if self.h_PCM_state is None:
            raise ValueError("TES_step: initial state not set. Use set_initial_state(T_PCM_eq_init) first.")

        if self.AS is None or self.wf != self.su.fluid:
            self.wf = self.su.fluid
            self.AS = CP.AbstractState("BICUBIC&HEOS", self.wf)

        Nx = self.Nx
        p_HTF = self.su.p # no pressure drop
        m_dot_tube = self.su.m_dot / self.n_tubes_parallel

        # HTF at the inlet of each mesh, and at the outlet of the last mesh (= outlet of the TES)
        h_HTF = np.zeros(Nx+1)
        T_HTF = np.zeros(Nx+1)
        self.Q_dot_local = np.zeros(Nx)
        self.h_PCM_new = np.zeros(Nx)
        self.T_PCM_new = np.zeros(Nx)

        h_HTF[0] = self.su.h
        self.AS.update(CP.HmassP_INPUTS, h_HTF[0], p_HTF)
        T_HTF[0] = self.AS.T()

        for i in range(Nx):

            if m_dot_tube <= 0:
                # no flow (e.g. everything in the bypass): no heat exchanged
                T_HTF_ex_mesh = T_HTF[i]
                h_HTF_ex_mesh = h_HTF[i]
                self.Q_dot_local[i] = 0

            else:
                # charge (+1) if the HTF is hotter than the PCM of this mesh, discharge (-1) otherwise
                if T_HTF[i] >= self.T_PCM_state[i]:
                    self.sign_mode = 1
                else:
                    self.sign_mode = -1


                T_HTF_ex_mesh, h_HTF_ex_mesh = self.solve_h_HTF_ex(h_HTF[i],
                                        self.T_PCM_state[i], p_HTF, m_dot_tube, self.L_mesh, self.n_subsegment_cell)
                
                # positive when charging
                # Q_dot_local is the Q_dot on only one mesh (of one tube branch)
                self.Q_dot_local[i] = m_dot_tube * (h_HTF[i] - h_HTF_ex_mesh)

            self.h_PCM_new[i] = self.h_PCM_state[i] + self.Q_dot_local[i]*self.dt/self.m_PCM_mesh
            self.T_PCM_new[i] = self.interpolate_T_PCM(self.h_PCM_new[i], self.h_profile, self.T_profile)

            h_HTF[i+1] = h_HTF_ex_mesh
            T_HTF[i+1] = T_HTF_ex_mesh

        self.h_HTF_nodes = h_HTF
        self.T_HTF_nodes = T_HTF
        # heat given by the HTF to the whole TES [W]
        self.Q_dot_tot = np.sum(self.Q_dot_local) * self.n_tubes_parallel

        self.update_connectors()
        self.solved = True


    def update_connectors(self):

        "Mass connectors"
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_p(self.su.p) # no pressure drop
        self.ex.set_h(self.h_HTF_nodes[-1])

        "Heat connector"
        self.Q.set_Q_dot(self.Q_dot_tot)


    def advance(self):
        """Accept the state computed by solve(): t_n -> t_n+1."""
        if not self.solved:
            raise ValueError("TES_step: solve() must be called before advance()")

        self.h_PCM_state = self.h_PCM_new.copy()
        self.T_PCM_state = self.T_PCM_new.copy()
        self.t = self.t + self.dt
        self.solved = False


    def print_results(self):
        print("=== TES Results ===")
        print(f"  - t: {self.t} [s]")
        print(f"  - Q_dot_tot: {self.Q_dot_tot} [W]")
        print(f"  - T_HTF_ex: {self.T_HTF_nodes[-1]} [K]")
        print(f"  - T_PCM (first, last mesh): {self.T_PCM_state[0]}, {self.T_PCM_state[-1]} [K]")
        print("=========================")

    def print_states_connectors(self):
        print("=== TES States ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], m_dot={self.ex.m_dot} [kg/s]")
        print(f"  - Q_dot: {self.Q.Q_dot} [W]")
        print("=========================")


# =============================================================================
# Level 2: memory of a TES_step (history, save / load, plots)
# =============================================================================
class TES_lattice:
    """
    Memory of a TES_step: history of all the time steps, save / load, plots (plot_all).
    The TES_step itself only knows its current state.

    Time convention (explicit Euler):
        row 0          : initial state. The HTF is filled with the inlet (no physical meaning),
                         so the HTF and the heat maps are plotted from the 2nd point.
        row n+1        : PCM at t_n+1, and HTF profile computed during the step n -> n+1 (with the PCM at t_n)
        Q_dot_local[n] : heat of the step n -> n+1 (last row = 0)

    After to_arrays():
        time [s], T_PCM_eq / h_PCM_eq (Nt x Nx), T_HTF / h_HTF (Nt x Nx+1: inlet of each mesh + outlet),
        T_HTF_ex, p_HTF, m_dot_in_TES, T_HTF_mix_out, Q_dot_local (one tube branch), Q_dot_tot [W], Q_tot, Q_cum [J]
    """

    def __init__(self, TES):
        self.TES = TES
        self.mode = None # 'charge' or 'discharge', for the titles (deduced from Q_tot if None)
        self.title = None # title of the figures; if None: built from the HTF inlet of the 1st step (see fig_title)
        self.bypass_on = False
        self.measured_inputs = False # True: HTF supply m_dot(t), T(t) imposed (e.g. measured), not constant
        self.skip_exp = skip_exp # only one experimental point out of skip_exp is plotted
        self.solving_time = None # s
        self.reset()


    def reset(self):
        self.time_list = []
        self.h_PCM_list = []
        self.T_PCM_list = []
        self.h_HTF_list = []
        self.T_HTF_list = []
        self.p_HTF_list = []
        self.m_dot_list = []
        self.T_HTF_mix_out_list = []
        self.Q_dot_local_list = []


    def record_initial_state(self):
        """To call after TES.set_initial_state() and TES.set_inputs() (HTF inlet at t_0)."""
        TES = self.TES
        self.reset()

        T_HTF_su = PropsSI("T", "H", TES.su.h, "P", TES.su.p, TES.su.fluid)

        self.time_list.append(TES.t)
        self.h_PCM_list.append(TES.h_PCM_state.copy())
        self.T_PCM_list.append(TES.T_PCM_state.copy())
        # HTF not computed yet: filled with the inlet (as in the old solve)
        self.h_HTF_list.append(np.ones(TES.Nx+1) * TES.su.h)
        self.T_HTF_list.append(np.ones(TES.Nx+1) * T_HTF_su)
        self.p_HTF_list.append(TES.su.p)
        self.m_dot_list.append(TES.su.m_dot)
        self.T_HTF_mix_out_list.append(T_HTF_su)


    def record_step(self, T_HTF_mix_out=None):
        """
        To call after TES.solve() and TES.advance().
        T_HTF_mix_out: HTF temperature after the mixing with the bypass [K] (default: outlet of the TES).
        """
        TES = self.TES

        if T_HTF_mix_out is None:
            T_HTF_mix_out = TES.T_HTF_nodes[-1]

        self.time_list.append(TES.t)
        self.h_PCM_list.append(TES.h_PCM_state.copy())
        self.T_PCM_list.append(TES.T_PCM_state.copy())
        self.h_HTF_list.append(TES.h_HTF_nodes.copy())
        self.T_HTF_list.append(TES.T_HTF_nodes.copy())
        self.p_HTF_list.append(TES.su.p)
        self.m_dot_list.append(TES.su.m_dot)
        self.T_HTF_mix_out_list.append(T_HTF_mix_out)
        self.Q_dot_local_list.append(TES.Q_dot_local.copy())


    def to_arrays(self):
        """Lists of the recorded steps -> arrays (Nt rows), as returned by the old solve."""
        Nx = self.TES.Nx

        self.time = np.array(self.time_list)
        self.h_PCM_eq = np.array(self.h_PCM_list)
        self.T_PCM_eq = np.array(self.T_PCM_list)
        self.h_HTF = np.array(self.h_HTF_list)
        self.T_HTF = np.array(self.T_HTF_list)
        self.p_HTF = np.array(self.p_HTF_list)
        self.m_dot_in_TES = np.array(self.m_dot_list)
        self.T_HTF_mix_out = np.array(self.T_HTF_mix_out_list)
        # Q_dot_local[n] is the step n -> n+1: last row = 0
        self.Q_dot_local = np.vstack([np.array(self.Q_dot_local_list).reshape(-1, Nx), np.zeros((1, Nx))])

        self.compute_totals()
        return self


    def compute_totals(self):
        # HTF temperature at the outlet of all the boxes
        self.T_HTF_ex = self.T_HTF[:,-1]
        # Heat duty of one simulated branch
        Q_dot_branch = np.sum(self.Q_dot_local, axis=1)
        # Heat duty of all the storage depending on time
        self.Q_dot_tot = Q_dot_branch * self.TES.n_tubes_parallel
        # Total energy exchanged by the storage
        self.Q_tot = np.sum(self.Q_dot_tot) * self.TES.dt
        # Cumulated energy
        self.Q_cum = np.cumsum(self.Q_dot_tot) * self.TES.dt


    # -------------------------------------------------------------------------
    # Save / load
    # -------------------------------------------------------------------------
    def save_results(self, result_folder, name_simu=""):
        """results.npz + simulation_parameters-<name_simu>.txt in result_folder."""
        TES = self.TES
        result_folder = Path(result_folder)

        np.savez_compressed(result_folder / "results.npz", time=self.time, m_dot_in_TES=self.m_dot_in_TES,
            T_HTF=self.T_HTF, h_HTF=self.h_HTF, T_HTF_mix_out=self.T_HTF_mix_out, p_HTF=self.p_HTF,
            T_PCM_eq=self.T_PCM_eq, h_PCM_eq=self.h_PCM_eq,
            Q_dot_local=self.Q_dot_local, Q_dot_tot=self.Q_dot_tot, Q_tot=self.Q_tot, Q_cum=self.Q_cum,
            tau=TES.dt, Nx=TES.Nx, n_box_parallel=TES.n_box_parallel, n_box_serie=TES.n_box_serie,
            n_pass_box=TES.n_pass_box, m_dot_tot=self.m_dot_in_TES[0], measured_inputs=self.measured_inputs)

        with open(result_folder / f"simulation_parameters-{name_simu}.txt", "w") as f:

            f.write("Model type = lumped (1D)\n")
            f.write(f"tau = {TES.dt} s\n")
            f.write(f"Nt = {len(self.time)}\n")
            f.write(f"Nx = {TES.Nx}\n")

            f.write("\n")
            f.write(f"n_box_parallel = {TES.n_box_parallel}\n")
            f.write(f"n_box_serie = {TES.n_box_serie}\n")
            f.write(f"n_pass_box = {TES.n_pass_box}\n")
            f.write(f"n_segment_box = {TES.n_segment_box}\n")
            f.write(f"n_subsegment_cell = {TES.n_subsegment_cell}\n")

            f.write("\n")
            if self.measured_inputs:
                f.write("HTF supply = measured (imposed m_dot(t), T_su(t))\n")
                f.write(f"T_su_HTF = {np.min(self.T_HTF[1:,0])-273.15:.2f} to {np.max(self.T_HTF[1:,0])-273.15:.2f} °C\n")
                f.write(f"m_dot = {np.min(self.m_dot_in_TES):.6f} to {np.max(self.m_dot_in_TES):.6f} kg/s\n")
            else:
                f.write(f"T_su_HTF = {self.T_HTF[1,0]-273.15:.2f} °C\n")
                f.write(f"m_dot_tot = {self.m_dot_in_TES[0]:.6f} kg/s\n")
            f.write(f"bypass_on = {self.bypass_on}\n")

            f.write("\n")
            f.write(f"p_HTF_init = {self.p_HTF[0]/1e5:.2f} bar\n")
            f.write(f"T_PCM_eq_init = {self.T_PCM_eq[0,0]-273.15:.2f} °C\n")

            f.write("\n")
            f.write(f"T_PCM = {TES.T_PCM-273.15:.2f} °C, Dh_latent = {TES.Dh_latent/1e3:.1f} kJ/kg, "
                    f"cp_sol = {TES.cp_sol:.0f} J/(kg.K), cp_liq = {TES.cp_liq:.0f} J/(kg.K), lambda_eq = {TES.lambda_eq} W/(m.K)\n")
            f.write(f"htc_correlation = {TES.htc_correlation}\n")

            if self.solving_time is not None:
                f.write("\n")
                f.write(f"Solving time = {self.solving_time:.2f} s\n")


    def load_results(self, file_path):
        """
        Load a results.npz (new or old lattice_TES_duCEA format). The TES_step must have the same
        parameters as the simulation (Nx, h_profile...). Old files: T_HTF has Nx columns (inlets of the meshes).
        """
        TES = self.TES
        data = np.load(file_path, allow_pickle=True)

        self.time = data["time"]
        self.m_dot_in_TES = data["m_dot_in_TES"]
        self.T_HTF = data["T_HTF"]
        self.T_HTF_mix_out = data["T_HTF_mix_out"]
        self.p_HTF = data["p_HTF"]
        self.T_PCM_eq = data["T_PCM_eq"]
        if "h_PCM_eq" in data:
            self.h_PCM_eq = data["h_PCM_eq"]
        else:
            self.h_PCM_eq = TES.interpolate_h_PCM(self.T_PCM_eq, TES.h_profile, TES.T_profile)

        self.Q_dot_local = data["Q_dot_local"]
        self.Q_dot_tot = data["Q_dot_tot"]
        self.Q_tot = data["Q_tot"]
        self.Q_cum = data["Q_cum"]
        self.T_HTF_ex = self.T_HTF[:,-1]

        TES.set_parameters(dt=data["tau"].item())
        self.measured_inputs = bool(data["measured_inputs"]) if "measured_inputs" in data else False
        # variable m_dot without measured inputs = bypass
        self.bypass_on = (not self.measured_inputs) and (not np.allclose(self.m_dot_in_TES, self.m_dot_in_TES[0]))
        return self


    # -------------------------------------------------------------------------
    # Plots
    # -------------------------------------------------------------------------
    def get_mode(self):
        if self.mode is not None:
            return self.mode
        if self.Q_tot >= 0:
            return "charge"
        return "discharge"


    def fig_title(self, extra=None):
        """
        self.title if given (e.g. with measured inputs, whose 1st step is not representative),
        else 'Charge – T_HTF,in = 135 °C – 700 kg/h' (simulation: no date).
        """
        if self.title is not None:
            return " – ".join([self.title, extra]) if extra else self.title
        parts = [self.get_mode().capitalize(),
                 rf"$T_{{HTF,in}}$ = {round(self.T_HTF[1,0]-273.15, 1):g} °C",
                 f"{self.m_dot_in_TES[0]*3600:.0f} kg/h"]
        if extra:
            parts.append(extra)
        return " – ".join(parts)


    def exp_col(self, exp, key):
        """One column of the experimental data (dict or DataFrame), one point out of self.skip_exp. None if absent."""
        if exp is None or key not in exp:
            return None
        return np.asarray(exp[key], dtype=float)[::self.skip_exp]


    def index_PCM(self):
        """Mesh index at 0 %, 50 % and 100 % of the length."""
        Nx = self.TES.Nx
        return {"PCM_0": round(0.0*(Nx-1)), "PCM_50": round(0.5*(Nx-1)), "PCM_100": round(1.0*(Nx-1))}


    def exp_T_PCM_keys(self):
        """
        Experimental PCM temperatures matching the model at 0 % / 100 % of the length:
        in charge the inlet (0 %) is the hottest, so max <-> 0 % and min <-> 100 %; the opposite in discharge.
        """
        if self.get_mode() == "charge":
            return {"PCM_0": "T_PCM_max_C", "PCM_50": "T_PCM_avg_C", "PCM_100": "T_PCM_min_C"}
        return {"PCM_0": "T_PCM_min_C", "PCM_50": "T_PCM_avg_C", "PCM_100": "T_PCM_max_C"}


    def plot_Q_vs_t(self, exp=None):
        """Q_dot_TES [kW] and Q_TES [kWh] (y2) vs time. > 0 in charge and in discharge."""
        sign = 1 if self.get_mode() == "charge" else -1

        fig, ax1 = new_axes()
        ax2 = add_y2(ax1, colors["Q_TES"])

        ax1.plot(self.time[:-1]/3600, sign*self.Q_dot_tot[:-1]/1e3, color=colors["Q_dot"], lw=line_width)
        ax2.plot(self.time/3600, sign*self.Q_cum/3.6e6, color=colors["Q_TES"], lw=line_width)

        t_exp = self.exp_col(exp, "t_h")
        if t_exp is not None:
            Q_dot_key = "Q_dot_smooth_kW" if "Q_dot_smooth_kW" in exp else "Q_dot_kW"
            ax1.scatter(t_exp, self.exp_col(exp, Q_dot_key), color=colors["Q_dot"], s=marker_size_exp, marker="x")
            ax2.scatter(t_exp, self.exp_col(exp, "Q_TES_kWh"), color=colors["Q_TES"], s=marker_size_exp, marker="x")
            # on ax2 (drawn over ax1), same place as the old text key
            model_exp_key(ax2, 0.97, 0.55)

        ax1.set_xlabel("Time [h]", fontsize=fontsize_label)
        ax1.set_ylabel(r"$\dot{Q}_{TES}$ [kW]", fontsize=fontsize_label)
        ax2.set_ylabel(r"$Q_{TES}$ [kWh]", fontsize=fontsize_label)
        color_y1(ax1, colors["Q_dot"])
        ax1.set_xlim(left=0)
        ax1.set_ylim(bottom=0)
        ax2.set_ylim(bottom=0)
        ax1.set_title(self.fig_title(), fontsize=fontsize_title, pad=14)
        fig.tight_layout()
        return fig


    def plot_T_vs_t(self, exp=None):
        """
        HTF inlet / outlet and PCM at 0, 50, 100 % of the length vs time.
        With bypass: HTF after the mixing, and m_dot_in_TES [kg/h] on y2.
        Labels stacked on the right, in the order of the curves (mean temperature over the simulation).
        """
        idx = self.index_PCM()
        t_h = self.time/3600

        fig, ax = new_axes()

        # (label, colour, mean temperature) of each curve, to stack the labels in the same order
        items = []

        ax.plot(t_h, self.T_HTF[:,0]-273.15, color=colors["HTF_in"], lw=line_width)
        items.append(("HTF inlet", colors["HTF_in"], np.mean(self.T_HTF[:,0])))
        ax.plot(t_h[1:], self.T_HTF_ex[1:]-273.15, color=colors["HTF_out"], lw=line_width)
        items.append(("HTF outlet", colors["HTF_out"], np.mean(self.T_HTF_ex[1:])))
        if self.bypass_on:
            ax.plot(t_h[1:], self.T_HTF_mix_out[1:]-273.15, color=colors["HTF_mix"], lw=line_width)
            items.append(("HTF mix outlet", colors["HTF_mix"], np.mean(self.T_HTF_mix_out[1:])))

        for key, label in (("PCM_0", "PCM 0 % length"), ("PCM_50", "PCM 50 % length"), ("PCM_100", "PCM 100 % length")):
            ax.plot(t_h, self.T_PCM_eq[:,idx[key]]-273.15, color=colors[key], lw=line_width)
            items.append((label, colors[key], np.mean(self.T_PCM_eq[:,idx[key]])))

        t_exp = self.exp_col(exp, "t_h")
        if t_exp is not None:
            ax.scatter(t_exp, self.exp_col(exp, "T_HTF_ex_C"), color=colors["HTF_out"], s=marker_size_exp, marker="x")
            for key, col in self.exp_T_PCM_keys().items():
                if col in exp:
                    ax.scatter(t_exp, self.exp_col(exp, col), color=colors[key], s=marker_size_exp, marker="x")

        items.sort(key=lambda item: item[2], reverse=True)
        stack = [(label, color) for label, color, _ in items]

        x_stack = 1.03
        if self.bypass_on:
            ax2 = add_y2(ax, colors["m_dot"])
            ax2.plot(t_h[1:], self.m_dot_in_TES[1:]*3600, color=colors["m_dot"], lw=line_width)
            ax2.set_ylabel(r"$\dot{m}_{in,TES}$ [kg/h]", fontsize=fontsize_label)
            ax2.set_ylim(bottom=0)
            x_stack = 1.24

        label_stack(ax, stack, x=x_stack, y_top=1.0, ha="left")
        if t_exp is not None:
            # below the stacked labels, one empty line in between
            model_exp_key(ax, x_stack, 1.0 - (len(stack) + 1)*dy_stack, ha="left")

        ax.set_xlabel("Time [h]", fontsize=fontsize_label)
        ax.set_ylabel(r"$T$ [°C]", fontsize=fontsize_label)
        ax.set_xlim(left=0)
        ax.set_title(self.fig_title(), fontsize=fontsize_title, pad=14)
        fig.tight_layout()
        fig.subplots_adjust(right=0.62 if self.bypass_on else 0.72)
        return fig


    def plot_heat_map(self):
        """T_PCM and T_HTF vs position and time, from the 2nd point (explicit Euler: row 0 = initialisation)."""
        L_mesh = self.TES.L_mesh

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=fig_size_heat_map, constrained_layout=True)

        t_h = self.time[1:]/3600
        # PCM: centre of the meshes; HTF: nodes (inlet of each mesh + outlet of the last one) [m]
        x_PCM = (np.arange(self.T_PCM_eq.shape[1]) + 0.5) * L_mesh
        x_HTF = np.arange(self.T_HTF.shape[1]) * L_mesh

        T1 = self.T_PCM_eq[1:, :] - 273.15
        T2 = self.T_HTF[1:, :] - 273.15

        vmin = np.floor(min(np.min(T1), np.min(T2)))
        vmax = np.ceil(max(np.max(T1), np.max(T2)))

        levels = np.linspace(vmin, vmax, 40)

        c1 = ax1.contourf(x_PCM, t_h, T1, levels=levels, cmap='turbo')
        c2 = ax2.contourf(x_HTF, t_h, T2, levels=levels, cmap='turbo')

        cbar = fig.colorbar(c1, ax=[ax1, ax2])
        cbar.set_label('$T$ [°C]', fontsize=fontsize_label)
        cbar.set_ticks(np.round(np.linspace(vmin, vmax, 6)).astype(int))
        cbar.ax.tick_params(labelsize=fontsize_tick)

        for ax, name in ((ax1, "$T_{PCM}$"), (ax2, "$T_{HTF}$")):
            ax.set_xlabel('$x$ [m]', fontsize=fontsize_label)
            ax.set_ylabel('Time [h]', fontsize=fontsize_label)
            ax.tick_params(labelsize=fontsize_tick)
            ax.set_title(name, fontsize=fontsize_title)

        fig.suptitle(self.fig_title(), fontsize=fontsize_title)
        return fig


    def plot_enthalpy_profile(self):
        """cp(T) and h(T) of the equivalent material, and the PCM states at 50 % of the length."""
        TES = self.TES
        i_50 = self.index_PCM()["PCM_50"]

        # only the temperature range of the simulation (+- 10 K), not the 350 K of the profile
        T_low = min(np.min(self.T_PCM_eq), np.min(self.T_HTF)) - 10
        T_high = max(np.max(self.T_PCM_eq), np.max(self.T_HTF)) + 10
        T_high = max(T_high, TES.T_PCM + 10)
        T_low = min(T_low, TES.T_PCM - 10)
        mask = (TES.T_profile >= T_low) & (TES.T_profile <= T_high)

        fig, ax1 = new_axes()
        ax2 = add_y2(ax1, colors["h"])

        ax1.plot(TES.T_profile[mask]-273.15, TES.cp_profile[mask]/1000, color=colors["cp"], lw=line_width)
        ax2.plot(TES.T_profile[mask]-273.15, TES.h_profile[mask]/1000, color=colors["h"], lw=line_width)
        ax2.scatter(self.T_PCM_eq[:,i_50]-273.15, self.h_PCM_eq[:,i_50]/1e3, color=colors["PCM_50"], s=80, zorder=3)

        label_stack(ax1, [("PCM state, 50 % length", colors["PCM_50"])], x=0.03, y_top=0.95, ha="left")

        ax1.set_xlabel("$T$ [°C]", fontsize=fontsize_label)
        ax1.set_ylabel(r"$c_{p,eq}$ [kJ/(kg·K)]", fontsize=fontsize_label)
        ax2.set_ylabel(r"$h_{eq}$ [kJ/kg]", fontsize=fontsize_label)
        color_y1(ax1, colors["cp"])
        ax1.set_title(self.fig_title(), fontsize=fontsize_title, pad=14)
        fig.tight_layout()
        return fig


    def plot_Q_vs_Q_TES(self, with_T_PCM=False, exp=None):
        """
        Q_dot_TES vs Q_TES (> 0 in charge and in discharge). with_T_PCM: PCM at 0, 50, 100 % of the length on y2,
        in shades of one colour (light -> dark = 0 -> 100 %).
        """
        sign = 1 if self.get_mode() == "charge" else -1
        Q_TES = sign*self.Q_cum[:-1]/3.6e6

        fig, ax = new_axes()
        ax.plot(Q_TES, sign*self.Q_dot_tot[:-1]/1e3, color=colors["Q_dot"], lw=line_width)

        has_exp = exp is not None and "Q_TES_kWh" in exp
        if has_exp:
            Q_dot_key = "Q_dot_smooth_kW" if "Q_dot_smooth_kW" in exp else "Q_dot_kW"
            ax.scatter(self.exp_col(exp, "Q_TES_kWh"), self.exp_col(exp, Q_dot_key), color=colors["Q_dot"],
                       s=marker_size_exp, marker="x")

        ax.set_xlabel(r"$Q_{TES}$ [kWh]", fontsize=fontsize_label)
        ax.set_ylabel(r"$\dot{Q}_{TES}$ [kW]", fontsize=fontsize_label)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

        if with_T_PCM:
            idx = self.index_PCM()
            color_y1(ax, colors["Q_dot"])
            ax2 = add_y2(ax, colors["T_PCM"])
            cols = shades(colors["T_PCM"], 3)
            for (key, label), c in zip((("PCM_0", "0 %"), ("PCM_50", "50 %"), ("PCM_100", "100 %")), cols):
                ax2.plot(Q_TES, self.T_PCM_eq[:-1,idx[key]]-273.15, color=c, lw=line_width)
                if has_exp and self.exp_T_PCM_keys()[key] in exp:
                    ax2.scatter(self.exp_col(exp, "Q_TES_kWh"), self.exp_col(exp, self.exp_T_PCM_keys()[key]),
                                color=c, s=marker_size_exp, marker="x")
            ax2.set_ylabel(r"$T_{PCM}$ [°C]", fontsize=fontsize_label)
            # key of the shades, stacked on the right (the curves end at the same point: no label at the end)
            stack = [(r"$T_{PCM}$ at length:", cols[-1])] + [(label, c) for label, c in zip(("0 %", "50 %", "100 %"), cols)]
            label_stack(ax, stack, x=1.24, y_top=1.0, ha="left")
            if has_exp:
                # below the stacked labels, one empty line in between
                model_exp_key(ax2, 1.24, 1.0 - (len(stack) + 1)*dy_stack, ha="left")

        elif has_exp:
            model_exp_key(ax, 0.97, 0.95)

        ax.set_title(self.fig_title(), fontsize=fontsize_title, pad=14)
        fig.tight_layout()
        if with_T_PCM:
            fig.subplots_adjust(right=0.62)
        return fig


    def plot_all(self, save_dir=None, name_simu="", exp=None):
        """
        All the figures in one command. Returns a dict of figures, saved in SVG in save_dir if given.
        exp: experimental data (dict or DataFrame, e.g. compute_derived(load_CEA_test(name)) of plots_tests_CEA.py)
        with columns t_h, Q_dot_kW (or Q_dot_smooth_kW), Q_TES_kWh, T_HTF_ex_C, T_PCM_min_C, T_PCM_avg_C, T_PCM_max_C.
        """
        figs = {
            "Q_vs_t": self.plot_Q_vs_t(exp=exp),
            "T_vs_t": self.plot_T_vs_t(exp=exp),
            "heat_map": self.plot_heat_map(),
            "enthalpy_profile": self.plot_enthalpy_profile(),
            "Q_vs_Q_TES": self.plot_Q_vs_Q_TES(exp=exp),
            "Q_and_T_PCM_vs_Q_TES": self.plot_Q_vs_Q_TES(with_T_PCM=True, exp=exp),
        }

        if save_dir is not None:
            for key, fig in figs.items():
                fig.savefig(Path(save_dir) / f"{key}-{name_simu}.svg", bbox_inches="tight")

        return figs


# =============================================================================
# Level 3: standalone simulation of the TES alone + sizing
# =============================================================================
class TES_simulation:
    """
    Standalone simulation of the TES alone (no cycle):
        simulate          : constant HTF inlet during Dt_ch_dis, with an optional bypass that keeps the HTF after
                            the mixing at T_target. Replaces the old solve(Dt_ch_dis, ...).
        simulate_measured : HTF supply m_dot(t), T(t) imposed (e.g. measured in a test), no bypass.
    Also contains the sizing of the TES (number of boxes in series / in parallel).

    The TES_step is in self.TES, its memory (TES_lattice) in self.history.
    """

    def __init__(self):
        self.TES = TES_step()
        self.history = TES_lattice(self.TES)
        self.mode = None
        self.sign_mode = None
        self.print_flag = False


    def set_parameters(self, **kwargs):
        self.TES.set_parameters(**kwargs)


    def simulate(self, Dt_ch_dis, tau, m_dot_tot, T_HTF_su, p_HTF, T_PCM_eq_init, fluid, T_target=None, bypass_on=False):
        """
        Dt_ch_dis: duration of the charge / discharge [s], tau: time step [s], m_dot_tot: HTF flow rate [kg/s],
        T_HTF_su: HTF inlet temperature [K], p_HTF: HTF pressure [Pa], T_PCM_eq_init: initial PCM temperature [K],
        T_target: HTF temperature after the mixing with the bypass [K] (only if bypass_on).
        Returns the history (TES_lattice).
        """
        t0_perf = perf_counter()
        TES = self.TES

        Nt = math.floor(Dt_ch_dis / tau) + 1
        TES.set_parameters(dt=tau)
        TES.set_initial_state(T_PCM_eq_init)

        # charge if the HTF is hotter than the PCM at the beginning
        if T_HTF_su > T_PCM_eq_init:
            self.mode = "charge"
            self.sign_mode = 1
        else:
            self.mode = "discharge"
            self.sign_mode = -1

        AS = CP.AbstractState("BICUBIC&HEOS", fluid)
        AS.update(CP.PT_INPUTS, p_HTF, T_HTF_su)
        h_HTF_su = AS.hmass()

        if bypass_on:
            AS.update(CP.PT_INPUTS, p_HTF, T_target)
            h_target = AS.hmass()

        m_dot_in_TES = m_dot_tot
        TES.set_inputs(fluid=fluid, h_su=h_HTF_su, P_su=p_HTF, m_dot=m_dot_in_TES)

        self.history.mode = self.mode
        self.history.bypass_on = bypass_on
        self.history.measured_inputs = False
        self.history.record_initial_state()

        # outlet of the TES at t_0: not computed yet, taken equal to the inlet (as in the old solve)
        h_TES_out = h_HTF_su

        for n in range(Nt-1):

            if bypass_on:

                h_HTF_in = h_HTF_su

                if abs(h_HTF_in-h_TES_out) < 1:
                    m_dot_in_TES = m_dot_tot
                else:
                    m_dot_target = m_dot_tot*(h_HTF_in-h_target)/(h_HTF_in-h_TES_out)
                    m_dot_target = np.clip(m_dot_target,0,m_dot_tot)
                    alpha = 0.5
                    m_dot_in_TES = alpha*m_dot_in_TES + (1-alpha)*m_dot_target

            TES.set_inputs(fluid=fluid, h_su=h_HTF_su, P_su=p_HTF, m_dot=m_dot_in_TES)
            TES.solve()
            TES.advance()

            h_TES_out = TES.h_HTF_nodes[-1]

            # calculating the resulting outlet temperature if bypass
            if bypass_on:
                h_HTF_mix_out = ((m_dot_tot-m_dot_in_TES)*h_HTF_su + m_dot_in_TES*h_TES_out) / m_dot_tot
                AS.update(CP.HmassP_INPUTS, h_HTF_mix_out, p_HTF)
                T_HTF_mix_out = AS.T()
            else:
                T_HTF_mix_out = None

            self.history.record_step(T_HTF_mix_out)

            if self.print_flag:
                print(f"Q_dot_tot_local[{n}]={TES.Q_dot_tot/1000:.3f} kW")

        self.history.to_arrays()
        self.history.solving_time = perf_counter() - t0_perf

        return self.history


    def simulate_measured(self, time, m_dot, T_HTF_su, tau, p_HTF, T_PCM_eq_init, fluid, Dt_ch_dis=None):
        """
        Simulation of the TES alone with imposed HTF supply conditions that vary in time (e.g. measured in a test).
        No bypass: all the flow goes through the TES.

        time: time of the data [s], from the start of the simulation (t = 0), increasing
        m_dot: HTF flow rate entering the whole TES at each time [kg/s]
        T_HTF_su: HTF supply temperature at each time [K]
        tau: time step [s]
        p_HTF: HTF pressure [Pa], one value or one value per time
        T_PCM_eq_init: initial temperature of the PCM + lattice [K], uniform (one value) or one value per mesh
        Dt_ch_dis: duration of the simulation [s] (default: up to the last time of the data)

        The supply conditions at t_n are linearly interpolated in the data, and used for the step t_n -> t_n+1
        (explicit Euler, as in simulate). The mode (charge / discharge) of the titles is deduced from Q_tot.
        Returns the history (TES_lattice).
        """
        t0_perf = perf_counter()
        TES = self.TES

        time = np.asarray(time, dtype=float)
        m_dot = np.asarray(m_dot, dtype=float)
        T_HTF_su = np.asarray(T_HTF_su, dtype=float)
        p_HTF = np.broadcast_to(np.asarray(p_HTF, dtype=float), time.shape)

        if Dt_ch_dis is None:
            Dt_ch_dis = time[-1]
        Nt = math.floor(Dt_ch_dis / tau) + 1

        AS = CP.AbstractState("BICUBIC&HEOS", fluid)

        def supply_at(t):
            """HTF supply at time t: h [J/kg], p [Pa], m_dot [kg/s]"""
            p_su = np.interp(t, time, p_HTF)
            AS.update(CP.PT_INPUTS, p_su, np.interp(t, time, T_HTF_su))
            return AS.hmass(), p_su, np.interp(t, time, m_dot)

        TES.set_parameters(dt=tau)
        TES.set_initial_state(T_PCM_eq_init)

        h_su, p_su, m_dot_su = supply_at(0)
        TES.set_inputs(fluid=fluid, h_su=h_su, P_su=p_su, m_dot=m_dot_su)

        self.mode = None
        self.sign_mode = None
        self.history.mode = None
        self.history.bypass_on = False
        self.history.measured_inputs = True
        self.history.record_initial_state()

        for n in range(Nt-1):

            h_su, p_su, m_dot_su = supply_at(n*tau)
            TES.set_inputs(fluid=fluid, h_su=h_su, P_su=p_su, m_dot=m_dot_su)
            TES.solve()
            TES.advance()

            self.history.record_step()

            if self.print_flag:
                print(f"Q_dot_tot_local[{n}]={TES.Q_dot_tot/1000:.3f} kW")

        self.history.to_arrays()
        self.history.solving_time = perf_counter() - t0_perf

        return self.history



    def design_generation(self, Q_TES_predesign, T_HTF_min, T_HTF_max):
        # we consider the maximum and minimum PCM temperatures being the HTF inlet and oulet extreme temperatures
        # coming from the steady-state analysis
        TES = self.TES

        # calculating the theoretical energy stored in 1 box:
        T_profile, cp_profile, h_profile = TES.enthalpy_profile_creation()

        h_min = TES.interpolate_h_PCM(T_HTF_min, h_profile, T_profile)
        h_max = TES.interpolate_h_PCM(T_HTF_max, h_profile, T_profile)
        Dh = h_max - h_min
        Q_box = TES.m_eq * Dh # enthalpy difference multiplied by TES mass

        N_box_th = math.ceil(Q_TES_predesign/Q_box)

        # N_box_tests = np.round(N_box_th*np.array([0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5])).astype(int)
        N_box_tests = np.round(N_box_th*np.array([1.55,1.55,1.60,1.65,1.7,1.75])).astype(int)
        k_values = np.array([0.4,0.5,0.6,0.7])
        configs = set()

        for N_box in N_box_tests:
            for k in k_values:
                N_par = max(1,round(k*np.sqrt(N_box)))
                N_ser = math.ceil(N_box/N_par)
                N_box_real = N_par*N_ser

                configs.add((N_box_real,N_par,N_ser))

        configs = sorted(configs)
        N_cases = len(configs)

        N_box_mat = np.zeros(N_cases,dtype=int)
        N_par_mat = np.zeros(N_cases,dtype=int)
        N_ser_mat = np.zeros(N_cases,dtype=int)

        for i,(N_box,N_par,N_ser) in enumerate(configs):

            N_box_mat[i] = N_box
            N_par_mat[i] = N_par
            N_ser_mat[i] = N_ser

        return N_box_mat, N_par_mat, N_ser_mat



    def TES_sizing(self, T_target, Q_TES_predesign, T_HTF_min, T_HTF_max, Dt_ch_dis, tau, m_dot_tot, T_HTF_init, p_HTF_init,
               T_PCM_eq_init, fluid, bypass_on):
        """
        Smallest configuration (n_box_parallel, n_box_serie) whose HTF outlet reaches T_target at the end of Dt_ch_dis.
        The other parameters (n_segment_box, n_pass_box, PCM...) are the ones already set.
        """

        # génération des configurations
        N_box_test, N_par, N_series = self.design_generation(Q_TES_predesign, T_HTF_min, T_HTF_max)

        target_achieved = np.full(len(N_box_test), False)

        target_margin = np.zeros(len(N_box_test))

        for i in range(len(N_box_test)):


            n_parallel_test = N_par[i]
            n_serie_test = N_series[i]
            print("n_parallel_test =", n_parallel_test)
            print("n_serie_test =", n_serie_test)

            self.set_parameters(n_box_parallel=n_parallel_test, n_box_serie=n_serie_test)
            history = self.simulate(Dt_ch_dis, tau, m_dot_tot, T_HTF_init, p_HTF_init, T_PCM_eq_init, fluid,
                                    T_target, bypass_on)

            target_margin[i] = self.sign_mode*(T_target - history.T_HTF_ex[-1])

        target_achieved = target_margin > 0

        valid = np.where(target_achieved)[0]

        if len(valid) == 0:
            best_attempt = np.argmax(target_margin)
            return {
                "N_parallel": N_par[best_attempt],
                "N_series": N_series[best_attempt],
                "N_box": N_box_test[best_attempt],
                "margin": target_margin[best_attempt],
                "index": best_attempt,
                "target_reached": False
                }

        opti = valid[np.argmin(N_box_test[valid])]

        best_config = {
            "N_parallel": N_par[opti],
            "N_series": N_series[opti],
            "N_box": N_box_test[opti],
            "margin": target_margin[opti],
            "index": opti,
            "target_reached": True
            }

        return best_config
