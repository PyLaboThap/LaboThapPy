# External imports
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import brentq
import warnings
import matplotlib.pyplot as plt
from CoolProp.Plots import PropertyPlot

# Internal imports
from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.heat_connector import HeatConnector

from labothappy.correlations.heat_exchanger.f_lmtd2 import f_lmtd2, F_shell_and_tube
# from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.hex_MB_correlations import check_correlation_name
from labothappy.correlations.pressure_drop.pressure_drop_distributor import pressure_drop
from labothappy.correlations.convection.heat_transfer_distributor import compute_htc, NEEDS_Q_FLUX
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction
from labothappy.correlations.properties.two_phase import compute_two_phase_density

HTC_TYPES   = ('correlation', 'user_defined')
HTC_REGIMES = ('single-phase', 'two-phase', 'supercritical')                                 # Keys of htc_corr
HTC_PHASES  = ('liquid', 'vapor', 'two-phase', 'vapor-wet', 'dryout', 'supercritical')    # Keys of htc_user
PHASE_TO_REGIME = {'liquid': 'single-phase', 'vapor': 'single-phase',
                   'two-phase': 'two-phase', 'vapor-wet': 'two-phase',
                   'supercritical': 'supercritical'}


T_EX_MIN_H = 218.0            # Lowest outlet temperature allowed on the hot side [K]
T_EX_MAX_C = 273.15 + 481.0   # Highest outlet temperature allowed on the cold side [K]


class _Side:
    """Settings and state of one fluid stream ('H' = hot and 'C' = cold).
    
    Enthalpy-ordered vectors (hvec, pvec, Tvec, ...) go from low to high
    enthalpy on both sides : index 0 is the cold inlet / hot outlet end.
    """

    def __init__(self, name):
        self.name = name
        self.label = 'hot' if name == 'H' else 'cold'
        # settings
        self.htc_type = None     # 'Correlation' or 'User-Defined'
        self.htc_corr = {}       # {'1P' : name, '2P' : name, 'SC' : name} correlation name for each regime
        self.htc_user = {}       # {'Liquid': value, ...} fixed htc for each phase
        self.dp_corr = {}        # {'1P' : name, '2P' : name, 'SC' : name} correlation name for each regime
        self.dp_user = 0.0       # Fixed total pressure drop in Pa
        # CoolProp states (re-created only if fluid or backend changes)
        self._AS_key = None # Remembers which backends and fluid the object were built on
        self.AS = self.AS_sat = self.AS_heos = None # AS = main property object, A_sat = only used for saturation properties; A_heos = a backup of the full slower equation of state
        # pressure-drop profile: cumulative fraction of DP vs. normalised
        # enthalpy measured from the side inlet (linear by default)
        self.dp = 0.0 # Total pressure drop
        self.dp_s = np.array([0.0, 1.0]) # Positions along the side
        self.dp_f = np.array([0.0, 1.0]) # Fraction of total pressure drop reached at each of these positions

    @property
    def two_phase_possible(self):
        return not (self.incomp or self.supercritical)

    def correlation(self, phase, kind): 
        """Correlation name for a cell phase ('supecritical' falls back to 'single-phase')."""
        corr = self.htc_corr if kind == 'htc' else self.dp_corr
        regime = PHASE_TO_REGIME[phase]

        # Vapor-wet cell: condensing film at the wall (two-phase htc), but superheated bulk flow (single-phase dp)
        if kind == 'dp' and phase == 'vapor-wet':
            regime = 'single-phase'

        name = corr.get(regime)  # Looks up the correlation name for the regime in the dictionary ie martin_holger_plate_HTC

        if name is None and regime == 'supercritical':   # No supercritical correlation given: use the single-phase one
            name = corr.get('single-phase')

        if name is None:
            what = 'heat transfer' if kind == 'htc' else 'pressure drop'
            raise ValueError(f"No '{regime}' {what} correlation given for the {self.label} side (cell phase: {phase}).")
        return name

    

class HexMBPlate(BaseComponent):
    """
    **Component**: Heat Exchanger

    **Model**: The model is based on the the work of Ian Bell "A generalized moving-boundary algorithm to predict the heat transfer
                    rate of counterflow heat exchangers for any phase configuration". 
    
    **Descritpion**:
    
        It is a moving-boundary model allowing for the precise estimation of its performance, pressure drops and fluid charges inside it (using a void fraction correlation for two-phase flow conditions). 
        For now, the method has been adapted for many geometries : 
        - Brazed-Plate heat exchangers
        - Shell-and-Tube heat exchangers
        - Cross-CounterFlow tube and fin heat exchangers

        The precision of the model and the geometry dependent charge estimation requires to know the geometry. 
        The required parameters are adapted to the used heat transfer and pressure drop correlations that can be chosen from the ones implemented.

    **Connectors**:
    
        su_H (MassConnector): Mass connector for the hot suction side.
        su_C (MassConnector): Mass connector for the cold suction side.

        ex_H (MassConnector): Mass connector for the hot exhaust side.
        ex_C (MassConnector): Mass connector for the cold exhaust side.

        Q_dot (HeatConnector): Heat connector for the heat transfer between the fluids

    **Parameters**:

    **Inputs**:
        P_su_H: Hot suction side pressure. [Pa]
        
        h_su_H: Hot suction side enthalpy. [J/kg]

        fluid_H: Hot suction side fluid. [-]

        m_dot_H: Hot suction side mass flowrate. [kg/s]

        P_su_C: Cold suction side pressure. [Pa]

        h_su_C: Cold suction side enthalpy. [J/kg]

        fluid_C: Cold suction side fluid. [-]

        m_dot_C: Cold suction side mass flowrate. [kg/s]
    
    **Ouputs**:
    
        h_ex_H: Hot exhaust side specific enthalpy. [J/kg]

        P_ex_H: Hot exhaust side pressure. [Pa]

        h_ex_C: Cold exhaust side specific enthalpy. [J/kg]

        P_ex_C: Cold exhaust side pressure. [Pa]

        Q: Heat transfer rate [W]

        charge_H : Hot fluid charge [kg]

        charge_C : Cold fluid charge [kg]
    
    """

    HEX_TYPES = ('Plate', 'Shell&Tube', 'Tube&Fins', 'PCHE')

    def __init__(self, hex_type):

        super().__init__()

        if hex_type not in self.HEX_TYPES:
            raise ValueError(f"Heat exchanger types implemented for this model are: {self.HEX_TYPES}.")
        self.hex_type = hex_type

        self.su_H = MassConnector()
        self.su_C = MassConnector()
        self.ex_H = MassConnector()
        self.ex_C = MassConnector()
        self.Q = HeatConnector()

        self.w = [None] # Length fraction at each cell

        self.eval = 0 # Diagnostic tool to count the calls to the objective function
        self.H = _Side('H')
        self.C = _Side('C')

    def get_required_inputs(self):
        #List of required inputs
        return['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']


    def get_required_parameters(self):
        """Names of the parameters the model needs."""
        general = ['flow_type',  # Type of flow (counter-current, co-current, ...)
                'htc_type',   # Type of heat transfer coefficient (user, correlations, ...)
                'dp_type',    # Type of pressure drop model
                'n_disc',     # Number of discretizations [-]
                'A_H', 'A_C', # Heat transfer area on the hot and cold side [m^2]
                'V_H', 'V_C', # Internal volume on the hot and cold side [m^3]
                ]

        return general + list(self.get_default_geometry().keys())


    def get_default_geometry(self):
        """Default values for the geometrical parameters."""
        if self.hex_type == 'Plate':
            return {
                # Neutral defaults: 0 means "effect neglected"
                'R_fouling': 0.0,  # Fouling factor [m^2K/W]
                'roughness': 0.0,  # Surface roughness [m]

                # No meaningful default: None means "not given"
                'k_wall': None,                 # Wall thermal conductivity [W/m/K]
                # 't_plate': None,                # Plate thickness [m]    -> COMPUTED!!
                'n_plates': None,            # Number of plates [-]
                'n_channels_H': None,        # Number of channels on the hot side [-]
                'n_channels_C': None,        # Number of channels on the cold side [-]
                'L_plate': None,             # Plate length [m]
                'W_plate': None,             # Plate width [m]
                'L_eff': None,               # Effective length of the plate [m]
                'plate_pitch': None,         # Plate pitch [m]
                'corrugation_amplitude': None,   # Corrugation amplitude [m]
                'corrugation_pitch': None,   # Corrugation pitch [m]
                'chevron_angle': None,       # Chevron angle [rad]
                # 'enlargement_factor': None,  # Enlargement factor [-]    -> COMPUTED!!!
                # 'A_cs_H': None,              # Cross-flow section of one channel, hot side [m^2]    -> COMPUTED !!!
                # 'A_cs_C': None,              # Cross-flow section of one channel, cold side [m^2]   -> COMPUTED !!!
                # 'D_h_H': None,               # Hydraulic diameter of one channel, hot side [m]      -> COMPUTED !!!
                # 'D_h_C': None,               # Hydraulic diameter of one channel, cold side [m]     -> COMPUTED !!!
            }
        else:
            return {}  # TODO: check from Basile's code


    def set_htc(self, htc_type="correlation", htc_corr_h=None, htc_corr_c=None, htc_user_h=None, htc_user_c=None):
        """
        htc_type : "user_defined" or "correlation"
        htc_corr_h, htc_corr_c : Correlation name for each regime on each side
        htc_user_h, htc_user_c : Fixed heat transfer coeffcient for eahc phase on each side
        """

        if htc_type not in ("correlation", "user_defined"):
            raise ValueError(f"htc_type must be 'correlation' or 'user_defined', get {htc_type!r}.")
        self.params['htc_type'] = htc_type

        for side, htc_corr, htc_user in ((self.H, htc_corr_h, htc_user_h), (self.C, htc_corr_c, htc_user_c)):
            side.htc_type = htc_type
            arg = f"htc_{'user' if htc_type == 'user_defined' else 'corr'}_{side.name.lower()}"  # Argument name, for error message

            if htc_type == 'user_defined':
                if htc_user is None:
                    raise ValueError(f"{arg} is required for the {side.label} side.")
                unknown = set(htc_user) - set(HTC_PHASES) # Checks if correctly written by user (liquid, vapor,...)
                if unknown:
                    raise ValueError(f"Unknown phase(s) {sorted(unknown)} in {arg}. Allowed: {HTC_PHASES}.")
                side.htc_user = dict(htc_user)

            else:
                if htc_corr is None or 'single-phase' not in htc_corr:
                    raise ValueError(f"{arg} with at least a 'single-phase' correlation is required for the {side.label} side.")
                unknown = set(htc_corr) - set(HTC_REGIMES)
                if unknown:
                    raise ValueError(f"Unknown regime(s) {sorted(unknown)} in {arg}. Allowed: {HTC_REGIMES}.")
                side.htc_corr = {regime: name for regime, name in htc_corr.items() if name is not None}  # Remove regimes given as None

    def set_dp(self, dp_type=None, dp_corr_h=None, dp_corr_c=None, dp_user_h=None, dp_user_c=None):
        """
        dp_type: None = No pressure drops, "user_defined" or "correlation_global" or "correlation_disc"

        """
        allowed = (None, "user_defined", "correlation_global", "correlation_disc")
        if dp_type not in allowed:
            raise ValueError(f"dp_type must be one of {allowed}.")

        self.params['dp_type'] = dp_type
        for side, corr, ud in ((self.H, dp_corr_h, dp_user_h), (self.C, dp_corr_c, dp_user_c)):
            side.dp_user = float(ud or 0.0)
            side.dp_corr = {}
            if dp_type in ("correlation_global", "correlation_disc"):
                side.dp_corr = {k: v for k, v in dict(corr).items() if v is not None}


    def _setup_geometry(self):
        """
        Side geometry computations
        """
        p = self.params
        uses_correlations = (p.get('htc_type') == 'correlation'
                            or p.get('dp_type') == 'correlation')
        
        if not uses_correlations:
            for side, su in ((self.H, self.su_H), (self.C, self.su_C)):
                s = side.name.upper()   # 'H' or 'C', to build the parameter names
                side.geom = {**p,
                            'A': p[f'A_{s}'],    # Heat transfer area [m^2]
                            'V': p[f'V_{s}'],    # Internal volume [m^3]
                            }
            return  # geometry not needed for user-defined htc and dp
    
        if self.hex_type == "Plate":
            # enlargement_factor = 
            p = self.params
            for side, su in ((self.H, self.su_H), (self.C, self.su_C)):
                s = side.name.upper()   # 'H' or 'C', to build the parameter names
                A_cs_channel = p['W_plate'] * p['corrugation_amplitude']
                enlargement_factor = 1/6 * (1+ np.sqrt(1 + (np.pi * p['corrugation_amplitude'] / p['corrugation_pitch'])**2) + 4 * np.sqrt(1 + (np.pi * p['corrugation_amplitude'] / p['corrugation_pitch'])**2) / 2)
                D_h = 2*p['corrugation_amplitude']/enlargement_factor
                side.geom = {**p,
                            'A': p[f'A_{s}'],    # Heat transfer area [m^2]
                            'V': p[f'V_{s}'],    # Internal volume [m^3]
                            'L_plate': p['L_plate'],    # Plate length [m]
                            'W_plate': p['W_plate'],    # Plate width [m]
                            'L_eff': p['L_eff'],    # Effective length of the plate [m]
                            'n_channels': p[f'n_channels_{s}'],    # Number of channels [-]
                            'plate_pitch': p['plate_pitch'],    # Plate pitch [m]
                            'corrugation_amplitude': p['corrugation_amplitude'],    # Corrugation amplitude [m]
                            'corrugation_pitch': p['corrugation_pitch'],    # Corrugation pitch [m]
                            'chevron_angle': p['chevron_angle'],    # Chevron angle [rad]
                            'A_cs_channel': A_cs_channel,    # Cross-flow section of one channel [m^2]
                            'D_h': D_h     # Hydraulic diameter of one channel [m]
                            }
        else:
            pass # To check from the code of Basile

    def _setup_fluids(self):
        """CoolProp states, inlet states, supercritical flags and saturation data."""
        for side, su in ((self.H, self.su_H), (self.C, self.su_C)):
            side.fluid = su.fluid
            side.incomp = su.AS.backend_name() == 'IncompressibleBackend' # Check if fluid is incompressible or not to avoid potential errors in the future
            if side.incomp:
                backend = "INCOMP"
            else:
                backend = "HEOS" if self.params.get('AS_Type') == 'HEOS' else "BICUBIC&HEOS"
            if side._AS_key != (backend, side.fluid): # Checks TO BE VERIFIED THESE LINES!!!
                side.AS = CP.AbstractState(backend, side.fluid)
                side.AS_sat = CP.AbstractState(backend, side.fluid)
                side.AS_heos = None if side.incomp else CP.AbstractState("HEOS", side.fluid)
                side._AS_key = (backend, side.fluid)

            side.m_dot, side.h_su, side.p_su = su.m_dot, su.h, su.p
            AS = self._state(side, side.h_su, side.p_su)
            side.T_su = AS.T()
            side.supercritical = (not side.incomp) and side.p_su >= AS.p_critical()
            if side.two_phase_possible:
                side.AS_sat.update(CP.PQ_INPUTS, side.p_su, 0)
                side.h_bub_su = side.AS_sat.hmass()
                side.AS_sat.update(CP.PQ_INPUTS, side.p_su, 1)
                side.h_dew_su = side.AS_sat.hmass()
                side.x_su = (side.h_su - side.h_bub_su) / (side.h_dew_su - side.h_bub_su)
            else:
                side.x_su = np.nan

    def _state(self, side, h, p):
        """Return an AbstractState set at (h, p). Falls back on the full EOS
        when the tabular backend fails (typically very close to saturation)."""
        try:
            side.AS.update(CP.HmassP_INPUTS, h, p)
            return side.AS
        except ValueError:
            if side.AS_heos is None:
                raise
            side.AS_heos.update(CP.HmassP_INPUTS, h, p)
            return side.AS_heos

    def _pressure_drop_correlation(self, side):
        """Pressure drop computed based on the correlations and the supply state."""
        if side.supercritical:
            phase = 'supercritical'
        elif side.two_phase_possible and 0.002 < side.x_su < 0.999:
            phase = 'two-phase'
        else:
            phase = 'liquid'

        return pressure_drop(side.AS, side.geom, side.G, side.correlation(phase, 'dp'),
                             m_dot=side.m_dot, T_wall=side.T_su)   # T_wall: rough estimate, no cells exist yet

    def _init_pressure_drops(self):
        dp_type = self.params.get('dp_type')
        for side in (self.H, self.C):
            if dp_type is None:
                side.dp = 0.0
            elif dp_type == "user_defined":
                side.dp = side.dp_user
            else:
                side.dp = self._pressure_drop_correlation(side)

    # =========================================================================
    # CELL BOUNDARIES
    # =========================================================================

    def _pressure(self, side, s):
        """Pressure at normalised enthalpy s measured from the side inlet."""
        if side.dp == 0.0:
            return np.full_like(np.asarray(s, dtype=float), side.p_su)
        return side.p_su - side.dp * np.interp(s, side.dp_s, side.dp_f)

    def _local_saturation_enthalpies(self, side, h_a, h_b):
        """Bubble and dew enthalpies at the local pressure (pressure drop
        included), strictly inside the enthalpy range ]h_a, h_b[ of the side."""
        if not side.two_phase_possible:
            return []
        out = []
        for Q_flag in (0, 1):
            h_sat = side.h_bub_su if Q_flag == 0 else side.h_dew_su
            if side.dp > 0.0:
                for _ in range(4):   # h_sat = h_sat(p(h_sat)): fixed point, converges in 2-3 steps
                    s = (h_sat - h_a) / (h_b - h_a) if side.name == 'C' else (h_b - h_sat) / (h_b - h_a)
                    side.AS_sat.update(CP.PQ_INPUTS, float(self._pressure(side, min(1.0, max(0.0, s)))), Q_flag)
                    h_sat = side.AS_sat.hmass()
            if h_a < h_sat < h_b:
                out.append(h_sat)
        return out

    def calculate_cell_boundaries(self, Q):
        """Enthalpy, pressure, temperature and quality at every cell boundary for
        the heat rate Q. Nodes are: n_disc+1 evenly spaced nodes, the phase-change
        points of both fluids, and their images on the other fluid (energy balance)."""
        H, C = self.H, self.C
        self.h_co = C.h_su + Q / C.m_dot
        self.h_ho = H.h_su - Q / H.m_dot

        n_nodes = max(int(self.params['n_disc']), 1) + 1
        q_regular = np.linspace(0.0, Q, n_nodes)
        q_sat = [C.m_dot * (h - C.h_su) for h in self._local_saturation_enthalpies(C, C.h_su, self.h_co)]
        q_sat += [H.m_dot * (h - self.h_ho) for h in self._local_saturation_enthalpies(H, self.h_ho, H.h_su)]

         # Only exact duplicates are merged: dropping regular nodes that are merely
        # close to a phase-change node changes the number of cells abruptly and
        # makes the residual jump; a tiny cell contributes a tiny w instead.
        if q_sat:
            tol = 1e-9 * Q
            keep = np.ones(n_nodes, dtype=bool)
            keep[1:-1] = np.min(np.abs(q_regular[1:-1, None] - np.array(q_sat)[None, :]), axis=1) > tol
            q_nodes = np.sort(np.concatenate([q_regular[keep], q_sat]))
        else:
            q_nodes = q_regular

        self.hvec_c = C.h_su + q_nodes / C.m_dot      # Numpy vectors with the cold enthalpies at each nodes
        self.hvec_h = self.h_ho + q_nodes / H.m_dot   # Numpy vectors with the hot enthalpies at each nodes
        self.hnorm_c = q_nodes / Q                    # Normalised position of each boundary, from 0 to 1.
        self.hnorm_h = self.hnorm_c.copy()
        self.Qvec_c = C.m_dot * np.diff(self.hvec_c)  # Heat exchanged in each cell.
        self.Qvec_h = H.m_dot * np.diff(self.hvec_h)

        # Pressures (cold inlet at s=0 / index 0, hot inlet at s=0 / index -1)
        self.pvec_c = self._pressure(C, self.hnorm_c)
        self.pvec_h = self._pressure(H, 1.0 - self.hnorm_h)

        for side, hvec, pvec in ((C, self.hvec_c, self.pvec_c), (H, self.hvec_h, self.pvec_h)):
            self._boundary_states(side, hvec, pvec)
        self.Tvec_c, self.svec_c, self.x_vec_c = C.Tvec, C.svec, C.x_flag
        self.Tvec_h, self.svec_h, self.x_vec_h = H.Tvec, H.svec, H.x_flag
        self.Tvec_sat_pure_c, self.Tvec_sat_pure_h = C.Tsat, H.Tsat

        self.DT_pinch = np.min(self.Tvec_h - self.Tvec_c)
        self.DT_ho_ci = self.Tvec_h[0] - self.Tvec_c[0]
        self.DT_hi_co = self.Tvec_h[-1] - self.Tvec_c[-1]

    def _boundary_states(self, side, hvec, pvec):
        n = len(hvec)
        side.hvec, side.pvec = hvec, pvec
        side.Tvec, side.svec = np.empty(n), np.empty(n)
        for i in range(n):
            AS = self._state(side, hvec[i], pvec[i])
            side.Tvec[i], side.svec[i] = AS.T(), AS.smass()

        side.Tsat = np.zeros(n)
        if not side.two_phase_possible:
            side.x = np.full(n, np.nan)
            side.x_flag = np.full(n, 3.0)
            return
        h_l, h_v, T_l, T_v = np.empty(n), np.empty(n), np.empty(n), np.empty(n)
        constant_p = np.all(pvec == pvec[0])
        for i in range(n):
            if constant_p and i > 0:
                h_l[i], h_v[i], T_l[i], T_v[i] = h_l[0], h_v[0], T_l[0], T_v[0]
                continue
            side.AS_sat.update(CP.PQ_INPUTS, pvec[i], 0)
            h_l[i], T_l[i] = side.AS_sat.hmass(), side.AS_sat.T()
            side.AS_sat.update(CP.PQ_INPUTS, pvec[i], 1)
            h_v[i], T_v[i] = side.AS_sat.hmass(), side.AS_sat.T()
        side.h_l, side.h_v = h_l, h_v
        side.Tsat = 0.5 * (T_l + T_v)
        side.x = (hvec - h_l) / (h_v - h_l)
        # Original convention: -2 subcooled, 2 superheated, 3 supercritical/incompressible
        side.x_flag = np.where(side.x < 0, -2.0, np.where(side.x > 1, 2.0, side.x))

    def _cell_states(self):
        """Phase, mean temperatures, pressures and qualities of every cell."""
        H,C = self.H, self.C
        Thi, Tho = self.Tvec_h[1:], self.Tvec_h[:-1]  # Tvec_h holds the hot fluid's temperature at each cell boundary
        Tci, Tco = self.Tvec_c[:-1], self.Tvec_c[1:]
        # Example:
        # Tvec_h = [40, 60, 80, 100]   # °C -> self.Tvec_h[1:] = [60, 80, 100] & self.Tvec_h[:-1] = [40, 60, 80]

        self.n_cells = n = len(Thi)
        T_wall = 0.25 * (Thi + Tho + Tci + Tco)  # /!\ WHY COMPUTED LIKE THIS????

        for side, T_in, T_out in ((H, Thi, Tho), (C, Tci, Tco)):
            side.h_mean = 0.5 * (side.hvec[1:] + side.hvec[:-1])    # Vectors containing the mean enthalpies over every cells
            side.p_mean = 0.5 * (side.pvec[1:] + side.pvec[:-1])
            side.T_mean = 0.5 * (T_in + T_out)
            side.T_wall = T_wall.copy()
            side.Tsat_mean = 0.5 * (side.Tsat[1:] + side.Tsat[:-1])
            if side.incomp:
                side.phases = np.full(n, 'liquid', dtype=object)
                side.x_mean = np.full(n, np.nan)
            elif side.supercritical:
                side.phases = np.full(n, 'supercritical', dtype=object)
                side.x_mean = np.full(n, np.nan)
            else:
                h_l_mean = 0.5 * (side.h_l[1:] + side.h_l[:-1])   # saturated liquid enthalpy, mean of the cell
                h_v_mean = 0.5 * (side.h_v[1:] + side.h_v[:-1])   # saturated vapour enthalpy, mean of the cell
                x_cell = (side.h_mean - h_l_mean) / (h_v_mean - h_l_mean)  # Only used to define the phase
                side.phases = np.where(x_cell < 0, 'liquid', np.where(x_cell > 1, 'vapor', 'two-phase')).astype(object)
                xc = np.clip(side.x, 0.0, 1.0)
                side.x_mean = 0.5 * (xc[1:] + xc[:-1])  # Physical quality passed to the two-phase correlations and the void fraction (must be between 0 and 1).

        # Hot vapour cooled by a wall below its saturation temperature: condensation
        # on the wall ('vapor-wet'); the wall temperature seen by the fluid is T_sat.
        if H.two_phase_possible:
            vap = H.phases == 'vapor'
            H.T_wall[vap] = np.maximum(T_wall[vap], H.Tsat_mean[vap])
            wet = vap & (T_wall <= H.Tsat_mean)
            H.phases[wet] = 'vapor-wet'
            H.x_mean[wet] = 1.0

        self.phases_h, self.phases_c = H.phases, C.phases

        # Counter-flow LMTD of every cell (also used by the Han boiling correlation)
        DTA = np.maximum(Thi - Tco, 1e-6)
        DTB = np.maximum(Tho - Tci, 1e-6)
        with np.errstate(divide='ignore', invalid='ignore'):
            lmtd = (DTA - DTB) / np.log(DTA / DTB)
        self.LMTD = np.where(np.abs(DTA - DTB) < 1e-9 * np.maximum(DTA, DTB), DTA, lmtd)

    def _correction_factors(self):
        """LMTD correction factor F per cell for non counter-flow arrangements."""
        n = self.n_cells
        self.F = np.ones(n)
        flow = self.params['flow_type']
        if flow == "CounterFlow" or (flow == 'Shell&Tube' and self.params.get('Tube_pass', 1) == 1):
            return
        Thi, Tho = self.Tvec_h[1:], self.Tvec_h[:-1]
        Tci, Tco = self.Tvec_c[:-1], self.Tvec_c[1:]
        for k in range(n):
            dTc, dTh = Tco[k] - Tci[k], Thi[k] - Tho[k]
            if dTc <= 1e-9 or dTh <= 1e-9:
                continue                                  # one side isothermal: F = 1
            R = dTh / dTc
            P = dTc / (Thi[k] - Tci[k])
            if flow == 'Shell&Tube':
                # (the original rule "P > 0.99 -> F = 0.01" was removed: it made the
                # residual jump by orders of magnitude where the exact F is still ~1)
                if R > 10:
                    F = 1.0
                elif self.params['Tube_pass'] % 2 == 0:
                    F = F_shell_and_tube(R, P, self.params['n_series'])
                    if not (np.isfinite(F) and F > 0):      # outside the formula's range
                        F = self._f_lmtd2(R, P)
                else:
                    F = self._f_lmtd2(R, P)
            else:
                F = self._f_lmtd2(R, P)
            self.F[k] = max(F, 0.0)

    def _f_lmtd2(self, R, P):
        """f_lmtd2 with two guards: C_r consistent with R (C_min/C_max of the
        cell), and R snapped to 1 when it is within round-off of 1 (the general
        formula of f_lmtd2 is 0/0 there and returns F = 0)."""
        if abs(R - 1.0) < 1e-6:
            R = 1.0
        F = f_lmtd2(R, P, self.params, min(R, 1.0 / R))
        return F if (np.isfinite(F) and 0.0 < F <= 1.0 + 1e-9) else 1.0

    def _compute_htc(self):
        """Heat transfer coefficient of every cell, hot side first (the cold
        side may need the hot-side value, e.g. Han boiling correlation)."""
        n = self.n_cells
        for side in (self.H, self.C):
            side.htc = np.empty(n)
            for k in range(n):
                phase = side.phases[k]

                "1) Fixed values given by the user"
                if side.htc_type == 'user_defined':
                    if phase not in side.htc_user:
                        raise ValueError(f"No fixed htc given for the '{phase}' phase on the {side.label} side "
                                         f"(cell {k}). Given: {list(side.htc_user)}.")
                    side.htc[k] = side.htc_user[phase]
                    continue

                "2) Correlation, evaluated at the mean state of the cell"
                AS = self._state(side, side.h_mean[k], side.p_mean[k])
                extra = dict(
                    T_wall = side.T_wall[k],                           # Wall temperature [K]
                    q_flux = self.Qvec_h[k] / (self.w[k] * side.geom['A']) if self.w[k] > 0 else None      # Heat flux if all cells had the same area [W/m^2] /!\ TO BE CHECKED!!
                )

                try:
                    side.htc[k] = compute_htc(AS, side.m_dot, side.geom, side.correlation(phase, 'htc'), **extra)
                except Exception as e:
                    raise RuntimeError(f"{side.label} side, cell {k} ({phase}, T = {side.T_mean[k]:.2f} K, "
                                       f"p = {side.p_mean[k]:.0f} Pa): {type(e).__name__}: {e}") from e

        self.htc_h, self.htc_c = self.H.htc, self.C.htc   # Kept on the component for plots and results
            
    def _compute_UA(self):
        """
        Available conductance of every cell [W/K], as if the cell covered the whole exchanger.

            1/UA = (1 + fouling) / (htc_h * A_h)  +  1 / (htc_c * A_c)

        The fouling factor is applied on the hot side only.
        Cell k really uses the fraction w[k] of the area, so w[k] = UA_req[k] / UA_avail[k].
        """
        H, C = self.H, self.C
        A_h, A_c = H.geom['A'], C.geom['A']    # Heat transfer area of each side [m^2]

        R_h = (1.0 + self.params.get('R_fouling', 0)) / (H.htc * A_h)   # Convection + fouling, hot side [K/W]
        R_c = 1.0 / (C.htc * A_c)                              # Convection, cold side [K/W]

        self.UA_avail = 1.0 / (R_h + R_c)

        # /!\ REMOVED THE MUTLIPASS VERSION!!!
        # /!\ IL MANQUE THE WALL CONDUCTION

    # =========================================================================
    # MAXIMUM HEAT RATE (PINCH ANALYSIS)
    # =========================================================================

    def _compute_Qmax(self):
        """External pinch, then internal pinch at the phase-change points, 
        then a final no-crossing check."""
        H, C = self.H, self.C
        p_ho, p_co = H.p_su - H.dp, C.p_su - C.dp

        "1) External pinch analysis"
        #     Computes the maximum heat transfer rate from external pinch analysis (Bell et al. 2025, eqs 4-6).

        #     Each fluid can at most leave at the inlet temperature of the other fluid:
        #         - hot side : cooled down to the cold inlet temperature
        #         - cold side :  heat up to the hot inlet temperature
        #     The smaller of the two heat transfer rate is the maximum possible heat transfer rate.

        def h_at(side, p, T):
            for dT in (0.0, 1e-3, -1e-3):
                try:
                    side.AS.update(CP.PT_INPUTS, p, T + dT)
                    return side.AS.hmass()
                except ValueError:
                    continue
            side.AS_heos.update(CP.PT_INPUTS, p, T)
            return side.AS_heos.hmass()

        Qmax_h = H.m_dot * (H.h_su - h_at(H, p_ho, max(C.T_su, T_EX_MIN_H )))
        Qmax_c = C.m_dot * (h_at(C, p_co, min(H.T_su, T_EX_MAX_C)) - C.h_su)
        Qmax = self.Qmax_ext = min(Qmax_h, Qmax_c)
        self.Qmax_int = None

        "2) Internal pinch analysis"
        #     Computes the maximum heat transfer rate from internal pinch analysis (Bell et al. 2025, eqs ???).
        self.calculate_cell_boundaries(Qmax)
        crossing = self.Tvec_c > self.Tvec_h + 1e-9
        if H.two_phase_possible and np.any(crossing):
            i = np.where(crossing & np.isclose(self.x_vec_h, 1.0, atol=1e-9))[0]
            if i.size:
                i = i[0]
                h_c_pinch = h_at(C, self.pvec_c[i], self.Tvec_h[i])
                Qmax = self.Qmax_int = C.m_dot * (h_c_pinch - C.h_su) + H.m_dot * (H.h_su - self.hvec_h[i])
                self.calculate_cell_boundaries(Qmax)
        crossing = self.Tvec_c > self.Tvec_h + 1e-9
        if C.two_phase_possible and np.any(crossing):
            i = np.where(crossing & np.isclose(self.x_vec_c, 0.0, atol=1e-9))[0]
            if i.size:
                i = i[0]
                h_h_pinch = h_at(H, self.pvec_h[i], self.Tvec_c[i])
                Qmax = self.Qmax_int = C.m_dot * (self.hvec_c[i] - C.h_su) + H.m_dot * (H.h_su - h_h_pinch)
                self.calculate_cell_boundaries(Qmax)

        # Safety net: a temperature crossing can remain (e.g. single-phase pinch
        # caused by a varying cp). The smallest temperature difference is
        # continuous in Q, so its zero is found with Brent's method.
        if self.DT_pinch < -1e-3:
            def dT_min(Q):
                self.calculate_cell_boundaries(Q)
                return self.DT_pinch
            Qmax = self.Qmax_int = brentq(dT_min, 1e-6 * Qmax, Qmax, xtol=1e-7 * Qmax)

        self.Qmax = Qmax
        return Qmax

    def _distributed_pressure_drops(self):
        """Cell-by-cell pressure drops for the current (converged) cells.
        Returns the new totals and updates the pressure profiles."""
        # /!\ ONLY USED IF DISTRIBUTED PRESSURE DROP IS SELECTED, OTHERWISE THE PRESSURE DROP IS COMPUTED ONCE AT THE BEGINNING OF THE SOLVE.
        frac = self.w/np.sum(self.w) # Length fraction of each cell, used to compute the pressure drop in each cell
        new = {}
        for side in (self.H, self.C):
            dp_cells = np.empty(self.n_cells)
            for k in range(self.n_cells):
                AS = self._state(side, side.h_mean[k], side.p_mean[k])
                extra = dict(
                    m_dot     = side.m_dot,
                    T_wall    = side.T_wall[k],
                    x         = side.x_mean[k],
                    Q         = self.Qvec_h[k],
                    q         = self.Qvec_h[k] * self.n_cells / side.geom['A'],
                    DT_lm     = self.F[k] * self.LMTD[k],
                    htc_other = self.H.htc[k] if side is self.C else np.nan,
                    h_min     = side.hvec[k],
                    h_max     = side.hvec[k + 1],
                )
                dp_full = pressure_drop(AS, side.geom, side.G, side.correlation(side.phases[k], 'dp'),
                                        m_dot=side.m_dot, T_wall=side.T_wall[k])   # Over the full flow length [Pa]
                dp_cells[k] = frac[k] * dp_full                                     # Over the cell length [Pa]
            side.dp_cells = dp_cells
            "Cumulative profile, measured from the supply of the side"
            # The cold fluid enters at index 0 and the hot fluid at the last index:
            # the hot-side arrays are reversed so that both start at their own supply.
            if side is self.C:
                s_su = self.hnorm_c         # Position of each boundary from the supply [-]
                dp_su = dp_cells            # Cell pressure drops, ordered from the supply [Pa]
            else:
                s_su = 1.0 - self.hnorm_h[::-1]
                dp_su = dp_cells[::-1]

            dp_total = float(np.sum(dp_cells))                      # Total pressure drop of the side [Pa]
            dp_cum   = np.concatenate([[0.0], np.cumsum(dp_su)])   # Pressure lost from the supply to each boundary [Pa]
            f_su     = dp_cum / dp_total if dp_total > 0.0 else s_su   # Fraction of dp_total lost at each boundary [-]

            new[side.name] = (dp_total, s_su, f_su)   # Not applied here: solve() checks convergence first

        self.dp_cells_h, self.dp_cells_c = self.H.dp_cells, self.C.dp_cells   # Kept for results and plots
        return new

    def _limit_dp(self, side, dp_total):
        """Limit the total pressure drop to 80% of the supply pressure."""
        dp_max = 0.8 * side.p_su
        if dp_total > dp_max:
            warnings.warn(f"{side.label} side pressure drop ({dp_total:.0f} Pa) above 80% of the inlet "
                                      f"pressure; limited to that value.", stacklevel=3)
            dp_total = dp_max
        return dp_total


    # =========================================================================
    # SOLVE
    # =========================================================================

    def _brent(self, a, b, rtol):
        Q = brentq(self._f, a, b, xtol=1e-9 * b, rtol=rtol)
        if self._last_Q != Q:
            self.objective_function(Q)
        return Q


    def solve(self, only_external=False, and_solve=True):
        
        "1) Checks and preparation"
        if not self.check_calculable():
            raise ValueError("Component not calculable: check the inputs.")

        # Some additional geometires need to be computed
        self._setup_geometry() # Pour l'instant rien, à remplir

        self.check_parametrized()

        self._setup_fluids()  # To set up the cold and hot side parameters and all

        if not (self.su_H.T - self.su_C.T > 1e-2 and self.su_H.m_dot > 0 and self.su_C.m_dot > 0):
            raise ValueError("Hot and cold temperatures seem to be reversed or a flow rate is not positive.")

        "2) Solver settings"

        rtol_dp = self.params.get('DP_rtol', 1e-3) # Relative tolerance for the pressure drops
        max_outer = int(self.params.get('max_outer_iter', 30)) # Maximum number of pressurre drop iterations
        dp_disc = self.params.get('dp_type') == "correlation_disc" # True when pressure drops are computed cell by cell.
        self._init_pressure_drops() # Gives a first total pressure drop
        Q = None
        self.outer_iterations = 0

        "3) Solve the heat rate, with distributed pressure drops if requested"

        for it in range(max_outer):
            self.outer_iterations = it +1
            self._compute_Qmax()
            if only_external or not and_solve: # Case where just want a first estimation of the maximum power
                    self.Q_dot = self.Qmax
                    return self.Qmax
            Q = self._solve_heat_rate(Q)
            if not dp_disc:
                break
            # Once the heat rate has been solved and the cells are known
            # recompute the pressure drop son both sides, and check if they have changed significantly. If so, repeat the process.
            new_dp_profile = self._distributed_pressure_drops()   # Dict that maps each side to its computed pressure-drop profile (dp_total, s_su, f_su)
            "Check convergence of the distributed pressure drops"
            dp_new = {}          # New values if each side, applied only if another pass is needed
            converged = True    # Becomes False as soon as one side still changes

            for side in (self.H, self.C):
                dp_total, s_su, f_su = new_dp_profile[side.name]
                dp_total = self._limit_dp(side, dp_total)      # At most 80% of the supply pressure

                tol = max(1.0, rtol_dp * dp_total)             # Absolute tolerance for the pressure drop convergence [Pa]
                if abs(dp_total - side.dp) > tol:
                    converged = False
                dp_new[side.name] = (dp_total, s_su, f_su)

            if converged:
                break                   # Q was solved with pressure drops that no longer change

            "Not converged: use the new pressure drops for the next pass"
            for side in (self.H, self.C):
                side.dp, side.dp_s, side.dp_f = dp_new[side.name]

        else:
            # Runs only if the for loop ended without 'break'
            warnings.warn(f"Distributed pressure drops not converged after {max_outer} passes.", stacklevel=2)

        "4) Finalize the solution"
        self._finalize(Q)
        return Q

    def _finalize(self, Q):
        H, C = self.H, self.C
        self.Q_dot = Q
        self.Q.set_Q_dot(Q)
        self.epsilon_th = Q / self.Qmax
        self.residual = 1.0 - float(np.sum(self.w))
        self.dp_h = float(self.pvec_h[-1] - self.pvec_h[0])
        self.dp_c = float(self.pvec_c[0] - self.pvec_c[-1])
        self.p_ho, self.p_co = self.pvec_h[0], self.pvec_c[-1]
        self.Avec_h = self.w * H.geom['A']   # Heat transfer area of each cell [m^2]
        self.Avec_c = self.w * C.geom['A']
        self._compute_charge()  # Computes the mass of fluid in each cell and the void fraction

        self.ex_H.reset()
        self.ex_H.set_fluid(H.fluid)
        self.ex_H.set_p(self.pvec_h[0])
        self.ex_H.set_h(self.hvec_h[0])
        self.ex_H.set_m_dot(H.m_dot)
        self.ex_C.reset()
        self.ex_C.set_fluid(C.fluid)
        self.ex_C.set_p(self.pvec_c[-1])
        self.ex_C.set_h(self.hvec_c[-1])
        self.ex_C.set_m_dot(C.m_dot)
        self.H_su, self.C_su, self.H_ex, self.C_ex = self.su_H, self.su_C, self.ex_H, self.ex_C
        self.solved = True

    def _compute_charge(self):
        """Fluid mass in every cell: cell volume times the mean of the boundary
        densities (void-fraction model in two-phase)."""
        void_fraction_model = self.params.get('void_fraction_model', 'Homogeneous')
        frac = self.w / np.sum(self.w)  # Length fraction of each cell, used to compute the pressure drop in each cell
        for side in (self.H, self.C):
            n = len(side.hvec)
            rho, void_fraction = np.empty(n), np.full(n, -1.0)
            for i in range(n):
                AS = self._state(side, side.hvec[i], side.pvec[i])
                x = side.x[i] if side.two_phase_possible else np.nan
                if side.two_phase_possible and 0.0 < x < 1.0:
                    side.AS_sat.update(CP.PQ_INPUTS, side.pvec[i], 0)
                    rho_l = side.AS_sat.rhomass()
                    side.AS_sat.update(CP.PQ_INPUTS, side.pvec[i], 1)
                    rho_v = side.AS_sat.rhomass()
                    D = side.geom.get('D')
                    A_cross = np.pi * D**2 / 4.0 if D is not None else 0.0
                    G = getattr(side, 'G', None)
                    m_dot = G * A_cross if G is not None else 0.0
                    void_fraction[i] = compute_void_fraction(AS, side.geom, m_dot, void_fraction_model=void_fraction_model) # /!\ Change from m_dot to G !!!
                    rho[i] = compute_two_phase_density(x=x, rho_l=rho_l, rho_g=rho_v, alpha=void_fraction[i])
                else:
                    rho[i] = AS.rhomass()
                    
            side.rhovec, side.void_fraction = rho, void_fraction
            side.Vvec = side.geom['V'] * frac
            side.charge_vec = side.Vvec * 0.5 * (rho[1:] + rho[:-1])  # Mass of fluid in each cell [kg]

        H, C = self.H, self.C
        self.charge_H = float(np.sum(H.charge_vec))
        self.charge_C = float(np.sum(C.charge_vec))
        self.void_fraction_H = float(np.mean(H.void_fraction[1:]))
        self.void_fraction_C = float(np.mean(C.void_fraction[1:]))



    def _solve_heat_rate(self, Q_guess=None):
        """Find Q such that sum(w) = 1, within ]0, Qmax[."""
        Qmax = self.Qmax
        rtol = self.params.get('Q_rtol', 1e-6)
        # Residual values are cached within one solve: brentq re-evaluates the
        # bracket ends, which were just computed to check the bracket.
        self._f_cache = {}

        def f(Q):
            if Q not in self._f_cache:
                self._f_cache[Q] = self.objective_function(Q)
            return self._f_cache[Q]
        
        self._f = f
        Q_hi = Qmax * (1.0 - 1e-5)

        # Warm start: narrow bracket around the previous solution
        if Q_guess is not None and Q_guess < Q_hi:
            a, b = 0.98 * Q_guess, min(1.02 * Q_guess, Q_hi)
            fa, fb = f(a), f(b)
            if fa > 0 > fb:
                return self._brent(a, b, rtol)

        f_hi = f(Q_hi)
        if f_hi >= 0.0:
            # The exchanger is large enough to reach the pinch-limited heat rate
            self.pinched = True
            if self._last_Q != Q_hi:
                self.objective_function(Q_hi)
            return Q_hi
        self.pinched = False
        for frac in (1e-2, 1e-4, 1e-6):
            Q_lo = frac * Qmax
            if f(Q_lo) > 0.0:
                return self._brent(Q_lo, Q_hi, rtol)
        raise RuntimeError("No heat rate found with sum(w) < 1, even at 1e-6*Qmax: "
                            "check the geometry and the heat transfer coefficients.")


    def objective_function(self, Q):
        """Return 1 - sum(w) for the heat rate Q (w_k: area fraction of cell k)."""
        self.calculate_cell_boundaries(Q)  # Determines the different cell boundaries
        self._cell_states()               # Computes the mean properties of each cells and the LMTDs
        self._correction_factors()        # Computes the correction factor F for the other cases than counter-flow heat exchanegers

        # Conductance each cell needs: does not depend on the htc, so computed here
        self.UA_req = self.Qvec_h / (self.F * self.LMTD)    # computes the conductance UA that eatch cell needs to transfer its heat (Qvec_h is an array)

        # If one of the htc correlation depends on the heat flux, the w for this Q must be computed based on an iteration
        needs_q_flux = any(
            side.htc_type == 'correlation' and side.correlation(phase, 'htc') in NEEDS_Q_FLUX
            for side in (self.H, self.C)
            for phase in side.phases
        )

        # Iteration loop until w converges
        self.w = np.full(self.n_cells, 1.0 / self.n_cells)
        for _ in range(20):
            self._compute_htc()
            self._compute_UA()                
            w_new = self.UA_req / self.UA_avail              # Ratio of the exchanger's area the cell needs
            converged = np.allclose(w_new, self.w, rtol=1e-4)
            self.w = w_new
            if converged or not needs_q_flux:
                break
        else:
            warnings.warn(f"Area fractions did not converge for Q = {Q:.1f} W")


        self._compute_htc()               # Computes the heat transfer coefficients for each cell for both streams
        self._compute_UA()                # Computes the conductance of every cell as if it covered the whole heat exchanger

        # if self._multipass -> passer pour le moment! c'est plus pour Basile tout ça
        self.w_sum = float(np.sum(self.w))
        self.eval += 1
        self._last_Q = Q
        return 1.0 - self.w_sum

    def plot_Ts_pair(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for side, svec, Tvec in (('H', self.svec_h, self.Tvec_h), ('C', self.svec_c, self.Tvec_c)):
                fluid = self.H.fluid if side == 'H' else self.C.fluid
                diagram = PropertyPlot(fluid, "TS", unit_system="EUR")
                plt.plot(0.001 * svec, Tvec - 273.15, 's-')
                diagram.calc_isolines()
                diagram.title(f"{'Hot' if side == 'H' else 'Cold'} side T-s diagram. Fluid: {fluid}")
                diagram.show()

    def plot_cells(self, fName='', dpi=400):
        plt.figure(figsize=(4, 3))
        plt.plot(self.hnorm_h, self.Tvec_h, 'rs-')
        plt.plot(self.hnorm_c, self.Tvec_c, 'bs-')
        plt.xlim(0, 1)
        plt.ylabel('T [K]')
        plt.xlabel(r'$\hat h$ [-]')
        plt.grid(True)
        plt.tight_layout(pad=0.2)
        if fName != '':
            plt.savefig(fName, dpi=dpi)
        plt.show()

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        for name in ('su_C', 'su_H', 'ex_C', 'ex_H'):
            c = getattr(self, name)
            print(f"  - {name}: fluid={c.fluid}, T={c.T}, p={c.p}, m_dot={c.m_dot}")
        print("======================")







    

            





