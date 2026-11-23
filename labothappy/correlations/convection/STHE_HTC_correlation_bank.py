# -*- coding: utf-8 -*-
"""
Shell-side heat transfer coefficient of shell & tube heat exchangers, in the
(AS, geom, cond) formalism of the BPHX correlations.

cond : dictionary of the cell conditions built by the model (_cond_dict). The shell correlations read
       m_dot (flow of the whole side), mu, mu_w (optional: without it mu_w = mu), Pr, k, cp.
geom : geometry of the side: all the heat exchanger parameters (Tube_OD, pitch_ratio, Shell_ID,
       central_spacing, tube_layout, ...). n_parallel (optional, default 1) divides the flow rate.
"""

import warnings

import numpy as np

from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props
from labothappy.correlations.convection.BPHE_HTC_correlation_bank import BPHE_cooper_boiling_HTC

#%% GEOMETRY HELPERS

def s_max_kern(Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout): # Maximum flow section (m**2)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    Shell_ID        : Shell inside diameter [m]
    central_spacing : Spacing between two baffles [m]
    tube_layout     : Angle of tube layout [°]
    
    Output
    ------
    
    s_max : Maximum flow area in the tube bank [m^2]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio # Tube Pitch
    s_max = (Shell_ID / p_T)*(p_T -  Tube_OD) * central_spacing

    return s_max


def s_max(Tube_OD, pitch_ratio, Shell_ID, D_OTL, central_spacing, tube_layout): # Maximum flow section (m**2)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    Shell_ID        : Shell inside diameter [m]
    D_OTL           : Outer tube limit diameter [m]
    central_spacing : Spacing between two baffles [m]
    tube_layout     : Angle of tube layout [°]: 0 square inline, 45 square staggered, 60 triangular staggered
    
    Output
    ------
    
    s_max : Maximum flow area in the tube bank [m^2]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio
    
    # TODO: the three effective pitch factors below are kept from the previous version and are not
    #       verified; the labels of 45 and 60 were swapped in the previous version of this function.
    if tube_layout == 0: # Square inline tube arrangement
    
        p_T_eff = p_T
    
    elif tube_layout == 45: # Square staggered tube arrangement
    
        p_T_eff = (2**(1/2))*p_T
    
    elif tube_layout == 60: # Triangular staggered tube arrangement
    
        p_T_eff = (2/3**(1/2))*p_T
    
    else:
        raise ValueError("Tube arrangement shall be either square inline (tube_layout = 0), "
                         "square staggered (tube_layout = 45) or triangular staggered (tube_layout = 60)")

    # (the bundle-to-shell gap term is not divided by p_T_eff: both terms are lengths)
    s_max = central_spacing * ((Shell_ID - D_OTL) + (D_OTL - Tube_OD)*(p_T - Tube_OD)/p_T_eff)
    
    return s_max


def s_L(Baffle_cut, Shell_ID, n_tubes, Tube_OD):
    """
    Input
    -----
    
    Baffle_cut      : Baffle cut [% of the shell diameter]
    Shell_ID        : Shell inside diameter [m]
    n_tubes         : Number of tubes [-]
    Tube_OD         : Tube outside diameter [m]
    
    Output
    ------
    
    s_L : Flow area in the longitudinal shell direction (baffle window) [m^2]:
          the fraction Baffle_cut of the shell section, minus the tubes in it (uniform tube density)
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    return (Baffle_cut/100)*(np.pi/4)*(Shell_ID**2 - n_tubes*Tube_OD**2)


def d_h(Tube_OD, pitch_ratio, tube_layout): # Hydraulic diameter (m)
    """
    Input
    -----
    
    Tube_OD         : Tube outside diameter [m]
    pitch_ratio     : Tube pitch_ratio (tube spacing/tube outside diameter) [-]
    tube_layout     : Angle of tube layout [°]
    
    Output
    ------
    
    D_h : Flow Hydraulic diameter in tube banks [m]
    
    Reference
    ---------
    https://cheguide.com/heat_exchanger_rating.html
    
    """
    
    p_T = Tube_OD * pitch_ratio # Tube Pitch

    if tube_layout == 0 or tube_layout == 45: # Squared tube arrangement
    
        D_h = 4*( (p_T**2 - np.pi*(Tube_OD**2 / 4)) / (np.pi*Tube_OD) )

    else:
    
        D_h = 2 * (3**(1/2) / np.pi) * (p_T**2 / Tube_OD) - Tube_OD
        
    return D_h

#%% KERN

def STHE_kern_shell_HTC(AS, geom, cond):
    """
    Kern shell-side heat transfer coefficient [W/(m^2 K)].

    cond : m_dot (flow of the whole side) [kg/s], mu, mu_w (optional) [Pa s], Pr, k [W/(m K)]
    geom : Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout; n_parallel (optional)

    Reference
    ---------
    Process Heat Transfer - D. Q. Kern
    """
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)        # flow rate per shell
    mu, Pr, k = cond['mu'], cond['Pr'], cond['k']
    mu_w = cond.get('mu_w', mu)

    S_T = s_max_kern(geom['Tube_OD'], geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])          # m^2
    D_hydro = d_h(geom['Tube_OD'], geom['pitch_ratio'], geom['tube_layout'])

    # Reynolds number on the equivalent diameter and the cross-flow mass velocity
    Re = (m_dot / S_T) * D_hydro / mu

    # McAdams (with baffles)
    Nu = 0.36 * Pr**(1/3) * Re**0.55 * (mu / mu_w)**0.14

    return Nu * k / D_hydro

#%% DONOHUE

def STHE_donohue_shell_HTC(AS, geom, cond):
    """
    Donohue shell-side heat transfer coefficient [W/(m^2 K)]. No viscosity correction.

    cond : m_dot (flow of the whole side) [kg/s], mu [Pa s], Pr, k [W/(m K)]
    geom : Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout, Baffle_cut [%], n_tubes;
           n_parallel (optional)

    Reference
    ---------
    Process Heat Transfer - D. Q. Kern
    """
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)        # flow rate per shell
    mu, Pr, k = cond['mu'], cond['Pr'], cond['k']
    d_o = geom['Tube_OD']

    # Cross-flow section (Kern) and longitudinal section (baffle window)
    S_T = s_max_kern(d_o, geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])          # m^2
    S_L = s_L(geom['Baffle_cut'], geom['Shell_ID'], geom['n_tubes'], d_o)   # m^2

    # Re = rho V_R d_o / mu with V_R = sqrt(V_L V_T) and V = m_dot / (rho S): the density cancels
    Re = m_dot * d_o / (mu * np.sqrt(S_T * S_L))

    Nu = 0.22 * Re**0.6 * Pr**0.33

    return Nu * k / d_o

#%% BELL-DELAWARE

def bell_delaware_coefs(layout_angle, Re):
    """

    Parameters
    ----------
    layout_angle : Tube layout angle
    Re : Reynolds number

    Returns
    -------
    a : coefficients
    b : coefficients
        
    Reference
    -------   
    https://cheguide.com/heat_exchanger_rating.html

    """
    # Define the breakpoints and corresponding 'a' and 'b' values
    breakpoints = [1, 10, 1e2, 1e3, 1e4, 1e5]
    
    if layout_angle == 0: # Squared Inline tube arrangement
        a_values = [
            [0.97, -0.667, 1.187, 0.37],
            [0.97, -0.667, 1.187, 0.37],
            [0.9, -0.631, 1.187, 0.37],
            [0.408, -0.46, 1.187, 0.37],
            [0.107, -0.266, 1.187, 0.37],
            [0.37, -0.395, 1.187, 0.37],
        ]
        
        b_values = [
            [35, -1, 6.3, 0.378],
            [35, -1, 6.3, 0.378],
            [32.1, -0.963, 6.3, 0.378],
            [6.09, -0.602, 6.3, 0.378],
            [0.0815, 0.022, 6.3, 0.378],
            [0.391, -0.148, 6.3, 0.378],
        ]
        
    elif layout_angle == 45: # Squared Staggered tube arrangement: 
        a_values = [
            [1.55, -0.667, 1.93, 0.5],
            [1.55, -0.667, 1.93, 0.5],
            [0.498, -0.656, 1.93, 0.5],
            [0.73, -0.5, 1.93, 0.5],
            [0.37, -0.396, 1.93, 0.5],
            [0.37, -0.396, 1.93, 0.5],
        ]   
        b_values = [
            [32, -1, 6.59, 0.52],
            [32, -1, 6.59, 0.52],
            [26.2, -0.913, 6.59, 0.52],
            [3.5, -0.476, 6.59, 0.52],
            [0.333, -0.136, 6.59, 0.52],
            [0.303, -0.126, 6.59, 0.52],
        ]
    
    elif layout_angle == 60: # Triangular Staggered tube arrangement
        a_values = [
            [1.4, -0.667, 1.45, 0.519],
            [1.4, -0.667, 1.45, 0.519],
            [1.36, -0.657, 1.45, 0.519],
            [0.593, -0.477, 1.45, 0.519],
            [0.321, -0.388, 1.45, 0.519],
            [0.321, -0.388, 1.45, 0.519],
        ]
        b_values = [
            [48, -1, 7, 0.5],
            [48, -1, 7, 0.5],
            [45.1, -0.973, 7, 0.5],
            [4.57, -0.476, 7, 0.5],
            [0.486, -0.152, 7, 0.5],
            [0.372, -0.123, 7, 0.5],
        ]
    
    else:
        raise ValueError("Tube arrangement shall be either square inline (tube_layout = 0), "
                         "square staggered (tube_layout = 45) or triangular staggered (tube_layout = 60)")

    # Interpolate on the log scale
    log_Re = np.log10(Re)
    log_breakpoints = np.log10(breakpoints)
    
    a = np.zeros_like(a_values[0])
    b = np.zeros_like(b_values[0])
    
    for i in range(len(a)):
        a[i] = np.interp(log_Re, log_breakpoints, [a_values[j][i] for j in range(len(a_values))])
        b[i] = np.interp(log_Re, log_breakpoints, [b_values[j][i] for j in range(len(b_values))])
    
    return a, b


def STHE_bell_delaware_shell_HTC(AS, geom, cond):
    """
    Bell-Delaware shell-side heat transfer coefficient [W/(m^2 K)].

    cond : m_dot (flow of the whole side) [kg/s], mu, mu_w (optional) [Pa s], cp [J/(kg K)], Pr
    geom : Tube_OD, pitch_ratio, Shell_ID, D_OTL, central_spacing, tube_layout, Baffle_cut [%],
           clear_TB, clear_BS, n_tubes, N_strips, Tube_L, Tubesheet_t, inlet_spacing, outlet_spacing;
           n_parallel (optional)

    References
    -------
    https://cheguide.com/heat_exchanger_rating.html

    + 

    https://fr.scribd.com/document/351507985/Shell-and-Tube-HX

    tube_layout : 
        30° or 60° (depends on definition) - Clean fluid (laminar + turbulent) (best for compactness)
        45° - Fouled + Laminar
        0 or 90° - Fouled + Turbulent (best for maintenance)
    """
    m_dot_shell = cond['m_dot'] / geom.get('n_parallel', 1)  # flow rate per shell
    mu, cp, Pr = cond['mu'], cond['cp'], cond['Pr']
    mu_w = cond.get('mu_w', mu)
    p = geom                                                  # geometry (all the heat exchanger parameters)

    "1) HTC"
    
    "1.1) Ideal HTC"

    # Effective sections and hydraulic diameter
    S_m = s_max(p['Tube_OD'], p['pitch_ratio'], p['Shell_ID'], p['D_OTL'], p['central_spacing'], p['tube_layout']) # m^2
    
    # Transversal to tube Flow speed and Reynolds number
    G_s = m_dot_shell/S_m
    
    Re = G_s*(p['Tube_OD']/mu)
        
    # Coefficients required for modelling tube bank flow
    a_vec, b_vec = bell_delaware_coefs(p['tube_layout'], Re)
    a = a_vec[2]/(1+0.14*Re**a_vec[3])
    j_i = a_vec[0]*(1.33/p['pitch_ratio'])**a * Re**a_vec[1]
        
    # Ideal Tube bank heat transfer coefficient
    h_id = j_i * cp * (m_dot_shell/S_m) * (Pr)**(-2/3) * (mu/mu_w)**0.14 
        
    "1.2) Correction factor for Baffle Window Flow : J_c"
    
    # Diameter effectively covered by the tube
    D_CTL = p['D_OTL'] - p['Tube_OD']
    
    theta_CTL = 2*np.arccos((p['Shell_ID'] * (1 - 2*p['Baffle_cut']/100))/D_CTL)
    F_w = (theta_CTL - np.sin(theta_CTL))/(2*np.pi)
    F_c = 1 - 2*F_w
    
    J_c = 0.55 + 0.72*F_c
    
    "1.3) Correction factor for Baffle Leakage : J_l"
    theta_ds = 2*np.arccos(1-2*(p['Baffle_cut']/100))
    
    # Clearances (tube to baffle and shell to baffle)
    S_tb = (np.pi/4)*((p['Tube_OD'] + p['clear_TB'])**2 - p['Tube_OD']**2)*p['n_tubes']*(1-F_w) # clear_TB - max - 0.794 mm in TEMA stdrads and between 0.397 mm 
    S_sb = 0.00436 * p['Shell_ID'] * p['clear_BS'] * (180/np.pi) * (2*np.pi - theta_ds) # - L_sb = 2 mm if D_s < 0.4m and = 1.6 + 0.004*D_s
    
    r_L = (S_sb + S_tb)/S_m
    r_S = S_sb /(S_sb + S_tb)

    J_l = 0.44*(1-r_S) + (1 - 0.44*(1-r_S))*np.e**(-2.2*r_L)
        
    "1.4) Correction factor for Bundle Bypass effects : J_b"
    
    # Tube pitch
    P_T = p['Tube_OD']*p['pitch_ratio']

    # Tube pitch relative to the flow    
    if p['tube_layout'] == 0: # Square inline tube arrangement
        P_p = P_T # for 90° layout
    elif p['tube_layout'] == 60: # Triangular staggered tube arrangement
        P_p = P_T * (3**(0.5)/2) 
    elif p['tube_layout'] == 45: # Square staggered tube arrangement
        P_p = P_T / (2**(0.5)) 
    else:
        raise ValueError(f"Unknown tube_layout {p['tube_layout']!r}: 0, 45 or 60.")
    
    # Number of effective tube rows crossed between two baffles
    N_tcc = (p['Shell_ID']/P_p)*(1 - 2*p['Baffle_cut']/100)
    
    # Number of sealing strips
    N_ss = p['N_strips']
    
    # Sealing strips ratio
    r_ss = N_ss/N_tcc

    if r_ss > 0.5:
        J_b = 1
    else:
        # Bypass flow area
        S_b = p['central_spacing'] * (p['Shell_ID'] - p['D_OTL'] - p['Tube_OD']/2) # TODO: check the - Tube_OD/2
        
        if Re < 100:
            C_j = 1.35
        else:
            C_j = 1.25
            
        J_b = np.e**(-C_j*(S_b / S_m)*(1 - (2*r_ss)**(1/3)))
    
    "1.5) Correction factor for adverse temperature gradient : J_r"
    
    # Effective number of baffles
    N_B = 1 + int((p['Tube_L'] - p['Tubesheet_t'] - p['inlet_spacing'] - p['outlet_spacing']) / p['central_spacing'])
    
    # Number of effective tube rows crossed in the main flow region ?
    N_tcw = (0.8/P_p)*(p['Shell_ID']*(p['Baffle_cut']/100) - ((p['Shell_ID']-(p['D_OTL']-p['Tube_OD']))/2))
            
    if Re >= 100:
        J_r = 1
    else:
        N_C = (N_tcw + N_tcc)*(1 + N_B)
        J_rl = (10/N_C)**0.18   

        if Re <= 20:
            J_r = J_rl
        else:
            J_r = J_rl + (20-Re)*(J_rl - 1)/80


    "1.6) Correction factor for unequal baffle spacing : J_s"
    
    if Re >= 100:
        n1 = 0.6
    else:
        n1 = 1/3
        
    J_s = ((N_B-1)+(p['inlet_spacing']/ p['central_spacing'])**(1-n1) + (p['outlet_spacing']/p['central_spacing'])**(1-n1))/((N_B-1)+(p['inlet_spacing']/p['central_spacing']) + (p['outlet_spacing']/p['central_spacing']))
    
    "1.7) Final Heat Transfer Coefficient"
    
    h = h_id*J_c*J_l*J_b*J_r*J_s

    return h

#%% SHELL-SIDE BOILING (pool boiling on the tubes, no bundle factor)

def STHE_mostinski_boiling_HTC(AS, geom, cond):
    """
    Nucleate pool boiling heat transfer coefficient on the shell side of a shell & tube
    evaporator using the Mostinski reduced-pressure corresponding states method.

    Parameters
    ----------
    AS : CoolProp abstract state object
        Fluid property provider.
    geom : dict
        Geometry dictionary (not used by this method, required for interface consistency).
    cond : dict
        Condition dictionary containing:
        - p : saturation pressure [Pa]
        - q : local heat flux on the outer tube surface [W/m^2]

    Returns
    -------
    h : float
        Boiling heat transfer coefficient [W/(m^2 K)].

    Reference
    ---------
    Mostinski, I.L. (1963) Application of the rule of corresponding states for the calculation
    of heat transfer and critical heat flux. Teploenergetika 4, 66 (English abstract:
    Br. Chem. Eng. 8(8), 586, 1963). As presented in Cao, E. (2009) Heat Transfer in Process
    Engineering, McGraw-Hill, and in Serth, R.W. (2014) Process Heat Transfer, 2nd ed.,
    Academic Press.

    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.boiling_nucleic.Montinsky (heat-flux form), which cites the references above. The
    constant 0.00417 requires the critical pressure in kPa (some sources write 3.7e-5 or
    3.596e-5 with p_crit in Pa).

    Validity
    --------
    Single tube / pool boiling; NO tube-bundle factor is applied (the bundle effect of kettle
    reboilers is not included: labothappy limitation).
    Independent of geometry, mass flux, and quality (geom, cond['x'] not used).
    """
    # 1) Read saturation pressure and local heat flux from conditions dictionary
    p_sat = cond['p']
    q = cond['q']

    # 2) Obtain critical pressure from the abstract state and compute reduced pressure
    p_crit = AS.p_critical()
    p_red = p_sat / p_crit

    # 3) Check that the reduced pressure is within the valid range ]0, 1[
    if not 0.0 < p_red < 1.0:
        raise ValueError(f"Mostinski: reduced pressure {p_red:.3f} outside ]0, 1[.")

    # 4) Check that the heat flux is finite and strictly positive
    if not (np.isfinite(q) and q > 0.0):
        raise ValueError(f"Mostinski: heat flux q = {q} must be finite and > 0.")

    # 5) Convert critical pressure from Pa to kPa
    p_crit_kPa = p_crit / 1000.0

    # 6) Compute the pressure factor F_p as a function of reduced pressure
    term_one = 1.8 * p_red**0.17
    term_two = 4.0 * p_red**1.2
    term_three = 10.0 * p_red**10
    F_p = term_one + term_two + term_three

    # 7) Compute the nucleate boiling heat transfer coefficient
    h = 0.00417 * p_crit_kPa**0.69 * q**0.7 * F_p

    # 8) Return the heat transfer coefficient
    return h
#%% TWO-PHASE : SHELL-SIDE CONDENSATION AND BUNDLE BOILING

def STHE_rose_kern_condensing_HTC(AS, geom, cond):
    """
    Mean condensation heat transfer coefficient on the outer surface of the tubes
    of a horizontal tube bundle (shell-side condensation, vapour in downward/cross
    flow) at cell conditions, using the single-tube Rose / Fujii / Nusselt
    correlations times the Kern inundation factor.

    Parameters
    ----------
    AS : CoolProp AbstractState
        CoolProp abstract state object for property evaluation.
    geom : dict
        Geometry dictionary containing:
        - Tube_OD : Tube outside diameter [m]
        - pitch_ratio : Tube pitch-to-diameter ratio [-]
        - Shell_ID : Shell inside diameter [m]
        - central_spacing : Spacing between baffles [m]
        - tube_layout : Angle of tube layout [deg] (0 or 45 for square, others for triangular)
        - N_rows : Number of tube rows in a column (optional) [-]
        - n_parallel : Number of parallel paths (optional, default 1) [-]
    cond : dict
        Condition dictionary containing:
        - m_dot : Flow rate of the whole shell side [kg/s]
        - p : Pressure [Pa]
        - x : Vapour quality [-]
        - T_w : Wall temperature [K]

    Returns
    -------
    h : float
        Mean condensation heat transfer coefficient [W/(m^2 K)] based on the
        outer tube surface area (pi * D * L), bundle-averaged (inundation included).

    References
    ----------
    - Bonneau, Josset, Melot, Auvity (2019), Synthesis of heat transfer correlations
      for condensation outside horizontal tube bundles, Nuclear Engineering and Design 349,
      doi 10.1016/j.nucengdes.2019.04.005, section 4 (recommendations); equation numbers
      below are those of Bonneau 2019.
    - Rose, J.W. (1984), Int. J. Heat Mass Transfer 27(1), 39-47.
    - Fujii et al. (1979), as given by Bonneau 2019 Eq. 70.
    - Nusselt, W. (1916), VDI-Z. 60, 541-546, 569-575.
    - Kern, D.Q. (1958), inundation factor N^(-1/6), as given by Bonneau 2019 Eq. 110.

    Validity
    --------
    - Fluids: steam and refrigerants (Bonneau 2019 database).
    - Rose single-tube correlations: F < 1e3 (gravity-controlled domain for F > 1e3).
    - See Bonneau 2019 section 4 for detailed validity limits.

    Definitions
    -----------
    - Nu and Re evaluated on the tube outer diameter (Tube_OD = D).
    - Liquid properties evaluated at saturation temperature via _sat_props.
    - dT = T_sat - T_w.
    - U = G * x / rho_v, where G is the mass flux on the Kern cross-flow area
      at the shell centre line (labothappy assumption).
    - Inundation (labothappy assumption): without geom['N_rows'], the tubes of each vertical
      column are counted from Shell_ID, the pitch and the layout (square 0: vertical pitch p,
      column spacing p; rotated square 45: sqrt(2) p, p / sqrt(2); triangular: sqrt(3) p, p / 2)
      and the Kern factor N^(-1/6) is averaged over the columns, weighted by their tube count.
    """
    # 1) Retrieve geometry and flow conditions
    m_dot = cond['m_dot']
    n_parallel = geom.get('n_parallel', 1)
    m_dot_shell = m_dot / n_parallel
    
    Tube_OD = geom['Tube_OD']
    pitch_ratio = geom['pitch_ratio']
    Shell_ID = geom['Shell_ID']
    central_spacing = geom['central_spacing']
    tube_layout = geom['tube_layout']
    
    S = s_max_kern(Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout)
    G = m_dot_shell / S
    
    x = cond['x']
    if x < 0.0 or x > 1.0:
        warnings.warn(f"Vapour quality x = {x} outside [0, 1], clipping to bounds.")
        x = np.clip(x, 0.0, 1.0)
        
    g = 9.81

    # 2) Obtain saturation properties at pressure p
    s = _sat_props(AS, cond['p'], cond)
    rho_l = s['rho_l']
    mu_l = s['mu_l']
    k_l = s['k_l']
    cp_l = s['cp_l']
    Pr_l = s['Pr_l']
    rho_v = s['rho_v']
    mu_v = s['mu_v']
    h_lv = s['h_lv']
    T_sat = s['T_sat']

    # 3) Temperature difference and check
    dT = T_sat - cond['T_w']
    if dT <= 0.0:
        raise ValueError(f"No condensation possible: dT = {dT} K <= 0.")

    # 4) Vapour velocity approaching the tubes
    U = G * x / rho_v

    # 5) Dimensionless groups (Bonneau nomenclature p. 2, liquid properties)
    Ja = cp_l * dT / h_lv
    H = Ja / Pr_l
    R = np.sqrt(rho_l * mu_l / (rho_v * mu_v))
    P = (Pr_l / Ja) * (rho_v / rho_l)
    Re_t = rho_l * U * Tube_OD / mu_l
    Fr = U**2 / (g * Tube_OD)
    
    if U <= 0.0:
        F = np.inf
    else:
        F = Pr_l / (Fr * Ja)

    # 6) Single-tube Nusselt number Nu_1 = h_1 * D / k_l
    D = Tube_OD
    if U <= 0.0 or F > 1000.0:
        # a) Nusselt Eq. 12 (gravity-controlled)
        Nu_1 = 0.728 * (rho_l**2 * g * D**3 * h_lv / (mu_l * k_l * dT))**0.25
    elif P > 10.0:
        # b) Fujii Eq. 70
        chi = 0.90 * (1.0 + 1.0 / (R * H))**(1.0 / 3.0)
        Nu_1 = (chi**4.0 + 0.276 * F)**0.25 * Re_t**0.5
    elif P <= F / 8.0:
        # c) Rose Eq. 90
        numerator = 0.9 * (1.0 + 1.0 / (R * H))**(1.0 / 3.0) + 0.728 * F**0.5
        denominator = (1.0 + 3.44 * F**0.5 + F)**0.25
        Nu_1 = (numerator / denominator) * Re_t**0.5
    else:
        # d) Rose Eq. 94
        pressure_term = (1.0 + 1.81 * P)**0.209
        numerator = 0.636 * pressure_term * (1.0 + 1.0 / (R * H))**(1.0 / 3.0) + 0.728 * F**0.5
        denominator = (1.0 + 3.51 * F**0.53 + F)**0.25
        Nu_1 = (numerator / denominator) * Re_t**0.5

    h_1 = Nu_1 * k_l / D

    # 7) Inundation factor (Kern Eq. 110)
    if 'N_rows' in geom:
        F_inund = geom['N_rows']**(-1.0 / 6.0)
    else:
        p_val = pitch_ratio * Tube_OD
        if tube_layout == 0:
            s_v = p_val
            s_c = p_val
        elif tube_layout == 45:
            s_v = np.sqrt(2.0) * p_val
            s_c = p_val / np.sqrt(2.0)
        else:
            s_v = np.sqrt(3.0) * p_val
            s_c = p_val / 2.0

        R_c = (Shell_ID - Tube_OD) / 2.0
        
        sum_numerator = 0.0
        sum_denominator = 0.0
        j = 0
        while j * s_c <= R_c:
            X = j * s_c
            n_j = np.floor(2.0 * np.sqrt(R_c**2 - X**2) / s_v) + 1.0
            if j == 0:
                sum_numerator = sum_numerator + n_j * n_j**(-1.0 / 6.0)
                sum_denominator = sum_denominator + n_j
            else:
                sum_numerator = sum_numerator + 2.0 * n_j * n_j**(-1.0 / 6.0)
                sum_denominator = sum_denominator + 2.0 * n_j
            j = j + 1
            
        F_inund = sum_numerator / sum_denominator

    # 8) Mean heat transfer coefficient over the bundle
    h = h_1 * F_inund

    return h


def STHE_shah_bundle_boiling_HTC(AS, geom, cond):
    """
    Local boiling heat transfer coefficient on the outer surface of tubes in a horizontal
    plain-tube bundle (shell-side flow boiling / kettle and flooded evaporators) [W/m2/K].

    Parameters
    ----------
    AS : CoolProp abstract state object
        Fluid property provider.
    geom : dict
        Geometry dictionary containing:
        - Tube_OD : tube outside diameter [m]
        - pitch_ratio : tube pitch / tube outside diameter [-]
        - Shell_ID : shell inside diameter [m]
        - central_spacing : central baffle spacing [m]
        - tube_layout : angle of tube layout [°] (0 or 45 = square, 60 = triangular)
        - n_parallel : optional, number of parallel shells [-] (default 1)
    cond : dict
        Condition dictionary containing:
        - m_dot : total mass flow rate on the shell side [kg/s]
        - p : saturation pressure [Pa]
        - x : local (cell mean) vapour quality [-]
        - q : local heat flux on the outer tube surface [W/m^2]

    Returns
    -------
    h : float
        Local boiling heat transfer coefficient on the nominal outer area [W/(m^2 K)].

    Reference
    ---------
    Shah, M.M. (2017), A correlation for heat transfer during boiling on bundles of horizontal
    plain and enhanced tubes, International Journal of Refrigeration 78, doi 10.1016/j.ijrefrig.
    2017.03.010, Eqs. 13-20 (modified correlation, section 3.2).
    Cooper, M.G. (1984), Advances in Heat Transfer 16, 157-239 (h_cooper).

    Validity Range (Shah 2017 Table 1)
    ---------------------------------
    p_red : 0.005 - 0.2866
    G : 0.15 - 1391 kg/(m^2 s)
    q : 1e3 - 1e6 W/m^2
    x : 0 - 0.98
    D : 3e-3 - 25.4e-3 m
    pitch_ratio : 1.17 - 2.0
    Re_L : 7 - 52117
    Database: 2173 data points, 12 fluids (water, ammonia, R-11, R-12, R-22, R-113,
    R-134a, R-236fa, R-410A, R-507A, isobutane, n-pentane; Shah 2017 Table 1).
    Mean absolute deviation 15.2 %.

    Definitions
    -----------
    Re_L = G * D / mu_l (all mass as liquid, tube OD).
    Bo and Fr are defined with G (labothappy: Kern cross-flow area at shell centre line).
    All fluid properties are evaluated at saturation temperature corresponding to pressure p.
    q is the heat flux on the outer tube surface.
    x is the local (cell mean) vapour quality.
    h_cooper is the Cooper (1984) original pool boiling heat transfer coefficient using M^-0.5
    and roughness term through the existing library function BPHE_cooper_boiling_HTC.
    """
    # 1) Shell-side mass flux and geometry definitions
    m_dot_shell = cond['m_dot'] / geom.get('n_parallel', 1)
    S = s_max_kern(geom['Tube_OD'], geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])
    G = m_dot_shell / S
    x = cond['x']
    D = geom['Tube_OD']
    g = 9.81

    # Check quality bounds
    if x < 0.0 or x > 1.0:
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: quality x = {x} outside [0, 1]. Clipped.")
        x = np.clip(x, 0.0, 1.0)

    # 2) Validate heat flux
    q = abs(cond['q'])
    if not np.isfinite(q) or q <= 0.0:
        raise ValueError(f"STHE_shah_bundle_boiling_HTC: heat flux q = {q} must be finite and > 0.")

    # 3) Compute Cooper pool boiling heat transfer coefficient before saturation properties
    cond_cooper = dict(cond)
    cond_cooper['q'] = q
    h_cooper = BPHE_cooper_boiling_HTC(AS, geom, cond_cooper)

    # 4) Obtain saturation properties and critical pressure
    p_critical = AS.p_critical()
    p_red = cond['p'] / p_critical
    sat_props = _sat_props(AS, cond['p'], cond)
    rho_l = sat_props['rho_l']
    mu_l = sat_props['mu_l']
    k_l = sat_props['k_l']
    Pr_l = sat_props['Pr_l']
    h_lv = sat_props['h_lv']

    # 5) Dimensionless numbers: Boiling number, Froude number, and modified Y parameter (Eq. 13)
    Bo = q / (G * h_lv)
    Fr = (G ** 2) / ((rho_l ** 2) * g * D)
    Y = Bo * (Fr ** 0.3)

    # 6) Plain tube factor
    F_pb = 1.0

    # Liquid-only Reynolds number (Eq. 20), also used in the validity check
    Re_L = G * D / mu_l

    # 7) High-Y branch (Eq. 14)
    if Y > 0.008:
        h = F_pb * h_cooper

    # 8) Low-Y branch (Y <= 0.008): maximum of Eq. 14, Eq. 20, Eq. 15, and Eq. 16
    else:
        # Eq. 14: h_14
        h_14 = F_pb * h_cooper

        # Liquid-only coefficient (Eq. 20)
        h_LT = 0.21 * (Re_L ** 0.62) * (Pr_l ** 0.4) * k_l / D

        # Eq. 15: phi0 components (Eqs. 17, 18, 19)
        phi_17 = 443.0 * (Bo ** 0.65) * F_pb
        phi_18 = 31.0 * (Bo ** 0.33) * F_pb
        phi_19 = 1.0
        phi0 = max(phi_17, phi_18, phi_19)
        h_15 = phi0 * h_LT

        # Eq. 16 (only if x > 0)
        if x > 0.0:
            Z = (((1.0 - x) / x) ** 0.8) * (p_red ** 0.4)
            phi = 2.3 / ((Z ** 0.08) * (Fr ** 0.22))
            h_16 = phi * h_LT
            h = max(h_14, h_15, h_16)
        else:
            h = max(h_14, h_15)

    # 9) Validity checks (Shah 2017 Table 1)
    if not (0.005 <= p_red <= 0.2866):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: reduced pressure p_red = {p_red:.4f} "
                      f"is outside validity range [0.005, 0.2866].")
    if not (0.15 <= G <= 1391.0):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: mass flux G = {G:.2f} kg/(m^2 s) "
                      f"is outside validity range [0.15, 1391].")
    if not (1e3 <= q <= 1e6):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: heat flux q = {q:.1f} W/m^2 "
                      f"is outside validity range [1e3, 1e6].")
    if not (0.0 <= x <= 0.98):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: quality x = {x:.3f} "
                      f"is outside validity range [0, 0.98].")
    if not (3e-3 <= D <= 25.4e-3):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: tube OD D = {D:.5f} m "
                      f"is outside validity range [3e-3, 25.4e-3].")
    if not (1.17 <= geom['pitch_ratio'] <= 2.0):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: pitch ratio = {geom['pitch_ratio']:.2f} "
                      f"is outside validity range [1.17, 2.0].")
    if not (7.0 <= Re_L <= 52117.0):
        warnings.warn(f"STHE_shah_bundle_boiling_HTC: liquid Reynolds number Re_L = {Re_L:.1f} "
                      f"is outside validity range [7, 52117].")

    # Return final heat transfer coefficient
    return h
