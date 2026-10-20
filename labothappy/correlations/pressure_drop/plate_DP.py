import numpy as np

def dp_martin_plate(G_ch, rho, mu, D_h, chevron_angle, L):
    """
    Single-phase frictional pressure drop in a chevron plate heat exchanger
    channel, after H. Martin, VDI Heat Atlas, chapter N6.

    Inputs
    ------
    G_ch          : mass flux in ONE channel [kg/m^2/s] (= m_dot_side / (n_channels_side*A_cs))
    rho           : density [kg/m^3]
    mu            : dynamic viscosity at bulk temperature [Pa.s]
    D_h           : hydraulic diameter of one channel, 2 * b / enlargement_factor [m]
    chevron_angle : corrugation angle measured FROM THE MAIN FLOW DIRECTION [rad]
    L             : flow length along the main flow direction [m]
                    (L_p for the whole plate, or w_i * L_p for a zone of a
                    moving boundary model)

    Output
    -------
    dp : frictional pressure drop over the length L [Pa]
         (port/distribution losses NOT included)
    """
    # Dimensionless number
    Re = G_ch * D_h / mu

    # Friction factors for phi = 0 (zeta0) and phi = 90 deg (zeta1_0)
    if Re < 2000:
        zeta0 = 64 / Re
        zeta1_0 = 597 / Re + 3.85
    else:
        zeta0 = (1.8 * np.log(Re) - 1.5) ** (-2)
        zeta1_0 = 39 / Re ** 0.289

    # Martin's fitting constants (a, b, c in the original paper)
    a_M, b_M, c_M = 3.8, 0.18, 0.36
    zeta1 = a_M * zeta1_0

    # Martin's equation 18: combined (Darcy-type) friction factor
    rhs = (np.cos(chevron_angle) / np.sqrt(b_M * np.tan(chevron_angle) + c_M * np.sin(chevron_angle) + zeta0 / np.cos(chevron_angle))
           + (1 - np.cos(chevron_angle)) / np.sqrt(zeta1))
    zeta = 1 / rhs**2

    # Pressure drop (Darcy form)
    # Equivalent to ACHP's: Hg * mu**2 * L / (rho * D_h**3), with Hg = zeta * Re**2 / 2
    dp = zeta * (L / D_h) * G_ch**2 / (2 * rho)
    return dp


def _zeta_martin(Re, chevron_angle):
    """
    Martin (VDI Heat Atlas N6) Darcy-type friction factor for a chevron plate
    channel. chevron_angle measured FROM THE MAIN FLOW DIRECTION [rad].
    """
    if Re < 2000:
        zeta0 = 64 / Re
        zeta1_0 = 597 / Re + 3.85
    else:
        zeta0 = (1.8 * np.log(Re) - 1.5) ** (-2)
        zeta1_0 = 39 / Re ** 0.289

    a_M, b_M, c_M = 3.8, 0.18, 0.36
    zeta1 = a_M * zeta1_0

    rhs = (np.cos(chevron_angle) / np.sqrt(b_M * np.tan(chevron_angle) + c_M * np.sin(chevron_angle) + zeta0 / np.cos(chevron_angle))
           + (1 - np.cos(chevron_angle)) / np.sqrt(zeta1))
    return 1 / rhs**2


def _fanning_tube(Re):
    """
    Smooth circular tube Fanning friction factor as used in ACHP's
    Lockhart-Martinelli implementation (laminar 16/Re, turbulent Blasius-type
    0.046 Re^-0.2, linear blend for 1000 < Re < 2000).
    """
    if Re < 1000:
        return 16 / Re
    elif Re > 2000:
        return 0.046 / Re ** 0.2
    else:
        w = (Re - 1000) / (2000 - 1000)
        return (1 - w) * 16 / Re + w * 0.046 / Re ** 0.2


def dp_lockhart_martinelli(G, x, rho_l, rho_v, mu_l, mu_v, D_h, L,
                                  C=None, chevron_angle=None):
    """
    Two-phase FRICTIONAL pressure drop with the Lockhart-Martinelli separated
    flow model, in the analytical form of Chisholm (1967).

    Original validity: adiabatic two-phase flow in horizontal, smooth circular
    tubes (Lockhart & Martinelli, 1949). Applying it to a plate heat exchanger
    is an extrapolation.

    Inputs
    ------
    G             : total mass flux in ONE channel/tube [kg/m^2/s]
    x             : vapour quality [-] (use the zone-average quality for a zone)
    rho_l, rho_v  : saturated liquid / vapour density [kg/m^3]
    mu_l, mu_v    : saturated liquid / vapour dynamic viscosity [Pa.s]
    D_h           : hydraulic diameter [m] (tube ID, or 2*b/enlargement_factor for a PHE)
    L             : flow length [m] (L_p, or w_i * L_p for a moving boundary zone)
    C             : Chisholm constant [-], optional. If None, chosen from the
                    laminar/turbulent regime of each phase (5, 10, 12 or 20)
    chevron_angle : corrugation angle FROM THE MAIN FLOW DIRECTION [rad], optional.
                    If None  -> smooth tube friction factors (original method).
                    If given -> Martin's PHE friction factor for each phase
                                (adaptation, C not refitted for plates).

    Output
    -------
    dp : frictional two-phase pressure drop over the length L [Pa]
         (acceleration, gravity and port losses NOT included)
    """
    x = float(x)

    # 1. Reynolds number of each phase flowing ALONE at its actual flow rate
    Re_l = G * (1 - x) * D_h / mu_l
    Re_v = G * x * D_h / mu_v

    # 2. Single-phase frictional pressure gradients [Pa/m]
    def dpdz_single(Re, G_phase, rho):
        if G_phase <= 0:
            return 0.0
        if chevron_angle is None:
            f_F = _fanning_tube(Re)                      # Fanning
            return 2 * f_F * G_phase**2 / (rho * D_h)
        else:
            zeta = _zeta_martin(Re, chevron_angle)       # Darcy
            return zeta * G_phase**2 / (2 * rho * D_h)

    dpdz_l = dpdz_single(Re_l, G * (1 - x), rho_l)
    dpdz_v = dpdz_single(Re_v, G * x, rho_v)

    # 3. Chisholm constant from the flow regime of each phase
    #    (1500 used as transition Re to ensure continuity, as in ACHP)
    if C is None:
        if Re_l > 1500 and Re_v > 1500:
            C = 20.0    # turbulent liquid - turbulent vapour
        elif Re_l < 1500 and Re_v > 1500:
            C = 12.0    # laminar liquid - turbulent vapour
        elif Re_l > 1500 and Re_v < 1500:
            C = 10.0    # turbulent liquid - laminar vapour
        else:
            C = 5.0     # laminar liquid - laminar vapour

    # 4. Two-phase gradient (Chisholm form, identical to phi_l^2*dpdz_l = phi_v^2*dpdz_v)
    dpdz = dpdz_l + C * np.sqrt(dpdz_l * dpdz_v) + dpdz_v

    dp = dpdz * L
    return dp