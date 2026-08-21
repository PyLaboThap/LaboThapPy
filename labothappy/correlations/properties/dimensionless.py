
import numpy as np

EPS = 1e-12

# def compute_reynolds(rho, v, d_hyd, mu):
#     """
#     Compute Reynolds number.
    
#     Re = ρ·v·d_hyd / η
    
#     Parameters
#     ----------
#     rho : float
#         Density [kg/m³]
#     v : float
#         Velocity [m/s]
#     d_hyd : float
#         Hydraulic diameter [m]
#     mu : float
#         Dynamic viscosity [Pa·s]
    
#     Returns
#     -------
#     float
#         Reynolds number [-]
#     """
#     return max(1.0, rho * v * d_hyd / max(EPS, mu))

# Changed to be able to take as inputs also G (to check if necessary!!)
def compute_reynolds(d_hyd, mu, rho=None, v=None, G=None):
    """
    Compute Reynolds number.

    Re = ρ·v·d_hyd / μ  =  G·d_hyd / μ   (since G = ρ·v)

    Parameters
    ----------
    d_hyd : float
        Hydraulic diameter [m]
    mu : float
        Dynamic viscosity [Pa·s]
    rho : float, optional
        Density [kg/m³] (used with v if G is not given)
    v : float, optional
        Velocity [m/s] (used with rho if G is not given)
    G : float, optional
        Mass flux [kg/(m²·s)] (overrides rho, v if given)

    Returns
    -------
    float
        Reynolds number [-]
    """
    if G is None:
        if rho is None or v is None:
            raise ValueError("Provide either G, or both rho and v.")
        G = rho * v

    return max(1.0, G * d_hyd / max(EPS, mu))



def compute_froude(G, L, rho):
    """
    Compute Froude number.

    Fr = G / (sqrt(g·d_hyd)·ρ)

    Parameters
    ----------
    G : float
        Mass flux [kg/(m²·s)]
    L : float
        Characteristic length [m]
    rho : float
        Density [kg/m³]
    
    Returns
    -------
    float
        Froude number [-]
    """
    g = 9.81  # gravitational acceleration [m/s²]
    Fr = G / (np.sqrt(g * L) * rho)
    return Fr

def compute_weber(G, L, sigma, rho):
    """
    Compute Weber number.

    We = G²·L / (σ·ρ)

    Parameters
    ----------
    G : float
        Mass flux [kg/(m²·s)]
    L : float
        Characteristic length [m]
    sigma : float
        Surface tension [N/m]
    rho : float
        Density [kg/m³]
    
    Returns
    -------
    float
        Weber number [-]
    """
    return (G**2 * L) / (sigma * rho)

        