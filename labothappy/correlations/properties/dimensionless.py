
import numpy as np

EPS = 1e-12

def compute_reynolds(rho, v, d_hyd, mu):
    """
    Compute Reynolds number.
    
    Re = ρ·v·d_hyd / η
    
    Parameters
    ----------
    rho : float
        Density [kg/m³]
    v : float
        Velocity [m/s]
    d_hyd : float
        Hydraulic diameter [m]
    mu : float
        Dynamic viscosity [Pa·s]
    
    Returns
    -------
    float
        Reynolds number [-]
    """
    return max(1.0, rho * v * d_hyd / max(EPS, mu))


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

        