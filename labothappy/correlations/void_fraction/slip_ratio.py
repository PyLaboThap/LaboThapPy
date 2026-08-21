"""
Slip ratio correlations

@author: elise.neven@uliege.be
"""

def slip_ratio_zivi(rho_l, rho_g):
    """
    Zivi slip ratio correlation.
    
    Formula: S = (ρ_l / ρ_g)^(1/3)
    
    Parameters
    ----------
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    
    Returns
    -------
    float
        Slip ratio S = v_g / v_l [-]

    References
    ----------
    S. M. Zivi. “Estimation of steady-state steam void-fraction by means of  the principle of minimum entropy production”. 
    In: Journal of Heat Transfer 86.2 (1964), pp. 247–251.
    """
    return (rho_l / max(EPS, rho_g)) ** (1.0 / 3.0)