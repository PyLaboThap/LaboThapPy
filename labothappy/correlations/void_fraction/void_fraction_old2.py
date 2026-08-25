"""
Void fraction correlations for two-phase flow.

Supports multiple slip ratio models:
- Homogeneous (S=1, default)
- Zivi (S = (ρ_l/ρ_g)^(1/3))
- Custom slip ratio
"""

EPS = 1e-12

#=======================================================================
# SLIP RATIO MODELS
#+======================================================================


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


def slip_ratio_homogeneous(rho_l, rho_g):
    """
    Homogeneous slip ratio (no slip).
    
    Formula: S = 1.0
    
    Both phases move at same velocity.
    """
    return 1.0

#=======================================================================
# VOID FRACTION COMPUTATION
#======================================================================

def compute_void_fraction(x, rho_l, rho_g, slip_model=None):
    """
    Compute void fraction using separated flow model with slip ratio.
    
    **Formula:**
    
        α = 1 / (1 + S · (ρ_g / ρ_l) · (1-x) / x)
    
    where S is the slip ratio (depends on model).
    
    Parameters
    ----------
    x : float
        Vapor quality [-]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    slip_model : str or float, optional
        Slip ratio model or custom value.
        
        - None or 'Homogeneous': S = 1.0 (DEFAULT)
        - 'Zivi': S = (ρ_l/ρ_g)^(1/3)
        - float: Use custom slip ratio value directly
    
    Returns
    -------
    alpha : float
        Void fraction (volume fraction of vapor) [0, 1]

    References
    ----------

    Examples
    --------
    >>> # Default: homogeneous model
    >>> alpha = compute_void_fraction(x=0.5, rho_l=500, rho_g=50)
    
    >>> # Zivi slip ratio
    >>> alpha = compute_void_fraction(x=0.5, rho_l=500, rho_g=50, 
    ...                                slip_model='Zivi')
    
    >>> # Custom slip ratio
    >>> alpha = compute_void_fraction(x=0.5, rho_l=500, rho_g=50, 
    ...                                slip_model=1.5)
    
    """
    
    # Safeguard quality
    x = max(EPS, min(1.0 - EPS, x))
    
    # ====================================================================
    # Determine slip ratio
    # ====================================================================
    
    if slip_model is None or slip_model == 'Homogeneous':
        S = slip_ratio_homogeneous(rho_l, rho_g)
    
    elif slip_model == 'Zivi':
        S = slip_ratio_zivi(rho_l, rho_g)
    
    elif isinstance(slip_model, (int, float)):
        # Custom slip ratio value
        S = float(slip_model)
    
    else:
        raise ValueError(
            f"Unknown slip_model: {slip_model}. "
            "Use None, 'Homogeneous', 'Zivi', or a numeric value."
        )
    
    # ====================================================================
    # Compute void fraction
    # ====================================================================
    
    # α = 1 / (1 + S · (ρ_g / ρ_l) · (1-x) / x)
    denominator = 1.0 + S * (rho_g / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)
    
    return alpha