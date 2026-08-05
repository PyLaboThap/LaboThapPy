
from correlations.properties.void_fraction import compute_void_fraction
EPS = 1e-12

def compute_two_phase_density(x, rho_l, rho_g, alpha=None):
    """
    Compute two-phase density.
    
    If alpha is provided: ρ_mix = 1 / (α/ρ_g + (1-α)/ρ_l)
    Otherwise: use homogeneous model
    
    Parameters
    ----------
    x : float
        Vapor quality [-]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    alpha : float, optional
        Void fraction [-]. If None, computed from homogeneous model.
    
    Returns
    -------
    float
        Mixture density [kg/m³]
    """
    x = max(EPS, min(1.0 - EPS, x))
    
    if alpha is None:
        alpha = compute_void_fraction(x, rho_l, rho_g, slip_model='Homogeneous')
    
    # ρ_mix = 1 / (α/ρ_g + (1-α)/ρ_l)
    return 1.0 / (alpha / rho_g + (1.0 - alpha) / rho_l)