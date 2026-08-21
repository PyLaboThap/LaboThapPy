"""
Void fraction correlations for two phase flows

@author: elise.neven@uliege.be
"""
import numpy as np

from labothappy.correlations.properties.dimensionless import compute_reynolds, compute_weber
from labothappy.toolbox.solvers.zero_brent import zero_brent


EPS = 1e-12
G = 9.81 # Gravitational acceleration

#======================================================================
# HOMOGENEOUS
#======================================================================

def void_fraction_homogeneous(rho_l, rho_g, x):
    # Compute slip ratio
    s = 1

    # Compute void fraction
    # α = 1 / (1 + S · (ρ_g / ρ_l) · (1-x) / x)
    denominator = 1.0 + s * (rho_g / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha
    
#======================================================================
# ZIVI
#======================================================================

def void_fraction_zivi(rho_l, rho_g, x):
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
    # Compute slip ratio
    s_zivi = (rho_l / max(EPS, rho_g)) ** (1.0 / 3.0)

    # Compute void fraction
    # α = 1 / (1 + S · (ρ_g / ρ_l) · (1-x) / x)
    denominator = 1.0 + s_zivi * (rho_g / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha


#======================================================================
# FAUKSE
#======================================================================

def void_fraction_fauske(rho_l, rho_v, x):
    # Compute slip ratio
    s_faukse = (rho_l / max(EPS, rho_v)) ** (1.0 / 2.0)

    # Compute void fraction
    # α = 1 / (1 + S · (ρ_g / ρ_l) · (1-x) / x)
    denominator = 1.0 + s_faukse * (rho_v / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)
    
    return alpha

#======================================================================
# PREMOLI
#======================================================================
def void_fratcion_premoli(rho_l, rho_g, mu_l, sigma, d_hyd, G):
    Re_l = compute_reynolds(d_hyd = d_hyd, mu = mu_l, G = G)
    We_l = compute_weber(G=G, sigma=sigma, L=d_hyd, rho=rho_l)

    alpha_h = void_fraction_homogeneous(rho_l, rho_g)
    y = alpha_h/(1-alpha_h)
    F1 = 1.578*(Re_l**-0.19)*(rho_l/rho_g)**0.22
    F2 = 0.0273*We_l*(Re_l^-0.51)*(rho_l/rho_g)^-0.08
    s_premoli = 1 + F1*np.sqrt(max(0, (y/(1+y*F2))-y*F2))

    denominator = 1.0 + s_premoli * (rho_g / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha

#======================================================================
# LOCKHART MARTINELLI
#======================================================================

# Remi Dickes version
def void_fraction_lockhart_Martinelli_v1(rho_g, rho_l, x, mu_g, mu_l):
    X_tt = (((1-x)/x)**0.9)*((mu_l/mu_g)**0.1)*((rho_g/rho_l)**0.5)

    if X_tt <= 10:
        alpha = min(1, max(0, (1+X_tt**0.8)**-0.378))
    else:
        alpha = min(1, max(0, 0.823-0.157*np.log(X_tt))) #/!\ = ln dans la formule de Rémis Dickes!!

    return alpha


# Cioncolinu and Thome version
def void_fraction_lockart_Martinelli_v2(rho_g, rho_l, mu_g, mu_l, G, x, d_hyd):

    Re_g = compute_reynolds(G = G*x, d_hyd = d_hyd, mu=mu_g)
    Re_l = compute_reynolds(G=G*(1-x), d_hyd = d_hyd, mu=mu_l)

    # Friction coefficients
    # /!\ a ajouter au friction factor library!!
    if Re_g < 1500:
        f_g = 16/Re_g
    else:
        f_g = 0.046/Re_g**0.2

    if Re_l < 1500:
        f_l = 16/Re_l
    else:
        f_l = 0.046/Re_l**0.2

    # Lockhart and Martinelli function (1949)
    X = (f_l/f_g)**(0.5) * ((1 - x)/x) * (rho_l/rho_g)**(0.5)

    # Void fraction
    alpha = 1/(1 + 0.28*X**(0.71))

    return alpha 

#======================================================================
# HUGHMARK
#======================================================================

# /!\ ALLER RELIRE LA FORMULE CAR ELLE N4EST PAS LA MEME ENTRE LA FORMULE ECRITE DANS LE CODE ET DANS LA THESE DE REMIS

def void_fraction_hughmark(x, rho_l, rho_v, mu_l, mu_v, d_hyd, G):

    def residual_void_fraction_hughmark(alpha_guess, x, alpha_h, rho_v, mu_v, mu_l, d_hyd, G):
        Z = (((d_hyd*G)/(mu_l+alpha_guess*(mu_v-mu_l)))**(1/6))*((((1/9.81/d_hyd)*(G*x/(rho_v*alpha_h*(1-alpha_h)))**2)**(1/8)))
        ln_Z = np.log(Z)
        p1 = -0.010060658854755
        p2 = 0.155594796014726
        p3 = -0.870912508715887
        p4 = 2.167004115373165
        p5 = -2.224608445535130
        ln_Kh = p1*ln_Z**4 + p2*ln_Z**3 + p3*ln_Z**2 + p4*ln_Z + p5
        Kh = np.exp(ln_Kh)
        alpha_new = Kh*alpha_h
        res= (alpha_guess-alpha_new)
        return alpha_new, res


    x_min = 0.001
    x_max = 0.99
    if x > x_min and x < x_max:
        x1 = x
    elif x < x_min:
        x = x_min
    elif x > x_max:
        x = x_max

    alpha_h = void_fraction_homogeneous(rho_l = rho_l, rho_g = rho_v, x =x)
    # OK jusque ici

    def f(alpha_guess):
            _, res = residual_void_fraction_hughmark(
                alpha_guess, x, alpha_h, rho_v, mu_v, mu_l, d_hyd, G
            )
            return res

    machep = np.finfo(float).eps
    t = 1e-10
    tol_f = 1e-8

    alpha, res_alpha = zero_brent(1e-6, 1 - 1e-6, machep, t, f, tol_f)

    return alpha


#======================================================================
# GRAHAM
#======================================================================
# /!\ for condensation only!
def void_fraction_graham(x, rho_v, G, D):

    Ft = (((x**3)*(G**2))/(9.81*(rho_v**2)*D*(1-q)))**0.5

    if Ft > 0.01032:
        alpha = min(1,max(0,1-np.exp(-1-0.3*np.log(Ft)-0.0328*(np.log(Ft))**2)))
    else:
        alpha = 0

    return alpha


#======================================================================
# AMRMAND TRESCHEV (K epsilon method)
#======================================================================

# /!\ Vérfier la plage de validité

def void_fraction_armand_treschev(x, rho_v, rho_l):
    alpha_h = void_fraction_homogeneous(rho_l = rho_l, rho_g=rho_v, x=x)
    K = 0.833 + 0.167 * x

    alpha = K*alpha_h
    return alpha


#======================================================================
# BANKOFF (K epsilon method)
#======================================================================

# /!\ Vérfier la plage de validité

def void_fraction_bankoff(P, x, rho_l, rho_v):

    alpha_h = void_fraction_homogeneous(rho_l=rho_l, rho_g=rho_v, x=x)
    K = 0.71 + 0.0145 * (P / 1e6)  # pressure term in MPa

    alpha = K*alpha_h

    return alpha


#======================================================================
# AXELSSON (Drift-flux correlations)
#======================================================================

def void_fraction_rouhani_axelsson(x, rho_l, rho_g, G, sigma):
    J_l = (1 - x) * G / rho_l
    J_g = x * G / rho_g

    C0 = 1 + 0.2 * (1 - x)
    V_drift = 1.18 * ((g * sigma * (rho_l - rho_g)) / (rho_l ** 2)) ** 0.25

    alpha = J_g / (C0 * (J_l + J_g) + V_drift)

    return alpha

#======================================================================
# DIX or CHEXAL model (Drift-flux correlations)
#======================================================================

def void_fraction_dix(x, rho_l, rho_g, G, g, sigma):
    J_l = (1 - x) * G / rho_l
    J_g = x * G / rho_g

    n = (rho_g / rho_l) ** 0.1
    C0 = J_g / (J_l + J_g) * (1 + (J_l / J_g) ** n)
    V_drift = 2.9 * ((g * sigma * (rho_l - rho_g)) / (rho_l ** 2)) ** 0.25

    alpha = J_g / (C0 * (J_l + J_g) + V_drift)

    return alpha

#======================================================================
# woldesemayat_ghajar model (Drift-flux correlations)
#======================================================================

def void_fraction_woldesemayat_ghajar(x, rho_l, rho_g, G, g, sigma, D, P, theta=np.pi/2):
    """
    Void fraction correlation by Woldesemayat and Ghajar (2007).
    
    Extends the DIX distribution parameter C0 with a drift velocity that
    accounts for channel diameter, inclination angle, and pressure.
    """

    J_l = (1 - x) * G / rho_l
    J_g = x * G / rho_g
    
    n = (rho_g / rho_l) ** 0.1
    C0 = J_g / (J_l + J_g) * (1 + (J_l / J_g) ** n)

    V_drift = (2.9 * ((g * sigma * D * (1 + np.cos(theta)) * (rho_l - rho_g)) / (rho_l ** 2)) ** 0.25
                   * (1.22 + 1.22 * np.sin(theta)) ** (101325 / P))

    alpha = J_g / (C0 * (J_l + J_g) + V_drift)

    return alpha

#======================================================================
# cioncolini_thome model
#======================================================================

def void_fraction_cioncolini_thome(x, rho_l, rho_g):
    """
    Void fraction correlation by Cioncolini and Thome (2012).

    Annular two-phase flow correlation, best suited for:
    0 < x < 1, 1e-3 < rho_g/rho_l < 1, 0.7 < eps < 1.

    Parameters
    ----------
    x : float
        Vapor quality [-]
    rho_l, rho_g : float
        Liquid and gas density [kg/m^3]

    Returns
    -------
    float
        Void fraction [-]
    """
    rho_ratio = rho_g / rho_l

    h = -2.129 + 3.129 * rho_ratio ** (-0.2186)
    n = 0.3487 + 0.6513 * rho_ratio ** 0.515

    return (h * x ** n) / (1 + (h - 1) * x ** n)
