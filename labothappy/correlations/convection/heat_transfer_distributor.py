"""
htc_distributor.py

Single entry point to compute a heat transfer coefficient from any correlation,
for any heat exchanger geometry:

    htc = compute(AS, m_dot, geom, correlation, q_flux=None, T_wall=None)

The model only provides the fluid state (AS, at the mean state of the zone),
the mass flow rate and the geometry of one side. Each adapter below translates
these into the specific inputs of its correlation, including the translation
from one geometry type to another when a correlation is used outside the
geometry it was made for.

author: Elise
"""

import CoolProp.CoolProp as CP

from labothappy.correlations.convection.plate_htc import htc_martin_plate_1phase, htc_cooper_pool_boiling, htc_longo_condensation

# Correlations whose htc depends on the heat flux
NEEDS_Q_FLUX = {'Cooper'}

# ===============================================================
# Entry point
# ===============================================================

def compute_htc(AS, m_dot, geom, correlation, q_flux=None, T_wall=None):
    """
    Heat transfer coefficient [W/m^2/K].

    AS          : CoolProp AbstractState, set at the MEAN state of the zone
    m_dot       : total mass flow rate on this side [kg/s]
    geom        : geometry of this side (PlateGeometry, TubeGeometry)
    correlation : name of the correlation (key of _ADAPTERS)
    q_flux      : heat flux through the wall [W/m^2], for boiling correlations
    T_wall      : wall temperature [K], optional
    """
    try:
        adapter = HTC_CORRELATIONS[correlation]
    except KeyError:
        raise ValueError(f"Unknown htc correlation '{correlation}'. "
                         f"Available: {list(HTC_CORRELATIONS)}") from None

    # Adapters may update AS (saturated or wall properties): save the zone
    # state and always restore it, so the caller's AS is never modified.
    p, h = AS.p(), AS.hmass()
    try:
        return adapter(AS, m_dot, geom, q_flux, T_wall)
    finally:
        AS.update(CP.HmassP_INPUTS, h, p)

# ===============================================================
# Adapters: all have the same signature (AS, m_dot, geom, q_flux, T_wall)
# ===============================================================

# ---------------- Single-phase ----------------

def _adapter_martin(AS, m_dot, geom, q_flux, T_wall):
    """Martin (VDI): chevron plates only. q_flux is not used."""
    if geom['chevron_angle'] is None:
        raise ValueError("Martin needs chevron_angle")

    # Single-phase properties at the current state
    mu = AS.viscosity()
    cp = AS.cpmass()
    k = AS.conductivity()

    # Wall viscosity (optional)
    if T_wall is None:
        mu_wall = None
    else:
        AS.update(CP.PT_INPUTS, AS.p(), T_wall)
        mu_wall = AS.viscosity()

    # Mass flux through ONE channel
    G_ch = m_dot / (geom['A_cs_channel'] * geom['n_channels'])

    return htc_martin_plate_1phase(
        G_ch=G_ch, mu=mu, cp=cp, k=k,
        D_h=geom['D_h'], chevron_angle=geom['chevron_angle'], mu_wall=mu_wall)

# ---------------- Evaporation ----------------

def _adapter_cooper(AS, m_dot, geom, q_flux, T_wall):
    """Cooper: nucleate pool boiling, independent of the geometry."""
    if q_flux is None:
            raise ValueError("Cooper needs q_flux.")

    return htc_cooper_pool_boiling(
        p=AS.p(), p_crit=AS.p_critical(), q_flux=q_flux,
        M_molar=AS.molar_mass(), roughness=getattr(geom, 'roughness', None))


# ---------------- Condensation ----------------

def _adapter_longo_cond(AS, m_dot, geom, q_flux, T_wall):
    """Longo: condensation in brazed plates. Uses D_h = 2b, without Phi."""
    if geom['corrugation_amplitude'] is None:
      raise ValueError("Longo condensation needs corrugation_amplitude.")

    x = AS.Q()                                   # read BEFORE _sat_props updates AS = x_avg over the cell
    p = AS.p()                                   # read before updating AS

    # Saturated liquid (bubble point for a mixture)
    AS.update(CP.PQ_INPUTS, p, 0.0)
    rho_L = AS.rhomass()
    mu_L = AS.viscosity()
    cp_L = AS.cpmass()
    k_L = AS.conductivity()

    # Saturated vapour (dew point for a mixture)
    AS.update(CP.PQ_INPUTS, p, 1.0)
    rho_V = AS.rhomass()
    # Mass flux through ONE channel
    G_ch = m_dot / (geom['A_cs_channel'] * geom['n_channels'])
    htc = htc_longo_condensation(x=x, G=G_ch, D_h=geom['D_h'], rho_L=rho_L, rho_V=rho_V, mu_L=mu_L, cp_L=cp_L, k_L=k_L)

    return htc


HTC_CORRELATIONS = {
    'Martin':     _adapter_martin,
    'Cooper':     _adapter_cooper,
    'Longo_cond': _adapter_longo_cond,
}
