# -*- coding: utf-8 -*-
"""
Single-phase frictional pressure drop of spiral plate heat exchanger channels (rectangular, curved), (AS, geom, cond)
formalism. Returns the pressure drop over the full flow length geom['L'] [Pa].

geom : Dh, R_curv, L, ch_b or ch_w (channel gap), H (plate height), K (roughness, default 0), theta
cond : G [kg/m2/s] mass flux of ONE channel, mu, rho or D

Straight rectangular duct (Shah & London, laminar, Darcy):
    f Re = 96 (1 - 1.3553 a + 1.9467 a^2 - 1.7012 a^3 + 0.9564 a^4 - 0.2537 a^5),   a = min(gap, H) / max(gap, H)
    turbulent: Churchill with Dh and K. Curved: Ito factors of curved_pipe_DP_correlation_bank (base f Re above).
No entrance / exit / turn losses (they matter in a spiral exchanger: not included).
"""
import numpy as np
from labothappy.correlations.pressure_drop.pipe_DP_correlation_bank import friction_factor_churchill
from labothappy.correlations.pressure_drop.curved_pipe_DP_correlation_bank import _curved_sp_dp, G_GRAVITY
from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import _D, _L, _need_cond
from labothappy.correlations.convection.SPHE_HTC_correlation_bank import _AR


def _fRe_rect(geom):
    a = _AR(geom)
    return 96.0 * (1 - 1.3553 * a + 1.9467 * a ** 2 - 1.7012 * a ** 3 + 0.9564 * a ** 4 - 0.2537 * a ** 5)


def SPHE_curved_rect_SP_DP(AS, geom, cond):
    """Curved rectangular channel (Ito curvature factors on the rectangular straight-duct friction)."""
    return _curved_sp_dp(geom, cond, f_lam_Re=_fRe_rect(geom))


def SPHE_rect_SP_DP(AS, geom, cond):
    """Straight rectangular duct (no curvature effect)."""
    _need_cond(cond, 'G', 'mu')
    D, L = _D(geom), _L(geom)
    rho = cond['rho'] if 'rho' in cond else cond['D']
    G = cond['G']
    Re = G * D / cond['mu']
    f = _fRe_rect(geom) / Re if Re < 2300 else friction_factor_churchill(geom.get('K', 0.0), D, Re)
    return f * (L / D) * G ** 2 / (2 * rho) + G_GRAVITY * rho * L * np.sin(np.radians(geom.get('theta', 0.0)))


def SPHE_picon_nunez_SP_DP(AS, geom, cond):
    """Picon-Nunez et al. (2007) after Minton (1970), your single_phase_spiral_plate_DP (SI-converted constants).
    geom : Dh, L, H (plate height), ch_w (channel gap); cond : G of ONE channel (m_dot = G ch_w H), mu, rho or D,
    mu_w (default mu). Laminar below Re_cr = 12730 curv^0.2 (needs R_curv), turbulent above."""
    _need_cond(cond, 'G', 'mu')
    from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import _R_curv
    D, L, R = _D(geom), _L(geom), _R_curv(geom)
    H, ch_w = geom['H'], geom['ch_w']
    mu, mu_w = cond['mu'], cond.get('mu_w', cond['mu'])
    rho = cond['rho'] if 'rho' in cond else cond['D']
    m_dot = cond['G'] * ch_w * H
    Re = cond['G'] * D / mu
    s = rho / 1000.0
    t1 = m_dot / (ch_w * H)
    if Re <= 12730.0 * (D / (2.0 * R)) ** 0.2:
        t2 = 1.8516 * mu ** 0.5 / (ch_w + 0.003175)
        dp = 0.11583 * L / s * t1 * (t2 * (mu_w / mu) ** 0.17 * (H / m_dot) ** 0.5 + 1.5 + 4.8768 / L)
    else:
        t2 = 0.56315 * mu ** (1.0 / 3.0) / (ch_w + 0.003175)
        dp = 5.9308e-4 * L / s * t1 ** 2 * (t2 * (H / m_dot) ** (1.0 / 3.0) + 1.5 + 4.8768 / L)
    return dp + G_GRAVITY * rho * L * np.sin(np.radians(geom.get('theta', 0.0)))

