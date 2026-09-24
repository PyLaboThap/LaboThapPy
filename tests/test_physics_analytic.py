# -*- coding: utf-8 -*-
"""
Analytic checks on the models whose answer is known in closed form.

Every assertion here is a statement of physics, not a recorded number. Nothing
in this file needs updating when a correlation is retuned or a solver is
replaced -- an energy balance is an energy balance. That is the property that
makes these tests worth having in a modelling library: they constrain the
result without freezing the implementation.

Four kinds of check are used:

* **Definitional identities.** Isentropic efficiency is a definition, so the
  model must reproduce it exactly from its own outputs.
* **Conservation.** Mass and energy across a component must balance.
* **The second law.** Entropy cannot fall across an adiabatic component; heat
  cannot cross from cold to hot.
* **Limiting cases.** At an efficiency of one, compression is isentropic. At an
  effectiveness of one, the heat exchanger reaches its thermodynamic maximum.

Fluids and conditions are chosen to sit well away from the saturation dome and
from CoolProp's range limits, so a failure means the model is wrong rather than
the property call.
"""

from __future__ import annotations

import math

import pytest

import CoolProp.CoolProp as CP

from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.expander.expander_csteff import ExpanderCstEff
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
from labothappy.component.pump.pump_csteff import PumpCstEff
from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic

# Superheated R134a: 20 K above saturation at 3 bar, compressed to 10 bar.
VAPOUR = dict(fluid="R134a", T_su=300.0, P_su=3.0e5, P_ex=10.0e5, m_dot=0.5)

# Subcooled water: 20 degC, pumped from 1 bar to 20 bar.
LIQUID = dict(fluid="Water", T_su=293.15, P_su=1.0e5, P_ex=20.0e5, m_dot=1.0)

REL_TOL = 1e-9


def _solved(component, **inputs):
    """Configure, solve, and assert the component believes it succeeded."""
    component.print_flag = 0
    component.set_inputs(**inputs)
    component.solve()
    assert component.solved is True, (
        "%s reported solved=False for a well-posed case; the model swallows "
        "the real exception inside solve()" % type(component).__name__
    )
    return component


def _h_isentropic(fluid, p_ex, s_su):
    state = CP.AbstractState("HEOS", fluid)
    state.update(CP.PSmass_INPUTS, p_ex, s_su)
    return state.hmass()


# ---------------------------------------------------------------------------
# Compressor
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("eta_is", [0.5, 0.7, 0.9])
def test_compressor_reproduces_its_own_isentropic_efficiency(eta_is):
    """eta_is = (h_ex_is - h_su) / (h_ex - h_su), by definition.

    Recomputing the isentropic enthalpy independently from CoolProp and
    inverting the definition must return the parameter that went in. This
    catches a sign error, a reciprocal, or an efficiency applied to the wrong
    enthalpy difference -- the three ways this formula is usually got wrong.
    """
    compressor = CompressorCstEff()
    compressor.set_parameters(eta_is=eta_is)
    _solved(compressor, P_su=VAPOUR["P_su"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_ex"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])

    h_ex_is = _h_isentropic(VAPOUR["fluid"], compressor.ex.p, compressor.su.s)
    recovered = (h_ex_is - compressor.su.h) / (compressor.ex.h - compressor.su.h)
    assert recovered == pytest.approx(eta_is, rel=1e-6)


def test_compressor_at_unit_efficiency_is_isentropic():
    """The limiting case: a perfect compressor produces no entropy."""
    compressor = CompressorCstEff()
    compressor.set_parameters(eta_is=1.0)
    _solved(compressor, P_su=VAPOUR["P_su"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_ex"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert compressor.ex.s == pytest.approx(compressor.su.s, rel=1e-6)


def test_compressor_is_irreversible_below_unit_efficiency():
    """The second law: a real compressor raises entropy and temperature."""
    compressor = CompressorCstEff()
    compressor.set_parameters(eta_is=0.7)
    _solved(compressor, P_su=VAPOUR["P_su"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_ex"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert compressor.ex.s > compressor.su.s
    assert compressor.ex.T > compressor.su.T
    assert compressor.ex.h > compressor.su.h


def test_compressor_power_is_mass_flow_times_specific_work():
    """Energy conservation on the work connector, and mass conservation."""
    compressor = CompressorCstEff()
    compressor.set_parameters(eta_is=0.7)
    _solved(compressor, P_su=VAPOUR["P_su"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_ex"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])

    specific_work = compressor.ex.h - compressor.su.h
    assert compressor.W.W_dot == pytest.approx(
        compressor.su.m_dot * specific_work, rel=REL_TOL
    )
    assert compressor.ex.m_dot == pytest.approx(compressor.su.m_dot, rel=REL_TOL)
    assert compressor.ex.fluid == compressor.su.fluid


def test_compressor_work_falls_as_efficiency_rises():
    """Monotonicity: a better compressor costs less work for the same duty."""
    work = []
    for eta_is in (0.5, 0.7, 0.9):
        compressor = CompressorCstEff()
        compressor.set_parameters(eta_is=eta_is)
        _solved(compressor, P_su=VAPOUR["P_su"], T_su=VAPOUR["T_su"],
                P_ex=VAPOUR["P_ex"], fluid=VAPOUR["fluid"],
                m_dot=VAPOUR["m_dot"])
        work.append(compressor.W.W_dot)
    assert work == sorted(work, reverse=True)


# ---------------------------------------------------------------------------
# Expander
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("eta_is", [0.5, 0.7, 0.9])
def test_expander_reproduces_its_own_isentropic_efficiency(eta_is):
    """eta_is = (h_su - h_ex) / (h_su - h_ex_is) for expansion.

    Note the inverted form relative to the compressor: the same parameter name
    means the reciprocal ratio. Testing both against their own definitions is
    what stops the two models drifting into using the same formula.
    """
    expander = ExpanderCstEff()
    expander.set_parameters(eta_is=eta_is)
    _solved(expander, P_su=VAPOUR["P_ex"], T_su=360.0,
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])

    h_ex_is = _h_isentropic(VAPOUR["fluid"], expander.ex.p, expander.su.s)
    recovered = (expander.su.h - expander.ex.h) / (expander.su.h - h_ex_is)
    assert recovered == pytest.approx(eta_is, rel=1e-6)


def test_expander_at_unit_efficiency_is_isentropic():
    expander = ExpanderCstEff()
    expander.set_parameters(eta_is=1.0)
    _solved(expander, P_su=VAPOUR["P_ex"], T_su=360.0,
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert expander.ex.s == pytest.approx(expander.su.s, rel=1e-6)


def test_expander_produces_work_and_entropy():
    """An expander must deliver positive power and cannot destroy entropy."""
    expander = ExpanderCstEff()
    expander.set_parameters(eta_is=0.7)
    _solved(expander, P_su=VAPOUR["P_ex"], T_su=360.0,
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert expander.ex.h < expander.su.h
    assert expander.ex.s > expander.su.s
    assert expander.W.W_dot > 0


def test_expander_work_rises_with_efficiency():
    work = []
    for eta_is in (0.5, 0.7, 0.9):
        expander = ExpanderCstEff()
        expander.set_parameters(eta_is=eta_is)
        _solved(expander, P_su=VAPOUR["P_ex"], T_su=360.0,
                P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"],
                m_dot=VAPOUR["m_dot"])
        work.append(expander.W.W_dot)
    assert work == sorted(work)


def test_isentropic_compression_then_expansion_returns_to_the_start():
    """A round trip at unit efficiency must recover the original state.

    This is the strongest available check that the two models are consistent
    with each other rather than merely each self-consistent: an error common to
    both formulas would survive the individual efficiency tests but not this.
    """
    compressor = CompressorCstEff()
    compressor.set_parameters(eta_is=1.0)
    _solved(compressor, P_su=VAPOUR["P_su"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_ex"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])

    expander = ExpanderCstEff()
    expander.set_parameters(eta_is=1.0)
    _solved(expander, P_su=compressor.ex.p, T_su=compressor.ex.T,
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])

    assert expander.ex.h == pytest.approx(compressor.su.h, rel=1e-5)
    assert expander.ex.T == pytest.approx(compressor.su.T, rel=1e-5)
    assert compressor.W.W_dot == pytest.approx(expander.W.W_dot, rel=1e-5)


# ---------------------------------------------------------------------------
# Pump
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("eta_is", [0.4, 0.6, 0.8])
def test_pump_reproduces_its_own_isentropic_efficiency(eta_is):
    pump = PumpCstEff()
    pump.set_parameters(eta_is=eta_is)
    _solved(pump, P_su=LIQUID["P_su"], T_su=LIQUID["T_su"],
            P_ex=LIQUID["P_ex"], fluid=LIQUID["fluid"], m_dot=LIQUID["m_dot"])

    h_ex_is = _h_isentropic(LIQUID["fluid"], pump.ex.p, pump.su.s)
    recovered = (h_ex_is - pump.su.h) / (pump.ex.h - pump.su.h)
    assert recovered == pytest.approx(eta_is, rel=1e-6)


def test_pump_work_is_close_to_the_incompressible_estimate():
    """v*dP/eta is the textbook approximation for a liquid pump.

    Water at 20 degC is very nearly incompressible, so the exact model must
    land within a few per cent of the hand calculation. This is the check that
    would catch a units error -- a factor of 1000 from kJ/kg, or a pressure
    read in bar.
    """
    pump = PumpCstEff()
    eta_is = 0.6
    pump.set_parameters(eta_is=eta_is)
    _solved(pump, P_su=LIQUID["P_su"], T_su=LIQUID["T_su"],
            P_ex=LIQUID["P_ex"], fluid=LIQUID["fluid"], m_dot=LIQUID["m_dot"])

    density = CP.PropsSI("D", "T", LIQUID["T_su"], "P", LIQUID["P_su"],
                         LIQUID["fluid"])
    estimate = (LIQUID["P_ex"] - LIQUID["P_su"]) / density / eta_is
    actual = pump.ex.h - pump.su.h
    assert actual == pytest.approx(estimate, rel=0.05)


# ---------------------------------------------------------------------------
# Isenthalpic valve
# ---------------------------------------------------------------------------


def test_valve_is_exactly_isenthalpic():
    """The defining property, and it should hold to the last bit.

    Enthalpy is copied rather than computed, so anything other than exact
    equality means the value made a round trip through another property.
    """
    valve = ValveIsenthalpic()
    _solved(valve, P_su=VAPOUR["P_ex"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert valve.ex.h == valve.su.h


def test_throttling_generates_entropy():
    """Throttling is irreversible: pressure falls, entropy rises."""
    valve = ValveIsenthalpic()
    _solved(valve, P_su=VAPOUR["P_ex"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert valve.ex.p < valve.su.p
    assert valve.ex.s > valve.su.s


def test_valve_conserves_mass():
    valve = ValveIsenthalpic()
    _solved(valve, P_su=VAPOUR["P_ex"], T_su=VAPOUR["T_su"],
            P_ex=VAPOUR["P_su"], fluid=VAPOUR["fluid"], m_dot=VAPOUR["m_dot"])
    assert valve.ex.m_dot == pytest.approx(valve.su.m_dot, rel=REL_TOL)


# ---------------------------------------------------------------------------
# Constant-effectiveness heat exchanger
# ---------------------------------------------------------------------------


HOT = dict(fluid_H="Water", T_su_H=360.0, P_su_H=5.0e5, m_dot_H=1.0)
COLD = dict(fluid_C="Water", T_su_C=300.0, P_su_C=5.0e5, m_dot_C=1.0)


def _heat_exchanger(eta):
    hx = HexCstEff()
    hx.set_parameters(eta=eta)
    inputs = dict(HOT)
    inputs.update(COLD)
    return _solved(hx, **inputs)


@pytest.mark.parametrize("eta", [0.3, 0.6, 0.9])
def test_heat_exchanger_conserves_energy(eta):
    """What the hot stream gives up, the cold stream takes in.

    With no declared losses this must balance to numerical precision, and it
    must equal the duty the model reports on its heat connector. Checking all
    three against each other catches a duty that is computed correctly but
    applied to only one side.
    """
    hx = _heat_exchanger(eta)
    released = hx.su_H.m_dot * (hx.su_H.h - hx.ex_H.h)
    absorbed = hx.ex_C.m_dot * (hx.ex_C.h - hx.su_C.h)
    assert released == pytest.approx(absorbed, rel=1e-9)
    assert hx.Q.Q_dot == pytest.approx(absorbed, rel=1e-9)


@pytest.mark.parametrize("eta", [0.3, 0.6, 0.9])
def test_heat_exchanger_respects_the_second_law(eta):
    """Heat flows hot to cold, and neither stream crosses the other inlet.

    The cold outlet cannot exceed the hot inlet temperature, and the hot outlet
    cannot fall below the cold inlet temperature. Violating either would mean
    the model had invented a temperature difference.
    """
    hx = _heat_exchanger(eta)
    assert hx.Q.Q_dot > 0
    assert hx.ex_C.T > hx.su_C.T
    assert hx.ex_H.T < hx.su_H.T
    assert hx.ex_C.T <= hx.su_H.T + 1e-6
    assert hx.ex_H.T >= hx.su_C.T - 1e-6


def test_heat_exchanger_duty_is_proportional_to_effectiveness():
    """Q = eta * Q_max, so halving eta must halve the duty exactly.

    Taking the ratio cancels Q_max, which means this holds whichever of the two
    Q_max branches the model took -- it is a statement about the model's
    structure, not about the property data.
    """
    duties = {eta: _heat_exchanger(eta).Q.Q_dot for eta in (0.3, 0.6, 0.9)}
    assert duties[0.6] / duties[0.3] == pytest.approx(2.0, rel=1e-6)
    assert duties[0.9] / duties[0.3] == pytest.approx(3.0, rel=1e-6)


def test_heat_exchanger_at_unit_effectiveness_pinches():
    """At eta = 1 one stream must reach the other inlet temperature.

    That is the definition of the thermodynamic maximum duty: the exchanger is
    infinitely large and the pinch closes to zero at one end.
    """
    hx = _heat_exchanger(1.0)
    pinch = min(
        abs(hx.ex_C.T - hx.su_H.T),
        abs(hx.ex_H.T - hx.su_C.T),
    )
    assert pinch < 1e-3, "residual pinch of %.3f K at unit effectiveness" % pinch


def test_heat_exchanger_conserves_mass_on_both_streams():
    hx = _heat_exchanger(0.6)
    assert hx.ex_H.m_dot == pytest.approx(hx.su_H.m_dot, rel=REL_TOL)
    assert hx.ex_C.m_dot == pytest.approx(hx.su_C.m_dot, rel=REL_TOL)
    assert hx.ex_H.fluid == hx.su_H.fluid
    assert hx.ex_C.fluid == hx.su_C.fluid


def test_heat_exchanger_total_entropy_rises():
    """The pair of streams together must produce entropy.

    Each stream on its own moves either way; only the sum is constrained. This
    is the check that a heat exchanger cannot be run backwards.
    """
    hx = _heat_exchanger(0.6)
    generated = (
        hx.ex_H.m_dot * (hx.ex_H.s - hx.su_H.s)
        + hx.ex_C.m_dot * (hx.ex_C.s - hx.su_C.s)
    )
    assert generated > 0
    assert math.isfinite(generated)
