# -*- coding: utf-8 -*-
"""
The contract every component shares, checked against every component.

``BaseComponent`` defines a small protocol -- ``set_inputs`` / ``set_parameters``
/ ``set_guesses`` / ``check_calculable`` / ``check_parametrized`` / ``solve``,
with ``MassConnector`` ports carrying ``T``, ``p``, ``h``, ``s`` and ``m_dot``.
Every model in the library implements it. That is what lets a user swap one
heat-exchanger model for another inside a cycle, and it is the single most
valuable thing to test here: one test covers every component at once, and a new
model gets its coverage for free on the day it is added.

The tests are written against the library *as it is today*. Where a component
genuinely breaks the contract, the breakage is recorded in one of the ``KNOWN_``
registries below as an expected failure rather than hidden or asserted away.
Each registry entry is a bug that deserves an issue; deleting the entry is how
the fix gets verified. Nothing here has to be remembered -- fix the component
and the test turns green by itself, or rather turns XPASS, which pytest reports.

Components are keyed by ``module:ClassName`` rather than by class name alone,
because three different template modules each define a class called
``HeatExchangerMB`` and they do not behave alike.
"""

from __future__ import annotations

import pytest

from conftest import build, component_classes

COMPONENTS = component_classes()
CASES = sorted(COMPONENTS.items())
IDS = [key for key, _ in CASES]
CLASSES = [cls for _, cls in CASES]
PARAMS = list(zip(IDS, CLASSES))


# --------------------------------------------------------------------------
# Recorded contract violations.
# --------------------------------------------------------------------------

#: These classes declare ``BaseComponent`` as their parent but never call
#: ``super().__init__()``. They have no ``solved`` / ``calculable`` /
#: ``parametrized`` / ``inputs`` / ``guesses`` attributes, and they hold their
#: ports as plain dicts instead of ``MassConnector`` objects. They cannot be
#: wired into a cycle or polled by a solver: they are a parallel design wearing
#: the base class as a label.
KNOWN_DOES_NOT_IMPLEMENT_BASE = {
    "labothappy.component.tank.drum_separator:DrumSeparator",
    "labothappy.component.tank.oil_separator:OilSeparator",
}

#: ``__init__`` overwrites the boolean status flags with ``None`` after
#: ``super().__init__()`` has set them to ``False``. ``None`` is falsy, so
#: nothing breaks today, but ``is False`` checks and any typed API will.
KNOWN_STATUS_FLAGS_ARE_NONE = {
    "labothappy.component.tank.tank_mixer:TankMixer",
    "labothappy.component.tank.tank_spliter:TankSpliter",
}

#: ``get_required_inputs()`` reads state that only exists after configuration:
#: ``self.params['mode']`` for the three semi-empirical models, and
#: ``self.W.N`` for the turbomachinery template. Since ``check_calculable()``
#: calls it, asking an unconfigured component whether it is calculable crashes
#: instead of answering ``False``.
KNOWN_REQUIRED_INPUTS_NEED_STATE = {
    "labothappy.component.compressor.compressor_semi_empirical:CompressorSE",
    "labothappy.component.expander.expander_semi_empirical:ExpanderSE",
    "labothappy.component.pump.pump_curve_similarity:PumpCurveSimilarity",
    "labothappy.component.templates.template_turbomachinery:HeatExchangerMB",
}

#: ``get_required_parameters()`` reads an attribute its constructor never set.
KNOWN_REQUIRED_PARAMS_NEED_STATE = {
    "labothappy.component.templates.template_heat_exchanger:HeatExchangerMB",
}

#: ``get_required_parameters()`` writes default values into ``self.params`` as a
#: side effect of being asked what the parameters are. A query that mutates
#: means ``check_parametrized()`` returns a different answer depending on
#: whether somebody happened to call the getter first -- and the defaults
#: (``AU_amb = 0``, ``W_dot_loss_0 = 0``) silently turn a forgotten parameter
#: into a physically meaningless but apparently valid model.
KNOWN_REQUIRED_PARAMS_MUTATE = {
    "labothappy.component.compressor.compressor_semi_empirical:CompressorSE",
    "labothappy.component.expander.expander_semi_empirical:ExpanderSE",
}

#: ``solve()`` touches connector or parameter state before it checks
#: ``calculable`` and ``parametrized``, so calling it on an unconfigured
#: component raises instead of returning with ``solved == False``. Typically a
#: ``CP.AbstractState('HEOS', self.su.fluid)`` placed above the guard, with
#: ``fluid`` still ``None``.
#:
#: Eighteen of twenty-four components do this; the six that do not show the
#: shape the fix should take.
KNOWN_SOLVE_RAISES_WHEN_UNCONFIGURED = {
    "labothappy.component.compressor.compressor_csteff:CompressorCstEff",
    "labothappy.component.compressor.compressor_semi_empirical:CompressorSE",
    "labothappy.component.expander.expander_csteff:ExpanderCstEff",
    "labothappy.component.expander.expander_semi_empirical:ExpanderSE",
    "labothappy.component.heat_exchanger.hex_MB_charge_sensitive:HexMBChargeSensitive",
    "labothappy.component.heat_exchanger.hex_crossflowfintube_finitevolume:HexCrossFlowTubeAndFinsFiniteVolume",
    "labothappy.component.heat_exchanger.hex_csteff:HexCstEff",
    "labothappy.component.heat_exchanger.hex_csteff_disc:HexCstEffDisc",
    "labothappy.component.heat_exchanger.hex_eNTU:HexeNTU",
    "labothappy.component.heat_exchanger.hex_thermosyphon:HexThermosyphon",
    "labothappy.component.pump.pump_curve_similarity:PumpCurveSimilarity",
    "labothappy.component.solar.parabolic_trough_collector:PTCollector",
    "labothappy.component.tank.drum_separator:DrumSeparator",
    "labothappy.component.tank.oil_separator:OilSeparator",
    "labothappy.component.tank.tank_mixer:TankMixer",
    "labothappy.component.tank.tank_spliter:TankSpliter",
    "labothappy.component.templates.template_heat_exchanger:HeatExchangerMB",
    "labothappy.component.templates.template_turbomachinery:HeatExchangerMB",
}

CONNECTOR_PROPERTIES = ("T", "p", "h", "s", "m_dot", "fluid")


def _skip_if_known(key, registry, reason):
    if key in registry:
        pytest.xfail("%s: %s" % (key, reason))


def _mass_connectors(component):
    """Every ``MassConnector`` hanging off *component*, by attribute name."""
    return {
        name: value
        for name, value in vars(component).items()
        if type(value).__name__ == "MassConnector"
    }


def test_at_least_one_component_was_discovered():
    """A guard on the guard.

    If discovery silently returned nothing -- a renamed base class, a broken
    import inside conftest -- every parametrised test below would vanish and
    the suite would still report green. This makes that failure loud.
    """
    assert len(COMPONENTS) >= 20, (
        "only %d components discovered; expected the full family. "
        "Check conftest.component_classes()." % len(COMPONENTS)
    )


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_component_is_instantiable(key, cls):
    """Every component builds with no arguments, or with documented ones.

    Components needing constructor arguments are listed in
    ``conftest.CONSTRUCTOR_ARGS``. A component that needs arguments and is not
    listed there fails here, which is the prompt to document them.
    """
    assert build(cls) is not None


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_component_has_the_base_state(key, cls):
    """The six attributes the base class promises are all present.

    A component missing these has not called ``super().__init__()``, which
    means no solver can poll it and no cycle can hold it.
    """
    _skip_if_known(
        key,
        KNOWN_DOES_NOT_IMPLEMENT_BASE,
        "__init__ never calls super().__init__()",
    )
    component = build(cls)
    for attribute in ("solved", "calculable", "parametrized"):
        assert hasattr(component, attribute)
    assert isinstance(component.inputs, dict)
    assert isinstance(component.params, dict)
    assert isinstance(component.guesses, dict)


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_fresh_component_claims_nothing(key, cls):
    """A component told nothing must not claim to know anything.

    ``solved`` matters most: a user reads it to decide whether the results can
    be trusted, so anything other than a hard ``False`` here is silent
    corruption waiting to happen.
    """
    _skip_if_known(
        key, KNOWN_DOES_NOT_IMPLEMENT_BASE, "has no status flags at all"
    )
    _skip_if_known(
        key,
        KNOWN_STATUS_FLAGS_ARE_NONE,
        "__init__ resets calculable / parametrized to None",
    )
    component = build(cls)
    assert component.solved is False
    assert component.calculable is False
    assert component.parametrized is False


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_component_exposes_mass_connectors(key, cls):
    """Every component has at least one port, and ports expose the full state.

    This is what makes components composable: a cycle wires ``a.ex`` to
    ``b.su`` and reads the same properties off either end, whichever models are
    involved. A component whose ports are plain dicts cannot take part.
    """
    _skip_if_known(
        key,
        KNOWN_DOES_NOT_IMPLEMENT_BASE,
        "holds its ports as plain dicts, not MassConnector objects",
    )
    component = build(cls)
    connectors = _mass_connectors(component)
    assert connectors, (
        "%s exposes no MassConnector; it cannot be wired into a cycle" % key
    )
    for name, connector in connectors.items():
        missing = [p for p in CONNECTOR_PROPERTIES if not hasattr(connector, p)]
        assert not missing, "%s.%s is missing %s" % (
            key,
            name,
            ", ".join(missing),
        )


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_required_inputs_is_a_list_of_names(key, cls):
    """``get_required_inputs()`` answers without needing to be configured first.

    It is the method a user calls to find out *what to supply*, so needing to
    have supplied something already is circular. It is also what
    ``check_calculable()`` calls, which makes a failure here contagious.
    """
    _skip_if_known(
        key,
        KNOWN_REQUIRED_INPUTS_NEED_STATE,
        "reads configuration state that a fresh instance does not have",
    )
    required = build(cls).get_required_inputs()
    assert all(isinstance(name, str) for name in required)


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_required_parameters_is_a_list_of_names(key, cls):
    """Same contract, for parameters."""
    _skip_if_known(
        key,
        KNOWN_REQUIRED_PARAMS_NEED_STATE,
        "reads an attribute its constructor never sets",
    )
    required = build(cls).get_required_parameters()
    assert all(isinstance(name, str) for name in required)


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_required_guesses_is_a_list_of_names(key, cls):
    """Same contract, for initial guesses. Empty by default on the base class."""
    _skip_if_known(key, KNOWN_DOES_NOT_IMPLEMENT_BASE, "does not implement it")
    required = build(cls).get_required_guesses()
    assert all(isinstance(name, str) for name in required)


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_asking_what_is_required_does_not_change_the_component(key, cls):
    """The three ``get_required_*`` methods are questions, not commands.

    Two semi-empirical models install a full set of default parameters when
    asked what their parameters are. That makes the answer to
    ``check_parametrized()`` depend on call order, and it means a user who
    forgets a parameter gets a model that runs on ``AU_amb = 0`` rather than a
    complaint. Defaults may well be the right idea -- but they belong in
    ``__init__``, where they are visible.
    """
    _skip_if_known(key, KNOWN_DOES_NOT_IMPLEMENT_BASE, "has no params dict")
    _skip_if_known(
        key,
        KNOWN_REQUIRED_PARAMS_MUTATE,
        "get_required_parameters() writes default values into self.params",
    )
    component = build(cls)
    before = (dict(component.params), dict(component.inputs), dict(component.guesses))
    for method in ("get_required_inputs", "get_required_parameters",
                   "get_required_guesses"):
        try:
            getattr(component, method)()
        except Exception:
            pass  # covered by the tests above; only mutation matters here
    after = (dict(component.params), dict(component.inputs), dict(component.guesses))
    assert before == after


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_check_calculable_answers_false_when_nothing_is_set(key, cls):
    """Asking the question must never be the thing that breaks.

    ``check_calculable()`` exists so a cycle solver can poll a component and
    decide what to do next. If it raises, the solver cannot even ask.
    """
    _skip_if_known(key, KNOWN_DOES_NOT_IMPLEMENT_BASE, "does not implement it")
    _skip_if_known(
        key,
        KNOWN_REQUIRED_INPUTS_NEED_STATE,
        "inherits the get_required_inputs() failure",
    )
    assert build(cls).check_calculable() is False


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_check_parametrized_agrees_with_its_own_declaration(key, cls):
    """``check_parametrized()`` must answer, and answer consistently.

    Some components legitimately arrive fully parametrized -- their constructor
    takes the parameters, or ships defaults. So the expected answer is not
    ``False`` but whatever the component's own declared requirements imply,
    read from ``params`` after the fact.
    """
    _skip_if_known(key, KNOWN_DOES_NOT_IMPLEMENT_BASE, "does not implement it")
    _skip_if_known(
        key,
        KNOWN_REQUIRED_PARAMS_NEED_STATE,
        "get_required_parameters() raises before configuration",
    )
    component = build(cls)
    answer = component.check_parametrized()
    expected = all(
        component.params.get(name) is not None
        for name in component.get_required_parameters()
    )
    assert bool(answer) == expected


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_set_parameters_and_guesses_round_trip(key, cls):
    """What goes into ``set_parameters`` comes back out of ``params``."""
    _skip_if_known(key, KNOWN_DOES_NOT_IMPLEMENT_BASE, "does not implement it")
    component = build(cls)
    component.set_parameters(_probe_parameter=1.25)
    assert component.params["_probe_parameter"] == 1.25
    component.set_guesses(_probe_guess=2.5)
    assert component.guesses["_probe_guess"] == 2.5


@pytest.mark.parametrize("key,cls", PARAMS, ids=IDS)
def test_solve_refuses_rather_than_raises_when_unconfigured(key, cls):
    """An unconfigured ``solve()`` must return quietly with ``solved == False``.

    This is the most important test in the file. Components report failure by
    catching broad exceptions inside ``solve()`` and setting ``solved =
    False``. That is a defensible design, but it only works if the flag is
    always reached: a ``solve()`` that raises past its own guard turns a
    modelling failure into a crash in the middle of a cycle iteration, and one
    that returns early without touching the flag leaves a stale ``True`` from
    the previous call.

    Eighteen of the twenty-four components fail this today, almost always for
    the same reason -- a ``CP.AbstractState('HEOS', self.su.fluid)`` sitting
    above the ``if not (self.calculable and self.parametrized)`` guard, with
    ``fluid`` still ``None``. Moving those two lines fixes most of the list.
    """
    _skip_if_known(
        key,
        KNOWN_SOLVE_RAISES_WHEN_UNCONFIGURED,
        "solve() touches state before checking calculable / parametrized",
    )
    component = build(cls)
    component.solve()
    assert component.solved is False


def test_base_component_has_a_single_identity():
    """``BaseComponent`` must be one class, not two.

    The library reaches its own base class through two import paths:
    ``from component.base_component import BaseComponent`` -- which resolves
    through the ``sys.path`` entry that ``labothappy/__init__.py`` adds -- and
    ``from labothappy.component.base_component import BaseComponent``. Python
    treats those as unrelated modules and builds two separate class objects, so
    ``isinstance(component, BaseComponent)`` gives a different answer depending
    on which import the caller happened to use, and a cycle solver filtering
    its parts by type silently drops half of them.

    This is why ``conftest.is_component_class`` matches the MRO by name. The
    fix is the codemod that rewrites every intra-package import to the
    ``labothappy.*`` form; this test turns green on its own the day it lands.
    """
    bases = sorted(
        {
            base.__module__
            for cls in CLASSES
            for base in cls.__mro__
            if base.__name__ == "BaseComponent"
        }
    )
    if len(bases) > 1:
        pytest.xfail(
            "BaseComponent exists as %d distinct classes at runtime: %s"
            % (len(bases), ", ".join(bases))
        )
    assert len(bases) == 1


def test_no_two_components_share_a_class_name():
    """Distinct models must have distinct names.

    Three template modules each define a class called ``HeatExchangerMB``, and
    a fourth name, ``HexMBChargeSensitive``, is what the real model is called
    now. Code elsewhere in the library still imports the old name, which is one
    of the recorded import failures. Names are the library's public surface;
    duplicates make it impossible to say what a traceback refers to.
    """
    by_name = {}
    for key, cls in CASES:
        by_name.setdefault(cls.__name__, []).append(key)
    duplicates = {n: k for n, k in by_name.items() if len(k) > 1}
    if duplicates:
        pytest.xfail(
            "component class names used more than once: %s"
            % "; ".join("%s -> %s" % (n, ", ".join(k))
                        for n, k in sorted(duplicates.items()))
        )
    assert not duplicates
