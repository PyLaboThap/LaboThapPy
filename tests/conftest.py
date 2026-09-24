# -*- coding: utf-8 -*-
"""
Shared configuration and helpers for the LaboThApPy test suite.

This module is imported by pytest before any test module, which makes it the
right place to do three things that every test depends on:

1. Force a non-interactive Matplotlib backend and neutralise ``show()``.
   Many LaboThApPy modules build figures at import time; without this the
   suite would block forever on a headless CI runner.

2. Import :mod:`labothappy` once, up front. The package ``__init__`` prepends
   the package directory to ``sys.path``, which is what makes the library's
   current ``from component.x import Y`` style imports resolve. Tests never
   rely on that bootstrap themselves -- they always name modules by their
   canonical ``labothappy.*`` dotted path -- so the suite keeps working
   unchanged once the bootstrap is removed.

3. Expose the discovery helpers (``iter_modules``, ``component_classes``) that
   the import and contract tests are parametrised over.
"""

from __future__ import annotations

import importlib
import importlib.util
import os
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# 1. Headless plotting. Must happen before any library module is imported.
# ---------------------------------------------------------------------------

import matplotlib  # noqa: E402

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt  # noqa: E402


def _noop(*_args, **_kwargs):
    return None


plt.show = _noop
plt.pause = _noop
matplotlib.figure.Figure.show = _noop

# ---------------------------------------------------------------------------
# 2. The library itself.
# ---------------------------------------------------------------------------

import labothappy  # noqa: E402

PACKAGE_ROOT = Path(labothappy.__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent

#: Directories that are never walked when discovering modules.
SKIP_DIRS = {"__pycache__", ".git", ".ipynb_checkpoints"}


# ---------------------------------------------------------------------------
# 3. Module discovery.
# ---------------------------------------------------------------------------


def iter_modules():
    """Yield ``(dotted_name, path)`` for every ``.py`` file under the package.

    ``dotted_name`` is always the canonical ``labothappy.*`` name, even when a
    path segment is not a legal Python identifier (the ``0D`` / ``1D`` sizing
    directories). Those names cannot be fed to :func:`importlib.import_module`,
    which is what :func:`import_module` below works around and what
    ``test_imports.test_package_names_are_legal_identifiers`` tracks.
    """
    for dirpath, dirnames, filenames in os.walk(PACKAGE_ROOT):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        for filename in sorted(filenames):
            if not filename.endswith(".py") or filename == "__init__.py":
                continue
            path = Path(dirpath) / filename
            rel = path.relative_to(REPO_ROOT)
            dotted = ".".join(rel.with_suffix("").parts)
            yield dotted, path


def is_dotted_importable(dotted):
    """True when every segment of *dotted* is a legal Python identifier."""
    return all(part.isidentifier() for part in dotted.split("."))


def import_module(dotted, path):
    """Import *dotted*, falling back to a by-path load for illegal names.

    Raises whatever the module raises. Callers catch ``BaseException`` because
    a handful of modules call bare ``exit()``, which raises ``SystemExit``.
    """
    if is_dotted_importable(dotted):
        return importlib.import_module(dotted)

    spec = importlib.util.spec_from_file_location(dotted, path)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise ImportError("no loader for %s" % path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[dotted] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(dotted, None)
        raise
    return module


# ---------------------------------------------------------------------------
# 4. Component discovery.
# ---------------------------------------------------------------------------

#: Components whose ``__init__`` takes required arguments. Each entry gives the
#: minimal kwargs needed to build a bare, unconfigured instance.
CONSTRUCTOR_ARGS = {
    "HexMBChargeSensitive": {"HTX_Type": "Plate"},
    "TankMixer": {"n_inlets": 2},
    "TankSpliter": {"outlet_repartition": [0.5, 0.5]},
}


def is_component_class(obj):
    """True for a concrete ``BaseComponent`` subclass.

    The check is by *name* along the MRO rather than ``issubclass``, because the
    library currently reaches ``BaseComponent`` through two different import
    paths (``component.base_component`` and
    ``labothappy.component.base_component``). Those produce two distinct class
    objects at runtime, so ``issubclass`` against either one sees only part of
    the family. ``test_component_contract`` tracks that duplication separately;
    discovery must not depend on it being fixed first.
    """
    return (
        isinstance(obj, type)
        and obj.__name__ != "BaseComponent"
        and any(base.__name__ == "BaseComponent" for base in obj.__mro__)
    )


def component_classes():
    """Return ``{"module:ClassName": class}`` for every importable component.

    Example scripts are excluded: they run full simulations at import and define
    no new components. Modules that fail to import are skipped here and covered
    by ``test_imports`` instead, so one broken module cannot silently erase the
    contract coverage of the others.
    """
    found = {}
    component_dir = PACKAGE_ROOT / "component"
    for dotted, path in iter_modules():
        if component_dir not in path.parents:
            continue
        if any(part.lower() == "examples" for part in path.parts):
            continue
        leaf = dotted.rsplit(".", 1)[-1]
        try:
            module = import_module(dotted, path)
        except BaseException:
            continue
        for name, obj in vars(module).items():
            # Only classes *defined* here, not ones imported from a sibling,
            # so each component is exercised exactly once.
            if is_component_class(obj) and obj.__module__.rsplit(".", 1)[-1] == leaf:
                found.setdefault("%s:%s" % (dotted, name), obj)
    return found


def build(cls):
    """Instantiate *cls* with the minimal constructor arguments it needs."""
    instance = cls(**CONSTRUCTOR_ARGS.get(cls.__name__, {}))
    instance.print_flag = 0
    return instance


# ---------------------------------------------------------------------------
# 5. Fixtures.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def props():
    """A cached CoolProp ``AbstractState`` factory.

    ``props(fluid)`` returns a ready ``AbstractState``. Reusing one per fluid
    keeps the analytic tests fast without leaking state between them, since
    every caller calls ``update`` before reading anything back.
    """
    import CoolProp.CoolProp as CP

    cache = {}

    def _get(fluid):
        if fluid not in cache:
            cache[fluid] = CP.AbstractState("HEOS", fluid)
        return cache[fluid]

    return _get


@pytest.fixture
def h_s_of_TP(props):
    """``h_s_of_TP(fluid, T, p)`` -> ``(h, s)``, SI units throughout."""
    import CoolProp.CoolProp as CP

    def _get(fluid, T, p):
        state = props(fluid)
        state.update(CP.PT_INPUTS, p, T)
        return state.hmass(), state.smass()

    return _get


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: takes more than a few seconds (example scripts)"
    )
