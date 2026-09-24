# -*- coding: utf-8 -*-
"""
Every module in the package must import, or be on the shrinking allowlist.

This is the cheapest test in the suite and the one that catches the most. An
import that raises means a module is unreachable to a user, to Sphinx autodoc,
and to any other module that depends on it -- so a broken import silently
removes a model from the library without anybody noticing.

The allowlist in ``known_import_failures.txt`` is a ratchet. A module that
fails and is not listed is a regression; a module that is listed but now
imports is a stale entry that must be deleted. Neither direction can be
ignored, so the list can only shrink.

Note on isolation: modules are imported into the running interpreter rather
than a subprocess. That is deliberate -- it is how a user experiences the
library, and it is what surfaces cross-module state leakage. It also means a
module already pulled in as a dependency of an earlier one is cached, which is
why the whole sweep costs far less than the sum of its parts.
"""

from __future__ import annotations

import traceback
from pathlib import Path

import pytest

from conftest import import_module, is_dotted_importable, iter_modules

ALLOWLIST_PATH = Path(__file__).parent / "known_import_failures.txt"


def _read_allowlist():
    """Return ``{module_name: comment}`` from the allowlist file."""
    entries = {}
    for raw in ALLOWLIST_PATH.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        module, _, comment = line.partition("#")
        entries[module.strip()] = comment.strip()
    return entries


KNOWN_FAILURES = _read_allowlist()
ALL_MODULES = sorted(iter_modules())
MODULE_NAMES = [dotted for dotted, _ in ALL_MODULES]

#: Modules that run a long simulation when imported. Together these four make
#: up roughly nine tenths of the sweep's three-minute runtime, so marking them
#: lets `pytest -m "not slow"` give a fast local loop while CI still runs
#: everything. Fixing the underlying problem -- example scripts that execute on
#: import rather than under `if __name__ == "__main__"` -- removes the need for
#: this list entirely.
SLOW_MODULES = {
    "labothappy.component.Examples.heat_exchanger.hex_thermosyphon_example",
    "labothappy.component.Examples.tank.tank_LV_separator_example",
    "labothappy.machine.examples.solver_comparison.HP_cycle_iterative",
    "labothappy.sizing.heat_exchanger.heat_pipe_HTX.sizing_HP_HTX",
}


def _params(modules):
    """Turn ``(dotted, path)`` pairs into pytest params, marking the slow ones."""
    return [
        pytest.param(
            dotted,
            path,
            id=dotted,
            marks=(pytest.mark.slow,) if dotted in SLOW_MODULES else (),
        )
        for dotted, path in modules
    ]


EXPECTED_OK = [m for m in ALL_MODULES if m[0] not in KNOWN_FAILURES]
EXPECTED_BAD = [m for m in ALL_MODULES if m[0] in KNOWN_FAILURES]


def _describe(exc):
    """A one-line-plus-traceback description, readable in CI output."""
    return "%s: %s\n\n%s" % (
        type(exc).__name__,
        exc,
        "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
    )


@pytest.mark.parametrize("dotted,path", _params(EXPECTED_OK))
def test_module_imports(dotted, path):
    """Every module not on the allowlist must import cleanly."""
    try:
        import_module(dotted, path)
    except BaseException as exc:  # noqa: BLE001 - some modules call exit()
        pytest.fail(
            "%s failed to import.\n\n"
            "If this is a new breakage, fix the module. If it is a known and "
            "accepted failure, add it to tests/known_import_failures.txt with "
            "the cause as a comment -- and open an issue for it.\n\n%s"
            % (dotted, _describe(exc))
        )


@pytest.mark.parametrize("dotted,path", _params(EXPECTED_BAD))
def test_known_failure_has_not_been_fixed_silently(dotted, path):
    """A listed module that now imports is a stale entry: delete the line.

    Without this half, the allowlist would slowly fill with entries nobody
    revisits, and the suite would report a green wall while coverage quietly
    rotted. Failing here is good news -- it means somebody fixed a module.
    """
    try:
        import_module(dotted, path)
    except BaseException:
        return  # still broken, as recorded
    pytest.fail(
        "%s now imports cleanly.\n\n"
        "Remove its line from tests/known_import_failures.txt so the module "
        "stays covered from now on. Recorded cause was: %s"
        % (dotted, KNOWN_FAILURES[dotted] or "(none given)")
    )


def test_allowlist_has_no_entries_for_missing_modules():
    """Guard against the allowlist drifting away from the tree.

    A line naming a module that no longer exists (renamed, moved, deleted)
    would otherwise sit there forever, silently excusing nothing.
    """
    stale = sorted(set(KNOWN_FAILURES) - set(MODULE_NAMES))
    assert not stale, (
        "tests/known_import_failures.txt names modules that do not exist:\n  "
        + "\n  ".join(stale)
    )


def test_allowlist_is_sorted_within_groups():
    """Keep the file mergeable.

    Several people will be deleting lines from this file at once. Sorted
    groups keep the diffs small and the conflicts rare.
    """
    groups = []
    current = []
    for raw in ALLOWLIST_PATH.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("#") or not line:
            if current:
                groups.append(current)
                current = []
            continue
        current.append(line.partition("#")[0].strip())
    if current:
        groups.append(current)
    for group in groups:
        assert group == sorted(group), (
            "allowlist group is out of order; sort it:\n  "
            + "\n  ".join(group)
        )


def test_package_names_are_legal_identifiers():
    """No directory or module may be unreachable by a dotted import.

    ``labothappy/sizing/turbomachinery/turbine/0D/`` and its siblings cannot
    appear in an ``import`` statement at all, because ``0D`` is not a Python
    identifier. The modules underneath are therefore dead to every consumer of
    the library, including Sphinx autodoc. Renaming the directories (``0D`` ->
    ``zero_d``, ``1D`` -> ``one_d``, or folding them into the model name) is
    the fix.

    This is marked xfail rather than skipped so it turns green by itself the
    day the rename lands.
    """
    illegal = sorted(d for d in MODULE_NAMES if not is_dotted_importable(d))
    if illegal:
        pytest.xfail(
            "%d modules live under a directory whose name is not a Python "
            "identifier and can never be imported:\n  %s"
            % (len(illegal), "\n  ".join(illegal))
        )


def test_no_module_shadows_a_stdlib_name():
    """A top-level module that shadows the standard library breaks users.

    The package directory is prepended to ``sys.path`` by
    ``labothappy/__init__.py``, so any module directly inside it competes with
    the stdlib for that name in every process that imports the library.
    """
    import sys

    package_root = Path(__import__("labothappy").__file__).resolve().parent
    stdlib = set(getattr(sys, "stdlib_module_names", ()))
    top_level = {
        p.stem
        for p in package_root.iterdir()
        if p.suffix == ".py" and p.stem != "__init__"
    } | {p.name for p in package_root.iterdir() if p.is_dir() and (p / "__init__.py").exists()}
    clashes = sorted(top_level & stdlib)
    assert not clashes, (
        "these names are put on sys.path by labothappy/__init__.py and shadow "
        "the standard library: " + ", ".join(clashes)
    )


def test_every_module_name_is_unique_case_insensitively():
    """Two modules differing only in case cannot coexist on Windows or macOS.

    Catching this in CI on Linux is the only way to stop it reaching a
    collaborator who cannot then check the repository out at all.
    """
    seen = {}
    for dotted in MODULE_NAMES:
        seen.setdefault(dotted.lower(), []).append(dotted)
    collisions = {k: v for k, v in seen.items() if len(v) > 1}
    assert not collisions, (
        "module names that collide on a case-insensitive filesystem:\n  "
        + "\n  ".join("%s -> %s" % (k, v) for k, v in sorted(collisions.items()))
    )
