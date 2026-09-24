# -*- coding: utf-8 -*-
"""
Put the Apache-2.0 header from LICENSE_HEADER.txt at the top of every source file.

Usage
-----

    python tools/add_license_headers.py            # report what would change
    python tools/add_license_headers.py --apply    # actually write the files
    python tools/add_license_headers.py --check    # exit 1 if any file lacks it

The default is a dry run. Nothing is written unless ``--apply`` is given.

Why this is a script and not a one-off edit
-------------------------------------------

The header goes on roughly two hundred files, which is a diff that conflicts
with every branch open at the time. Running it is therefore a scheduling
decision, not a coding one: merge the outstanding branches first, run this on a
quiet tree, and merge the result the same day. Keeping it as a script means the
decision can be taken later, and that the same header can be applied again to
files added afterwards.

What it is careful about
------------------------

* **Placement.** A shebang line and a PEP 263 coding declaration must stay in
  the first two lines of the file or Python stops honouring them. The header is
  inserted after those, before the module docstring.
* **Idempotence.** A file that already contains the header is left alone, so
  the script can be re-run safely and used as a CI check.
* **Encoding.** Files are read as UTF-8 and written back as UTF-8 with the line
  ending they already had. A file that is not valid UTF-8 is reported and
  skipped rather than guessed at -- the header contains accented characters and
  writing it into a mis-decoded file would corrupt the source.
"""

from __future__ import annotations

import argparse
import io
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
HEADER_PATH = REPO_ROOT / "LICENSE_HEADER.txt"

#: Directories never touched: third-party or generated content, and the build
#: artifacts that should not be in the repository in the first place.
SKIP_DIRS = {
    "__pycache__",
    ".git",
    ".ipynb_checkpoints",
    "build",
    "_build",
    "dist",
    "docs",
    "LaboThapPy.egg-info",
}

#: Trees that do get the header.
INCLUDE_ROOTS = ("labothappy", "tests", "tools")

#: A file already carrying the header is recognised by this phrase. Matching on
#: the licence name rather than the whole block means a file whose header was
#: hand-edited (a different year, an extra contributor line) is still left
#: alone instead of gaining a second copy.
SENTINEL = "Licensed under the Apache License, Version 2.0"

CODING_RE = re.compile(rb"^[ \t\f]*#.*?coding[:=][ \t]*([-_.a-zA-Z0-9]+)")


def read_header():
    """The header as a list of ``#``-prefixed lines, without a trailing blank."""
    text = HEADER_PATH.read_text(encoding="utf-8").rstrip("\n")
    lines = []
    for line in text.split("\n"):
        lines.append("# " + line if line.strip() else "#")
    return lines


def iter_source_files():
    """Every ``.py`` file under the included roots, in a stable order."""
    for root_name in INCLUDE_ROOTS:
        root = REPO_ROOT / root_name
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.py")):
            if any(part in SKIP_DIRS for part in path.relative_to(REPO_ROOT).parts):
                continue
            # This script quotes the sentinel phrase, so it would always look
            # as though it already had the header. Exclude it explicitly.
            if path.resolve() == Path(__file__).resolve():
                continue
            yield path


def split_preamble(lines):
    """Return ``(preamble, rest)``.

    The preamble is the shebang and the PEP 263 coding declaration, which have
    to stay within the first two lines of the file. Everything else, including
    the module docstring, belongs after the header.
    """
    index = 0
    if index < len(lines) and lines[index].startswith("#!"):
        index += 1
    # PEP 263 allows the coding line on line 1 or line 2.
    if index < len(lines) and index < 2:
        if CODING_RE.match(lines[index].encode("utf-8", "replace")):
            index += 1
    return lines[:index], lines[index:]


def newline_of(raw_text):
    """The line ending the file already uses, so the diff stays minimal."""
    if "\r\n" in raw_text:
        return "\r\n"
    if "\r" in raw_text:
        return "\r"
    return "\n"


def process(path, header_lines, apply_changes):
    """Return one of ``"has-header"``, ``"would-add"``, ``"added"``, ``"skipped"``."""
    try:
        raw = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return "skipped"

    if SENTINEL in raw:
        return "has-header"

    newline = newline_of(raw)
    lines = raw.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    preamble, rest = split_preamble(lines)

    # One blank line between the header and whatever follows, and no run of
    # blank lines left over from the file's own leading whitespace.
    while rest and not rest[0].strip():
        rest.pop(0)

    new_lines = preamble + header_lines + [""] + rest
    new_text = newline.join(new_lines)

    if apply_changes:
        with io.open(path, "w", encoding="utf-8", newline="") as handle:
            handle.write(new_text)
        return "added"
    return "would-add"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument(
        "--apply", action="store_true", help="write the files (default: dry run)"
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit 1 if any file is missing the header; writes nothing",
    )
    args = parser.parse_args(argv)

    if not HEADER_PATH.exists():
        print("missing %s" % HEADER_PATH, file=sys.stderr)
        return 2

    header_lines = read_header()
    counts = {"has-header": 0, "would-add": 0, "added": 0, "skipped": 0}
    missing = []
    skipped = []

    for path in iter_source_files():
        result = process(path, header_lines, args.apply and not args.check)
        counts[result] += 1
        rel = path.relative_to(REPO_ROOT).as_posix()
        if result in ("would-add", "added"):
            missing.append(rel)
        elif result == "skipped":
            skipped.append(rel)

    print("%d files already carry the header" % counts["has-header"])
    if args.apply and not args.check:
        print("%d files updated" % counts["added"])
    else:
        print("%d files would be updated" % counts["would-add"])
    if skipped:
        print("\n%d files skipped -- not valid UTF-8, fix the encoding first:"
              % len(skipped))
        for rel in skipped:
            print("  " + rel)

    if args.check and (missing or skipped):
        print("\nRun: python tools/add_license_headers.py --apply")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
