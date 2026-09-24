# LaboThapPy: A Thermodynamic Modelling Library

**LaboThapPy** is an open-source Python library that consolidates a wide range of thermodynamic models developed over years of research at the University of Liège. Its primary purpose is to provide a standardized framework for detailed thermodynamic system and component modeling.

## Why LaboThapPy?

At the University of Liège's Thermodynamic Laboratory, we noticed a recurring issue:  
- Identical or similar models were being recreated across different software and codes.  
- Some models were lost when researchers left.  

**LaboThapPy** was created to address these challenges by:  
- Establishing a shared repository of thermodynamic models.  
- Providing a standardized framework for model implementation.  

## Authors

This library was initially developed by:  
- **Elise Neven** (University of Liège) – [elise.neven@uliege.be](mailto:elise.neven@uliege.be)  
- **Basile Chaudoir** (University of Liège) – [Basile.Chaudoir@uliege.be](mailto:Basile.Chaudoir@uliege.be)  

## Documentation

Detailed documentation is available on [ReadTheDocs](https://labothappy.readthedocs.io).

## Reporting Issues

If you encounter any issues or have suggestions, please report them through the [GitHub Issues page](../../issues).

## Installation

```bash
git clone https://github.com/PyLaboThap/PyLaboThap.git
cd PyLaboThap
pip install -e .
```

Python 3.10 or later. `pip install -e ".[optimization]"` adds the
particle-swarm sizing routines; `pip install -e ".[dev]"` adds everything,
including the test suite and the documentation toolchain.

## Running the tests

```bash
pip install -e ".[test]"
pytest                      # everything, about three minutes
pytest -m "not slow"        # the same minus four long-running modules, ~1 min
```

The suite has three parts:

| File | What it checks |
| --- | --- |
| `tests/test_imports.py` | Every module in the package imports, or is listed in `tests/known_import_failures.txt` |
| `tests/test_component_contract.py` | Every component implements the `BaseComponent` protocol |
| `tests/test_physics_analytic.py` | Energy and mass balances, the second law, and closed-form identities for the analytic models |

`tests/known_import_failures.txt` is a **ratchet**. A module that fails to
import and is not listed is a new breakage; a module that is listed but now
imports cleanly is a stale entry and the test tells you to delete the line. The
list can therefore only get shorter.

Contract violations are recorded the same way, as `xfail` entries in
`tests/test_component_contract.py`. Each one names a real bug. Fixing the bug
turns the test from `xfail` into `XPASS`, which is reported — so nothing needs
to be remembered, and nothing is quietly suppressed.
