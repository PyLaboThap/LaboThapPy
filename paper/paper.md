---
title: 'LaboThApPy: a multi-fidelity Python library of thermodynamic component and cycle models'
tags:
  - Python
  - thermodynamics
  - energy systems
  - organic Rankine cycle
  - heat pump
  - heat exchanger
authors:
  # TODO before submission: confirm author order, and add an `orcid:` line to
  # every author who has one -- JOSS asks for them and they cannot be guessed.
  # The field is left out rather than filled with a placeholder, because an
  # ORCID that fails its checksum breaks the Open Journals PDF build.
  # The full contributor list, with affiliations as the project records them,
  # is in AUTHORS.txt at the repository root.
  - name: Elise Neven
    affiliation: 1
  - name: Basile Chaudoir
    affiliation: 1
  - name: Mattéo Hauglustaine
    corresponding: true
    affiliation: "2, 3"
  - name: Alanis Zeoli
    affiliation: 1
  - name: Titouan Janod
    affiliation: 1
  - name: Samuel Gendebien
    affiliation: 1
  - name: Marie Peeters
    affiliation: 1
  - name: Andres Hernandes
    affiliation: 1
affiliations:
  - name: Thermodynamics Laboratory, University of Liège, Belgium
    index: 1
  - name: Institute of Mechanics, Materials and Civil Engineering, UCLouvain, Belgium
    index: 2
  - name: Thermal Engineering and Combustion Unit, University of Mons, Belgium
    index: 3
date: 24 September 2026
bibliography: paper.bib
---

# Summary

`LaboThApPy` is an open-source Python library of steady-state thermodynamic
models for the components of energy conversion systems — compressors,
expanders, pumps, heat exchangers, valves, separators, tanks and solar
collectors — together with the machinery to assemble them into complete cycles
such as organic Rankine cycles (ORC) and heat pumps. Fluid properties come from
`CoolProp` [@bell2014coolprop], the numerical work rests on `NumPy`
[@harris2020numpy] and `SciPy` [@virtanen2020scipy], and every component can
draw its own temperature-entropy diagram through `Matplotlib`
[@hunter2007matplotlib].

The organising idea is *multi-fidelity modelling behind one interface*. Every
component, whatever its internal complexity, is a Python object implementing the
same short protocol: the user sets inlet states through typed connectors, sets
model parameters, and calls `solve()`; outlet states appear on the outlet
connectors. A compressor described by a single constant isentropic efficiency
and a compressor described by the semi-empirical model of @lemort2009scroll —
with its leakage area, supply and exhaust heat transfer conductances and
electromechanical losses — are interchangeable at the point of use. The same
holds for heat exchangers, which range from a constant effectiveness or a
constant pinch, through an $\varepsilon$-NTU formulation, to a charge-sensitive
moving-boundary model [@bell2015movingboundary] specialised for plate,
shell-and-tube, tube-and-fin and printed-circuit geometries, and a
finite-volume cross-flow finned-tube model.

That interchangeability is the point. A researcher can lay out a cycle with the
cheapest model of every component, confirm that the architecture makes sense,
and then replace only the components that matter for the question at hand with
detailed ones, without rewriting the cycle. The library also ships sizing
routines — mean-line design for axial and radial turbomachinery, and
optimisation-driven heat exchanger sizing — so that a component can be designed
and then simulated off-design inside the same framework.

# Statement of need

Thermodynamic component models are written again and again. In a research
group, a semi-empirical expander model is implemented for one thesis in EES,
re-implemented for the next in MATLAB, and re-implemented once more in Python;
each version is validated against the same experimental campaign, and each is
lost when its author leaves. `LaboThApPy` began as an answer to exactly this
problem at the Thermodynamics Laboratory of the University of Liège: a single
place where a model, once validated, stays available, documented and runnable.

Existing tools do not close this gap. Commercial process simulators are closed,
expensive, and difficult to extend with a new correlation. Modelica libraries
such as `ThermoCycle` [@quoilin2014thermocycle] are powerful for dynamic
simulation but require a Modelica toolchain and an equation-based mindset that
is a poor fit for a quick parametric study. Within Python, `TESPy`
[@witte2020tespy] solves thermal engineering *networks* very effectively, but
its strength is the network solver rather than a catalogue of alternative
component physics; `IDAES` [@lee2021idaes] is a full process-systems
optimisation framework whose weight is hard to justify for a single ORC study.

`LaboThApPy` occupies the space between a handful of ad-hoc scripts and a
process-systems framework. It is aimed at researchers and engineers working on
ORC, heat pump, refrigeration and waste-heat-recovery systems
[@quoilin2013orcsurvey], who need models at a specific and often unusual
fidelity, who want to read and modify the physics rather than treat it as a
black box, and who want the component they build for one study to be reusable
in the next. Each model carries a docstring stating its assumptions, its
connectors, its parameters and its inputs, so that the modelling choices behind
a result are legible in the source.

The library is used in ongoing doctoral research on ORC systems, heat pumps and
thermally integrated energy storage across the three partner universities, and
it is the reference implementation for the component models developed in those
projects.

# Acknowledgements

The authors are indebted to Pr. Vincent Lemort, Pr. Francesco Contino and
Pr. Ward De Paepe for their guidance.

<!-- TODO before submission: add funding sources, grant numbers and project
     names, as a sentence here. JOSS expects financial support to be declared.
     Delete this comment once done. -->

# References
