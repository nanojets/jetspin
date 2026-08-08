# Changelog

This file summarizes the main user-visible changes between JETSPIN releases.

## Version 1.22 — May 2017

Version 1.22 extends the physical models available in JETSPIN and improves
simulation control, output, and numerical robustness.

### New features

- Added time-dependent external electric fields, including rectangular,
  RC-circuit, and orthogonal rotating waveforms.
- Added solvent evaporation modelling, including the evolution of polymer
  fraction, viscosity, jet radius, and related thermodynamic parameters.
- Added stochastic forcing through a configurable random-noise model.
- Added controlled job-time and close-time directives, allowing simulations
  to stop cleanly before an external execution-time limit is reached.
- Extended deterministic integrators to account for aerodynamic drag.
- Added support for recording evaporation-related quantities in simulation
  output and restart data.

### Examples and documentation

- Added example 7, demonstrating an evaporating jet under an orthogonal
  rotating electric field.
- Updated the user manual to version 1.22.
- Added the LaTeX manual sources, figures, bibliography, reproducible local
  build command, and automatic PDF build through GitHub Actions.
- Added a numerical smoke-test suite covering all seven bundled examples.

### Fixes and maintenance

- Fixed Coulomb-force handling used by evaporation simulations.
- Fixed removed-bead data output and related evaporation bookkeeping.
- Improved restart handling.
- Added GNU Fortran compatibility flags for legacy MPI argument conventions
  on both older and current compiler releases.
- Fixed out-of-bounds argument passing in the one-dimensional RK2 and RK4
  integrators, detected by the new runtime-checking smoke tests.
- Retained the gfortran compatibility corrections for the parser and I/O
  modules.

## Version 1.21 — July 2016

Version 1.21 introduced:

- dynamic mesh refinement;
- a multiple-timestep algorithm for long-range Coulomb forces;
- the Kelvin–Voigt rheological model;
- nanoparticle and mass-impurity modelling;
- Lorentz-force support;
- PDB/PSF topology output and additional trajectory tools;
- examples 5 and 6 for dynamic refinement and quasi-Newtonian fluids.
