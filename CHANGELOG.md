# Changelog

This file summarizes the main user-visible changes between JETSPIN releases.

## Version 2.0-beta.1 — October 2026 (pre-release)

Version 2.0 adds a single-GPU OpenACC backend, makes stochastic and
restarted runs reproducible across builds, MPI ranks and restarts, and adds
the validation suites that check these properties. This first pre-release,
2.0-beta.1, is tagged `v2.0-beta.1` on the `development` branch; `master`
stays at 1.22 until the final 2.0 is merged into it (tag `v2.0`). Numerical
results must still be validated for each application. The detailed record
of the work and of its validation is `docs/STATE.md`.

### New features

- Single-GPU OpenACC build `nvfortran-openacc` (NVIDIA HPC SDK; compute
  capability and CUDA toolkit set with `GPUCC` and `CUDA_VERSION`, see
  `docs/introduction/compiling.md`). Three-dimensional Euler, RK2 and RK4
  runs (Maxwell or Kelvin–Voigt model) and Platen runs (Maxwell model), with
  or without evaporation and air drag, take one device step once the nozzle
  bead index reaches 100: the state stays on the GPU, nozzle insertion, collector removal and
  capacity growth run on the device, and accepted refinement events
  interpolate on the device. Below 100 beads, and for the options not yet
  ported (among them multiple-step Coulomb, refinement with the RK schemes,
  time-dependent fields, the Lorentz force, one-dimensional systems and
  MPI), the OpenACC build runs the CPU code. Development builds check the device path
  against the host equations (`nvfortran-openacc-host`, force and Coulomb
  oracles, Akima and refinement-assembly comparisons).
- Sequential Gaussian pool for the stochastic (Platen) integrator, the
  default in every build: a run reads the same noise on CPU and GPU and for
  any number of MPI ranks. Its size is set with `noise pool`.
- Exact restart: `save.dat` is written in double precision with a versioned
  state record (generator state, pool position, refinement and tag
  counters, device-step state); a restarted run continues the uninterrupted
  one exactly.
- `dynamic refinement capacity` sets the bead-capacity reserve and growth
  increment of refining runs; a refinement threshold below twenty times the
  resolution raises an advisory warning (108), since a fine threshold with a
  short refinement interval can thin the cross section over many events.
- Test Cases 9–25: fixed and dynamic topologies, evaporation, Kelvin–Voigt,
  refinement and collector removal, and two runs grown from a single bead
  (Tests 24 and 25), with CPU and GPU references.
- Test suites: numerical regression (`tests/regression`, serial and two MPI
  ranks against versioned baselines), restart (`tests/restart`),
  refinement (`tests/refinement`), dynamic evaporation
  (`tests/performance/dynamic`), besides the smoke tests.
- Markdown documentation in `docs/` alongside the LaTeX manual, whose
  sources are split by chapter.

### Behaviour changes

- A bead that reaches the collector is held there in every run: it
  discharges on the grounded electrode, stops moving and leaves the Coulomb
  sums. `removing yes` only deletes the collected beads. Before, a run
  without removal let its beads cross the collector plane, keeping their
  charge and their Coulomb interactions.
- In a run without removal, dynamic refinement leaves the fiber deposited on
  the collector unchanged and remeshes only the free jet.
- A run makes `final time` divided by `timestep` steps, rounded up, with
  every compiler. Before, compilers that round divisions exactly (GFortran,
  NVFORTRAN 25.5) made one step more when the final time was a multiple of
  the timestep.
- `traj.xyz`, `frame%06d.xyz` and `frame%06d.pdb` give the coordinates in
  cm times the rescale factor (`print xyz rescalexyz`, `print pdb
  rescalepdb`, default 1). Before, `traj.xyz` was written in reduced units
  without `rescalexyz`, and the frame files applied the length unit twice
  with a rescale factor.
- `primary cutoff` is read in cm also without multiple step, where it
  truncates the one-dimensional Coulomb sum; before, that comparison was
  made in reduced units.
- Example 4 (Platen) reads the Gaussian pool; its regression baselines were
  regenerated.
- Test Case 25 was redefined with a 50 % initial polymer fraction.

### Fixes

- Dynamic refinement: the endpoint volume of a remeshed jet, the floor of
  the evaporated volume at an accepted event, and the Akima refit of the
  cross section (positivity, despiking, monotonicity).
- Restart: the evaporation state, the capacity reserve, and the position in
  the Gaussian pool are restored.
- The loop timer and `JETSPIN_PROFILE` use 64-bit clock counts (NVFORTRAN
  builds wrapped after 2147 s); after a restart the throughput counts only
  the steps of the restarted run.
- The input summary reported `dragvel` as on by default; it is off.
- Multiple-step Coulomb sums (since 1.21 unless stated):
  - the per-rank arrays were one element short, so a jet that filled its
    arrays wrote past them (wrong forces or a crash);
  - in one-dimensional runs the second full evaluation truncated the sum at
    the primary cutoff while the list rebuild did not, so the far field was
    extrapolated with the wrong sign and grew over the interval;
  - a new list requested just after an update (insertion, removal, remesh,
    larger arrays) was dropped until the next scheduled update; a removal
    did not request one, nor (with the collector rule of this version) a
    bead frozen at the collector without removal;
  - with evaporation the second full evaluation used the non-evaporative
    masses (since 1.22);
  - once per interval one force evaluation kept the forces of the previous
    call, or took the direct sum when `erms` was printed, so printing `erms`
    changed the trajectory;
  - the maximum-displacement test skipped the first bead of each rank, so
    serial and parallel runs could rebuild the list at different steps.
- `job time` is checked against the elapsed time since the start of the
  run; it used the processor time of rank 0, which can fall behind the
  limit of a batch system.

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
