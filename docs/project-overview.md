# Project overview

This page records the context needed to understand and maintain JETSPIN
without duplicating the detailed scientific manual.

## Project status

- The latest stable release is JETSPIN 1.22, tagged as `v1.22`.
- `master` represents the stable line.
- `development` is the integration branch for ongoing work.
- The next release is JETSPIN 2.0. Its first pre-release, 2.0-beta.1, is
  tagged `v2.0-beta.1` on `development` (the manual, the program banner and
  the [changelog](../CHANGELOG.md) use this name); later pre-releases are
  numbered beta.2, ... or rc.1, ..., and the final 2.0 is merged into
  `master` and tagged `v2.0`.
- The implementation is research software and numerical results must be
  validated for each scientific application.

## Purpose and physical model

JETSPIN simulates electrospinning by discretizing an electrified polymer jet
into charged beads connected by viscoelastic elements. The code supports
one- and three-dimensional systems and includes:

- Maxwell and Herschel–Bulkley constitutive behaviour;
- Kelvin–Voigt rheology;
- external static and time-dependent electric fields;
- Coulomb, viscoelastic, surface-tension, gravitational, and aerodynamic
  forces;
- stochastic forcing;
- bead insertion and removal;
- multiple-timestep Coulomb calculations;
- dynamic mesh refinement;
- Yarin 2001 solvent evaporation and concentration-dependent rheological
  solidification;
- a JETSPIN extension coupling evaporation to Kelvin–Voigt rheology and
  dynamic refinement.

The [scientific PDF manual](../manual/manual.pdf) is authoritative for the
equations, nondimensionalization, assumptions, citations, and limitations.

## Implementation map

The code is written in modular Fortran. Important entry points include:

| Path | Responsibility |
| --- | --- |
| `source/main.f90` | Program lifecycle and main integration loop |
| `source/nanojet_mod.f90` | Jet state, parameters, allocation, and bead data |
| `source/eom_mod.f90` | Standard equations of motion |
| `source/eom_ev_mod.f90` | Equations of motion with evaporation |
| `source/driver_eom_mod.f90` | Selection and dispatch of equation variants |
| `source/integrator_mod.f90` | General Euler, Heun, RK4, and Platen integration |
| `source/integrator_kv_ev_mod.f90` | Kelvin–Voigt/evaporation integration |
| `source/device_step_mod.f90` | OpenACC device step: gate, engagement, stage evaluator, Euler/RK2/RK4 and Platen steps |
| `source/openacc_accelerator_mod.f90` | Module `accelerator_mod` for the OpenACC targets: device kernels and data management |
| `source/host_accelerator_mod.f90` | Module `accelerator_mod` for the CPU targets: stubs (selected by the Makefile through `ACCOBJ`) |
| `source/coulomb_force_mod.f90` | Coulomb sums (direct, mirror charges, multiple step) and their offload |
| `source/viscoelastic_force_mod.f90`, `source/electric_field_mod.f90` | Viscoelastic forces; static and time-dependent external fields |
| `source/dynamic_refinement_mod.f90` | Adaptive jet remeshing |
| `source/fit_mod.f90` | Host and OpenACC Akima coefficient/interpolation kernels |
| `source/io_mod.f90` | Input parsing, output, restart, and reporting |
| `source/statistic_mod.f90` | Runtime and statistical observables |
| `source/utility_mod.f90` | Random-number generator and sequential Gaussian pool |
| `source/profiling_mod.f90` | Optional region profiler (`JETSPIN_PROFILE`) |
| `source/parse_mod.f90`, `source/error_mod.f90` | Input-string parsing; numbered warnings and errors |
| `source/serial_version_mod.f90` | Serial communication abstraction |
| `source/parallel_version_mod.f90` | MPI replicated-data implementation |

For the communication model, bead partitioning, collectives, scalability,
and MPI maintenance invariants, see the dedicated
[MPI parallelization guide](introduction/parallelization.md).
For the moving active interval, capacity growth, compaction, refinement, and
workspace-resizing invariants, see
[dynamic allocation and bead indexing](introduction/dynamic-allocation.md).
The adaptive remeshing algorithm, conservation rules, anchors, and
refinement-specific validation are covered by the
[dynamic-refinement guide](introduction/dynamic-refinement.md).
Stochastic stream ownership, Gaussian broadcasts, and reproducibility limits
are documented in
[random numbers and MPI reproducibility](introduction/random-numbers.md).

## Repository map

| Directory | Contents |
| --- | --- |
| `source/` | Fortran source code |
| `build/` | Program Makefile |
| `execute/` | Default runtime working directory |
| `examples/` | Reference, regression, and GPU-development `input.dat` cases |
| `tests/smoke/` | Serial, runtime-debug, and MPI smoke suite |
| `tests/regression/` | Eight-case numerical baselines and serial/MPI comparisons |
| `tests/evaporation/` | Yarin and Kelvin–Voigt evaporation checks |
| `tests/refinement/` | Anchored Akima-remeshing, repeated growth, and collector-removal checks |
| `tests/restart/` | Exact-restart check: restarted runs against uninterrupted ones (CPU, GPU, MPI) |
| `tests/performance/` | CPU/GPU benchmark records (Tests 9-13, 20) and the dynamic evaporation validation (Tests 16, 17) |
| `docs/` | Concise operational Markdown documentation |
| `manual/` | Scientific LaTeX sources and tracked PDF |
| `tools/` | Trajectory conversion and auxiliary tools |

## Operational invariants

These behaviours are important when changing code or documentation:

- Input is free-format and case-insensitive.
- The runtime input file must be named `input.dat`.
- The final input record must be `finish`.
- Dimensional input and output use the centimetre–gram–second system.
- Programs read and write files in their current working directory.
- Some historical input spellings, including `evaporation umidity`, are
  retained for compatibility and must not be silently renamed.
- The build Makefile is invoked from `source/` through
  `make -C source -f ../build/Makefile <target>`; it is not copied.
- The stable executable name is `main.x`.
- Jet arrays use a moving active interval inside a larger zero-based
  allocation; `inpjet` is not necessarily zero and `mxnpjet` is not the bead
  count.
- A bead that reaches the collector is frozen there in every run (no
  motion, no Coulomb interaction); `removing yes` only deletes collected
  beads. Without removal, dynamic refinement leaves the deposited fiber
  unchanged.
- A run makes `final time` divided by `timestep` steps, rounded up (to
  within a relative 1e-12), with every compiler.
- The OpenACC build reproduces the CPU build of the same compiler byte for
  byte until the device step engages (the nozzle bead index reaches 100).

## Validation expectations

Before merging numerical or build changes, run checks proportional to their
scope. The standard pre-merge set is:

```sh
tests/smoke/run.sh serial
tests/smoke/run.sh debug
tests/smoke/run.sh mpi
tests/regression/run.sh
tests/restart/run.sh gfortran     # default backend: nvfortran
make -C manual
```

On a GPU node add `tests/regression/run.sh openacc` (its cases stay below
the device gate, so it checks the OpenACC build's host path),
`tests/restart/run.sh openacc`,
`tests/performance/dynamic/validate_evaporation.sh` and the refinement
runners of `tests/refinement/`; the GPU kernels are exercised by those and by
Test Cases 9-25.

The serial and debug suites exercise all eight examples plus evaporation,
Kelvin–Voigt, and dynamic-refinement regression paths. MPI mode builds the
parallel implementation and runs a representative case on two processes.
The numerical regression suite compares 1,000-step serial runs of all eight
cases with versioned baselines, then compares two-rank MPI output with the
matching serial run. Its baselines, tolerances, stochastic coverage, and
failure workflow are described in the
[numerical-regression guide](introduction/numerical-regression.md). GitHub
Actions mirrors the smoke checks.

## Documentation policy

JETSPIN deliberately uses two complementary documentation formats:

- `docs/` is the primary source for concise operational guidance that should
  render well on GitHub.
- `manual/` is the primary source for scientific theory, full equations,
  bibliography, numbered references, and archival PDF output.

Avoid copying long scientific passages or complete directive tables into
both formats. Markdown pages should summarize, link to executable examples,
and refer readers to the relevant LaTeX source or PDF section. When a command,
file name, branch convention, or runtime behaviour changes, update the
Markdown guide. When a model, equation, parameter definition, or scientific
claim changes, update the LaTeX manual and rebuild `manual/manual.pdf`.

The root `README.md` remains a concise landing page and should link to these
sources instead of duplicating them.
