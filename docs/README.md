# JETSPIN documentation

This directory provides a GitHub-friendly guide to building, running, and
configuring JETSPIN. The documentation is organized into short pages that
mirror the sections of the LaTeX manual.

The [PDF user manual](../manual/manual.pdf) remains the authoritative
scientific reference for the mathematical model, equations, bibliography,
and complete input/output tables.

For repository architecture, development status, validation expectations,
and documentation-maintenance rules, see the
[project overview](project-overview.md).

## Getting started

1. [Compile JETSPIN](introduction/compiling.md)
2. [Understand MPI parallelization](introduction/parallelization.md)
3. [Follow the OpenACC GPU port](introduction/openacc.md)
4. [Understand dynamic allocation and bead indexing](introduction/dynamic-allocation.md)
5. [Understand dynamic refinement](introduction/dynamic-refinement.md)
6. [Understand random numbers and MPI reproducibility](introduction/random-numbers.md)
7. [Run numerical regression tests](introduction/numerical-regression.md)
8. [Run JETSPIN](introduction/running.md)
9. [Prepare an input file](data/input.md)
10. [Read the output files](data/output.md)

## Example simulations

The [example index](examples/README.md) groups the cases and their
validation records.

- [Test Case 1: one-dimensional reference case](examples/test-1.md)
- [Test Case 2: one-dimensional bead insertion](examples/test-2.md)
- [Test Case 3: three-dimensional PVP electrospinning](examples/test-3.md)
- [Test Case 4: gas counterflow](examples/test-4.md)
- [Test Case 5: dynamic refinement](examples/test-5.md)
- [Test Case 6: Kelvin–Voigt fluid](examples/test-6.md)
- [Test Case 7: evaporation and rotating electric field](examples/test-7.md)
- [Test Case 8: Yarin 2001 reference parameters](examples/test-8.md)
- [Test Case 9: 1,000-bead CPU/GPU benchmark](examples/test-9.md)
- [Test Case 10: Euler CPU/GPU benchmark](examples/test-10.md)
- [Test Case 11: RK2 CPU/GPU benchmark](examples/test-11.md)
- [Test Case 12: Platen stochastic CPU/GPU benchmark](examples/test-12.md)
- [Test Case 13: dynamic-topology CPU/GPU benchmark](examples/test-13.md)
- [Test Case 14: larger bounded dynamic topology](examples/test-14.md)
- [Test Case 15: forced array reallocation](examples/test-15.md)
- [Test Case 16: dynamic Maxwell GPU evaporation](examples/test-16.md)
- [Test Case 17: dynamic Kelvin–Voigt GPU evaporation](examples/test-17.md)
- [Test Case 18: high-resolution dynamic Maxwell evaporation](examples/test-18.md)
- [Test Case 19: high-resolution dynamic Kelvin–Voigt evaporation](examples/test-19.md)
- [Test Case 20: stochastic Platen Maxwell evaporation](examples/test-20.md)
- [Test Case 21: anchored dynamic refinement](examples/test-21.md)
- [Test Case 22: repeated refinement and capacity growth](examples/test-22.md)
- [Test Case 23: refinement with collector removal](examples/test-23.md)
- [Test Case 24: long run to a stationary bead count](examples/test-24.md)
- [Test Case 25: evaporative counterpart of Test 24](examples/test-25.md)

## Test suites

- [Smoke tests](../tests/smoke/README.md): Examples 1-8, serial, runtime
  checks and MPI (also run by GitHub Actions)
- [Numerical regression](../tests/regression/README.md): versioned baselines,
  two-rank MPI comparison, OpenACC mode
- [Exact restart](../tests/restart/README.md): restarted runs against
  uninterrupted ones, CPU, GPU and MPI
- [Dynamic refinement](../tests/refinement/README.md): Tests 21-23 with their
  oracles and comparison builds
- [Evaporation](../tests/evaporation/README.md): Yarin 2001 and Kelvin–Voigt
  evaporation checks
- CPU/GPU records: [Test 9](../tests/performance/test9/README.md),
  [Tests 10-12 and 20](../tests/performance/integrators/README.md),
  [dynamic topology and evaporation, Tests 13-19](../tests/performance/dynamic/README.md),
  and their [build provenance](../tests/performance/BUILD-PROVENANCE.md)

The development log, with the evidence behind every change of the GPU port,
is [`STATE.md`](STATE.md).

## Investigation notes

- [Refinement robustness under repeated remeshing](refinement-robustness-investigation.md)

## Scientific model

The mathematical-model pages will be added progressively. For now, consult
the [PDF manual](../manual/manual.pdf) or the modular
[LaTeX sources](../manual/).

## Documentation policy

Markdown pages prioritize concise operational guidance and links to working
examples. The LaTeX sources retain the full scientific discussion. Avoid
maintaining the same long passage in both formats; use the
[project overview](project-overview.md) to decide where a change belongs.
