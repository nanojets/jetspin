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

- [Test Case 1: one-dimensional reference case](examples/test-1.md)
- [Test Case 2: one-dimensional bead insertion](examples/test-2.md)
- [Test Case 3: three-dimensional PVP electrospinning](examples/test-3.md)
- [Test Case 4: gas counterflow](examples/test-4.md)
- [Test Case 5: dynamic refinement](examples/test-5.md)
- [Test Case 6: Kelvin–Voigt fluid](examples/test-6.md)
- [Test Case 7: evaporation and rotating electric field](examples/test-7.md)
- [Test Case 8: Yarin 2001 reference parameters](examples/test-8.md)
- [Test Case 9: 1,000-bead CPU/GPU benchmark](examples/test-9.md)
- [Test Case 16: dynamic Maxwell GPU evaporation](examples/test-16.md)
- [Test Case 17: dynamic Kelvin–Voigt GPU evaporation](examples/test-17.md)

## Scientific model

The mathematical-model pages will be added progressively. For now, consult
the [PDF manual](../manual/manual.pdf) or the modular
[LaTeX sources](../manual/).

## Documentation policy

Markdown pages prioritize concise operational guidance and links to working
examples. The LaTeX sources retain the full scientific discussion. Avoid
maintaining the same long passage in both formats; use the
[project overview](project-overview.md) to decide where a change belongs.
