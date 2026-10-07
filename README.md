# JETSPIN

JETSPIN is an open-source Fortran simulator for the dynamics of electrified
viscoelastic jets in electrospinning. It represents a jet as a sequence of
charged beads connected by viscoelastic elements and supports serial, MPI,
and single-GPU OpenACC execution (NVIDIA HPC SDK).

The latest stable release is
[**JETSPIN 1.22**](https://github.com/nanojets/jetspin/releases/tag/v1.22).
Ongoing work is integrated on the `development` branch before it reaches
`master`; its first pre-release,
[**JETSPIN 2.0-beta.1**](https://github.com/nanojets/jetspin/releases/tag/v2.0-beta.1),
adds the single-GPU OpenACC build, reproducible stochastic and restarted
runs and the validation suites (see the [changelog](CHANGELOG.md)).

## Capabilities

JETSPIN includes one- and three-dimensional jet models, Maxwell,
Herschel–Bulkley, and Kelvin–Voigt rheology, aerodynamic drag, dynamic
refinement, time-dependent electric fields, stochastic forcing, and the
Yarin 2001 solvent-evaporation and rheological-solidification model.

## Quick start

From the repository root, compile the serial GFortran version:

```sh
make -C source -f ../build/Makefile gfortran
```

Then select an example and run it:

```sh
cp examples/input-1/input.dat execute/input.dat
cd execute
./main.x
```

See the [compilation guide](docs/introduction/compiling.md) and
[running guide](docs/introduction/running.md) for other compilers, MPI, and
runtime details.

## Documentation

- [Operational documentation](docs/README.md): build, execution, input,
  output, and worked examples in GitHub-friendly Markdown.
- [Project overview](docs/project-overview.md): architecture, repository
  conventions, maintenance context, and documentation policy.
- [Scientific user manual](manual/manual.pdf): authoritative equations,
  algorithms, references, and complete parameter tables.
- [LaTeX manual sources](manual/): modular source for the scientific manual.
- [Changelog](CHANGELOG.md): user-visible changes by release.

The Markdown guide intentionally avoids reproducing the complete scientific
manual. Operational guidance belongs in `docs/`; mathematical derivations
and the formal reference belong in `manual/`.

The Markdown documentation in `docs/`, the recent revisions of the LaTeX
manual, and the development records were written with the assistance of an
artificial-intelligence (AI) agent, working under the direction of the
JETSPIN developers.

## Validation

Run the local numerical smoke suite with:

```sh
tests/smoke/run.sh serial
```

Runtime-checking and MPI variants are also available:

```sh
tests/smoke/run.sh debug
tests/smoke/run.sh mpi
```

The same modes run in GitHub Actions, together with compilation of the
LaTeX manual.

The [numerical regression suite](docs/introduction/numerical-regression.md)
(`tests/regression/run.sh`) compares 1,000-step runs of Examples 1-8 with
versioned baselines and two MPI ranks with one, and `tests/restart/run.sh`
checks that a restarted run continues the uninterrupted one exactly.

The OpenACC build is checked by further suites: the restart check
(`tests/restart/run.sh openacc`), the refinement runners
(`tests/refinement/`), the dynamic evaporation validation
(`tests/performance/dynamic/validate_evaporation.sh`), and the CPU/GPU
records of Test Cases 9-25. The regression suite's `openacc` mode stays below
the device gate (at most 25 beads) and checks the OpenACC build's host path.

Every three-dimensional Euler, RK2 or RK4 run (Maxwell or Kelvin–Voigt) and
every Platen run (Maxwell), with or without evaporation, air drag, insertion,
removal and, for Platen, dynamic refinement, takes one device step once the
nozzle bead index reaches 100: the state stays on the GPU and an ordinary
step returns only a 20-byte topology record. Below that, and for the options
not yet ported, the OpenACC build runs the CPU code and reproduces the CPU
build byte for byte. See the
[OpenACC porting status](docs/introduction/openacc.md) for the coverage, the
numerical constraints, and the remaining work.

The test cases cover the device step from fixed 1,000-bead workloads to
growing jets:
[Test Cases 9-12](docs/examples/test-9.md) (RK4, Euler, RK2 and Platen with a
fixed bead set), [13-15](docs/examples/test-13.md) (insertion, removal and
capacity growth), [16-19](docs/examples/test-16.md) (Maxwell and
Kelvin–Voigt evaporation), [20](docs/examples/test-20.md) (stochastic Platen
with evaporation), [21-23](docs/examples/test-21.md) (dynamic Akima
refinement, capacity growth and collector removal). [Test Case
24](docs/examples/test-24.md) is a long production run that grows a jet from
a single nozzle bead until nozzle insertion balances collector removal, and
[Test Case 25](docs/examples/test-25.md) is the same input with Yarin
evaporation and a 50 % initial polymer fraction; its page explains why
Yarin's 6 % makes the dried jet unstable in this configuration.

## Citation

If JETSPIN contributes to published work, please cite:

> M. Lauricella, G. Pontrelli, I. Coluzza, D. Pisignano, and S. Succi,
> “JETSPIN: A specific-purpose open-source software for simulations of
> nanofiber electrospinning,” *Computer Physics Communications* **197**
> (2015), 227–238.

The models implemented in JETSPIN are reviewed, in the broader context of
electrospinning and solution blowing, in the following works, which may
also be cited:

> M. Lauricella, S. Succi, E. Zussman, D. Pisignano, and A. L. Yarin,
> “Models of polymer solutions in electrified jets and solution blowing,”
> *Reviews of Modern Physics* **92** (2020), 035004.

> A. L. Yarin, F. Pierini, E. Zussman, and M. Lauricella, “Modelling of
> nanofiber formation processes,” in *Materials and Electro-mechanical and
> Biomedical Devices Based on Nanofibers*, Springer (2024), 237–326.

The evaporation implementation follows A. L. Yarin, S. Koombhongse, and
D. H. Reneker, *Journal of Applied Physics* **89** (2001), 3018–3026.
See the scientific manual for model scope and attribution.

## License and disclaimer

JETSPIN is distributed under the
[Open Software License 3.0](LICENSE.md). It is experimental research
software and is provided without warranty; users are responsible for
validating results for their applications.
