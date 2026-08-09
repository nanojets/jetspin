# JETSPIN

JETSPIN is an open-source Fortran simulator for the dynamics of electrified
viscoelastic jets in electrospinning. It represents a jet as a sequence of
charged beads connected by viscoelastic elements and supports serial and MPI
execution.

The latest stable release is
[**JETSPIN 1.22**](https://github.com/nanojets/jetspin/releases/tag/v1.22).
Ongoing work is integrated on the `development` branch before it reaches
`master`.

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

The longer [Test Case 9 benchmark](docs/examples/test-9.md) is intentionally
excluded from the smoke and numerical-regression suites. It provides a fixed
1,000-bead workload for CPU/GPU performance measurements and optional
region-level profiling during OpenACC development. The current single-GPU
milestone offloads direct Coulomb interactions and, for the fixed Test 9
configuration, the three-dimensional equation-of-motion assembly including
its local curvature calculation. Test 9 keeps its state and RK4 scratch data
resident on the device across timesteps. See the
[OpenACC porting status](docs/introduction/openacc.md) for supported
configurations, numerical constraints, and remaining work.

## Citation

If JETSPIN contributes to published work, please cite:

> M. Lauricella, G. Pontrelli, I. Coluzza, D. Pisignano, and S. Succi,
> “JETSPIN: A specific-purpose open-source software for electrospinning
> simulations of nanofibers,” *Computer Physics Communications* **197**
> (2015), 227–238.

The evaporation implementation follows A. L. Yarin, S. Koombhongse, and
D. H. Reneker, *Journal of Applied Physics* **89** (2001), 3018–3026.
See the scientific manual for model scope and attribution.

## License and disclaimer

JETSPIN is distributed under the
[Open Software License 3.0](LICENSE.md). It is experimental research
software and is provided without warranty; users are responsible for
validating results for their applications.
