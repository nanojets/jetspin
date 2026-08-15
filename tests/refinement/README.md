# Dynamic-refinement validation

`run.sh` builds an isolated executable, runs Test Case 21, and checks the
accepted Akima-remeshing event. The checks cover initial anchor creation,
increased resolution, preservation of existing anchor positions, velocities,
stress, and both radii, strict path ordering, conservation of mass, charge,
reference volume, and evaporated volume, finite statistics, and a non-zero
final refinement count.

Run it from the repository root:

```sh
tests/refinement/run.sh
```

The default is the GFortran CPU reference. The same event-level checks cover
the NVFORTRAN CPU, native OpenACC, and development host-force-oracle paths:

```sh
tests/refinement/run.sh nvfortran
tests/refinement/run.sh openacc
tests/refinement/run.sh force-oracle
```

Append `capacity-growth` to any backend to reduce only the developer reserve
from 100 to 50 entries. The accepted Akima mesh must then exceed the original
capacity, which validates topology/evaporation release and rebind, Platen
scratch reallocation, and Gaussian-history stride resizing:

```sh
tests/refinement/run.sh nvfortran capacity-growth
tests/refinement/run.sh openacc capacity-growth
tests/refinement/run.sh force-oracle capacity-growth
```

The NVIDIA variants require an NVHPC environment; `GPUCC` defaults to `80`.
The native GPU trajectory is compared through refinement invariants rather
than binary identity because target-centric Coulomb accumulation changes the
floating-point summation order. The force oracle is for diagnosis only.

The test takes roughly half a minute on the development host. Set
`JETSPIN_REFINEMENT_KEEP=1` to retain its temporary build and output directory.

## Repeated-refinement stress test

`run_test22.sh` executes the longer Test Case 22 workload. It requires exactly
three accepted remeshing events and three capacity increases, and checks every
event for anchor preservation, ordered coordinates, separate reference and
evaporated-volume conservation, mass/charge conservation, finite output, and
a consistent repacked Gaussian history.

This is intentionally a same-compiler NVFORTRAN CPU/GPU test. The chaotic,
stochastic trajectory is not expected to accept the same number of events at
the same final time with a different compiler; the single-event Test 21 runner
continues to provide the GFortran refinement check.

```sh
tests/refinement/run_test22.sh nvfortran
tests/refinement/run_test22.sh openacc
tests/refinement/run_test22.sh force-oracle
```

The runner sets the developer-only initial reserve and refinement growth
increment to 20. This deliberately causes repeated rebinds; production values
remain 100. For accelerator builds the runner additionally requires exactly
three `OpenACC refinement capacity rebind` diagnostics. Its default timeout is
180 seconds and can be changed with `JETSPIN_REFINEMENT_TIMEOUT`.

## Refinement and collector-removal test

`run_test23.sh` extends the repeated-growth workload to the collector. It
requires exactly three accepted remeshes and capacity changes, plus at least
one device topology removal before and after the final remesh. The same
anchor, conservation, ordering, Gaussian-history, and finite-output checks
are applied at every refinement event.

```sh
tests/refinement/run_test23.sh nvfortran
tests/refinement/run_test23.sh openacc
tests/refinement/run_test23.sh force-oracle
tests/refinement/run_test23.sh host-akima
tests/refinement/run_test23.sh akima-compare
```

This is also a same-compiler NVFORTRAN CPU/GPU test. The exact later removal
schedule is intentionally not compared because the different direct-Coulomb
summation order is amplified by the stochastic bending trajectory near the
collector. The runner's default timeout is 240 seconds to accommodate the
development force oracle.

The `host-akima` backend keeps the historical CPU spline as an oracle while
the rest of the accelerator lifecycle remains enabled. The `akima-compare`
backend evaluates both spline implementations for every field at every event.
It requires 33 comparison records (11 fields times 3 events), finite metrics,
coefficient relative error at most `1e-12`, interpolated absolute error at
most `1e-9`, and interpolated relative error at most `1e-12`.
