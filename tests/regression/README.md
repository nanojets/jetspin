# Numerical regression tests

This suite runs shortened, deterministic versions of all eight examples for
1,000 integration steps. It compares serial output against a versioned
baseline and a two-rank MPI run against the matching serial run.

Each case retains its original timestep. The generated final time is exactly
1,000 timesteps and output is sampled every 200 timesteps, giving six rows
including the initial state.

Regression Case 4 starts with 25 beads. This deliberately exceeds the
10-bead minimum MPI chunk on two ranks and exercises stochastic Platen draws
on both ranks rather than merely launching an idle second process.

The inputs retain their model-specific directives, but use a common
`statout.dat` schema without CPU timings. `nstep`, bead count `n`, and
refinement count `nref` are compared exactly. Continuous observables use:

```text
abs(actual - reference) <= atol + rtol * abs(reference)
```

Run the suite from the repository root:

```sh
tests/regression/run.sh gfortran
```

The same cases can be compiled entirely with NVIDIA HPC SDK. Load the
`nvhpc/24.3` module first so both `nvfortran` and `mpif90` select NVIDIA
Fortran:

```sh
tests/regression/run.sh nvfortran
```

This uses compiler-specific baselines under `baselines/nvfortran/` and then
compares the two-rank NVFORTRAN/MPI result with the matching serial result.
The compiler-specific baseline is necessary because Fortran compilers may
provide different `random_number` sequences.

Compare the OpenACC executable with the matching NVFORTRAN CPU trajectory on
an accessible NVIDIA GPU with:

```sh
tests/regression/run.sh openacc
```

| Comparison | Relative tolerance | Absolute tolerance |
| --- | ---: | ---: |
| Serial versus baseline | `1e-7` | `1e-10` |
| Two-rank MPI versus serial | `1e-6` | `1e-9` |
| OpenACC versus NVFORTRAN CPU, cases 1--7 | `1e-6` | `1e-9` |
| OpenACC versus NVFORTRAN CPU, case 8 (Maxwell evaporation) | `3e-2` | `1e-8` |

Override the defaults without editing the test:

```sh
JETSPIN_REGRESSION_RTOL=1e-8 \
JETSPIN_REGRESSION_ATOL=1e-11 \
JETSPIN_MPI_REGRESSION_RTOL=1e-7 \
JETSPIN_MPI_REGRESSION_ATOL=1e-10 \
JETSPIN_OPENACC_REGRESSION_RTOL=1e-7 \
JETSPIN_OPENACC_REGRESSION_ATOL=1e-10 \
tests/regression/run.sh nvfortran
```

Case 8 has a separate default because the evaporating Maxwell trajectory is
chaotic to floating-point perturbations caused by the GPU Coulomb reduction
order. Override it independently when tighter or looser acceptance is needed:

```sh
JETSPIN_OPENACC_EVAPORATION_RTOL=3e-2 \
JETSPIN_OPENACC_EVAPORATION_ATOL=1e-8 \
tests/regression/run.sh openacc
```

This is an intentional numerical tolerance, not a baseline update: the GPU
and host Coulomb forces agree to approximately `1e-9` per component before
trajectory-level error amplification. The dynamic event counters `n`, `curn`,
and `curc` are omitted from the OpenACC case-8 comparison because a tiny
floating-point perturbation can move an insertion/removal event by one output
interval; continuous observables remain checked.

Regenerate baselines intentionally after an accepted numerical change with:

```sh
tests/regression/run.sh update
tests/regression/run.sh nvfortran update
```

Baseline updates must be reviewed as scientific results. Do not update them
merely to make an unexplained regression pass. Set
`JETSPIN_REGRESSION_KEEP=1` to retain temporary runs for investigation, and
`JETSPIN_REGRESSION_TIMEOUT` to change the 30-second per-case timeout. For
diagnosis, `JETSPIN_REGRESSION_FIRST_CASE` and
`JETSPIN_REGRESSION_LAST_CASE` select an inclusive subset of cases.
On systems with multiple MPI installations, select matching GFortran/OpenMPI
wrappers explicitly with `JETSPIN_MPIFC` and `MPIEXEC`.
