# Compiling JETSPIN

JETSPIN requires a Fortran compiler. Parallel builds also require an MPI
implementation. Run all commands below from the repository root; the
Makefile does not need to be copied into `source/`.

## Serial build with GFortran

```sh
make -C source -f ../build/Makefile gfortran
```

Object and module files are created in `source/`. The resulting executable
is written to `execute/main.x`.

## OpenMPI build with GFortran

```sh
make -C source -f ../build/Makefile gfortran-mpi
```

Run `mpif90 --showme:command` to confirm that the MPI wrapper uses GFortran
when selecting this target.

## NVIDIA HPC SDK builds

The shell initialization adds the NVIDIA module directory. Start from a clean
module environment and load the installed SDK with:

```sh
module purge
module use /opt/nvidia/hpc_sdk/modulefiles
module load nvhpc/24.3
```

The installed module is named `nvhpc/24.3`, not `hpcsdk`. Confirm the compiler
and MPI wrapper selection with:

```sh
nvfortran --version
mpif90 --showme:command
```

The second command must report `nvfortran`.

Build the serial CPU executable with NVIDIA Fortran, without OpenACC
offloading:

```sh
make -C source -f ../build/Makefile nvfortran
```

Build the MPI CPU executable using the HPC SDK `mpif90` wrapper, again without
OpenACC offloading:

```sh
make -C source -f ../build/Makefile nvfortran-mpi
```

Build with OpenACC offloading. The defaults target compute capability 8.0
(including NVIDIA A30 and A100) and the CUDA 12.3 toolkit shipped with HPC
SDK 24.3. The target also selects `nofma` to preserve the numerical behaviour
of the nearly straight three-point curvature calculation:

```sh
make -C source -f ../build/Makefile nvfortran-openacc
```

The OpenACC GPU and host-validation targets expose the standard `_OPENACC`
preprocessing macro. Imports and calls that initialize or select accelerator
code are enclosed by this macro. CPU targets preprocess the same source
without defining it, so those calls are absent from the compiled CPU
translation units rather than being executed as no-op runtime branches.

Select another GPU compute capability or CUDA toolkit through Make variables.
Pass the numeric capability without the `cc` prefix; the Makefile constructs
the NVFORTRAN option:

```sh
make -C source -f ../build/Makefile nvfortran-openacc \
  GPUCC=90 CUDA_VERSION=12.3
```

HPC SDK 25.5 ships the CUDA 11.8 and 12.9 toolkits, not 12.3: build with
`CUDA_VERSION=12.9`. Such a build also runs with a driver of the CUDA 12.3
generation (checked on A30 GPUs, driver 545). NVFORTRAN 24.3 at `-O3`
rewrites a division `x/y` as `x*(1/y)` (`-Mrecip-div`; `-Mno-recip-div`
turns it off), while 25.5, like GFortran, rounds divisions exactly. The two
compilers therefore give results that differ in the last digits, and
trajectories that are sensitive to them (insertion and removal thresholds,
bending instability) separate. Compare an OpenACC run with a CPU build of
the same compiler version.

Use `nvaccelinfo` and `nvidia-smi` to check that the compiler runtime can see
the selected GPU before running the resulting executable.

When timing a GPU run, bind it to the NUMA node of its GPU (`nvidia-smi topo
-m` lists the affinity), for example `numactl --cpunodebind=N --membind=N
./main.x` with `N` that node: placement of the host process alone changes
latency-bound runs by up to 8 %.

For compiler and numerical validation without an accessible NVIDIA device:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-host
```

This executes the OpenACC code path on the host. It does not measure GPU
performance. See [OpenACC porting status](openacc.md) for the implemented
kernel, data movement, limitations, and validation expectations.

Two model-independent development targets isolate force-kernel numerical
differences in the device step: Euler, RK2, RK4, and stochastic Platen,
Maxwell and Kelvin–Voigt, with or without evaporation, fixed bead sets and
jets with insertion and refinement:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-force-oracle GPUCC=80
make -C source -f ../build/Makefile nvfortran-openacc-coulomb-oracle GPUCC=80
```

`nvfortran-openacc-force-oracle` enables
`JETSPIN_DEV_HOST_FORCE_ORACLE`. At every force evaluation of the device
step it downloads the current stage, evaluates the complete trusted CPU
force equations, and uploads the derivatives. Integration updates, dynamic
topology and the tail of the Platen step (`accelerator_platen_update`,
`accelerator_platen_end_step`, with the final stress derivative) stay on
the GPU; since milestone M3 (2026-10-06) the tail is not oracled, while the
oracle builds evaluated the final stress derivative on the host before.

`nvfortran-openacc-coulomb-oracle` enables only
`JETSPIN_DEV_HOST_COULOMB_ORACLE`. It downloads the state required by the direct
Coulomb sum, evaluates that sum on the CPU in its established order, uploads
`ycf`, and keeps every other force term on the GPU. The two macros are not
aliases: the first is a complete force oracle and the second isolates only
Coulomb accumulation. Their deliberate per-stage transfers make both targets
unsuitable for production or performance measurements.

Two additional development targets isolate dynamic-refinement interpolation:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-host-akima GPUCC=80
make -C source -f ../build/Makefile nvfortran-openacc-compare-akima GPUCC=80
```

`nvfortran-openacc-host-akima` enables `JETSPIN_DEV_HOST_AKIMA` and retains
the historical host coefficient/interpolation path while the surrounding
OpenACC refinement lifecycle remains active. `nvfortran-openacc-compare-akima`
enables `JETSPIN_COMPARE_AKIMA`: it computes the trusted host result followed
by the device result and reports field-wise maximum coefficient and value
differences. These are numerical-validation builds, not performance builds.

The former `nvfortran-openacc-dynamic-platen` target, which opened the
persistent device path to non-evaporative Platen runs growing from a single
bead, was removed on 2026-10-05: `nvfortran-openacc` now does this itself
(see [OpenACC](openacc.md)).

## Runtime environment variables

The OpenACC build reads these variables at run time (see
[OpenACC](openacc.md) for the device step they refer to); the last four
also act in the CPU builds. A switch is on when its value is `1`.

| Variable | Effect |
| --- | --- |
| `JETSPIN_OPENACC_COULOMB_MIN_BEADS=<n>` | Fewest active beads for which a run outside the device step offloads a Coulomb sum, one-dimensional (non-evaporative) or three-dimensional; default 128, 0 offloads every call. Inside the device step every sum runs on the GPU. |
| `JETSPIN_OPENACC_SYNC` | Runs the Platen device step of a run with insertion synchronously instead of on one asynchronous queue. The RK device steps and fixed bead sets are always synchronous. |
| `JETSPIN_OPENACC_DISABLE_PERSISTENT` | Keeps the device step closed for every run, fixed bead sets included: the OpenACC build runs the CPU build's code, offloading only Coulomb sums of at least `JETSPIN_OPENACC_COULOMB_MIN_BEADS` beads. The Gaussian-pool decision is unchanged. |
| `JETSPIN_OPENACC_DISABLE_EOM` | Keeps the device step closed, as the previous one, since 2026-10-07; before, the device stage skipped the equations of motion and integrated stale derivatives. A device EOM kernel that refuses a stage now stops the run with error 22. |
| `JETSPIN_OPENACC_DISABLE_COULOMB` | Keeps the one-dimensional Coulomb sum on the host; no effect on the three-dimensional sums. |
| `JETSPIN_PROFILE` | Prints at the end of the run the wall time and number of calls of the timed parts of the loop (integrator, Coulomb sums, insertion, removal, statistics, output, restart). Also `yes`, `true`, `on`. |
| `JETSPIN_TOPOLOGY_SNAPSHOT` | Writes the active state at every insertion and removal to `topology-state.dat` (downloaded from the device first in the OpenACC build). |
| `JETSPIN_REFINEMENT_INITIAL_RESERVE=<n>` | Spare capacity allocated at the start of a run with insertion and tagged beads (refinement) that starts with at least 100 beads; default the allocation increment, 100. A small value reaches capacity growth in minutes instead of hours. |
| `JETSPIN_REFINEMENT_GROWTH_INCREMENT=<n>` | Capacity added when an accepted refinement event outgrows the arrays; default the allocation increment. |

Development macros are passed through `FPPFLAGS_EXTRA` ([below](#make-variables));
the oracle and comparison targets set theirs. To compute the Coulomb sums on
the host while the rest of the step stays on the device, use the
`nvfortran-openacc-coulomb-oracle` target (the macro
`JETSPIN_DISABLE_COULOMB_EVAP`, which did this for the evaporative 3-D sum
only and read a stale host state inside the device step, was removed on
2026-10-07).

## Make variables

Pass these on the `make` command line, as `GPUCC` above. Make also reads
those with a default set by `?=` (all but `EX` and `BINROOT`) from the
environment, so an unrelated `CUDA_VERSION` exported by another tool
reaches the build.

| Variable | Default | Used by |
| --- | --- | --- |
| `EX` | `main.x` | every target: name of the executable |
| `BINROOT` | `../execute`, relative to `source/` | every target: directory that receives the executable |
| `MPIFC` | `mpif90` | `gfortran-mpi`, `gfortran-mpidebugger`, `nvfortran-mpi`: MPI compiler wrapper |
| `MPI_COMPAT_FLAG` | `-fallow-argument-mismatch` | `gfortran-mpi`, `gfortran-mpidebugger`; a GFortran older than 10 needs `-Wno-argument-mismatch` instead, as `tests/regression/run.sh` selects |
| `NVFORTRAN` | `nvfortran` | `nvfortran` and the `nvfortran-openacc*` targets |
| `GPUCC` | `80` | the OpenACC GPU targets (not `nvfortran-openacc-host`): compute capability without `cc` |
| `CUDA_VERSION` | `12.3` | the OpenACC GPU targets; `12.9` with HPC SDK 25.5 |
| `FPPFLAGS_EXTRA` | empty | `nvfortran`, `nvfortran-mpi` and the `nvfortran-openacc*` targets only: extra preprocessor flags (development macros) |

The `gfortran`, `gfortran-debugger`, `intel*` and `cygwin*` targets name
their compilers directly.

## Other targets

| Target | Purpose |
| --- | --- |
| `gfortran` | Optimized serial GFortran build |
| `gfortran-mpi` | Optimized OpenMPI/GFortran build |
| `gfortran-debugger` | Serial build with runtime checks |
| `gfortran-mpidebugger` | MPI build with runtime checks |
| `intel` | Optimized serial Intel Fortran build |
| `intel-mpi` | Intel MPI build |
| `intel-openmpi` | Intel Fortran with OpenMPI |
| `intel-debugger` | Serial Intel Fortran build with runtime checks and traceback |
| `intel-mpidebugger` | Intel MPI build with runtime checks and traceback |
| `cygwin`, `cygwin-mpi` | Windows/Cygwin builds |
| `nvfortran` | Optimized serial CPU build with NVIDIA Fortran |
| `nvfortran-mpi` | MPI CPU build with the HPC SDK NVFORTRAN wrapper |
| `nvfortran-openacc` | Single-GPU OpenACC build; configurable with `GPUCC` and `CUDA_VERSION` |
| `nvfortran-openacc-host` | OpenACC code-path validation on the host |
| `nvfortran-openacc-force-oracle` | Development-only complete host-force oracle for all supported integrators, with or without evaporation |
| `nvfortran-openacc-coulomb-oracle` | Development-only host direct-Coulomb oracle, with or without evaporation |
| `nvfortran-openacc-host-akima` | Development-only historical host-Akima oracle inside the OpenACC refinement path |
| `nvfortran-openacc-compare-akima` | Development-only host/device Akima A/B comparison |
| `nvfortran-openacc-compare-refinement` | Development-only host/device A/B comparison of the accepted-event volume, mass, and charge reconstruction (`JETSPIN_COMPARE_REFINEMENT_ASSEMBLY`) |
| `help` | Display available targets |
| `clean` | Remove objects and module files from `source/` |

For example, clean a build with:

```sh
make -C source -f ../build/Makefile clean
```

The complete build rules are in [`build/Makefile`](../../build/Makefile).
See also the corresponding [LaTeX section](../../manual/compiling.tex).
Benchmark-specific compiler provenance and exact effective flags are recorded
in [`tests/performance/BUILD-PROVENANCE.md`](../../tests/performance/BUILD-PROVENANCE.md).
