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

Use `nvaccelinfo` and `nvidia-smi` to check that the compiler runtime can see
the selected GPU before running the resulting executable.

For compiler and numerical validation without an accessible NVIDIA device:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-host
```

This executes the OpenACC code path on the host. It does not measure GPU
performance. See [OpenACC porting status](openacc.md) for the implemented
kernel, data movement, limitations, and validation expectations.

Two model-independent development targets isolate force-kernel numerical
differences in the supported serial three-dimensional Euler, RK2, RK4, and
fixed-topology stochastic Platen paths, with or without evaporation:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-force-oracle GPUCC=80
make -C source -f ../build/Makefile nvfortran-openacc-coulomb-oracle GPUCC=80
```

`nvfortran-openacc-force-oracle` enables
`JETSPIN_DEV_HOST_FORCE_ORACLE`. It downloads the current RK stage, evaluates
the complete trusted CPU force equations, uploads the derivatives, and leaves
integration updates and dynamic topology on the GPU. The same interface
covers non-evaporative simulations, Maxwell and Kelvin–Voigt deterministic
evaporation, and Maxwell stochastic Platen evaporation.

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
| `cygwin`, `cygwin-mpi` | Windows/Cygwin builds |
| `nvfortran` | Optimized serial CPU build with NVIDIA Fortran |
| `nvfortran-mpi` | MPI CPU build with the HPC SDK NVFORTRAN wrapper |
| `nvfortran-openacc` | Single-GPU OpenACC build; configurable with `GPUCC` and `CUDA_VERSION` |
| `nvfortran-openacc-host` | OpenACC code-path validation on the host |
| `nvfortran-openacc-force-oracle` | Development-only complete host-force oracle for all supported integrators, with or without evaporation |
| `nvfortran-openacc-coulomb-oracle` | Development-only host direct-Coulomb oracle, with or without evaporation |
| `nvfortran-openacc-host-akima` | Development-only historical host-Akima oracle inside the OpenACC refinement path |
| `nvfortran-openacc-compare-akima` | Development-only host/device Akima A/B comparison |
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
