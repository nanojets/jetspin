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
