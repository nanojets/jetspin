# NVIDIA benchmark build provenance

The numerical and timing records for Tests 9--13 are tied to the following
NVHPC environment. Record these values with every future replacement of a
benchmark baseline.

```text
Compiler: nvfortran 24.3-0 64-bit target on x86-64 Linux -tp znver3
GPU: NVIDIA A30, compute capability 8.0
CUDA toolkit selected by Makefile: 12.3
GPU device flags: -gpu=cc80,cuda12.3,nofma
```

Initialize the modules with:

```sh
module purge
module use /opt/nvidia/hpc_sdk/modulefiles
module use --append "$HOME/modulefiles"
module load nvhpc/24.3
nvfortran --version
nvidia-smi -L
```

The CPU-side NVFORTRAN reference path used for paired benchmarks is the
OpenACC host-validation target:

```text
Compile flags: -O2 -acc=host -Minfo=accel -Mpreprocess
Link flags:    -O2 -acc=host
Target:        nvfortran-openacc-host
```

The GPU path is:

```text
Compile flags: -O3 -acc=gpu -gpu=cc80,cuda12.3,nofma -Minfo=accel -Mpreprocess
Link flags:    -O3 -acc=gpu -gpu=cc80,cuda12.3,nofma
Target:        nvfortran-openacc
```

Build commands from the repository root:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-host
make -C source -f ../build/Makefile nvfortran-openacc \
  GPUCC=80 CUDA_VERSION=12.3
```

`-tp znver3` is the target reported by this NVFORTRAN installation; it is
recorded even though it is supplied by the compiler environment rather than
explicitly by the repository Makefile.

## Records of Tests 9-12 regenerated on 2026-10-07

Since 2026-10-07 a bead that reaches the collector is frozen there in every
run, and the leading bead of the fixed 1,000-bead jets starts on the
collector: the records of Tests 9-12 (`test9/` and `integrators/`) were
regenerated with the working tree of that date, on the standard inputs:

```text
Compiler: nvfortran 24.3-0 64-bit target on x86-64 Linux -tp icelake-server
GPU: NVIDIA A30, compute capability 8.0, driver of the CUDA 12.3 generation
CPU reference target: nvfortran       (-O3 -Mpreprocess)
GPU target:           nvfortran-openacc (-O3 -acc=gpu -gpu=cc80,cuda12.3,nofma)
```

The CPU and A30 records are identical, row for row, and NVFORTRAN 25.5
(`CUDA_VERSION=12.9`) gives the same rows. The timings quoted in the README
files were measured with the earlier records and builds.

## NVHPC 25.5

NVHPC 25.5 is supported since 2026-10-07; 24.3 remains the reference
compiler of these records. 25.5 ships the CUDA 11.8 and 12.9 toolkits, so
build its OpenACC targets with `CUDA_VERSION=12.9` (the Makefile default is
`12.3`):

```sh
make -C source -f ../build/Makefile nvfortran-openacc \
  GPUCC=80 CUDA_VERSION=12.9
```

25.5 rounds divisions exactly, while 24.3 at `-O3` multiplies by the
reciprocal (`-Mrecip-div`), so the two compilers' trajectories differ at
roundoff level and chaotic runs separate. Compare a new record only with
records of the same compiler. On Tests 24 and 25, 25.5 was also about 6 %
slower than 24.3.
