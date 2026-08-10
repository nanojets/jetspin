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
