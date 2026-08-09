# OpenACC GPU porting status

JETSPIN's NVIDIA GPU port is being developed incrementally with explicit
OpenACC data regions. The normal GFortran, Intel, and MPI targets continue to
select the original CPU implementation.

## Implemented milestone

The `nvfortran-openacc` target currently accelerates the direct Coulomb
summation when all of these conditions hold:

- execution is serial (`mxrank == 1`);
- multiple-step Coulomb summation is disabled;
- evaporation is disabled;
- the system is either one- or three-dimensional.

The build defaults to `GPUCC=80` and `CUDA_VERSION=12.3`, producing
`-gpu=cc80,cuda12.3`. This covers the NVIDIA A30 and A100. Select a different
architecture without the `cc` prefix, for example:

```sh
make -C source -f ../build/Makefile nvfortran-openacc \
  GPUCC=90 CUDA_VERSION=12.3
```

Unsupported configurations retain the existing CPU path. In particular, MPI,
the multiple-step/neighbour-list algorithm, and the evaporation-specific
Coulomb routine are not GPU kernels yet.

The accelerator kernel assigns one target bead to each parallel iteration.
That iteration visits every other active bead and writes only its target
force. This target-centric formulation is race-free and needs no atomics. It
evaluates each physical pair twice instead of sharing a pair contribution as
the CPU implementation does, so its floating-point accumulation order is
different.

The implementation preserves the existing model conventions, including:

- charge and mass scaling;
- the softening cross section associated with the higher-index bead;
- frozen-bead exclusion;
- the one-dimensional distance cutoff;
- mirror-charge contributions and their three-dimensional cutoff.

## Explicit data region

Coordinates, bead properties, cross sections, and the Coulomb force array are
named explicitly in each OpenACC data region. No managed/unified-memory build
mode is used.

At this milestone, the time integrators and dynamic bead operations still run
on the CPU. Coordinates can change at every integrator stage and the computed
force is immediately consumed by host code. Each Coulomb evaluation therefore
copies current inputs to the device and the result back to the host. Keeping
these arrays permanently resident now would not remove those transfers and
would make host reallocation easier to mishandle.

This call-scoped data region also makes the current implementation safe when
insertion, removal, or dynamic refinement changes `mxnpjet`: the next call
maps the newly allocated host arrays and their new capacity.

## Numerical validation

Compare the OpenACC path with a CPU reference built by the same `nvfortran`
version. Different Fortran compilers may use different pseudo-random-number
sequences, so a GFortran trajectory is not a suitable direct reference for a
stochastic NVIDIA build.

The direct CPU algorithm accumulates one pair into both beads, whereas the
race-free accelerator algorithm accumulates independently by target bead.
Their results should therefore be compared numerically, not as binary files.
The initial development checks use the regression-suite criterion with
`rtol=1e-6` and `atol=1e-9`. Run the complete CPU/GPU comparison with:

```sh
tests/regression/run.sh openacc
```

The initial 1,000-step validation on an NVIDIA A30 passed all eight cases;
the written observables matched the NVFORTRAN CPU results at their output
precision.

`nvfortran-openacc-host` is useful for checking the accelerated control path
and dynamic allocation without a GPU. It does not replace the real-device
regression above and cannot establish GPU performance.

## Next porting stages

1. Port the evaporation-specific direct Coulomb kernel.
2. Move the equation-of-motion and integrator stages into explicit device
   regions, keeping the main bead state resident across timesteps.
3. Add explicit device teardown/recreation hooks around capacity changes and
   synchronize only topology metadata and requested output fields.
4. Port the local Akima coefficient loops, replace the interpolation interval
   scan with a GPU-suitable search, and then address dynamic refinement.
5. Evaluate one-GPU-per-rank MPI execution only after the single-GPU numerical
   path is stable.

The intended steady state is a persistent device-resident simulation with
host transfers for output, checkpoints, and topology changes—not a transfer
of the complete jet at every timestep.
