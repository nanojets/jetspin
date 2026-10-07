# Test Case 9 numerical records

This directory preserves the numerical outputs associated with the initial
Test Case 9 CPU/GPU performance measurement. It is deliberately separate from
the normal regression suite because Test 9 takes much longer than Cases 1--8.

## Provenance

The exact compiler, module setup, and CPU/GPU flags are recorded in
[`../BUILD-PROVENANCE.md`](../BUILD-PROVENANCE.md). Use that record before
regenerating these files.

- source state: the repository state committed together with these records;
- source parent before the benchmark work: `53463bb4d05259da9c81bc11b79154be0ab263c9`;
- date: 2026-08-09 (records regenerated on 2026-10-07, see below);
- input: `examples/input-9/input.dat`;
- input SHA-256: `ec46d68f7850c03b28b0d1ae28f0d13d1621a103cf655f9c667f8e3e25fd5953`;
- compiler: NVFORTRAN 24.3-0, target `x86-64 Linux -tp znver3`;
- accelerator: NVIDIA A30, compute capability 8.0;
- accelerator build: `nvfortran-openacc GPUCC=80 CUDA_VERSION=12.3`;
- integration: 1,000 steps with 1,000 fixed beads;
- output sampling: initial state and every 200 steps.

The files are full sampled trajectories, not merely their final rows:

- `baseline-nvfortran-cpu.statout`;
- `baseline-openacc-a30.statout`.

The original files (2026-08-09) retained small CPU/GPU floating-point
differences. They were regenerated on 2026-10-07, when the leading bead,
which starts on the collector, became frozen there (see
`../BUILD-PROVENANCE.md`): the CPU (`nvfortran` target) and A30 records are
now identical. Their SHA-256 values are:

```text
d48ec723d6735e38cc2883dccf8cc522c5fd635b32c4db69f63945da9d260284  baseline-nvfortran-cpu.statout
d48ec723d6735e38cc2883dccf8cc522c5fd635b32c4db69f63945da9d260284  baseline-openacc-a30.statout
```

## Comparing future results

Pass a new CPU result, a new GPU result, or both:

```sh
sh tests/performance/test9/compare.sh /path/to/cpu/statout.dat
sh tests/performance/test9/compare.sh /path/to/cpu/statout.dat \\
    /path/to/gpu/statout.dat
```

The defaults are `rtol=1e-6` and `atol=1e-9`, matching the OpenACC
regression criterion. Override them only for diagnostic work:

```sh
JETSPIN_TEST9_RTOL=1e-7 JETSPIN_TEST9_ATOL=1e-10 \\
    sh tests/performance/test9/compare.sh cpu.statout gpu.statout
```

When both paths are supplied, the script checks the new CPU result against
the CPU baseline, the new GPU result against the A30 baseline, and the new GPU
result against the new CPU result. Performance times are intentionally not
stored in the numerical files and must be recorded separately because they
depend on hardware and system load.

Do not replace these records merely to make a comparison pass. If an
intentional numerical change requires new reference data, review every
observable and commit the updated records together with the responsible code
change.

## Development-oracle comparison

`nvfortran-openacc-force-oracle` evaluates every RK4 force stage through the
trusted CPU equations; `nvfortran-openacc-coulomb-oracle` moves only direct
Coulomb accumulation to the host. In the current NVFORTRAN 24.3/A30 check,
the standard GPU, complete-force oracle, and Coulomb-only oracle had worst
normalized differences of `7.97e-5`, `3.00e-7`, and `8.52e-5` from the CPU
record. This diagnostic indicates that direct-Coulomb summation order is not
the dominant Test 9 RK4 difference. Oracle runs are intentionally excluded
from performance comparisons.
