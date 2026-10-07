# Fixed 1,000-bead integrator benchmarks

The compiler version, module initialization, and exact CPU/GPU flags for
these records are maintained in
[`../BUILD-PROVENANCE.md`](../BUILD-PROVENANCE.md).

Tests 10 and 11 reuse the fixed 1,000-bead geometry and physics of Test 9,
but select the explicit Euler and second-order Heun integrators respectively.
Each case executes 1,000 steps and samples the trajectory every 200 steps.
Test 12 uses the same fixed geometry with the Platen stochastic integrator and
the stochastic air-drag parameters from Test 4. Its complete Gaussian history
is generated on the CPU before loop timing and copied to the GPU once.
Test 20 is the corresponding 100-step Maxwell evaporation Platen check.

The versioned files record a paired run built with NVFORTRAN 24.3 and run on
one NVIDIA A30 (`cc80`, CUDA 12.3, `nofma`). The measured integration-loop
times were:

| Case | NVFORTRAN CPU | OpenACC A30 | Speedup |
| --- | ---: | ---: | ---: |
| Test 10, Euler | `9.787760 s` | `0.741555 s` | `13.20x` |
| Test 11, RK2/Heun | `19.464682 s` | `1.254380 s` | `15.52x` |
| Test 12, Platen | `29.050473 s` | `1.874810 s` | `15.50x` |
| Test 20, Maxwell evaporation Platen | `3.177595 s` | `0.224058 s` | `14.18x` |

Both CPU/GPU comparisons pass with `rtol=1e-6` and `atol=1e-9`. The worst
normalized differences were `4.00e-7` for Euler and `8.99e-7` for RK2.
The Platen outputs agree exactly in all written columns.
For Test 20, the standard GPU path and both oracle paths also agree exactly
with the fresh NVFORTRAN CPU output in all six rows and fourteen columns.
These timings describe this development system and are not portable
performance guarantees.

Since 2026-10-05 the non-evaporative Coulomb kernel gives each target bead one
gang, reduces over the sources across its vector lanes, and no longer uploads
array descriptors at every call. On one A30 bound to its NUMA node, against the
previous kernel on the same node, the loop times are `0.122 s` against
`0.589 s` (Test 10), `0.223 s` against `1.078 s` (Test 11), and `0.354 s`
against `1.659 s` (Test 12); Test 20 uses the
evaporative kernel and is unchanged. All three still pass `compare.sh`
against the CPU records (worst normalized differences `8e-7`, `9e-7`, and
`0`).

The development-only `nvfortran-openacc-force-oracle` and
`nvfortran-openacc-coulomb-oracle` targets use the same interfaces here as in
the evaporation tests. For Euler their worst normalized differences from the
NVFORTRAN CPU record were `7.00e-7` and `5.00e-7`; for RK2 they were
`6.00e-10` and `1.00e-7`. Both targets add per-stage host/device transfers and
must not be used for the timing table above.
The same two oracle interfaces cover Platen with and without evaporation;
Tests 12 and 20 both have a maximum written-output difference of zero.

On 2026-10-07 the records of Tests 10-12 were regenerated: a bead that
reaches the collector is now frozen there, and the leading bead of these jets
starts on the collector (before, it crossed the plane by about 6e-4 cm in
1,000 steps). With the common device step (Euler 10, RK2 14, Platen 13 kernel
launches per step) the A30 records equal the CPU records row for row, and
NVFORTRAN 25.5 gives the same rows; provenance in
[`../BUILD-PROVENANCE.md`](../BUILD-PROVENANCE.md). The timings and
normalized differences above were measured with the earlier records. Test 20
has no versioned record: its check compares the CPU and A30 outputs of the
same build.

Compare new outputs against the records and against each other with:

```sh
tests/performance/integrators/compare.sh euler CPU/statout.dat GPU/statout.dat
tests/performance/integrators/compare.sh rk2 CPU/statout.dat GPU/statout.dat
tests/performance/integrators/compare.sh platen CPU/statout.dat GPU/statout.dat
```

Set `RTOL` or `ATOL` in the environment to test a different tolerance.
