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

The development-only `nvfortran-openacc-force-oracle` and
`nvfortran-openacc-coulomb-oracle` targets use the same interfaces here as in
the evaporation tests. For Euler their worst normalized differences from the
NVFORTRAN CPU record were `7.00e-7` and `5.00e-7`; for RK2 they were
`6.00e-10` and `1.00e-7`. Both targets add per-stage host/device transfers and
must not be used for the timing table above.
The same two oracle interfaces cover Platen with and without evaporation;
Tests 12 and 20 both have a maximum written-output difference of zero.

Compare new outputs against the records and against each other with:

```sh
tests/performance/integrators/compare.sh euler CPU/statout.dat GPU/statout.dat
tests/performance/integrators/compare.sh rk2 CPU/statout.dat GPU/statout.dat
tests/performance/integrators/compare.sh platen CPU/statout.dat GPU/statout.dat
```

Set `RTOL` or `ATOL` in the environment to test a different tolerance.
