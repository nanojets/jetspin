# Euler and RK2 1,000-bead benchmarks

Tests 10 and 11 reuse the fixed 1,000-bead geometry and physics of Test 9,
but select the explicit Euler and second-order Heun integrators respectively.
Each case executes 1,000 steps and samples the trajectory every 200 steps.

The versioned files record a paired run built with NVFORTRAN 24.3 and run on
one NVIDIA A30 (`cc80`, CUDA 12.3, `nofma`). The measured integration-loop
times were:

| Case | NVFORTRAN CPU | OpenACC A30 | Speedup |
| --- | ---: | ---: | ---: |
| Test 10, Euler | `9.787760 s` | `0.741555 s` | `13.20x` |
| Test 11, RK2/Heun | `19.464682 s` | `1.254380 s` | `15.52x` |

Both CPU/GPU comparisons pass with `rtol=1e-6` and `atol=1e-9`. The worst
normalized differences were `4.00e-7` for Euler and `8.99e-7` for RK2.
These timings describe this development system and are not portable
performance guarantees.

Compare new outputs against the records and against each other with:

```sh
tests/performance/integrators/compare.sh euler CPU/statout.dat GPU/statout.dat
tests/performance/integrators/compare.sh rk2 CPU/statout.dat GPU/statout.dat
```

Set `RTOL` or `ATOL` in the environment to test a different tolerance.
