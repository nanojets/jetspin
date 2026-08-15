# Test Case 22 input

Test Case 22 is the repeated-refinement stress counterpart of Test Case 21.
It retains Maxwell rheology, Yarin evaporation, stochastic Platen air drag,
nozzle insertion, a pre-extended 400-element jet, and 0.10 cm numerical anchor
spacing. The timestep is reduced to `2.5e-8 s`, and the 15,800-step run accepts
three Akima remeshing events.

The input itself does not change production allocation policy. The dedicated
runner supplies developer-only 20-entry initial and growth reserves so every
accepted event exceeds the current capacity. Normal JETSPIN runs continue to
reserve and grow by 100 entries.

The three-event contract is a same-compiler NVFORTRAN CPU/GPU comparison.
Use Test Case 21 for the portable single-event GFortran refinement check.

Run the complete validation with:

```sh
tests/refinement/run_test22.sh nvfortran
tests/refinement/run_test22.sh openacc
tests/refinement/run_test22.sh force-oracle
```

In the OpenACC build, Akima coefficient construction and interpolation run on
the GPU. Target-mesh preparation, conservation, host reallocation, and the
persistent-data rebind remain rare event-level host work.
