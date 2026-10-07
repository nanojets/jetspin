# Test Case 20: stochastic Platen Maxwell evaporation

Test Case 20 provides a fixed-topology CPU/GPU baseline for the Maxwell
evaporative equations with the stochastic Platen integrator (`system 4`,
`integrator 4`). It starts with 1,000 beads over 16 cm and runs 100 steps.

The case combines the stochastic air-drag settings from Test 4 with the
evaporation settings used by Test 16. Insertion, removal, and dynamic
refinement are disabled so that differences isolate the three Platen EOM
evaluations, evaporation-dependent mass scaling, and Gaussian-noise indexing.

Run it from the executable directory with:

```sh
cp examples/input-20/input.dat execute/input.dat
(cd execute && ./main.x)
```

For the MPI baseline, use two ranks and compare the resulting `statout.dat`
against the serial output with the numerical comparison tool.

The OpenACC path pre-generates a pool of 600,600 Gaussian doubles on the CPU,
each step's slice in the draw order of the per-step block, and copies it to the
device once before loop timing; CPU and GPU outputs were bit-identical to the
records made before this pool became the default layout
([random numbers](../introduction/random-numbers.md)), until the collector
rule of 2026-10-07 held the leading bead at 16 cm. Its three Platen drift evaluations,
positive/negative predictors, stochastic velocity update, Heun position,
volume and stress updates, and statistics then execute inside persistent data
regions. No noise or stage array is transferred during an ordinary timestep.

Since 2026-10-06 the step is the common Platen device step
(`device_platen_step` in `source/device_step_mod.f90`) with the fused
two-kernel tail of the dynamic runs: 12 kernel launches per step instead of
19, and the same output, byte for byte. Since 2026-10-07 a thirteenth kernel,
`accelerator_freeze_at_collector`, holds the leading bead on the collector
(see [Test 9](test-9.md)); CPU and A30 rows remain identical.

An NVFORTRAN 24.3 comparison on one NVIDIA A30, with a build of 2026-08-14,
measured `3.177595 s` on the CPU and `0.224058 s` with OpenACC for the
integration loop, a `14.18x` speedup. The timing predates the evaporative
Coulomb kernel rewrite (2026-09-30) and the common device step, and has not
been measured again; Test 20 has no versioned record. The standard GPU output and both development oracle outputs—the
complete host-force oracle and the host Coulomb-only oracle—match all six CPU
rows and fourteen columns exactly at the written precision.
