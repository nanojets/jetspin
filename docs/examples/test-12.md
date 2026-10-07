# Test Case 12: Platen stochastic CPU/GPU benchmark

Test Case 12 combines the fixed 1,000-bead geometry and constant axial field
of Test 9 with the Platen integrator and stochastic air-drag parameters of Test
4. Insertion, removal, evaporation, dynamic refinement, and multiple-step
Coulomb summation remain disabled. It executes 1,000 steps and writes six
statistical samples.

Before loop timing, rank 0 generates the Gaussian pool, each step's slice in
the draw order of the per-step block (bead, component, draw). The pool
contains 6,006,000 double precision values (48,048,000 bytes). CPU
integration reads the same values; the OpenACC build copies them to the device
once during initialization. No noise generation or noise transfer occurs
inside the temporal loop. A run that consumes the whole pool (`noise pool`,
default 100,000,000 values) reuses it cyclically; Test 12 is below this limit
and never wraps. The outputs were bit-identical to the records made before
the pool became the default layout
([random numbers](../introduction/random-numbers.md)), until the collector
rule of 2026-10-07 held the leading bead at 16 cm (records regenerated).

The persistent GPU path performs the three Platen force evaluations,
positive/negative predictors, stochastic velocity update, Heun position and
stress updates, and statistics on the device. Host synchronization follows
the same output and restart rules as Tests 9--11.

Since 2026-10-06 the step is the common Platen device step
(`device_platen_step` in `source/device_step_mod.f90`) with the fused
two-kernel tail of the dynamic runs: 12 kernel launches per step instead of
18, and the same output, byte for byte. Since 2026-10-07 a thirteenth kernel,
`accelerator_freeze_at_collector`, holds the leading bead on the collector
(see [Test 9](test-9.md)); CPU and A30 rows remain identical.

Both development oracle targets cover the same Platen sequence. The complete
force oracle evaluates every trusted host force or partial stress evaluation;
the Coulomb-only oracle replaces only the direct sum. With either target, all
six statistical rows and fourteen columns match the NVFORTRAN CPU record
exactly. These oracle builds deliberately transfer stage data and are not
performance configurations.

An initial NVFORTRAN 24.3 comparison on one NVIDIA A30 measured `29.050473 s`
on CPU and `1.874810 s` with OpenACC, a `15.50x` speedup. CPU and GPU outputs
were identical in all six rows and fourteen columns at the written precision.
With the Coulomb kernel of 2026-10-05 (one gang per target, see
[Test 9](test-9.md)) and before the fused tail of 2026-10-06, the A30 loop
took `0.354 s` against `1.659 s` with the previous kernel, the process bound
to the GPU's NUMA node.

- [Input file](../../examples/input-12/input.dat)
- [Versioned numerical records](../../tests/performance/integrators/README.md)
