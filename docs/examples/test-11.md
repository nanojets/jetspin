# Test Case 11: RK2 CPU/GPU benchmark

Test Case 11 is the second-order Heun counterpart of the fixed 1,000-bead Test
9 benchmark. It retains the same geometry, material parameters, constant
axial electric field, and 1,000 steps, while selecting integrator `2`.

Since 2026-10-06 the OpenACC build runs the common device step
(`device_rk_step` in `source/device_step_mod.f90`) from the first step, with
14 kernel launches per step since 2026-10-07, one of which freezes the beads
that reach the collector (13 on 2026-10-06); it keeps both equation-of-motion stages, the
intermediate state, and the primary state on the device. The final Heun update
is fused with the path-length and maximum-stress reductions. Scheduled scalar
output transfers only the selected bead and reduced quantities; full state is
synchronized only when an output format or restart requires it.

An initial NVFORTRAN 24.3 comparison on one NVIDIA A30 measured `19.464682 s`
on CPU and `1.254380 s` with OpenACC, a `15.52x` speedup. The six saved samples
pass the CPU/GPU comparison with `rtol=1e-6` and `atol=1e-9`; the worst
normalized difference is `8.99e-7`. Since 2026-10-07 the leading bead, which
starts on the collector, is frozen there (see [Test 9](test-9.md)): the A30
rows equal the CPU's, and the versioned records were regenerated. With the
Coulomb kernel of 2026-10-05 (one gang per target, see [Test 9](test-9.md))
the A30 loop takes `0.223 s` against `1.078 s` with the previous kernel, the
process bound to the GPU's NUMA node.

The common development targets `nvfortran-openacc-force-oracle` and
`nvfortran-openacc-coulomb-oracle` also cover both non-evaporative RK2 stages.
Against the same CPU record their worst normalized differences were
`6.00e-10` and `1.00e-7`, respectively. They are numerical diagnostics, not
performance builds.

- [Input file](../../examples/input-11/input.dat)
- [Versioned numerical records](../../tests/performance/integrators/README.md)
