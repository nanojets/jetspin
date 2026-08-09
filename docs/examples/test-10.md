# Test Case 10: Euler CPU/GPU benchmark

Test Case 10 is the Euler-integrator counterpart of the fixed 1,000-bead Test
9 benchmark. It retains the same geometry, material parameters, constant
axial electric field, and 1,000 steps, while selecting integrator `1`.

The persistent OpenACC path keeps the primary state, force arrays, and Euler
derivatives on the device. Its final state update is fused with the path-length
and maximum-stress reductions. Scheduled scalar output transfers only the
selected bead and reduced quantities; full state is synchronized for output
formats that require it and for restart data.

An initial NVFORTRAN 24.3 comparison on one NVIDIA A30 measured `9.787760 s`
on CPU and `0.741555 s` with OpenACC, a `13.20x` speedup. The six saved samples
pass the CPU/GPU comparison with `rtol=1e-6` and `atol=1e-9`; the worst
normalized difference is `4.00e-7`.

- [Input file](../../examples/input-10/input.dat)
- [Versioned numerical records](../../tests/performance/integrators/README.md)
