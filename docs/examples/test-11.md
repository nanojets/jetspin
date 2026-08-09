# Test Case 11: RK2 CPU/GPU benchmark

Test Case 11 is the second-order Heun counterpart of the fixed 1,000-bead Test
9 benchmark. It retains the same geometry, material parameters, constant
axial electric field, and 1,000 steps, while selecting integrator `2`.

The persistent OpenACC path keeps both equation-of-motion stages, the
intermediate state, and the primary state on the device. The final Heun update
is fused with the path-length and maximum-stress reductions. Scheduled scalar
output transfers only the selected bead and reduced quantities; full state is
synchronized only when an output format or restart requires it.

An initial NVFORTRAN 24.3 comparison on one NVIDIA A30 measured `19.464682 s`
on CPU and `1.254380 s` with OpenACC, a `15.52x` speedup. The six saved samples
pass the CPU/GPU comparison with `rtol=1e-6` and `atol=1e-9`; the worst
normalized difference is `8.99e-7`.

- [Input file](../../examples/input-11/input.dat)
- [Versioned numerical records](../../tests/performance/integrators/README.md)
