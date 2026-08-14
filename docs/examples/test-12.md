# Test Case 12: Platen stochastic CPU/GPU benchmark

Test Case 12 combines the fixed 1,000-bead geometry and constant axial field
of Test 9 with the Platen integrator and stochastic air-drag parameters of Test
4. Insertion, removal, evaporation, dynamic refinement, and multiple-step
Coulomb summation remain disabled. It executes 1,000 steps and writes six
statistical samples.

Before loop timing, rank 0 generates the complete Gaussian history in a fixed
`step, bead, component, draw` order. The history contains 6,006,000 double
precision values (48,048,000 bytes). CPU integration reads this same layout;
the OpenACC build copies it to the device once during initialization. No noise
generation or noise transfer occurs inside the temporal loop.
The implementation caps a pre-generated history at 100,000,000 values; longer
fixed runs reuse it cyclically. Test 12 is below this limit and never wraps.

The persistent GPU path performs the three Platen force evaluations,
positive/negative predictors, stochastic velocity update, Heun position and
stress updates, and statistics on the device. Host synchronization follows
the same output and restart rules as Tests 9--11.

Both development oracle targets cover the same Platen sequence. The complete
force oracle evaluates every trusted host force or partial stress evaluation;
the Coulomb-only oracle replaces only the direct sum. With either target, all
six statistical rows and fourteen columns match the NVFORTRAN CPU record
exactly. These oracle builds deliberately transfer stage data and are not
performance configurations.

An initial NVFORTRAN 24.3 comparison on one NVIDIA A30 measured `29.050473 s`
on CPU and `1.874810 s` with OpenACC, a `15.50x` speedup. CPU and GPU outputs
were identical in all six rows and fourteen columns at the written precision.

- [Input file](../../examples/input-12/input.dat)
- [Versioned numerical records](../../tests/performance/integrators/README.md)
