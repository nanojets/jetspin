# Test Case 16: dynamic topology with evaporation

Test Case 16 is the Maxwell reference for the GPU dynamic-topology path
with the Yarin evaporation model. It uses the Maxwell RK4 integrator, starts
with 100 beads, enables both nozzle insertion and collector removal, and keeps
dynamic Akima refinement disabled. The small initial capacity deliberately
forces host reallocations.

The complete Maxwell RK4 chain is device-resident: all four force stages,
three intermediate state constructions, the final weighted update, direct
Coulomb summation, evaporation, charge smoothing/restoration, and statistics
run on the GPU. The normal build performs no host/device array transfer between
RK stages.

The NVIDIA A30 run reproduces the CPU topology totals: 111 additions, 122
removals, two reallocations, and 89 active beads after 1,000 steps. It is also
bit-for-bit identical to the preceding saved GPU trajectory, which guards the
porting increments against algorithmic regressions. A CPU/GPU pointwise
trajectory comparison is intentionally not the acceptance criterion: the
different direct-Coulomb accumulation order creates roundoff-level changes
that are amplified by the physical bending instability and by topology
threshold crossings.

A transfer audit with all development macros disabled shows only topology
decision scalars on ordinary timesteps. Bead records are transferred on actual
insertion/removal and output events; complete active arrays are synchronized
only for the two capacity reallocations and the final checkpoint. The
`nvfortran-openacc-host-forces` target remains available solely to isolate
forces during development and must not be used for performance measurements.

- [Input file](../../examples/input-16/input.dat)
- [Input-file notes](../../examples/input-16/README.md)
