# Test Case 19: high-resolution dynamic Kelvin–Voigt evaporation

Test Case 19 is the Kelvin–Voigt counterpart of Test 18. It uses the
Kelvin–Voigt RK4 integrator with evaporation, nozzle insertion, and collector
removal enabled. The initial jet length is 16 cm and `points 800`, giving a
nominal spacing of 0.02 cm (200 micrometres) between adjacent nodes. Dynamic
refinement is disabled and the case is serial.

This case is a high-resolution GPU porting and performance probe. Because the
direct Coulomb force is sensitive to the bead spacing and the bending dynamics,
CPU/GPU trajectories should be compared with a tolerance rather than as
bitwise-identical states.

- [Input file](../../examples/input-19/input.dat)
- [Test 18 Maxwell counterpart](test-18.md)

The NVFORTRAN CPU reference completed in about 14.5 s with 500 additions,
76 removals, five reallocations, and 1224 active beads. Since 2026-10-07 the
NVFORTRAN 24.3/A30 standard run gives the same events, step for step, and
statistics rows within `1e-5` relatively. Before, it gave 498 additions, 77
removals and 1221 active beads (both oracles the same): the device path
inserted one step later whenever the arrays had to grow, which was taken for
threshold-sensitive topology after roundoff amplification. NVFORTRAN 25.5,
CPU and A30 alike, gives 500 additions, 61 removals and 1239 active beads:
at step 62 a bead stops 8e-6 of a length unit (relative 1.8e-8) short of the
collector instead of reaching it, and the removals that follow differ. The
two compilers round divisions differently (`docs/introduction/compiling.md`);
the reference values are those of NVFORTRAN 24.3. This
high-resolution case remains a porting probe rather than a strict pointwise
regression.
