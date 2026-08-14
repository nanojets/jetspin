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
76 removals, five reallocations, and 1224 active beads. A fresh NVFORTRAN
24.3/A30 standard run completed normally with 498 additions, 77 removals,
five reallocations, and 1221 active beads. Both the complete-force oracle and
the Coulomb-only oracle produced the same aggregate GPU totals. The remaining
three-bead final-count difference is consistent with threshold-sensitive
topology after roundoff amplification; this high-resolution case remains a
porting probe rather than a strict pointwise regression.
