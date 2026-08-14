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
76 removals, five reallocations, and 1224 active beads. An exploratory A30
run completed in about 2.7 s with 499 additions, 898 removals, three
reallocations, and 401 active beads. The topology difference is expected for
this roundoff-sensitive high-resolution probe and is not currently a strict
regression failure.
