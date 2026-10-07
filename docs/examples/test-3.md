# Test Case 3: three-dimensional PVP electrospinning

This case models a polyvinylpyrrolidone solution in a conventional
three-dimensional electrospinning setup. It uses a 16 cm collector distance,
a 0.005 cm jet radius at the nozzle, a rotating nozzle perturbation, bead
insertion, and the Maxwell rheological model.

The material and operating parameters are model parameters taken from the
input file, which describes a PVP solution (1300 kDa, ethanol:water 17:3 v:v,
about 2.5 wt%) under operating conditions typical of PVP electrospinning
(9 kV, collector at 16 cm). The charge density, 4.4e4 statC/cm³ (1.47e-2 C/L), and the
elastic modulus, 5e4 g/(cm s²) (5000 Pa), are those of the input file. The
example is not compared with an experiment.

- [Input file](../../examples/input-3/input.dat)
- [LaTeX manual section](../../manual/test3.tex)

## OpenACC build

The OpenACC build runs the CPU build's code, byte for byte, until the jet
arrays reach index 100 (`npjet>=100`), and then the common device step
(RK2, Maxwell, no air drag; `device_rk_step` in
`source/device_step_mod.f90`, see [OpenACC](../introduction/openacc.md)).
With seed 317 it engages at step 1,114,161. Over 1e7 steps (2026-10-07, A30
against the NVFORTRAN CPU build) the axial position of the leading bead `x`,
the nozzle velocity `vn`, the bead count `n` and the currents are identical
in every printed row, the worst off-axis distance `yz` difference is 6.7e-5,
and every insertion and removal falls at the CPU's step up to step 2,364,190
(1,125,560 before the insertion fix of 2026-10-07, see
[Test 15](test-15.md)); both runs make 878 insertions and 767 removals.
At about 110 beads the A30 is slower than one CPU core: 2179 s against
about 1800 s for the 1e7 steps, because the RK step is still synchronous
and not fused (milestone M4 of `docs/STATE.md`).
