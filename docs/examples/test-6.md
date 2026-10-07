# Test Case 6: Kelvin–Voigt fluid

This three-dimensional simulation activates the Kelvin–Voigt constitutive
model with `kvfluid yes`. It provides the reference case for the
quasi-Newtonian/Kelvin–Voigt integration path.

Coulomb forces use the multiple-step algorithm (`multiple step every 100`).
Over 1,000,000 steps it agrees with the direct sum within 9e-8 (axial
position) and 4.4e-6 (off-axis distance), with the same 146 insertions and
35 removals (one a step later); two MPI ranks agree with the serial run
within 6e-6 (2026-10-07). Before the multiple-step corrections of
2026-10-07 the run departed from the direct sum by up to 25 % in the
off-axis distance, with one more insertion, and, once the jet filled its
arrays, the per-rank arrays one element short made a GFortran build stop
near step 668,500 (out-of-bounds write; NVFORTRAN overwrote memory
silently). The OpenACC build runs this case with the CPU code.

- [Input file](../../examples/input-6/input.dat)
- [LaTeX manual section](../../manual/test6.tex)
