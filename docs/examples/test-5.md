# Test Case 5: dynamic refinement

This three-dimensional case activates adaptive remeshing with
`dynamic refinement yes`. The refinement interval and segment-length
threshold are specified in the input file.

Coulomb forces use the multiple-step algorithm (`multiple step every 100`).
Over 1,000,000 steps it agrees with the direct sum (`multiple step no`)
within 6e-8 in every printed row, with the same 88 insertions; two MPI ranks
agree with the serial run within 2e-8 (NVFORTRAN 24.3 and GFortran,
2026-10-07, after the multiple-step corrections listed in `CHANGELOG.md`).
The OpenACC build runs this case with the CPU code.

- [Input file](../../examples/input-5/input.dat)
- [LaTeX manual section](../../manual/test5.tex)
