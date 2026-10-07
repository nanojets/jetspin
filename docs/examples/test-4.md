# Test Case 4: gas counterflow

This three-dimensional case enables aerodynamic drag (`system 4` and
`airdrag yes`) with a gas velocity opposed to the main jet direction. It
uses the stochastic Platen integrator required by the aerodynamic model.

- [Input file](../../examples/input-4/input.dat)
- [LaTeX manual section](../../manual/test4.tex)

Since 2026-10-06 this case reads its stochastic forcing from the
pre-generated Gaussian pool, in the CPU build as in the OpenACC one, like the
other Platen runs of the device step's models
([random numbers](../introduction/random-numbers.md)); before, it drew a
Gaussian block at every step. The noise realization therefore differs from
earlier records, and the regression baselines of case 4 were regenerated.
The OpenACC build runs the CPU code, byte for byte, until the jet has 100
beads, and the device step above that
([OpenACC](../introduction/openacc.md)): the jet reaches 100 beads after
about 1.12 million steps (0.011 s) and the collector at about 0.015 s, and
then carries about 115 beads.

Seed ensembles over the first 5 million steps (2026-10-07, means over
0.03-0.05 s) show that the change of noise layout leaves the physics
unchanged: 8 CPU seeds with the pool and 8 with the former step-by-step
noise agree within one standard error or so for every observable (off-axis
distance 2.959 cm, 114.6 beads, 21.0 degrees, path length 177.5 cm), and
24 A30 seeds agree with 24 CPU seeds in the same way. The noise matters
little here: the off-axis distance varies by about 8e-4 cm between seeds.
