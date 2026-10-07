# Test Case 1: one-dimensional reference case

This is the basic one-dimensional benchmark. It uses the dimensionless
parameters $Q=12$, $V=2$ (printed as `Phi`), and $F_{ve}=12$ in the
definition of Reneker et al. (1727.5 with JETSPIN's own definition; the
`Fve` line of the output prints both), from established discrete-jet
reference calculations.

In the OpenACC build this one-dimensional case runs the CPU build's code on
the host and launches no kernel: one-dimensional systems do not take the
device step, and the jet keeps its single bead.

- [Input file](../../examples/input-1/input.dat)
- [LaTeX manual section](../../manual/test1.tex)
