# Test Case 2: one-dimensional bead insertion

This case keeps the material parameters, jet radius at the nozzle, collector
distance, and external potential of Test Case 1, but starts from an initial
segment of 0.1 cm (0.319 cm in Test Case 1), which sets the length scale:
$Q\simeq122.1$, $V\simeq6.38$ (printed as `Phi`), and $F_{ve}\simeq122.1$
in the definition of Reneker et al. (1727.5 with JETSPIN's own definition).
It enables injection of new beads with `inserting yes` and their removal at
the collector with `removing yes`, and uses the fourth-order Runge-Kutta
integrator (`integrator 3`) with a timestep of 1e-7 s up to a final time of
2 s.

In the OpenACC build this one-dimensional case runs the CPU build's code,
since one-dimensional systems do not take the device step; its Coulomb sum
is computed on the GPU, with copies at every call, only while the jet has at
least 128 active beads.

- [Input file](../../examples/input-2/input.dat)
- [LaTeX manual section](../../manual/test2.tex)
