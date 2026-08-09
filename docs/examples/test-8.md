# Test Case 8: Yarin 2001 reference parameters

This reference case uses the 6 wt% aqueous-PEO material and ambient
parameters reported by Yarin, Koombhongse, and Reneker (2001). It exercises
the default vapour-diffusivity correlation, the evaporation cutoff, and
concentration-dependent viscosity and relaxation time.

Important regression targets include:

$$
\frac{V}{V_0}=\frac{0.06}{0.9}=0.0666666667,
\qquad c_p=0.9.
$$

At the cutoff, the reference material-property ratios are approximately
$\mu/\mu_0=43.978335$, $\theta/\theta_0=15$, and $G/G_0=2.931889$.

- [Input file](../../examples/input-8/input.dat)
- [Case notes](../../examples/input-8/README.md)
- [LaTeX manual section](../../manual/test8.tex)
