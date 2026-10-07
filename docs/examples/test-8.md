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

## OpenACC build

Since 2026-10-06 this case (RK4, evaporation, no air drag) takes the common
device step (`device_rk_step` in `source/device_step_mod.f90`, see
[OpenACC](../introduction/openacc.md)). The OpenACC build runs the CPU
build's code, byte for byte, until the jet arrays reach index 100
(`npjet>=100`), and engages at step 13,699, with 89 active beads. Over the
full 1e6 steps (2026-10-07, up to about 7,000 beads) the A30 takes 521 s
against 1957 s on one CPU core. Insertions fall at the CPU's steps up to
step 26,122 and removals up to step 34,783 (before the insertion fix of
2026-10-07 the events agreed only up to step 13,839, see
[Test 15](test-15.md)); then the two runs separate chaotically, with events
a few tens of steps apart either way, and both make 7,090 insertions and 130
removals and end with the same bead count. The restart check
(`tests/restart/run.sh`) restarts this case at step 14,000, after an
insertion has compacted the arrays below 100 beads, and the restarted run
resumes on the device as the uninterrupted one does.

NVFORTRAN 25.5 and 24.3 round divisions differently: over the first 20,000
steps their runs separate after the 128th topology event, with the same
totals. The reference values are those of NVFORTRAN 24.3.
