# Test Case 23 input

Test Case 23 extends the repeated-refinement Test Case 22 workload with
collector removal. It uses Maxwell rheology, Yarin evaporation, stochastic
Platen air drag, nozzle insertion, a pre-extended 400-element jet, and
0.10 cm numerical anchor spacing.

The initial jet occupies the first 8 cm of a shorter 12 cm domain. The
external potential is reduced in the same proportion, so the physical
potential gradient is unchanged from Test Case 22. During the 16,000-step
run, the jet undergoes three accepted Akima remeshes and reaches the
collector. Beads are removed both before and after the third remesh.

The dedicated runner supplies developer-only 20-entry initial and growth
reserves. Every accepted event therefore exceeds the current capacity and
exercises release, host reallocation, Gaussian-history repacking, and device
rebind while collector removal changes the lower active bound. Production
runs continue to reserve and grow by 100 entries.

Run the complete validation with:

```sh
tests/refinement/run_test23.sh nvfortran
tests/refinement/run_test23.sh openacc
tests/refinement/run_test23.sh force-oracle
tests/refinement/run_test23.sh host-akima
tests/refinement/run_test23.sh akima-compare
```

The normal OpenACC path computes Akima coefficients and interpolation on the
GPU. `host-akima` retains the historical host spline as an oracle, while
`akima-compare` evaluates both implementations field by field at all three
events. Target-mesh preparation and conservation remain on the host.
