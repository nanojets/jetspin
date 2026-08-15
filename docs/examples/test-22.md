# Test Case 22: repeated dynamic refinement and capacity growth

Test Case 22 extends the hybrid Test Case 21 path from one Akima event to
three. It retains Maxwell rheology, Yarin evaporation, stochastic Platen air
drag, nozzle insertion, and 81 active numerical anchor beads. The initial jet
contains 400 elements over the first 8 cm of a 16 cm domain.

The test runs 15,800 steps at `2.5e-8 s`. Developer-only environment overrides
set both the initial reserve and each refinement growth increment to 20. This
forces every accepted remesh to exceed the current capacity. Production runs
continue to use the normal 100-entry `incnpjet` policy.

The event-count contract is intentionally restricted to NVFORTRAN CPU and GPU
builds. A different compiler can perturb this stochastic bending trajectory
enough to accept a different number of marginal remeshes within the same final
time; Test Case 21 remains the portable single-event GFortran check.

For each event the checker requires:

- a denser and strictly ordered target mesh;
- an unchanged active-anchor count and unchanged anchor fields;
- separate conservation of reference and evaporated volume;
- conservation of mass and charge;
- a consistent Gaussian-history stride and positive retained cycle; and
- finite output with a final `nref` of three.

The validated NVFORTRAN 24.3 results are:

| Path | Event steps | Active elements | Capacity sequence | Final elements |
| --- | --- | --- | --- | ---: |
| CPU | 14,301; 14,944; 15,719 | 413→455; 455→498; 499→536 | 420→477→520→558 | 536 |
| Native A30 | 14,301; 14,944; 15,701 | 413→455; 455→498; 499→536 | 420→477→520→558 | 536 |
| A30 complete-force oracle | 14,301; 14,944; 15,719 | 413→455; 455→498; 499→536 | 420→477→520→558 | 536 |

All anchor-field differences are zero at printed precision. Reference volume,
evaporated volume, mass, and charge are conserved at approximately `1e-16`.
The native CPU/GPU comparison passes with `rtol=2.5e-2`; its largest relative
difference is about 2.38 percent in final `vz`. The force oracle reproduces
the CPU topology and passes with `rtol=6e-7`. The 18-step shift in the third
native event is therefore attributed to amplification of the different direct
Coulomb summation order rather than to the reallocation lifecycle.

The normal OpenACC path computes Akima coefficients and spline interpolation
on the GPU. An A30 transfer audit found four complete state downloads: the
three accepted refinement events and final shutdown. The event path performs
exactly three
Gaussian-history uploads, topology rebinds, and evaporation-state rebinds.
Unsuccessful threshold scans return only three reduction scalars, so ordinary
timesteps do not download the complete jet state.

Run:

```sh
tests/refinement/run_test22.sh nvfortran
tests/refinement/run_test22.sh openacc
tests/refinement/run_test22.sh force-oracle
```

- [Input file](../../examples/input-22/input.dat)
- [Input-file notes](../../examples/input-22/README.md)
