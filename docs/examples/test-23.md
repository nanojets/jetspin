# Test Case 23: refinement interleaved with collector removal

Test Case 23 extends the repeated-refinement Test Case 22 path through the
collector. It retains Maxwell rheology, Yarin evaporation, stochastic Platen
air drag, nozzle insertion, and 0.10 cm numerical anchor spacing. The initial
jet contains 400 elements over the first 8 cm of a 12 cm domain.

The applied potential is scaled from Test Case 22 to preserve the same
physical potential gradient. The test runs 16,000 steps at `2.5e-8 s`.
Developer-only environment overrides set both the initial reserve and each
refinement growth increment to 20 entries, forcing every accepted remesh to
replace the current allocation.

The validation requires three accepted Akima events, three capacity changes,
and at least one collector removal both before and after the third event. It
checks every remesh for:

- a denser, strictly ordered target mesh;
- unchanged active-anchor fields;
- separate conservation of reference and evaporated volume;
- conservation of mass and charge;
- a consistent repacked Gaussian-history stride; and
- finite output with a final `nref` of three.

With NVFORTRAN 24.3, the validated results are:

| Path | Refinement steps | First removal | Removals | Final elements |
| --- | --- | ---: | ---: | ---: |
| CPU | 14,301; 14,944; 15,718 | 15,647 | 10 | 527 |
| Native A30 | 14,301; 14,944; 15,700 | 15,647 | 4 | 532 |
| A30 complete-force oracle | 14,301; 14,944; 15,717 | 15,647 | 7 | 530 |

All three paths preserve 81 active anchors at each accepted event and grow
capacity from 420 to 477, then 520, and finally 557. The native trajectory is
not expected to reproduce the exact later removal schedule: direct Coulomb
accumulation uses a different floating-point summation order on the GPU, and
the stochastic bending trajectory amplifies that difference near the
collector. The event-level topology and conservation contracts therefore
validate this test rather than binary trajectory identity.

The GPU removal primitive updates the active lower bound and clears only the
removed bead. It does not download the full jet. Full-state communication is
reserved for accepted refinement events and final shutdown. At each accepted
event the host constructs the normalized target mesh, the GPU computes Akima
coefficients and all 11 field interpolations, and the host applies conservation
and anchor rules before the persistent state is rebound.

An A30 transfer audit records four complete state downloads: the three
accepted Akima events and final shutdown. It also records exactly three
Gaussian-history uploads, topology rebinds, and evaporation-state rebinds.
Each timestep returns one four-byte removal flag; an accepted removal returns
only the collected point fields and clears that one device slot.

The Akima A/B build evaluates the historical host routine and the new device
routine at every event. Across 33 field/event comparisons, the largest
coefficient relative error is `2.48e-15`, the largest interpolated absolute
error is `7.28e-12`, and the largest interpolated relative error is
`3.32e-15`. Mass and charge density are identical at printed precision. An
independent normal-device run and host-Akima-oracle run produce byte-identical
`statout.dat` files and the same event/removal topology.

Run:

```sh
tests/refinement/run_test23.sh nvfortran
tests/refinement/run_test23.sh openacc
tests/refinement/run_test23.sh force-oracle
tests/refinement/run_test23.sh host-akima
tests/refinement/run_test23.sh akima-compare
```

- [Input file](../../examples/input-23/input.dat)
- [Input-file notes](../../examples/input-23/README.md)
