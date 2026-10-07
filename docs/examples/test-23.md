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
- a single Gaussian pool allocation that capacity growth does not rebuild; and
- finite output with a final `nref` of three.

With NVFORTRAN 24.3, the validated results are:

| Path | Refinement steps | First removal | Removals | Final elements |
| --- | --- | ---: | ---: | ---: |
| CPU | 14,268; 14,971; 15,691 | 15,656 | 9 | 526 |
| Native A30 | 14,268; 14,971; 15,691 | 15,656 | 9 | 526 |
| A30 complete-force oracle | 14,268; 14,971; 15,691 | 15,656 | 9 | 526 |

NVFORTRAN 25.5 reproduces the CPU row and its removal schedule exactly.
All three paths preserve 81 active anchors at each accepted event and grow
capacity from 420 to 481, 517, and 556. The native path and the force oracle
reproduce the CPU events and removal schedule step for step; all 81
statistics rows agree with the CPU within `8.0e-10` (native) and `6.1e-8`
(oracle) relatively, including the three rows written after the leading bead
has reached the collector. The Coulomb-only oracle, which evaluates only the
direct sum on the host, writes a statistics file byte-identical to the native
run's.

Two of the persistent-path defects fixed on 2026-09-30 (see
[OpenACC](../introduction/openacc.md)) showed here. The device kept the full
charge on the bead being inserted at the nozzle, so the native events moved
(14,256; 14,976; 15,780, with ten removals); and the Platen kernels let a bead
frozen at the collector keep moving, so the oracles, which keep those kernels
on the device, differed from the CPU by up to 1.4 percent once the leading
bead had reached the collector. The device kernels still evaluate in a
different floating-point order; the event-level topology and conservation
contracts therefore remain the acceptance criteria rather than binary
trajectory identity.

The GPU removal primitive updates the active lower bound and clears only the
removed bead. It does not download the full jet. Full-state communication is
reserved for accepted refinement events and final shutdown. At each accepted
event the host walks the mass boundaries and constructs the normalized target
mesh; the GPU computes Akima coefficients and all 11 field interpolations,
reconstructs bead and evaporated volumes, applies the reference-volume
conservation rescale, and converts mass and charge densities back to
per-bead quantities; the host keeps the anchor bookkeeping and final
assembly, the conservative radius-floor and evaporation-limit corrections,
and the invariant checks before the persistent state is rebound.

An A30 transfer audit records four complete state downloads: the three
accepted Akima events and final shutdown. It also records exactly three
topology rebinds and three evaporation-state rebinds; the Gaussian pool is
uploaded once at startup.
Each timestep returns one 20-byte topology record carrying the insertion and
removal decisions (until 2026-10-01, seven separate scalar transfers); an
accepted removal returns only the collected point fields and clears that one
device slot. With the scan result, an ordinary step moves two small records
in total, against thirteen transfers before.

The Akima A/B build evaluates the historical host routine and the new device
routine at every event. Across 33 field/event comparisons, the largest
coefficient relative error is `1.03e-14`, the largest interpolated absolute
error is `1.82e-12`, and the largest interpolated relative error is
`1.89e-15`. Mass and charge density are identical at printed precision. An
independent normal-device run and host-Akima-oracle run produce byte-identical
`statout.dat` files and the same event/removal topology.

The `refinement-compare` backend (target
`nvfortran-openacc-compare-refinement`) evaluates the host and device
reconstruction of volume, mass, charge, and evaporated volume at every event;
all three events agree at roundoff (worst absolute `8.9e-16`, worst relative
`1.4e-16`).

Run:

```sh
tests/refinement/run_test23.sh nvfortran
tests/refinement/run_test23.sh openacc
tests/refinement/run_test23.sh force-oracle
tests/refinement/run_test23.sh host-akima
tests/refinement/run_test23.sh akima-compare
tests/refinement/run_test23.sh refinement-compare
```

- [Input file](../../examples/input-23/input.dat)
- [Input-file notes](../../examples/input-23/README.md)
