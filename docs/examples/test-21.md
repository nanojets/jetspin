# Test Case 21: anchored dynamic refinement with Platen evaporation

Test Case 21 validates the hybrid OpenACC port of dynamic Akima refinement. It
combines:

- Maxwell rheology (`system 4`);
- stochastic Platen integration (`integrator 4`);
- the Yarin evaporation model;
- aerodynamic drag and its Gaussian forcing;
- nozzle insertion; and
- dynamic refinement with 0.10 cm anchor spacing and target length.

The case starts from 400 elements at 0.02 cm resolution, distributed from the
nozzle to 8 cm in a 16 cm nozzle-to-collector domain. This half-span
initialization leaves room for field-driven stretching before collection.
Interior anchors are initialized throughout the pre-extended jet every
0.10 cm. They are ordinary material beads used as fixed interpolation knots
during a remeshing event; no nanoparticle mass or other variable-mass feature
is enabled.

The reduced charge density limits sensitivity of the direct Coulomb sum, while
the increased external potential makes the short validation run reach the
refinement threshold. These are deliberate test accelerators rather than a
reference physical configuration. `dynamic refinement every 5.d-5` prevents
an immediate second remesh and leaves a post-refinement trajectory for
validation.

The port keeps Maxwell/Platen integration, Gaussian-history access, insertion,
and the repeated refinement-threshold scan on the GPU. Each eligible scan
reduces only three values and returns them to the host in one transfer: the
last over-threshold segment index, total path length, and the nozzle
correction used to estimate the target mesh size. When the complete historical acceptance criterion can be satisfied, the
active state is downloaded once. The host walks the mass boundaries and
prepares normalized source and target coordinates; Akima slopes, tangents,
cubic coefficients, and all 11 field interpolations then execute on the GPU,
which also reconstructs bead and evaporated volumes, applies the
reference-volume conservation rescale, and converts mass and charge densities
back to per-bead quantities. The anchor bookkeeping and final state assembly,
the conservative radius-floor and evaporation-limit corrections
(`enforce_radius_floor_conservative`, `enforce_evlim_conservative`), and the
invariant checks remain on the host. The remeshed state is uploaded once and
persistent device integration resumes.

The initial 400-element allocation reserves 100 additional slots. The current
single refinement therefore changes active bounds without reallocating a
mapped array. A separate `capacity-growth` mode reduces only this developer
reserve to 50 entries. The accepted mesh then exceeds the original capacity:
JETSPIN releases the old persistent mappings before host reallocation,
resizes the Platen and Coulomb workspaces, and binds the new arrays once; the
Gaussian pool stays mapped unchanged. The Akima GPU kernels are unchanged in both
modes.

Current event-level results are:

| Build | Event step | Active elements across event | Final elements |
| --- | ---: | ---: | ---: |
| GFortran CPU | 7,136 | 413 to 452 | 453 |
| NVFORTRAN 24.3 CPU | 7,131 | 412 to 453 | 455 |
| NVFORTRAN 24.3 native A30 | 7,131 | 412 to 453 | 455 |
| NVFORTRAN 24.3 A30 with complete-force oracle | 7,131 | 412 to 453 | 455 |

All paths start with 79 interior anchors and end with `nref=1`. Existing anchor
position, velocity, stress, reference radius, and evaporation radius are
preserved at the interpolation knots to the checker's `1e-12` tolerance.
Reference and evaporated-volume conservation errors remain at order `1e-16`,
and no capacity reallocation occurs.

The capacity-growth results are:

| Build | Event step | Active elements across event | Final elements | Capacity |
| --- | ---: | ---: | ---: | ---: |
| NVFORTRAN 24.3 CPU | 7,131 | 412 to 453 | 455 | 450 to 555 |
| NVFORTRAN 24.3 native A30 | 7,131 | 412 to 453 | 455 | 450 to 555 |
| NVFORTRAN 24.3 A30 with complete-force oracle | 7,131 | 412 to 453 | 455 | 450 to 555 |

The capacity-growth `statout.dat` is byte-identical to the standard-mode one,
on the CPU and on the A30: the Gaussian pool does not depend on capacity, so
the reserve size no longer changes the noise.

The aerodynamic noise comes from the 100,000,000-value Gaussian pool
(800,000,000 bytes), generated once before the loop and copied to the device
once, in both modes. No random values or noise arrays cross the host/device
boundary inside the timestep loop, and no Gaussian draw is performed there.
The independent Langevin `noise yes` term is intentionally not enabled.

The native A30 run and the complete-force oracle both reproduce the
NVFORTRAN CPU event and topology exactly; their written statistics differ
from the CPU by at most `2.4e-10` and `1.4e-7` relatively. CPU and GPU read
the same Gaussian pool; what remains is the different floating-point
evaluation order of the device kernels. Before the persistent-path fixes of
2026-09-30 (see [OpenACC](../introduction/openacc.md)) the A30 kept the full
charge on the bead being inserted at the nozzle; its fourth insertion then
fell one step earlier, which offset the pool position, and the event moved
to step 7,125.

An `NVCOMPILER_ACC_NOTIFY=2` audit of the native A30 run observed 6,132
eligible threshold scans, each returning its three results in one transfer
(until 2026-10-01, three uploads and three downloads). The insertion and
removal decisions return as one 20-byte record per step (before, two uploads
and four downloads).
The full active state was downloaded at the accepted refinement event and at
final shutdown, and the remeshed state was uploaded once. Within that rare
event, source and target coordinates are uploaded once, each of the 11 source
fields is uploaded once, and each interpolated field is downloaded once. There
is no full-state transfer on an ordinary timestep. Repeating the audit in
capacity-growth mode again found exactly those two full downloads, one
topology/evaporation rebind, and one Platen/Coulomb workspace recreation at
the accepted event; the Gaussian pool is uploaded once at startup and never
transferred again.

Run the case and its event-level checker with:

```sh
tests/refinement/run.sh
tests/refinement/run.sh nvfortran
tests/refinement/run.sh openacc
tests/refinement/run.sh force-oracle
tests/refinement/run.sh nvfortran capacity-growth
tests/refinement/run.sh openacc capacity-growth
tests/refinement/run.sh force-oracle capacity-growth
```

To run it manually:

```sh
cp examples/input-21/input.dat execute/input.dat
(cd execute && ./main.x > run.log 2>&1)
python3 tests/refinement/check_test21.py \
  --run-log execute/run.log --statout execute/statout.dat
```

The NVIDIA variants require NVHPC. `GPUCC=80` is the default used for the
NVIDIA A30.

- [Input file](../../examples/input-21/input.dat)
- [Input-file notes](../../examples/input-21/README.md)
