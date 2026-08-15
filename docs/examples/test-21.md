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
reduces only three values to the host: the last over-threshold segment index,
total path length, and the nozzle correction used to estimate the target mesh
size. When the complete historical acceptance criterion can be satisfied, the
active state is downloaded once. The host prepares normalized source and
target coordinates; Akima slopes, tangents, cubic coefficients, and all 11
field interpolations then execute on the GPU. Conservation, anchor checks, and
the final state assembly remain on the host. The remeshed state is uploaded
once and persistent device integration resumes.

The initial 400-element allocation reserves 100 additional slots. The current
single refinement therefore changes active bounds without reallocating a
mapped array. A separate `capacity-growth` mode reduces only this developer
reserve to 50 entries. The accepted mesh then exceeds the original capacity:
JETSPIN releases the old persistent mappings before host reallocation,
resizes the Platen and Coulomb workspaces, repacks the Gaussian-history stride,
and binds the new arrays once. The Akima GPU kernels are unchanged in both
modes.

Current event-level results are:

| Build | Event step | Active elements across event | Final elements |
| --- | ---: | ---: | ---: |
| GFortran CPU | 7,134 | 413 to 454 | 455 |
| NVFORTRAN 24.3 CPU | 7,118 | 412 to 460 | 462 |
| NVFORTRAN 24.3 native A30 | 7,122 | 412 to 459 | 461 |
| NVFORTRAN 24.3 A30 with complete-force oracle | 7,118 | 412 to 460 | 462 |

All paths start with 79 interior anchors and end with `nref=1`. Existing anchor
position, velocity, stress, reference radius, and evaporation radius are
preserved at the interpolation knots to the checker's `1e-12` tolerance.
Reference and evaporated-volume conservation errors remain at order `1e-16`,
and no capacity reallocation occurs.

The capacity-growth results are:

| Build | Event step | Active elements across event | Final elements | Capacity |
| --- | ---: | ---: | ---: | ---: |
| NVFORTRAN 24.3 CPU | 7,156 | 413 to 455 | 456 | 450 to 557 |
| NVFORTRAN 24.3 native A30 | 7,155 | 413 to 455 | 456 | 450 to 557 |
| NVFORTRAN 24.3 A30 with complete-force oracle | 7,156 | 413 to 455 | 456 | 450 to 557 |

CPU and force oracle have identical topology; their printed values agree
within `3e-7` relatively. The native GPU also has identical event and final
topology and differs from the CPU statistics by at most `0.8%`. This tighter
agreement than the standard case is incidental to the changed Gaussian
storage stride and must not be interpreted as a new global tolerance.

The aerodynamic Gaussian history is generated for all 500 reserved indices
before loop timing and copied to the device once. It contains 24,048,000
doubles (192,384,000 bytes) for this 8,000-step run. No random values or noise
arrays cross the host/device boundary inside the timestep loop. The independent
Langevin `noise yes` term is intentionally not enabled.

In capacity-growth mode the history grows from 450 to 557 indexed entries.
All values for existing indices are retained and new indices are generated
once, producing 26,784,000 doubles (214,272,000 bytes). No Gaussian draw is
performed in the timestep loop.

The native GPU and CPU trajectories are not binary-identical because the
target-centric direct-Coulomb kernel changes the accumulation order and the
bending instability amplifies the resulting floating-point perturbation. The
complete-force oracle reproduces the NVFORTRAN CPU event and topology exactly;
their final written statistics differ by at most `3.6e-6` relatively. The
native A30 path is therefore accepted through event topology, conservation,
anchor invariants, finite output, and transfer behaviour.

An `NVCOMPILER_ACC_NOTIFY=2` audit of the native A30 run observed 6,123
eligible threshold scans, each exchanging only the three reduction scalars.
The full active state was downloaded at the accepted refinement event and at
final shutdown, and the remeshed state was uploaded once. Within that rare
event, source and target coordinates are uploaded once, each of the 11 source
fields is uploaded once, and each interpolated field is downloaded once. There
is no full-state transfer on an ordinary timestep. Repeating the audit in
capacity-growth mode again found exactly those two full downloads, one
topology/evaporation rebind, one resized-history upload, and one Platen/Coulomb
workspace recreation at the accepted event.

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
