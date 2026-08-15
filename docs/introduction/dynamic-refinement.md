# Dynamic refinement

JETSPIN's dynamic-refinement procedure prevents the Lagrangian jet
discretization from becoming excessively coarse as the filament stretches.
When an element exceeds a configured length threshold, part of the jet is
remeshed and its state is reconstructed on a denser parametrization.

The scientific formulation is described in the
[LaTeX manual](../../manual/refinement.tex). This page concentrates on the
runtime workflow and the invariants needed to maintain the implementation in
[`dynamic_refinement_mod.f90`](../../source/dynamic_refinement_mod.f90).

## Configuration

A minimal input fragment is:

```text
dynamic refinement yes
dynamic refinement every 1.d-3
dynamic refinement threshold 0.4d0
```

The supported controls are:

| Directive | Purpose |
| --- | --- |
| `dynamic refinement yes` | Enable refinement and bead tagging |
| `dynamic refinement no` | Disable refinement |
| `dynamic refinement every <time>` | Minimum interval between accepted refinements |
| `dynamic refinement threshold <length>` | Maximum target element length |
| `dynamic refinement start <time>` | Do not refine before this time |
| `dynamic refinement anchor <length>` | Spacing used to preserve anchor beads as interpolation knots |

Times and lengths use the normal input units: seconds and centimetres before
internal nondimensionalization. Both `every` and `threshold` are required
when refinement is enabled. The threshold must be at least twice the base
resolution; values below five times the resolution are raised to that safer
minimum. If internal defaults are used, the code derives them from the base
resolution and emits a warning.

[`examples/input-5/input.dat`](../../examples/input-5/input.dat) is the
historical refinement case. [Test Case 21](../examples/test-21.md) combines
pre-extended anchors, insertion, Maxwell evaporation, and stochastic Platen
integration in a focused remeshing validation.

## Runtime workflow

The integrator calls `driver_dynamic_refinement` before advancing the next
time step. The driver follows this sequence:

```text
check enabled state and start time
                 |
                 v
wait for the configured check interval
                 |
                 v
measure element lengths along the active jet
                 |
                 v
estimate whether remeshing adds useful resolution
                 |
                 v
construct an anchored target mesh
                 |
                 v
Akima-interpolate continuous bead fields
                 |
                 v
reconstruct and renormalize conserved quantities
                 |
                 v
update bounds, allocation flags, and dependent workspaces
```

When bead insertion is active, refinement is deferred until at least ten
active beads exist. A threshold crossing alone is also insufficient: the
estimated fitted mesh must contain more points than the current active mesh.
After the configured `every` interval has elapsed, the threshold is tested on
each subsequent timestep until a useful refinement is accepted. The interval
counter is reset only at that accepted event. Thus `every` is a minimum delay
between remeshes, not a permanently sparse sampling period. The counter
`irefinementdone`, available as the `nref` observable, increases only when a
refinement is accepted.

## Path parametrization and anchors

`compute_length_path` assigns each bead a normalized path coordinate `jetpt`.
The remesher uses this coordinate, rather than an individual Cartesian axis,
so the same algorithm works for one- and three-dimensional trajectories.

The refined region extends from the collected/nozzle side through
`irefbeadstart`, the last element found above the threshold. Beads marked in
`jetbd` divide it into anchored segments. Each segment receives a number of
new intervals equal to the ceiling of its arc length divided by the requested
`dynamic refinement threshold`. The base insertion `resolution` does not
define this target mesh. Anchor locations are copied into the target mesh
exactly; the unrefined tail is retained rather than interpolated.

Anchor beads reduce interpolation drift and preserve selected material
locations, but their array indices may still change when the active interval
is rebased. They are ordinary jet material unless the separate variable-mass
feature is explicitly enabled; an anchor flag alone does not represent a
nanoparticle or freeze a bead between timesteps.

For the historical one-bead startup, anchors continue to be assigned as new
material is inserted at the nozzle. A pre-extended initial jet now receives
interior anchors before its first timestep, at the configured anchor spacing.
The two path endpoints are not tagged because they already occur in both
meshes; treating an endpoint as an interior anchor would create a zero-length
anchored segment.

## Interpolated and reconstructed state

Akima interpolation is applied to the continuous fields on the new path
mesh, including:

- coordinates and velocities;
- viscoelastic stress;
- mass and charge densities;
- cross-sectional radius;
- evaporation cross-sectional radius when evaporation is enabled.

Stress, mass, charge, and radii are forced non-negative after interpolation.
Discrete flags such as bead anchors and breakup markers are mapped separately
rather than treated as continuous scalar data. Collected-bead flags are
rebuilt from the remeshed coordinates.

Before interpolation, mass and charge are converted to densities using the
reference material volume. After remeshing they are converted back to bead
quantities. This avoids interpreting an old per-bead amount as though it
belonged unchanged to a new bead volume.

## Volume conservation and evaporation

Reference volume `jetvl` is reconstructed from the refined element length
and interpolated cross section. The refined portion is then renormalized so
its total reference volume equals the pre-refinement total.

With evaporation, instantaneous post-evaporation volume `jetve` follows the
same procedure independently, using the evaporation radius `jetce`. Thus the
algorithm preserves the two totals separately:

- `jetvl` describes reference material volume;
- `jetve` describes current volume after solvent loss.

They must not be merged, substituted for one another, or normalized with a
shared factor. Their local distribution may change slightly through spline
interpolation even though each total is conserved.

## OpenACC path

The serial Maxwell/Platen configuration in Test Case 21 keeps integration,
insertion, and the jet state in a persistent OpenACC data region. Once the
minimum refinement interval has elapsed, a device reduction returns only the
last over-threshold segment index, total path length, and the nozzle-length
correction. These scalars reproduce the inexpensive part of the historical
acceptance test without downloading bead arrays.

Only when that test predicts a denser target mesh does JETSPIN download the
active topology and evaporation state. The host computes normalized path
coordinates and target knots. Akima then follows four accelerator phases:

1. one independent secant slope per source interval;
2. serial extrapolation of the four endpoint slopes;
3. one independent Akima tangent and one cubic coefficient set per knot; and
4. one independent interval search and interpolation per target knot.

The endpoint extrapolation is constant work; the three knot-oriented phases
are OpenACC parallel loops. They contain no reduction or order-dependent sum.
Source and target coordinates are uploaded once per accepted event. Each of
the 11 source fields is then uploaded once and its interpolated result is
downloaded once.

The interpolated cross-section radius then feeds a second device kernel that
reconstructs bead volume/evaporation volume, rescales each for
reference-volume conservation, and converts the interpolated mass/charge
densities back to per-bead quantities -- the same volume-conservation and
mass/charge-conversion rules described above, executed on the device instead
of the host. Only the data-dependent bookkeeping that builds the normalized
target mesh in the first place (the mass-boundary walk, `jetptc`/`jetbdc`
assembly) and the anchor save/restore/final-assembly steps remain host-side:
they are a small, rare, inherently sequential scan, not a per-bead parallel
operation. After those host checks, the remeshed state is uploaded and device
integration resumes.

The normal Test 21 allocation reserves one `incnpjet` block and does not need
to grow. The capacity-growth validation deliberately reduces that reserve to
50 entries. When the accepted mesh exceeds the old capacity, JETSPIN detaches
the topology and evaporation mappings before the host allocations are
replaced. It then resizes the Platen workspaces, semantically repacks the
pre-generated Gaussian history for the new stride, generates values only for
new bead indices, and binds the completed host mesh and resized history once.
The same device Akima kernels operate before and after capacity growth.

Test Case 22 validates repeated use of this lifecycle. Its physical input is
derived from Test 21, while developer-only environment overrides reduce the
initial reserve and refinement growth increment to 20 entries. Three accepted
events must therefore perform three independent mapping releases, device Akima
remeshes, Gaussian-history repacks, workspace reallocations, and device
rebinds. Normal runs still reserve and grow by `incnpjet`, currently 100.

Test Case 23 adds collector removal without changing this hybrid boundary.
The device topology primitive advances the lower active bound and clears the
removed bead without downloading the full state. Accepted refinement events
still perform the only intermediate complete-state downloads. The test
requires removal on both sides of the third remesh, proving that a rebased
active interval survives the Akima and device-rebind lifecycle.

Two development builds validate this boundary. `JETSPIN_DEV_HOST_AKIMA`
retains the historical host coefficient/interpolation path as an oracle.
`JETSPIN_COMPARE_AKIMA` executes the host reference and the GPU path at each
accepted event and reports coefficient and interpolated-value errors for all
11 fields. Test 23 performs 33 such comparisons; its largest relative errors
are `2.48e-15` for coefficients and `3.32e-15` for values. The largest
absolute value difference is `7.28e-12`.

A third development build, `JETSPIN_COMPARE_REFINEMENT_ASSEMBLY`, applies the
same host-then-device oracle pattern to the volume/conservation/density
kernel: back up the pre-reconstruction state, evaluate the trusted host
reference, restore that state, run the device kernel so it stays
authoritative, then compare. All three Test 23 events agree with the host
reference at or near roundoff (worst absolute `7.1e-15`, worst relative
`4.05e-16`). Building this oracle caught two real ordering mistakes before
they could reach the standard path: an early version invoked the device
kernel a second time after the comparison helper already ran it, silently
reapplying the density-to-mass conversion; a later version evaluated the host
reference from the already-device-converted state instead of the shared
pre-reconstruction input. Both corrupted the event's mass/charge invariant by
tens of percent -- exactly the class of error this A/B methodology exists to
catch before a device kernel is trusted.

## Allocation and MPI interaction

Refinement can increase the number of beads beyond `mxnpjet`. In that case
`define_akima_bounds` grows capacity to the required upper bound plus the
normal `incnpjet` reserve. The remesher may also compact and rebase the active
interval without growing capacity.

The resulting flags have different meanings:

- `doreorder` records that bead indices were remapped;
- `doallocate` records that capacity-dependent consumers may require larger
  arrays;
- `lneighlistdo` requests a new neighbour list for multiple-timestep Coulomb
  calculations.

The fundamental jet state is replicated on every MPI rank. Interpolation
uses root-generated service data and collective communication to reproduce
the remeshed arrays. All ranks must leave refinement with identical active
state, bounds, and feature flags. Rank-local chunks and cached integrator
workspaces must then match the new `inpjet`, `npjet`, `mxnpjet`, and
`mxchunk` values.

See [dynamic allocation and bead indexing](dynamic-allocation.md) and
[MPI parallelization](parallelization.md) for the two supporting memory
models.

## Backup and recovery

The module provides grow-on-demand backup arrays through `store_backup` and
`restore_backup`. The saved state includes active bounds, simulation time,
fundamental bead fields, reference volume, and `jetve` when evaporation is
enabled. Restore broadcasts the recovered state so every rank resumes from
the same configuration.

Any new fundamental bead field that affects the equations must be considered
for both backup and restore, in addition to interpolation and allocation.

## Maintenance invariants

Changes to refinement should preserve these rules:

- the target path coordinates remain ordered and valid for Akima fitting;
- anchored and unrefined regions retain their intended physical locations;
- mass, charge, and volume are treated consistently as densities or
  per-bead quantities at each stage;
- `jetvl` and `jetve` are conserved and normalized independently;
- continuous fields are interpolated, while logical or identity fields are
  explicitly remapped;
- optional arrays are handled symmetrically on every allocation path;
- inactive array entries are cleared before collective summation;
- capacity growth updates MPI chunk sizing and all cached workspaces;
- code never assumes `inpjet == 0` or stable bead array indices;
- restart and backup formats are updated when new persistent state is added.

## Validation

The standard smoke suite exercises the historical refinement example and a
three-way Kelvin--Voigt, evaporation, and refinement path. The latter forces
both an Akima remeshing event and growth beyond the initial bead capacity.
Test Case 21 adds the pre-extended anchored Maxwell/Platen path and checks the
event-level remeshing invariants.
Test Case 22 extends that check to three consecutive events and three capacity
increases.
Test Case 23 combines those three events with collector removal before and
after the final remesh.

Run:

```sh
tests/smoke/run.sh serial
tests/smoke/run.sh debug
tests/smoke/run.sh mpi
tests/refinement/run.sh
tests/refinement/run.sh openacc
tests/refinement/run.sh force-oracle
tests/refinement/run_test22.sh nvfortran
tests/refinement/run_test22.sh openacc
tests/refinement/run_test22.sh force-oracle
tests/refinement/run_test23.sh nvfortran
tests/refinement/run_test23.sh openacc
tests/refinement/run_test23.sh force-oracle
tests/refinement/run_test23.sh host-akima
tests/refinement/run_test23.sh akima-compare
tests/refinement/run_test23.sh refinement-compare
```

For numerical changes, also compare conservation totals immediately before
and after refinement, verify that `nref` becomes non-zero, check for ordered
path coordinates and finite state, and compare serial and MPI trajectories
within an explicitly chosen tolerance.
