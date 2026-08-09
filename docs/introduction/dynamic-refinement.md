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
| `dynamic refinement every <time>` | Interval between refinement checks |
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
reference refinement case.

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
The counter `irefinementdone`, available as the `nref` observable, increases
only when a refinement is accepted.

## Path parametrization and anchors

`compute_length_path` assigns each bead a normalized path coordinate `jetpt`.
The remesher uses this coordinate, rather than an individual Cartesian axis,
so the same algorithm works for one- and three-dimensional trajectories.

The refined region extends from the collected/nozzle side through
`irefbeadstart`, the last element found above the threshold. Beads marked in
`jetbd` divide it into anchored segments. Each segment receives a number of
new intervals based on its fraction of total path length and the requested
threshold. Anchor locations are copied into the target mesh exactly; the
unrefined tail is retained rather than interpolated.

Anchor beads reduce interpolation drift and preserve selected material
locations, but their array indices may still change when the active interval
is rebased.

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

Run:

```sh
tests/smoke/run.sh serial
tests/smoke/run.sh debug
tests/smoke/run.sh mpi
```

For numerical changes, also compare conservation totals immediately before
and after refinement, verify that `nref` becomes non-zero, check for ordered
path coordinates and finite state, and compare serial and MPI trajectories
within an explicitly chosen tolerance.
