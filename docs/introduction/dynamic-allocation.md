# Dynamic allocation and bead indexing

JETSPIN stores the jet in dynamically allocated, zero-based arrays. Its
memory management is unusual because the allocated capacity, the active bead
interval, and the MPI work allocation are three different concepts. Code must
not assume that the first active bead is at index zero or that the last active
index is the allocated upper bound.

The main implementation is in
[`nanojet_mod.f90`](../../source/nanojet_mod.f90), with refinement-specific
reallocation in
[`dynamic_refinement_mod.f90`](../../source/dynamic_refinement_mod.f90).

## Core indices and flags

| Name | Meaning |
| --- | --- |
| `inpjet` | Index of the first active bead |
| `npjet` | Index of the last active bead/current active upper endpoint |
| `mxnpjet` | Allocated upper bound of the fundamental jet arrays |
| `incnpjet` | Capacity-growth increment; `dynamic refinement capacity <i>` input directive, default 100 beads |
| `reallocation_increment` | Capacity increase applied by `reallocate_jet`; currently 100 beads |
| `mxchunk` | Maximum size of a rank-local MPI work chunk |
| `doallocate` | Signals that capacity-dependent work arrays may need resizing |
| `doreorder` | Signals that bead indices have been rebased or remapped |

The fundamental arrays are allocated with bounds `0:mxnpjet`, while the live
state normally occupies only `inpjet:npjet`:

```text
0                 inpjet                 npjet              mxnpjet
| reclaimed slots |------ active jet ------| spare capacity |
```

Consequently, `mxnpjet + 1` is the number of allocated elements, whereas the
active interval contains `npjet - inpjet + 1` elements when both endpoints
are included by the operation in question.

## Initial allocation

`allocate_jet` establishes the problem dimension, initializes capacity in
blocks of `incnpjet`, updates `mxchunk`, and allocates the fundamental state.
Coordinates, velocities, stress, mass, charge, radius, reference volume, and
related bead data all follow the same zero-based capacity.

Several arrays exist only when their corresponding feature is active:

- `jetfm` for multiple-timestep calculations;
- `jetce` and `jetve` for evaporation;
- `jetbd` for variable mass or bead tagging;
- `jetbr` for breakup;
- `jetlb` for bead tracking, allocated only on rank 0.

Allocation, copying, remeshing, restart, and deallocation must remain
symmetric for these optional fields.

## Insertion, removal, and compaction

Bead removal normally does not shift arrays or shrink their capacity. A bead
is marked as collected, its state is cleared where appropriate, and `inpjet`
advances. This leaves reusable storage before the active interval.

Bead insertion advances `npjet`. If it passes `mxnpjet`, `reallocate_jet` is
called. Despite its name, this routine performs two related operations:

1. it compacts and rebases the active data toward index zero, normally
   retaining one bead immediately before `inpjet` for neighbour-dependent
   calculations;
2. it increases `mxnpjet` by `incnpjet` only when compaction does not recover
   enough space.

Therefore `reallocate_jet` can reorder every bead-aligned array without
increasing the physical allocation. It copies data through shared service
buffers, clears or reallocates the destination arrays, restores the shifted
active interval, updates `inpjet` and `npjet`, and sets `doreorder`. When the
capacity changes it also sets `doallocate` and recomputes `mxchunk`.

This distinction is important for output and tracking code: a stable physical
bead cannot be identified permanently by its current array index.

## Dynamic refinement

Dynamic refinement is a separate topology-changing path. It estimates the
new remeshed upper bound and, when necessary, sets capacity to that bound plus
`incnpjet`. Akima interpolation then reconstructs the bead-aligned fields in
a rebased interval.

Production runs set `incnpjet` with the `dynamic refinement capacity <i>`
input directive (positive integer, default 100 beads if omitted); it applies
uniformly to the initial reserve and every later growth increment. Test Case
22 additionally overrides only the refinement growth increment through the
developer environment variable `JETSPIN_REFINEMENT_GROWTH_INCREMENT`,
independently of the input value, to stress-test several release/rebind
cycles in a short run.

The full workflow and its conservation rules are documented in the
[dynamic-refinement guide](dynamic-refinement.md).

Refinement uses its own grow-on-demand service and backup arrays. It may
change both the bead count and index mapping, so it propagates `doreorder` and
combines its local allocation decision with `doallocate`. When evaporation is
enabled, reference volume `jetvl` and instantaneous post-evaporation volume
`jetve` are distinct state and both must survive remeshing, backup, restore,
and MPI synchronization.

In the single-GPU Maxwell/Platen path, an accepted refinement first downloads
the active state for host target-mesh preparation. Akima coefficient
construction and field interpolation then execute on the GPU. If the target
mesh exceeds `mxnpjet`, the old topology and evaporation mappings are deleted
before any host allocation changes. The new arrays, capacity-dependent Platen
workspace, and indexed Gaussian history are then rebound once. No stale device
address is retained across the host reallocation.

Test Case 23 verifies the same replacement while the active lower bound moves.
Collector removal is handled by the existing device topology primitive: it
advances `inpjet`, clears the removed slot, and returns only event metadata.
It does not require the complete jet to be copied to the host.

## Relationship with MPI workspaces

Every MPI rank holds the complete fundamental arrays through `mxnpjet`, but
many derivative and force work arrays are rank-local and sized through
`mxchunk`. After a topology or capacity change:

1. `set_mxchunk(mxnpjet)` updates the possible local workspace size when
   capacity changes;
2. `set_chunk(inpjet, npjet)` redistributes the current active interval;
3. integrators and force modules resize cached workspaces when signalled.

The Kelvin--Voigt evaporation integrator explicitly tracks both the last
allocated `mxnpjet` and `mxchunk`, because dynamic refinement can invalidate
either family of arrays. See the
[MPI parallelization guide](parallelization.md) for the replicated-data and
collective-communication model.

## Maintenance checklist

When adding a bead-aligned state variable, inspect every relevant lifecycle
path rather than changing only `allocate_jet`:

- initial allocation, initialization, clearing, and final deallocation;
- bead insertion, collection, and removal;
- compaction in `reallocate_jet`;
- dynamic-refinement interpolation or discrete remapping;
- refinement backup and restore;
- restart input and output;
- MPI broadcasts or reductions;
- integrator stage, derivative, and force workspaces;
- optional-feature allocation/deallocation symmetry.

Preserve these invariants:

- arrays remain zero-based unless an algorithm explicitly documents another
  convention;
- loops use the correct active interval and never assume `inpjet == 0`;
- allocated capacity is not treated as the active bead count;
- inactive entries supplied to an `MPI_SUM` operation are zero;
- `set_mxchunk` follows a capacity change and `set_chunk` follows an active
  interval change;
- `doallocate` reaches every capacity-dependent consumer;
- `doreorder` reaches code whose identifiers or output depend on bead indices.

Relevant validation includes serial and runtime-checking smoke tests, the
dynamic-refinement cases, Kelvin--Voigt evaporation with refinement, and an
MPI comparison against the serial result.
