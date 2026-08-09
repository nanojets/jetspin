# MPI parallelization

JETSPIN uses a **replicated-data MPI strategy with bead-level work
partitioning**. It is not a spatial domain-decomposition code.

Every MPI rank stores the complete fundamental state of the jet, including
the active bead coordinates, velocities, stresses, masses, charges, reference
volumes, and—when enabled—post-evaporation volumes. Expensive calculations
are divided among ranks, and collective operations reconstruct a consistent
full state on every rank.

The main communication abstraction is implemented in
[`parallel_version_mod.f90`](../../source/parallel_version_mod.f90). The
serial build substitutes [`serial_version_mod.f90`](../../source/serial_version_mod.f90),
which exposes the same interface without MPI communication.

## Execution model

At startup, JETSPIN initializes MPI, obtains the rank identifier `idrank` and
the number of ranks `mxrank`, and allocates the domain bookkeeping arrays.
For each active jet configuration, `set_chunk` assigns bead indices to the
ranks.

The normal integration path follows this pattern:

```text
complete jet state replicated on every rank
                  │
                  ▼
       assign a subset of beads to each rank
                  │
                  ▼
       calculate local derivative contributions
                  │
                  ▼
       MPI_Allreduce partial result arrays
                  │
                  ▼
complete updated jet state replicated on every rank
```

This design makes the complete geometry immediately available for non-local
interactions, especially the long-range Coulomb force.

## Contiguous bead chunks

`set_chunk(inpjet, npjet)` normally divides the active bead interval into
contiguous ranges. Each rank receives:

- `mystart:myend`, the bead interval owned for computation;
- `memystart:memyend`, an extended interval that includes up to one adjacent
  bead at each boundary for calculations involving neighbouring elements.

The chunk size is balanced across ranks, subject to the historical minimum
chunk size `nchunkmin = 10`. When the jet contains too few beads to give each
rank a useful chunk, some higher-numbered ranks may receive no active beads.
MPI efficiency therefore depends strongly on having enough beads per rank.

Temporary work arrays can use `mxchunk`, which is derived from the maximum
jet capacity and the number of ranks. Fundamental jet arrays remain fully
replicated.

The time integrators in [`integrator_mod.f90`](../../source/integrator_mod.f90)
and [`integrator_kv_ev_mod.f90`](../../source/integrator_kv_ev_mod.f90)
typically iterate over `mystart:myend`. Their partial stage arrays are then
combined with `sum_world_darr`.

## Cyclic Coulomb-force distribution

The direct Coulomb calculation uses a cyclic distribution for its outer bead
loop:

```fortran
do ipoint = inpjet + idrank, npjet, mxrank
```

Each rank therefore handles every `mxrank`-th outer-loop bead while reading
the complete replicated jet geometry. The partial force arrays are summed
with `MPI_Allreduce`. This cyclic assignment generally balances the
triangular pair-interaction workload better than a simple contiguous block.

The implementation is in
[`coulomb_force_mod.f90`](../../source/coulomb_force_mod.f90).

## Collective communication

The MPI version wraps the main collective operations behind routines in
`version_mod`:

| Wrapper family | MPI operation | Purpose |
| --- | --- | --- |
| `sum_world_*` | `MPI_Allreduce` with `MPI_SUM` | Assemble partial arrays and scalar contributions on every rank |
| `min_world_*`, `max_world_*` | `MPI_Allreduce` | Compute global extrema |
| `and_world_*`, `or_world_*` | `MPI_Allreduce` | Combine logical decisions |
| `bcast_world_*` | `MPI_Bcast` from rank 0 | Replicate root-generated state or decisions |
| `get_sync_world` | MPI barrier | Explicit synchronization |

Many integrator stages perform several full-array reductions. Communication
cost can consequently become significant when the bead count per rank is
small.

## Dynamic topology and refinement

Bead insertion, removal, and dynamic refinement change the active index
range. Code that changes `npjet` or the maximum allocated capacity must keep
the domain metadata and local workspaces synchronized:

1. update or reallocate the replicated jet state;
2. call `set_chunk` for the new active bead range;
3. call `set_mxchunk` when the maximum capacity changes;
4. resize local integrator or force workspaces when required;
5. synchronize root-generated arrays before continuing collective stages.

Dynamic refinement performs several operations on rank 0 and broadcasts the
reconstructed arrays afterward. When evaporation is active, both the
reference volume `jetvl` and instantaneous volume `jetve` must be preserved
and communicated. The relevant implementation is in
[`dynamic_refinement_mod.f90`](../../source/dynamic_refinement_mod.f90).
The underlying capacity, compaction, and index-rebasing rules are documented
in [dynamic allocation and bead indexing](dynamic-allocation.md). The Akima
remeshing sequence and conservation rules are documented separately in the
[dynamic-refinement guide](dynamic-refinement.md).

## Scalability characteristics

The replicated strategy has deliberate trade-offs:

### Advantages

- simple access to the complete jet geometry on every rank;
- straightforward implementation of long-range pair interactions;
- no explicit halo exchange for the fundamental jet arrays;
- the same high-level algorithms can use serial and MPI communication
  backends.

### Limitations

- every rank stores the complete fundamental jet state;
- full-array `MPI_Allreduce` calls occur repeatedly during multi-stage
  integration;
- the direct Coulomb calculation remains quadratic in the number of beads;
- ranks may be idle for small systems;
- performance is limited by collective latency when there are too few beads
  per rank.

The scientific manual historically recommends parallel execution only when
the simulation has roughly 50 or more beads per rank. This is a practical
guideline rather than a hard program requirement and should be confirmed by
benchmarking on the target machine.

## Maintenance invariants

Changes to MPI-sensitive code should preserve the following rules:

- Every rank must call collective operations in the same order and with the
  same element count.
- Arrays assembled with `MPI_SUM` must contain zero outside the locally
  computed contribution, unless the routine intentionally supplies a
  separate output buffer.
- Fundamental state must be identical on all ranks before a calculation that
  reads non-local bead data.
- Random values that affect replicated state must be generated in a global
  order and addressed independently of rank ownership; see
  [random numbers and MPI reproducibility](random-numbers.md).
- Changes to bead count or capacity must update chunk metadata and all local
  workspaces.
- New replicated state variables must be included in restart, refinement,
  backup/restore, and broadcast paths where applicable.
- The serial and MPI versions of `version_mod` must retain compatible public
  interfaces.
- Root-only I/O or topology changes must be followed by the required
  broadcasts before distributed computation resumes.

## Building and testing MPI

Build the OpenMPI/GFortran version from the repository root:

```sh
make -C source -f ../build/Makefile gfortran-mpi
```

Run the MPI smoke test with:

```sh
tests/smoke/run.sh mpi
```

The smoke test builds the MPI executable and runs a representative case on
two processes. Numerical or communication changes should also be compared
against serial results when practical.
