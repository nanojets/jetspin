# Random numbers and MPI reproducibility

JETSPIN uses pseudo-random numbers for nozzle perturbations, stochastic
forcing, and optional variable-mass insertion. A fixed input `seed` makes a
run repeatable only when every execution consumes the random stream in the
same logical order. MPI work partitioning must therefore not determine which
random value is assigned to a physical bead.

The generator and shared Gaussian buffer are implemented in
[`utility_mod.f90`](../../source/utility_mod.f90). Their use by the stochastic
Platen integrators is in
[`integrator_mod.f90`](../../source/integrator_mod.f90).

## Rank-independent Gaussian blocks

Before each stochastic Platen step, rank 0 generates a block indexed by:

```text
global bead index, Cartesian component, Gaussian draw number
```

One-dimensional systems reserve two Gaussian values per active bead. A
three-dimensional system reserves two values for each of its three
components. Rank 0 advances the random stream in global bead order and then
broadcasts the complete block. Every rank reads values using the global bead
index rather than a rank-local loop counter.

```text
rank 0: generate one globally ordered Gaussian block
                           |
                           v
                  broadcast the block
                           |
          +----------------+----------------+
          v                v                v
       rank 0           rank 1           rank N
   read owned beads  read owned beads  read owned beads
```

This makes the stochastic assignment independent of `mystart`, `myend`, and
the number of MPI ranks. The allocated buffer follows `mxnpjet` and is reused
until jet capacity changes. Its maximum storage is six double-precision
values per allocated bead, rather than a very large simulation-wide random
archive.

## Fixed GPU benchmark histories

Tests 12 and 20 deliberately use a different storage policy because their
topology and duration are fixed. Before timing, rank 0 generates every
timestep block in `step, bead, component, draw` order. CPU and GPU paths
address this same history; the OpenACC build copies it to the device once.
Test 12 stores 6,006,000 doubles for 1,000 steps, while Test 20 stores 600,600
doubles for 100 steps. Test 21 uses the same pre-generated policy with dynamic
topology: its standard 500-entry capacity stores 24,048,000 doubles for 8,000
steps. This guarantees indexed Gaussian assignments without a per-step
host/device transfer.

Test 22 applies the dynamic policy three times in one run. Each capacity
increase repacks retained values into the new stride and generates only the
new bead-index slots before a single device upload.

The pre-generated history is capped at 100,000,000 double precision values
(about 763 MiB). If a fixed simulation needs more timesteps than fit in that
limit, timestep addressing wraps to the first stored block. CPU and GPU remain
reproducible, but the Gaussian sequence then repeats periodically and should
not be interpreted as independent noise beyond that period.

## Sequential pool for the dynamic Platen fork

The layout above indexes the pre-generated history as a four-dimensional array
`(step, bead, component, draw)` whose bead stride is hard-wired to
`mxnpjet + 1`. That stride is the reason a capacity change forces a full host
repack plus a device remap. It also means the covered timestep count *shrinks*
as capacity grows, because the 100,000,000-value cap is divided by a
per-step size that follows allocated capacity rather than live beads: one
observed run went from 74,404 to 51,440 to 36,791 covered steps without ever
reading the added slots.

The `nvfortran-openacc-dynamic-platen` build (macro `JETSPIN_GPU_DYNAMIC_PLATEN`,
see [compiling](compiling.md)) replaces this with a single flat random
sequence, allocated once and never remapped. Each timestep reserves
`(active beads) * 6` consecutive values:

```text
begin_gaussian_history_step(mystart, myend)
    base   <- cursor            (start of this step's slice)
    window <- myend - mystart + 1
    first  <- mystart           (bead index that maps to offset zero)
    cursor <- mod(cursor + window*6, values)
```

`begin_gaussian_history_step` is called once per step ahead of the `systype`
dispatch in `platen()` and `platen_ev()`, so it also covers the persistent
device branch that returns before the host stage code. A read resolves to

```text
index = mod(base + (ipoint - first)
            + window*((icomponent - 1) + 3*(idraw - 1)), values)
```

The cursor walks forward and wraps only once the entire sequence has been
consumed, and reservation and read share the same modulo, so a slice
straddling the end of the pool is served correctly instead of discarding its
tail. `resize_gaussian_history` becomes a no-op under the macro: the remap
hazard is removed structurally rather than patched at each call site.

Because consumption tracks live beads instead of reserved capacity, the pool
also lasts longer. At 351 active beads with capacity 452 it reserves 2,106
values per step (about 47,500 steps before wrapping) against 2,718 values per
step (36,791 steps) for the stride layout.

Both device kernels, `accelerator_platen_velocity` and
`accelerator_platen_evap_velocity`, take three extra arguments carrying the
slice (`base`, `window`, `values`) and compute the index through one shared
expression covering both layouts. For the historical layout the index is
always below the array size, so the `mod` is an identity and macro-off results
are unchanged; this is confirmed by regression comparison against the previous
commit.

The pool consumes a different noise ordering from the stride layout, so
macro-on Platen trajectories are statistically equivalent but not bit-identical
to macro-off ones. Stored references therefore do not need regenerating unless
the pool layout is ever promoted to the default, at which point the contract
described in this document is what changes.

## Other random paths

Nozzle perturbation draws are generated on rank 0 and broadcast as scalar
values. Any future random path must follow the same rule or use the shared
globally indexed buffer. Calling `gauss()` independently inside a rank-local
bead loop is not MPI reproducible.

The intrinsic generator is initialized separately on each rank for legacy
compatibility, with a rank-dependent seed. This is safe only when non-root
streams do not determine replicated physical state. Rank 0 has the same seed
as the serial build and is the authoritative generator for shared stochastic
data.

## Dynamic topology

Insertion, removal, compaction, and refinement can change bead indices and
capacity. The ordinary per-step buffer is prepared after chunk assignment and
therefore follows the current global indices. For Tests 21--23, capacity
growth of the pre-generated history instead repacks every retained step into
the new bead stride, preserves values for existing indices, and generates one
extension for newly available indices. In an OpenACC run, the old history
mapping is deleted before the host allocation changes and the resized history
is copied to the device once.

Test Case 23 is also a debugging guard: CPU and accelerator comparisons must
use the same pre-generated-history policy before a difference in collector
removal timing is attributed to topology code. Falling back to on-the-fly CPU
draws changes the stochastic trajectory before the first removal.

The current scheme associates a draw with the bead's array index during that
step. It guarantees serial/MPI agreement for the same topology evolution; it
does not define a permanent stochastic identity for a bead across a remeshing
event. See [dynamic allocation and bead indexing](dynamic-allocation.md).

Capacity can grow through two independent routes, and the stride layout must
be resized on both: dynamic refinement, and insertion overflow in
`reallocate_jet`. Only the first called `resize_gaussian_history` until this
was corrected; a run that grew through insertion alone kept a history sized
for the previous `mxnpjet` and indexed past its end, which on the GPU
surfaces as `CUDA_ERROR_ILLEGAL_ADDRESS` inside
`accelerator_platen_velocity`. `reallocate_jet` now calls
`resize_gaussian_history(mxnpjet)` whenever it sets `doallocate`.

## Restart limitation

The historical restart format does not serialize the Fortran intrinsic
random-generator state. A restarted stochastic run can therefore continue
from the saved physical state but is not guaranteed to reproduce the exact
uninterrupted random sequence. Exact stochastic restart would require saving
and restoring `random_seed(get=...)` state, with an explicit restart-format
version and compatibility policy.

The sequential pool adds a second, distinct restart gap. Its cursor is
genuinely stateful — it is the cumulative sum of active bead counts, not
derivable from the step number — and it is not written to the checkpoint, so
a restarted macro-on run resumes the noise sequence from offset zero. This is
harmless statistically but breaks bit-exact restart reproducibility, and would
have to be addressed before the pool layout could become the default.

## Maintenance invariants

- Only rank 0 generates random values that affect replicated physical state.
- Shared random values are broadcast before any rank consumes them.
- Random values inside distributed loops are addressed by global bead index,
  component, and draw role, never by a local consumption counter.
- All ranks prepare the same buffer shape and enter its broadcast in the same
  collective order.
- A new stochastic integrator must document how many draws it reserves per
  bead and stage.
- Optional branches must not shift later assignments differently on
  different ranks.
- Capacity changes must resize the buffer before indexed access, on every
  route that can change capacity, not only dynamic refinement.
- A pre-generated history must be allocated by a decision that does not depend
  on the bead count, because that decision is taken once before the timestep
  loop, when a jet growing from a single bead has not yet reached any
  size-based eligibility threshold.
- Restart reproducibility must not be claimed until generator state is part
  of the restart format.

## Regression coverage

The numerical regression suite modifies Test Case 4 to start with 25 beads.
On two ranks this exceeds the ten-bead minimum chunk and forces both ranks to
integrate stochastic beads. The MPI trajectory is compared with the matching
serial trajectory using the configured numerical tolerance.

Run this focused check with:

```sh
JETSPIN_REGRESSION_FIRST_CASE=4 \
JETSPIN_REGRESSION_LAST_CASE=4 \
tests/regression/run.sh
```
