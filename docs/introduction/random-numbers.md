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

## Pre-generated Gaussian pool

Stochastic Platen runs that the OpenACC build can execute on the device
(`system 4`, `integrator 4`, serial) read their noise from a pre-generated
Gaussian pool instead of the per-step block above. Rank 0 generates the pool
once, before the timestep loop, and the OpenACC build copies it to the device
once. No random number is drawn inside the loop, and CPU and GPU runs read the
same values. The pool is used by:

- the fixed 1,000-bead Platen benchmarks, Tests 12 and 20;
- the dynamic Platen path (insertion, collector removal, dynamic
  refinement), evaporative (Tests 21--23 and 25) and, since 2026-10-05,
  non-evaporative (Test 24), whatever the initial bead count.

The pool is one flat sequence, allocated once and never remapped. Each
timestep reserves `(active beads) * 6` consecutive values:

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

in the host function `gaussian_history_value` and in the device kernels
`accelerator_platen_velocity` and `accelerator_platen_evap_velocity`.
Reservation and read share the same modulo, so a slice straddling the end of
the pool is served correctly.

How the pool is sized and filled depends on the topology:

- **Fixed topology** (Tests 12 and 20): the pool holds `nsteps * window * 6`
  values, 6,006,000 for Test 12 and 600,600 for Test 20, and each step's slice
  is filled in the draw order of the per-step block (bead, component, draw).
  A serial run reading the pool therefore uses exactly the numbers that an MPI
  or non-history run draws step by step; Tests 12 and 20 are bit-identical to
  their records from before the pool became the default.
- **Dynamic topology** (Tests 21--23 and 25): the bead count changes during
  the run, so the pool takes its whole size at once (100,000,000 values by
  default) and is filled in index order. Capacity growth needs no rebuild and
  no device remap: `resize_gaussian_history` is a no-op.

### The noise repeats after the pool is consumed

The cursor wraps once the pool has been consumed, so **the Gaussian sequence
repeats periodically**. The period is

```text
period (timesteps) = values / (6 * active beads)
```

With the default 100,000,000 values, a jet of 270 active beads repeats its
noise every about 62,000 timesteps; a 100-million-step Test 25 run goes
through the pool about 1,600 times. The noise is not independent beyond one
period, and the user must judge whether that is acceptable for the quantity
being studied.

The run log reports the pool once, before the loop:

```text
Gaussian history pool: values=100000000 covers 165016 steps at 101 beads
```

For a fixed topology the bead count is the active window and the step count
is exact. For a dynamic topology it is the initial capacity (here Test 25:
one active bead and 100 reserved slots), and the period at any other active
count follows from the formula above.

The lever is the pool size, set in `input.dat`:

```text
noise pool 4.d8
```

`noise pool` gives the number of pre-generated values, between 1,000,000 and
2,000,000,000 (default 100,000,000). Each value costs 8 bytes on the host and,
in an OpenACC build, 8 bytes again on the device: 4.d8 values need about
3 GiB on each side, and rank 0 needs a time proportional to the pool size to
generate them before the loop (about 5 s for the default pool with the
NVFORTRAN CPU build).
A larger pool lengthens the period proportionally; it does not change the
statistics of the noise within a period. The directive has no effect on runs
that do not use the pool.

Until 2026-09-30 the default layout indexed a four-dimensional history
`(step, bead, component, draw)` with a bead stride of `mxnpjet + 1`: every
capacity change forced a host repack and a device remap, and the covered cycle
shrank as capacity grew (74,404, then 51,440, then 36,791 steps in one
observed run). The sequential pool was then available only in the
`nvfortran-openacc-dynamic-platen` build. The dynamic evaporative Platen path
also used the history only when the jet already had 100 beads at the one-shot
allocation decision, so a single-bead start never used it; that decision is
now size-independent.

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
therefore follows the current global indices. The pre-generated pool does not
depend on capacity at all: each step reads the next `(active beads) * 6`
values, so capacity growth leaves it untouched on the host and on the device.

Test Case 23 is also a debugging guard: CPU and accelerator comparisons must
use the same pre-generated-history policy before a difference in collector
removal timing is attributed to topology code. Falling back to on-the-fly CPU
draws changes the stochastic trajectory before the first removal.

The current scheme associates a draw with the bead's array index during that
step. It guarantees serial/MPI agreement for the same topology evolution; it
does not define a permanent stochastic identity for a bead across a remeshing
event. See [dynamic allocation and bead indexing](dynamic-allocation.md).

Capacity can grow through two independent routes, dynamic refinement and
insertion overflow in `reallocate_jet`; both still call
`resize_gaussian_history`, which is a no-op for the pool. With the former
stride layout, a run that grew through insertion alone once kept a history
sized for the previous `mxnpjet` and indexed past its end
(`CUDA_ERROR_ILLEGAL_ADDRESS` inside `accelerator_platen_velocity`).

## Restart

Since 2026-10-06 the restart file `save.dat` ends with a versioned block
(restart state version 1, after the bead records) that holds the two random
states of a run, together with the counters of the dynamic refinement and of
the anchor tagging:

- the cursor of the pre-generated pool and the pool size. The cursor is the
  cumulative sum of the active bead counts and cannot be derived from the
  step number. A restarted run regenerates the pool from the seed of
  `input.dat` and then moves the cursor to the saved position; a different
  `noise pool` size prints a warning and the cursor is taken modulo the new
  size.
- the state of the intrinsic generator (`random_seed(get=...)`), from which
  the runs without a pool draw their noise step by step. It is restored
  after the pool has been drawn from the input seed. Its layout belongs to
  the compiler's runtime: when the size differs (an executable built with
  another compiler) a warning is printed and the generator is not restored.

With the same input file and executable, a restarted run therefore continues
the uninterrupted one exactly (`tests/restart/run.sh`, Tests 24 and 25).
Files written before 2026-10-06 have no such block; they are still read,
with a warning, and the noise then restarts from the beginning of the pool
and from the input seed.

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
- Capacity changes must resize the per-step buffer before indexed access, on
  every route that can change capacity, not only dynamic refinement.
- Every integrator path that reads the pre-generated pool must call
  `begin_gaussian_history_step` exactly once per timestep, before any read.
- A pre-generated history must be allocated by a decision that does not depend
  on the bead count, because that decision is taken once before the timestep
  loop, when a jet growing from a single bead has not yet reached any
  size-based eligibility threshold.
- Every random state that a run advances (pool cursor, intrinsic generator)
  must be written to the restart file; restart state version 1 holds both. A
  new state means a new version, read only when present.

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
