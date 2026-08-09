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

## Fixed GPU benchmark history

Test 12 deliberately uses a different storage policy because its topology and
duration are fixed. Before timing, rank 0 generates all 1,000 timestep blocks
in `step, bead, component, draw` order. CPU and GPU paths address this same
history; the OpenACC build copies its 6,006,000 doubles to the device once.
This guarantees identical Gaussian assignments without a per-step host/device
transfer. General Platen simulations with dynamic topology continue to use the
per-step block described above.

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
capacity. A fresh Gaussian block is prepared after chunk assignment for the
current integration step, so it uses the current global indices consistently
on all ranks. Capacity growth automatically reallocates the shared block.

The current scheme associates a draw with the bead's array index during that
step. It guarantees serial/MPI agreement for the same topology evolution; it
does not define a permanent stochastic identity for a bead across a remeshing
event. See [dynamic allocation and bead indexing](dynamic-allocation.md).

## Restart limitation

The historical restart format does not serialize the Fortran intrinsic
random-generator state. A restarted stochastic run can therefore continue
from the saved physical state but is not guaranteed to reproduce the exact
uninterrupted random sequence. Exact stochastic restart would require saving
and restoring `random_seed(get=...)` state, with an explicit restart-format
version and compatibility policy.

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
- Capacity changes must resize the buffer before indexed access.
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
