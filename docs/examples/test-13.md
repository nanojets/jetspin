# Test Case 13: dynamic-topology CPU/GPU benchmark

Test 13 starts with 1,024 beads uniformly distributed from nozzle to collector
and enables both `inserting yes` and `removing yes`. It retains the constant
axial field and material model of Test 9 but uses an intentionally artificial
2,000 cm/s nozzle velocity so that a short 1,000-step run exercises topology.

The CPU build and the OpenACC device step make 13 insertions and 13
removals at the same steps, finish with 1,024 active beads, and remain
between 1,023 and 1,025 active beads; the original call-scoped and
persistent paths that preceded the device step made the same 13 insertions
and 13 removals. The program prints a structured `Topology event:` record
whenever either boundary changes the topology.

The second dynamic milestone reserves capacity for 1,280 beads before the
OpenACC mapping and keeps RK4, Coulomb, and EOM arrays resident while
`inpjet` and `npjet` change. Collector detection and clamping execute on the
device. Nozzle distance checks, release of the blocked bead, initialization
of the new record, the `npjet` update, and the placement and charge
smoothing of the blocked bead before every force evaluation also execute on
the device. The host receives one 20-byte topology record each step and the
two new tail records when an insertion actually occurs. Removed records are
synchronized individually for removal output and then cleared on both host
and device. This run stays within its reserve; capacity growth, which
reallocates the arrays on the host and rebuilds the device mapping, is
exercised by [Test 15](test-15.md). Since 2026-10-06 the path is the common
device step of all RK runs (`device_rk_step` in
`source/device_step_mod.f90`), engaged from the first step.

The original call-scoped NVFORTRAN 24.3/A30 comparison measured `40.659696 s`
on CPU and `3.834931 s` with OpenACC and reproduced the CPU topology stream.
With host-side insertion the bounded persistent A30 path took `2.829376 s`.
The final fully device-side insertion path takes `2.623558 s`
(`381.162 steps/s`); those timings predate the 2026-10-01 changes below.
With the Coulomb kernel of 2026-10-05 (one gang per target, see
[Test 9](test-9.md)) Tests 9-14 run 4.2 to 4.8 times faster on the A30 than
with the previous kernel.

The persistent A30 run reproduces the CPU topology stream (26 events) step
for step, and its statistics pass the paired `rtol=3e-4` comparison with the
NVFORTRAN CPU run (worst normalized difference 0.94). Until 2026-10-01 its
insertions fell progressively earlier (step 593 instead of 600 for the
seventh). The persistent path placed the blocked nozzle bead only when it was
created and smoothed its charge on the stale host copy, so the device Coulomb
sum saw a full charge at the nozzle. That drift had been attributed to GPU
RK4 rounding; the diagnosis that the EOM itself agrees with the CPU to about
machine precision was correct, the attribution of the remaining difference
was not.

`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` keeps the device step closed: the
OpenACC build then runs the CPU build's code, with the three-dimensional
Coulomb sums offloaded (with copies) from 128 active beads.

Validation tools are in `tests/performance/dynamic/`: `topology-events.txt`
is the reference event stream (the other `topology-events*.txt` files keep
the historical streams of the persistent paths);
`compare.sh CPU/run.log GPU/run.log CPU/statout.dat GPU/statout.dat` checks
the `Topology event:` lines of both logs against it and compares the two
statistics files with `rtol=3e-4` and `atol=1e-6` (`RTOL`, `ATOL` override
them); `compare_state.py` compares two `topology-state.dat` files written
with `JETSPIN_TOPOLOGY_SNAPSHOT=1`.

- [Input file](../../examples/input-13/input.dat)
- [Comparison record](../../tests/performance/dynamic/README.md)
