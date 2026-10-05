# Test Case 13: dynamic-topology CPU/GPU benchmark

Test 13 starts with 1,024 beads uniformly distributed from nozzle to collector
and enables both `inserting yes` and `removing yes`. It retains the constant
axial field and material model of Test 9 but uses an intentionally artificial
2,000 cm/s nozzle velocity so that a short 1,000-step run exercises topology.

Both the original call-scoped and persistent runs produce 13 insertions and
13 removals, finish with 1,024 active beads, and remain between 1,023 and
1,025 active beads. The program prints a structured `Topology event:` record
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
and device. Capacity growth and general compaction remain future work.

The original call-scoped NVFORTRAN 24.3/A30 comparison measured `40.659696 s`
on CPU and `3.834931 s` with OpenACC and reproduced the CPU topology stream.
With host-side insertion the bounded persistent A30 path took `2.829376 s`.
The final fully device-side insertion path takes `2.623558 s`
(`381.162 steps/s`); those timings predate the 2026-10-01 changes below.

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

`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` restores the call-scoped diagnostic
path.

- [Input file](../../examples/input-13/input.dat)
- [Comparison record](../../tests/performance/dynamic/README.md)
