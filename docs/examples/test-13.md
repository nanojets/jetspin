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
of the new record, and the `npjet` update also execute on the device. The host
receives only topology scalars each step and the two new tail records when an
insertion actually occurs. Removed records are synchronized individually for
removal output and then cleared on both host and device. Capacity growth and
general compaction remain future work.

The original call-scoped NVFORTRAN 24.3/A30 comparison measured `40.659696 s`
on CPU and `3.834931 s` with OpenACC and reproduced the CPU topology stream.
With host-side insertion the bounded persistent A30 path took `2.829376 s`.
The final fully device-side insertion path takes `2.623558 s`
(`381.162 steps/s`).
Moving the RK4 intermediate updates from host to GPU changes floating-point
rounding; the insertion threshold amplifies this into progressively earlier
insertion events, although event counts, removals, active bounds, and final
bead count remain unchanged. The old and new streams are stored separately.

A full-state event-by-event diagnosis confirmed exact bead indices, frozen
flags, masses, charges, volumes, and event metadata. The first CPU/device
difference occurs in a transverse quantity close to zero at step 40. The same
OpenACC EOM executed on the CPU agrees with the original CPU implementation to
about machine precision, identifying device floating-point evaluation of the
EOM/curvature kernel—not insertion or removal—as the source subsequently
amplified by the dynamic trajectory.

`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` restores the call-scoped diagnostic
path. On the same A30 it exactly restores the previous 26 event timesteps,
which isolates the changed topology timing to GPU RK4 arithmetic rather than
to removal or Coulomb indexing. Moving the insertion distance test itself to
the GPU shifts one later threshold crossing from step 910 to step 909; counts,
bounds, and final topology remain unchanged.

- [Input file](../../examples/input-13/input.dat)
- [Comparison record](../../tests/performance/dynamic/README.md)
