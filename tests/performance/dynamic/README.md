# Dynamic-topology benchmark

Test 13 initializes 1,024 beads across the nozzle--collector distance and uses
an intentionally high 2,000 cm/s injection speed to force topology activity
within 1,000 RK4 steps. The reference event sequence contains 13 insertions
and 13 removals, ends with 1,024 active beads, and never falls below 1,023.

The current milestone reserves 1,280 bead slots and keeps the RK4, Coulomb,
EOM, state, and scratch arrays persistently mapped while the active bounds
change. Removal detection and collector clamping run on the device and return
one decision scalar. Nozzle distance checks, blocked-bead release, record
initialization, and the `npjet` update also run on the device. Only topology
scalars return every step; the two new tail records are downloaded when an
insertion actually occurs. Removed records are downloaded individually for
removal output and cleared on both sides.
Reallocation beyond the reserved capacity and general compaction are not yet
device-resident.

The device-insertion persistent A30 path completes in `2.623558 s`
(`381.162 steps/s`), compared with `2.829376 s` for persistent RK4 with
host-side insertion. It
retains 13 additions, 13 removals, and 1,024 final active beads. Its insertion
events occur progressively earlier than in the call-scoped baseline because
the threshold is sensitive to GPU RK4 rounding. The two acceptance streams
are therefore stored separately:

- `topology-events.txt`: previous call-scoped/CPU topology baseline;
- `topology-events-persistent-host-insertion-a30.txt`: intermediate
  persistent RK4 baseline with host-side insertion;
- `topology-events-persistent-a30.txt`: current device-insertion A30 baseline.

Setting `JETSPIN_OPENACC_DISABLE_PERSISTENT=1` restores the call-scoped path
and exactly reproduces the previous event stream on the same A30. Replacing
only the device cross-section calculation with the historical host routine
does not change the persistent stream. Together these controls isolate the
timing change to GPU RK4 arithmetic, not topology or Coulomb indexing.
The insertion distance test on the GPU moves only the later step-910 crossing
to step 909 relative to the intermediate host-insertion baseline.

A paired NVFORTRAN 24.3/A30 run measured `40.659696 s` on CPU and `3.834931 s`
with OpenACC. The topology event streams were identical. Dynamic curvature
amplifies small floating-point differences, so the paired observable check
uses `rtol=3e-4` and `atol=1e-6`; this is deliberately separate from the
normal `1e-6` regression criterion.

## Full-state diagnosis

Set `JETSPIN_TOPOLOGY_SNAPSHOT=1` to write `topology-state.dat`. At every
event it records every active bead's index, frozen flag, position, stress,
velocity, mass, charge, and volume. Compare two files with:

```sh
tests/performance/dynamic/compare_state.py CPU/topology-state.dat \
  GPU/topology-state.dat --rtol 3e-4 --atol 1e-6
```

The diagnostic comparison found the first device difference at step 40 in a
transverse quantity near machine zero. Indices, flags, masses, charges, and
volumes remained bit-for-bit identical at every event. Running the OpenACC
kernel on the CPU and comparing it with the original CPU EOM gave differences
of only about `3e-16`, which rules out a material algorithm mismatch in the
accelerator EOM. Runtime isolation also showed that GPU EOM/curvature is the
main source; GPU Coulomb with CPU EOM diverges later and much less.

For diagnosis, `JETSPIN_OPENACC_DISABLE_EOM=1` and
`JETSPIN_OPENACC_DISABLE_COULOMB=1` independently force those components back
to their CPU implementations. `JETSPIN_OPENACC_DISABLE_PERSISTENT=1` restores
the call-scoped Test 13 control path. These switches are disabled by default.

Compare a new paired run with:

```sh
tests/performance/dynamic/compare.sh \
  CPU/run.log GPU/run.log CPU/statout.dat GPU/statout.dat
```
