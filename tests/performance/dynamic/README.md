# Dynamic-topology benchmark

Test 13 initializes 1,024 beads across the nozzle--collector distance and uses
an intentionally high 2,000 cm/s injection speed to force topology activity
within 1,000 RK4 steps. The reference event sequence contains 13 insertions
and 13 removals, ends with 1,024 active beads, and never falls below 1,023.

The current milestone keeps insertion, removal, reallocation, and compaction
on the host. OpenACC Coulomb and eligible EOM calls use call-scoped mappings;
this is the numerical baseline for a later persistent dynamic-data port.

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
to their CPU implementations. These switches are disabled by default.

Compare a new paired run with:

```sh
tests/performance/dynamic/compare.sh \
  CPU/run.log GPU/run.log CPU/statout.dat GPU/statout.dat
```
