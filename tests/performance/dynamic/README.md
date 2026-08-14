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

## Maxwell evaporation topology

Test 16 starts from 100 beads and combines Maxwell rheology, Yarin evaporation,
insertion, removal, and forced capacity growth. The stored input selects RK4;
variants changing only `integrator` to `1` or `2` validate Euler and RK2. The
standard OpenACC path keeps each integrator's one, two, or four stages and its
intermediate state on the device. The 1,000-step CPU and A30 acceptance totals
for all three integrators are 111 additions, 122 removals, two reallocations,
and 89 active beads.

With `NVCOMPILER_ACC_NOTIFY=2`, the normal build shows no jet-state, Coulomb,
force, stress, or derivative-array transfer between stages for Euler, RK2, or
RK4. Ordinary steps exchange only topology decision scalars. New/removed
records, selected statistical samples, capacity rebinds, and the final
checkpoint account for the remaining data traffic. The complete-force and
Coulomb-only oracle targets are diagnostic exceptions and deliberately copy
their required stage data.

Three-step CPU/GPU `statout.dat` and XYZ geometry comparisons for Euler, RK2,
and RK4 are identical before the first topology event. The CPU crosses the
first insertion threshold at step 4 and both standard GPU and oracle paths
at step 5. The bending instability then amplifies roundoff, so aggregate
topology, pre-event agreement, and the transfer audit are the acceptance
criteria.

## Kelvin--Voigt evaporation topology

Test 17 repeats the Test 16 workload with Kelvin--Voigt rheology. The stored
input selects RK4; variants changing only `integrator` to `1` or `2` validate
Euler and RK2. Every force stage, concentration-dependent evaporation stress,
direct Coulomb sum, state update, and statistics reduction remains on the
device. The accepted 1,000-step CPU and A30 totals for all three integrators
are 111 additions, 122 removals, two reallocations, and 89 active beads.

The standard transfer policy is the same as Test 16: no jet-state, force,
stress, or RK-derivative array moves between stages. Only topology decisions
cross every timestep; event records, output samples, capacity rebinds, and the
final checkpoint cause data movement. The development-only
`nvfortran-openacc-force-oracle` build transfers stage state and derivatives to
evaluate the trusted CPU equations for either rheology.
`nvfortran-openacc-coulomb-oracle` transfers only the state needed by the host
direct sum and uploads its force array. Both are numerical-isolation tools.

Three-step CPU/GPU trajectories for Euler and RK2 are identical at
`rtol=1e-12`, `atol=1e-13` before the first topology event, but their insertion
threshold occurs at step 4 on the CPU and step 5 on the GPU.
The bending instability then amplifies the roundoff difference, so topology
totals and the transfer audit are the acceptance criteria rather than a strict
pointwise trajectory comparison.

## Reproducible Maxwell/Kelvin--Voigt validation

Run all six deterministic dynamic-evaporation combinations with:

```sh
module use /opt/nvidia/hpc_sdk/modulefiles
module load nvhpc/24.3
tests/performance/dynamic/validate_evaporation.sh
```

The script builds isolated NVFORTRAN CPU and standard OpenACC executables,
then runs Tests 16 and 17 with Euler, RK2, and RK4. It checks the accepted
1,000-step topology totals and performs strict three-step pre-event CPU/GPU
comparisons. Test 16 additionally requires byte-identical XYZ geometry. The
standard GPU target is built without either diagnostic oracle macro.

Set `GPUCC` and `CUDA_VERSION` to select another NVIDIA target. Set
`JETSPIN_DYNAMIC_EVAP_KEEP=1` to retain build logs, inputs, and outputs in the
reported temporary directory.
