# OpenACC GPU porting status

JETSPIN's NVIDIA GPU port is being developed incrementally with explicit
OpenACC data regions. The normal GFortran, Intel, and MPI targets continue to
select the original CPU implementation.

## Implemented milestones

The `nvfortran-openacc` target currently accelerates the direct Coulomb
summation when all of these conditions hold:

- execution is serial (`mxrank == 1`);
- multiple-step Coulomb summation is disabled;
- evaporation is enabled for the serial 3D path;
- the system is either one- or three-dimensional.

The build defaults to `GPUCC=80` and `CUDA_VERSION=12.3`, producing
`-gpu=cc80,cuda12.3,nofma`. This covers the NVIDIA A30 and A100. `nofma`
prevents floating-point contraction in the sensitive three-point curvature
calculation. Select a different architecture without the `cc` prefix, for
example:

```sh
make -C source -f ../build/Makefile nvfortran-openacc \
  GPUCC=90 CUDA_VERSION=12.3
```

The GPU and OpenACC-host targets make the standard `_OPENACC` preprocessing
macro available. Accelerator imports, initialization, Coulomb dispatch, and
EOM dispatch are compiled only when this macro is present. Normal CPU and MPI
targets do not contain references to those accelerator entry points.

Unsupported configurations retain the existing CPU path. In particular, MPI
and the multiple-step/neighbour-list algorithm are not GPU kernels yet.

The accelerator kernel assigns one target bead to each parallel iteration.
That iteration visits every other active bead and writes only its target
force. This target-centric formulation is race-free and needs no atomics. It
evaluates each physical pair twice instead of sharing a pair contribution as
the CPU implementation does, so its floating-point accumulation order is
different.

The implementation preserves the existing model conventions, including:

- charge and mass scaling;
- the softening cross section associated with the higher-index bead;
- frozen-bead exclusion;
- the one-dimensional distance cutoff;
- mirror-charge contributions and their three-dimensional cutoff.

The three-dimensional equation-of-motion force assembly is offloaded for the
fixed Tests 9--12 and 20 configurations and for the bounded dynamic
evaporation paths in Tests 16 and 17. One kernel is launched per integrator
stage. Its local
three-point curvature calculation is device-side: an iteration reads the
current bead and its two neighbours. No host curvature array is built or
transferred. Disabling fused multiply-add preserves the accepted trajectory
for the initially straight geometry. Configurations outside the explicitly
validated gates use the original CPU EOM path.

## Explicit data region

Coordinates, bead properties, cross sections, and the Coulomb force array are
named explicitly in each OpenACC data region. No managed/unified-memory build
mode is used.

Configurations outside the explicitly validated persistent gates in Tests
9--13, 16, 17, 20, and 21 continue to use separate call-scoped Coulomb and EOM
data regions.

Test 13 now records the bounded persistent dynamic-topology milestone. It
preallocates 1,280 slots, keeps RK4 and force data resident as `inpjet` and
`npjet` change, and performs collector detection and clamping on the device.
The host receives one removal decision and, when needed, the removed record.
Nozzle insertion, including threshold checks, blocked-bead release, record
initialization, and `npjet` update, also runs on the device. The host receives
topology scalars and, on an actual insertion, downloads only the injected
mass and charge required by the host statistics. Test 15 additionally
exercises small-capacity teardown, host
reallocation, and device remapping. General device-side compaction remains
outside this path. The
persistent A30 execution retains all 26 events and the same final topology,
but GPU RK4 rounding moves insertion threshold crossings progressively
earlier. The call-scoped and persistent A30 streams are versioned separately;
`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` restores the former exactly.
An optional full-state snapshot at every topology event verified exact bead
metadata and properties. Component-isolation switches showed that the growing
trajectory difference originates primarily in device EOM/curvature arithmetic;
the OpenACC EOM run on the host matches the original CPU path near machine
precision. See the [Test 13 record](../examples/test-13.md).

Tests 16 and 17 complete the Maxwell and Kelvin--Voigt evaporation paths for
Euler, RK2, and RK4. Their one, two, or four force/stress evaluations,
intermediate states, final state update, topology operations, and capacity
rebinding use the same persistent-data strategy. CPU and A30 executions of all
three integrators and both rheologies report 111 additions, 122 removals, two
reallocations, and 89 active beads. Three-step pre-event CPU/GPU comparisons
are identical at `rtol=1e-12` and `atol=1e-13`; Maxwell XYZ geometry is also
byte-identical at written precision. The first insertion threshold occurs at
step 4 on the CPU and step 5 on the GPU, so aggregate topology and transfer
behaviour, rather than pointwise trajectory identity, form the dynamic
acceptance criteria.

Tests 9--12 use an explicit persistent-data path. The primary jet state,
static bead properties, Coulomb force, EOM derivatives, and integrator scratch
arrays are mapped once and remain resident across all timesteps.
Each Coulomb stage computes its cross sections on the device; EOM consumes the
device Coulomb force directly, and all integrator updates execute on the device.

For the fixed stochastic Platen Tests 12 and 20, the complete Gaussian history
is generated on the CPU before loop timing in a fixed
step/bead/component/draw order. Test 12 requires 6,006,000 doubles
(48,048,000 bytes); the 100-step Test 20 requires 600,600 doubles. The OpenACC
build transfers this history once during initialization and indexes it on the
device; the CPU path indexes the same layout. No random generation or noise
transfer occurs inside the measured loop.
The history is capped at 100,000,000 doubles (about 763 MiB). Longer fixed
runs wrap to its beginning, preserving CPU/GPU reproducibility while making
the noise periodic after the stored interval.
The host no longer receives stage intermediates or Coulomb forces.

The fixed-topology Maxwell Platen evaporation path used by Test 20 also keeps
its three drift evaluations, stochastic velocity update, Heun position,
volume and stress updates, and statistics on the device. Its state and
Gaussian history therefore follow the same persistent-data policy as Test 12.

Test 21 extends Maxwell Platen evaporation to insertion and dynamic refinement.
Its standard 24,048,000-double Gaussian history covers the reserved 500-bead
capacity and is uploaded once. Eligible refinement checks return only three
reduction scalars. At an accepted event, the complete active state is
downloaded once for host target-mesh preparation. Akima slopes, tangents,
cubic coefficients, and the 11 field interpolations then run on the GPU. The
resulting cross-section radius is then used, still on the GPU, to reconstruct
bead volume and evaporation volume, rescale each for reference-volume
conservation, and convert the interpolated mass/charge densities back to
per-bead quantities. Host code retains the data-dependent, rare (a few events
per run) bookkeeping instead: normalized target-mesh/anchor-mesh construction,
anchor state save/restore, and the final mesh assembly.
If the resulting mesh exceeds capacity, the old
topology, evaporation state, Platen scratch arrays, Coulomb workspace, and
Gaussian history mappings are released before their host allocations change.
Existing indexed Gaussian values are preserved during stride repacking; only
new capacity slots are generated. The resized state and history are rebound
once. An A30 transfer audit found no complete state transfer on an ordinary
timestep. The one full-state download at final shutdown is independent of
refinement.

Test 22 repeats this lifecycle three times in 15,800 steps. The native A30 and
NVFORTRAN CPU runs retain the same active counts across all events and finish
with 536 elements. The first two event steps are identical; the third occurs
at step 15,701 on the A30 and 15,719 on the CPU because of the target-centric
Coulomb accumulation order. Their 80-row statistics compare within 2.5
percent. The complete-force oracle reproduces all CPU event steps and compares
within `6e-7` relatively.

An A30 transfer audit records four complete state downloads: one for each of
the three accepted refinement events and one at final shutdown. It also records
exactly three Gaussian-history uploads, three topology rebinds, and three
evaporation-state rebinds. No complete state is downloaded during an ordinary
timestep or an unsuccessful threshold scan. Across the three Akima events,
coordinates are uploaded six times (source and target once per event), source
fields 33 times, and interpolated results are downloaded 33 times. These
kilobyte-scale payloads are confined to accepted events.

Test 23 combines this repeated lifecycle with collector removal. The Maxwell
Platen persistent-path eligibility includes `removing yes`; the established
device topology primitive advances the active lower bound and clears the
collected bead. Removal itself transfers only point/event data, while accepted
Akima events retain the event-only target-mesh/anchor-bookkeeping boundary
described above.
The A30 audit again finds four complete state downloads (three refinement
events plus shutdown) and three history/topology/evaporation rebinds. The
ordinary removal check returns one four-byte control scalar; each of the four
accepted removals downloads and clears only one bead.

The standard build's CPU/GPU `statout.dat` agreement (`rtol=3e-2`) holds for
79 of 81 rows; the last two diverge once the jet reaches the `x=12 cm`
collector and the two builds fork onto different removal schedules (CPU ten
removals/527 final active beads, GPU four/532). Rebuilding Test 23 with the
existing narrow `nvfortran-openacc-coulomb-oracle` target -- direct Coulomb
sum on the host, everything else including topology and the reconstruction
kernel above still on the device -- confirms direct-Coulomb summation order
as the dominant source: its third event lands at the CPU's own step (15718),
its first seven removals match the CPU step-for-step, and its final active
count (530) matches the complete-force oracle exactly. Agreement then holds
for 80 of 81 rows, leaving only the last row diverging. This is a diagnostic
confirmation of an already-documented sensitivity, not a change to the
standard build's accepted topology.

The Yarin evaporation rate, the evaporation-specific direct Coulomb force, and
the concentration-dependent constitutive updates are now available in the
serial 3D OpenACC path. Maxwell stress uses the same neighbour tangent and
relative velocity convention as the CPU implementation. For Test 16, the
Maxwell Euler, RK2, and RK4 paths execute their one, two, or four force
evaluations and all intermediate/final state updates on the device. Charge
smoothing/restoration, nozzle geometry, evaporation cross sections, and the
final statistics reduction are device-resident as well. Neither state nor
derivative arrays are transferred between stages.
Kelvin--Voigt stress uses the corresponding relative acceleration and the full
product-rule derivative of the concentration-dependent viscosity and modulus.
These stress kernels are explicit parallel regions; no per-bead host round
trip is added.

The per-step path-length and maximum-stress reductions are fused with each
integrator's final state-update kernel. Their scalar accumulators also remain device
resident. A small follow-up kernel preserves the CPU rule that the last bead
index wins when several beads share the maximum stress. For ordinary
statistical output, only the seven state values of the selected bead and four
scalar accumulators are downloaded. XYZ, PDB, binary trajectory, and periodic
restart events request a complete state explicitly. With the standard Test 9
input, the five scheduled samples transfer 420 bytes in total and the final
restart performs the only 56,056-byte full-state download. The accelerator
records the last synchronized timestep so coincident output and restart events
never duplicate a transfer. Test 13 uses its separate bounded dynamic transfer
policy described above.

## Numerical validation

Compare the OpenACC path with a CPU reference built by the same `nvfortran`
version. Different Fortran compilers may use different pseudo-random-number
sequences, so a GFortran trajectory is not a suitable direct reference for a
stochastic NVIDIA build.

The direct CPU algorithm accumulates one pair into both beads, whereas the
race-free accelerator algorithm accumulates independently by target bead.
Their results should therefore be compared numerically, not as binary files.
The initial development checks use the regression-suite criterion with
`rtol=1e-6` and `atol=1e-9`. Run the complete CPU/GPU comparison with:

```sh
tests/regression/run.sh openacc
```

The standard 1,000-step regression validation on an NVIDIA A30 passed all
eight normal regression cases; the written observables matched the NVFORTRAN
CPU results at their output precision. The separate 1,000-bead RK4, Euler,
RK2, and Platen trajectories also pass their A30 baselines with `rtol=1e-6` and
`atol=1e-9`; see the [benchmark index](../examples/README.md) for their records.

`nvfortran-openacc-host` is useful for checking the accelerated control path
and dynamic allocation without a GPU. It does not replace the real-device
regression above and cannot establish GPU performance.

For a complete force-stage oracle, use the development-only target:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-force-oracle GPUCC=80
```

This enables the single `JETSPIN_DEV_HOST_FORCE_ORACLE` interface. At every
supported Euler, RK2, RK4, or fixed-topology Platen force evaluation, with or
without evaporation, the current state is downloaded, the complete trusted
CPU force equations (including direct Coulomb) are evaluated, and the
derivatives are uploaded. Maxwell and Kelvin--Voigt deterministic evaporation,
non-evaporative Platen, and Maxwell evaporative Platen share this interface.
Integration updates and dynamic topology remain on the GPU. The target is a
correctness diagnostic only; it is unsuitable for performance measurements or
production runs.

To isolate only direct-Coulomb accumulation, use:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-coulomb-oracle GPUCC=80
```

This enables only `JETSPIN_DEV_HOST_COULOMB_ORACLE`. It downloads only the
state required by the trusted direct sum, evaluates Coulomb on the host at
every supported force evaluation, restores the nozzle charge when evaporation
is active, and uploads `ycf`; all other force terms and integration updates
stay on the device. It has the same meaning for Euler, RK2, RK4, and Platen,
for non-evaporative simulations and for the supported evaporation rheologies.
It is distinct from the complete force oracle and likewise adds per-stage PCIe
transfers. Additional preprocessor switches can be passed through
`FPPFLAGS_EXTRA`.

The unified oracle checks include the fixed 1,000-bead stochastic cases.
Tests 12 and 20 match their NVFORTRAN CPU statistical outputs exactly in all
six rows and fourteen columns with both the complete-force oracle and the
Coulomb-only oracle. These checks also cover the positive and negative Platen
predictors and the partial Heun position, evaporation, and stress evaluations.

The bounded dynamic-topology path is enabled for Maxwell Euler, RK2, and RK4
evaporation (Test 16), including insertion, removal, and capacity rebinds.
All three A30 runs reproduce the CPU event totals (111 insertions, 122 removals,
two reallocations, 89 active beads). Runtime transfer audits of the standard
build (no diagnostic macros) confirm that no state, force, stress, or
derivative array crosses the PCIe boundary between stages. Each timestep
exchanges only topology decision scalars. Insertion, removal, statistical
output, capacity growth, and the final checkpoint transfer only the records
required by those events. Full active-state transfers occur at the two
reallocations and at the final checkpoint. CPU/GPU trajectories are not
expected to be pointwise identical after the first topology threshold because
target-centric Coulomb accumulation and the subsequent bending instability
amplify floating-point ordering differences; topology totals, pre-event
agreement, and transfer behaviour are the acceptance criteria.

Test 17 applies the same transfer policy to all deterministic Kelvin--Voigt
evaporation integrators. Runtime transfer audits of Euler, RK2, and RK4 confirm
the same absence of per-stage array transfers. The historical CPU
Kelvin--Voigt evaporation EOM omits air drag and lift even when the input
enables air drag; the device implementation preserves that established
semantics.

Test 21 validates the hybrid dynamic-refinement path. The native A30 event is
accepted through anchor preservation, ordered path coordinates, topology, and
separate reference/evaporated-volume conservation rather than binary
trajectory identity. With NVFORTRAN 24.3, the CPU and complete-force-oracle
runs accept the same step and final topology; their final written statistics
agree within `3.6e-6` relatively. The native target-centric Coulomb path shifts
the accepted event by four steps, as expected for a bending-sensitive
trajectory.

Test 23 validates the Akima device kernels directly. A comparison build runs
the historical host spline and the accelerator spline for 11 fields at each
of three remeshes. The maximum coefficient relative difference is
`2.48e-15`; interpolated values differ by at most `7.28e-12` absolutely and
`3.32e-15` relatively. The calculation has no reduction: source slopes,
interior tangents, cubic coefficients, and target interpolations are
independent, with
only the constant-size endpoint extrapolation executed serially.

A second comparison build validates the bead volume/evaporation-volume
reconstruction, conservation rescale, and density-to-quantity conversion that
follow the Akima interpolation. It backs up the pre-reconstruction state,
evaluates the host reference, restores that state, then runs the device
kernel so it remains authoritative, and reports the maximum difference
against the host reference at each event. All three Test 23 events agree at
or near roundoff (worst absolute `7.1e-15`, worst relative `4.05e-16`), and
the standard build reproduces the same event/removal/final-active topology as
before this kernel existed. Building this oracle exposed two real ordering
bugs — the device kernel invoked twice, and the host reference evaluated from
already-converted state — each of which corrupted the event's mass/charge
invariant by tens of percent before being fixed; this is the reason every new
device stage in this project must be validated with an A/B oracle rather than
trusted from inspection alone.

## Next porting stages

1. Investigate packing the maximum stress and bead index into one deterministic
   reduction so that its two follow-up kernels can also be removed.
2. Extend the refinement-capacity lifecycle beyond the currently validated
   serial Maxwell/Platen combination when additional GPU model combinations
   are enabled; current transfers occur only on resize events.
3. The remaining host-side work at an accepted event is the data-dependent
   mass-boundary walk, target-mesh/anchor-mesh construction, and anchor
   save/restore/final-assembly bookkeeping. It is a small, rare (a few events
   per run), inherently sequential scan rather than a per-bead parallel
   operation, so it is not a priority target; the next candidate is instead
   determining whether the accepted-event full-state round trip can be
   removed without duplicating that bookkeeping's model logic.
4. Evaluate one-GPU-per-rank MPI execution only after the single-GPU numerical
   path is stable.

The intended steady state is a persistent device-resident simulation with
host transfers for output, checkpoints, and topology changes—not a transfer
of the complete jet at every timestep.
