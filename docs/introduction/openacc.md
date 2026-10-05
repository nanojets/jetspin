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

A three-dimensional run that has not engaged a persistent device path
offloads the Coulomb sum only while the jet has at least 128 active beads
(`JETSPIN_OPENACC_COULOMB_MIN_BEADS`; 0 offloads every call, as before
2026-10-05). Such a call copies the jet to the device and the forces back,
which costs more than the host sum on a shorter jet; the
[Test 25 time budget](../examples/test-25.md#phase-1-the-coulomb-sum-on-the-host)
gives the measured crossover. A persistent run always sums on the device.

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
different. The evaporative 3D kernel (`accelerator_coulomb_evap_3d`)
additionally spreads the sources of each target over the vector lanes of one
gang and combines them with a reduction, for the pair and the mirror terms
alike; the previous one-thread-per-target layout left the device almost idle
for jets of a few hundred beads.

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

## Dynamic-topology Platen fork

The persistent gates listed above all require a bead count that a realistic
run does not have when the decision is taken. `prepare_integrator_random_history`
runs once, before the timestep loop, and every eligibility test it consults
requires `npjet >= 100`; a jet started from a single bead therefore never
allocates its Gaussian history, and since `allocated(gaussianhistory)` is a
required condition in the activation gate, the persistent path can never
engage for the rest of the run however large the jet grows. This was confirmed
by profiling: on `examples/input-24`, the only kernel `nsys` recorded on the
device was the Coulomb summation, with the integrator running on the host.

The `nvfortran-openacc-dynamic-platen` target (see
[compiling](compiling.md)) addresses this for the non-evaporative Platen
integrator only. Its eligibility logic lives in
[`openacc_dynamic_platen_mod.f90`](../../source/openacc_dynamic_platen_mod.f90),
which splits the decision in two: `dynamic_platen_accelerator_configured` tests
only size-independent model configuration and is what the one-shot history
allocation consults, while `dynamic_platen_accelerator_eligible` adds the
`npjet >= 100` and `mxnpjet > npjet` thresholds and gates the per-step
activation. The history is the sequential pool described in
[random numbers](random-numbers.md), the default layout for every build since
2026-09-30, which removes the remap that capacity growth used to force.

With the macro on, a step in the persistent branch runs entirely on the
device: charge smoothing and the Coulomb/electric driver four times, four EOM
stages, then the Platen predictor, velocity, position, stress-statistics, and
non-inserted-position kernels, with no per-step host transfer. The host is
reached only for topology events and output.

A single-bead run using `input-24` physics without evaporation engaged at step
1,379,105 and passed four capacity-growth, reset, and re-engage cycles at
`npjet` 223, 353, 496, and 576, completing its full six million steps with
`Program closed correctly`, no NaN, and no device errors. The first two of
those cycles are exactly where earlier builds failed, before the Coulomb
mapping reset above was corrected. The jet bridges the collector distance and
finishes with its bead count oscillating around 570 rather than drifting, so
the balanced insertion/removal regime is reached on the device from a
single-bead start. That run predates the nozzle-insertion fix of 2026-10-01
described below and is not a reference for the physics: an NVFORTRAN CPU run
of the same input has about 512 active beads between steps 4 and 5 million.

The same allocation defect affected `platen_ev`: a 400-bead start engaged at
step 1 while a single-bead start never engaged, its history still unallocated
once the jet had grown to 262 beads. Since 2026-09-30 the one-shot decision
uses the size-independent `dynamic_evaporative_platen_configured`, and the
per-step gate keeps `dynamic_evaporative_platen_eligible` with its size
thresholds: a single-bead Test 25 engages the persistent path when it first
exceeds 100 beads (step 1,446,413, 120 beads, in the reference run). The
decision does not depend on `JETSPIN_OPENACC_DISABLE_PERSISTENT`, so CPU,
persistent and non-persistent runs read the same noise. `rk4sys_ev` has not
been examined for the same defect.

## Persistent Coulomb mapping reset

A persistent run maps the Coulomb force array to the device once and keeps it
there. When a topology change forces the host array to be reallocated, that
mapping has to be torn down first, or the device retains a pointer to storage
the host no longer owns.

`reset_coulomb_accelerator` used to clear its bookkeeping flags without ever
issuing the matching `exit data delete`, so the mapping was orphaned rather
than removed. The consequence depends on where the host allocator places the
new array: if it lands elsewhere, the next `enter data` fails with a
*partially present* error; if it reuses the same address, OpenACC finds the
stale entry still valid and silently reuses the old device buffer. The second
case is the dangerous one, because it produces wrong numbers with no
diagnostic — a jet that freezes at a fixed bead count, or NaN a few hundred
steps later.

The routine now takes the force array as an optional argument and, when it is
supplied, deletes `coulcrossec` and the array before clearing its flags. The
one caller that omitted it, in `rk4sys`, was corrected to pass it. This is a
defect of the shipped `nvfortran-openacc` target, not of any development fork:
`examples/input-15` on the previous commit produced NaN from its first printed
line with the bead count frozen at 103, against `x = 35.99` and 211 beads on
CPU; after the fix the GPU result matches the CPU one to eight significant
digits. `examples/input-13` and `input-14` hide it at their shipped 1,000-step
length, which performs no reallocation at all, and reproduce it when extended
to 25,000 steps.

Any future device-resident array must therefore pair its `enter data` with an
`exit data delete` on every path that can reallocate the host storage.
Bookkeeping flags alone do not unmap anything.

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
Nozzle insertion, including threshold checks, blocked-bead release, record
initialization, and `npjet` update, also runs on the device, and so do the
placement and the charge smoothing of the blocked bead before every force
evaluation. One serial kernel decides release, insertion, capacity
exhaustion, and removal; the outcome returns to the host as one 20-byte
record per step. On an actual insertion the host downloads only the injected
mass and charge required by the host statistics, and on a removal only the
removed record. Test 15 additionally exercises small-capacity teardown, host
reallocation, and device remapping. General device-side compaction remains
outside this path.

The persistent A30 run reproduces all 26 CPU topology events step for step.
Until 2026-10-01 its insertions fell progressively earlier (step 593 instead
of 600 for the seventh): the device placed the blocked bead only when it was
created and never smoothed its charge, so the device Coulomb sum saw a full
charge at the nozzle. That drift had been attributed to GPU RK4 rounding.
`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` still restores the call-scoped path.
An optional full-state snapshot at every topology event verifies bead
metadata and properties. See the [Test 13 record](../examples/test-13.md).

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

For the fixed stochastic Platen Tests 12 and 20, the Gaussian pool is
generated on the CPU before loop timing, each step's slice in the draw order
of the per-step block. Test 12 requires 6,006,000 doubles (48,048,000 bytes);
the 100-step Test 20 requires 600,600 doubles. The OpenACC build transfers the
pool once during initialization and indexes it on the device; the CPU path
reads the same values. No random generation or noise transfer occurs inside
the measured loop. The pool holds at most `noise pool` values (default
100,000,000, about 763 MiB); a run that consumes it wraps to its beginning,
preserving CPU/GPU reproducibility while making the noise periodic (see
[random numbers](random-numbers.md#the-noise-repeats-after-the-pool-is-consumed)).
The host no longer receives stage intermediates or Coulomb forces.

The fixed-topology Maxwell Platen evaporation path used by Test 20 also keeps
its three drift evaluations, stochastic velocity update, Heun position,
volume and stress updates, and statistics on the device. Its state and
Gaussian history therefore follow the same persistent-data policy as Test 12.

Test 21 extends Maxwell Platen evaporation to insertion and dynamic refinement.
Its Gaussian pool (100,000,000 doubles by default) is uploaded once and is
not affected by capacity growth. Eligible refinement checks return only their
three reduction results, in one transfer. At an accepted event, the complete active state is
downloaded once for host target-mesh preparation. Akima slopes, tangents,
cubic coefficients, and the 11 field interpolations then run on the GPU. The
resulting cross-section radius is then used, still on the GPU, to reconstruct
bead volume and evaporation volume, rescale each for reference-volume
conservation, and convert the interpolated mass/charge densities back to
per-bead quantities. Host code retains the data-dependent, rare (a few events
per run) bookkeeping instead: normalized target-mesh/anchor-mesh construction,
anchor state save/restore, and the final mesh assembly.
If the resulting mesh exceeds capacity, the old
topology, evaporation state, Platen scratch arrays, and Coulomb workspace
mappings are released before their host allocations change, and the resized
state is rebound once. The Gaussian pool stays mapped unchanged. An A30 transfer audit found no complete state transfer on an ordinary
timestep. The one full-state download at final shutdown is independent of
refinement.

Test 22 repeats this lifecycle three times in 15,800 steps. The native A30
run and the complete-force oracle reproduce all CPU event steps, topologies,
and capacities, and all 80 statistics rows agree with the NVFORTRAN CPU run
within `8.0e-10` and `3.2e-8` relatively. CPU and GPU read the same sequential
Gaussian pool; what remains is the different floating-point evaluation order
of the device kernels. Before the persistent-path fixes described in the
Test 25 section below, the native build moved its fourth nozzle insertion one
step earlier, which offset the pool position, and every later event
differed.

An A30 transfer audit records four complete state downloads: one for each of
the three accepted refinement events and one at final shutdown. It also records
three topology rebinds and three evaporation-state rebinds; the Gaussian pool
is uploaded once at startup and never transferred again. No complete state is downloaded during an ordinary
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
events plus shutdown), three topology/evaporation rebinds, and a single
Gaussian-pool upload at startup. Each
ordinary step returns one 20-byte topology record, which carries the insertion
and removal decisions together (until 2026-10-01, two uploads and five
downloads of separate scalars); each accepted
removal downloads and clears only one bead.

The standard build reproduces the CPU events (14,268, 14,971, 15,691), the
nine removals, and the 526 final active beads step for step, and all 81
statistics rows agree with the CPU within `8.0e-10` relatively, including
the rows written after the leading bead has reached the collector. The
complete-force oracle agrees within `6.1e-8`. Rebuilding Test 23 with the
narrow `nvfortran-openacc-coulomb-oracle` target -- direct Coulomb sum on the
host, everything else including topology and the reconstruction kernel above
still on the device -- gives a statistics file byte-identical to the native
run's: at the printed precision the device Coulomb kernel adds nothing to the
remaining roundoff-level difference.

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
returns only one 20-byte topology record. Insertion, removal, statistical
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
trajectory identity. With NVFORTRAN 24.3, the CPU, complete-force-oracle, and
native A30 runs accept the same step (7,131) and final topology; the written
statistics of the oracle and of the native run agree with the CPU within
`1.4e-7` and `2.4e-10` relatively.

Test 23 validates the Akima device kernels directly. A comparison build runs
the historical host spline and the accelerator spline for 11 fields at each
of three remeshes. The maximum coefficient relative difference is
`1.03e-14`; interpolated values differ by at most `1.82e-12` absolutely and
`1.89e-15` relatively. The calculation has no reduction: source slopes,
interior tangents, cubic coefficients, and target interpolations are
independent, with
only the constant-size endpoint extrapolation executed serially.

A second comparison build validates the bead volume/evaporation-volume
reconstruction, conservation rescale, and density-to-quantity conversion that
follow the Akima interpolation. It backs up the pre-reconstruction state,
evaluates the host reference, restores that state, then runs the device
kernel so it remains authoritative, and reports the maximum difference
against the host reference at each event. All three Test 23 events agree at
or near roundoff (worst absolute `8.9e-16`, worst relative `1.4e-16`), and
the standard build reproduces the same event/removal/final-active topology as
before this kernel existed. Building this oracle exposed two real ordering
bugs — the device kernel invoked twice, and the host reference evaluated from
already-converted state — each of which corrupted the event's mass/charge
invariant by tens of percent before being fixed; this is the reason every new
device stage in this project must be validated with an A/B oracle rather than
trusted from inspection alone.

## Dynamic evaporative Platen from a single bead (Test 25)

Test 25 grows the jet from one bead, so the device sees tens to a few hundred
beads. Per-step wall times (NVFORTRAN 24.3, one A30, 2.5 million steps,
median per bead-count band, from the `cpu` statistics key):

| Active beads | CPU | GPU, old Coulomb kernel | GPU, current kernel | Current kernel, persistent path |
| --- | ---: | ---: | ---: | ---: |
| 1-30 | 33 us | 281 us | 253 us | 248 us |
| 60-100 | 280 us | 503 us | 407 us | 390 us |
| 100-150 | 542 us | 641 us | 490 us | 458 us |
| 150-200 | 849 us | 799 us | 568 us | 454 us |
| 250-300 | 1899 us | 1155 us | 784 us | 447 us |
| Whole run | 1522 s | 1557 s | 1188 s | 931 s |

These columns predate the fixes of 2026-09-30 and 2026-10-01. Current
numbers, measured with each process bound to its GPU's NUMA node, are in the
[Test 25 time budget](../examples/test-25.md#where-the-a30-run-spends-its-time),
which also lists what runs on the CPU and on the GPU before and after the
persistent path engages.

The persistent column was measured with a scratch version of the fix that is
now in the code (see the dynamic-topology Platen fork above). Before it, the
Gaussian history of `platen_ev` was allocated only if the jet already had 100
beads when `prepare_integrator_random_history` ran, which is never true for a
single-bead start, so the persistent path never engaged and only the Coulomb
kernel ran on the device. The path now engages at step 1,446,413 (120 beads),
right after an accepted refinement event. Over the 5 million steps of the
Test 25 reference the A30 takes about 1570 s and the NVFORTRAN CPU build
5893 s.

Nsight Systems profiles of the current build explain the remaining cost:

- non-persistent path (up to 95 beads in Test 25): three Coulomb calls per
  Platen step, each with 13 host-to-device copies, 2 device-to-host copies,
  6 stream synchronizations, and 2 kernels. A call cost 78 us of wall time
  (`JETSPIN_PROFILE=1`), against 18 us for the same sum on the host, and the
  device was busy for about 35 us of it. Since 2026-10-05 the sum stays on
  the host below 128 active beads, and phase 1 of Test 25 takes 156 s
  instead of 430 s;
- persistent path (about 215 beads): 31 kernel launches, 59 stream
  synchronizations, one 20-byte device-to-host copy and no upload per step.
  The kernels take 220 us of the 399 us step, 66 us of them in the three
  Coulomb sums (at about 273 beads, with collector removal: 34 launches,
  65 synchronizations, 236 us of kernels in a 402 us step). The step time is flat at about 400 us from 120 to 300 beads,
  so it is set by launch and synchronization overhead rather than by the
  bead count. Before the Coulomb kernel rewrite of 2026-09-30 the three
  Coulomb calls took 88 us each at about 170 beads, 72 % of the device time.

Below about 100 beads the CPU build is also faster than the persistent path;
there the A30 run now does the same work on the host.

### Persistent-path defects found by a seed ensemble (2026-09-30)

Once the size-independent history decision let the persistent path engage,
a seed ensemble of Test 25 (five NVFORTRAN CPU seeds, three A30 seeds, all
with the sequential Gaussian pool) showed a systematic difference. Up to the
engagement step all runs agree; afterwards the A30 jet grows longer and
denser: between 4 and 5 million steps the mean active count was 297 on every
A30 seed and 268 on the CPU, and the path length 123 cm against 112 cm, with
a seed-to-seed spread below 0.5 percent on both builds. Three defects of the
persistent evaporative Platen path were responsible:

- **Nozzle charge smoothing.** While a newly inserted bead is still being
  released at the nozzle, its charge is scaled by a smooth cut-off before
  every Coulomb evaluation and restored afterwards. `smooth_charge` and
  `restore_charge` dispatched to their device versions only for
  `systype 3`; a `systype 4` (stochastic Maxwell) persistent run smoothed
  the stale host copy, and the device Coulomb sum saw the bead with its full
  charge. This caused the difference before the jet reaches the collector.
  With the dispatch extended to `systype 4`, the A30 run of seed 317 follows
  the CPU trajectory through 2 million steps: the same bead count at every
  sample and a path-length difference of at most `3e-8` relatively, the same
  as before engagement, against `5e-4` without the fix at 1.54 million
  steps.
- **Lead-bead curvature after removal.** Once beads have been collected,
  `eom3_ev` and `eom4_ev` give the lead bead the surface-tension and lift
  terms computed with the last collected bead. The Maxwell evaporative device
  stage did not request them (`collector_curvature`), unlike the
  Kelvin--Voigt stage.
- **Frozen and inserting beads in the Platen kernels.** The predictor and
  velocity kernels gave the stochastic force to beads frozen at the
  collector and to the bead still being inserted, and the position kernel
  moved frozen beads with their velocity; `eom4_ev` and `eom4_pos_ev` keep
  them fixed. This is why the Test 23 complete-force oracle, which keeps these
  kernels on the device, differed from the CPU by 1.4 percent once the lead
  bead had reached the collector; with the fix all 81 rows agree within
  `6.1e-8`.

The first defect was invisible to the oracles, which smooth the charge on
the host, and to Tests 16 and 17, which use `systype 3`. Tests 21--23 run
with a reduced charge density and were accepted on event contracts, so the
native difference was attributed to the Coulomb summation order. With the
fixes their native runs reproduce the CPU events, capacities, and removals,
and every statistics row agrees within `8.0e-10`. Over the full 5 million
steps of Test 25 the fixed A30 run gives the CPU stationary values (268 active
beads, 112 cm path length).

The non-evaporative dynamic device paths had the same omissions, fixed on
2026-10-01. The `systype 3` dynamic RK4 path (Tests 13--15) also placed the
blocked nozzle bead only on the stale host arrays, so on the device it stayed
where it had been created; the persistent non-evaporative Platen branch
smoothed the charge once per step on the host copy and placed the bead only
at the end of the step. Placement (`place_inserting_bead`), smoothing, and
restoring now bracket every force evaluation on the device, the
collector-curvature terms are requested, and the non-evaporative Platen
kernels treat frozen and inserting beads as `eom4` does. Tests 13 and 14
now reproduce the CPU topology stream step for step (Test 13 used to insert
progressively earlier, step 593 instead of 600 for the seventh insertion).
For the dynamic Platen fork the reference is the same build with
`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` (host integration, same Gaussian
pool): on Test 24 physics from a single bead, 2 million steps, the fixed fork
follows it within `1.2e-8` in path length through 1.8 million steps and
`7e-7` at 2 million, with the same bead count in all 100 frames and 86 of 88
insertion steps identical (the last two one step apart). The unfixed fork
left it immediately after engaging at step 1,379,105 (`5e-4` at 1.48 million
steps, first bead-count difference at 1.42 million) and ended with 92
insertions and 282 active beads against 88 and 276.

## Next porting stages

1. For the dynamic evaporative Platen path: run the persistent step on one
   asynchronous queue with a single synchronization per step. An ordinary
   step now returns only two small records (topology decisions, and the
   refinement scan when it runs); fuse each force stage into fewer kernels
   afterwards if needed.
2. Investigate packing the maximum stress and bead index into one deterministic
   reduction so that its two follow-up kernels can also be removed.
3. Extend the refinement-capacity lifecycle beyond the currently validated
   serial Maxwell/Platen combination when additional GPU model combinations
   are enabled; current transfers occur only on resize events.
4. The remaining host-side work at an accepted event is the data-dependent
   mass-boundary walk, target-mesh/anchor-mesh construction, and anchor
   save/restore/final-assembly bookkeeping. It is a small, rare (a few events
   per run), inherently sequential scan rather than a per-bead parallel
   operation, so it is not a priority target; the next candidate is instead
   determining whether the accepted-event full-state round trip can be
   removed without duplicating that bookkeeping's model logic.
5. Evaluate one-GPU-per-rank MPI execution only after the single-GPU numerical
   path is stable.

The intended steady state is a persistent device-resident simulation with
host transfers for output, checkpoints, and topology changes—not a transfer
of the complete jet at every timestep.
