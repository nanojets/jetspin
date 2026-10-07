# OpenACC GPU porting status

JETSPIN's NVIDIA GPU port is being developed incrementally with explicit
OpenACC data regions. The normal GFortran, Intel, and MPI targets continue to
select the original CPU implementation.

## Device strategy and its coverage

Tests 24 and 25 set the strategy that every run is to follow (work plan in
`docs/STATE.md`, from 2026-10-06):

1. a size-independent configuration test decides once, before the loop,
   whether the run can use the device path (and whether it reads the
   Gaussian pool); a per-step gate opens the device path when `npjet`, the
   index of the nozzle bead, reaches 100 and keeps it open (`npjet` counts
   the beads in the arrays, collected ones included until a compaction
   removes them, not the active ones);
2. below the gate the OpenACC build runs the CPU build's code, Coulomb sums
   included, and reproduces the CPU build byte for byte;
3. above it the state stays on the device and an ordinary step returns one
   20-byte topology record;
4. one asynchronous queue, one wait per step, small kernels fused;
5. stochastic forcing from the sequential Gaussian pool;
6. exact restart.

Euler, RK2, RK4 and Platen runs follow it through one device step
(`source/device_step_mod.f90`, milestones M1-M3): one set of model
conditions (`device_step_supported`, which for the Platen scheme also
decides the Gaussian pool, in every build and for any number of ranks), one
configuration test that adds the data layout (`device_step_configured`),
one sticky gate at `npjet` = 100 (`device_step_eligible`), one engagement
routine that maps the state and prints `OpenACC device step engaged
(scheme, model) at step N with M active beads` (`resumed`, after a restart
from a state on which the gate had opened: the restart file records it,
restart state version 2, since a compaction can leave fewer than 100 beads
in the arrays of an engaged run), and one stage evaluator for
the Maxwell and Kelvin-Voigt models, each with or without evaporation, with
or without air drag. The RK schemes combine its stages with their update
kernels (`device_rk_step`); the Platen scheme (`device_platen_step`) adds the
stochastic terms, one fused preparation kernel per force evaluation and the
two-kernel tail (`accelerator_platen_update`, `accelerator_platen_end_step`),
with or without insertion. The integrators of `integrator_mod` and
`integrator_kv_ev_mod` contain only the CPU build's code and dispatch to it
above the gate. Still outside it (`device_step_supported`,
`device_step_configured`): multiple-step Coulomb, refinement and the other
options that tag beads with the RK schemes, time-dependent external fields,
Lorentz force, upper potential, drag velocity, breakup, bead tracking
(`print binary ... style 4`), a fixed bead set with `removing yes`, 1D
systems and MPI runs; they run the CPU code in the OpenACC build. The RK
schemes take the device step on system 3, the Platen scheme on system 4.

Coverage measured on 2026-10-06 (A30, NVHPC 24.3, `NVCOMPILER_ACC_NOTIFY`;
launches and transfers per step over the test input, or over the first 200
steps of the examples):

| Case | Model | Device path | Launches / step | Transfers / step |
| --- | --- | --- | ---: | ---: |
| Ex 1, 2 | 1D, RK2 / RK4 | host, identical to the CPU build (the 1D Coulomb sum ran on the device at every step until M1) | 0 | 0 |
| Ex 3, 5, 6, 7 | 3D, from one bead | host below 100 beads, identical to the CPU build; Example 3 takes the device step above it, Examples 5 and 6 need multiple-step Coulomb, Example 7 a time-dependent field | 0 | 0 |
| Ex 4 | Platen, no refinement | host below 100 beads, identical to the CPU build, both reading the Gaussian pool since M3 (regression baselines regenerated); the device step above it | 0 | 0 |
| Ex 8 | RK4 with evaporation, no air drag | host, then the device step from step 13,699 (`npjet` 100, 89 active beads); identical to the CPU build up to it | 0 (first 200 steps) | 0 |
| Test 9 / 10 / 11 | RK4 / Euler / RK2, fixed 1000 beads | device step from step 1; one more launch since 2026-10-07 (freezing at the collector) | 22 / 10 / 14 | 0.2 |
| Test 12 | Platen, fixed 1000 beads | device step from step 1, fused tail (18 launches before M3), freezing at the collector (since 2026-10-07) | 13 | 0.2 |
| Tests 13-15 | RK4 with insertion | device step from step 1 | 28-30 | 1.5-1.6 |
| Tests 16, 18 | RK4 with evaporation | device step from step 1 | 29 / 30 | 4.3 / 22 |
| Tests 17, 19 | RK4 Kelvin-Voigt with evaporation | device step from step 1 | 33 / 34 | 4.3 / 4.2 |
| Test 20 | Platen with evaporation, fixed 1000 beads | device step from step 1, fused tail (19 launches before M3), freezing at the collector (since 2026-10-07) | 13 | 1.8 |
| Tests 21-25 | Platen with refinement | device step, asynchronous queue | 12-14 | about 2 |

The RK steps and the fixed-geometry Platen step are still synchronous, and
the RK steps are not fused (one launch per kernel): milestone M4.

Since 2026-10-07 a bead whose insertion needs larger arrays is inserted in
the step that grows them, as in the CPU build. Before, the device topology
path grew the arrays and inserted the bead one step later. A jet started
from one bead fills its initial arrays (100 beads) just when the gate opens,
so such runs left the CPU trajectory at their first insertion on the device
(Examples 3, 4, 6, 8), the others at their first reallocation (Tests 15-19,
24); now they separate only through roundoff, much later (`docs/STATE.md`,
item 3 of the pre-merge checks).

## Implemented milestones

Inside the device step the whole step runs on the GPU, the Coulomb sums
included. Outside it the `nvfortran-openacc` target runs the CPU build's
code and offloads only the direct Coulomb summation, when all of these
conditions hold:

- multiple-step Coulomb summation is disabled;
- the sum is the non-evaporative one-dimensional one (since 2026-10-06) or a
  three-dimensional one, with or without evaporation;
- the jet has at least 128 active beads
  (`JETSPIN_OPENACC_COULOMB_MIN_BEADS`; 0 offloads every call, as before
  2026-10-05).

Such a call copies the jet to the device and the forces back, which costs
more than the host sum on a shorter jet; the
[Test 25 time budget](../examples/test-25.md#phase-1-the-coulomb-sum-on-the-host)
gives the measured crossover.

Runtime switches (all of them, with the developer ones, in
[Compiling](compiling.md#runtime-environment-variables)):

- `JETSPIN_OPENACC_DISABLE_PERSISTENT=1` keeps the device step closed for
  every run, fixed bead sets included: the OpenACC build runs the CPU
  build's code, offloading only the Coulomb sums above; the Gaussian-pool
  decision does not change.
- `JETSPIN_OPENACC_DISABLE_EOM=1` does the same since 2026-10-07. Before, it
  made the device stage skip the equations of motion and integrate stale
  derivatives. A device EOM kernel that refuses a stage now stops the run
  with error 22 (`the OpenACC device step cannot evaluate the equations of
  motion of this run`); `device_step_supported` admits only what the
  kernels evaluate, so this is a guard, not an expected outcome.
- `JETSPIN_OPENACC_DISABLE_COULOMB=1` keeps the one-dimensional sum on the
  host; it has never affected the three-dimensional sums.
- `JETSPIN_OPENACC_SYNC=1` runs the Platen device step of a run with
  insertion synchronously; the RK device steps and fixed bead sets are
  always synchronous (milestone M4).

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
device-step dispatch are compiled only when this macro is present. Normal CPU and MPI
targets do not contain references to those accelerator entry points.

Unsupported configurations retain the existing CPU path. In particular, MPI
and the multiple-step/neighbour-list algorithm are not GPU kernels yet.

The accelerator kernel assigns one target bead to each parallel iteration.
That iteration visits every other active bead and writes only its target
force. This target-centric formulation is race-free and needs no atomics. It
evaluates each physical pair twice instead of sharing a pair contribution as
the CPU implementation does, so its floating-point accumulation order is
different. The three-dimensional kernels, evaporative
(`accelerator_coulomb_evap_3d`, since 2026-09-30) and non-evaporative
(`compute_coulomelec_openacc_3d`, since 2026-10-05), additionally spread the
sources of each target over the vector lanes of one gang and combine them with
a reduction, for the pair and the mirror terms alike; the previous
one-thread-per-target layout left the device almost idle for jets of a few
hundred beads. The non-evaporative kernel also takes its arrays as
assumed-shape dummies; before, a data region around allocatable dummies made
the runtime upload nine array descriptors, and wait, at every call. With both
changes the fixed 1,000-bead Tests 9--12 and the dynamic Tests 13 and 14 run
4.2 to 4.8 times faster on an A30 and still pass their recorded
comparisons.

The implementation preserves the existing model conventions, including:

- charge and mass scaling;
- the softening cross section associated with the higher-index bead;
- frozen-bead exclusion;
- the one-dimensional distance cutoff;
- mirror-charge contributions and their three-dimensional cutoff.

The three-dimensional equation-of-motion force assembly runs on the device
only inside the device step: one kernel per force evaluation
(`accelerator_eom3_stage`, which for an evaporating Maxwell jet also
computes the stress and the evaporation rate), followed for the
Kelvin-Voigt model by its stress kernel (`accelerator_kv_stress_3d`). Its local
three-point curvature calculation is device-side: an iteration reads the
current bead and its two neighbours. No host curvature array is built or
transferred. Disabling fused multiply-add preserves the accepted trajectory
for the initially straight geometry. Outside the device step the CPU
build's EOM code runs. Until milestones M1-M3 (2026-10-06) separate gates
shaped on the tests offloaded it (fixed Tests 9--12 and 20, dynamic Tests
13--25), and the non-evaporative Euler, RK2 and RK4 host paths offloaded
the EOM stage whenever air drag was on.

## Non-evaporative dynamic Platen (Test 24)

Since 2026-10-06 (milestone M3) the gates and persistent branches described
in this section and the next ones (`dynamic_platen_configured`,
`dynamic_platen_eligible`, `dynamic_evaporative_platen_configured`,
`dynamic_evaporative_platen_eligible`, `fixed_accelerator_geometry`,
`fixed_evaporative_platen_eligible`) are replaced by the device step above;
the step itself is unchanged (Tests 24 and 25 byte-identical over 2.1
million steps), and it now also serves runs without refinement (Example 4)
and the fixed bead sets of Tests 12 and 20, which lost their separate,
unfused kernels. The text below is the history of the Platen paths.

Since 2026-10-05 the non-evaporative stochastic Platen run growing from a
single bead (Test 24) follows the evaporative one (Test 25) in every build.
`dynamic_platen_configured` (size-independent) decides once, before the loop,
that the run reads the sequential Gaussian pool, on the CPU as on the GPU;
`dynamic_platen_eligible` adds `npjet >= 100`, `mxnpjet > npjet`, and
`JETSPIN_OPENACC_DISABLE_PERSISTENT` and opens the persistent device path in
`platen()`. The step then runs on one asynchronous queue with the same fused
small kernels as the evaporative step. Until then this path existed only as
the `nvfortran-openacc-dynamic-platen` build option (macro
`JETSPIN_GPU_DYNAMIC_PLATEN`, module `openacc_dynamic_platen_mod.f90`, both
removed); the standard builds drew Test 24's noise step by step and
integrated it on the host with only the Coulomb sums on the GPU. The history
below describes that option.

The persistent gates listed above all require a bead count that a realistic
run does not have when the decision is taken. `prepare_integrator_random_history`
runs once, before the timestep loop, and every eligibility test it consults
requires `npjet >= 100`; a jet started from a single bead therefore never
allocates its Gaussian history, and since `allocated(gaussianhistory)` is a
required condition in the activation gate, the persistent path can never
engage for the rest of the run however large the jet grows. This was confirmed
by profiling: on `examples/input-24`, the only kernel `nsys` recorded on the
device was the Coulomb summation, with the integrator running on the host.

The `nvfortran-openacc-dynamic-platen` target addressed this for the non-
evaporative Platen integrator only. Its eligibility logic, now in
`integrator_mod.f90`, splits the decision in two: `dynamic_platen_configured`
tests only size-independent model configuration and is what the one-shot
history allocation consults, while `dynamic_platen_eligible` adds the `npjet
>= 100` and `mxnpjet > npjet` thresholds and gates the per-step activation.
The history is the sequential pool described in
[random numbers](random-numbers.md), the default layout for every build since 2026-09-30, which
removes the remap that capacity growth used to force.

A step in the persistent branch runs entirely on the device: three charge
smoothings, Coulomb sums, and EOM stages, the Platen predictor, velocity, and
position updates, the stress derivative at the new state, and the stress
statistics, with no per-step host transfer; since 2026-10-05 the small
kernels are fused as in the evaporative step. Since 2026-10-06 the velocity
and position updates share one kernel (`accelerator_platen_update`, also used
by the evaporative step), and the stress derivative at the new state is
evaluated inside the end-of-step kernel: a fourth EOM stage computed all the
derivatives for that one alone. The host is reached only for
topology events and output.

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
persistent and non-persistent runs read the same noise.

The evaporative Euler, RK2 and RK4 integrators (`eulsys_ev`, `rk2sys_ev`,
`rk4sys_ev`, system 3) do not share the defect: they draw no Gaussian pool,
and their gate, `evaporative_dynamic_accelerator_eligible`, had no size
threshold at all, so a jet grown from one bead ran on the device from the
first step, at about 430 us per step with a few beads against 7 us on the
host (checked on 2026-10-06). Since then the gate requires 100 beads, as the
Platen gates, and stays open once it has opened, because these integrators
test it at every step and `reallocate_jet` can lower `npjet`. Below the gate
`rk4sys_ev` now runs the CPU build's code: until then its host path also
recomputed the Maxwell evaporative stress and the stage updates on the
device, copying the arrays at every call (600 us per step). A single-bead
evaporative RK4 run is byte-identical to the CPU build up to the
engagement. Tests 16-19 start with at least 100 beads and engage at the
first step, as before.

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
CPU; after the fix the GPU result matched the CPU one to about eight
significant digits (every row is identical since the collector rule of
2026-10-07). `examples/input-13` and `input-14` hide it at their shipped 1,000-step
length, which performs no reallocation at all, and reproduce it when extended
to 25,000 steps.

Any future device-resident array must therefore pair its `enter data` with an
`exit data delete` on every path that can reallocate the host storage.
Bookkeeping flags alone do not unmap anything.

## Explicit data region

Coordinates, bead properties, cross sections, and the Coulomb force array are
named explicitly in each OpenACC data region. No managed/unified-memory build
mode is used.

Outside the device step no jet array stays on the device: a Coulomb sum
offloaded there (at least 128 active beads) maps the jet and the force
array for that call only. Until milestones M1-M3 (2026-10-06) the runs
outside the per-test persistent gates of Tests 9--13, 16, 17, 20 and 21 used
such call-scoped Coulomb and EOM data regions.

Test 13 recorded the first persistent dynamic-topology path, behind a gate
shaped on its geometry (at least 1,000 beads and 1,280 slots) until M1. Its
mechanisms are those of the device step with insertion today: the RK and
force data stay resident as `inpjet` and `npjet` change, and collector
detection and clamping run on the device.
Nozzle insertion, including threshold checks, blocked-bead release, record
initialization, and `npjet` update, also runs on the device, and so do the
placement and the charge smoothing of the blocked bead before every force
evaluation. One serial kernel decides release, insertion, capacity
exhaustion, and removal; the outcome returns to the host as one 20-byte
record per step. On an actual insertion the host downloads only the injected
mass and charge required by the host statistics, and on a removal only the
removed record. Test 15 additionally exercises small-capacity teardown, host
reallocation, and device remapping. Compaction and capacity growth remain
host operations, after which the device step rebuilds its workspace
(`device_step_request_reset`).

The A30 run of Test 13 reproduces all 26 CPU topology events step for step.
Until 2026-10-01 its insertions fell progressively earlier (step 593 instead
of 600 for the seventh): the device placed the blocked bead only when it was
created and never smoothed its charge, so the device Coulomb sum saw a full
charge at the nozzle. That drift had been attributed to GPU RK4 rounding.
`JETSPIN_OPENACC_DISABLE_PERSISTENT=1` keeps the device step closed, so the
run executes the CPU build's code with Coulomb sums of at least 128 active
beads offloaded. `JETSPIN_TOPOLOGY_SNAPSHOT=1` writes the full active state
at every insertion and removal to `topology-state.dat`, to verify bead
metadata and properties. See the [Test 13 record](../examples/test-13.md).

Tests 16 and 17 complete the Maxwell and Kelvin--Voigt evaporation paths for
Euler, RK2, and RK4. Their one, two, or four force/stress evaluations,
intermediate states, final state update, topology operations, and capacity
rebinding use the same persistent-data strategy. CPU and A30 executions of all
three integrators and both rheologies report 111 additions, 122 removals, two
reallocations, and 89 active beads. Three-step pre-event CPU/GPU comparisons
are identical at `rtol=1e-12` and `atol=1e-13`; Maxwell XYZ geometry is also
byte-identical at written precision. Until 2026-10-07 the first insertion
came at step 4 on the CPU and step 5 on the GPU: the device topology path
inserted a bead one step late whenever the arrays had to grow for it. Since
the fix Tests 15, 16, 17 and 19 give the CPU's events step for step and
Test 18 its totals; pointwise trajectories still separate at roundoff
level, so events, pre-event agreement and transfer behaviour form the
dynamic acceptance criteria.

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
never duplicate a transfer. Since M1-M3 the runs with insertion take the
same device step and add the 20-byte topology record described above.

## Numerical validation

Compare the OpenACC path with a CPU reference built by the same `nvfortran`
version. Different Fortran compilers may use different pseudo-random-number
sequences, so a GFortran trajectory is not a suitable direct reference for a
stochastic NVIDIA build.

NVHPC 25.5, built with `CUDA_VERSION=12.9`, passes the same checks as 24.3
(2026-10-07: regression, restart, evaporation and refinement suites, every
oracle, Examples 1-8 and Tests 9-25), and its OpenACC runs give the events of
its own CPU build. Its results differ from 24.3's in the last digits, since
24.3 at `-O3` divides by multiplying with the reciprocal; trajectories that
amplify roundoff separate (Example 8, Tests 18, 19). Tests 24 and 25 take
about 6 % longer on the A30 with 25.5 (host phase 9-16 %, device phase 3-5 %
per step), so 24.3 remains the reference compiler.

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
CPU results at their output precision. These cases never exceed 25 beads, so
since 2026-10-06 the OpenACC build runs the CPU build's code on them (below
the device-step gate and the 128-bead Coulomb threshold): the comparison
checks the OpenACC build's host path. The GPU kernels are validated by
Tests 9--25, `tests/performance/dynamic/validate_evaporation.sh` (Tests 16
and 17, Euler, RK2 and RK4), the refinement runners in `tests/refinement/`
and the restart check (`tests/restart/run.sh openacc`). The separate 1,000-bead RK4, Euler,
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
force evaluation of the device step (Euler, RK2, RK4 and Platen; Maxwell and
Kelvin--Voigt, with or without evaporation; fixed bead sets and jets with
insertion, with or without refinement) the current stage is downloaded, the
complete trusted CPU force equations (including direct Coulomb) are
evaluated, and the derivatives are uploaded. Integration updates, dynamic
topology and the two-kernel tail of the Platen step
(`accelerator_platen_update`, `accelerator_platen_end_step`, final stress
derivative included) remain on the GPU: since M3 the tail is not oracled
(the oracle builds evaluated the final stress derivative on the host
before). The target is a
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
Tests 12 and 20 matched their NVFORTRAN CPU statistical outputs exactly in all
six rows and fourteen columns with both the complete-force oracle and the
Coulomb-only oracle. Before M3 these checks also covered parts of the
Platen tail, the final stress evaluation among them, which the oracle builds
then evaluated on the host; since M3 they cover the three force evaluations
of the step (at the initial state and at the two predicted states), and the
Platen force oracles of Tests 12 and 20--22 give the device step's events.

The device step covers Maxwell Euler, RK2, and RK4
evaporation (Test 16), including insertion, removal, and capacity rebinds.
All three A30 runs reproduce the CPU event totals (111 insertions, 122 removals,
two reallocations, 89 active beads), and since 2026-10-07 every event at the
CPU's step. Runtime transfer audits of the standard
build (no diagnostic macros) confirm that no state, force, stress, or
derivative array crosses the PCIe boundary between stages. Each timestep
returns only one 20-byte topology record. Insertion, removal, statistical
output, capacity growth, and the final checkpoint transfer only the records
required by those events. Full active-state transfers occur at the two
reallocations and at the final checkpoint. Until 2026-10-07 the device
inserted one step late when the arrays grew (first insertion at step 5
against 4), which was ascribed to the Coulomb summation order. CPU/GPU
trajectories are still not pointwise identical, because target-centric
Coulomb accumulation and the subsequent bending instability amplify
floating-point ordering differences in the transverse components; event
steps, pre-event agreement, and transfer behaviour are the acceptance
criteria.

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

### Development note: diagnostic strategies for new device kernels

Three strategies have been used to check a device kernel against the trusted
host code; a new kernel should get at least one of them.

- **Substitution oracle** (`nvfortran-openacc-force-oracle`,
  `nvfortran-openacc-coulomb-oracle`): the host result replaces the device
  result at every call and the run goes on with it. A trajectory or event
  difference between the oracle and the native run then isolates the
  kernel. It changes the run and transfers data at every call.
- **A/B comparison build** (`nvfortran-openacc-compare-akima`,
  `nvfortran-openacc-compare-refinement`): both results are computed, the
  device one stays authoritative, and the build prints the largest
  difference at each call. The run is unchanged, so the difference measures
  the kernel itself.
- **In-line comparison at run time**: until 2026-10-07 the variable
  `JETSPIN_COULOMB_DIAGNOSTIC=1` did an A/B comparison of the evaporative
  three-dimensional Coulomb sum inside the standard build. After each device
  sum it downloaded the state and the forces, recomputed the sum with a
  plain host loop, and printed `COULOMB_DIAGNOSTIC maxdiff=` with the bead
  and component of the largest difference (and both values above 1e-10). It
  served while the evaporative kernel was written (2026-08) and was removed
  because it covered one path only (no non-evaporative or one-dimensional
  sum), its download assumed arrays mapped by the device step and stopped
  the run outside it, and the environment lookup sat in the production path.
  The removed in-line Maxwell stage comparisons of `rk4sys_ev`
  (`JETSPIN_COMPARE_MAXWELL_*`, until M1) were of the same kind. The macro
  `JETSPIN_DISABLE_COULOMB_EVAP` (until 2026-10-07) was a partial
  substitution oracle: it kept the evaporative 3-D Coulomb sum on the host,
  but without refreshing the host state from the device, so it gave wrong
  results inside the device step; the Coulomb-oracle build replaces it.

If such a check is needed again, write it as an A/B comparison build: a
macro (as `JETSPIN_COMPARE_AKIMA`) and a Make target, so that the standard
build does not carry it; cover every path of the quantity compared (for
Coulomb: evaporative and non-evaporative, 1-D and 3-D, inside and outside the
device step); download with `update self(...) if_present`, since outside
the device step the data are already on the host; and print one summary per
call or per output interval, not per bead.

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
Test 25 reference the time-integration loop takes about 690 s on the A30
and 5863 s with the NVFORTRAN CPU build.

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
  Since 2026-10-05 the dynamic evaporative step runs on one asynchronous
  queue (below): at about 219 beads it made 31 launches, one stream
  synchronization, and one 20-byte download per step, took about 172 us,
  and kept the device busy for about 84 % of it. With its small kernels
  fused (below) it makes 15 to 18 launches and takes about 150 us.

The queue (`accelerator_queue`) carries the charge smoothing and restoring,
the placement of the inserting bead, the Coulomb, EOM, and stress kernels,
the Platen updates, the statistics, and the topology check of the Platen
device step of a run with insertion, evaporative or not (since 2026-10-05
for the evaporative step, then for the non-evaporative one). The topology
record is read back on the
queue and is the only wait of an ordinary step; every routine that moves
data between host and device waits on the queue first, so output, removals,
refinement events, and capacity resets see the completed device state. The
RK device steps and the fixed bead sets (Tests 12 and 20) are still
synchronous (milestone M4), and so are the oracle builds and runs with
`JETSPIN_OPENACC_SYNC=1`. The
5-million-step Test 25 run and Tests 21--23 are byte-identical with and
without the queue.

The Platen device step fuses its small kernels. Before each force evaluation one
serial kernel (`accelerator_platen_stage_prep`) restores the charge smoothed
for the previous evaluation, smooths it again, and places the inserting
bead; `accelerator_eom3_stage` with `fev_evap` computes the evaporation rate
and the Maxwell stress in its own loop instead of a second kernel; and one
single-gang kernel (`accelerator_platen_end_step`, after the per-bead
`accelerator_platen_update`) does the stress
update with its statistics, the placement, the topology decisions, the
freezing, and the step's statistics, after which `accelerator_topology_check`
only reads the record back. The serial pieces (`device_smooth_charge`,
`device_place_inserting_bead`, `device_topology_decide`) are
`!$acc routine seq` routines shared with the separate smoothing, placement
and topology kernels of the RK device step, which the Platen step of the
oracle builds also uses for smoothing and placement; the unfused Platen
kernels were removed in M3. The results do not change.

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
at the end of the step. Placement, smoothing, and restoring now bracket
every force evaluation on the device (placement by `place_inserting_bead`
until M3; now by `device_place_inserting_bead`, called by
`accelerator_compute_posnoinserted_3d` in the RK step and by
`accelerator_platen_stage_prep` and `accelerator_platen_end_step` in the
Platen step), the
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

1. For the dynamic Platen paths the tail of the step is two kernels since
   2026-10-06: the per-bead updates and the single-gang end of the step,
   which also evaluates the final stress. Merging them would need the
   device routines of the end of the step with explicit-shape arrays (NVHPC
   24.3 and 25.5 both fail otherwise, `docs/STATE.md`), would save one
   launch, and would make the updates grow with the bead count on one
   multiprocessor. Beyond that the step time is in the Coulomb and force
   kernels themselves.
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
