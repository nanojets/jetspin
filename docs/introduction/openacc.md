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
fixed Tests 9--12 configuration and for the bounded dynamic evaporation paths
in Tests 16 and 17. One kernel is launched per integrator stage. Its local
three-point curvature calculation is device-side: an iteration reads the
current bead and its two neighbours. No host curvature array is built or
transferred. Disabling fused multiply-add preserves the accepted trajectory
for the initially straight geometry. Configurations outside the explicitly
validated gates use the original CPU EOM path.

## Explicit data region

Coordinates, bead properties, cross sections, and the Coulomb force array are
named explicitly in each OpenACC data region. No managed/unified-memory build
mode is used.

Configurations outside Tests 9--13 continue to use separate call-scoped
Coulomb and EOM data regions.

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

For stochastic Platen integration, the complete Gaussian history is generated
on the CPU before loop timing in a fixed step/bead/component/draw order. Test
12 requires 6,006,000 doubles (48,048,000 bytes). The OpenACC build transfers
this history once during initialization and indexes it on the device; the CPU
path indexes the same layout. No random generation or noise transfer occurs
inside the measured loop.
The history is capped at 100,000,000 doubles (about 763 MiB). Longer fixed
runs wrap to its beginning, preserving CPU/GPU reproducibility while making
the noise periodic after the stored interval.
The host no longer receives stage intermediates or Coulomb forces.

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

For Test 16 Maxwell force isolation, use the development-only target:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-host-forces GPUCC=80
```

This keeps the OpenACC executable but enables both
`JETSPIN_DEV_HOST_COULOMB_EVAP` and `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY`.
The evaporative Coulomb and Maxwell geometric forces are evaluated by the
trusted host routines and copied back to the device at each call. It is a
diagnostic build only; it is intentionally unsuitable for performance
measurements or production runs. Additional preprocessor switches can be
passed through `FPPFLAGS_EXTRA`.

For the equivalent Test 17 Kelvin--Voigt isolation, use:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-kv-host-forces GPUCC=80
```

This target enables `JETSPIN_DEV_HOST_COULOMB_EVAP` and
`JETSPIN_DEV_HOST_KV_FORCES`. Each complete Kelvin--Voigt force stage is
evaluated with the trusted CPU routines, while RK state updates remain on the
GPU. It exists only to distinguish force-kernel differences from integrator
and topology effects.

For development-only diagnosis of CPU/GPU Coulomb differences, compile the
evaporative Coulomb module with `-DJETSPIN_DEV_HOST_COULOMB_EVAP`. This macro
disables the device Coulomb kernel, evaluates the full evaporative Coulomb sum
on the host at every call, and updates the device copy of `ycf` before the
next device force stage. It is intentionally not enabled by the standard
Makefile targets because it adds a host/device transfer at every RK stage and
is unsuitable for performance measurements.

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

## Next porting stages

1. Investigate packing the maximum stress and bead index into one deterministic
   reduction so that its two follow-up kernels can also be removed.
2. Move capacity growth itself to a device-aware allocator if reallocation
   frequency becomes significant; current transfers occur only on resize
   events, not on ordinary timesteps.
3. Port the local Akima coefficient loops, replace the interpolation interval
   scan with a GPU-suitable search, and then address dynamic refinement.
4. Evaluate one-GPU-per-rank MPI execution only after the single-GPU numerical
   path is stable.

The intended steady state is a persistent device-resident simulation with
host transfers for output, checkpoints, and topology changes—not a transfer
of the complete jet at every timestep.
