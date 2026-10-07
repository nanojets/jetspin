# Test Case 9: 1,000-bead CPU/GPU benchmark

Test Case 9 is a performance-oriented three-dimensional Maxwell simulation.
It initializes 1,000 beads uniformly over the complete 16 cm distance from
the nozzle to the collector and executes 1,000 integration steps.

The material, aerodynamic, gravity, and perturbation parameters are derived
from Test Case 7. Two coupled features are deliberately removed:

- the orthogonal rotating field is replaced by a constant axial electric
  potential;
- evaporation is disabled so this case isolates the non-evaporative force and
  integration paths.

Insertion, removal, multiple-step Coulomb summation, and dynamic refinement
are disabled, keeping the workload fixed at 1,000 beads. The leading bead
starts on the collector; since 2026-10-07 it is frozen there, like any bead
that reaches the collector, with or without removal (before, it crossed the
collector plane by about 6e-4 cm over the 1,000 steps). The leading bead of
Tests 10-12, 15 and 20 also starts on the collector and is frozen in the same
way. The fixed workload makes the case
suitable for comparing the direct O(N²) Coulomb implementations without
topology or evaporation coupling.

Build and run CPU and GPU executables with:

```sh
module purge
module use /opt/nvidia/hpc_sdk/modulefiles
module load nvhpc/24.3

make -C source -f ../build/Makefile clean
make -C source -f ../build/Makefile nvfortran EX=jetspin-cpu.x

make -C source -f ../build/Makefile clean
make -C source -f ../build/Makefile nvfortran-openacc EX=jetspin-gpu.x
```

Run each executable from a separate directory containing a copy of the Test 9
input. One reproducible layout is:

```sh
mkdir -p /tmp/jetspin-test9/cpu /tmp/jetspin-test9/gpu
cp execute/jetspin-cpu.x examples/input-9/input.dat /tmp/jetspin-test9/cpu/
cp execute/jetspin-gpu.x examples/input-9/input.dat /tmp/jetspin-test9/gpu/

cd /tmp/jetspin-test9/cpu
./jetspin-cpu.x | tee run.log

cd /tmp/jetspin-test9/gpu
./jetspin-gpu.x | tee run.log
```

Compare the reported `Time-integration loop wall time`, which excludes
simulation preparation and accelerator initialization. Run the two programs
sequentially on an otherwise idle system, and record the compiler, GPU, and
build options with every benchmark result.

- [Input file](../../examples/input-9/input.dat)
- [LaTeX manual section](../../manual/test9.tex)

## Initial reference measurement

An initial run with NVFORTRAN 24.3 and one NVIDIA A30 produced:

| Backend | Loop wall time | Throughput |
| --- | ---: | ---: |
| NVFORTRAN CPU | `40.659856 s` | `24.594 steps/s` |
| OpenACC A30 | `3.785055 s` | `264.197 steps/s` |

The measured speedup was `10.742x`. Both trajectories retained exactly 1,000
beads and passed the numerical comparison (`rtol=1e-6`, `atol=1e-9`), with a
maximum observed absolute difference of `1.7e-15`. These timings are an
initial reference, not a portable performance guarantee.

## Development oracles

The OpenACC diagnostic targets use the same interface as the evaporation
tests:

```sh
make -C source -f ../build/Makefile nvfortran-openacc-force-oracle GPUCC=80
make -C source -f ../build/Makefile nvfortran-openacc-coulomb-oracle GPUCC=80
```

The complete-force oracle evaluates all four RK4 force stages through the
trusted CPU equations and uploads their derivatives. The Coulomb-only oracle
moves only the four direct sums to the host. In an A30 check made before
2026-10-07 against the versioned NVFORTRAN CPU output, the standard GPU path, complete-force
oracle, and Coulomb-only oracle had worst normalized differences of
`7.97e-5`, `3.00e-7`, and `8.52e-5`, respectively. This shows that the Test 9
RK4 discrepancy is not dominated by Coulomb accumulation order. Both oracle
builds deliberately transfer data at every stage and must not be benchmarked.

## Subroutine profiling

Enable the optional internal profiler without rebuilding:

```sh
JETSPIN_PROFILE=1 ./main.x
```

The terminal report lists `Integrator total`, `Coulomb (nested)`, the bead
topology operations (`Add bead`, `Remove bead`, `Erase bead`, `Breakup
check`), `Statistics`, `Scheduled output`, `Restart output`, `EOM evaluation
(nested)`, and `RK update (nested)` (`source/profiling_mod.f90`). The three
nested timers are inside `Integrator total`; they must not be added to the
integrator time. The EOM and RK update regions (`prof_eom`,
`prof_rk_update`) are instrumented only in the CPU code of `eulsys`, `rk2sys`
and `rk4sys` (`source/integrator_mod.f90`); the OpenACC device step
(`source/device_step_mod.f90`) has no such regions, so a run that takes it,
as Test 9 does from its first step, reports zero for them and, within the
integrator, times only `Coulomb (nested)`. The EOM and RK figures below were
measured with the accelerator paths of the porting milestones, which ran
inside those integrators. Profiling is disabled by default to avoid adding
clock calls to production runs.

The report is printed on the terminal only and contains cumulative time,
percentage of the complete temporal loop, and call count for each selected
region. It instruments these regions explicitly; it is not an automatic
profile of every Fortran subroutine. For example, the fourth-order Runge--Kutta
integrator performs four Coulomb evaluations per step, so this case reports
4,000 Coulomb calls.

Use the same profiler setting on both sides of a timing comparison because
clock reads add a small overhead. In MPI execution the detailed region table
is local to rank 0 and is not a maximum-over-ranks profile. The synchronized
complete-loop timer remains the appropriate value for MPI scaling.

An initial profiled Test 9 run produced:

| Region | NVFORTRAN CPU | OpenACC A30 | A30 share of loop |
| --- | ---: | ---: | ---: |
| Complete temporal loop | `39.813621 s` | `3.320607 s` | `100%` |
| Integrator total | `39.804322 s` | `3.314476 s` | `99.82%` |
| Coulomb, nested | `38.232985 s` | `2.164045 s` | `65.17%` |
| Statistics | `0.006663 s` | `0.004526 s` | `0.14%` |

The measured Coulomb-region speedup was approximately `17.67x`, while the
complete-loop speedup was approximately `11.99x`. Later porting milestones
should record the same table to show whether time moves from Coulomb into the
remaining host-side integrator work or data transfers.

A subsequent profiler refinement separated the remaining A30 integrator
time into explicitly instrumented regions:

| Region | A30 time | Share of loop | Calls |
| --- | ---: | ---: | ---: |
| Coulomb, nested | `2.026551 s` | `56.28%` | 4,000 |
| EOM evaluation, nested | `1.420609 s` | `39.45%` | 4,000 |
| RK update, nested | `0.014712 s` | `0.41%` | 3,000 |

This measurement identifies the three-dimensional equation-of-motion chain,
not the inexpensive RK vector updates, as the next relevant accelerator
target. Absolute timings vary between runs; the numerical trajectory matched
the versioned A30 baseline exactly at the saved output precision.

## Device EOM and curvature milestone

This section and the following ones up to the selective-output milestone
record the porting milestones of August 2026, which the common device step
(below) has superseded; their normalized differences refer to the records of
that time.

The first equation-of-motion port uses one explicit OpenACC kernel per RK4
stage. The Test 9 force assembly and its local three-point curvature
construction execute entirely on the device. Each bead reads only its own
coordinates and those of its two neighbours. No curvature array is constructed
or transferred by the host.

The initially straight jet makes the circumcentre calculation sensitive to
floating-point contraction. The OpenACC target therefore uses NVFORTRAN's
`nofma` GPU option. Without it, fused multiply-add changed the trajectory
beyond the accepted tolerance; with it, the device calculation preserves the
versioned reference while remaining parallel.

An updated paired measurement compiled both executables from the same source
with NVFORTRAN 24.3 and ran them sequentially with profiling enabled:

| Region | NVFORTRAN CPU | OpenACC A30 | Speedup | A30 share | Calls |
| --- | ---: | ---: | ---: | ---: | ---: |
| Complete temporal loop | `38.616119 s` | `3.159748 s` | `12.22x` | `100%` | 1 |
| Integrator total | `38.609060 s` | `3.152703 s` | `12.25x` | `99.78%` | 1,000 |
| Coulomb, nested | `37.059132 s` | `2.250507 s` | `16.47x` | `71.22%` | 4,000 |
| EOM evaluation, nested | `1.398234 s` | `0.736914 s` | `1.90x` | `23.32%` | 4,000 |
| RK update, nested | `0.018556 s` | `0.030343 s` | `0.61x` | `0.96%` | 4,000 |
| Statistics | `0.005216 s` | `0.004980 s` | `1.05x` | `0.16%` | 1,000 |

The complete temporal loop is `12.22x` faster. Direct Coulomb evaluation
obtains the largest regional speedup and remains the dominant A30 cost. The
device EOM, including curvature, is `1.90x` faster. RK updates are slower in
this call-scoped implementation but account for less than one percent of the
GPU loop. Their call count is 4,000 because all four RK stages are
instrumented separately.

The full six-snapshot, fourteen-column trajectory passed against the A30
baseline with `rtol=1e-6` and `atol=1e-9`; its worst normalized difference was
`7.82e-5`. That initial EOM kernel was deliberately limited to the fixed Test 9
configuration and used call-scoped data transfers. Unsupported configurations
retained the CPU EOM path. Absolute timings vary with system load; preserve the
compiler, profiler setting, and execution order when repeating the comparison.

## Persistent-data milestone

The next milestone retains the complete Test 9 state, Coulomb force, EOM
derivatives, and RK4 scratch arrays on the A30 across timesteps. Cross-section
calculation and all four RK updates also execute on the device. Coulomb output
is consumed directly by EOM and never crosses back to the host.

The CPU statistics of that time required the primary state after every
timestep, so the final RK stage still updated seven arrays on the host once
per step. This is 56,056 bytes per step for the 1,001 array entries (indices
0-1,000) of the 1,000-bead jet; no stage intermediate was transferred. A representative profiled run produced:

| Region | Call-scoped A30 | Persistent A30 | Change |
| --- | ---: | ---: | ---: |
| Complete temporal loop | `3.159748 s` | `2.451978 s` | `1.29x` faster |
| Integrator total | `3.152703 s` | `2.445458 s` | `1.29x` faster |
| Coulomb, nested | `2.250507 s` | `2.006970 s` | `1.12x` faster |
| EOM evaluation, nested | `0.736914 s` | `0.112892 s` | `6.53x` faster |
| RK update, nested | `0.030343 s` | `0.191860 s` | includes host synchronization |

The complete loop is `15.75x` faster than the paired NVFORTRAN CPU run of
`38.616119 s`. The persistent EOM time demonstrates the benefit of eliminating
four sets of input/output mappings per step. RK time increases because its
timer now includes the single full-state device-to-host synchronization; it
should not be interpreted as slower RK arithmetic.

The persistent trajectory passes the unchanged A30 baseline with `rtol=1e-6`
and `atol=1e-9`. Its worst normalized difference is `7.8e-5`, and it differs
from the preceding call-scoped GPU trajectory by at most `2.0e-7` normalized.

## Device-resident statistics milestone

Path length and maximum stress are now accumulated on the device after each
step. The statistic counters remain device resident, and the complete primary
state is copied to the host only at the five scheduled 200-step samples. The
final restart detects that the step-1000 snapshot is already current and does
not download it again. The former 56,056-byte transfer on every timestep is
therefore absent. `NV_ACC_NOTIFY=2` confirms exactly 35 array downloads and 20
scalar-statistic downloads: five synchronization events in total.

A representative A30 run took `2.468580 s` for the complete loop. Statistics
took `0.056676 s`; this is mostly the launch overhead of the small reduction
kernels rather than data transfer. The trajectory passed the unchanged A30
baseline with `rtol=1e-6` and `atol=1e-9`; the worst normalized difference was
`7.97e-5`. This measurement is the reference for the following fusion step.

## Fused RK-statistics milestone

The path-length and maximum-stress reductions are subsequently fused with the
final RK4 state update. Each thread reconstructs the updated position of its
next bead from the same RK coefficients, avoiding a neighbour race without an
extra state-update kernel. Two small follow-up kernels retain the deterministic
last-index rule for equal maximum stresses.

In a representative A30 run, the statistics region decreased from
`0.056676 s` to `0.036373 s`, and the complete loop took `2.457896 s`. An
extended comparison that also printed `lp`, `mxst`, and `mxsx` matched the
preceding commit exactly for maximum stress and its position. Path length
differed by at most `2e-8 cm` because of reduction ordering. The standard Test
9 trajectory still passed the A30 baseline with a worst normalized difference
of `7.97e-5`.

## Selective-output milestone

The output path no longer downloads all seven state arrays for ordinary
statistics. At each configured `print time`, it transfers only the seven state
values of the selected bead plus four accumulated statistics. The interval is
derived from `print time / timestep`; 200 steps is specific to the standard
Test 9 input and is not hardcoded.

`NV_ACC_NOTIFY=2` reported 280 bytes for five point samples, 140 bytes for the
associated scalar statistics, and one 56,056-byte full-state download for the
final restart: 56,476 bytes in total. The preceding implementation transferred
280,420 bytes. A diagnostic input with XYZ output enabled correctly reverted
to complete snapshots at the requested XYZ cadence and produced a valid
trajectory file. The standard numerical comparison retained the worst
normalized difference of `7.97e-5`.

## Coulomb kernel with one gang per target (2026-10-05)

The non-evaporative three-dimensional Coulomb kernel now spreads the sources
of each target bead over the vector lanes of one gang and combines them with a
reduction, as the evaporative kernel does since 2026-09-30, instead of giving
each target one thread, and takes its arrays as assumed-shape dummies, which
removes nine descriptor uploads and a wait per call. On one A30, with the
process bound to the GPU's NUMA node, the Test 9 loop takes `0.508 s` against
`2.121 s` with the previous kernel on the same node (NVFORTRAN CPU
`24.94 s`). Only the summation order changes: the output still passes the CPU
and A30 records with a worst normalized difference of `7.99e-5`.

## One device step and the collector rule (2026-10-06/07)

Since 2026-10-06 Tests 9-11 run the common device step (`device_rk_step` in
`source/device_step_mod.f90`) from the first step, synchronous and not fused
yet (milestone M4 of `docs/STATE.md`): 22, 10 and 14 kernel launches per
step for RK4, Euler and RK2 since 2026-10-07, one of which,
`accelerator_freeze_at_collector`, freezes the beads that reach the
collector (21, 9 and 13 on 2026-10-06). Since 2026-10-07 a bead that reaches the
collector is frozen there: the leading bead of this jet, which starts on the
collector, stays at 16 cm (before, it crossed the plane by 6e-4 cm in 1,000
steps, `x = 16.00064` cm at the last row). With it the A30 `statout.dat` rows
equal the CPU's in every printed row; the normalized differences quoted in
the sections above predate this change.

## Versioned numerical records

The CPU and A30 OpenACC `statout.dat` records are stored under
[`tests/performance/test9/`](../../tests/performance/test9/). That directory
records their provenance and provides a comparison command for future porting
milestones. It remains separate from the automatic regression matrix so that
the 1,000-bead benchmark does not slow routine validation. The records were
regenerated on 2026-10-07 for the collector rule (NVFORTRAN 24.3, CPU and A30
rows identical; NVFORTRAN 25.5 gives the same rows).
