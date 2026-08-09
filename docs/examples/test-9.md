# Test Case 9: 1,000-bead CPU/GPU benchmark

Test Case 9 is a performance-oriented three-dimensional Maxwell simulation.
It initializes 1,000 beads uniformly over the complete 16 cm distance from
the nozzle to the collector and executes 1,000 integration steps.

The material, aerodynamic, gravity, and perturbation parameters are derived
from Test Case 7. Two coupled features are deliberately removed:

- the orthogonal rotating field is replaced by a constant axial electric
  potential;
- evaporation is disabled because its direct Coulomb routine has not yet
  been ported to OpenACC.

Insertion, removal, multiple-step Coulomb summation, and dynamic refinement
are disabled, keeping the workload fixed at 1,000 beads. This makes the case
suitable for comparing the current direct O(N²) Coulomb implementations
without topology changes or CPU fallback from the evaporation model.

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

## Subroutine profiling

Enable the optional internal profiler without rebuilding:

```sh
JETSPIN_PROFILE=1 ./main.x
```

The terminal report separates the integrator, Coulomb calculation, bead
topology operations, statistics, scheduled output, and restart output. The
Coulomb timer is nested inside `Integrator total`; it must not be added to the
integrator time. Profiling is disabled by default to avoid adding clock calls
to production runs.

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
retain the CPU EOM path. Absolute timings vary with system load; preserve the
compiler, profiler setting, and execution order when repeating the comparison.

## Persistent-data milestone

The next milestone retains the complete Test 9 state, Coulomb force, EOM
derivatives, and RK4 scratch arrays on the A30 across timesteps. Cross-section
calculation and all four RK updates also execute on the device. Coulomb output
is consumed directly by EOM and never crosses back to the host.

The existing CPU statistics require the primary state after every timestep,
so the final RK stage still updates seven arrays on the host once per step.
This is about 56 kB per step for 1,001 beads; no stage intermediate is
transferred. A representative profiled run produced:

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

## Versioned numerical records

The CPU and A30 OpenACC `statout.dat` files from this initial measurement are
stored under [`tests/performance/test9/`](../../tests/performance/test9/).
That directory records their provenance and provides a comparison command for
future porting milestones. It remains separate from the automatic regression
matrix so that the 1,000-bead benchmark does not slow routine validation.
