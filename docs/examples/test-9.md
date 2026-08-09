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

## Versioned numerical records

The CPU and A30 OpenACC `statout.dat` files from this initial measurement are
stored under [`tests/performance/test9/`](../../tests/performance/test9/).
That directory records their provenance and provides a comparison command for
future porting milestones. It remains separate from the automatic regression
matrix so that the 1,000-bead benchmark does not slow routine validation.
