# Test Case 25: evaporative counterpart of Test 24

Test Case 25 is Test Case 24 with Yarin evaporation enabled. The two
`input.dat` files differ by exactly two directives:

```text
input-24:  evaporation no      evaporation polymer frac 0.06d0 (unused)
input-25:  evaporation yes     evaporation polymer frac 0.50d0
```

Everything else is shared: canonical JETSPIN electrostatics
(`density charge 44000`, `collector distance 16`,
`external potential 30.02076857`), Maxwell rheology (`viscosity 20`,
`elastic modulus 50000`), stochastic Platen air drag, nozzle insertion and
collector removal from Test 23, a single-nozzle-bead start from Example 3,
integration to `final time 0.5d0` s at `timestep 5.d-9` s, and Example 5's
refinement cadence (`dynamic refinement threshold 0.4` cm,
`dynamic refinement every 1.d-3` s).

The evaporation law is Yarin's (2001) unchanged: `bconstant 7`,
`mconstant 0.1`, `tconstant 1`, relative humidity 0.165 at 293.15 K. Only the
initial polymer mass fraction is raised from Yarin's 6 % to 50 %.

## What evaporation does in this case

With `polymer frac 0.50` the solidification cutoff `cp = 0.9` is reached at
`V/V0 = cp0/0.9 = 0.556`. Every bead collected in the reference run has
reached it, so at the collector the jet is 90 % polymer and 10 % solvent by
mass, and 89 % of the initial solvent has evaporated. At the cutoff the Yarin
law gives

| Quantity | Nozzle | Collector | Ratio |
| --- | ---: | ---: | ---: |
| Polymer mass fraction | 0.50 | 0.90 | 1.8 |
| Element volume and mass | `V0` | `0.556 V0` | 0.556 |
| Viscosity | 20 P | 49.7 P | 2.49 |
| Relaxation time | 0.4 ms | 0.72 ms | 1.80 |
| Elastic modulus | 5.0e4 | 6.9e4 | 1.38 |

## Why the polymer fraction is 0.50

With Yarin's own 6 % polymer fraction this configuration does not survive.
It stops on the code's numerical instability check at step 1,078,685 on an
NVIDIA A30 (bead 30, `x = 12.68` cm) and at step 1,067,695 with the
NVFORTRAN CPU build. The cause is the parameter regime, not the numerics:

- In the model an evaporating element keeps its charge while its mass and
  cross-section fall with its volume. With `cp0 = 0.06` the cutoff is
  `V/V0 = 0.067`, so the charge-to-mass ratio grows fifteenfold, and once the
  jet is thin it dries in well under a millisecond.
- With `tconstant 1` (relaxation time proportional to `cp`) the elastic
  modulus grows only 2.9 times while the cross-section shrinks fifteenfold.
  The elastic resistance of a dried element, `G·A`, is therefore about five
  times weaker than that of the fresh jet.
- The dried jet head is then pushed away by Coulomb repulsion: its off-axis
  distance jumps from 3 to 20 cm between steps 1.00 and 1.04 million. The next
  accepted refinement event inserts over 160 beads with segments of about
  `1e-5` cm, and a stress NaN follows within a few thousand steps.

Diagnostic variants of the same input (1.5 million steps, NVFORTRAN CPU)
confirm that the outcome follows the elastic resistance of the dried jet, not
the viscosity increase alone:

| Variant, `cp0 = 0.06` | Viscosity, modulus at cutoff | `G·A` at cutoff | Outcome |
| --- | --- | ---: | --- |
| `bconstant 0` | ×1, ÷15 | 1/225 | head ejected to 2100 cm off-axis; NaN at 1,016,381 |
| `bconstant 0`, `tconstant 0` | ×1, ×1 | 1/15 | NaN at 1,019,317 |
| Yarin (`bconstant 7`, `tconstant 1`) | ×44, ×2.9 | 1/5 | NaN at 1,067,695 |
| `tconstant 0`, modulus ×2, ×5, ×10 | same factor for both | 0.13 – 0.67 | NaN between 1.017 and 1.066 million |
| `tconstant 0`, modulus ×15 | ×15, ×15 | 1 | reaches the collector, NaN at 1,430,536 |
| `tconstant 0`, modulus ×30 or ×44 | same factor for both | 2 – 2.9 | clean to 1.5 million |
| `evaporation no` (Test 24) | — | 1 | clean |

With `evaporation umidity 1.0` (no net evaporation) the evaporative path
reproduces Test 24 bit for bit when `bconstant 0` and `tconstant 0`, and to
`2.6e-8` relatively with the Yarin law (the recomputed `cp` equals `cp0` to
roundoff), so the evaporative code path itself adds nothing spurious.

Raising the initial polymer fraction bounds the mass loss, and therefore the
growth of the charge-to-mass ratio, without touching the Yarin law. With
`tconstant 0` and `bconstant` chosen so that viscosity and modulus only double
at the cutoff, polymer fractions from 0.20 to 0.60 (0.20 at relative
humidity 0.165, the others at both 0.165 and 0.9) all ran 5 million steps
cleanly. The humidity changes only how fast the cutoff is reached: every
collected bead had reached it in all cases.

Test 25 therefore changes only the polymer fraction and keeps all three Yarin
coefficients, which carry the literature justification: `B` and `m` were
fitted by Yarin et al. to the envelope cone of a 6 % aqueous PEO jet, and the
relaxation-time law is their stated rheological assumption. With the Yarin
law, 0.50 and 0.577 both run 5 million steps cleanly; 0.577 would limit the
viscosity increase to exactly two (modulus ×1.28), and 0.50 was adopted as
the round value.

Yarin's fit was obtained for a much more viscous solution than JETSPIN's
canonical rheology (`mu0 = 1e4` P, `theta0 = 10` ms, `G0 = 1e6`, nozzle radius
150 µm against 20 P, 0.4 ms, `5e4` and 50 µm here), so the 6 % value is not
directly transferable to this configuration.

## Reference results

NVFORTRAN 24.3 CPU build, truncated at 5,000,000 steps (5 % of the target),
`Program closed correctly`, no NaN (`final time 2.5d-2` and `print time 1.d-4`
in a copy of the input):

- the jet first reaches the collector at step 2,120,839;
- between steps 4 and 5 million: collector velocity 2535 cm/s, off-axis
  distance of the leading bead 2.8 cm (2.5 – 3.2), cone angle `angl` 19.7°
  (17.4 – 22.4°), 242 – 295 active beads, path length 112 cm;
- fibre radius at the collector (`rc`, computed from the evaporated volume)
  2.8 µm (2.3 – 3.0 µm), against a nozzle radius of 50 µm;
- the envelope cone reconstructed from `traj.xyz` over the same interval
  opens from the nozzle with a local half-angle growing from about 8° to 12°
  and reaches a radius of 2.9 cm (95th percentile) at the collector;
- 20 accepted refinement events, 221 topology additions, 922 removals, no
  array reallocations.

The noise realization hardly matters: four more CPU seeds (318 – 321) give
mean active counts between 267.8 and 268.4, a path length of 111.8 cm, and an
off-axis distance of 2.78 cm over the same window.

One NVIDIA A30 runs the same 5 million steps in about 1570 s, against
5893 s for the CPU build (see [where the time goes](#where-the-a30-run-spends-its-time)).
It reads the same Gaussian pool and reproduces the CPU run byte for byte up
to step 1,446,413, where it engages the persistent path; afterwards it stays
on the CPU trajectory, and with seed 317 the active bead count first differs
at step 4.64 million.
Between steps 4 and 5 million it gives the same values as the CPU: collector
velocity 2535 cm/s, off-axis distance 2.8 cm, cone angle 19.71° (CPU 19.74°), 242 – 295
active beads, path length 112 cm, fibre radius 2.8 µm, and a 2.9 cm envelope
radius at the collector, with 221 additions and 922 removals, as on the CPU.
With the build of 2026-09-30, seeds 318 and 319 likewise followed their CPU
runs to 3.2 and 2.7 million steps and gave 267.9 and 268.2 active beads and a
111.8 cm path length. Before the
persistent-path defects found on 2026-09-30 were fixed (see
[OpenACC](../introduction/openacc.md)), the A30 gave 297 active beads and a
123 cm path length on every seed.

An NVFORTRAN CPU run of Test 24 over the same interval has a 3.8 cm off-axis
distance and a 27° cone:
evaporation and the stiffening it causes make the bending loops smaller, as
reported by Yarin et al.

A full-length reference does not yet exist.

## Where the A30 run spends its time

An A30 run of Test 25 has two phases. The persistent device path of
`platen_ev` needs at least 100 beads (`npjet>=100`), so it stays closed while
the jet grows from the single nozzle bead. With seed 317 it opens at step
1,446,413, when an accepted refinement event takes the jet from 95 to 120
active beads. The next refinement event that grows the capacity (step
2,046,413, from 222 to 375 beads) resets the device state, and the path
engages again in the same step.

| Work in each timestep | Phase 1: steps 1 – 1,446,412, 1 – 95 active beads | Phase 2: steps 1,446,413 – 5,000,000, 120 – 300 active beads |
| --- | --- | --- |
| Platen predictor, velocity, position, and stress updates | CPU | GPU |
| Forces other than Coulomb (viscoelastic, surface tension, evaporation, air drag, external field) and the noise read from the Gaussian pool | CPU | GPU |
| Coulomb sum, three per step | CPU, below 128 active beads (until 2026-10-05 the GPU: each call uploaded 13 arrays, ran a cross-section and a Coulomb kernel, and downloaded 2 arrays) | GPU, on the device-resident state |
| Nozzle charge smoothing and restoring, placement of the inserting bead | CPU | GPU |
| Insertion, removal, and collector-freezing decisions | CPU | GPU; one 20-byte record returns to the CPU |
| Removal bookkeeping (counters, collected-bead statistics) | CPU | CPU; the removed bead is downloaded and its entries uploaded again |
| Path-length and maximum-stress statistics | CPU | accumulated on the GPU, downloaded at print steps |
| Refinement threshold scan, from `1.d-3` s after the last accepted event until the next one | CPU | GPU; three values return to the CPU |
| Accepted refinement event (3 in phase 1, 17 in phase 2) | CPU | Akima interpolation and reconstruction on the GPU; mass-boundary walk, target mesh, and conservation on the CPU, then one upload |
| `statout.dat`, `traj.xyz`, and `save.dat` every 20,000 steps | CPU | CPU, after one download of the full state |

The CPU generates the Gaussian pool before the loop and copies it to the GPU
once. In this run every refinement scan of phase 2 was accepted at its first
step, so an ordinary phase-2 step moves only the 20-byte record. Phase 1 now
runs entirely on the CPU and gives the same `statout.dat`, byte for byte, as
the CPU build.

### Measured times

A30 runs with each process bound to the CPU cores and memory of its GPU's
NUMA node (`numactl --cpunodebind --membind`), seed 317, 5 million steps,
`print time 1.d-4` (NVHPC 24.3); phase times come from the wall-clock time
of the printed lines. The CPU runs are one bound run stopped at the end of
phase 1 and one unbound full run:

| Build | Phase 1 | Phase 2 | Time-integration loop |
| --- | ---: | ---: | ---: |
| A30, current | 156 s (108 µs/step) | 1418 s (399 µs/step) | 1574 s |
| A30, 2026-10-01 | 429 s, 430 s (298 µs/step) | 1410 s, 1425 s (399 µs/step) | 1840 s, 1856 s |
| A30, before the 2026-10-01 transfer reduction | 431 s, 433 s (299 µs/step) | 1585 s, 1581 s (445 µs/step) | 2016 s, 2014 s |
| CPU | 166 s bound, 156 s unbound (114 and 108 µs/step) | 5708 s (1606 µs/step) | 5863 s |

The CPU loop time comes from the timestamps: with NVFORTRAN the code's own
`Time-integration loop wall time` wrapped after 2147 s until the 64-bit fix of
2026-10-05 and reported 1568 s for that run.

The builds of 2026-09-30 and 2026-10-01 run the same code in phase 1 and give
byte-identical `traj.xyz`. On this path they differ only in the topology
check, which now returns one 20-byte record instead of two uploads and five
downloads per step, and in the refinement scan: the 2026-10-01 build saves
46 µs per phase-2 step and 8.3 % of the whole run.

Median time per step over the 20,000-step print intervals, by active-bead
count:

| Active beads | CPU | A30, current | A30, 2026-10-01 | A30, before 2026-10-01 |
| --- | ---: | ---: | ---: | ---: |
| 1 – 30 | 31 µs | 30 µs | 247 µs | 249 µs |
| 30 – 60 | 91 µs | 92 µs | 293 µs | 294 µs |
| 60 – 100 | 258 µs | 258 µs | 388 µs | 389 µs |
| 100 – 150 | 502 µs | 398 µs | 401 µs | 459 µs |
| 150 – 200 | 787 µs | 396 µs | 401 µs | 455 µs |
| 200 – 250 | 1513 µs | 400 µs | 399 µs | 448 µs |
| 250 – 300 | 1833 µs | 403 µs | 402 µs | 447 µs |

### Phase 1: the Coulomb sum on the host

Until 2026-10-05 phase 1 offloaded the three Coulomb sums of each step and
was 2.6 – 2.8 times slower than the CPU build. The `JETSPIN_PROFILE=1`
timers of a run stopped at step 1,446,411 showed why: the 4,339,233 Coulomb
calls took 340 s, 78 µs each, 79 % of the loop, while the rest of the step
cost 62 µs on the host. The CPU build spends 79 s on the same calls (18 µs
each) and 59 µs per step on the rest. Nsight Systems at about 94 beads
showed what a call cost: the two kernels ran for 3.4 and 9.2 µs and the 15
copies for 23 µs, so the device was busy for less than half of the call; the
rest was the launch, copy-enqueue, and six stream synchronizations of each
call.

A run that has not engaged a persistent path now computes a 3-D Coulomb sum
on the host while the jet has fewer than 128 active beads
(`JETSPIN_OPENACC_COULOMB_MIN_BEADS`; 0 offloads every call, as before). A
persistent run always sums on the device. The threshold comes from runs of
the first 2.1 million steps of Tests 25 and 24 with the persistent path kept
closed (`JETSPIN_OPENACC_DISABLE_PERSISTENT=1`), the Coulomb sum either always
offloaded or always on the host. One step makes three Coulomb calls, so a
third of the difference of the median step times is the extra cost of one
offloaded call:

| Active beads | Test 25, evaporative | Test 24, non-evaporative |
| --- | ---: | ---: |
| 1 – 20 | +72 µs | +66 µs |
| 40 – 60 | +64 µs | +66 µs |
| 80 – 100 | +33 µs | +45 µs |
| 120 – 140 | −20 µs | +12 µs |
| 160 – 200 | −95 µs | −41 µs |
| 250 – 300 | −381 µs | −286 µs |

Offloading pays above about 110 beads with evaporation and 150 without; 128
lies between the two crossovers. With it, phase 1 of Test 25 takes 156 s
instead of 430 s, the time of the CPU build, and the first 2.1 million steps
of Test 24, which never engages a persistent path, take 586 s instead of
854 s.

### Phase 2

In phase 2 the step time stays at about 400 µs from 120 to 300 beads, while
the CPU build grows from 0.5 to 1.8 ms. Nsight Systems at about 215 beads
counts per step 31 kernel launches, 59 stream synchronizations, one 20-byte
download, and no upload; the kernels take 220 µs, of which 66 µs are the
three Coulomb sums with their cross sections and 89 µs the three force
stages and the stress updates. At about 273 beads, with collector removal
active, the counts are 34 launches and 65 synchronizations, and the kernels
take 236 µs (73 µs Coulomb). The device is therefore busy for 55 – 60 % of
the step, and the remainder is the launch and synchronization latency of
the small kernels issued one after the other. Issuing them on one
asynchronous queue with a single synchronization per step targets that
idle time.

### Process placement

Process placement matters at this scale. Without binding, the operating
system may run the host process on the other socket: an unbound phase-1 run,
seen on a core of the other socket, took 464 s instead of 431 s, and an
unbound comparison of the 2026-09-30 and 2026-10-01 builds gave 1932 and
1935 s against 2001 and 2002 s (3.4 %) instead of 8.3 %. Timing comparisons
must therefore bind each run to its GPU's NUMA node.

## Use

Track `n`, `yz`, and `angl` to follow the stationary regime; adding `visc`,
`gc`, or `evrc` to `printstat list` reports the viscosity, modulus, and
evaporated volume fraction of the collected jet. The case is excluded from the
smoke and regression matrices because of its length.

- [Input file](../../examples/input-25/input.dat)
- [Test 24, the non-evaporative twin](test-24.md)
- [Refinement robustness investigation](../refinement-robustness-investigation.md)
- [Dynamic refinement](../introduction/dynamic-refinement.md)
