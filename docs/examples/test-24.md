# Test Case 24: long run to a stationary bead count

Test Case 24 is a long production run whose purpose is to reach a
statistically stationary active-bead count: once the jet has bridged the
16 cm nozzle-to-collector gap with a developed bending instability, insertion
at the nozzle and removal at the collector should balance and the bead count
should oscillate around a mean instead of drifting. Unlike Tests 21--23 it is
not meant to finish in a smoke-test timeframe: it integrates
`final time 0.5d0` s at `timestep 5.d-9` s, about 100 million steps.

It is assembled from three existing sources rather than invented from scratch:

- **Canonical JETSPIN electrostatics** — `density charge 44000`,
  `collector distance 16`, `external potential 30.02076857`, the set shared
  with Examples 3 and 6 and with Tests 10, 15 and 16--20.
- **Rheology, air drag and topology from Test Case 23** — Maxwell rheology,
  stochastic Platen air drag, nozzle insertion, collector removal, and the
  same `dynamic refinement anchor 0.10` cm.
- **Initialization and print directives from Example 3** — the jet starts
  from a single nozzle bead (`points 1`) instead of the pre-extended
  400-element mesh that Tests 21--23 use to validate the refinement
  algorithm.

The single-bead start is what makes this case distinct: every other validated
refinement test begins from an already-extended jet, so the young-jet growth
regime had never been exercised.

## This case is non-evaporative by definition

Test 23, the rheology parent, uses Yarin evaporation; Test 24 does not, and
its evaporation directives are read but unused. The evaporative configuration
is [Test Case 25](test-25.md), whose `input.dat` differs from this one by two
lines: `evaporation yes` and `evaporation polymer frac 0.50d0`.

With Yarin's own 6 % polymer fraction the evaporative configuration stops on a
stress NaN near step 1.07 million, because the dried jet keeps its charge
while losing fifteen times its mass and cross-section. With 50 % it runs
stably; [Test Case 25](test-25.md) carries that analysis.

## The refinement cadence is load-bearing

The cadence is Example 5's — `dynamic refinement threshold 0.4` cm and
`dynamic refinement every 1.d-3` s — not Test 23's tighter `0.10` cm /
`1.d-5` s, and this is a deliberate choice. With the tighter values this
input accumulates dozens of accepted refinement events in under 200,000
steps, each re-fitting the cross-section radius from the previous event's own
output. That compounds a nonphysical thinning far beyond the
order-of-magnitude reduction electrospinning genuinely produces, collapsing
bead mass toward zero until the resulting acceleration overflows the stress
equation.

`warning(108)` fires whenever a refinement threshold is set below twenty
times the base discretization resolution. It is informational, not enforced:
Example 5 and Tests 21--23 are themselves below that ratio and remain valid
short-window references.

## Related work

Long single-bead-start runs with repeated remeshing exposed three
refinement-robustness defects, two fixed and one contained through the cadence
choice above; they are recorded in
[the refinement robustness note](../refinement-robustness-investigation.md).

## Reference results

Since 2026-10-05 Test 24 follows Test 25: CPU and GPU builds read the same
sequential Gaussian pool, and the GPU runs the persistent device path once the
jet holds 100 beads. The values below come from these builds (NVFORTRAN 24.3,
5,000,000 steps, `final time 2.5d-2` and `print time 1.d-4` in a copy of the
input, 5 % of the target duration).

NVFORTRAN CPU build, seed 317, `Program closed correctly`:

- the first bead reaches the collector at step 2,541,113;
- between steps 4 and 5 million: collector velocity 1963 cm/s, off-axis
  distance of the leading bead 3.79 cm, cone angle `angl` 26.6°, 471 – 553
  active beads (mean 511.4), path length 211.4 cm, radius at the collector
  (`rc`) 3.1 µm;
- 219 topology additions, 1,041 removals, one array reallocation.

Four more CPU seeds (318 – 321) give mean active counts between 511.3 and
512.0, a path length of 211.3 – 211.4 cm, and an off-axis distance of 3.79 cm
over the same window: the stationary regime does not depend on the noise
realization. Each CPU run takes about 14,000 s.

One NVIDIA A30 runs the same 5 million steps in about 620 s. It reproduces the
CPU run byte for byte up to step 1,379,105, where it engages the persistent
path, and then follows the CPU trajectory: with seed 317 the active bead count
first differs at step 2.66 million (2.18 million with seed 318). Between
steps 4 and 5 million the five A30 seeds (317 – 321) give 511.1 – 512.0
active beads (471 – 553), a path length of 211.3 – 211.4 cm, an off-axis
distance of 3.78 – 3.79 cm, a cone angle of 26.6 – 26.7°, and a collector
velocity of 1963 – 1964 cm/s, with 219 additions and 1,042 – 1,045 removals.

Before 2026-10-05 the CPU and the standard GPU builds drew this case's noise
step by step, and an earlier CPU run of the same seed gave, over the same
window, a path length of 212 cm, collector velocity 1962 cm/s, `yz` = 3.8 cm,
26.8 degrees, and about 512 active beads: the noise layout does not change
the stationary regime. The persistent path then existed only in a
development build; its longest run (6,000,000 steps, about 570 active beads at
the end) predates the device-side charge smoothing and placement of the
inserting bead (2026-10-01) and halved `lp`, and is not a reference.

A full-length reference does not yet exist; at the current A30 speed it takes
about 4 hours (see below).

## Where the A30 run spends its time

The A30 run has the two phases of Test 25, with the same split of work between
CPU and GPU (see [Test 25](test-25.md#where-the-a30-run-spends-its-time)):
until step 1,379,104 (up to 95 active beads) the CPU integrates the whole
step, Coulomb sums included; from step 1,379,105, where an accepted refinement
event takes the jet to 121 beads, the step runs on the
device-resident state, on one asynchronous queue with the fused small
kernels, and an ordinary step returns one 20-byte topology record. The
non-evaporative step computes its stress derivative at the new state with a
fourth EOM stage instead of the evaporative stress kernel.

Seed 317, process bound to the GPU's NUMA node: the loop takes 621 s, 114 s
in phase 1 (83 µs/step) and 507 s in phase 2 (140 µs/step); seed 318 takes
616 s. The standard GPU build of the same morning, which kept this case on the
host with only the Coulomb sums on the GPU, took 4546 s. Median time per step
over the 20,000-step print intervals:

| Active beads | CPU | A30, current | A30, before 2026-10-05 |
| --- | ---: | ---: | ---: |
| 1 – 30 | 19 µs | 21 µs | 27 µs |
| 30 – 60 | 59 µs | 65 µs | 74 µs |
| 60 – 100 | 148 µs | 159 µs | 245 µs |
| 100 – 150 | 405 µs | 113 µs | 444 µs |
| 150 – 200 | 689 µs | 116 µs | 602 µs |
| 200 – 300 | 1156 µs | 123 µs | 737 µs |
| 300 – 400 | 2515 µs | 132 µs | 1049 µs |
| 400 – 500 | 4086 µs | 145 µs | 1373 µs |
| 500 – 600 | 5037 µs | 152 µs | 1467 µs |

Four changes of 2026-10-05 produce the difference: the persistent path
itself (until then a development build), the asynchronous queue with the
fused kernels shared with Test 25, the non-evaporative Coulomb kernel, which
now gives each target bead one gang and reduces over the sources across its
vector lanes, and the removal of the nine array-descriptor uploads that kernel
made, with a wait, at every call. Over the first 2.1 million steps phase 2
averaged 575, 347, 175, and 119 µs per step after each change in turn. With
about 500 beads the A30 step costs about 150 µs against 4 – 5 ms on the
CPU. Nsight Systems at about 540 beads counts per step 15 kernel launches, one
stream synchronization, the 20-byte download, and no upload; the kernels take
122 µs, half of it in the three Coulomb sums.

## Use

Track the printed `n` column over time to see the stationary oscillation. The
case is excluded from the smoke and regression matrices because of its length.

- [Input file](../../examples/input-24/input.dat)
- [Refinement robustness investigation](../refinement-robustness-investigation.md)
- [Test 25, the evaporative counterpart](test-25.md)
- [Test 23, its rheology parent](test-23.md)
- [Dynamic refinement](../introduction/dynamic-refinement.md)
