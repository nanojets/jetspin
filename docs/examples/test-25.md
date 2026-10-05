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

One NVIDIA A30 runs the same 5 million steps in about 1940 s, against
5893 s for the CPU build. It engages the persistent path at
step 1,446,263, reads the same Gaussian pool, and stays on the CPU trajectory:
with seed 317 the first removal falls at the same step, and the active bead
count first differs at step 3.26 million.
Between steps 4 and 5 million it gives the same values as the CPU: collector
velocity 2535 cm/s, off-axis distance 2.8 cm, cone angle 19.75° (CPU 19.74°), 242 – 295
active beads, path length 112 cm, fibre radius 2.8 µm, and a 2.9 cm envelope
radius at the collector, with 221 additions and 925 removals. Seeds 318 and
319 likewise follow their CPU runs to 3.2 and 2.7 million steps and give
267.9 and 268.2 active beads and a 111.8 cm path length. Before the
persistent-path defects found on 2026-09-30 were fixed (see
[OpenACC](../introduction/openacc.md)), the A30 gave 297 active beads and a
123 cm path length on every seed.

An NVFORTRAN CPU run of Test 24 over the same interval has a 3.8 cm off-axis
distance and a 27° cone:
evaporation and the stiffening it causes make the bending loops smaller, as
reported by Yarin et al.

A full-length reference does not yet exist.

## Use

Track `n`, `yz`, and `angl` to follow the stationary regime; adding `visc`,
`gc`, or `evrc` to `printstat list` reports the viscosity, modulus, and
evaporated volume fraction of the collected jet. The case is excluded from the
smoke and regression matrices because of its length.

- [Input file](../../examples/input-25/input.dat)
- [Test 24, the non-evaporative twin](test-24.md)
- [Refinement robustness investigation](../refinement-robustness-investigation.md)
- [Dynamic refinement](../introduction/dynamic-refinement.md)
