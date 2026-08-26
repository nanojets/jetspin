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

Test 23, the rheology parent, uses Yarin evaporation; Test 24 does not. That
is now a definition rather than a temporary state: the evaporative
configuration is [Test Case 25](test-25.md), whose `input.dat` differs from
this one by exactly one line.

The split exists because evaporation does not survive this configuration, and
keeping the two as separate inputs makes the comparison reproducible with a
one-keyword difference:

```text
Test 24 (evaporation no)   -> 6,000,000 steps, Program closed correctly
Test 25 (evaporation yes)  -> ERROR - numerical instability at nstep 1078685
```

Test 24 is therefore the case to use for the stationary-bead-count study it
was designed for; Test 25 is the reproducer for the open evaporation defect,
and carries its analysis.

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

## Related open work

Long single-bead-start runs with repeated remeshing exposed three
refinement-robustness defects, two fixed and one contained through the cadence
choice above; they are recorded in
[the refinement robustness note](../refinement-robustness-investigation.md).
The remaining open item bearing on this case is the evaporation instability
that keeps the evaporative configuration in [Test Case 25](test-25.md).

## Reference results

The longest clean run is 6,000,000 steps on an NVIDIA A30,
`Program closed correctly`, zero NaN and zero device errors, ending with
x = 16.00 cm, collector velocity 1986.5 cm/s, `yz` = 4.09 cm, 28.7 degrees,
path length 111.4 cm, 571 active beads, 291 topology additions and 1,616
removals. The bead count oscillates around 570 over the last two million steps
rather than drifting, so this run does demonstrate the balanced
insertion/removal regime the case was designed to reach — at 6% of the
target duration, not over the full run.

A full-length reference does not yet exist. The developed bending instability
carries several hundred active beads, so the direct Coulomb cost per step is
substantial and a 100-million-step run is correspondingly expensive.

## Use

Track the printed `n` column over time to see the stationary oscillation. The
case is excluded from the smoke and regression matrices because of its length.

- [Input file](../../examples/input-24/input.dat)
- [Refinement robustness investigation](../refinement-robustness-investigation.md)
- [Test 25, the evaporative counterpart](test-25.md)
- [Test 23, its rheology parent](test-23.md)
- [Dynamic refinement](../introduction/dynamic-refinement.md)
