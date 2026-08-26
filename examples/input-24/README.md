# Test Case 24 input

Test Case 24 is a long production run intended to reach a statistically
stationary active-bead count: the number of beads inserted at the nozzle and
the number removed at the collector should balance out and oscillate around
a mean value once the jet has fully bridged the 16 cm nozzle-to-collector
domain with a developed bending instability.

It combines three sources:

- Canonical JETSPIN electrostatics: `density charge 44000`,
  `collector distance 16`, `external potential 30.02076857`. This is the
  parameter set used by Examples 3 and 6 and by Test Cases 10, 15 and 16-20.
- Rheology, air drag, insertion and removal from Test Case 23: Maxwell
  rheology, stochastic Platen air drag, nozzle insertion, collector removal.
- Initialization and print directives from Example 3: the jet starts from a
  single nozzle bead (`points 1`) instead of the pre-extended 400-element
  mesh used to validate the refinement algorithm, and integrates for
  `final time 0.5d0` s at `timestep 5.d-9` s (~100 million steps) instead of
  Test Case 23's short 16,000-step validation window.

## This case is non-evaporative by definition

Test Case 23, from which the rheology is taken, uses Yarin evaporation. This
input does not. The evaporative configuration is Test Case 25
(`../input-25/`), whose input differs from this one by exactly one line.

The split is deliberate. With evaporation enabled this configuration reaches
the code's own numerical instability check at step 1,078,685 with a stress NaN
at bead 30, roughly 1% into its target, while the non-evaporative one runs
cleanly for six million steps. Keeping the two as separate inputs makes that
controlled comparison reproducible with a one-keyword difference. The failure
occurs on the ordinary host integrator path, so it is not an accelerator
artifact.

## Refinement cadence

The dynamic-refinement cadence is Example 5's (`dynamic refinement threshold
0.4` cm, `dynamic refinement every 1.d-3` s), not Test Case 23's tighter
`0.10` cm / `1.d-5` s. The anchor spacing (`dynamic refinement anchor 0.10`
cm) is unchanged from Test 23. This choice is deliberate and load-bearing:
with the tighter cadence, growing this jet from a single bead accumulates
dozens of accepted refinement events in under 200,000 steps, each re-fitting
the cross-section radius from the previous event's own output, which compounds
a nonphysical thinning far beyond anything electrospinning genuinely produces
and eventually crashes the integrator. See
`docs/refinement-robustness-investigation.md` for the full investigation, and
`error_mod.f90` warning 108, emitted automatically when a refinement threshold
is set below twenty times the base discretization resolution.

## Running it

This is a genuinely long run: unlike Tests 21-23, it is not meant to complete
in a smoke-test timeframe. The developed bending instability carries several
hundred active beads, so the direct Coulomb cost per step is substantial.
Track the printed `n` column over time to see whether it settles into a
stationary oscillation once insertion and collector removal balance out.
