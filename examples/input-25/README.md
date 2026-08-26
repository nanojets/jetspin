# Test Case 25 input

Test Case 25 is the evaporative counterpart of Test Case 24. Every directive
is identical to `examples/input-24/input.dat` except one: `evaporation yes`
instead of `evaporation no`.

Everything else is shared with Test 24:

- Canonical JETSPIN electrostatics: `density charge 44000`,
  `collector distance 16`, `external potential 30.02076857`.
- Maxwell rheology, stochastic Platen air drag, nozzle insertion and
  collector removal from Test Case 23, with `dynamic refinement anchor 0.10`
  cm.
- Single-nozzle-bead start (`points 1`), integration to `final time 0.5d0` s
  at `timestep 5.d-9` s, and print directives from Example 3.
- Example 5's refinement cadence: `dynamic refinement threshold 0.4` cm and
  `dynamic refinement every 1.d-3` s.

## Why the pair exists

Evaporation does not currently survive this configuration, and the two inputs
are kept side by side so that the comparison stays reproducible with a
one-keyword difference:

```text
input-24 (evaporation no)   -> 6,000,000 steps, Program closed correctly,
                               zero NaN, 571 active beads, bending developed
input-25 (evaporation yes)  -> ERROR - numerical instability at nstep 1078685
                               stress NaN at bead 30, x = 12.68 cm, npjet = 262
```

That is about 1% into the 1e8-step target. The failure occurs on the ordinary
host integrator path — the failing run never engaged the persistent
accelerator path — so it is not an accelerator artifact. Bending develops
normally up to the failure: at step 800,000 the jet is at off-axis distance
1.66 cm, 27 degrees, 48 active beads.

## Status

This input **does not currently run to completion**, and that is its purpose:
it is the reproducer for an open evaporation instability. Do not use it as a
validation or performance reference until that issue is resolved.

The stress-NaN signature matches the two evaporation/refinement defects that
Test 24's investigation already fixed, each by a different route, so the first
hypothesis to test is whether the cross-event cross-section compounding
resurfaces under genuine bending, where the field being re-fit at each
accepted refinement event is far less smooth than in the configuration where
the refinement-cadence workaround was originally validated.
