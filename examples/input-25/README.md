# Test Case 25 input

Test Case 25 is the evaporative counterpart of Test Case 24. Every directive
is identical to `examples/input-24/input.dat` except two: `evaporation yes`
instead of `evaporation no`, and `evaporation polymer frac 0.50d0` instead of
the unused `0.06d0`.

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

## Evaporation parameters

The Yarin (2001) law is used unchanged: `bconstant 7`, `mconstant 0.1`,
`tconstant 1`, relative humidity 0.165 at 293.15 K. The initial polymer mass
fraction is 0.50 instead of Yarin's 0.06. At the solidification cutoff
(`cp = 0.9`, reached by every collected bead) the element keeps 56 % of its
volume, the viscosity grows 2.49 times, the relaxation time 1.80 times, and
the elastic modulus 1.38 times.

With the 6 % fraction the dried jet keeps its charge while losing fifteen
times its mass and cross-section; its head is ejected laterally and the run
stops with a stress NaN near step 1.07 million. The analysis, the diagnostic
variants, and the choice of 0.50 are in `docs/examples/test-25.md`.

## Running it

Like Test 24, this is a long run and is not meant to complete in a
smoke-test timeframe. The reference run, truncated at 5 million steps,
reaches the collector at step 2.1 million and then holds a stationary regime
with about 270 active beads, a 20° cone, and a 2.8 µm fibre radius at the
collector.
