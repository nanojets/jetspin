# Test Case 25: evaporative counterpart of Test 24

Test Case 25 is Test Case 24 with evaporation enabled, and nothing else
changed. The two `input.dat` files differ by exactly one line:

```text
input-24:  evaporation no
input-25:  evaporation yes
```

Everything else is shared: canonical JETSPIN electrostatics
(`density charge 44000`, `collector distance 16`,
`external potential 30.02076857`), Maxwell rheology, stochastic Platen air
drag, nozzle insertion and collector removal from Test 23, a single-nozzle-bead
start from Example 3, integration to `final time 0.5d0` s at
`timestep 5.d-9` s, and Example 5's refinement cadence
(`dynamic refinement threshold 0.4` cm, `dynamic refinement every 1.d-3` s).

## This case does not currently complete

That is the point of it. Test 24 was originally meant to carry Test 23's Yarin
evaporation along with its rheology, but evaporation does not survive the
corrected electrostatics. Rather than silently dropping it, the two
configurations are kept as separate inputs so the failure stays reproducible
with a single-keyword difference:

| | Test 24 (`evaporation no`) | Test 25 (`evaporation yes`) |
| --- | --- | --- |
| outcome | 6,000,000 steps, `Program closed correctly` | `ERROR - numerical instability` |
| stopped at | — | step 1,078,685 |
| failure | none, zero NaN | stress NaN at bead 30, x = 12.68 cm |
| beads at stop | 571 (oscillating, stationary) | 262 |

Step 1,078,685 is about 1% into the 1e8-step target.

Two facts narrow the search. The failing run never engaged the persistent
accelerator path — `persistent=F` at every probe point — so this is the
ordinary host integrator and not an accelerator artifact. And the bending
instability is developing normally right up to the failure: at step 800,000
the jet is at off-axis distance 1.66 cm, bending angle 27 degrees, with 48
active beads.

The stress-NaN signature is the same one that the `cp`/`evlim` defect and the
cross-section thinning defect each produced by different routes during Test
24's investigation. The first hypothesis to test is therefore whether the
cross-event compounding of the cross-section radius resurfaces under genuine
bending: the refinement-cadence workaround that contains it was validated on a
configuration whose cross-section field was far smoother than it is here.

## Use

Use this case as the reproducer for the open evaporation instability described
above. It is not a validation or performance reference and is excluded from
the smoke and regression matrices.

- [Input file](../../examples/input-25/input.dat)
- [Test 24, the non-evaporative twin](test-24.md)
- [Refinement robustness investigation](../refinement-robustness-investigation.md)
- [Dynamic refinement](../introduction/dynamic-refinement.md)
