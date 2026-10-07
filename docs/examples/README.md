# Example simulations

JETSPIN includes eight numerical reference cases (Examples 1-8) and
seventeen extended validation or performance cases (Test Cases 9-25). Each
case can be run by copying its `input.dat` next to a compiled `main.x`. The
smoke suite runs Examples 1-8, each shortened to 1,000 steps, and checks that
they close correctly with finite, evolving output (the `mpi` mode runs
Example 1 on two ranks):

```sh
tests/smoke/run.sh serial
```

The regression suite (`tests/regression/run.sh`) compares the same eight
examples, also shortened to 1,000 steps, with versioned baselines.

Test Cases 9-25 are not part of the smoke and regression matrices. Some
have their own runners and records: Test 9 under
`tests/performance/test9/`, Tests 10-12 and 20 under
`tests/performance/integrators/`, Test 13 and the paired CPU/GPU checks of
Tests 16 and 17 under `tests/performance/dynamic/`, and Tests 21-23 under
`tests/refinement/`; Tests 12, 13, 16, 17, 24 and 25 are also cases of the
restart check (`tests/restart/run.sh`). The others are run by hand as their
pages describe.

Test Cases 9-12 run fixed 1,000-bead direct Coulomb workloads that compare
the RK4, Euler, RK2, and stochastic Platen CPU/GPU implementations.
Test Case 13 is the separate dynamic-topology baseline with active nozzle
insertion and collector removal. Test Case 14 extends it to 1,500 initial
beads with bounded persistent capacity.
Test Case 15 forces `reallocate_jet` with insertion and without removal; in
the OpenACC build it checks the device capacity growth, which since
2026-10-07 inserts the bead in the same step as the CPU build.
Test Case 16 validates the complete device-resident Maxwell Euler, RK2, and
RK4 chains with dynamic insertion, removal, reallocation, and evaporation
enabled.
Test Case 17 validates the complete device-resident Kelvin–Voigt Euler, RK2,
and RK4 chains with the same dynamic evaporation and capacity-growth workload.
Test Case 18 is the high-resolution Maxwell evaporation performance probe
with 800 points over 16 cm.
Test Case 19 is the corresponding high-resolution Kelvin–Voigt evaporation
performance probe.
Test Case 20 is the fixed-topology stochastic Platen Maxwell evaporation
CPU/GPU validation case.
Test Case 21 is the anchored dynamic-refinement Maxwell/Platen evaporation
validation for the persistent OpenACC path. Akima coefficient construction and
field interpolation execute on the GPU at an accepted remeshing event; target
mesh construction and conservation remain host event work.
Test Case 22 repeats that event three times and deliberately forces a capacity
increase at every remesh to stress persistent-data release/rebind.
Test Case 23 shortens the collector distance so device-side removal occurs
both before and after the final remesh.
Test Case 24 is a long production run, about 100 million steps, that grows a
jet from a single nozzle bead until insertion and collector removal balance
into a stationary bead count. It is non-evaporative.
Test Case 25 is Test Case 24 with Yarin evaporation enabled and a 50 %
initial polymer fraction; with Yarin's 6 % the dried jet becomes unstable, as
analysed on its page. It is equally long.

## Pages

- [Test Case 1: one-dimensional reference case](test-1.md)
- [Test Case 2: one-dimensional bead insertion](test-2.md)
- [Test Case 3: three-dimensional PVP electrospinning](test-3.md)
- [Test Case 4: gas counterflow](test-4.md)
- [Test Case 5: dynamic refinement](test-5.md)
- [Test Case 6: Kelvin–Voigt fluid](test-6.md)
- [Test Case 7: evaporation and rotating electric field](test-7.md)
- [Test Case 8: Yarin 2001 reference parameters](test-8.md)
- [Test Case 9: 1,000-bead CPU/GPU benchmark](test-9.md)
- [Test Case 10: Euler CPU/GPU benchmark](test-10.md)
- [Test Case 11: RK2 CPU/GPU benchmark](test-11.md)
- [Test Case 12: Platen stochastic CPU/GPU benchmark](test-12.md)
- [Test Case 13: dynamic-topology CPU/GPU benchmark](test-13.md)
- [Test Case 14: larger bounded dynamic topology](test-14.md)
- [Test Case 15: forced array reallocation](test-15.md)
- [Test Case 16: dynamic topology with evaporation](test-16.md)
- [Test Case 17: Kelvin–Voigt dynamic topology with evaporation](test-17.md)
- [Test Case 18: high-resolution dynamic Maxwell evaporation](test-18.md)
- [Test Case 19: high-resolution dynamic Kelvin–Voigt evaporation](test-19.md)
- [Test Case 20: stochastic Platen Maxwell evaporation](test-20.md)
- [Test Case 21: anchored dynamic refinement with Platen evaporation](test-21.md)
- [Test Case 22: repeated dynamic refinement and capacity growth](test-22.md)
- [Test Case 23: refinement interleaved with collector removal](test-23.md)
- [Test Case 24: long run to a stationary bead count](test-24.md)
- [Test Case 25: evaporative counterpart of Test 24](test-25.md)
