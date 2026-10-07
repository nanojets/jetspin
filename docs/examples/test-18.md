# Test Case 18: high-resolution dynamic Maxwell evaporation

Test Case 18 is the high-resolution counterpart of Test 16. It uses the
Maxwell RK4 integrator with evaporation, nozzle insertion, and collector
removal enabled. The initial jet length is 16 cm and `points 800`, giving a
nominal spacing of 0.02 cm (200 micrometres) between adjacent nodes. As in
Test 16, dynamic refinement is disabled and the case is serial.

Because the initial capacity already contains 800 intervals, this case is
intended to measure the dynamic topology path at a finer spatial resolution,
not to force frequent capacity reallocations. Record topology events and the
final numerical state before using it as a GPU regression baseline.

- [Input file](../../examples/input-18/input.dat)
- [Test 16 comparison](test-16.md)

The NVFORTRAN CPU reference completed in about 14.2 s with 500 additions,
899 removals, three reallocations, and 401 active beads. Since 2026-10-07 the
NVFORTRAN 24.3/A30 standard run gives the same totals, every insertion at
the CPU's step; one removal at step 809 is decided in a different event
line, and the transverse components, of order `1e-12` cm, differ by up to
their own size. Before, the device path inserted one step later whenever
the arrays had to grow, and the A30 run gave 499 additions, 898 removals and
401 active beads (both oracles 499, 899 and 400). The case remains a
performance/porting probe rather than a strict pointwise regression case.
