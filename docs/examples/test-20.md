# Test Case 20: stochastic Platen Maxwell evaporation

Test Case 20 provides a fixed-topology CPU/MPI baseline for the Maxwell
evaporative equations with the stochastic Platen integrator (`system 4`,
`integrator 4`). It starts with 1,000 beads over 16 cm and runs 100 steps.

The case combines the stochastic air-drag settings from Test 4 with the
evaporation settings used by Test 16. Insertion, removal, and dynamic
refinement are disabled so that differences isolate the three Platen EOM
evaluations, evaporation-dependent mass scaling, and Gaussian-noise indexing.

Run it from the executable directory with:

```sh
cp examples/input-20/input.dat execute/input.dat
(cd execute && ./main.x)
```

For the MPI baseline, use two ranks and compare the resulting `statout.dat`
against the serial output with the numerical comparison tool. GPU support for
evaporative Platen is not enabled yet; this case therefore establishes the
CPU/MPI reference before that port.
