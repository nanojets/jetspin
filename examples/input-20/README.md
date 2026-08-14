# Test Case 20: stochastic Platen Maxwell evaporation

Test Case 20 is a short, fixed-topology baseline for the evaporative Maxwell
equations integrated with the stochastic Platen scheme. It uses 1,000 beads,
the stochastic air-drag parameters of Test 4, and the evaporation parameters
used by the dynamic Maxwell cases.

Insertion, removal, and dynamic refinement are disabled deliberately. This
keeps the comparison focused on the three Platen force evaluations, the
evaporation-dependent mass scaling, and the rank-independent Gaussian block.
The case is intended to be run in serial and with two MPI ranks before the
evaporative Platen GPU port is enabled.
