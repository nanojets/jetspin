# Test Case 14: larger bounded dynamic topology

Test 14 starts with 1,500 beads and reserves 1,756 slots before persistent
OpenACC mapping. Insertion, removal, RK4, Coulomb, and EOM remain device
resident; no reallocation occurs during the run. It stresses the bounded
capacity path before device teardown and remapping are implemented for true
capacity growth.

- [Input file](../../examples/input-14/input.dat)
