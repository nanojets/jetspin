# Test Case 14: larger bounded dynamic topology

Test 14 starts with 1,500 beads and reserves 1,756 slots before persistent
OpenACC mapping. Insertion, removal, RK4, Coulomb, and EOM remain device
resident; no reallocation occurs during the run. It stresses the bounded
capacity path and serves as the no-overflow control for Test 15's explicit
teardown and remapping case.

- [Input file](../../examples/input-14/input.dat)
