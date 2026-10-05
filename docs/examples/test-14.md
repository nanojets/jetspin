# Test Case 14: larger bounded dynamic topology

Test 14 starts with 1,500 beads and reserves 1,756 slots before persistent
OpenACC mapping. Insertion, removal, RK4, Coulomb, and EOM remain device
resident; no reallocation occurs during the run. It stresses the bounded
capacity path and serves as the no-overflow control for Test 15's explicit
teardown and remapping case.

On an NVIDIA A30 (NVFORTRAN 24.3) the run reproduces the CPU topology stream
step for step: 37 events, 18 additions, 19 removals, 1,499 final active beads.
Nine of the eleven statistics rows agree with the CPU within `3e-4`
relatively; the largest difference, `4.8e-3`, is in the transverse velocity.
Before the 2026-10-01 fix of the nozzle charge smoothing and of the placement
of the blocked bead on the device, the A30 made 19 additions and finished
with 1,500 beads.

- [Input file](../../examples/input-14/input.dat)
