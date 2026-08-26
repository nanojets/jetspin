# Test Case 15: forced array reallocation

Test 15 is a dynamic allocation test. It starts with 100 beads and enables
insertion while the initial allocation has no spare capacity. The first
successful insertion exceeds the allocation and must call `reallocate_jet`.
The nozzle velocity is set to 200,000 cm/s to produce repeated insertions in
the short run.
The program reports `Array reallocations:` at shutdown; the test acceptance
condition is a value greater than zero and at least one topology addition.
Each subsequent reallocation grows the capacity by 100 slots, while the
initial allocation remains governed by the normal 100-slot increment.

The test enters the small persistent GPU topology gate and exercises the same
device-to-host synchronization, capacity growth, and device remapping path
used by the larger dynamic benchmark, but with a deliberately small capacity.
CPU and GPU trajectories are not required to match bead-for-bead; acceptance
checks the successful completion, topology additions, and positive
reallocation count.

That relaxed criterion was originally attributed to the high insertion speed
amplifying GPU rounding. That attribution was wrong. The real cause was
`reset_coulomb_accelerator` clearing its bookkeeping flags without issuing the
matching `exit data delete`, so a reallocation left the device holding a
mapping to host storage that no longer existed. The GPU run produced NaN from
its first printed line with the bead count frozen at 103, against 211 beads
and `x = 35.99` on CPU. After the fix the two agree to eight significant
digits, so this test is in fact a regression guard for that defect and its
acceptance criterion could be tightened. See
[OpenACC](../introduction/openacc.md).

- [Input file](../../examples/input-15/input.dat)
