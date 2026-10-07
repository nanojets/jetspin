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

In the OpenACC build the run takes the common device step from the first
step (100 beads) and exercises the device-to-host synchronization, capacity
growth, and device remapping path used by the larger dynamic benchmarks, with
a deliberately small capacity. The acceptance check (completion, topology
additions, positive reallocation count) does not require bead-for-bead
identity, although the two builds now agree row for row (below).

That relaxed criterion was originally attributed to the high insertion speed
amplifying GPU rounding. That attribution was wrong. The real cause was
`reset_coulomb_accelerator` clearing its bookkeeping flags without issuing the
matching `exit data delete`, so a reallocation left the device holding a
mapping to host storage that no longer existed. The GPU run produced NaN from
its first printed line with the bead count frozen at 103, against 211 beads
and `x = 35.99` on CPU. After the fix the two runs had the same bead count
and `x` equal to about eight significant digits (35.99200794 and 35.99201023
cm); that comparison predates the collector rule below, and since the rule
all rows are identical. This test is in fact a regression guard for that defect and
its acceptance criterion could be tightened. See
[OpenACC](../introduction/openacc.md).

On an NVIDIA A30 (NVFORTRAN 24.3) both builds make 111 additions and two
reallocations and finish with 211 beads, every insertion at the same step,
and all eleven statistics rows are identical (within `4.3e-5` before the
collector rule below; up to 14 percent in the transverse velocity before the
2026-10-01 nozzle-insertion fix).

The test does not remove beads. Since 2026-10-07 a bead that reaches the
collector is frozen there in any case (it discharges on the grounded
electrode), so the leading bead stays at `x = 16` cm. Before, a run without
removal froze nothing: the leading bead crossed the collector plane keeping
its charge and its Coulomb interactions, and reached `x = 35.99` cm in 1,000
steps.
Until 2026-10-07 two insertions fell one step later on the A30, the first
one (step 5 against 4) and the one at the second reallocation (906 against
905): when the capacity was exhausted, the device path grew it and inserted
at the next step, whereas the CPU reallocates and inserts in the same step.
The device path now decides the step again after the growth.

- [Input file](../../examples/input-15/input.dat)
