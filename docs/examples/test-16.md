# Test Case 16: dynamic topology with evaporation

Test Case 16 is the Maxwell reference for the GPU dynamic-topology path
with the Yarin evaporation model. It starts with 100 beads, enables both
nozzle insertion and collector removal, and keeps dynamic Akima refinement
disabled. The stored input selects RK4; changing only `integrator` to `1` or
`2` selects Euler or RK2. The small initial capacity deliberately forces host
reallocations.

All three deterministic Maxwell chains are device-resident. Euler executes
one force/stress stage, RK2 executes two, and RK4 executes four. Their
intermediate and final state updates, direct Coulomb summation, evaporation,
charge smoothing/restoration, and statistics run on the GPU. The normal build
performs no host/device state, force, or derivative-array transfer between
stages. The run takes the common device step (`source/device_step_mod.f90`)
from step 1, since it starts with 100 beads (sticky gate at `npjet` = 100);
the step is still synchronous and not yet fused (milestone M4 of
`docs/STATE.md`). Tests 17-19 take it from step 1 too.

For Euler, RK2, and RK4, NVIDIA A30 and CPU runs reproduce the same topology
totals: 111 additions, 122 removals, two reallocations, and 89 active beads
after 1,000 steps. Three-step pre-event `statout.dat` comparisons have zero
difference at `rtol=1e-12`, `atol=1e-13`; the complete written XYZ geometry is
also byte-identical. Since 2026-10-07 the 111 insertions and 122 removals
occur at the same 222 steps on both, for each integrator; before, the first insertion came at step
5 on the GPU against 4 on the CPU, and so did the insertion at the second
reallocation. The arrays are full at those insertions, and the device path
grew them and inserted one step later (not, as first assumed, an effect of
the Coulomb summation order). A full pointwise trajectory comparison is not
the acceptance criterion: the transverse components, of order `1e-12` cm,
differ by up to their own size, while `vx` agrees within `1e-5`.

A transfer audit with all development macros disabled shows only topology
decision scalars on ordinary timesteps; since 2026-10-01 they return as one
20-byte record per step. Bead records are transferred on actual
insertion/removal and output events; complete active arrays are synchronized
only for the two capacity reallocations and the final checkpoint. The
rheology-independent `nvfortran-openacc-force-oracle` target isolates complete
force stages during development. The narrower
`nvfortran-openacc-coulomb-oracle` target evaluates only the direct Coulomb
sum on the host. Both reproduce the accepted 1,000-step topology totals for
Euler, RK2, and RK4 while leaving state updates on the GPU. Their deliberate
per-stage transfers make them unsuitable for performance measurements.

Run the complete Maxwell and Kelvin--Voigt deterministic validation with:

```sh
tests/performance/dynamic/validate_evaporation.sh
```

The script builds independent NVFORTRAN CPU and standard OpenACC executables;
the GPU build does not enable either diagnostic oracle macro.

- [Input file](../../examples/input-16/input.dat)
- [Input-file notes](../../examples/input-16/README.md)
