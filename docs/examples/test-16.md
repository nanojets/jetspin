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
stages.

For Euler, RK2, and RK4, NVIDIA A30 and CPU runs reproduce the same topology
totals: 111 additions, 122 removals, two reallocations, and 89 active beads
after 1,000 steps. Three-step pre-event `statout.dat` comparisons have zero
difference at `rtol=1e-12`, `atol=1e-13`; the complete written XYZ geometry is
also byte-identical. The first insertion occurs at step 4 on the CPU and step
5 on the GPU. A full pointwise trajectory comparison after that threshold is
intentionally not the acceptance criterion: the different direct-Coulomb
accumulation order creates roundoff-level changes that are amplified by the
physical bending instability and later topology crossings.

A transfer audit with all development macros disabled shows only topology
decision scalars on ordinary timesteps. Bead records are transferred on actual
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
