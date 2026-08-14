# Test Case 17: Kelvin–Voigt dynamic topology with evaporation

Test Case 17 repeats Test Case 16 with `kvfluid yes`, using the same 100-bead
capacity-stress configuration and the same insertion/removal and evaporation
parameters. Its stored input selects RK4; changing only `integrator` to `1` or
`2` exercises the equivalent Euler or RK2 path. All three deterministic
Kelvin–Voigt evaporation integrators are device-resident for this configuration.

On the GPU, Euler executes one force/stress stage, RK2 executes two, and RK4
executes four. Their intermediate and final updates, charge
smoothing/restoration, direct Coulomb force, evaporation, and final statistics
reductions execute without transferring jet or derivative arrays between
stages. Per-step host/device traffic is limited to topology decision scalars.
Actual insertions, removals, scheduled output, capacity growth, and the final
checkpoint transfer only their required data; the two reallocations
necessarily synchronize and rebind the active state.

For Euler, RK2, and RK4, the accepted 1,000-step CPU and NVIDIA A30 runs all
produce 111 additions, 122 removals, two reallocations, and 89 active beads.
Three-step pre-event CPU/GPU comparisons are identical at `rtol=1e-12` and
`atol=1e-13`. The first insertion threshold is crossed at step 4 on the CPU
and step 5 on the GPU. Pointwise trajectories then diverge because the
race-free target-centric GPU Coulomb sum uses a different accumulation order
and the bending instability amplifies the roundoff-level perturbation. Equal
topology totals, pre-event agreement, clean completion, and the absence of
per-stage array transfers are therefore the acceptance criteria.

For numerical isolation only, the
`nvfortran-openacc-kv-host-forces` target evaluates every Kelvin–Voigt force
stage through the trusted CPU equations and uploads the derivatives to the
device. It must not be used for performance measurements. The standard target
uses only device force kernels.

The established CPU Kelvin–Voigt evaporation equation omits aerodynamic drag
and lift even when `airdrag yes` appears in the input. The accelerator follows
the same semantics rather than silently changing the physical model.

The paired CPU/GPU validation for all three deterministic integrators is
automated by:

```sh
tests/performance/dynamic/validate_evaporation.sh
```

- [Input file](../../examples/input-17/input.dat)
- [Input-file notes](../../examples/input-17/README.md)
