# Test Case 17: Kelvin–Voigt dynamic topology with evaporation

Test Case 17 is the Kelvin–Voigt counterpart of Test Case 16. It uses the
same 100-bead, insertion/removal, evaporation, and forced-reallocation setup,
with `kvfluid yes` selecting the Kelvin–Voigt evaporation integrator.

The stored input selects RK4. Setting `integrator 1` or `integrator 2` creates
the corresponding Euler or RK2 validation without changing any other
parameter. The `nvfortran-openacc` path keeps all force evaluations,
intermediate/final updates, and Kelvin–Voigt evaporation stress on the GPU for
all three deterministic integrators. Jet and derivative arrays remain resident
between stages. Ordinary timesteps exchange only the small insertion/removal
control state; bead records move for actual topology, output, checkpoint, or
capacity-growth events.

The accepted 1,000-step CPU and A30 executions report 111 additions, 122
removals, two reallocations, and 89 final active beads for Euler, RK2, and RK4.
Their trajectories are not required to agree point by point after topology
thresholds are crossed:
the target-centric GPU Coulomb sum has a different floating-point order and
the bending instability amplifies roundoff-level perturbations.

The historical CPU Kelvin–Voigt evaporation equation does not add the air-drag
or lift terms even though this input retains `airdrag yes`; the GPU path
deliberately preserves that established behaviour.
