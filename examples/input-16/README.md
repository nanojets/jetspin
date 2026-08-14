# Test Case 16: dynamic topology with evaporation

Test Case 16 validates the CPU and OpenACC GPU implementations of dynamic
insertion, removal, and capacity growth with the Yarin evaporation model.
It starts with a deliberately small capacity (`points 100`), enables both
`inserting yes` and `removing yes`, and uses the Maxwell RK4 integrator.

The case is intentionally serial and does not enable dynamic refinement, MPI,
or breakup. In the OpenACC build, the four Maxwell RK4 stages and evaporation
state remain on the device while insertion, removal, and capacity rebinding
are active. The 1,000-step acceptance totals are 111 additions, 122 removals,
two reallocations, and 89 final active beads. CPU and GPU trajectories may
separate after roundoff-level direct-Coulomb differences are amplified, so use
these topology totals and the saved same-build GPU result for validation.
