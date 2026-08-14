# Example simulations

JETSPIN includes eight numerical reference cases and five performance
benchmarks. Each case can be run by copying its
`input.dat` next to a compiled `main.x`. For quick validation of all cases,
use the repository smoke suite:

```sh
tests/smoke/run.sh serial
```

Test Cases 9--12 are intentionally excluded from the normal short smoke and
regression matrices: their fixed 1,000-bead direct Coulomb workloads compare
the RK4, Euler, RK2, and stochastic Platen CPU/GPU implementations.
Test Case 13 is the separate dynamic-topology baseline with active nozzle
insertion and collector removal. Test Case 14 extends it to 1,500 initial
beads with bounded persistent capacity.
Test Case 15 is a CPU-side forced `reallocate_jet` check.
Test Case 16 validates the complete device-resident Maxwell Euler, RK2, and
RK4 chains with dynamic insertion, removal, reallocation, and evaporation
enabled.
Test Case 17 validates the complete device-resident Kelvin–Voigt Euler, RK2,
and RK4 chains with the same dynamic evaporation and capacity-growth workload.
Test Case 18 is the high-resolution Maxwell evaporation performance probe
with 800 points over 16 cm.
Test Case 19 is the corresponding high-resolution Kelvin–Voigt evaporation
performance probe.
Test Case 20 is a fixed-topology stochastic Platen Maxwell evaporation
baseline for serial/MPI verification.

See the individual pages in this directory for the purpose of each case.
