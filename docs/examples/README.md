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

See the individual pages in this directory for the purpose of each case.
