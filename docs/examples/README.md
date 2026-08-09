# Example simulations

JETSPIN includes eight numerical reference cases and three performance
benchmarks. Each case can be run by copying its
`input.dat` next to a compiled `main.x`. For quick validation of all cases,
use the repository smoke suite:

```sh
tests/smoke/run.sh serial
```

Test Cases 9--11 are intentionally excluded from the normal short smoke and
regression matrices: their fixed 1,000-bead direct Coulomb workloads compare
the RK4, Euler, and RK2 CPU/GPU implementations respectively.

See the individual pages in this directory for the purpose of each case.
