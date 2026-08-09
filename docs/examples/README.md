# Example simulations

JETSPIN includes eight numerical reference cases and one performance
benchmark. Each case can be run by copying its
`input.dat` next to a compiled `main.x`. For quick validation of all cases,
use the repository smoke suite:

```sh
tests/smoke/run.sh serial
```

Test Case 9 is intentionally excluded from the normal short smoke and
regression matrices: its fixed 1,000-bead direct Coulomb workload is designed
for CPU/GPU performance measurements.

See the individual pages in this directory for the purpose of each case.
