# Numerical smoke tests

Run the complete smoke suite from the repository root:

```sh
tests/smoke/run.sh
```

The script builds the serial executable with `gfortran` in a temporary
directory, then runs shortened copies of all seven example inputs. Each
simulation executes 1,000 integration steps while retaining its original
case-specific directives.

A case passes when JETSPIN closes normally, produces no error or non-finite
value (`NaN` or `Inf`), and writes at least two numerical rows to
`statout.dat`. Original example inputs and repository files are not modified.

Set `JETSPIN_SMOKE_TIMEOUT` to change the per-case timeout from its default
of 30 seconds.
