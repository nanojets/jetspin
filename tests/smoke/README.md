# Numerical smoke tests

Run the complete smoke suite from the repository root:

```sh
tests/smoke/run.sh
```

Runtime-checking and MPI variants are also available:

```sh
tests/smoke/run.sh debug
tests/smoke/run.sh mpi
```

The script builds the serial executable with `gfortran` in a temporary
directory, then runs shortened copies of all eight example inputs. Each
simulation executes 1,000 integration steps while retaining its original
case-specific directives.

A case passes when JETSPIN closes normally, produces no error or non-finite
value (`NaN` or `Inf`), writes at least two numerical rows to `statout.dat`,
and evolves at least one observable. The suite also confirms that insertion,
dynamic refinement, rotating-field, and evaporation directives are recognized
in their corresponding examples. Original inputs and repository files are not
modified. MPI mode runs one representative case on two processes.

Set `JETSPIN_SMOKE_TIMEOUT` to change the per-case timeout from its default
of 30 seconds.
