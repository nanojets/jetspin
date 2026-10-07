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

The script builds the executable in a temporary directory (Make target
`gfortran`; `gfortran-debugger` in debug mode, `gfortran-mpi` in MPI mode),
then runs shortened copies of all eight example inputs. Each simulation
executes 1,000 integration steps while retaining its original case-specific
directives. MPI mode runs Example 1 only, on two processes.

A case passes when JETSPIN closes normally, produces no error or non-finite
value (`NaN` or `Inf`), writes at least two numerical rows to `statout.dat`,
and evolves at least one observable. The suite also confirms that insertion,
dynamic refinement, rotating-field, and evaporation directives are recognized
in their corresponding examples. Original inputs and repository files are not
modified.

Serial and debug modes add the evaporation checks of `tests/evaporation/`:
`check_yarin2001.py` on Example 8 (including the diffusivity printed in its
run log), `check_kv_evaporation.py`, three 100-step Kelvin-Voigt evaporation
runs of Example 8 (`kvfluid yes`, integrators 1-3), and a 6,000-step run of
Example 5 with Kelvin-Voigt, evaporation, RK4 and a lowered refinement
threshold, which must accept a refinement and grow past 100 beads.

Set `JETSPIN_SMOKE_TIMEOUT` to change the per-run timeout from its default
of 30 seconds, and `JETSPIN_SMOKE_KEEP=1` to keep the temporary build and
run directories.
