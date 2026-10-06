# Running JETSPIN

JETSPIN reads a case-insensitive file named `input.dat` from its current
working directory. Output files are written to that same directory.

After compiling, copy an example input into `execute/`:

```sh
cp examples/input-1/input.dat execute/input.dat
cd execute
```

Run a serial executable with:

```sh
./main.x
```

Run an MPI executable on four processes with:

```sh
mpirun -np 4 ./main.x
```

The program prints selected observables to the terminal and writes only the
outputs requested by directives in `input.dat`. Restart state is stored in
the binary `save.dat` file.

To restart, copy `save.dat` to `restart.dat` in the run directory and add
`restart yes` to `input.dat`. The run resumes from the saved step and time
and continues up to `final time`; to extend a finished run, increase
`final time`. On restart `statout.dat` and `traj.xyz` are appended (and
created if missing), while the other outputs, including `save.dat`, are
overwritten. `save.dat` is written every `restart dump` steps (100,000 by
default) and at the end of the run.

Since 2026-10-06 `save.dat` holds the bead state in double precision,
followed by the run state that the beads do not carry: the position in the
Gaussian pool, the state of the random generator, and the counters of the
dynamic refinement and of the anchor tagging (see
[random numbers](random-numbers.md#restart)). With the same `input.dat` and
executable, a restarted run continues the uninterrupted one exactly, with or
without evaporation; `tests/restart/run.sh [nvfortran|openacc]` checks this
for Tests 24 and 25. The statistics accumulated since the last print are not
saved, so the printed rows are all reproduced when `restart dump` is a
multiple of the print interval; `restart reset` has no further effect.
Restart files written before that date are still read, in single precision
and without the run state (the noise then restarts from the beginning of the
pool), and a run with evaporation could not be restarted at all
("restart file is corrupted").

See [Input files](../data/input.md), [Output files](../data/output.md), and
the corresponding [LaTeX section](../../manual/running.tex).
