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
`final time`. On restart `statout.dat` and `traj.xyz` are appended, while
the other outputs, including `save.dat`, are overwritten. Stochastic runs
resume from the saved physical state but not from the saved random-generator
state (see [random numbers](random-numbers.md#restart-limitation)).

See [Input files](../data/input.md), [Output files](../data/output.md), and
the corresponding [LaTeX section](../../manual/running.tex).
