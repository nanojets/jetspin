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

See [Input files](../data/input.md), [Output files](../data/output.md), and
the corresponding [LaTeX section](../../manual/running.tex).
