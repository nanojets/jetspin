# Compiling JETSPIN

JETSPIN requires a Fortran compiler. Parallel builds also require an MPI
implementation. Run all commands below from the repository root; the
Makefile does not need to be copied into `source/`.

## Serial build with GFortran

```sh
make -C source -f ../build/Makefile gfortran
```

Object and module files are created in `source/`. The resulting executable
is written to `execute/main.x`.

## OpenMPI build with GFortran

```sh
make -C source -f ../build/Makefile gfortran-mpi
```

Run `mpif90 --showme:command` to confirm that the MPI wrapper uses GFortran
when selecting this target.

## Other targets

| Target | Purpose |
| --- | --- |
| `gfortran` | Optimized serial GFortran build |
| `gfortran-mpi` | Optimized OpenMPI/GFortran build |
| `gfortran-debugger` | Serial build with runtime checks |
| `gfortran-mpidebugger` | MPI build with runtime checks |
| `intel` | Optimized serial Intel Fortran build |
| `intel-mpi` | Intel MPI build |
| `intel-openmpi` | Intel Fortran with OpenMPI |
| `cygwin`, `cygwin-mpi` | Windows/Cygwin builds |
| `help` | Display available targets |
| `clean` | Remove objects and module files from `source/` |

For example, clean a build with:

```sh
make -C source -f ../build/Makefile clean
```

The complete build rules are in [`build/Makefile`](../../build/Makefile).
See also the corresponding [LaTeX section](../../manual/compiling.tex).
