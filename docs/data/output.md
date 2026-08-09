# Output files

JETSPIN writes output in the directory from which `main.x` is executed.
Files are created only when their corresponding input directives are active.

At normal termination, JETSPIN also reports two performance measurements on
the terminal:

```text
Time-integration loop wall time:       0.078000 s
Time-integration throughput:          12833.333 steps/s
```

The wall-clock timer starts immediately before the temporal integration loop
and stops immediately after it. Input parsing, initial allocation, restart
loading, output-file setup, and other simulation preparation are excluded.
Work performed inside the loop—including integration, bead topology changes,
statistics, and scheduled output—is included. MPI timing is synchronized and
reports the elapsed time of the slowest participating rank.

Set `JETSPIN_PROFILE=1` to print an optional subroutine timing breakdown after
the loop timer. It includes integrator, Coulomb, bead topology, statistics,
scheduled output, and restart regions. The Coulomb region is nested within
the integrator total and must not be added to it. Profiling is disabled by
default and produces terminal output only; it does not create a profiling
file.

The profiler instruments selected code regions rather than discovering and
timing every Fortran subroutine automatically. For every region, the report
shows its cumulative wall time, its percentage of the complete temporal-loop
time, and its call count. Region times are inclusive: in particular, the
integrator measurement includes its Coulomb calls. In an MPI run the detailed
breakdown is the measurement from rank 0; it is not reduced to the maximum
time over all ranks. Use the synchronized complete-loop wall time for MPI
scaling measurements.

Clock reads add a small amount of overhead. Compare two performance runs with
profiling either enabled in both or disabled in both.

| Output | Description |
| --- | --- |
| `statdat.dat` or `statout.dat` | Time-dependent and statistical observables |
| `traj.xyz` | A time-ordered XYZ trajectory with a fixed bead count |
| `frame%06d.xyz` | Individual XYZ geometry frames |
| `frame%06d.pdb` | Individual PDB geometry frames |
| `frame%06d.psf` | Bond topology accompanying PDB output |
| `save.dat` | Binary restart and accumulated-statistics state |

Use `print list` to select terminal observables and `printstat list` to
select columns written to the statistics file. Common keys include:

| Key | Meaning |
| --- | --- |
| `t`, `ts` | Unscaled and scaled time |
| `x`, `y`, `z` | Coordinates of the bead farthest from the nozzle |
| `vx`, `vy`, `vz` | Velocity components |
| `n` | Current number of beads |
| `rc` | Jet radius at the collector |
| `curn`, `curc` | Current at the nozzle and collector |
| `mxst`, `mxsx` | Maximum stress and its axial position |

XYZ and PDB coordinates are expressed in centimetres multiplied by the
configured output rescaling factor. The full list of observables is in the
[PDF manual](../../manual/manual.pdf) and
[`manual/output.tex`](../../manual/output.tex).
