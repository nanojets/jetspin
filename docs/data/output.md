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
reports the elapsed time of the slowest participating rank. The throughput
divides the steps of this run by that time (since 2026-10-07; before, a
restarted run also counted the steps made before the restart). The clocks
use 64-bit counts: NVFORTRAN builds wrapped them after 2147 s until
2026-10-05.

Set `JETSPIN_PROFILE=1` to print an optional subroutine timing breakdown after
the loop timer. It includes integrator, Coulomb, bead topology (add, remove,
erase bead), breakup check, statistics, scheduled output, and restart
regions. The Coulomb region is nested within
the integrator total and must not be added to it. Profiling is disabled by
default and produces terminal output only; it does not create a profiling
file.

The profiler instruments selected code regions rather than discovering and
timing every Fortran subroutine automatically. For every region, the report
shows its cumulative wall time, its percentage of the complete temporal-loop
time, and its call count. Region times are inclusive: in particular, the
integrator measurement includes its Coulomb, EOM-evaluation, and RK-update
regions. The EOM-evaluation and RK-update timers exist only in the host
Euler, RK2 and RK4 integrators without evaporation: they read zero calls in
Platen, evaporative, Kelvin–Voigt and device-step runs. In an MPI run the detailed breakdown is the measurement from rank 0;
it is not reduced to the maximum time over all ranks. Use the synchronized
complete-loop wall time for MPI scaling measurements.

Clock reads add a small amount of overhead. Compare two performance runs with
profiling either enabled in both or disabled in both.

| Output | Description |
| --- | --- |
| `statout.dat` | Time-dependent and statistical observables selected with `printstat list` (text: in a new run a record at step 0, then one per `print time`; each record is the step, format `i10`, then one `g20.10` column per key, under a three-line `#` header) |
| `statdat.dat` | Developer-only binary dump of all observables, written only with `printstat binary` |
| `traj.xyz` | A time-ordered XYZ trajectory with a fixed bead count (`print xyz maxnum`, from the leading bead; missing beads at the origin) |
| `frame%06d.xyz` | Individual XYZ geometry frames, with all the beads |
| `frame%06d.pdb` | Individual PDB geometry frames, with all the beads; numbered by step divided by the frame interval, so the first frame of a new run is `frame000001` (a restarted run continues the numbering) |
| `frame%06d.psf` | Bond topology accompanying PDB output (`vmd -psf frame000001.psf -pdb frame000001.pdb`) |
| `save.dat` | Binary restart state, written every `restart dump` steps and at the end of the run: bead records in double precision and the run state (pool cursor, random generator, refinement counters, device gate); see [running](../introduction/running.md) |
| `traj.dat`, `bead.dat` | Developer-only binary trajectory (`print binary`, reduced units after a parameter header) and record of the removed beads (`print binary removed`, with `removing yes`). [`tools/traj2dcd.f90`](../../tools/traj2dcd.f90) converts `traj.dat` (any style) to `trajout.xyz` and `trajout.dcd` in cm, for `vmd trajout.xyz -dcd trajout.dcd` |
| `topology-state.dat` | With `JETSPIN_TOPOLOGY_SNAPSHOT=1`: at every insertion or removal a line `event step inpjet npjet active add`, then one row per active bead (index, frozen flag, position, stress, velocity, mass, charge, volume; with evaporation also the evaporated volume and the concentration) |

A restarted run appends to `statout.dat` (after a new three-line `#` header
block), `traj.xyz` and the developer binaries (`statdat.dat`, `traj.dat`,
`bead.dat`), creating them if they are missing.

Use `print list` to select terminal observables and `printstat list` to
select columns written to the statistics file. Common keys include:

| Key | Meaning |
| --- | --- |
| `t`, `ts` | Unscaled and scaled time |
| `x`, `y`, `z` | Coordinates of the leading bead (`inpjet`), the bead farthest from the nozzle along the jet |
| `vx`, `vy`, `vz` | Velocity components of the leading bead |
| `n` | Current number of beads |
| `rc` | Jet radius at the collector |
| `curn`, `curc` | Current at the nozzle and collector |
| `mxst`, `mxsx` | Maximum stress along the jet during the print interval, and its x position |
| `nref` | Number of accepted dynamic-refinement events |
| `angl` | Bending (cone) angle of the farthest bead, in degrees |
| `cpu`, `cpur`, `cpue` | Processor time (`cpu_time` of rank 0, not the wall clock of `job time`): of the last print interval (the first record: since the program started); estimated remaining, `cpu` times all the records from t=0 minus those of this run (it overestimates after a restart); elapsed since the first record of this run |
| `nms` | Mean number of steps per neighbour-list build during the print interval (zero without multiple step, or when the list was not built in the interval; until 2026-10-07 the latter case divided by zero) |
| `erms` | Once per multiple-step interval (last step before a scheduled update) the extrapolated outer Coulomb forces are compared with those of the direct sum; the modulus of the difference of the outer Coulomb accelerations, averaged over the beads and over the checks of the print interval, in cm s^-2 (until 2026-10-07 a reduced acceleration times `chargescale**2/lengthscale**2`, labelled `erms (dyne)`). Computed only when `erms` is selected; zero without multiple step. Until 2026-10-07 the check step used the direct sum as force when `erms` was printed, and stale forces otherwise; now it does not change the forces |

A bead that reaches the collector is frozen there at x = h and stays the
leading bead until it is removed: once the jet has reached the collector,
the leading-bead keys describe the first deposited bead without
`removing yes`, and the last bead that reached the collector with it. The
collector keys (`curc`, `vc`, `svc`, `mfc`, `rc`, `rrr`, `visc`, `gc`,
`evrc`, `emfc`) are accumulated only when beads are removed, so they are
zero in a print interval without removals and always zero without
`removing yes`; `curn` and `mfn` are accumulated at each insertion and need
`inserting yes`.

## Terminal log lines

Besides the observables of `print list`, rank 0 prints these lines:

| Line | When |
| --- | --- |
| `Gaussian history pool: values=V covers S steps at N beads` | Before the loop, in a run that reads the sequential Gaussian pool ([random numbers](../introduction/random-numbers.md)) |
| `OpenACC device step engaged (scheme, model) at step N with M active beads` | OpenACC build, once, when the device step engages; `resumed (...)` after a restart of an engaged run ([OpenACC](../introduction/openacc.md)) |
| `Topology event: step=N add=T/F remove=R active=M` | Every step with an insertion or a removal, in every build |
| `Dynamic refinement event: ...` and the `geometry check`, `thin-bead classification`, `geometry range`, `invariants`, `anchor fields`, `conserved amounts` lines | Every accepted refinement event ([dynamic refinement](../introduction/dynamic-refinement.md)) |
| `Akima endpoint diagnostic: field=...` | Seven lines per accepted refinement event (endpoint slopes and extrema of the fitted fields) |
| `Dynamic refinement capacity: old= new=`, `OpenACC refinement capacity rebind: old= new=` | An accepted event that grows the arrays (the second in OpenACC builds) |
| `Initial dynamic-refinement anchors: ...` | Start of a refining run |
| `Numerical instability ...` | Before error 14, with the failing bead, its neighbour and the context |
| `Rattao blowup diagnostic: ...` | Evaporative runs, for every bead and stage whose relaxation-time ratio leaves [1e-2, 1e2] (a sign of an unstable run) |
| `Topology additions:`, `Topology removals:`, `Array reallocations:`, `Topology active beads:` | End of the run; the totals include the steps before a restart, since they are restored from `save.dat` |
| `Program closed correctly` | Normal end |
| `CPU time = T seconds on N CPU` (or `CPUs`) | Normal end: processor time of rank 0 since the program started (not the wall clock), N the number of MPI ranks |

The `print list` records have the `statout.dat` format (step `i10`, then
`g20.10` columns); their three-line `#` header is repeated every 100
terminal records.

XYZ and PDB coordinates (`traj.xyz`, `frame%06d.xyz`, `frame%06d.pdb`) are
expressed in centimetres multiplied by `print xyz rescalexyz` or
`print pdb rescalepdb` (default 1). Until 2026-10-07 `traj.xyz` was in
reduced (internal) length units unless `rescalexyz` was given, and the frame
files applied the length unit twice when a rescale factor was given. The
coordinates are written with `f10.5` (XYZ) and `f8.2` (PDB): a scale factor
may be needed to resolve a thin jet in the PDB files. The full list of observables is in the
[PDF manual](../../manual/manual.pdf) and
[`manual/output.tex`](../../manual/output.tex).
