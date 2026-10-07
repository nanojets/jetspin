# Exact restart check

`run.sh` checks that a restarted run continues the uninterrupted one exactly.
For each case it runs three copies of the case's input, with `restart dump`
at the restart step and shorter print and trajectory intervals:

- run A stops at the restart step and writes `save.dat`;
- run B reads it as `restart.dat` (`restart yes`) and continues;
- run C runs the same steps without interruption.

`check_restart.py` then requires every `statout.dat` row, every terminal row
and every `traj.xyz` frame written by B after the restart step to equal C's,
character for character.

| Case | Input | Steps (restart at) | What it exercises |
| --- | --- | --- | --- |
| `t24`, `t25` | Tests 24, 25 | 200,000 (100,000) | Platen with insertion, removal and refinement, without and with evaporation: pool cursor, refinement and anchor counters |
| `t12` | Test 12 | 1,000 (500) | Platen on a fixed 1,000-bead jet |
| `t13` | Test 13 | 1,000 (500) | RK4 with insertion and removal |
| `t16` | Test 16 | 1,000 (500) | RK4 with evaporation |
| `t17` | Test 17 | 1,000 (500) | RK4, Kelvin-Voigt with evaporation |
| `ex4` | Example 4 | 20,000 (10,000) | Platen with insertion, without refinement, from one bead |
| `ex8` | Example 8 | 20,000 (14,000) | RK4 with evaporation from one bead, restarted after a compaction |

The final times of `t24` and `t25` are exact multiples of the timestep
(`5.d-4` and `1.d-3` with `5.d-9`); the others are half a step short of the
step count (for instance `4.995d-6` with a timestep of `1.d-8` for 500
steps). Since 2026-10-07 a run ends at the first step whose time reaches the
final time to within a relative `1e-12` (`integration_last_step`), so both
forms give the listed counts with every compiler. Test 12 sizes its pool by
the final time, which also checks that a fixed jet extended by a restart
resumes its pool.

In the OpenACC build Tests 12, 13, 16 and 17 run on the device step from
the first step, so the restart maps a saved state onto the device; Tests 24
and 25 and Example 4 stay below 100 beads, on the CPU build's code. Example 8
engages the device step at step 13,699, when the arrays reach index 100
with 89 beads active; at step 13,840 an insertion compacts the arrays to 91
beads, below the gate, which stays open, and they reach 100 again at step
15,109. The restart at step 14,000 falls in that window: the restart file
records that the gate had opened (restart state version 2), so the
restarted run resumes on the device as the uninterrupted run does; with a
version 1 file it continued on the host up to step 15,109.

```sh
tests/restart/run.sh nvfortran          # NVFORTRAN CPU build, all cases
tests/restart/run.sh openacc            # OpenACC build
tests/restart/run.sh gfortran           # GFortran
tests/restart/run.sh gfortran-mpi       # two ranks ($MPIEXEC, $JETSPIN_MPIFC)
tests/restart/run.sh nvfortran-mpi      # two ranks, NVHPC MPI
tests/restart/run.sh openacc ex8 t17    # selected cases
```

Only `GPUCC` (default `80`) is passed to `make`; `CUDA_VERSION` comes from
the environment through the Makefile default (`12.3`), so export
`CUDA_VERSION=12.9` before building with NVHPC 25.5.

Each backend takes about two to three minutes for all cases (the MPI ones
run B and C one after the other, within the two ranks).

`JETSPIN_RESTART_KEEP=1` keeps the work directory, `JETSPIN_RESTART_TIMEOUT`
(default 600 s) bounds each run and `JETSPIN_RESTART_EXE` uses an existing
executable instead of building one. The device path of Tests 24 and 25 (from
about 1.4 million steps) was checked by hand on an A30 by restarting both
tests at steps 1,600,000 and 2,600,000 (see `docs/STATE.md`).
