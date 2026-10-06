# Exact restart check

`run.sh` checks that a restarted run continues the uninterrupted one exactly,
for the dynamic stochastic Platen runs of Test 24 (no evaporation) and
Test 25 (evaporation). For each test it builds the requested executable and
runs three copies of the test's input, with a 20,000-step print interval and
`restart dump 100000`:

- run A stops at step 100,000 and writes `save.dat`;
- run B reads it as `restart.dat` (`restart yes`) and stops at step 200,000;
- run C runs the same 200,000 steps without interruption.

`check_restart.py` then requires every `statout.dat` row, every terminal row
and every `traj.xyz` frame written by B after step 100,000 to equal C's,
character for character. This exercises the double-precision bead records,
the Gaussian-pool cursor, and the refinement and anchor counters stored in
`save.dat` since 2026-10-06; before that date Test 25 could not be restarted
("restart file is corrupted") and Test 24 restarted with a different noise
sequence.

```sh
tests/restart/run.sh nvfortran   # NVFORTRAN CPU build, about one minute
tests/restart/run.sh openacc     # OpenACC build
```

The OpenACC check stays on the host path (fewer than 100 beads). The
persistent device path was checked by hand on an A30 by restarting both
tests at steps 1,600,000 and 2,600,000 (see `docs/STATE.md`).
`JETSPIN_RESTART_KEEP=1` keeps the work directory and
`JETSPIN_RESTART_TIMEOUT` (default 600 s) bounds each run.
