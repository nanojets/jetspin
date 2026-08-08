# Yarin-2001 evaporation regression checks

`check_yarin2001.py` protects the parameter mapping and the main evaporation/rheology formulas used by Test Case 8.

It checks the Test 8 input against fixed numerical targets:

- `theta0 = mu0/G0 = 0.01 s`;
- atmospheric-pressure fallback `D_a = 0.242001758074 cm^2/s` at `293.15 K`;
- Yarin solidification cutoff `V/V0 = cp0/0.9 = 0.066666666667`;
- `cp = 0.9` at the cutoff;
- `mu/mu0 = 43.9783350021`;
- `theta/theta0 = 15`;
- `G/G0 = 2.93188900014`.

When the Fortran sources are supplied, the checker also verifies that the production code still contains the Yarin concentration, viscosity, relaxation-time and cutoff expressions, together with the corrected Seaver atmospheric-pressure diffusivity fallback.  When a JETSPIN `run.log` is supplied it compares the diffusivity actually initialized by the executable with the analytical value above.

Run the formula/source checks directly from the repository root with:

```sh
python3 tests/evaporation/check_yarin2001.py \
    --input examples/input-8/input.dat \
    --nanojet-source source/nanojet_mod.f90 \
    --eom-source source/eom_ev_mod.f90
```

The normal smoke-test driver runs Test Case 8 in both serial and debug builds and adds the run-time diffusivity check automatically:

```sh
tests/smoke/run.sh serial
tests/smoke/run.sh debug
```

The MPI smoke suite remains intentionally short and currently runs only Test Case 1.  The GitHub Actions `Smoke tests` workflow is triggered by changes under `tests/evaporation/` as well as by source and example changes.

The checker is deliberately a regression guard, not an independent reimplementation of the full Yarin trajectory.  Test Case 8 itself exercises the coupled JETSPIN equations, while the fixed values above provide stable diagnostics for the evaporation cutoff and concentration-dependent rheology.
