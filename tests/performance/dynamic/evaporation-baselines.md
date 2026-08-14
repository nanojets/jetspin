# Dynamic evaporation CPU baselines

These baselines use the 100-bead Test Cases 16 and 17 with a serial CPU build.
The stored inputs use RK4, Yarin evaporation, insertion, removal, and a
deliberately small capacity that forces reallocation. Variants changing only
`integrator` to `1` or `2` validate Euler and RK2.

| Case | Rheology | Additions | Removals | Reallocations | Final active beads |
|---|---|---:|---:|---:|---:|
| 16 | Maxwell | 111 | 122 | 2 | 89 |
| 17 | Kelvin–Voigt | 111 | 122 | 2 | 89 |

Both runs completed with `Program closed correctly` and no NaN or infinity
diagnostics. The complete Maxwell and Kelvin--Voigt GPU paths reproduce these
aggregate totals for Euler, RK2, and RK4. Three-step pre-event comparisons for
all three integrators are exact at `rtol=1e-12`, `atol=1e-13`; their written
Maxwell XYZ geometries are byte-identical between the same NVFORTRAN CPU and
A30 builds. Preserve the complete `statout.dat` and topology-event stream with
the exact compiler and flags used for any new reference build. Set
`JETSPIN_TOPOLOGY_SNAPSHOT=1` to write
`topology-state.dat`; its event records include `jetve` and `jetce` as the final
two columns. Compare same-path snapshots with:

```sh
python3 tests/performance/dynamic/compare_state.py \
  cpu/topology-state.dat gpu/topology-state.dat \
  --rtol 1e-6 --atol 1e-12
```

Because direct Coulomb accumulation order differs on the CPU and GPU and the
bending instability amplifies tiny perturbations, a full pointwise CPU/GPU
trajectory is not the Test 16 or Test 17 acceptance criterion. Use completion,
topology totals, transfer audits, and comparison with a saved same-build GPU
reference instead.
