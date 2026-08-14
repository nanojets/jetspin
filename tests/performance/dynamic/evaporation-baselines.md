# Dynamic evaporation CPU baselines

These baselines use the 100-bead Test Cases 16 and 17 with the serial
`gfortran` target. Both runs use RK4, Yarin evaporation, insertion, removal,
and a deliberately small capacity that forces reallocation.

| Case | Rheology | Additions | Removals | Reallocations | Final active beads |
|---|---|---:|---:|---:|---:|
| 16 | Maxwell | 111 | 122 | 2 | 89 |
| 17 | Kelvin–Voigt | 111 | 122 | 2 | 89 |

Both runs completed with `Program closed correctly` and no NaN or infinity
diagnostics. The complete Maxwell RK4 GPU path in Test 16 reproduces these
aggregate totals and its final sampled output is bit-for-bit identical to the
preceding saved OpenACC result. Preserve the complete `statout.dat` and
topology-event stream with the exact compiler and flags used for any new
reference build. Set `JETSPIN_TOPOLOGY_SNAPSHOT=1` to write
`topology-state.dat`; its event records include `jetve` and `jetce` as the final
two columns. Compare same-path snapshots with:

```sh
python3 tests/performance/dynamic/compare_state.py \
  cpu/topology-state.dat gpu/topology-state.dat \
  --rtol 1e-6 --atol 1e-12
```

Because direct Coulomb accumulation order differs on the CPU and GPU and the
bending instability amplifies tiny perturbations, a full pointwise CPU/GPU
trajectory is not the Test 16 acceptance criterion. Use completion, topology
totals, and comparison with the saved same-build GPU reference instead.
