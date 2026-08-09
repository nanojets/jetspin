# Test Case 13: dynamic-topology CPU/GPU benchmark

Test 13 starts with 1,024 beads uniformly distributed from nozzle to collector
and enables both `inserting yes` and `removing yes`. It retains the constant
axial field and material model of Test 9 but uses an intentionally artificial
2,000 cm/s nozzle velocity so that a short 1,000-step run exercises topology.

The reference run produces 13 insertions and 13 removals. CPU and OpenACC runs
perform every event at the same timestep, finish with 1,024 active beads, and
remain between 1,023 and 1,025 active beads. The program prints a structured
`Topology event:` record whenever either boundary changes the topology.

This first milestone is a baseline rather than a persistent dynamic GPU port:
insertion, removal, capacity growth, and compaction still execute on the host,
while accelerator force calls use call-scoped mappings. It establishes the
event sequence that future device-resident topology work must preserve.

An initial NVFORTRAN 24.3/A30 comparison measured `40.659696 s` on CPU and
`3.834931 s` with OpenACC. The topology streams are exact. The sampled physical
observables pass with `rtol=3e-4` and `atol=1e-6`; the looser dedicated
tolerance accounts for curvature amplification after repeated topology events.

- [Input file](../../examples/input-13/input.dat)
- [Comparison record](../../tests/performance/dynamic/README.md)
