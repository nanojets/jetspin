# Test Case 21: anchored dynamic refinement with Platen evaporation

Test Case 21 validates the first hybrid OpenACC port of dynamic Akima
refinement. It combines Maxwell rheology, Yarin evaporation, stochastic
Platen integration, aerodynamic drag, nozzle insertion, and numerical anchor
beads. Integration and the repeated threshold scan remain on the GPU. When a
refinement event is accepted, target-mesh construction remains on the host,
while Akima coefficient construction and spline interpolation run on the GPU.

The initial jet contains 400 elements at 0.02 cm resolution and occupies only
the first 8 cm of the 16 cm nozzle-to-collector distance. The remaining span
allows the electric field to stretch the jet. Interior anchors are initialized
every 0.10 cm and retained as interpolation knots; they carry the same material
model as every other bead and do not represent nanoparticles.

The high potential and reduced charge density are deliberate test
accelerators. This input is a numerical validation workload, not a recommended
physical parameter set. Its accepted GFortran run performs one refinement,
preserves all pre-existing anchor coordinates, conserves reference and
evaporated volumes to roundoff, and ends with `nref = 1`.

Run and validate the complete case with:

```sh
tests/refinement/run.sh
tests/refinement/run.sh openacc
```

The normal input reserves 100 capacity entries. The test runner's optional
`capacity-growth` mode uses a developer-only 50-entry reserve so this same
physical input exercises release, host reallocation, and persistent-device
rebind during its accepted refinement event.
