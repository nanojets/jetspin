# Output files

JETSPIN writes output in the directory from which `main.x` is executed.
Files are created only when their corresponding input directives are active.

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
