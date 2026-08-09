# Input files

Every simulation requires a file named `input.dat` in the directory where
`main.x` is launched. Input is free-format and case-insensitive. Blank lines
and records beginning with `#` are ignored, and the final record must be
`finish`.

All dimensional quantities use the centimetre–gram–second (CGS) system.

## Minimal structure

```text
system 3
integrator 3
timestep 1.d-6
final time 1.d-3

initial length 0.02d0
nozzle cross 0.005d0
density mass 0.84d0
density charge 44000.d0
viscosity 20.d0
elastic modulus 50000.d0
collector distance 16.d0

finish
```

Use an existing [`examples/input-*`](../../examples/) case as the starting
point for a new simulation.

## Main directive groups

| Group | Representative directives |
| --- | --- |
| Model | `system`, `integrator`, `timestep`, `final time` |
| Initial jet | `initial length`, `points`, `resolution`, `nozzle cross` |
| Material | `density mass`, `density charge`, `viscosity`, `elastic modulus`, `surface tension` |
| Injection | `inserting`, `removing`, `nozzle velocity`, `nozzle stress` |
| Perturbation | `perturbation yes`, `perturbation frequency`, `perturbation amplitude` |
| Electric field | `external potential`, `external potential type`, `external potential vector` |
| Aerodynamics | `airdrag yes`, `airdrag airdensity`, `airdrag airviscosity`, `airdrag airvelocity` |
| Rheology | `hbfluid yes`, `kvfluid yes`, `yield stress` |
| Refinement | `dynamic refinement yes`, `dynamic refinement every`, `dynamic refinement threshold` |
| Evaporation | `evaporation yes`, `evaporation polymer frac`, `evaporation temperature`, `evaporation umidity` |
| Output | `print time`, `print list`, `printstat list`, `print xyz`, `print pdb frame` |
| Restart | `restart yes`, `restart dump`, `restart reset` |

Some historical keywords use spellings retained for compatibility, such as
`evaporation umidity`. Copy them exactly as documented.

The authoritative directive table is in the [PDF manual](../../manual/manual.pdf)
and its [LaTeX source](../../manual/input.tex). Evaporation-specific input is
described in [`evaporation-input.tex`](../../manual/evaporation-input.tex).
