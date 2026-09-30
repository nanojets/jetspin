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
| Injection | `inserting`, `removing`, `nozzle velocity`, `nozzle stress`, `dragvel yes` |
| Perturbation | `perturbation yes`, `perturbation frequency`, `perturbation amplitude` |
| Electric field | `external potential`, `external potential type`, `external potential vector` |
| Aerodynamics | `airdrag yes`, `airdrag airdensity`, `airdrag airviscosity`, `airdrag airvelocity` |
| Stochastic noise | `noise yes`, `noise variance`, `noise diffusivity`, `noise pool` |
| Rheology | `hbfluid yes`, `kvfluid yes`, `yield stress` |
| Refinement | `dynamic refinement yes`, `dynamic refinement every`, `dynamic refinement threshold` |
| Evaporation | `evaporation yes`, `evaporation polymer frac`, `evaporation temperature`, `evaporation umidity` |
| Output | `print time`, `print list`, `printstat list`, `print xyz`, `print pdb frame` |
| Restart | `restart yes`, `restart dump`, `restart reset` |

Some historical keywords use spellings retained for compatibility, such as
`evaporation umidity`. Copy them exactly as documented.

`noise pool n` sets the number of pre-generated Gaussian values read by the
stochastic Platen integrator when it uses the pre-generated pool (fixed
1,000-bead benchmarks and the dynamic evaporative Platen path), between
1,000,000 and 2,000,000,000; the default is 100,000,000. The pool is read
cyclically, so the noise repeats after `n / (6 * active beads)` timesteps:
about 62,000 steps for 270 beads with the default. Raising `n` lengthens the
period at 8 bytes per value on the host and again on the GPU. Whether the
period is long enough is the user's choice; see
[random numbers](../introduction/random-numbers.md#the-noise-repeats-after-the-pool-is-consumed).

## How records are read

The parser is `read_input` in [`io_mod.f90`](../../source/io_mod.f90), with
the string helpers in [`parse_mod.f90`](../../source/parse_mod.f90). Knowing
its rules avoids silent misreadings:

- Each line is read up to 150 characters; anything beyond is lost. Leading
  blanks are removed and the whole line is converted to lower case.
- Lines starting with `#` or `!`, and blank lines, are skipped. There are no
  end-of-line comments: text after the values is still part of the record.
- A keyword is recognized when it matches the **beginning of any word** in
  the line, not necessarily the first one. Top-level directives are tested
  in a fixed order and the first match wins; sub-keywords (`yes`, `every`,
  `freq`, ...) are then tested the same way.
- A numeric value is the **first number** found scanning the line from the
  left (`1.d-3`, `5e-9`, `-2`, `.5` are all accepted). Directives with
  several values, such as `external potential vector`, take the following
  numbers in order.
- An unknown top-level directive stops the run with
  `ERROR - unknown directive in input file`. An unknown sub-keyword prints
  warning 61 and the run stops after reading (`ERROR - incomplete input
  file`). Reaching the end of the file without `finish` also stops the run;
  lines after `finish` are ignored.

Consequences worth remembering:

- Abbreviations are equivalent: `dynam refin every` and `perturb freq`, used
  by the manual and several examples, are the same as
  `dynamic refinement every` and `perturbation frequency`.
- An invented keyword can silently match another directive:
  `timesteps 15000` is read as `timestep 15000`, a 15000 s integration step.
  The run length is set only by `final time`.
- Keep comments on their own `#` line. A trailing note such as
  `viscosity 20.d0 initial value` is taken by the `initial` directive, which
  is tested before `viscosity`, and fails there.

`dragvel yes` gives each newly inserted nozzle bead the dragging velocity
described in the manual's jet-insertion section. It is off by default: an
inserted bead then starts with the nozzle velocity along the jet axis.

The authoritative directive table is in the [PDF manual](../../manual/manual.pdf)
and its [LaTeX source](../../manual/input.tex). Evaporation-specific input is
described in [`evaporation-input.tex`](../../manual/evaporation-input.tex).
