# Input files

Every simulation requires a file named `input.dat` in the directory where
`main.x` is launched. Input is free-format and case-insensitive. Blank lines
and records beginning with `#` are ignored, and the final record must be
`finish`.

All dimensional quantities use the centimetre–gram–second (CGS) system
(`units real` changes only the charge density and the potentials, below).

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
| Model | `system`, `integrator`, `timestep`, `final time`, `seed`, `units` |
| Initial jet | `initial length`, `points`, `resolution`, `nozzle cross`, `collector distance` |
| Material | `density mass`, `density charge`, `viscosity`, `elastic modulus`, `surface tension` |
| Injection | `inserting`, `removing`, `nozzle velocity`, `nozzle stress`, `dragvel yes` |
| Perturbation | `perturbation yes`, `perturbation frequency`, `perturbation amplitude` |
| Electric field, gravity | `external potential`, `external potential type`, `external potential vector`, `external potential freq`, `external potential phase`, `external potential time`, `gravity yes` |
| Aerodynamics | `airdrag yes`, `airdrag airdensity`, `airdrag airviscosity`, `airdrag airvelocity`, `airdrag amplitude` |
| Stochastic noise | `noise yes`, `noise variance`, `noise diffusivity`, `noise pool` |
| Rheology | `hbfluid yes`, `hbfluid consistency`, `hbfluid index`, `kvfluid yes`, `yield stress` |
| Coulomb multiple step | `multiple step yes`, `multiple step every`, `primary cutoff`, `maximum displ` |
| Refinement | `dynamic refinement yes`, `every`, `threshold`, `start`, `anchor`, `capacity` (see the [dynamic-refinement guide](../introduction/dynamic-refinement.md#configuration)) |
| Evaporation | `evaporation yes`, `evaporation polymer frac`, `evaporation temperature`, `evaporation umidity`, ... (see [`evaporation-input.tex`](../../manual/evaporation-input.tex)) |
| Output | `print time`, `print list`, `printstat list`, `print xyz`, `print xyz frame`, `print xyz maxnum`, `print xyz rescalexyz`, `print pdb frame`, `print pdb rescalepdb`, `print binary`, `printstat binary` |
| Restart, job control | `restart yes`, `restart dump`, `restart reset`, `job time`, `close time` |

Defaults and conditions that the names do not show:

| Directive | Meaning, default and conditions |
| --- | --- |
| `seed <i>` | Seed of the random generator; default 1 |
| `units cgs`, `units real` | Default `cgs`. `real` reads `density charge` in C/L and the `external potential` values (scalar and vector) in V and converts them; every other quantity stays CGS |
| `resolution <f>`, `points <i>` | Exactly one is required (warning 7 if neither, 22 if both); the resulting step must not be smaller than `nozzle cross` (warning 62, error 8) |
| `elastic modulus`, `initial length`, `perturb ampl` | Checked at input: elastic modulus at least 1e-4 (warning 88), initial length not smaller than `nozzle cross` (warnings 42-44), perturbation amplitude not larger than the initial length (warnings 45-47) |
| `system <i>`, `integrator <i>` | Not range-checked at input: an unsupported value or combination (e.g. `integrator 4` with `kvfluid yes`) stops the run at the first step with `ERROR - ftype is wrong!` (system) or `ERROR - ktype error.` (integrator) |
| `collector distance <f>` | Required (warning 16): distance of the collector from the nozzle along x |
| `gravity yes` | Adds gravity along x; default no |
| `external potential freq`, `phase`, `time` | Waveform frequency in 1/s (required when `type` is not 0, warning 71), phase in rad of the rotating y and z components of `type 3` (ignored by types 1 and 2; default 0), RC time constant in s (required with `type 2`, warning 78) |
| `external potential vector <f> <f> <f>` | Potentials of the x, y, z field components (required with `type 3`, warning 90, and accepted only with it, warning 92): each component is the value divided by `collector distance`, y and z rotating at `freq`; in V with `units real`. With a vector the scalar `external potential` is not used for the field, while the output key `v` still prints the scalar times the waveform |
| `airdrag airdensity`, `airdrag airviscosity` | Required with `airdrag yes` (warnings 38, 39). The deterministic drag works with every integrator; `system 1` and `kvfluid yes` have no air-drag term and ignore it |
| `airdrag airvelocity <f>` | Velocity of the air along x (the jet axis), subtracted from the bead velocity in the drag term; default 0 (warning 49 with `airdrag yes`) |
| `airdrag amplitude <f>` | Amplitude of the random air-drag force: only with `system 4` and `integrator 4` (warnings 35, 34), and required there with `airdrag yes` (warning 33) |
| `hbfluid consistency`, `hbfluid index` | Required with `hbfluid yes` (warnings 40, 41); without it both are 1 |
| `multiple step every <n>` | Required with `multiple step yes` (warning 76), and n > 2 (warning 75). The neighbour list is rebuilt every n steps, and at the next Coulomb evaluation after an insertion, a removal, a remesh, a bead frozen at the collector (with or without removal) or an array growth (until 2026-10-07 a request arriving just after an update was dropped until the next scheduled one) |
| `primary cutoff <f>` | Distance in cm, with or without multiple step. With `multiple step yes` it is the neighbour-list radius and is required (warning 74); every direct sum used by the algorithm is then complete. Without multiple step it only truncates the 1-D (`system 1`) direct Coulomb sum (and the image terms of the developer-only `mirror`); the 3-D direct sum ignores it. Default: no cutoff. Until 2026-10-07 it was converted from cm only with multiple step, and it truncated the 1-D direct sum also under multiple step, which corrupted the extrapolation |
| `maximum displ <f>` | Only with `multiple step yes`: the neighbour list is also rebuilt when at least two beads (any beads, on any rank) have moved more than f/2 since the last rebuild. Default off. Keep `multiple step every` times `timestep` short compared with the motion over the cutoff (see the manual, multiple time step) |
| `print time <f>` | Terminal and `statout.dat` records every f seconds; default 1000 timesteps (warning 59) |
| `print xyz <f>`, `print xyz frame <f>`, `print xyz maxnum <i>`, `print xyz rescalexyz <f>` | `traj.xyz` every f seconds; one `frame%06d.xyz` every f seconds; beads written to `traj.xyz` only (default 100; the frame files hold all the beads); coordinate scale factor (default 1). `print pdb frame <f>` writes `frame%06d.pdb` and `frame%06d.psf` every f seconds, and `print pdb rescalepdb <f>` is the PDB scale factor (default 1). All these coordinates are in cm times the scale factor (until 2026-10-07 `traj.xyz` was in reduced units without `rescalexyz`, and the frame files applied the length unit twice with a scale factor) |
| `print binary <f> [style <i>]`, `print binary removed`, `printstat binary` | Developer binaries `traj.dat` (every f seconds, reduced units, layout `style` 1-9, default 1; warning 63 for a style outside 1-9, 77 for a style not allowed with the active options), `bead.dat` (beads removed at the collector, needs `removing yes`) and `statdat.dat` (all observables at every record); see [output files](output.md) |
| `restart dump <i>` | `save.dat` every i steps (default 100000); it is also written at the end of the run |
| `job time <t> [m, h or d]`, `job time indef`, `close time <t>` | Once less than `close time` seconds (default 0) are left of `job time` (seconds, or minutes, hours, days; default unlimited), the loop stops at the end of the step, the output files are closed and `save.dat` is written as at the final time; no extra statistics record is forced. Write the unit as a separate word (`job time 12 h`): a `d` attached to the number is read as an exponent mark, so `job time 2d` is 2 seconds. The clock is the elapsed (wall-clock) time since the program started, as the batch system counts it (until 2026-10-07 the processor time of rank 0, which does not count waits and could fall behind) |

The directives `lorentz`, `magnetic field`, `variable mass`, `wall`,
`mirror`, `breakup`, `ultimate strength` and `print pdb tagbead`, and the
`type` option of `dragvel yes`, are read only by developer builds
(`ldevelopers` in [`nanojet_mod.f90`](../../source/nanojet_mod.f90), `.false.`
in released builds). A released build stops at the first seven with `ERROR -
unknown directive in input file`, stops after reading at `print pdb tagbead`
(warning 61, `ERROR - incomplete input file`), and reads `dragvel yes type
<i>` as `dragvel yes`.

A bead that reaches the collector is frozen there, with or without `removing
yes`: it stops and leaves the Coulomb sums, as a discharged bead on the
grounded electrode. `removing yes` also deletes the collected beads from the
arrays (the leading one when the next reaches the collector); without it they
stay on the collector as deposited fiber, which the dynamic refinement
ignores. Until 2026-10-07 a run without
removal froze nothing, and its beads crossed the collector plane with their
charge.

Some historical keywords use spellings retained for compatibility, such as
`evaporation umidity`. Copy them exactly as documented.

`noise pool n` sets the number of pre-generated Gaussian values in the pool
read by every stochastic Platen run (`system 4`, `integrator 4`) whose
options the device step supports (`device_step_supported` in
[`device_step_mod.f90`](../../source/device_step_mod.f90): Tests 12, 20-25
and Example 4), in every build; the other runs draw their noise step by
step. `n` lies between 1,000,000 and 2,000,000,000, default 100,000,000; a
value outside prints warning 109 and the run stops after reading (`ERROR -
incomplete input file`). The pool is read
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
  `ERROR - unknown directive in input file`. For most directives an unknown
  sub-keyword prints warning 61 and the run stops after reading (`ERROR -
  incomplete input file`). `primary`, `maximum`, `multiple`, `restart`,
  `dragvel`, `evaporation`, and `dynam` without `refin` ignore it silently
  instead, and `external potential <word> f`, `print xyz <word> f` and
  `print binary <word> f` are read as `external potential f`, `print xyz f`
  and `print binary f`. An unknown key in `print list` or `printstat list`
  prints warning 48 or 58, then warning 60 with the valid keys, and the run
  stops after reading. Reaching the end of the file without `finish` also
  stops the run; lines after `finish` are ignored.

Consequences worth remembering:

- Abbreviations are equivalent: `dynam refin every` and `perturb freq`, used
  by the manual and several examples, are the same as
  `dynamic refinement every` and `perturbation frequency`.
- An invented keyword can silently match another directive:
  `timesteps 15000` is read as `timestep 15000`, a 15000 s integration step.
  The run length is set only by `final time`: the run makes `final time`
  divided by `timestep` steps, rounded up, and exactly that many when the
  final time is a multiple of the timestep. Until 2026-10-07 such a run made
  one step more with the compilers that round divisions exactly (GFortran,
  NVFORTRAN 25.5) than with NVFORTRAN 24.3, because the end test compared
  times already divided by the time unit.
- `print xyz rescale 1000` is read as `print xyz 1000`, a `traj.xyz` frame
  every 1000 s: the scale directive is `print xyz rescalexyz`.
- Keep comments on their own `#` line. A trailing note such as
  `viscosity 20.d0 initial value` is taken by the `initial` directive, which
  is tested before `viscosity`, and fails there.

`dragvel yes` gives each newly inserted nozzle bead the dragging velocity
described in the manual's jet-insertion section. It is off by default: an
inserted bead then starts with the nozzle velocity along the jet axis. Until
2026-10-07 the input summary printed `velocity drag yes by default` when the
directive was omitted, although the default was already off; it now prints
`velocity drag no by default`.

The authoritative directive table is in the [PDF manual](../../manual/manual.pdf)
and its [LaTeX source](../../manual/input.tex). Evaporation-specific input is
described in [`evaporation-input.tex`](../../manual/evaporation-input.tex).
