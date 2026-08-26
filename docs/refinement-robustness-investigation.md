# Refinement robustness under repeated remeshing (2026-08-15/16)

This note records an investigation into the numerical robustness of the
dynamic-refinement path when a jet grows from a single nozzle bead through
many tens of consecutive accepted refinement events. It is kept separate from
any example case: the configuration used as a probe has since changed, while
the three code fixes and the cadence guidance the investigation produced are
in the shipped code and still apply.

Every previously validated refinement test starts from a pre-extended jet, so
this regime -- a young jet remeshed repeatedly as it grows -- had never been
exercised before. Three distinct problems surfaced; two were fixed, and the
third was contained through a cadence choice rather than resolved at the
mechanism level.

See also [dynamic refinement](introduction/dynamic-refinement.md).

- The probe configuration was a long production run intended to reach a
  statistically stationary active-bead count: Maxwell rheology with Yarin
  evaporation and a tight refinement cadence, combined with a
  single-nozzle-bead startup and a long integration window. It is the first
  case that grows from one bead through many tens of accepted dynamic-
  refinement events rather than a short validation window, and it still does
  not run to completion. This entry documents everything found and fixed in
  this session; the case remains untracked and unresolved as a checkpoint.
- GFortran CPU is bit-identical between `master` and `development` for this
  input: both crash with `ERROR - numerical instability` at the same step,
  confirming the failure is not a `development`-branch GPU-porting
  regression but a pre-existing property of the shared CPU refinement model.
- **Root cause 1 (evaporation path, fixed and committed).** The accepted-
  event reconstruction of `jetve` (post-evaporation volume) in
  `reconstruct_refinement_state_host` (`dynamic_refinement_mod.f90`) did not
  enforce the `evlim` floor that every ordinary per-timestep integrator path
  already enforces (`jetve(i)/jetvl(i)>=evlim`). An Akima undershoot on the
  evaporation-radius field let this ratio collapse to about 0.00916 (against
  a floor of `evlim=cp0/(1-evsolvlim)~=0.0667` for this input), i.e. roughly
  7x below its own per-bead floor, uniformly across nearly the whole active
  jet at once. This inflates the polymer mass fraction `cp=cp0*jetvl/jetve`
  to about 6.55 (655%, matching the number already recorded in
  issue 2 below) and, in `eom_ev_mod.f90`'s
  `eom4_ev` (and the analogous `eom1_ev`/`eom3_ev`/`eom4_pos_ev` branches),
  the evaporation-corrected viscosity ratio `ratmu` explodes accordingly,
  driving the Maxwell stress derivative `fst=(1/rattao)*consistency*ratmu*(...)`
  into the thousands on that same event and then, integrated over roughly
  1700 further ordinary timesteps, into a genuine double-precision overflow.
  Fixed with `enforce_evlim_conservative` (new subroutine,
  `dynamic_refinement_mod.f90`), called right after the existing
  evaporated-volume conservation rescale: a water-filling algorithm clamps
  any bead below `evlim*jetvl(i)` up to that floor, then takes the resulting
  deficit back from the still-compliant beads in proportion to their own
  surplus above their floor, iterating in case that step creates a new
  violator. This restores the floor while preserving the segment's
  evaporated-volume total to roundoff -- the same invariant the existing
  `evaporation_volume_relative_difference` check already validates -- rather
  than the previously-reverted naive clamp, which broke that conservation.
  If the floor is intrinsically infeasible against the segment total (no
  redistribution can satisfy it), the code now aborts with the new
  `error(20)` rather than silently violating either constraint. Confirmed
  by direct instrumentation (`JETSPIN_DEV`-style print, guarded on
  `rattao<1e-2 .or. rattao>1e2`, added at every `fst=(1/rattao)*(...)`
  assignment site in `eom_ev_mod.f90`): with the fix, this diagnostic never
  fires again for this input. Validated with zero regressions on Tests
  21/22/23 (GFortran).
- **Root cause 2 (evaporation-independent, mitigated but not eliminated).**
  Even with evaporation entirely disabled (`evaporation no`), the same class
  of crash reproduces: `jetcr` (cross-section radius), fit independently by
  Akima and then squared into `jetvl=length*pi*jetcr**2` (and, via
  `convert_from_density`, into `jetms`), can collapse across many
  *successive* accepted refinement events to a physically nonsensical scale,
  eventually driving `jetms` toward zero. In the non-evaporative `eom4`
  branch (`eom_mod.f90`) this makes `Fvet=Fve/jetms(ipoint)`,
  `attt=att/jetms(ipoint)` explode, and the resulting huge acceleration
  integrates into huge velocity within a handful of steps, which the stress
  equation's `consistency*(beadvelup/beadlenup)` term then overflows into a
  `st_nan=T` detection matching the evaporative case's signature. Confirmed
  this is *not* physically expected breakup: the user's domain expectation
  is that this jet reaches the collector intact (electrospinning genuinely
  thins a fibre by orders of magnitude -- up to about 1/1000 of the nozzle
  radius over the full nozzle-to-collector run is plausible -- but this
  specific collapse is faster/deeper than that and does not occur at all
  when Akima refinement is not exercised). `breaking_mod.f90` already
  implements a physically-motivated breakup detector
  (`condition_breakup_1`: root-find a local cubic through `jetcr` for an
  imminent radius-zero singularity) and a mesh-compaction handler
  (`clean_breakup`), but the `breakup yes` input keyword is gated behind the
  compile-time `parameter :: ldevelopers=.false.` in `nanojet_mod.f90` (which
  also gates several other experimental features), so it was never available
  to this input; it was not enabled during this investigation because
  flipping `ldevelopers` unlocks that broader, unvalidated developer-mode
  surface, not just breakup, and the user's diagnosis is that this is an
  interpolation-oscillation defect rather than a case for that mechanism.
- Diagnosis method for root cause 2: `Dynamic refinement thin-bead
  classification` (added to `dynamic_refinement_mod.f90`'s
  `reconstruct_refinement_state_host`) shows the per-event minimum-radius
  bead is *never* an anchor (`jetbd`), and for most of the run sits within a
  few beads of the fitted segment's collector-side endpoint (`jptinit`,
  which by construction is never tagged as an anchor -- see
  `initialize_tagged_beads` in `nanojet_mod.f90`). In the last 1-3 events
  before each crash its minimum-radius location instead jumps tens of beads
  inward and its value drops roughly an order of magnitude *between two
  consecutive events*, i.e. each event re-fits from the previous event's
  already-slightly-degraded state and can make it measurably worse -- a
  cross-event compounding drift, not a single-event artifact. Explicitly
  ruled out as contributors, with direct evidence: `mass_density`/
  `charge_density` are exactly constant throughout this input (uniform
  `density mass`/`density charge`, nothing in this model varies them
  spatially) and never show the pattern; position/segment length
  (`min_segment_length_cm`) fluctuates non-monotonically and shrinks only
  about 13x over the same 21 events where radius shrinks about 160x and
  never once increases, so the defect is isolated to the cross-section
  radius/area reconstruction, not a general Akima weakness across all 11
  interpolated fields; stress and velocity do eventually show huge
  excursions, but only in the single event immediately before each crash,
  after radius has already been collapsing smoothly and monotonically for
  many prior events -- consistent with them being a late-stage Newton's-law
  consequence (tiny mass -> huge acceleration -> huge velocity -> stress
  overflow) rather than an independent source.
- Five layered, incremental defenses were implemented in
  `dynamic_refinement_mod.f90`/`fit_mod.f90`, each validated with zero
  regressions on `tests/refinement/run.sh gfortran standard` (Test 21) and
  direct GFortran runs of `examples/input-22`/`examples/input-23` (clean
  completion, matching topology). None of them, individually or combined,
  fully eliminates the crash; each measurably changes (not monotonically
  improves, given the stochastic Platen trajectory's sensitivity to any code
  change) how far the run gets before it. In order of application inside
  the `radius_area`/`evap_radius_area` fit blocks:
  1. **Fit `ln(pi*jetcr**2)` (log cross-section area), not the raw radius.**
     Reference volume is linear in area but quadratic in radius, so a radius
     undershoot is squared into the reconstructed volume/mass; area removes
     that amplification. Fitting its *logarithm* additionally guarantees the
     recovered area is strictly positive by construction (`dexp` of anything
     is positive -- a raw area fit could and did go measurably negative
     before `dabs()`, confirmed directly: `raw_fit_min=-1.014459E-06`) and
     linearises the roughly exponential thinning trend, which is much better
     conditioned for Akima's tangent estimate than the raw area's steeper
     curvature. `jetcr` is recovered as `sqrt(exp(fit)/pi)`. The source
     array passed to `fit_akima` is guarded with
     `max(dabs(jetcr(:)),1.d-300)` before `dlog` to avoid a spurious
     `log(0)` (`IEEE_DIVIDE_BY_ZERO`) on unused array slots outside
     `inpjet:npjet`; this is cosmetic (the run already completed correctly
     without it) but avoids an unnecessary FP flag.
  2. **`despike_median_filter`** (new subroutine): a branchless 3-point
     median filter over the fitted target segment, iterated up to 3 passes,
     skipping anchor points entirely (an earlier version that did not skip
     anchors broke the `anchor_displacement`/radius-exactness invariant and
     was caught by Test 21 -- always exclude `jetbd` from any post-hoc
     filter on these fields).
  3. **`clamp_to_local_source_range`** (new subroutine): bounds each
     non-anchor target point to the range spanned by its two immediate
     bracketing pre-event source points (a cheap, local approximation of a
     shape-preserving/monotone-interpolation guarantee that classic Akima
     does not provide). Confirmed by direct measurement that this clamp
     cannot by itself catch a drift that already compounded through earlier
     events, because by the time of a later event that drift is already
     baked into the "legitimate" local source data too.
  4. **`clamp_to_anchor_envelope`** (new subroutine, safety factor `1.d2`):
     bounds interior non-anchor points between two bracketing *anchors*
     instead of immediate source neighbours, since anchors are restored
     exactly every event (`radius_max_difference_cm=0`) and so should be
     drift-immune. Found empirically insufficient in this run because an
     anchor itself can already be tagged from an already-partially-degraded
     bead (anchors are frozen at tagging time, not validated against
     anything), so trusting an anchor as ground truth does not always hold.
  5. **`limit_akima_tangents_monotone`** (new subroutine in `fit_mod.f90`,
     operating on the module-level `mak`/`tak` arrays): a Fritsch-Carlson
     sufficient-condition limiter on the Akima knot tangents -- zero the
     tangent at any knot where the two adjacent secant slopes disagree in
     sign (a local extremum), then rescale each segment's two endpoint
     tangents, iterated a few passes, to satisfy the classic
     circle-of-radius-3 monotonicity bound. Opted into via a new optional
     `lmonotone` argument on `setup_akima`, enabled only for
     `field_name` in `{radius_area, evap_radius_area, mass_density,
     charge_density}` (never position/velocity/stress, which may have
     genuine local extrema). This measurably delayed the crash the most of
     any single change so far, but the compounding drift still eventually
     recurs -- even a genuinely monotone-consistent fit, re-run from its own
     prior output across ~20-30 events, appears to accumulate some drift,
     which may be a resampling/requantization effect distinct from classic
     Akima overshoot and not yet isolated.
  6. **`enforce_radius_floor_conservative`** (new subroutine, water-filling
     identical in structure to the evlim fix, floor
     `minimum_bead_radius_cm=1.d-7` cm, i.e. 1 nm): found empirically too
     permissive to matter here -- the measured collapsing radii (order
     5-40 nm) never actually cross an absolute 1 nm floor, so this call
     never fires for this input. Electrospinning legitimately reaches
     sub-micron fibre radii, so an *absolute* floor is the wrong kind of
     check; a floor relative to the jet's own current/local scale was
     discussed but not implemented, superseded by the anchor-envelope and
     monotonicity approaches above. Kept as a defense-in-depth against a
     genuine sign-crossing that the other layers miss, and as the
     `error(21)` infeasibility guard.
- New error codes added to `error_mod.f90`: `case(20)` (evlim floor
  infeasible against the segment's evaporated-volume total) and `case(21)`
  (radius floor infeasible against the segment's reference-volume total).
  Both abort with a clear message rather than silently breaking either the
  floor or the conservation invariant.
- **GPU/OpenACC scope of all fixes above: originally host-only; now ported,
  see the dedicated section below.** This investigation itself, per the
  user's direction, targeted GFortran CPU only, and at the time of writing
  three of the six defenses (`enforce_evlim_conservative`,
  `enforce_radius_floor_conservative`, `limit_akima_tangents_monotone`)
  were reachable only from the host paths. They have since been ported to
  their OpenACC counterparts; see "GPU/OpenACC port of the three host-only
  fixes" below for what changed and how it was validated.
- Guarded, low-noise diagnostic instrumentation was added and left in the
  tracked source (development-only, no effect on normal output unless a
  threshold is crossed): `Rattao blowup diagnostic` in `eom_ev_mod.f90`
  (fires when `rattao<1e-2 .or. rattao>1e2`, reporting `yvl`, `yve`,
  `cmass`, `rattao`, `fst` for the offending bead); `Akima endpoint
  diagnostic` in `fit_mod.f90` (fires for `field_name` in
  `{radius_area, evap_radius_area, mass_density, charge_density, stress,
  vx, vy, vz}`, reporting the endpoint segment slopes/tangent and the raw
  fitted value/min/max before any of the filters above run); `Dynamic
  refinement thin-bead classification` in `dynamic_refinement_mod.f90`
  (reports whether the event's minimum-radius bead is an anchor and its
  distance from the segment endpoint); and, in `integrator_mod.f90`'s
  instability-detection block, `Numerical instability last-good neighbor`/
  `domain extent` (reports the physical `x` position in cm of the last
  still-finite neighbour and of both the collector-side (`inpjet`) and
  nozzle-side (`npjet`) active-mesh extremes). These were essential for
  distinguishing "near the nozzle/head" from "near the collector/tail" and
  for confirming which field's interpolation was implicated; keep them for
  any future continuation of this investigation.
- Confirmed collector-reaching state at the point of failure differs by
  configuration: with evaporation enabled, the jet does reach the collector and 34 beads have already been
  removed there when the crash occurs (`collector_side_bead` at
  `x_cm=12.000000` exactly); with evaporation disabled, it does not (stops
  at about `x=8.97` of 12 cm, zero removals).
- Status: the probe configuration still does not run to completion under any of
  the above combinations. This is accepted as an incremental checkpoint, not
  a resolution. Suggested next directions, not started: investigate whether
  the residual drift is a cumulative resampling/requantization effect
  intrinsic to re-fitting the same evolving quantity across dozens of
  events (rather than a per-event Akima quality issue, which the above five
  layers already target); consider a true monotone/shape-preserving cubic
  Hermite (PCHIP) replacement rather than a post-hoc limiter on top of
  Akima; consider whether a region should stop being re-fit after some
  number of consecutive events without becoming stable.

### Resolution: refinement cadence, not interpolation quality, was the practical driver

- Comparing against `examples/input-5` (the historical single-nozzle-bead
  refinement case) showed it uses a much coarser cadence than the probe:
  `dynamic refinement threshold 0.4` cm (20x the 0.02 cm base resolution)
  and `dynamic refinement every 1.d-3` s, versus the probe's `0.10` cm
  (5x resolution) and `1.d-5` s. The tight cadence lets a jet
  growing from one bead accumulate several dozen accepted refinement events
  within under `2e5` timesteps; each event re-fits the cross-section from
  the previous event's own output, which is exactly the precondition for
  the cross-event compounding drift documented above.
- Re-ran the probe configuration (with evaporation) with only
  `threshold` changed to `0.4d0` and `every` changed to `1.d-3`, all five
  cross-section-fit defenses and the `evlim` fix still active. **This run
  reached full completion**: all `1e8` (100 million) target timesteps,
  `t=0.5` s, 499 accepted refinement events, `Program closed correctly`,
  zero `ERROR` lines. (An earlier check of this same run at `6.5e7` steps,
  before completion, was reported elsewhere as "interrupted externally"; the run in
  fact continued unattended and completed on its own while later GPU-port
  work was underway. Those earlier notes have been superseded and should
  be read as describing an intermediate checkpoint, not the final state.)
  Aggregate result: `Topology additions=17871`, `Topology removals=27320`,
  `Array reallocations=0`, wall time `5055 s` (`~19782 steps/s`).
  The jet first bridges the full 12 cm nozzle-to-collector distance and
  reaches its first collector removal at step `95528` -- under 0.1% into
  the run. From that point on the active-bead count is already in its
  final stationary regime: the very first accepted refinement event
  (`step=249999`) and the very last (`step=99849999`, essentially the
  entire `1e8`-step span apart) report the identical
  `active_before=17 active_after=36` pattern, with every event in between
  reporting either `17/36` or `18/37`. Between two consecutive refinement
  events the active count settles back down to `17`-`18` from the
  post-event `36`-`37`, in a repeating sawtooth that never drifts across
  the full run. `min_radius_cm` measured at every one of the 499 events
  stayed at its initial value (`~5.40e-4` cm) for the entire run -- no
  collapse, gradual or sudden. This is the first time this input has run
  to completion and demonstrated, over its full intended duration rather
  than an extrapolated partial run, the stationary active-bead behavior it
  was designed to reach.
- The guidance that came out of this is `dynamic refinement threshold 0.4d0`
  with `dynamic refinement every 1.d-3` for long single-bead-start runs,
  leaving the anchor spacing unchanged.
- Added a general advisory, `warning(108)` in `error_mod.f90`, raised from
  `set_refinement_threshold()` (`dynamic_refinement_mod.f90`) whenever the
  final (configured or auto-derived) `refinementthreshold` is below 20
  times the base `resolution`. It is advisory only and does not alter
  behavior: Example 5 itself (5x resolution) and Tests 21-23 all use
  thresholds below this recommendation and remain validated short-window
  references where the compounding drift was never exercised long enough
  to appear. Documented in `manual/input.tex`'s directive table and
  `manual/refinement.tex` Sec. "Numerical robustness of the cross-section
  reconstruction", and in `docs/introduction/dynamic-refinement.md`.
- This is a practical workaround for the probe's specific parameters, not
  a fix for the underlying cross-event compounding mechanism itself, which
  remains as documented above (still suspected to be a cumulative
  resampling/requantization effect, not yet isolated further). An input
  that genuinely needs both a tight refinement threshold and a long,
  many-event integration window would still be expected to hit this
  failure mode; the five interpolation-side defenses reduce its severity
  but do not eliminate it, per the extensive earlier testing in this
  section.
