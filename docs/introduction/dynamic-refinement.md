# Dynamic refinement

JETSPIN's dynamic-refinement procedure prevents the Lagrangian jet
discretization from becoming excessively coarse as the filament stretches.
When an element exceeds a configured length threshold, part of the jet is
remeshed and its state is reconstructed on a denser parametrization.

The scientific formulation is described in the
[LaTeX manual](../../manual/refinement.tex). This page concentrates on the
runtime workflow and the invariants needed to maintain the implementation in
[`dynamic_refinement_mod.f90`](../../source/dynamic_refinement_mod.f90).

## Configuration

A minimal input fragment is:

```text
dynamic refinement yes
dynamic refinement every 1.d-3
dynamic refinement threshold 0.4d0
```

The supported controls are:

| Directive | Purpose |
| --- | --- |
| `dynamic refinement yes` | Enable refinement and bead tagging |
| `dynamic refinement no` | Disable refinement |
| `dynamic refinement every <time>` | Minimum interval between accepted refinements |
| `dynamic refinement threshold <length>` | Maximum target element length |
| `dynamic refinement start <time>` | Do not refine before this time |
| `dynamic refinement anchor <length>` | Spacing used to preserve anchor beads as interpolation knots |
| `dynamic refinement capacity <i>` | Bead-count capacity reserve/growth increment `incnpjet` (positive integer, default 100) |

Times and lengths use the normal input units: seconds and centimetres before
internal nondimensionalization. Both `every` and `threshold` are required
when refinement is enabled. The threshold must be at least twice the base
resolution; values below five times the resolution are raised to that safer
minimum. If internal defaults are used, the code derives them from the base
resolution and emits a warning.

`anchor` is optional. If omitted, it defaults to five times the base
resolution (warning 81). If given explicitly but set below one resolution,
it is likewise raised to five times the resolution (warning 85) -- the same
floor value used for the default, just reached from an invalid explicit
input instead. There is no upper bound: an anchor spacing coarser than the
refinement threshold is accepted without warning.

A threshold below twenty times the base resolution triggers a further
advisory (warning 108): combined with a short `every` interval on a case
that accumulates many consecutive accepted events (e.g. a jet growing from a
single bead over a long run), this cadence was found to compound a
nonphysical cross-section thinning across events -- see [Numerical
robustness of the cross-section fit](#numerical-robustness-of-the-cross-section-fit)
below. The warning is advisory only: the historical Example 5 and Tests
21-23 references all use thresholds below this recommendation and remain
validated short-window cases.

[`examples/input-5/input.dat`](../../examples/input-5/input.dat) is the
historical refinement case. [Test Case 21](../examples/test-21.md) combines
pre-extended anchors, insertion, Maxwell evaporation, and stochastic Platen
integration in a focused remeshing validation.

## Runtime workflow

The integrator calls `driver_dynamic_refinement` before advancing the next
time step. The driver follows this sequence:

```text
check enabled state and start time
                 |
                 v
wait for the configured check interval
                 |
                 v
measure element lengths along the active jet
                 |
                 v
estimate whether remeshing adds useful resolution
                 |
                 v
construct an anchored target mesh
                 |
                 v
Akima-interpolate continuous bead fields
                 |
                 v
reconstruct and renormalize conserved quantities
                 |
                 v
update bounds, allocation flags, and dependent workspaces
```

When bead insertion is active, refinement is deferred until at least ten
active beads exist. A threshold crossing alone is also insufficient: the
estimated fitted mesh must contain more points than the current active mesh.
After the configured `every` interval has elapsed, the threshold is tested on
each subsequent timestep until a useful refinement is accepted. The interval
counter is reset only at that accepted event. Thus `every` is a minimum delay
between remeshes, not a permanently sparse sampling period. The counter
`irefinementdone`, available as the `nref` observable, increases only when a
refinement is accepted.

## Path parametrization and anchors

`compute_length_path` assigns each bead a normalized path coordinate `jetpt`.
The remesher uses this coordinate, rather than an individual Cartesian axis,
so the same algorithm works for one- and three-dimensional trajectories.

The refined region extends from the collected/nozzle side through
`irefbeadstart`, the last element found above the threshold. Beads marked in
`jetbd` divide it into anchored segments. Each segment receives a number of
new intervals equal to the ceiling of its arc length divided by the requested
`dynamic refinement threshold`. The base insertion `resolution` does not
define this target mesh. Anchor locations are copied into the target mesh
exactly; the unrefined tail is retained rather than interpolated.

Anchor beads reduce interpolation drift and preserve selected material
locations, but their array indices may still change when the active interval
is rebased. They are ordinary jet material unless the separate variable-mass
feature is explicitly enabled; an anchor flag alone does not represent a
nanoparticle or freeze a bead between timesteps.

For the historical one-bead startup, anchors continue to be assigned as new
material is inserted at the nozzle. A pre-extended initial jet now receives
interior anchors before its first timestep, at the configured anchor spacing.
The two path endpoints are not tagged because they already occur in both
meshes; treating an endpoint as an interior anchor would create a zero-length
anchored segment.

## Interpolated and reconstructed state

Akima interpolation is applied to the continuous fields on the new path
mesh, including:

- coordinates and velocities;
- viscoelastic stress;
- mass and charge densities;
- cross-sectional radius;
- evaporation cross-sectional radius when evaporation is enabled.

Stress, mass, charge, and radii are forced non-negative after interpolation.
Discrete flags such as bead anchors and breakup markers are mapped separately
rather than treated as continuous scalar data. Collected-bead flags are
rebuilt from the remeshed coordinates.

Before interpolation, mass and charge are converted to densities using the
reference material volume. After remeshing they are converted back to bead
quantities. This avoids interpreting an old per-bead amount as though it
belonged unchanged to a new bead volume.

## Volume conservation and evaporation

Reference volume `jetvl` is reconstructed from the refined element length
and interpolated cross section. The refined portion is then renormalized so
its total reference volume equals the pre-refinement total.

With evaporation, instantaneous post-evaporation volume `jetve` follows the
same procedure independently, using the evaporation radius `jetce`. Thus the
algorithm preserves the two totals separately:

- `jetvl` describes reference material volume;
- `jetve` describes current volume after solvent loss.

They must not be merged, substituted for one another, or normalized with a
shared factor. Their local distribution may change slightly through spline
interpolation even though each total is conserved.

After the evaporated-volume rescale, `enforce_evlim_conservative` restores
the same `evlim` floor (`jetve(i)/jetvl(i)>=evlim`) that every ordinary
per-timestep integrator path already enforces, but which the reconstruction
above did not check before 2026-08-16: an Akima undershoot on the
evaporation-radius field could otherwise leave a bead's evaporated volume
far below its own floor (observed: to about 1/7 of it, uniformly across
nearly the whole active jet in one event), inflating the polymer mass
fraction `cp` well past 1 and, in the Maxwell evaporative stress equation,
driving the stress derivative toward overflow. The fix is a water-filling
redistribution: beads below their floor are clamped up to it, and the
resulting deficit is taken back from the still-compliant beads in
proportion to their surplus, iterated in case that step creates a new
violator, so the segment's evaporated-volume total is preserved to
roundoff rather than broken (as a naive clamp alone would). If the floor is
intrinsically infeasible against the segment total, the code aborts
(`error(20)`) instead of silently violating either constraint.

## Numerical robustness of the cross-section fit

The cross-section radius `jetcr` (and the evaporation radius `jetce`) is the
one field, among the 11 interpolated per accepted event, that is squared
back into a conserved quantity (`jetvl=length*pi*jetcr**2`, and via
`jetms=density*jetvl`, into bead mass). A long single-nozzle-bead-start run
exercising several dozen consecutive accepted refinement events
(`examples/input-24`, untracked, see its local `STATUS.md`) exposed that
this specific reconstruction can compound a nonphysical thinning across
*successive* events, well beyond anything electrospinning's genuine
order-of-magnitude fibre thinning would produce, eventually collapsing bead
mass toward zero and overflowing the stress equation -- reproducing with
evaporation on or off, so it is independent of the `evlim` fix above.
Position/velocity/stress do not show the same defect: segment length
(derived from interpolated positions) fluctuates non-monotonically and
shrinks only modestly over the same run where radius shrinks relentlessly
and monotonically by two orders of magnitude more, and stress/velocity only
show large excursions in the single event immediately preceding a crash,
consistent with a late-stage consequence of the mass collapse (via Newton's
second law) rather than an independent defect in their own interpolation.

Five layered defenses were added to the `radius_area`/`evap_radius_area`
handling in `fit_jet_akima` (`dynamic_refinement_mod.f90`) and `fit_mod.f90`,
each validated against Tests 21-23 with zero regressions, none of which by
itself has eliminated the failure on this specific long run:

1. Fit `ln(pi*jetcr**2)` (log cross-section area), not the raw radius or
   area. Volume is quadratic in radius but linear in area, so fitting area
   removes one amplification step; its logarithm additionally guarantees a
   strictly positive recovered value by construction (a raw area fit was
   observed to go measurably negative before an unprincipled `dabs()`) and
   linearises the roughly exponential thinning trend, which is much better
   conditioned for Akima's tangent estimate.
2. `despike_median_filter`: a 3-point median filter over the fitted target
   segment, skipping anchor points (which must remain exact), removing an
   isolated single-event interpolation spike.
3. `clamp_to_local_source_range`: bounds each non-anchor target point to the
   range spanned by its two immediate bracketing pre-event source points --
   a cheap local approximation of a shape-preserving interpolation
   guarantee that classic Akima does not provide. Confirmed unable, by
   itself, to catch a drift already compounded through earlier events,
   since that drift is by then baked into the "legitimate" source data too.
4. `clamp_to_anchor_envelope`: bounds interior points between two
   bracketing *anchors* instead, with a generous (100x) margin, since
   anchors are restored exactly every event. Found insufficient in
   isolation because an anchor can itself already be tagged from an
   already-degraded bead -- anchors are frozen at tagging time, not
   validated against anything.
5. `limit_akima_tangents_monotone` (`fit_mod.f90`): a Fritsch-Carlson
   sufficient-condition monotonicity limiter on the Akima knot tangents,
   opted into via `setup_akima`'s optional `lmonotone` argument only for
   `field_name` in `{radius_area, evap_radius_area, mass_density,
   charge_density}` (never position/velocity/stress, which may have
   genuine local extrema). This delayed the failure the most of any single
   change tried, but the compounding drift still eventually recurs across
   enough events -- possibly a cumulative resampling/requantization effect
   distinct from classic Akima overshoot, not yet isolated.

An absolute-floor variant (`enforce_radius_floor_conservative`, same
water-filling structure as the `evlim` fix, default 1 nm) is also present
but was found empirically too permissive for this input: the measured
collapsing radii (order 5-40 nm) never actually cross an absolute
nanometre-scale floor, since electrospinning legitimately reaches sub-micron
fibre radii. It remains as a defense-in-depth against a genuine sign
crossing the other layers might miss, and as the `error(21)` infeasibility
guard, but a floor relative to the jet's own current scale (rather than an
absolute physical constant) would be needed to actually engage here, and
was not implemented. See `docs/STATE.md` for the full investigation,
including the confirmed-ruled-out alternative explanations and the
guarded diagnostic instrumentation kept in the tracked source for any
future continuation.

The practical fix for this specific long run turned out to be its
refinement **cadence**, not something the five interpolation-side defenses
above could fully absorb by themselves. The tight `0.10` cm threshold
(5x resolution) combined with a short `1.d-5` s `every` interval let a jet
growing from a single bead accumulate several dozen accepted events within
under `2e5` timesteps -- each re-fitting the cross-section from the
previous event's own output. Reusing Example 5's coarser, already-validated
cadence instead (`threshold 0.4` cm, 20x resolution; `every 1.d-3` s) let
the same input run past `6.5e7` of its `1e8` timesteps (324 accepted
events) with the minimum bead radius measured at each event staying
constant at its initial value the entire time, and the active-bead count
settled into the stationary insertion/removal oscillation the case was
designed to reach. This is the cadence now shipped with
`examples/input-24`.

**GPU/OpenACC scope: host-only.** This whole investigation targeted
GFortran CPU only. Of the fixes above, the log-area transform and the
three despiking/local-source/anchor-envelope calls apply on both backends
(they wrap `fit_akima` in `fit_jet_akima` regardless of which backend it
dispatches to), but `enforce_evlim_conservative`,
`enforce_radius_floor_conservative`, and `limit_akima_tangents_monotone`
are reachable only from the host paths
(`reconstruct_refinement_state_host`, `setup_akima`) and were not ported to
their OpenACC counterparts (`accelerator_reconstruct_refinement_state`,
`fit_akima_accelerator`). A GPU build would still reproduce the `evlim`/`cp`
blowup and would miss the monotonicity limiter's share of the cadence
mitigation.

## OpenACC path

The serial Maxwell/Platen configuration in Test Case 21 keeps integration,
insertion, and the jet state in a persistent OpenACC data region. Once the
minimum refinement interval has elapsed, a device reduction returns only the
last over-threshold segment index, total path length, and the nozzle-length
correction. These scalars reproduce the inexpensive part of the historical
acceptance test without downloading bead arrays.

Only when that test predicts a denser target mesh does JETSPIN download the
active topology and evaporation state. The host computes normalized path
coordinates and target knots. Akima then follows four accelerator phases:

1. one independent secant slope per source interval;
2. serial extrapolation of the four endpoint slopes;
3. one independent Akima tangent and one cubic coefficient set per knot; and
4. one independent interval search and interpolation per target knot.

The endpoint extrapolation is constant work; the three knot-oriented phases
are OpenACC parallel loops. They contain no reduction or order-dependent sum.
Source and target coordinates are uploaded once per accepted event. Each of
the 11 source fields is then uploaded once and its interpolated result is
downloaded once.

The interpolated cross-section radius then feeds a second device kernel that
reconstructs bead volume/evaporation volume, rescales each for
reference-volume conservation, and converts the interpolated mass/charge
densities back to per-bead quantities -- the same volume-conservation and
mass/charge-conversion rules described above, executed on the device instead
of the host. Only the data-dependent bookkeeping that builds the normalized
target mesh in the first place (the mass-boundary walk, `jetptc`/`jetbdc`
assembly) and the anchor save/restore/final-assembly steps remain host-side:
they are a small, rare, inherently sequential scan, not a per-bead parallel
operation. After those host checks, the remeshed state is uploaded and device
integration resumes.

The normal Test 21 allocation reserves one `incnpjet` block and does not need
to grow. `incnpjet` is set by the `dynamic refinement capacity <i>` input
directive (positive integer, default 100 if omitted) and applies uniformly to
the initial reserve and every later growth increment. The capacity-growth
validation instead uses the developer-only `JETSPIN_REFINEMENT_INITIAL_RESERVE`
environment override to reduce that reserve to 50 entries, independently of
the input value. When the accepted mesh exceeds the old capacity, JETSPIN detaches
the topology and evaporation mappings before the host allocations are
replaced. It then resizes the Platen workspaces, semantically repacks the
pre-generated Gaussian history for the new stride, generates values only for
new bead indices, and binds the completed host mesh and resized history once.
The same device Akima kernels operate before and after capacity growth.

Test Case 22 validates repeated use of this lifecycle. Its physical input is
derived from Test 21, while developer-only environment overrides reduce the
initial reserve and refinement growth increment to 20 entries. Three accepted
events must therefore perform three independent mapping releases, device Akima
remeshes, Gaussian-history repacks, workspace reallocations, and device
rebinds. Normal runs still reserve and grow by `incnpjet`, which defaults to
100 and is otherwise set by the `dynamic refinement capacity <i>` input
directive.

Test Case 23 adds collector removal without changing this hybrid boundary.
The device topology primitive advances the lower active bound and clears the
removed bead without downloading the full state. Accepted refinement events
still perform the only intermediate complete-state downloads. The test
requires removal on both sides of the third remesh, proving that a rebased
active interval survives the Akima and device-rebind lifecycle.

Two development builds validate this boundary. `JETSPIN_DEV_HOST_AKIMA`
retains the historical host coefficient/interpolation path as an oracle.
`JETSPIN_COMPARE_AKIMA` executes the host reference and the GPU path at each
accepted event and reports coefficient and interpolated-value errors for all
11 fields. Test 23 performs 33 such comparisons; its largest relative errors
are `2.48e-15` for coefficients and `3.32e-15` for values. The largest
absolute value difference is `7.28e-12`.

A third development build, `JETSPIN_COMPARE_REFINEMENT_ASSEMBLY`, applies the
same host-then-device oracle pattern to the volume/conservation/density
kernel: back up the pre-reconstruction state, evaluate the trusted host
reference, restore that state, run the device kernel so it stays
authoritative, then compare. All three Test 23 events agree with the host
reference at or near roundoff (worst absolute `7.1e-15`, worst relative
`4.05e-16`). Building this oracle caught two real ordering mistakes before
they could reach the standard path: an early version invoked the device
kernel a second time after the comparison helper already ran it, silently
reapplying the density-to-mass conversion; a later version evaluated the host
reference from the already-device-converted state instead of the shared
pre-reconstruction input. Both corrupted the event's mass/charge invariant by
tens of percent -- exactly the class of error this A/B methodology exists to
catch before a device kernel is trusted.

## Allocation and MPI interaction

Refinement can increase the number of beads beyond `mxnpjet`. In that case
`define_akima_bounds` grows capacity to the required upper bound plus the
normal `incnpjet` reserve. The remesher may also compact and rebase the active
interval without growing capacity.

The resulting flags have different meanings:

- `doreorder` records that bead indices were remapped;
- `doallocate` records that capacity-dependent consumers may require larger
  arrays;
- `lneighlistdo` requests a new neighbour list for multiple-timestep Coulomb
  calculations.

The fundamental jet state is replicated on every MPI rank. Interpolation
uses root-generated service data and collective communication to reproduce
the remeshed arrays. All ranks must leave refinement with identical active
state, bounds, and feature flags. Rank-local chunks and cached integrator
workspaces must then match the new `inpjet`, `npjet`, `mxnpjet`, and
`mxchunk` values.

See [dynamic allocation and bead indexing](dynamic-allocation.md) and
[MPI parallelization](parallelization.md) for the two supporting memory
models.

## Backup and recovery

The module provides grow-on-demand backup arrays through `store_backup` and
`restore_backup`. The saved state includes active bounds, simulation time,
fundamental bead fields, reference volume, and `jetve` when evaporation is
enabled. Restore broadcasts the recovered state so every rank resumes from
the same configuration.

Any new fundamental bead field that affects the equations must be considered
for both backup and restore, in addition to interpolation and allocation.

## Maintenance invariants

Changes to refinement should preserve these rules:

- the target path coordinates remain ordered and valid for Akima fitting;
- anchored and unrefined regions retain their intended physical locations;
- mass, charge, and volume are treated consistently as densities or
  per-bead quantities at each stage;
- `jetvl` and `jetve` are conserved and normalized independently;
- continuous fields are interpolated, while logical or identity fields are
  explicitly remapped;
- optional arrays are handled symmetrically on every allocation path;
- inactive array entries are cleared before collective summation;
- capacity growth updates MPI chunk sizing and all cached workspaces;
- code never assumes `inpjet == 0` or stable bead array indices;
- restart and backup formats are updated when new persistent state is added.

## Validation

The standard smoke suite exercises the historical refinement example and a
three-way Kelvin--Voigt, evaporation, and refinement path. The latter forces
both an Akima remeshing event and growth beyond the initial bead capacity.
Test Case 21 adds the pre-extended anchored Maxwell/Platen path and checks the
event-level remeshing invariants.
Test Case 22 extends that check to three consecutive events and three capacity
increases.
Test Case 23 combines those three events with collector removal before and
after the final remesh.

Run:

```sh
tests/smoke/run.sh serial
tests/smoke/run.sh debug
tests/smoke/run.sh mpi
tests/refinement/run.sh
tests/refinement/run.sh openacc
tests/refinement/run.sh force-oracle
tests/refinement/run_test22.sh nvfortran
tests/refinement/run_test22.sh openacc
tests/refinement/run_test22.sh force-oracle
tests/refinement/run_test23.sh nvfortran
tests/refinement/run_test23.sh openacc
tests/refinement/run_test23.sh force-oracle
tests/refinement/run_test23.sh host-akima
tests/refinement/run_test23.sh akima-compare
tests/refinement/run_test23.sh refinement-compare
```

For numerical changes, also compare conservation totals immediately before
and after refinement, verify that `nref` becomes non-zero, check for ordered
path coordinates and finite state, and compare serial and MPI trajectories
within an explicitly chosen tolerance.

The five cross-section-fit defenses in the previous section were each
validated the same way -- `tests/refinement/run.sh gfortran standard` and
direct GFortran runs of `examples/input-22`/`examples/input-23` -- with zero
regressions, but they do not yet make the long single-nozzle-bead-start
`examples/input-24` run to completion; see `docs/STATE.md` for the open
investigation.
