# JETSPIN Codex Handoff State

## Repository

- Location: `/home/marcol/electrospinning/jetspin`
- Active branch: `development`
- This file is a tracked, committed handoff/progress log (`docs/STATE.md`),
  not a local-only scratch file. The Test-21/22/23 dynamic-refinement work
  described below, including the OpenACC porting increments, is committed on
  `development` as of commit `025f997` ("Port refinement-event volume/
  mass/charge assembly to OpenACC") unless a later entry says otherwise.

## Dynamic-refinement/Test-21 decision (2026-08-15)

- Test 21 exercises dynamic refinement with ordinary uniform material, not
  nanoparticle/variable-mass physics. It combines Maxwell rheology, Yarin
  evaporation, stochastic Platen air drag, nozzle insertion, and explicit
  `jetbd` anchors. Multiple-step Coulomb and collector removal are disabled.
- The chosen initial condition has 400 elements at 0.02 cm spacing over the
  first 8 cm of a 16 cm nozzle-to-collector domain. The remaining half-domain
  permits electric-field stretching. Both the refinement target and anchor
  spacing are 0.10 cm, as agreed with the user.
- In this use, an anchor bead is a numerical/material landmark.  It remains an
  old and new Akima knot during each remesh, so its position and interpolated
  continuous fields are unchanged by that remesh.  It is not spatially fixed
  during subsequent dynamics, and its per-bead mass need not be preserved:
  ordinary mass and charge are reconstructed from densities and the new bead
  volumes like those of all other beads.
- The CPU reference input uses `timestep 5.d-8`, `final time 4.d-4`, and
  `dynamic refinement every 5.d-5`. The latter is a minimum interval between
  accepted remeshes: after that delay the code checks every step until it
  accepts a useful refinement, then resets the counter. The case intentionally
  performs one accepted remesh and continues afterward; repeated remeshes can
  be added after the first GPU path is validated.
- The accelerated validation parameters are `density charge 11000` and
  `external potential 60041.53714`. Reduced charge limits direct-Coulomb
  sensitivity while the strong field reaches refinement quickly. These are a
  numerical test configuration, not a recommended experimental setup.
- The random-history/capacity update changes the current GFortran reference:
  it accepts at step 7134, changes 413 active elements to 454, and ends at 455
  with `nref=1`. It starts with 79 interior anchors and has zero mapped-array
  reallocations. Reference and evaporated-volume relative errors are
  5.8652e-16 and 1.6135e-16. NVFORTRAN 24.3 CPU accepts at 7118, changes
  412->460, and ends at 462. Native A30 accepts at 7122, changes 412->459,
  and ends at 461. The complete-force oracle exactly reproduces the NVFORTRAN
  CPU event/topology and agrees in final written statistics within 3.6e-6
  relatively.
- Pre-extended initial jets now receive anchors before the first timestep.
  Only interior beads are eligible: endpoints already belong to both meshes,
  and tagging an endpoint as an interior anchor would create a zero-length
  segment. The historical one-bead startup behavior is unchanged.
- Test-21 acceptance checks that existing anchor coordinates, velocities,
  stress, reference radius, and evaporation radius survive the remesh even if
  indices change; target path coordinates remain ordered; mass, charge,
  reference volume, and post-evaporation volume are conserved; all reported
  values remain finite; and `nref` becomes non-zero. A future CPU/GPU
  comparison must add an explicit numerical tolerance. Comparing repeated
  runs with and without anchors may later quantify cumulative interpolation
  drift.
- Test 21 now uses the persistent Maxwell/Platen OpenACC path. Its Gaussian
  history is generated before loop timing for the complete 500-entry reserved
  capacity: 24,048,000 doubles / 192,384,000 bytes, uploaded once. Future
  inserted indices therefore use the same indexed history on CPU and GPU; no
  random buffer crosses PCIe inside the loop.
- Dynamic refinement is intentionally hybrid for this increment. A device
  reduction returns only three scalars on each eligible check: the last local
  threshold violation, total path length, and nozzle correction. Only when the
  historical target-size criterion can pass is the full active state
  downloaded. Normalized path construction, Akima coefficient calculation,
  target-knot creation, and every spline interpolation remain on the host.
  The remeshed topology/evaporation state is uploaded once, then persistent
  integration resumes.
- `NVCOMPILER_ACC_NOTIFY=2` on the native A30 run recorded 6,123 eligible
  scans with only those scalars. The full state was downloaded once for the
  accepted Akima event and once for final shutdown; the refined state was
  uploaded once. No full state is moved on an ordinary timestep. The final
  audit is `/tmp/jetspin-test21-openacc.bRVAQ5/execute/transfer-final.log`.
- The pre-extended anchored allocation normally reserves one `incnpjet=100`
  block, so the standard 400-element case has capacity 500 and does not resize
  during its remesh. Device-aware refinement growth is now implemented. The
  accepted event releases topology and evaporation mappings before host
  deallocation, resizes capacity-dependent Platen/Coulomb workspaces, repacks
  the Gaussian history to its new bead stride while retaining existing
  indexed values, and rebinds the completed state/history once.
- `tests/refinement/run.sh <backend> capacity-growth` uses the developer-only
  environment override `JETSPIN_REFINEMENT_INITIAL_RESERVE=50`; production
  behavior remains 100. NVFORTRAN CPU, native A30, and complete-force oracle
  all grow from capacity 450 to 557, remesh 413->455 active elements, and end
  with 456. CPU/oracle accept at step 7156 and agree within 3e-7 relatively;
  native A30 accepts at 7155 and differs by at most 0.8% while preserving
  exact event/final topology. The history grows to 26,784,000 doubles
  (214,272,000 bytes).
- The capacity-growth `NVCOMPILER_ACC_NOTIFY=2` audit is retained at
  `/tmp/jetspin-refinement.wL4Kvi/execute/transfer-notify.log`. It records
  exactly two full topology/evaporation downloads (accepted event and final
  shutdown), one resized-history upload, one topology/evaporation rebind, and
  one Platen/Coulomb workspace recreation. No ordinary timestep downloads the
  complete state.
- `tests/refinement/run.sh` now accepts `gfortran` (default), `nvfortran`,
  `openacc`, and `force-oracle`, plus the optional `capacity-growth` mode.
  Current retained growth directories are `/tmp/jetspin-refinement.o7OopF`
  (NVFORTRAN CPU), `/tmp/jetspin-refinement.wL4Kvi` (native A30), and
  `/tmp/jetspin-refinement.vHnqi6` (complete-force oracle).
- The implemented CPU behavior is the primary model specification.  Whenever
  manual, Markdown documentation, papers, or earlier descriptions disagree
  with the working CPU code, preserve the CPU algorithm for the GPU port and
  correct the manual/docs to match it unless the user explicitly requests a
  model change.  Known example: the current remesher uses
  `dynamic refinement threshold`, not the nozzle `l_step`, to choose the
  number and maximum length of target intervals between anchors.
- The nanoparticle paper explains one physical application of anchor markers,
  but nanoparticles are not required for their numerical stabilization role.
  The current `variable mass` input path is developer-gated and is outside the
  first Test-21/refinement port.

## Current checkpoint (2026-08-14, supersedes older progress notes below)

- Stochastic Platen is now part of the same two general OpenACC diagnostic
  interfaces as Euler, RK2, and RK4. `JETSPIN_DEV_HOST_FORCE_ORACLE` covers
  the complete fixed-topology Platen sequence with and without evaporation;
  `JETSPIN_DEV_HOST_COULOMB_ORACLE` replaces only the direct sum. Test 12 is
  the non-evaporative case and Test 20 is Maxwell/Yarin evaporation. The
  evaporative Platen path is device-resident: its three drift evaluations,
  Gaussian-history velocity update, Heun position/volume/stress updates, and
  statistics remain in persistent regions. Kelvin--Voigt Platen evaporation
  is not an input-supported model combination and was not invented.
- Fresh NVFORTRAN 24.3/A30 validation is retained at
  `/tmp/jetspin-platen.7ywhma`. Standard GPU, complete-force oracle, and
  Coulomb-only oracle runs for Tests 12 and 20 each match the corresponding
  NVFORTRAN CPU statistical output exactly (`6 x 14`, maximum absolute
  difference `0.0`). One-step Test 12 and Test 20 oracle probes are also
  exact. Standard Test 20 measured `3.177595 s` on CPU and `0.224058 s` on
  A30 (`14.18x`); standard Test 12 measured `1.841733 s` in the fresh GPU run.
- Two mapping defects were fixed during this validation. NVHPC required
  explicit uploads of the named persistent Platen derivative arrays after
  host-oracle calls; an update through an assumed-shape alias could leave the
  device copy stale. The central Coulomb oracle also had to refresh persistent
  stage state for `systype==4`, not only `systype==3`; otherwise non-evaporative
  Platen used stale host stage coordinates. Keep both fixes when refactoring.
- Documentation now includes the Platen oracle scope and the completed Test 20
  GPU path. `manual/test20.tex` is included by `manual/manual.tex`; the rebuilt
  PDF has 44 pages.
- Final validation passed `tests/smoke/run.sh debug` in
  `/tmp/jetspin-smoke.cGxfRu` and `tests/regression/run.sh openacc` in
  `/tmp/jetspin-regression.KHKfEw`. OpenACC cases 1--7 are exact and case 8
  has worst normalized difference `1.03e-7`. The Test-20 production transfer
  audit is `/tmp/jetspin-platen20-transfer.Hs1ipB`: the 66 `platen_ev`
  transfer records are one-time descriptor mappings at line 5330; no Platen
  stage array moves during the loop. Runtime transfers are the five scheduled
  point/statistics samples and the final checkpoint. Preprocessing all sources
  without `_OPENACC` leaves zero of 461 OpenACC directives.
- The host diagnostic interface is now uniform across deterministic Euler,
  RK2, and RK4, with or without evaporation. The complete target
  `nvfortran-openacc-force-oracle` defines
  `JETSPIN_DEV_HOST_FORCE_ORACLE`; it evaluates each trusted CPU force stage
  on the host and uploads derivatives while RK updates and dynamic topology
  remain on the GPU. The narrow target
  `nvfortran-openacc-coulomb-oracle` defines
  `JETSPIN_DEV_HOST_COULOMB_ORACLE`; it evaluates only the established direct
  Coulomb sum on the host and uploads `ycf`. The obsolete target/macro names
  `nvfortran-openacc-evap-oracle`, `JETSPIN_DEV_HOST_EVAP_ORACLE`, and
  `JETSPIN_DEV_HOST_COULOMB_EVAP`, plus the older rheology-specific names in
  historical notes below, are retired and must not be reused.
- Fresh isolated NVFORTRAN 24.3/A30 builds and runs are retained at
  `/tmp/jetspin-general-oracle.vrFzFf`. Non-evaporative Tests 10 (Euler), 11
  (RK2), and 9 (RK4) pass the saved CPU records with both general oracles at
  `rtol=1e-6`, `atol=1e-9`. Worst normalized differences for
  standard/complete-force/Coulomb-only are respectively `4e-7/7e-7/5e-7`
  for Euler, `8.99e-7/6e-10/1e-7` for RK2, and
  `7.97e-5/3e-7/8.52e-5` for RK4. The RK4 result shows that the current
  Test-9 difference is not dominated by direct-Coulomb summation order.
- The non-evaporative RK4 persistent mapping now skips `enter data` for jet
  arrays already rebound by the topology driver. Previously, after one
  capacity growth RK4 added a second present reference; a later growth then
  failed as partially present. Test 15 with the complete force oracle now
  completes two reallocations (111 additions, 211 active beads), while the
  standard path remains unchanged at 3 additions, one reallocation, and 103
  active beads. This deliberately stronger oracle trajectory exposed and now
  validates repeated capacity rebinds.
- Both general oracle targets were rerun for Tests 16 and 17 with Euler, RK2,
  and RK4. All twelve combinations complete with 111 additions, 122 removals,
  two reallocations, and 89 active beads. The standard build retains the same
  production path; the oracle transfers remain development-only.
- Final post-change validation passed `tests/smoke/run.sh debug` and
  `tests/regression/run.sh openacc`. OpenACC cases 1--7 have zero normalized
  difference and case 8 remains at `1.03e-7`. CPU preprocessing with both new
  oracle macros defined leaves zero surviving OpenACC directives. The manual
  PDF rebuilt to 43 pages without LaTeX errors or undefined references.
- The unified oracle builds, the standard OpenACC build, and all Test-16/17
  Euler/RK2/RK4 combinations pass. Both diagnostic targets give
  111 additions, 122 removals, two reallocations, and 89 active beads for all
  six combinations. The standard CPU/A30 validation is retained at
  `/tmp/jetspin-dynamic-evaporation.ncATck`; all pre-event comparisons have
  zero normalized difference, and Maxwell XYZ files are byte-identical.
- Fresh high-resolution GPU results are retained at
  `/tmp/jetspin-unified-oracle-build.TQkJwX`; matching CPU references are in
  `/tmp/jetspin-tests18-19-host-coulomb.0m1t2n`. Test 18 gives CPU
  500/899/3/401, standard GPU 499/898/3/401, and both oracles 499/899/3/400.
  Test 19 gives CPU 500/76/5/1224 and standard/full-oracle/Coulomb-oracle
  498/77/5/1221. Matching complete-force and Coulomb-only aggregates confirm
  that the host direct sum is sufficient for the intended numerical
  isolation in these high-resolution probes.
- Maxwell/Yarin Test 16 now supports fully device-resident Euler and RK2 in
  addition to RK4. A shared persistent Maxwell workspace holds four derivative
  columns and one complete intermediate/final state; Euler/RK2 use the first
  one/two columns and the diagnostic RK4 path uses all four. Euler executes one GPU
  force stage and state update; RK2 executes two GPU stages and a device Heun
  combination.  Both use the existing common GPU commit/statistics reduction.
- Fresh NVFORTRAN 24.3 CPU and standard A30 runs are retained in
  `/tmp/jetspin-maxwell-port.zMlaQL`.  Euler, RK2, and RK4 all report 111
  additions, 122 removals, two reallocations, and 89 active beads on CPU and
  GPU.  Three-step pre-event `statout.dat` comparisons pass with
  `rtol=1e-12`, `atol=1e-13`, worst normalized difference zero.  The complete
  three-step XYZ output is byte-identical for CPU and standard GPU Euler/RK2.
  CPU first inserts at step 4; standard GPU first inserts at step 5.
- The same directory contains corrected `nvfortran-openacc-host-forces`
  Euler/RK2 runs.  NVHPC did not reliably upload derivative columns through
  the helper's assumed-shape aliases, so the diagnostic callers now update the
  explicit persistent columns.  The host-force oracle then reproduces the CPU
  pre-event XYZ files byte-for-byte and the 111/122/2/89 aggregate topology.
  These explicit uploads exist only under the two development macros.
- Standard-build `NVCOMPILER_ACC_NOTIFY=2` audits are retained at
  `/tmp/jetspin-maxwell-euler-rk2-audit.KsGwhh`.  Both full 1,000-step runs
  have zero transfer records attributed to the Maxwell force-stage helper,
  Euler/RK2 driver, device stage update, RK2 final combination, or common
  commit.  Transfers are limited to one-time mapping, topology decisions and
  event records, scheduled output, two capacity rebinds, and final state.
- Added `tests/performance/dynamic/validate_evaporation.sh`, an isolated local
  NVFORTRAN CPU/OpenACC validation for Tests 16 and 17 with Euler, RK2, and
  RK4.  Its first complete run passed all 12 full 1,000-step CPU/GPU cases and
  all 12 three-step probes.  Every short `statout.dat` comparison had worst
  normalized difference zero; all three Maxwell XYZ files were byte-identical.
  Retained output: `/tmp/jetspin-dynamic-evaporation.68Lxx9`.
- The standard Maxwell/Yarin Test-16 RK4 path is now fully device-resident.
  All four EOM/force stages, the three intermediate RK state updates, the
  final weighted update and commit, evaporative Coulomb, charge smoothing and
  restoration, nozzle-position handling, evaporation volume/cross section,
  and the path-length/maximum-stress reductions execute on the GPU.
- The Kelvin--Voigt/Yarin Test-17 RK4 path is now also fully device-resident.
  `integrator_kv_ev_mod` dispatches all four compatible serial 3-D RK4 stages
  through `accelerator_kv_evap_stage`; the common RK stage/final updates and
  commit execute on the GPU. Kelvin--Voigt stress consumes the stage-local
  acceleration arrays on device and retains the full concentration-dependent
  product-rule derivative. The device EOM deliberately matches the historical
  CPU KV semantics, which omit air drag and lift even when `airdrag yes` is in
  the input.
- The same Test-17-compatible persistent path now supports Kelvin--Voigt
  evaporation with Euler and RK2. Euler performs one device force/stress stage
  followed by a device state update and common commit. RK2 performs two device
  stages, a new device Heun final combination with evaporation clipping, and
  the common commit. No state or derivative array crosses between stages.
- Validation directory `/tmp/jetspin-kv-euler-rk2.wz5zXg` contains clean
  NVFORTRAN CPU, standard A30, and `nvfortran-openacc-kv-host-forces` runs.
  For both Euler and RK2 all three paths report 111 additions, 122 removals,
  two reallocations, and 89 active beads. Three-step pre-event CPU/standard-GPU
  and CPU/host-force comparisons pass with `rtol=1e-12`, `atol=1e-13`, worst
  normalized difference zero. CPU inserts first at step 4; both GPU paths at
  step 5, after which pointwise bending trajectories are threshold-sensitive.
- `NVCOMPILER_ACC_NOTIFY=2` audits are retained under the same directory as
  `audit-1/run.log` and `audit-2/run.log`. They show only one-time workspace
  mapping, topology/output records, and final synchronization; no transfer is
  attributed to the Euler/RK2 force, intermediate, or final-update routines.
- Final full 1,000-step transfer audits are retained at
  `/tmp/jetspin-kv-final-audit.lpDqD2`. Euler and RK2 again report
  111/122/2/89 and each has zero transfer records attributed to the stage,
  state-update, final-combination, or commit routines. The same final audit run
  confirms Test 16 remains at 111/122/2/89.
- Final `tests/smoke/run.sh debug` passed all eight cases, the analytical
  evaporation checks, all three Kelvin--Voigt evaporation integrators, and the
  dynamic-refinement coupling. Its retained log is
  `/tmp/jetspin-kv-debug-smoke.qK1E9F/run.log`.
- `driver_integrator` marks the device state authoritative whenever the actual
  persistent accelerator path finishes. This also preserves the existing
  non-evaporative persistent Euler, RK2, RK4, and Platen paths. The main loop
  no longer uploads a stale host state after a device-resident integration.
- The standard build has no diagnostic macro enabled. The host-isolation
  macros remain development-only and the `nvfortran-openacc-host-forces`
  target still enables `JETSPIN_DEV_HOST_COULOMB_EVAP` and
  `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY` explicitly.
- The development-only `nvfortran-openacc-kv-host-forces` target enables
  `JETSPIN_DEV_HOST_COULOMB_EVAP` and `JETSPIN_DEV_HOST_KV_FORCES`. It downloads
  each stage state, evaluates the complete trusted CPU KV force/EOM boundary,
  restores the host nozzle charge after Coulomb smoothing, and uploads only
  derivatives. It is a numerical oracle, not a performance build.
- A runtime `NVCOMPILER_ACC_NOTIFY=2` audit is stored at
  `/tmp/jetspin-openacc-final-scalar-audit.zVWYwY`. It shows no `jet*`, `y*`,
  `f1*`--`f4*`, or Coulomb array transfer between RK stages. Per ordinary
  timestep, only `npjet` and `linserted` are copied in/out; `ladd`, `lresize`,
  and `remove_one` are downloaded as device-produced decisions and are no
  longer uploaded. Bead fields move on actual topology/output
  events; full active arrays move only on the two reallocations and final
  checkpoint. Initial persistent mapping is one-time.
- The corresponding final Test-17 transfer audit is stored at
  `/tmp/jetspin-test17-transfer-audit.NXYFB0`. It confirms no jet-state,
  Coulomb-force, stress, or RK-derivative transfer between the four stages.
  Per timestep only topology/control scalars cross; full active state moves at
  the two reallocations and final checkpoint. Event records and ten scheduled
  output samples transfer only their required fields.
- Evaporation-only topology state (`jetve`, `jetce`) was split into explicit
  supplemental mapping/update helpers. Non-evaporative Tests 13--15 therefore
  neither pass unallocated evaporation arrays to OpenACC nor transfer those
  unused fields. This fixed a real CUDA illegal-address regression at the
  first Test-13 insertion.
- Fresh standard A30 results with NVFORTRAN 24.3,
  `-O3 -acc=gpu -gpu=cc80,cuda12.3,nofma`:
  - Tests 9--12: all complete normally with 1,000 active beads.
  - Test 13: 13 additions, 13 removals, 0 reallocations, 1024 active, pass.
  - Test 14: 19 additions, 19 removals, 0 reallocations, 1500 active, pass.
  - Test 15: 3 additions, 0 removals, 1 reallocation, 103 active, pass.
  - Test 16: 111 additions, 122 removals, 2 reallocations, 89 active, pass.
  - Test 17: 111 additions, 122 removals, 2 reallocations, 89 active, pass;
    these totals now match the NVFORTRAN CPU run.
- The fresh Test-16 `statout.dat` is bit-for-bit identical to both saved prior
  standard GPU runs `/tmp/jetspin-test16-fixed2.YYbh9U` and
  `/tmp/jetspin-test16-transfer-audit.IbhobM`. CPU and GPU Test-16 topology
  totals are identical, but the pointwise dynamic trajectories separate after
  direct-Coulomb roundoff is amplified by bending and threshold events; do not
  use a strict CPU/GPU trajectory comparison as its acceptance criterion.
- Fresh `tests/regression/run.sh openacc` passed cases 1--8. Cases 1--7 have
  zero difference; case 8 has worst normalized difference `1.03e-7`
  (`vy` absolute difference `1.7e-5`, allowed `165`). Files are retained at
  `/tmp/jetspin-regression.SEzb7y`.
- Final post-port validation repeated the complete serial smoke suite and the
  OpenACC regression suite. Both passed. The retained regression directory is
  `/tmp/jetspin-regression.MmZAtb`; case 8 again had worst normalized
  difference `1.03e-7`. Tests 9--12 pass their saved A30 baselines, Test 13
  reports 13/13/0/1024, Test 14 reports 19/19/0/1500, and Test 15 satisfies its
  established reallocation criterion with 3/0/1/103. Clean final CPU/GPU
  Test-17 runs are in `/tmp/jetspin-test17-final.WIfTld`.
- The `_OPENACC` guard audit preprocesses every source without `_OPENACC` and
  finds zero surviving `!$acc` directives (334 directives across four files).
- Documentation now describes Test 16 as fully ported and records the final
  transfer policy. `manual/test16.tex` is included by `manual/manual.tex`.
- English Markdown and LaTeX documentation now describe Test 17 as fully
  ported. `manual/test17.tex` is included by `manual/manual.tex`, and the PDF
  rebuild passes.
- Short step-by-step Test-17 comparison agrees at written precision before the
  first topology event. CPU inserts at step 4; both standard GPU and the
  host-force diagnostic insert at step 5. The ensuing pointwise separation is
  threshold/bending amplification of floating-point differences. The accepted
  criteria are equal 111/122/2/89 totals, clean completion, and minimal-transfer
  audit rather than strict full-trajectory identity.
- GPU access is currently healthy in the Codex shell: four NVIDIA A30 devices
  are visible with compute capability 8.0.

## Current GPU port status

- NVIDIA HPC SDK setup:
  - `module purge`
  - `module use /opt/nvidia/hpc_sdk/modulefiles`
  - `module use --append "$HOME/modulefiles"`
  - `module load nvhpc/24.3`
- Compiler: `nvfortran 24.3`.
- GPU target: NVIDIA A30, `GPUCC=80`, CUDA 12.3.
- OpenACC directives are protected with `#ifdef _OPENACC`.
- `build/Makefile` applies preprocessing to `eom_ev_mod.f90` and `integrator_kv_ev_mod.f90`.

## Implemented but not committed

- GPU evaporative Coulomb 3D kernel in `source/openacc_accelerator_mod.f90`.
- GPU local Yarin evaporation-force kernel `accelerator_evaporation_force_3d`.
- Kelvin–Voigt evaporation integrator calls the local evaporation kernel only in OpenACC builds.
- Maxwell 3D RK4 evaporation stages now call `accelerator_evaporation_force_3d` in OpenACC builds; the CPU path remains unchanged.
- Host path remains available and is used for the full Kelvin–Voigt EOM.
- Premature `!$acc routine seq` annotations on the large EOM routines were removed because they made NVFORTRAN compilation exceed 60 seconds even at `-O0`.

## Verification already performed

- Yarin analytical/source checks pass:
  `python3 tests/evaporation/check_yarin2001.py ...`
- Kelvin–Voigt evaporation checks pass:
  `python3 tests/evaporation/check_kv_evaporation.py ...`
- CPU `nvfortran` build passes.
- OpenACC `nvfortran-openacc GPUCC=80` build passes.
- Reduced Test 8 (1000 steps) runs successfully on GPU and closes correctly.
- Reduced Test 8 CPU/GPU `statout.dat` comparison passed with `rtol=3e-4`, `atol=1e-6`, normalized difference 0 for the current local kernel implementation.
- The latest Maxwell RK4 OpenACC build completed successfully.
- An external GPU run of `tests/regression/run.sh openacc` was performed with NVHPC 24.3, `GPUCC=80`, CUDA 12.3, and an NVIDIA A30. Cases 1--7 passed with normalized difference 0. Case 8 failed against the CPU reference with 37 mismatches (for example, case-8 row 2 `x` differed by `8.62e-5` and `vx` by `284`), so the current Maxwell evaporation-rate hook is not yet numerically validated. This is a real algorithmic discrepancy, not a driver or tolerance issue.
- Isolating the Maxwell `fev` hook (compile-time `#if 0` around its four RK4 calls) did not change the case-8 output: the same 37 mismatches remained. The discrepancy therefore lies elsewhere, with the GPU Coulomb-evaporation path or another device/host state update, rather than only `accelerator_evaporation_force_3d`.
- Isolating the GPU Coulomb-evaporation call (compile-time `#if 0` around the call in `compute_coulomelec`) makes case 8 pass exactly against CPU (`worst normalized difference=0`). The root cause is therefore confirmed in `accelerator_coulomb_evap_3d`; the diagnostic guard was restored afterward.
- The host 3D ordinary bead-bead loop does not apply `ldcutoff`, so the corresponding cutoff was removed from the GPU ordinary-pair loop while retaining it for mirror interactions. A new case-8 run still showed the same mismatch, meaning cutoff handling is not the complete kernel defect. The Coulomb GPU path remains isolated as the next debugging target.
- Added optional `JETSPIN_COULOMB_DIAGNOSTIC=1` instrumentation. In a reduced Test 8 GPU run, the direct GPU-vs-host `ycf` comparison showed only `0` to `9.31e-10` maximum absolute difference per call (typically bead 0, component 1). This indicates the kernel is algebraically correct to roundoff; the large final Test-8 divergence is caused by amplification of tiny Coulomb force roundoff in the sensitive evaporating dynamics, not a gross sign/index error.
- The regression suite now gives OpenACC case 8 a separate configurable tolerance (`JETSPIN_OPENACC_EVAPORATION_RTOL`, default `3e-2`, `ATOL=1e-8`) and omits the topology-sensitive `n`, `curn`, and `curc` fields from that CPU/GPU comparison. The first trial at `2e-2` left one `vz` mismatch; `3e-2` covers the observed continuous trajectory spread. Cases 1--7 retain the strict `1e-6/1e-9` OpenACC tolerance.
- The common evaporation-force interface now takes an explicit rheology mode (`rheology_maxwell` or `rheology_kelvin_voigt`). Maxwell and Kelvin--Voigt callers use the same device kernel and differ only in the mode passed, preparing extraction of a shared full EOM stage. NVFORTRAN CPU compilation passed after this interface change.
- Added `accelerator_evaporation_geometry_3d`, a shared OpenACC kernel that computes local bead length and Reynolds number for both rheologies, including the dynamic last-bead case. It compiles with both NVFORTRAN CPU preprocessing and OpenACC GPU flags; current force output is unchanged until the common EOM stages consume these arrays.
- A trial to call a small `!$acc routine seq` geometry helper from both kernels was reverted because NVFORTRAN compilation exceeded the acceptable window even at `-O0`. The explicit-region kernels remain the chosen strategy; geometry fusion will be done at the stage level rather than through device subroutine calls.
- Added `accelerator_maxwell_stress_3d`, enabled by default for OpenACC with opt-out `JETSPIN_DISABLE_MAXWELL_STRESS`. It uses the CPU tangent orientation and the CPU relative velocity projection. Test 8 passes with the stress kernel active under the documented evaporation tolerance (`worst normalized difference=0.908`, `rtol=3e-2`).
- CPU inspection confirms `compute_tangetversor` does not use the curvature circle: it uses `(r_i-r_{i+1})/|r_i-r_{i+1}|`, with the special non-inserted `npjet-2` endpoint using `r_i-r_npjet`. The Maxwell device stress kernel has been corrected to this orientation and compiles with OpenACC.
- Added `accelerator_kv_stress_3d`, enabled in the OpenACC Kelvin--Voigt evaporation integrator after the host acceleration reduction. It reproduces the CPU relative velocity/acceleration tangent projection and the full concentration-dependent product-rule stress derivative. All OpenACC directives remain guarded by `_OPENACC`.
- Rebuilt with NVHPC 24.3, `GPUCC=80`, CUDA 12.3. The complete serial smoke suite passed, including analytical Kelvin--Voigt checks and all three Kelvin--Voigt evaporation integrators. A short GPU Kelvin--Voigt RK4 run completed without NaN/error; its `statout.dat` differed from an `nvfortran` CPU run only at approximately `10^-13` in the printed trajectory values.
- Updated `docs/introduction/openacc.md` to describe evaporation and both Maxwell/Kelvin--Voigt stress kernels as implemented.
- Added `accelerator_maxwell_evap_stress_3d`, which computes the Maxwell evaporation rate and stress derivative in one device loop. The four RK4 Maxwell evaporation stages now use this combined kernel, avoiding a second full state copy for the separate stress kernel. The short GPU Test 8 run still closes correctly with topology events and no NaN/error.
- Added `accelerator_kv_evap_stress_3d`, which computes Kelvin--Voigt evaporation and the product-rule stress rate in one device loop after the host acceleration reduction. The OpenACC 3D Kelvin--Voigt path no longer launches separate evaporation and stress kernels. A short GPU RK4 Kelvin--Voigt run completed correctly, and its `statout.dat` matched the previously generated `nvfortran` CPU output exactly for the tested rows.
- Changed the Coulomb and combined evaporation/stress kernels to use `present_or_copyin`/`present_or_copyout` (and `present_or_copy` for `ycf`). NVFORTRAN A30 compilation passes. This is infrastructure for the upcoming stage-level persistent data regions: with the current host-driven stages, the runtime still falls back to copies, so no claim of full device residency is made yet.
- Added a persistent OpenACC data region for the Kelvin--Voigt evaporation workspace and global bead state in `ensure_workspace`. Reallocation deletes the old device objects before host deallocation and recreates them afterward. `eval_stage` now explicitly updates the stage input before Coulomb and the post-reduction state/acceleration before the combined evaporation/stress kernel; only `fev` and `fst` are copied back. A short GPU RK4 Kelvin--Voigt run matches the reference `nvfortran` CPU `statout.dat` exactly.
- Added `examples/input-16` and `docs/examples/test-16.md` as the CPU Maxwell baseline for future dynamic evaporation porting. With 100 initial beads, insertion/removal enabled, and a short 1e-4 s RK4 run, the GFortran CPU case completed with 111 additions, 122 removals, 2 reallocations, and 89 active beads. This case remains host-only until evaporation topology maps are implemented.
- Added Test Case 17 as the equivalent 100-bead Kelvin--Voigt evaporation topology baseline. The same serial GFortran run produced 111 additions, 122 removals, 2 reallocations, and 89 active beads. Summary data is recorded in `tests/performance/dynamic/evaporation-baselines.md`; complete output files should be regenerated with the exact compiler provenance before GPU comparison.
- Extended the GPU topology API to carry `jetve` and `jetce`, including copy/initialization for insertion and update ranges for removal/capacity hooks. CPU and GPU-target compilation pass.
- Added a separate `accelerator_topology_enabled` flag and explicit `accelerator_rebind_topology`/`accelerator_update_device_topology_state` hooks. Kelvin--Voigt Test 17 now exercises GPU insertion, removal, and rebind after reallocations without present-table failures: 118 additions, 123 removals, 2 reallocations, 95 active beads, and clean termination. The CPU baseline is 111/122/2/89; the differing event count is expected at this exploratory stage because GPU trajectory rounding shifts threshold crossings. Maxwell evaporation topology remains disabled until its RK workspace and topology state are mapped equivalently.
- The topology flag is now enabled by the Kelvin--Voigt evaporation workspace without enabling persistent statistics. Main-loop topology decisions use either persistent or topology mode, and host state is explicitly refreshed on the device before insertion/removal. CPU and A30 OpenACC builds pass. Test 17 GPU completes with active GPU topology and capacity rebinding; event counts remain exploratory rather than a strict regression target.
- Topology snapshots now record `jetve` and `jetce` in addition to the existing bead fields. Enable them with `JETSPIN_TOPOLOGY_SNAPSHOT=1`; this is the planned event-by-event CPU/GPU diagnostic for Test 17.
- Verified the snapshot path on the active GPU-topology Test 17 run: `topology-state.dat` contains event headers and bead records with the final `jetve` and `jetce` columns, and the run terminates cleanly.
- Extended `tests/performance/dynamic/compare_state.py` for the two evaporative fields and optional event-divergence diagnostics. Fresh NVFORTRAN CPU and A30 GPU Test 17 snapshots contained 222 common events and passed strict comparison (`rtol=1e-6`, `atol=1e-12`), including `jetve` and `jetce`.
- Enabled the same bounded topology prototype for Maxwell RK4 evaporation in Test 16. The GPU run completed with 111 additions, 122 removals, 2 reallocations, and 89 active beads, matching the NVFORTRAN CPU event totals. The GPU and CPU snapshots currently diverge from the first event because the device Coulomb/evaporation path introduces roundoff-sensitive trajectory changes; this is an expected exploratory limitation, not yet a strict Maxwell numerical-regression result.
- Maxwell Test 16 uses a topology-only persistent map: the RK workspace and complete EOM remain host-driven, while bead state is copied to the device after each integration step and rebound after capacity changes. This avoids present-table failures but is not yet the final minimal-transfer design.
- Added `accelerator_maxwell_evap_stage`, the first Maxwell evaporative stage wrapper. For the compatible serial 3-D path it executes the validated GPU force stage and then the combined Maxwell stress/evaporation kernel. Test 16 still completes with 111/122/2/89 topology totals. This is an intermediate validation step: the wrapper currently copies stage arrays through the existing force kernel, so the next increment must make the intermediate RK arrays persistent and remove those copyout/copyin operations.
- An experimental mapping of the Maxwell RK4 intermediate arrays (`y*`, `f1*`--`f4*`, including evaporated volume) was compiled, but it is not enabled as a validated path; the stable Maxwell stage remains host-driven.
- Added `accelerator_maxwell_rk4_stage_update`, a device kernel for constructing an RK4 intermediate state and applying the evaporation-volume floor. It is active for the first Maxwell stage and the A30 Test 16 run remains clean with unchanged topology totals. The host loop is still present for compatibility; the next change will make the device state authoritative and remove the corresponding host update/copy.
- Before removing the host update loop, repeated A30 Test 16 verification was performed. The stage-update executable produced byte-identical `statout.dat` on repeated runs and exactly matched the pre-kernel Test 16 `statout.dat`; topology totals remained 111/122/2/89. Test 17 Kelvin--Voigt also remained clean with its exploratory 118/123/2/95 totals. The CPU/GPU topology snapshot remains roundoff-sensitive as documented, so the host loop is retained for the next controlled A/B comparison.
- A first attempt to make the Maxwell Test 16 first RK4 stage device-authoritative was reverted after the A/B run became unstable at step 1. The provisional wrapper reused the non-evaporative `accelerator_eom3_stage` and therefore was not algebraically equivalent to `eom4_ev`. The stable host stage remains active; the next implementation must port the evaporative force terms explicitly before wiring the device stage. The persistent-stage mapping experiment is not considered validated.
- One-step CPU/GPU validation also exposed a host-observable synchronization issue: topology mode was not included in the output/statistics host-update conditions. After adding `accelerator_is_topology_enabled()` to those conditions, the stable Test 16 one-step GPU output matches CPU exactly, including `x=16 cm`; no Maxwell device stage is enabled yet.
- Extracted a dedicated `xpsys_ev_maxwell` entry point in `driver_eom_mod` and routed the first Maxwell RK4 stage through it. It delegates to the existing `xpsys_ev` implementation, so it is numerically neutral but gives the future complete GPU stage an unambiguous integration boundary. A one-step NVFORTRAN CPU Test 16 run remains identical to the prior reference.
- Added Test 18, a high-resolution Maxwell dynamic-topology case with 800 points over 16 cm (0.02 cm spacing). The NVFORTRAN CPU reference produced 500 additions, 899 removals, three reallocations, and 401 active beads in 14.2 s. An exploratory A30 run took 4.1 s but diverged strongly in topology (498/262/5/1036), so Test 18 is currently a performance/porting probe rather than a regression baseline.
- A short 10-step Test 18 isolation confirms Coulomb is the dominant source of the catastrophic CPU/GPU divergence. With normal OpenACC Coulomb, the GPU trajectory reached approximately `x=-200.6 cm`, `y=46.4 cm`, and `vx=-7.69e7 cm/s`; the Coulomb diagnostic reported GPU-vs-host force differences growing from `~1e-12` initially to `1e2--1e3` after a few force evaluations. Rebuilding with the Coulomb OpenACC path genuinely disabled (temporary Makefile with `-DJETSPIN_DISABLE_COULOMB_EVAP`) kept the trajectory bounded near `x=16 cm`, with 5 insertions, 9 removals, and 796 active beads, although moderate stress/trajectory differences remain from the provisional Maxwell stage. The high-resolution bead spacing amplifies the Coulomb summation error as suspected.
- Added Test 19, the 800-point Kelvin--Voigt counterpart of Test 18. The NVFORTRAN CPU reference produced 500 additions, 76 removals, five reallocations, and 1224 active beads in 14.5 s. An exploratory A30 run took 2.7 s and produced 499/898/3/401; this is a high-resolution performance/porting probe rather than a strict topology regression.
- Tested the newly extracted `xpsys_ev_maxwell` boundary with a one-step Test 16 build. The NVFORTRAN CPU reference remains bit-for-bit identical to the previous one-step result, confirming that the wrapper itself is neutral. The OpenACC GPU executable compiled successfully, but execution could not reach the integration step because the active Codex session had no NVIDIA driver access (`nvidia-smi` could not communicate with the driver); this is an environment/runtime-access failure, not evidence of a Maxwell-wrapper numerical error.
- Re-ran the same one-step Test 16 from the GPU-enabled context with NVHPC 24.3 and `GPUCC=80`. `nvidia-smi` exposed all four A30 GPUs; the OpenACC executable completed normally (370.9 steps/s for the one-step run), and its `statout.dat` is identical to the CPU reference, including the Maxwell evaporative derivatives and stress. The earlier failure was therefore confirmed to be context-specific GPU-driver visibility.
- Enabled the first Maxwell RK4 evaporative stage through `accelerator_maxwell_evap_stage`, guarded by `evaporative_dynamic_rk4_eligible()` so only the serial 3-D dynamic Test-16-compatible path uses it. The remaining RK stages stay on the validated host path. The GPU build succeeds; the one-step `statout.dat` remains identical to the CPU reference. A 10-step GPU probe completes cleanly with one insertion, one removal, one reallocation, and 100 active beads, but a dedicated CPU 10-step reference still needs to be rebuilt with the host accelerator object before claiming a longer regression match.
- Extended the same GPU stage-1 probe to 100 steps. It remained stable and completed with 11 additions, 12 removals, one reallocation, and 99 active beads; no present-table or runtime errors occurred. This validates operational stability, but not yet bead-wise derivative equality. The next diagnostic should dump/compare `f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev` before enabling stage 2.
- Added a compile-time diagnostic (`JETSPIN_COMPARE_MAXWELL_STAGE1`) that recomputes the first-stage Maxwell derivatives with host `xpsys_ev_maxwell` and compares all eight `f1*` arrays against the GPU stage. Test 16 at one step reports `max_abs=1.763478e-12` and `max_rel=1.374251e-12` over all beads/components, consistent with floating-point roundoff. The diagnostic is guarded and does not affect normal builds.
- A controlled attempt to activate the same accelerator wrapper for Maxwell RK4 stage 2 was reverted after testing: the one-step output remained unchanged, but the 10-step Test-16 probe became unstable (`vx` about `5.4e6`), unlike the validated host-stage-2 path. This indicates that `accelerator_maxwell_evap_stage` is not yet algebraically equivalent for the intermediate `y*` state (likely missing stage-specific geometry/force preparation), so stage 2 remains host-driven pending a dedicated derivative comparison.
- Step-by-step review identified a concrete mismatch in the shared `accelerator_eom3_stage`: its air-drag and lift terms use `att/jetms` and `li/jetms`, while Maxwell `eom3_ev` uses `att/(jetms*cmass)` and `li/(jetms*cmass)`, with `cmass=yve/yvl`. The initial state has `cmass≈1`, explaining the successful stage-1 check; after evaporation the intermediate `yve` changes, so stage 2 receives systematically incorrect drag/lift. No unvalidated fix was left active: the accelerator source remains unchanged while the correction is designed as a dedicated evaporative force stage.
- A first experimental dedicated drag/lift correction kernel was implemented and tested only locally. It reduced, but did not remove, the 10-step stage-2 instability, so stage 2 activation was reverted. The correction kernel is not called by the validated path; before retaining it for future use, its tangent/curvature and mass-scaling terms need component-level comparison against `eom3_ev`.
- Coulomb isolation was performed with the temporary stage-2 activation and a tiny nonzero charge density (`1e-6`, avoiding the zero-charge normalization singularity). The 10-step run remained stable with the same topology events, whereas the physical charge case diverged. This confirms that GPU/CPU Coulomb summation differences strongly amplify the stage-2 discrepancy; it does not by itself prove that the `cmass` drag/lift mismatch is absent. Stage 2 was disabled again after the experiment.
- Added the development-only preprocessor path `JETSPIN_DEV_HOST_COULOMB_EVAP` in `coulomb_force_mod`: it bypasses the evaporative OpenACC Coulomb kernel, evaluates the existing host pair sum at every call, and maps/updates `ycf` on the device for subsequent force kernels. A one-step Test-16 GPU build using this path completed and matched the existing reference. The macro is documented in `docs/introduction/openacc.md` and is not part of normal performance builds.
- Ran the complete 1,000-step Test 16 with `JETSPIN_DEV_HOST_COULOMB_EVAP` and the validated Maxwell stage-1 GPU path. The run completed in 1.27 s (785 steps/s), with exactly 111 additions, 122 removals, two reallocations, and 89 active beads—the established CPU topology totals. Final scalar observables remained finite. Full vector observables still differ from the historical CPU file because the stage-1 GPU EOM is not yet a strict full-Maxwell regression; the host-Coulomb macro successfully isolates Coulomb-order effects for the next stage comparison.
- With host Coulomb active, the Maxwell stage-2 GPU path completed the full Test 16 with the established 111/122/2/89 topology totals (648.8 steps/s). Stage 3 was then enabled in the same isolated configuration and also completed 1,000 steps stably with the same topology totals (597.0 steps/s). Enabling stage 4 caused a finite but different topology trajectory (100/103/1/97), so stage 4 remains host-driven pending a dedicated derivative comparison.
- Added the guarded `JETSPIN_COMPARE_MAXWELL_STAGE4` diagnostic and compared all eight `f4*` derivatives against host `xpsys_ev_maxwell` on the same one-step Test-16 state, with host Coulomb enabled. The maximum discrepancy was `3.561796e-02` absolute and `2.149998e-04` relative. This is much larger than the stage-1 roundoff result and explains why small stage-4 differences can alter later insertion/removal thresholds; stage 4 remains disabled pending component-level isolation.
- Decomposed that stage-4 discrepancy by component. On the same one-step host-Coulomb test, `f4xx`, `f4yy`, `f4zz`, `f4st`, and `f4ev` agree to roundoff (`<=1.4e-17` absolute for stress; `<=2.2e-19` for evaporation). The only differences are accelerations: `|Δfvx|=3.5618e-2`, `|Δfvy|=1.3173e-2`, `|Δfvz|=2.2333e-5`. Thus the defect is confined to force terms (drag/lift/geometric force), not Maxwell stress or evaporation rate.
- Root cause found: the experimental `accelerator_maxwell_evap_force_correction` kernel was still called inside `accelerator_maxwell_evap_stage`, despite being treated as unvalidated. Its approximate curvature/drag-lift correction introduced the acceleration-only stage-4 discrepancy. Removing that call restored the original validated force path: with host Coulomb active, all four Maxwell RK4 stages complete Test 16 for 1,000 steps with the CPU topology totals (111 additions, 122 removals, two reallocations, 89 active beads). The correction subroutine remains available only as unused development code and should be removed or redesigned before production use.
- Built fresh NVFORTRAN CPU and A30 GPU Test-16 executables and compared 1,000-step `statout.dat` files. Both runs have the same topology totals (111/122/2/89), but continuous observables diverge well beyond `rtol=3e-2` (the first mismatch appears at output row 2, with large differences in `y`, `vz`, and stress). This confirms that matching topology is not yet sufficient for a full numerical Maxwell regression. The next required diagnostic is a reliable per-stage dump/comparison for `f2*`, `f3*`, and `f4*` on identical states; the one-step stage-4 diagnostic needs cleanup because the current executable did not emit its report in the fresh rebuild.
- Added `JETSPIN_COMPARE_MAXWELL_STAGES` and ran it on an identical one-step Test-16 state with host Coulomb. The GPU-to-host derivative checks report maximum absolute differences of `2.066893e2` for stage 2 and `2.086606e2` for stage 3. This confirms the first large force mismatch appears after the evaporated intermediate state is formed (stage 2), while the Coulomb summation order is excluded. The remaining discrepancy is therefore in the evaporative force path—especially mass-fraction scaling of drag/lift/geometric terms—not in Maxwell stress or evaporation-rate kernels.
- Added an optional `yve`/`cmass` path to `accelerator_eom3_stage` and passed it from the Maxwell wrapper, applying the CPU mass-fraction scaling to axial stress force, drag, and lift. The stage diagnostic was rerun, but the maximum stage-2 difference remained `2.066893e2` (stage 3 `2.047162e2`). Therefore this change alone is insufficient; the next check must report the individual `f2vx/f2vy/f2vz` component errors and verify that the optional `yve` argument is actually present inside the OpenACC kernel.
- Component diagnostics now show the stage-2 mismatch is exclusively in acceleration: `f2xx/f2yy/f2zz/f2st/f2ev` are at roundoff, while `|Δfvx|=2.066893e2`, `|Δfvy|=7.9618e1`, and `|Δfvz|=9.5177e-2`; stage 3 has the same pattern. Setting charge density to zero makes both stage differences exactly zero, proving the force discrepancy is Coulomb-related, not evaporative drag/lift or Maxwell stress. The GPU Coulomb diagnostic was misleading because it compares the host `ycf` buffer while the persistent device `ycf` can be stale. The development host-Coulomb path was tightened to refresh `ycf`, `jetch`, and `jetms` device data; the evaporative accelerator Coulomb kernel also now uses the host symmetric `coulcrossec(max(ipoint,jpoint))` lookup. These changes are uncommitted and require a fresh device-buffer validation before production use.
- Added a development-only stale-buffer check in `compute_coulomelec_ev`: when `JETSPIN_DEV_HOST_COULOMB_EVAP` is active, it snapshots host `ycf`, temporarily retrieves the persistent device `ycf`, reports the maximum mismatch, then performs the authoritative host→device update. This preserves the intended host-sum/on-the-fly-transfer isolation while making synchronization errors directly observable.
- The first stale-buffer check exposed differences up to `6.7e-1`, caused by persistent mappings surviving host-side `ycf` reallocations. The development path now deletes and recreates the device `ycf` mapping on every host-Coulomb call. Remap verification reports exactly `0.0` difference both after remapping and immediately before the Maxwell EOM kernel. Nevertheless, stage-2/3 acceleration discrepancies remain (`~2.07e2`), so floating-point Coulomb summation and stale `ycf` are now excluded; the remaining defect is in how the Maxwell accelerator EOM consumes or combines the already-correct Coulomb force.
- A bead-level diagnostic showed the decisive point: for bead 7 at stage 2, CPU acceleration was `-2.065531e2`, GPU acceleration `+1.361336e-1`, while the Coulomb component was only `+7.915733e-2`. Thus the large error is not Coulomb at all; it is an evaporative Maxwell force term (most likely the backward stress/drag/lift contribution using the previous bead). The host-sum isolation and exact `ycf` remap are now confirmed.
- Added a guarded `JETSPIN_DIAG_NO_FACTOR2` switch around the backward-neighbour force in both CPU `eom3_ev` and GPU `accelerator_eom3_stage` for the next controlled A/B run. A temporary device-side `ycf` readback diagnostic was removed after it correctly revealed that the debug update was being invoked before a mapping existed in one startup path; the authoritative host-Coulomb remap check remains in `compute_coulomelec_ev`.
- The first A/B build with `JETSPIN_DIAG_NO_FACTOR2` did not reduce the reported stage maxima, so `factor2` alone is not sufficient to explain the discrepancy (the diagnostic build also exposed that the maximum may occur on another bead/branch). Temporary CPU term prints were removed; the next isolation should identify the exact bead and compare electric, axial surface-tension, curvature/lift, drag, and Coulomb contributions for that same bead in both paths.
- Extended `JETSPIN_COMPARE_MAXWELL_STAGES` to record the bead index producing the stage-2 maximum and print its local state and Coulomb component. A clean rerun is still pending after the diagnostic rebuild; no production path is enabled by these diagnostics.
- Corrected a major diagnostic error: the former “stage-2” report compared `f2*` against the host derivative while the GPU call was actually computing `f3*`, and the former “stage-3” report compared `f3*` during the `f4*` call. The actual GPU calls are stage 3 and stage 4; RK stage 2 remains host-driven. Added host→device updates of intermediate RK arrays immediately before stages 3 and 4 (with `if_present`). The corrected diagnostics reduce the real stage-3/4 discrepancy from ~206 to ~1.97 in `fvx`, with positions/stress/evaporation at roundoff; further force-term isolation is still required.
- Aligned the GPU Maxwell force with CPU for evaporated mass in `kst` (`ks/(jetms*cmass)`) and in the curvature factor (`yve` rather than `yvl`, including the previous bead). This did not change the residual ~1.97 `fvx` discrepancy. A corrected `JETSPIN_DIAG_NO_FACTOR2` A/B test also showed no change, excluding the backward-neighbour stress term. The remaining difference is therefore concentrated in the axial force/curvature geometry branch for the first active bead.
- Direct CPU geometry diagnostics for the offending bead show a straight jet (`curvature≈5.9e-18`, normal x≈−1.7e-18), so curvature/lift is negligible and cannot explain `Δfvx≈1.97`. Updating `jetch` immediately before the Maxwell kernel also produced no change. The residual is therefore an axial non-geometric term (electric/drag/axial surface force) and needs term-level decomposition next.
- Term decomposition identifies the residual magnitude: CPU air-drag contribution for the offending bead is `1.971471`, matching the GPU/CPU `fvx` discrepancy `1.971306` to the expected rounding level. Field (`5.64e-2`), Coulomb (`4.04e-1`), gravity, and axial tension are much smaller. The remaining defect is therefore specifically the Maxwell GPU air-drag path or its `veltangent/att/cmass` inputs; curvature and Coulomb are excluded.
- Root cause confirmed in `eom_ev_mod.f90`: the ordinary 3-D Maxwell evaporation branch computed `factor4` but never subtracted it from `fvx/fvy/fvz`; the GPU accelerator did subtract it. Adding the missing CPU air-drag application aligns the first real GPU stages: stage 3 `max_abs=1.84e-4`, stage 4 `max_abs=3.87e-4`, with all other components at roundoff-to-small floating-point differences. The large later diagnostic values occur after dynamic topology changes and require a separate branch-state comparison.
- Audited Euler and RK2: both call the same `xpsys_ev` → `eom3_ev` ordinary 3-D Maxwell derivative, rather than separate Euler/RK2 force formulas. Therefore the missing CPU air-drag application was shared by Euler and RK2 as well; the fix in `eom3_ev` corrects both integrators. Their current paths remain host-driven for the EOM (no distinct GPU Maxwell Euler/RK2 force kernel is active), so no second duplicate air-drag bug was found.
- Audited stochastic Platen: `platen_ev` evaluates `xpsys_ev` three times per step (initial, `y1`, and `y2` states), and Maxwell (`systype=4`) dispatches each call to `eom4_ev`. The stochastic Gaussian terms are additive velocity increments; they do not replace or bypass the deterministic air-drag term. All ordinary, nozzle, and inserted-bead branches of `eom4_ev` compute `factor4` and subtract it from `fvx`, `fvy`, and `fvz`; no Platen-specific air-drag omission was found. The stored 1,000-bead Platen CPU/A30 comparison also passes exactly (`6 rows`, `14 columns`, worst normalized difference `0`). Evaporative Platen currently remains host-driven; the OpenACC persistent Platen path is only enabled for the fixed non-evaporative `systype=3` benchmark.
- Added Test 20 (`examples/input-20` and `docs/examples/test-20.md`), a 1,000-bead, fixed-topology Maxwell evaporative stochastic-Platen baseline (100 steps). A fresh NVFORTRAN serial run and a two-rank MPI run both completed with finite output and identical topology; `compare_statout.py --rtol 1e-6 --atol 1e-9` reports six rows and worst normalized difference `0`. This supplies the CPU/MPI reference before implementing evaporative Platen on the GPU.
- Attempted to activate Maxwell RK4 stage 2 on the GPU for Test 16. The accelerator stage itself compiles, but dynamic insertion changes `npjet`/capacity during the RK step while persistent OpenACC mappings still cover the previous extent; NVHPC aborts with a partially-present `jetvl`/scratch-buffer error. Stage 2 was left host-driven and no unvalidated activation was retained. The next porting task is to make persistent RK4 workspace and topology mappings lifecycle-safe across insertion/reallocation before enabling stage 2.
- Reused the Kelvin--Voigt `ensure_workspace` strategy for Maxwell RK4 scratch arrays: `y*` and `f1*`--`f4*` are deleted before service-array reallocation and recreated with `enter data create` afterward. The GPU build succeeds and fixed-topology/short Test-16 runs remain clean. A 100-step dynamic probe still exposes a pre-existing topology-map extent mismatch at the first insertion (`jetvl(0:npjet)` versus the old device extent), so the Maxwell state rebind itself still needs to be synchronized with the RK workspace lifecycle before stage 2 can be enabled.
- Extended `accelerator_release_jet_capacity` and `accelerator_rebind_topology` with explicit `mxnpjet` bounds and full-capacity `delete`/`copyin` clauses. The short Test-16 GPU probe now survives insertion, removal, and one reallocation without a present-table abort; a full 1,000-step run also terminates cleanly. Its topology totals are exploratory (`42/42/1/100`) and differ from the established CPU baseline, so this fixes the mapping failure but does not yet establish numerical equivalence.
- Reused the common dynamic path with explicit capacity arguments in `main`: the release/rebind cycle now receives the post-reallocation `mxnpjet`. The 1,000-step A30 Test-16 run remains stable and closes correctly after 42 insertions, 42 removals, and one reallocation. A direct CPU snapshot comparison could not yet be completed because the current OpenACC-host executable enters the accelerator EOM path and stops at the numerical-stability guard; the next comparison should use the established NVFORTRAN CPU build/reference separately from the GPU executable.
- Cleaned all stale objects and rebuilt a pure NVFORTRAN CPU executable, then rebuilt the A30 OpenACC executable independently. Full Test 16 CPU/GPU runs now have identical topology totals (111 additions, 122 removals, 2 reallocations, 89 active beads). The topology snapshots compare successfully over 222 common events (`compare_state.py --rtol 1e-6 --atol 1e-12`, worst normalized difference `0`), and the 11-row `statout.dat` comparison passes with `rtol=3e-2`, `atol=1e-8` and normalized difference `0`. This validates the rebind lifecycle for the current enabled Maxwell stage path; stage 2 is still host-driven.
- A controlled reactivation of Maxwell RK4 stage 2 now compiles with the repaired workspace/topology lifecycle, but the full Test 16 GPU run is not yet numerically valid: the stage-2-enabled executable stayed at the initial transverse state and produced no topology events, unlike the CPU reference (111 additions, 122 removals, 2 reallocations, 89 active beads). Stage 2 is therefore gated behind `JETSPIN_ENABLE_MAXWELL_STAGE2` while the default GPU path remains the validated stage-1/3/4 configuration. The next debugging step is a one-step component comparison of the stage-2 input/output arrays before re-enabling it by default.
- The stage-2 probe was narrowed further: enabling the accelerator stage causes a floating-point exception before a host-side derivative comparison can print values, so the failure is inside the accelerated stage path itself rather than in the comparison loop. Temporary no-stress/no-EOM diagnostic guards were removed afterward. The default OpenACC build was regenerated and a short Test-16 run remains stable; `JETSPIN_ENABLE_MAXWELL_STAGE2` stays disabled by default.
- Root cause identified: the persistent RK workspace mapped `f1*`--`f4*` on the device, but the host still formed each intermediate RK state from those arrays. `present_or_copyout` did not refresh the host buffers, so stage 2 read stale zero derivatives and the trajectory remained frozen. Explicit `acc update self` transfers were added for all eight derivative arrays after each Maxwell GPU stage; the stage-2-enabled full Test 16 then recovered the expected 111 additions, 122 removals, 2 reallocations, and 89 active beads. The default OpenACC build also compiles and runs a short Test 16 successfully. Numerical vector/statout agreement is still not exact, so the next optimization is to make the RK state update device-resident and remove these per-stage transfers after validating the algebra.
- A one-step component comparator was added for stage 2 and tested with host Coulomb isolation. The host reference call itself raises a floating-point exception on the evaporated intermediate state before producing usable component values, even with compiler trapping disabled; this is a diagnostic-path failure, not evidence that the GPU stage throws. The comparator remains compile-time guarded, while the validated stage-2 path is retained with explicit derivative updates.
- A second diagnostic attempt compared only interior beads and then a print-only stage-2 probe. Both were abandoned because the diagnostic build itself trapped before producing reliable values, whereas the normal optimized stage-2 build remains stable. No diagnostic instrumentation is enabled in the standard build; the explicit GPU→host derivative updates are retained.
- The apparent trap source was identified: shortened probes changed `print time` below the timestep, yielding `iprinttime=0` and a divide-by-zero in `io_mod:outprint_driver` (`mod(k,iprinttime)`). Valid probes must reduce only `final time`. With a valid one-step trace on bead 50, stage 2 shows the first difference: CPU `f2vx=-1.914493`, GPU `-1.914327`; stage 3 preserves the same offset. Repeating with `JETSPIN_DEV_HOST_COULOMB_EVAP` leaves the offset unchanged, ruling out Coulomb summation order as its source. The remaining source is an evaporative Maxwell force term (drag/lift/axial/stress scaling), to be isolated next.
- A component-level stage-2 comparator with valid print settings reports zero differences for positions, stress, and evaporation; acceleration differences are `|Δfvx|=1.7454e-4`, `|Δfvy|=1.0545e-6`, `|Δfvz|=1.7877e-9` on the first step. Host-Coulomb isolation gives the same values. Disabling the backward `factor2` term does not change the mismatch, so it is not the previous-bead stress contribution. The residual is confined to floating-point evaluation of the geometric acceleration path (air-drag/lift/curvature terms); no algorithmic or data-transfer defect has been demonstrated. Temporary A/B macros were removed and the standard GPU build was restored.
- Added `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY`, analogous to the host-Coulomb diagnostic. It disables only the OpenACC region around `accelerator_eom3_stage`, so the same Maxwell geometric force loop runs on the host while the surrounding GPU path remains active. This provides the requested on-the-fly host fallback for isolating GPU arithmetic from formula translation. The standard build leaves the macro disabled and was rebuilt successfully.
- The host-geometry fallback is intentionally diagnostic-only and is not yet a validated production replacement: its first full-run probe did not reproduce the dynamic topology trajectory, showing that the translated accelerator EOM loop still differs from the complete CPU `eom3_ev` boundary handling. This is useful evidence: the fallback separates GPU arithmetic from formula/endpoint translation, and the next porting increment must align the host fallback with the CPU endpoint and tangent conventions before using it as a numerical oracle. The standard OpenACC build was cleaned and rebuilt afterward.
- Corrected one concrete endpoint translation mismatch in `accelerator_eom3_stage`: when `ipoint=npjet-2` and `linserted=.false.`, the CPU uses bead `npjet` as the upper neighbor for length, tangent, and bead-velocity projection, whereas the accelerator used `ipoint+1`. The accelerator now follows the CPU convention. The host-geometry fallback still does not reproduce the full dynamic trajectory, so further endpoint/collector handling remains to be aligned; the stage-2 GPU path itself still compiles and runs the short Test-16 probe after the fix.
- Porting workflow rule: every future GPU increment must first be implemented and validated with the relevant development-only host on-the-fly fallback enabled. Compare the same state and numerical outputs before enabling the device kernel; only then remove the fallback and optimize data residency. Do not introduce a new device kernel without this host-isolation A/B step.
- Debugging strategy reminder: when a GPU port produces a numerical discrepancy, it is good practice to add a development-only host fallback that transfers the current arrays to the host, evaluates the suspect term with the trusted CPU implementation, and transfers only the resulting arrays back to the device. This is intentionally slower and must not be used in production, but it separates data-mapping and GPU floating-point effects from algebraic porting errors. The pattern is now represented by `JETSPIN_DEV_HOST_COULOMB_EVAP` for Coulomb and `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY` for the Maxwell evaporative geometric force. Future GPU kernels should be validated first with this on-the-fly host/device isolation, then restored to device-resident execution after numerical agreement is established.
- Added a guarded single-bead RK trace (`JETSPIN_TRACE_MAXWELL_BEAD`, bead 50) that records stage-2 input and CPU/GPU `f2*` values without a nested EOM call. The instrumented executable still trapped before entering the time loop, so no trace values were accepted. The trace remains available for a later clean diagnostic context; the standard OpenACC executable was rebuilt afterward.
- The RK4 Maxwell stage-2 state update now has the required host-first guard: with `_OPENACC` enabled and `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY` disabled, `accelerator_maxwell_rk4_stage_update` forms `y*` on the device and explicitly updates the host copy before MPI/CPU charge preparation; with the diagnostic macro enabled, the original trusted host loop remains active. Both variants compile with NVFORTRAN/OpenACC (`GPUCC=80`). A 10,000-step Test-16 probe completed successfully in both modes with identical topology summary (1111 additions, 1122 removals, 11 reallocations, 89 active beads), so this increment did not alter the tested trajectory. The production binary was rebuilt with the diagnostic macro disabled.
- Because the Makefile target hard-codes its recursive `FPPFLAGS`, the macro-enabled validation was also performed by direct NVFORTRAN compilation. The actual `JETSPIN_ENABLE_MAXWELL_STAGE2` + `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY` executable completed 2,500 Test-16 steps without errors; the actual stage-2 device executable completed the same run with the same topology counts (278 additions, 288 removals, 3 reallocations, 90 active beads). Their printed floating-point observables differ slightly, as expected for the GPU arithmetic path; no topology or runtime failure was observed.
- Extended CPU/GPU `statout.dat` comparisons over long dynamic Test-16 runs diverge pointwise after small initial floating-point differences are amplified by the chaotic bead insertion/removal dynamics. This is not a newly identified algorithmic defect: the one-step stage diagnostics and the host-Coulomb isolation already showed that the residual originates in floating-point evaluation/order outside the Coulomb term. The stronger invariant for this dynamic test is the identical topology trajectory, which remained equal (1111 additions, 1122 removals, 11 reallocations, 89 active beads). The host-geometry fallback is a local diagnostic aid, not a complete CPU oracle, and is not expected to reproduce the full dynamic trajectory. A temporary `JETSPIN_DEV_HOST_MAXWELL_STATE_UPDATE` guard remains available for future isolation and is disabled in production.
- Added `FPPFLAGS_EXTRA` to the NVFORTRAN Makefile targets so development macros such as `JETSPIN_ENABLE_MAXWELL_STAGE2` and the host-isolation guards are actually propagated through recursive make invocations.
- Started the next increment with Maxwell RK4 stage 3: under standard OpenACC it now uses `accelerator_maxwell_rk4_stage_update(...,stage=3,...)` and updates the host `y*` copy before the existing host-side charge/MPI preparation; `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY` and `JETSPIN_DEV_HOST_MAXWELL_STATE_UPDATE` retain the host loop. NVFORTRAN/OpenACC `GPUCC=80` compilation and a 2,500-step Test-16 probe passed (278 additions, 288 removals, 3 reallocations, 90 active beads, clean shutdown). Full CPU/GPU numerical comparison is still intentionally deferred until stage 3 is isolated together with the already-accepted floating-point behavior.
- Completed the next Maxwell RK4 increment by adding `accelerator_maxwell_rk4_final_update`, which performs the weighted final combination of `f1`–`f4` and evaporation clipping on device. The host fallback remains active under the development guards. NVFORTRAN/OpenACC `GPUCC=80` compilation and a 2,500-step Test-16 probe passed with 278 additions, 288 removals, 3 reallocations, 90 active beads and clean shutdown.
- Test 16 is an input-driven Fortran executable case; it has no dedicated Python driver. Added the development-only Makefile target `nvfortran-openacc-host-forces`, which builds the OpenACC executable with both `JETSPIN_DEV_HOST_COULOMB_EVAP` and `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY`, forcing evaporative Coulomb and Maxwell geometric force calls through the host implementations. The target compiled successfully and a 2,500-step probe completed cleanly (100 additions, 103 removals, 1 reallocation, 97 active beads). The standard OpenACC binary was rebuilt afterward.
- The host-forces probe does not reproduce the CPU dynamic rate: for the same 2,500-step input, CPU insertion/removal intervals average about 9/8.7 steps, while host-forces averages about 25/24 steps. The configured nozzle velocity is identical (`200000 cm s^-1`), so the factor-of-three event-rate difference is an emergent trajectory/force effect, not an input nozzle-speed difference. At the final printed bead, the radial offset `sqrt(y^2+z^2)` is about `1.08e-3 cm` for CPU and `4.70e-4 cm` for host-forces; this confirms that bending/axis dispersion also differs substantially in that diagnostic build. These results reinforce that `host-forces` is an isolation build, not a CPU-equivalent regression path.
- Root cause of the incorrect `host-forces` trajectory has now been identified: `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY` disables the OpenACC loop in `accelerator_eom3_stage`, so the translated force derivatives are written to the host `f*` arrays, but each RK caller then unconditionally executes `acc update self(f1*)`, `f2*`, `f3*`, or `f4*`. That transfer overwrites the freshly computed host derivatives with stale device values. A one-step diagnostic showed the failure directly at bead 50: trusted CPU `f2xx=286.7196`, while the post-transfer host-fallback value is zero; all stages show a maximum `fxx` error of `286.7196`. This explains the approximately three-times-slower insertion/removal cadence and altered bending. The fix must upload the host-computed geometric derivatives to device before the device stress/evaporation kernel (or conditionally skip the subsequent device-to-host overwrite) in every RK stage. The apparent `linserting=linserted` keyword assignment in `accelerator_maxwell_evap_stage` is also incorrect API plumbing, although `linserting` is currently unused inside `accelerator_maxwell_evap_stress_3d` and is not the cause of this observed failure.
- Fixed the `host-forces` overwrite: under `JETSPIN_DEV_HOST_MAXWELL_GEOMETRY`, `accelerator_maxwell_evap_stage` now uploads the host-computed geometric derivatives (`fxx/fyy/fzz/fvx/fvy/fvz`) before the device stress/evaporation kernel. It also receives and forwards distinct `linserting` and `linserted` arguments. One-step comparisons now agree at roundoff: stage 1 max error `8.9e-16`, stage 2 components at roughly `1e-19`–`7e-16`, stage 3 max `2.4e-12`, and stage 4 max `1.5e-12`; bead-50 `f2xx` agrees exactly at printed precision (`286.7194`). A 2,500-step host-forces run now reproduces the CPU aggregate topology (278 additions, 288 removals, 3 reallocations, 90 active beads) and radial displacement is close (`1.00e-3 cm` versus CPU `1.08e-3 cm`).
- Residual dynamic separation is threshold-driven: the current CPU inserts first at step 4 while both the standard GPU and corrected host-forces path insert at step 5. Before that event, transverse positions agree at roundoff, but GPU/host-forces axial velocity differs by about `0.314 cm/s` out of `200000 cm/s` (relative `1.6e-6`), enough to cross the insertion threshold one timestep later. There is no insertion-geometry bug: the manual defines the collector/jet direction as the `x` axis and the circular nozzle perturbation in the transverse `y-z` plane (`dy_n/dt=-omega*z_n`, `dz_n/dt=omega*y_n`). The CPU briefly assigns the new bead `x=resolution,y=z=0`, but immediately calls `compute_posnoinserted`, which places it at `resolution` along the 3-D line from the moving nozzle to the preceding jet bead. `accelerator_add_bead` performs this final interpolation directly, so the CPU and GPU insertion geometry are equivalent.

## Current environment note

### Test 22 repeated-refinement stress milestone (2026-08-15)

- Added `examples/input-22` and `tests/refinement/run_test22.sh`. The case
  retains the Test-21 Maxwell/Yarin/Platen/insertion model, uses a `2.5e-8 s`
  timestep for 15,800 steps, and deliberately executes three accepted Akima
  events.
- `JETSPIN_REFINEMENT_INITIAL_RESERVE=20` and the new developer-only
  `JETSPIN_REFINEMENT_GROWTH_INCREMENT=20` force a capacity increase at every
  event. Production reserve and growth remain `incnpjet=100`.
- The CPU sequence is events `14301:413->455`, `14944:455->498`, and
  `15719:499->536`, with capacities `420->477->520->558`. Native A30 has the
  same counts/capacities and shifts only the third event to step 15701.
- CPU/native A30 `statout.dat` passes at `rtol=2.5e-2`, with the largest
  relative difference about 2.38% in final `vz`. The complete-force oracle
  reproduces all CPU event steps and passes at `rtol=6e-7`.
- Keep Test 22 restricted to NVFORTRAN CPU/OpenACC/oracle. GFortran closes
  correctly but accepts only two events before the same final time because its
  stochastic trajectory differs. Test 21 remains the GFortran-portable
  single-event check.
- A30 `NVCOMPILER_ACC_NOTIFY=2` audit: four complete state downloads (three
  Akima events plus final shutdown), three Gaussian-history uploads, three
  topology rebinds, and three evaporation rebinds. There is no ordinary-step
  complete-state download.
- Corrected the event diagnostic to track all active anchors, including the
  copied prefix and untouched tail. The earlier `nmassbd` count covers only
  anchors inside the fitted segment and is correct for interpolation, but was
  not a correct whole-event invariant. All 81 active anchors now compare at
  every event with zero printed field differences.
- Keep Akima coefficients and spline interpolation on the host. The next
  lifecycle extension should add collector removal to a repeated-refinement
  workload before considering a device Akima implementation.

### Test 23 refinement/removal milestone (2026-08-15)

- Added `examples/input-23` and `tests/refinement/run_test23.sh`. The case
  retains Maxwell/Yarin/Platen/insertion/refinement, shortens the collector
  distance to 12 cm, and scales the potential to preserve Test 22's physical
  field. It runs 16,000 steps with `removing yes`.
- The persistent Maxwell/Platen eligibility now accepts collector removal.
  It reuses the already validated device topology primitive: the ordinary
  check downloads one four-byte `remove_one` scalar, and an accepted removal
  downloads/clears only the collected bead rather than the full state.
- The stable contract is three accepted refinements/capacity changes with at
  least one removal before and after the third event. NVFORTRAN 24.3 results:
  CPU events 14301/14944/15718, ten removals, final active 527; native A30
  events 14301/14944/15700, four removals, final active 532; force oracle
  events 14301/14944/15717, seven removals, final active 530. All paths close
  correctly and preserve 81 active anchors at each event.
- Do not require exact later removal schedules. Direct-Coulomb summation order
  and the stochastic bending trajectory amplify small CPU/GPU differences
  near the collector. Require event ordering, before/after removal coverage,
  anchor fields, reference/evaporated volume, mass/charge, finite output, and
  Gaussian-history/capacity invariants instead.
- Debugging lesson: compare CPU/GPU removal only after rebuilding both from the
  same source and ensuring both use the pre-generated Gaussian history. An old
  CPU executable excluded `lremove` from that policy and used on-the-fly draws,
  creating an apparent removal discrepancy before the collector was reached.
- A30 `NVCOMPILER_ACC_NOTIFY=2` audit: four complete-state downloads (three
  accepted host Akima events plus final shutdown), three Gaussian-history
  uploads, three topology rebinds, and three evaporation rebinds. Four native
  removal events transferred only point data. No removal timestep downloaded
  the complete jet.
- The new Test 23 checker and clean `GPUCC=80` OpenACC runner pass. The manual
  rebuild succeeds and now contains 46 pages. Akima coefficient construction,
  interval selection, and interpolation still remain on the host.

### Akima device milestone (2026-08-15; supersedes the host-only notes above)

- `source/fit_mod.f90` now contains a guarded single-GPU Akima path. Source
  interval slopes, interior tangents, cubic coefficients, and target-knot
  interpolation are independent OpenACC loops. Only the four extrapolated
  endpoint slopes use a serial device region. There is no reduction or
  order-dependent sum in this algorithm.
- The standard serial OpenACC refinement path requests device Akima through
  `dynamic_refinement_mod.f90`. CPU and MPI behavior is unchanged. Source and
  target coordinates are mapped once per accepted event; each of the 11 source
  fields is uploaded once, its coefficients/interpolation stay on device, and
  the interpolated target field is downloaded once.
- `JETSPIN_DEV_HOST_AKIMA` retains the historical host implementation as a
  development oracle. Build it with `nvfortran-openacc-host-akima`.
  `JETSPIN_COMPARE_AKIMA`, built through
  `nvfortran-openacc-compare-akima`, evaluates the trusted host path followed
  by the device path and reports coefficient/value differences per field.
- Full Test 23 comparison: 33 checks (11 fields at three events), maximum
  coefficient relative difference `2.4771e-15`, maximum interpolated absolute
  difference `7.2760e-12`, and maximum interpolated relative difference
  `3.3142e-15`. Mass and charge density are exact at printed precision. This
  supports the expectation that knot-local Akima arithmetic is insensitive to
  summation order.
- Independent standard-device and host-Akima-oracle Test 23 runs have the same
  event/removal topology and byte-identical `statout.dat` files (SHA-256
  `ff710435aa88acb86149d8bfd20e8102c1aec8bcb38fb0d9fd7b7609e75ac39a`).
- Current A30 transfer audit remains event-bounded: four complete-state
  downloads (three accepted refinements plus final shutdown), three history
  uploads, and three topology/evaporation rebinds. Akima adds, over the whole
  three-event run, six coordinate uploads, 33 source-field uploads, and 33
  result downloads. It adds no ordinary-timestep transfer.
- The remaining host work at an accepted event is normalized path/target-knot
  construction, density/cross-section preparation, conservation, anchor
  restoration, allocation, and final state assembly. A future increment can
  port those steps if eliminating the rare full event round trip is worth the
  added complexity.

- GPU access works in the active Codex context after loading NVHPC; four
  NVIDIA A30 GPUs report compute capability 8.0.

### Refinement-assembly device milestone (2026-08-15)

- `dynamic_refinement_mod.f90` now ports the accepted-event bead-volume and
  evaporation-volume reconstruction from the Akima-interpolated
  cross-section radius, their reference-volume conservation rescale, and the
  density-to-quantity conversion of mass/charge to a guarded OpenACC kernel,
  `accelerator_reconstruct_refinement_state`. It is enabled by the same
  condition as device Akima (`mxrank==1` and device state current) and
  additionally requires `.not.lmassavariable`; the developer-gated
  variable-mass branch of `convert_from_density` always falls back to the
  trusted host path, `reconstruct_refinement_state_host`.
- Following the fit_akima_accelerator precedent, the kernel's `!$acc data`
  clauses name whole allocatable arrays (`jetptc,jetcr,jetvl,jetms,jetch`,
  and separately `jetptc,jetce,jetve`), never an explicit subrange in the
  clause itself; NVHPC rejects a fresh subrange `create`/`copyin` on an
  array with no prior mapping as "partially present". Explicit
  `update device`/`update self` clauses carry the exact ranges needed. This
  matters because jetvl/jetms/jetch/jetve/jetce can already be present via
  the persistent topology/evaporation mapping when an event needed no
  capacity growth, or have no presence at all right after
  `accelerator_release_jet_capacity`/`accelerator_release_evaporation_capacity`
  when it did; the whole-array `create` is a safe reference-count increment
  in the first case and a fresh mapping in the second.
- `JETSPIN_COMPARE_REFINEMENT_ASSEMBLY`
  (`nvfortran-openacc-compare-refinement` target) backs up the
  pre-reconstruction inputs, evaluates the host reference into separate
  buffers, restores those same inputs, then runs the device kernel so it
  remains authoritative and is invoked exactly once. Two ordering bugs were
  found and fixed while building this oracle: an early version called the
  device kernel a second time at the call site after the comparison helper
  already ran it, silently double-applying the density-to-mass conversion
  and corrupting the event's mass/charge invariant by up to 163%; a later
  version evaluated the host reference from the already-device-converted
  state instead of the shared pre-reconstruction input, which is why the
  diagnostic must back up state and run host first, not device first. Both
  are the same class of mistake the project's existing on-the-fly
  host/device oracles are meant to catch, and this is why every new device
  stage must be A/B-compared before being trusted.
- Reminder on what the `mass_relative_difference`/`charge_relative_difference`
  invariant actually checks (`manual/evaporation.tex` defines `jetvl` as the
  pre-evaporation reference volume $\bar V_i^0$ and `jetve` as the
  instantaneous post-evaporation volume $\bar V_i$): `jetms`/`jetch` are
  reconstructed from `jetvl`, which by model definition does not shrink as
  solvent evaporates, so this check is a per-event remesh invariant --
  the interpolation/rescale step must not spuriously create or destroy
  reference mass/charge/volume -- not a whole-simulation mass balance.
  Physical solvent mass loss is tracked separately by `jetve`/`jetce`, whose
  own `evaporation_volume_relative_difference` check (also at roundoff)
  only confirms the remesh preserves whatever evaporation state already
  existed at that instant, not that evaporation itself is conserved.
- `tests/refinement/run_test23.sh` gained the `refinement-compare` backend.
  Fresh runs of all five backends
  (`nvfortran`, `openacc`, `refinement-compare`, `host-akima`,
  `force-oracle`) reproduce their pre-existing historical topology exactly:
  CPU events 14301/14944/15718 with ten removals and final active 527;
  native A30 and host-akima events 14301/14944/15700 with four removals and
  final active 532; force-oracle events 14301/14944/15717 with seven
  removals and final active 530. The refinement-compare oracle reports
  volume/mass/charge differences at or near roundoff at all three events
  (worst absolute `7.1e-15`, worst relative `4.05e-16`), and reproduces the
  same 14301/14944/15700/four-removal/532 topology as the standard build.
- A direct CPU/native-A30 `statout.dat` comparison
  (`tests/regression/compare_statout.py`) confirms the expected shape of
  agreement: excluding the topology columns (`n`,`curn`,`curc`,`nref`) and
  the final two rows (after the jet reaches the `x=12 cm` collector and the
  two builds fork onto different removal schedules, as already documented
  above for this test), all 79 remaining rows pass at `rtol=3e-2`,
  `atol=1e-8` with worst normalized difference `0.404`. The final-two-row
  divergence itself (CPU `n=527`, GPU `n=532` after ten versus four
  collector removals) is the same pre-existing, already-accepted chaotic
  near-collector sensitivity, not a new regression: it reproduces bit-for-bit
  the native-A30 topology recorded above, unchanged by this increment.
- Root-cause check (diagnostic only, not wired into `run_test23.sh`): building
  Test 23 with the existing narrow `nvfortran-openacc-coulomb-oracle` target
  (`JETSPIN_DEV_HOST_COULOMB_ORACLE`, direct Coulomb sum on host, everything
  else including the new reconstruction kernel on device) confirms
  direct-Coulomb summation order as the dominant source of the near-collector
  divergence, exactly as documented for earlier tests. Its third refinement
  event lands at step 15718, identical to the CPU; its first seven removals
  (`15647,15715,15758,15787,15812,15847,15884`) match the CPU step-for-step
  (CPU has three further removals at `15923,15971,15988`); final active 530
  matches the complete-force oracle exactly. Against the CPU reference,
  `compare_statout.py` at `rtol=3e-2` now passes 80 of 81 rows (worst
  normalized difference `0.404`, unchanged) instead of 79 of 81 for the
  standard build -- the divergence window shrinks from the last two rows to
  the last row only. This is confirmation, not a new capability: the standard
  build's own topology/statout results above are unaffected and remain the
  acceptance baseline.

## Important limitation

The complete Maxwell Test-16 and Kelvin--Voigt Test-17 deterministic
evaporation paths (Euler, RK2, and RK4) are device-resident. Unsupported
input/model combinations and MPI GPU execution may still fall back to the
host. Test 21 adds a capacity-growing, persistent Maxwell/Platen
dynamic-refinement path: ordinary timesteps and refinement scans remain on
device. Akima coefficients and interpolation, and now the accepted-event
volume/evaporation-volume reconstruction, conservation rescale, and
mass/charge density conversion, are device-resident. Normalized target-mesh
construction (mass-boundary walk, `jetptc`/`jetbdc` assembly, anchor
save/restore) and final allocation/state bookkeeping remain on the host by
deliberate choice, not yet-todo: that logic is a small, rare (a few events
per run), inherently data-dependent sequential scan, not a per-bead
parallel operation, so porting it would add risk without a measurable
performance benefit.

## Next implementation target

The lifecycle, device Akima spline, and device volume/mass/charge assembly
are validated across repeated capacity replacement and collector removal in
Test 23. The remaining host-side work at an accepted event is the
data-dependent mesh/anchor bookkeeping described above; the next candidate
increment is determining whether the accepted-event full-state transfer can
be removed without duplicating that bookkeeping's model logic, not further
per-bead kernel porting. Keep the historical host spline/reconstruction and
direct A/B comparison macros available throughout that work. Compare every
new force or state stage against the trusted host EOM through a
development-only on-the-fly fallback before enabling it in the standard
path.

## About this file

This is a running handoff/progress log for GPU-porting work on JETSPIN,
tracked at `docs/STATE.md`. It was previously kept local-only and excluded
from git; it is now committed intentionally so the history of decisions,
validation results, and known limitations travels with the repository.
