
#ifdef JETSPIN_GPU_DYNAMIC_PLATEN
module openacc_dynamic_platen_mod

!***********************************************************************
!
!     JETSPIN module extending the OpenACC persistent execution path of
!     the non-evaporative stochastic Platen integrator (integrator 4,
!     system 4, evaporation no) to dynamic bead topologies: ongoing
!     nozzle insertion together with dynamic mesh refinement, instead of
!     the fixed 1,000-bead, insertion-complete benchmark geometry
!     required by fixed_accelerator_geometry in integrator_mod.f90.
!
!     This module is a fork, in the sense requested for this refactor:
!     it lives entirely in a new file and is compiled in and consulted
!     only when JETSPIN_GPU_DYNAMIC_PLATEN is defined, so the original
!     fixed-geometry persistent path in integrator_mod.f90 is left
!     completely unmodified for default builds (macro undefined).
!
!     The eligibility conditions mirror those already validated in
!     integrator_mod.f90 for the sibling dynamic-topology persistent
!     paths: dynamic_rk4_accelerator_eligible (non-evaporative RK4) and
!     dynamic_evaporative_platen_eligible (evaporative Platen). The
!     integrator==4 check present in those siblings is intentionally
!     omitted here: this function is only ever consulted from inside
!     the platen() subroutine, which driver_integrator only calls when
!     integrator==4 already holds, and adding that check would force a
!     circular module dependency (integrator is a public variable of
!     integrator_mod itself, not of nanojet_mod).
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     last modification August 2026
!
!***********************************************************************

 use nanojet_mod, only : systype,levaporation,lKVfluid,npjet,mxnpjet, &
                         inpjet,linserting,lmultiplestep,lairdrag, &
                         lflorentz,luppot,ldragvel,typemass,ltrackbeads, &
                         ltagbeads,lbreakup
 use version_mod, only : mxrank,mystart,myend
 use electric_field_mod, only : nfieldtype
 use dynamic_refinement_mod, only : lrefinement

 implicit none

 private

 public :: dynamic_platen_accelerator_eligible
 public :: dynamic_platen_accelerator_configured

contains

! Escape hatch for A/B benchmarking against the non-persistent path,
! mirroring the existing dynamic eligibility functions in integrator_mod.f90.
 logical function persistent_disabled_by_env()
  implicit none
  character(len=16) :: disable_persistent
  disable_persistent=''
  call get_environment_variable('JETSPIN_OPENACC_DISABLE_PERSISTENT', &
   disable_persistent)
  persistent_disabled_by_env=trim(disable_persistent)=='1'
 end function persistent_disabled_by_env

! Every condition this fork requires EXCEPT how many beads exist right now.
! A run's static setup (system/integrator choice, evaporation off, airdrag
! on, dynamic refinement on, single rank, ...) is fixed from the first line
! of the input file and does not depend on how far the jet has grown yet.
 logical function dynamic_platen_configured_common()
  implicit none
  dynamic_platen_configured_common=systype.eq.4 .and. &
   .not.levaporation .and. .not.lKVfluid .and. mxrank.eq.1 .and. &
   mystart.eq.inpjet .and. myend.eq.npjet .and. linserting .and. &
   .not.lmultiplestep .and. lairdrag .and. .not.lflorentz .and. &
   .not.luppot .and. nfieldtype.eq.0 .and. .not.ldragvel .and. &
   typemass.eq.0 .and. .not.ltrackbeads .and. ltagbeads .and. &
   .not.lbreakup .and. lrefinement
 end function dynamic_platen_configured_common

! Used ONLY to decide, once, before the timestep loop starts, whether to
! pre-generate the Gaussian pool at all (see
! prepare_integrator_random_history in integrator_mod.f90). At that point
! the jet may still be a single bead (npjet==1), so this deliberately does
! NOT gate on npjet/mxnpjet; the pool does not depend on capacity. It does
! not honour JETSPIN_OPENACC_DISABLE_PERSISTENT either (since 2026-10-01):
! that switch only keeps the per-step gate below closed, so persistent and
! non-persistent runs of this build read the same noise.
 logical function dynamic_platen_accelerator_configured()
  implicit none
  dynamic_platen_accelerator_configured=dynamic_platen_configured_common()
 end function dynamic_platen_accelerator_configured

! Used to decide, every step, whether to actually switch platen() into the
! persistent/fused GPU execution path. Here npjet>=100 and mxnpjet>npjet
! DO matter: there is no benefit to mapping device state for a handful of
! beads, and some spare reserved capacity should exist before mapping.
 logical function dynamic_platen_accelerator_eligible()
  implicit none
  if(persistent_disabled_by_env())then
    dynamic_platen_accelerator_eligible=.false.
    return
  endif
  dynamic_platen_accelerator_eligible=dynamic_platen_configured_common() &
   .and. npjet>=100 .and. mxnpjet>npjet
 end function dynamic_platen_accelerator_eligible

end module openacc_dynamic_platen_mod
#endif
