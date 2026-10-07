 module device_step_mod

!***********************************************************************
!
!     JETSPIN module for the device-resident timestep of the OpenACC
!     build (work plan of 2026-10-06 in docs/STATE.md).
!
!     The strategy of Tests 24 and 25 for every run: a size-independent
!     configuration test (device_step_configured), a gate that opens when
!     the jet has device_min_beads beads and stays open
!     (device_step_eligible), the CPU build's code below it, and above it
!     one step on the device-resident state.  The step is written once:
!     a stage evaluator (nozzle charge smoothing and inserting bead,
!     Coulomb sums, model kernel) and the stage combinations of the
!     scheme, with the evaporated volume as an option.  Milestone M1:
!     Euler, RK2 and RK4 for the 3-D Maxwell model with and without
!     evaporation; it replaces the persistent branches of eulsys, rk2sys,
!     rk4sys and rk4sys_ev and the eulsys/rk2sys/rk4sys_maxwell_ev_device
!     drivers of integrator_mod.  Milestone M2: the Kelvin-Voigt model
!     with and without evaporation (the former with the device drivers of
!     integrator_kv_ev_mod until then), and runs without air drag.
!     Milestone M3: the Platen scheme, Maxwell with and without
!     evaporation, with or without insertion and refinement, formerly the
!     persistent branches of platen and platen_ev with their own gates
!     (fixed 1000-bead geometry, dynamic runs with refinement only).
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!
!***********************************************************************

 use version_mod,       only : mystart,myend,mxchunk,mxrank,idrank
 use nanojet_mod,       only : mxnpjet,npjet,inpjet,systype, &
                         dresolution,thresolution,ivelocity,istress, &
                         imassa,icharge,ivolume,collector_h=>h, &
                         airdragamp,noisediff,noisefric, &
                         jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetst, &
                         jetms,jetch,jetvl,jetve,jetce,jetfr,evlim, &
                         cp0,Bev,mev,tev, &
                         evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity, &
                         resolution,linserted,liniperturb,lairdrag, &
                         lflorentz,luppot,pfreq,consistency,findex, &
                         yieldstress,att,fve,gr,ks,li,v,velext,linserting, &
                         lremove,lmultiplestep,lKVfluid,levaporation, &
                         ldragvel,typemass,ltrackbeads,ltagbeads,lbreakup, &
                         compute_posnoinserted
 use electric_field_mod, only : nfieldtype
 use error_mod,         only : error
 use utility_mod,       only : gaussianhistory,gaussianhistorybase, &
                         gaussianhistorywindow,gaussianhistoryvalues
 use coulomb_force_mod, only : smooth_charge,restore_charge, &
                         coulforce,compute_coulomelec_driver, &
                         set_coulomb_accelerator_persistent, &
                         reset_coulomb_accelerator
#ifdef _OPENACC
 use accelerator_mod,   only : accelerator_eom3_stage, &
                         accelerator_maxwell_evap_stage, &
                         accelerator_kv_stage, &
                         accelerator_maxwell_rk4_stage_update, &
                         accelerator_evap_rk2_final_update, &
                         accelerator_maxwell_rk4_final_update, &
                         accelerator_maxwell_commit_state, &
                         accelerator_compute_posnoinserted_3d, &
                         accelerator_freeze_at_collector, &
                         accelerator_mark_device_state, &
                         accelerator_platen_predict, &
                         accelerator_platen_stage_prep, &
                         accelerator_platen_update, &
                         accelerator_platen_end_step, &
                         accelerator_set_async, &
                         accelerator_set_persistent, &
                         accelerator_set_topology_enabled, &
                         accelerator_is_topology_enabled
 use statistic_mod,     only : counterlpath,ncounterlpath,maxstress, &
                         maxstressposx
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
 use driver_eom_mod,    only : xpsys,xpsys_ev,xpsys_ev_maxwell, &
                         xpsys_KV_pos_v,xpsys_KV_st
 use eom_ev_mod,        only : eom3_KV_pos_v_ev,eom3_KV_st_ev
#endif
#endif

 implicit none

 private

 integer, parameter, public :: scheme_euler=1,scheme_rk2=2,scheme_rk4=3
 integer, parameter, public :: scheme_platen=4
! Value of npjet (the index of the nozzle bead) at which the device step
! engages; the gate then stays open for the rest of the run.
 integer, parameter, public :: device_min_beads=100

 public :: device_step_supported
 public :: device_step_configured
 public :: device_step_eligible
 public :: device_step_request_reset
 public :: device_step_engaged
 public :: device_step_restore_engaged
#ifdef _OPENACC
 public :: device_rk_step
 public :: device_platen_step
#endif

! engaged: the gate has opened once, and stays open (npjet can fall below
! device_min_beads when reallocate_jet compacts the arrays); part of the
! restart state.  restored: engaged was read from the restart file and the
! step has not run since.
 logical, save :: engaged=.false.
 logical, save :: restored=.false.
 logical, save :: reset_requested=.false.
 logical, save :: workspace=.false.
 logical, save :: workspace_mapped=.false.
 integer, save :: workspace_mxnpjet=-1
 integer, save :: workspace_mxchunk=-1
! Stage derivatives f(s), by bead offset, and the stage states: y for the
! RK schemes and the first Platen prediction, z for the second.  The
! evaporated-volume arrays exist without evaporation too, unused, so that
! one set of kernels serves both models.  One allocation per array and
! stage: with the stages as columns of two-dimensional arrays (the first
! version of this module) the evaporative Platen step of Test 25 took
! about 4 us more per step, on the host side only.
 type stage_derivatives
   double precision, allocatable :: xx(:),yy(:),zz(:),st(:)
   double precision, allocatable :: vx(:),vy(:),vz(:),ev(:)
 end type stage_derivatives
 type(stage_derivatives), save :: f(4)
 double precision, allocatable, save :: yxx(:),yyy(:),yzz(:),yst(:)
 double precision, allocatable, save :: yvx(:),yvy(:),yvz(:),yev(:)
 double precision, allocatable, save :: zxx(:),zyy(:),zzz(:),zst(:)
 double precision, allocatable, save :: zvx(:),zvy(:),zvz(:),zev(:)
! The oracle builds exchange forces with the host between kernels: no
! asynchronous queue, no fused preparation.
#if defined(JETSPIN_DEV_HOST_FORCE_ORACLE) || defined(JETSPIN_DEV_HOST_COULOMB_ORACLE) || defined(JETSPIN_DEV_HOST_COULOMB_ACTIVE)
 logical, parameter :: host_oracle=.true.
#else
 logical, parameter :: host_oracle=.false.
#endif

 contains

 logical function device_step_supported(scheme)

!***********************************************************************
!
!     The model and option conditions of the device step, whatever the
!     data layout: the schemes and models it implements, and the options
!     its kernels do not cover yet (work plan, milestone M5).  The Platen
!     scheme runs the 3-D stochastic system (systype 4), with refinement
!     and tagged beads; the RK schemes systype 3, without tagged beads
!     (refinement with RK, milestone M5).  For the Platen scheme it also
!     decides the Gaussian pool, once before the loop, in every build and
!     for any number of ranks (prepare_integrator_random_history), so that
!     CPU, GPU, serial and MPI runs read the same noise.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: scheme

  device_step_supported=.not.lmultiplestep .and. .not.lflorentz .and. &
   .not.luppot .and. nfieldtype==0 .and. .not.ldragvel .and. &
   typemass==0 .and. .not.ltrackbeads .and. &
   .not.lbreakup .and. (linserting .or. .not.lremove)
  if(scheme==scheme_platen)then
    device_step_supported=device_step_supported .and. systype==4 .and. &
     .not.lKVfluid
  else
    device_step_supported=device_step_supported .and. systype==3 .and. &
     .not.ltagbeads
  endif

 end function device_step_supported

 logical function device_step_configured(scheme)

!***********************************************************************
!
!     Every condition of the device step except the bead count: the
!     supported models and options, one rank holding the whole jet.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: scheme

  device_step_configured=device_step_supported(scheme) .and. &
   mxrank==1 .and. mystart==inpjet .and. myend==npjet

 end function device_step_configured

 logical function device_step_eligible(scheme)

!***********************************************************************
!
!     The per-step gate: the configuration test and device_min_beads
!     beads, or an earlier engagement (also one before a restart); the
!     Platen scheme also needs the Gaussian pool.
!     JETSPIN_OPENACC_DISABLE_PERSISTENT=1 or JETSPIN_OPENACC_DISABLE_EOM=1
!     keeps the run on the CPU code.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: scheme

  character(len=16) :: disable
  logical, save :: checked=.false.,disabled=.false.

  device_step_eligible=.false.
#ifdef _OPENACC
! Read once: the gate is evaluated at every step, and the environment
! lookup was a measurable part of the host work between two steps.
! Until 2026-10-07 JETSPIN_OPENACC_DISABLE_EOM made the device stage skip
! the equations of motion and integrate stale derivatives.
  if(.not.checked)then
    disable=''
    call get_environment_variable('JETSPIN_OPENACC_DISABLE_PERSISTENT',disable)
    disabled=trim(disable)=='1'
    disable=''
    call get_environment_variable('JETSPIN_OPENACC_DISABLE_EOM',disable)
    disabled=disabled .or. trim(disable)=='1'
    checked=.true.
  endif
  if(disabled)return
  if(scheme==scheme_platen .and. .not.allocated(gaussianhistory))return
  device_step_eligible=device_step_configured(scheme) .and. &
   (npjet>=device_min_beads .or. engaged)
#endif

 end function device_step_eligible

 subroutine device_step_request_reset()

!***********************************************************************
!
!     A capacity change (reallocate_jet, refinement) invalidates the
!     device workspace: rebuild it at the next device step.
!
!***********************************************************************

  implicit none

  reset_requested=.true.

 end subroutine device_step_request_reset

 logical function device_step_engaged()

!***********************************************************************
!
!     Whether the gate has opened, for the restart file (restart state
!     version 2).  Always false in the CPU builds, unless a restart file
!     written by the OpenACC build said otherwise.
!
!***********************************************************************

  implicit none

  device_step_engaged=engaged

 end function device_step_engaged

 subroutine device_step_restore_engaged()

!***********************************************************************
!
!     A restarted run whose gate had opened keeps it open, as the
!     uninterrupted run does: after a compaction the jet can hold fewer
!     than device_min_beads beads.
!
!***********************************************************************

  implicit none

  engaged=.true.
  restored=.true.

 end subroutine device_step_restore_engaged

#ifdef _OPENACC
 subroutine ensure_device_step(scheme,k)

!***********************************************************************
!
!     Map the jet arrays (unless the topology driver already did) and
!     the workspace, and mark the run as device-resident.  Rebuilt after
!     a capacity change.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: scheme,k

  character(len=40) :: label
  integer :: s

  if(workspace .and. workspace_mxnpjet>=mxnpjet .and. &
   workspace_mxchunk>=mxchunk .and. .not.reset_requested)then
    call accelerator_set_topology_enabled(.true.)
    call accelerator_set_persistent(.true.)
    call set_coulomb_accelerator_persistent(.true.)
    return
  endif

  if(workspace)then
! Complete the queued work (and leave the asynchronous queue) before the
! workspace goes.
    call accelerator_set_async(.false.)
    if(workspace_mapped)then
      do s=1,4
!$acc exit data delete(f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy, &
!$acc& f(s)%vz,f(s)%ev)
      enddo
!$acc exit data delete(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
!$acc& zxx,zyy,zzz,zst,zvx,zvy,zvz,zev)
      workspace_mapped=.false.
    endif
    do s=1,4
      deallocate(f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy,f(s)%vz, &
       f(s)%ev)
    enddo
    deallocate(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
    deallocate(zxx,zyy,zzz,zst,zvx,zvy,zvz,zev)
    call reset_coulomb_accelerator(coulforce)
  endif

  do s=1,4
    allocate(f(s)%xx(0:mxchunk),f(s)%yy(0:mxchunk),f(s)%zz(0:mxchunk))
    allocate(f(s)%st(0:mxchunk),f(s)%vx(0:mxchunk),f(s)%vy(0:mxchunk))
    allocate(f(s)%vz(0:mxchunk),f(s)%ev(0:mxchunk))
    f(s)%ev(:)=0.d0
  enddo
  allocate(yxx(0:mxnpjet),yyy(0:mxnpjet),yzz(0:mxnpjet),yst(0:mxnpjet))
  allocate(yvx(0:mxnpjet),yvy(0:mxnpjet),yvz(0:mxnpjet),yev(0:mxnpjet))
  allocate(zxx(0:mxnpjet),zyy(0:mxnpjet),zzz(0:mxnpjet),zst(0:mxnpjet))
  allocate(zvx(0:mxnpjet),zvy(0:mxnpjet),zvz(0:mxnpjet),zev(0:mxnpjet))
  yev(:)=0.d0
  zev(:)=0.d0

  if(.not.accelerator_is_topology_enabled())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet),jetzz(0:mxnpjet), &
!$acc& jetst(0:mxnpjet),jetvx(0:mxnpjet),jetvy(0:mxnpjet),jetvz(0:mxnpjet), &
!$acc& jetvl(0:mxnpjet),jetms(0:mxnpjet),jetch(0:mxnpjet),jetfr(0:mxnpjet))
    if(levaporation)then
!$acc enter data copyin(jetve(0:mxnpjet),jetce(0:mxnpjet))
    endif
  endif
  do s=1,4
!$acc enter data copyin(f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy, &
!$acc& f(s)%vz,f(s)%ev)
  enddo
!$acc enter data copyin(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
!$acc& zxx,zyy,zzz,zst,zvx,zvy,zvz,zev)

  workspace=.true.
  workspace_mapped=.true.
  workspace_mxnpjet=mxnpjet
  workspace_mxchunk=mxchunk
  reset_requested=.false.
  call accelerator_set_topology_enabled(.true.)
  call accelerator_set_persistent(.true.)
  call set_coulomb_accelerator_persistent(.true.)
! One asynchronous queue and one wait per step, for the topology record
! (dynamic Platen runs; the RK schemes and fixed bead sets follow in M4).
  if(scheme==scheme_platen .and. linserting .and. .not.host_oracle) &
   call accelerator_set_async(.true.)

  if(.not.engaged .or. restored)then
    if(idrank==0)then
      select case(scheme)
      case(scheme_euler)
        label='Euler'
      case(scheme_rk2)
        label='RK2'
      case(scheme_rk4)
        label='RK4'
      case default
        label='Platen'
      end select
      if(lKVfluid .and. levaporation)then
        label=trim(label)//', Kelvin-Voigt with evaporation'
      elseif(lKVfluid)then
        label=trim(label)//', Kelvin-Voigt'
      elseif(levaporation)then
        label=trim(label)//', Maxwell with evaporation'
      else
        label=trim(label)//', Maxwell'
      endif
      if(restored)then
        write(6,'(3a,i0,a,i0,a)')'OpenACC device step resumed (', &
         trim(label),') at step ',k,' with ',npjet-inpjet, &
         ' active beads, engaged before the restart'
      else
        write(6,'(3a,i0,a,i0,a)')'OpenACC device step engaged (', &
         trim(label),') at step ',k,' with ',npjet-inpjet,' active beads'
      endif
    endif
    engaged=.true.
    restored=.false.
  endif

 end subroutine ensure_device_step

 subroutine device_stage(tstage,k,s,xs,ys,zs,ss,vxs,vys,vzs,ves, &
   stochastic,prep)

!***********************************************************************
!
!     Derivatives of stage s from the stage state (xs ... vzs, and ves
!     with evaporation): smoothed nozzle charge and inserting bead,
!     Coulomb sums, model kernel, charge restored.  stochastic: the
!     deterministic part of the stochastic system (Platen: friction, no
!     stress force on a compressed element).  prep (Platen, fused step):
!     1 for the first force evaluation of the step, 2 for the later ones,
!     one kernel that restores the charge smoothed for the previous
!     evaluation, smooths it for this one and places the inserting bead;
!     the end of the step restores the last.  Absent or 0: smoothing,
!     placement and restoring as separate kernels.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: k,s
  double precision, intent(in) :: tstage
  double precision, allocatable, intent(inout) :: xs(:),ys(:),zs(:)
  double precision, allocatable, intent(inout) :: ss(:),vxs(:),vys(:),vzs(:)
  double precision, allocatable, intent(inout) :: ves(:)
  logical, intent(in), optional :: stochastic
  integer, intent(in), optional :: prep

  logical :: ok,stoc
  integer :: mode
  double precision :: friction
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
  integer :: ipoint,j,nactive
  double precision :: fstocx,fstocy,fstocz
  double precision, allocatable :: ax(:),ay(:),az(:)
#endif

  stoc=.false.
  if(present(stochastic))stoc=stochastic
  friction=0.d0
  if(stoc)friction=noisefric
  mode=0
  if(present(prep) .and. .not.host_oracle)mode=prep

#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
! Development oracle: the trusted CPU equations of motion on the
! downloaded stage, derivatives uploaded; the state updates stay on the
! device.
  nactive=myend-mystart
  if(lKVfluid)then
!$acc update self(xs(0:npjet),ys(0:npjet),zs(0:npjet),ss(0:npjet), &
!$acc& vxs(0:npjet),vys(0:npjet),vzs(0:npjet), &
!$acc& jetvl(0:npjet),jetms(0:npjet),jetch(0:npjet),jetfr(0:npjet)) if_present
    if(levaporation)then
!$acc update self(ves(0:npjet)) if_present
    endif
    allocate(ax(0:mxnpjet),ay(0:mxnpjet),az(0:mxnpjet))
    ax(:)=0.d0; ay(:)=0.d0; az(:)=0.d0
    call smooth_charge(xs,ys,zs)
    call compute_posnoinserted(xs,ys,zs)
    if(levaporation)then
      call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
    else
      call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs)
    endif
    j=0
    do ipoint=mystart,myend
      if(levaporation)then
        call eom3_KV_pos_v_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,f(s)%xx(j),f(s)%yy(j),f(s)%zz(j),ax(ipoint),ay(ipoint), &
         az(ipoint),f(s)%ev(j),tstage,k)
      else
        call xpsys_KV_pos_v(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl, &
         coulforce,f(s)%xx(j),f(s)%yy(j),f(s)%zz(j),ax(ipoint),ay(ipoint), &
         az(ipoint),tstage,k)
      endif
      j=j+1
    enddo
    j=0
    do ipoint=mystart,myend
      if(levaporation)then
        call eom3_KV_st_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,ax,ay,az,f(s)%ev(j),f(s)%st(j),tstage,k)
      else
        call xpsys_KV_st(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl, &
         coulforce,ax,ay,az,f(s)%st(j),tstage,k)
      endif
      f(s)%vx(j)=ax(ipoint); f(s)%vy(j)=ay(ipoint); f(s)%vz(j)=az(ipoint)
      j=j+1
    enddo
    call restore_charge()
    deallocate(ax,ay,az)
  elseif(levaporation)then
!$acc update self(xs(0:npjet),ys(0:npjet),zs(0:npjet),ss(0:npjet), &
!$acc& vxs(0:npjet),vys(0:npjet),vzs(0:npjet),ves(0:npjet), &
!$acc& jetvl(0:npjet),jetms(0:npjet),jetch(0:npjet),jetfr(0:npjet)) if_present
    call smooth_charge(xs,ys,zs)
    call compute_posnoinserted(xs,ys,zs)
    call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
    j=0
    do ipoint=mystart,myend
      if(stoc)then
        call xpsys_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves,coulforce, &
         f(s)%xx(j),f(s)%yy(j),f(s)%zz(j),f(s)%st(j),f(s)%vx(j),f(s)%vy(j),f(s)%vz(j), &
         f(s)%ev(j),tstage,k,fstocx,fstocy,fstocz)
      else
        call xpsys_ev_maxwell(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,f(s)%xx(j),f(s)%yy(j),f(s)%zz(j),f(s)%st(j),f(s)%vx(j),f(s)%vy(j), &
         f(s)%vz(j),f(s)%ev(j),tstage,k)
      endif
      j=j+1
    enddo
    call restore_charge()
  else
    call smooth_charge(xs,ys,zs)
    call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
     xs,ys,zs)
    call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs)
!$acc update self(xs(0:npjet),ys(0:npjet),zs(0:npjet),ss(0:npjet), &
!$acc& vxs(0:npjet),vys(0:npjet),vzs(0:npjet),jetvl(0:npjet), &
!$acc& coulforce(0:npjet,1:3)) if_present
    j=0
    do ipoint=mystart,myend
      if(stoc)then
        call xpsys(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,coulforce, &
         f(s)%xx(j),f(s)%yy(j),f(s)%zz(j),f(s)%st(j),f(s)%vx(j),f(s)%vy(j),f(s)%vz(j), &
         tstage,k,fstocx,fstocy,fstocz)
      else
        call xpsys(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,coulforce, &
         f(s)%xx(j),f(s)%yy(j),f(s)%zz(j),f(s)%st(j),f(s)%vx(j),f(s)%vy(j),f(s)%vz(j), &
         tstage,k)
      endif
      j=j+1
    enddo
    call restore_charge()
  endif
!$acc update device(f(s)%xx(0:nactive),f(s)%yy(0:nactive),f(s)%zz(0:nactive), &
!$acc& f(s)%st(0:nactive),f(s)%vx(0:nactive),f(s)%vy(0:nactive), &
!$acc& f(s)%vz(0:nactive),f(s)%ev(0:nactive))
#else
  if(mode==0)then
    call smooth_charge(xs,ys,zs)
    call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
     xs,ys,zs)
  elseif(.not.linserted)then
    call accelerator_platen_stage_prep(npjet,mode==2,thresolution, &
     dresolution,resolution,xs,ys,zs,jetch)
  endif
  if(levaporation)then
    call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
  else
    call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs)
  endif
  if(lKVfluid)then
    call accelerator_kv_stage(mystart,myend,npjet,xs,ys,zs,ss, &
     vxs,vys,vzs,jetvl,ves,coulforce,jetms,jetch,jetfr, &
     f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy,f(s)%vz, &
     f(s)%ev,linserting,linserted,liniperturb,lairdrag,lflorentz,luppot, &
     nfieldtype,pfreq,consistency,findex,yieldstress,att,fve,gr,ks,li,v, &
     velext,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev, &
     tev,evlim,evaporative=levaporation)
  elseif(levaporation)then
    call accelerator_maxwell_evap_stage(mystart,myend,npjet,xs,ys,zs,ss, &
     vxs,vys,vzs,jetvl,ves,coulforce,jetms,jetch,jetfr, &
     f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy,f(s)%vz, &
     f(s)%ev,linserting,linserted,liniperturb,lairdrag,lflorentz,luppot, &
     nfieldtype,pfreq,consistency,findex,yieldstress,att,fve,gr,ks,li,v, &
     velext,stoc,friction,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity, &
     cp0,Bev,mev,tev)
  else
! collector_curvature: the lead bead after removals has its surface
! tension and lift with the last collected bead, as in eom3 (the other
! models and the former rk4sys device path had it; Euler and RK2 not).
    ok=accelerator_eom3_stage(mystart,myend,npjet,xs,ys,zs,ss,vxs,vys, &
     vzs,jetvl,coulforce,jetms,jetch,jetfr,f(s)%xx,f(s)%yy,f(s)%zz, &
     f(s)%st,f(s)%vx,f(s)%vy,f(s)%vz,linserted,liniperturb,lairdrag, &
     lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
     att,fve,gr,ks,li,v,velext,stoc,friction,collector_curvature=.true.)
! device_step_supported admits only what the kernel evaluates; a refusal
! would leave the stage derivatives stale.
    if(.not.ok)call error(22)
  endif
  if(mode==0)call restore_charge()
#endif

 end subroutine device_stage

 subroutine stage_update(h,mode,s)

!***********************************************************************
!
!     Stage state from the step's initial state and the derivatives of
!     stage s: y = x + 0.5 h f (mode 1) or y = x + h f (mode 3).
!
!***********************************************************************

  implicit none

  double precision, intent(in) :: h
  integer, intent(in) :: mode,s

  if(levaporation)then
    call accelerator_maxwell_rk4_stage_update(mystart,myend,h,mode, &
     jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
     f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy,f(s)%vz, &
     f(s)%ev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
  else
    call accelerator_maxwell_rk4_stage_update(mystart,myend,h,mode, &
     jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetvl,jetvl, &
     f(s)%xx,f(s)%yy,f(s)%zz,f(s)%st,f(s)%vx,f(s)%vy,f(s)%vz, &
     f(s)%ev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim,evaporative=.false.)
  endif

 end subroutine stage_update

 subroutine device_rk_step(scheme,timesub,h,k)

!***********************************************************************
!
!     One Euler, RK2 or RK4 step on the device-resident state, then the
!     new state committed with the path-length and maximum-stress
!     statistics, and the inserting bead placed.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: scheme,k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h

  call ensure_device_step(scheme,k)

  select case(scheme)
  case(scheme_euler)
    call first_stage()
    call stage_update(h,3,1)
  case(scheme_rk2)
    call first_stage()
    call stage_update(h,3,1)
    call device_stage(timesub+h,k,2,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
    if(levaporation)then
      call accelerator_evap_rk2_final_update(mystart,myend,h, &
       jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
       f(1)%xx,f(1)%yy,f(1)%zz,f(1)%st,f(1)%vx,f(1)%vy,f(1)%vz, &
       f(1)%ev,f(2)%xx,f(2)%yy,f(2)%zz,f(2)%st,f(2)%vx,f(2)%vy, &
       f(2)%vz,f(2)%ev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
    else
      call accelerator_evap_rk2_final_update(mystart,myend,h, &
       jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetvl,jetvl, &
       f(1)%xx,f(1)%yy,f(1)%zz,f(1)%st,f(1)%vx,f(1)%vy,f(1)%vz, &
       f(1)%ev,f(2)%xx,f(2)%yy,f(2)%zz,f(2)%st,f(2)%vx,f(2)%vy, &
       f(2)%vz,f(2)%ev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim, &
       evaporative=.false.)
    endif
  case default
    call first_stage()
    call stage_update(h,1,1)
    call device_stage(timesub+0.5d0*h,k,2,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
    call stage_update(h,1,2)
    call device_stage(timesub+0.5d0*h,k,3,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
    call stage_update(h,3,3)
    call device_stage(timesub+h,k,4,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
    if(levaporation)then
      call accelerator_maxwell_rk4_final_update(mystart,myend,h, &
       jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
       f(1)%xx,f(1)%yy,f(1)%zz,f(1)%st,f(1)%vx,f(1)%vy,f(1)%vz, &
       f(1)%ev,f(2)%xx,f(2)%yy,f(2)%zz,f(2)%st,f(2)%vx,f(2)%vy, &
       f(2)%vz,f(2)%ev,f(3)%xx,f(3)%yy,f(3)%zz,f(3)%st,f(3)%vx, &
       f(3)%vy,f(3)%vz,f(3)%ev,f(4)%xx,f(4)%yy,f(4)%zz,f(4)%st, &
       f(4)%vx,f(4)%vy,f(4)%vz,f(4)%ev, &
       yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
    else
      call accelerator_maxwell_rk4_final_update(mystart,myend,h, &
       jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetvl,jetvl, &
       f(1)%xx,f(1)%yy,f(1)%zz,f(1)%st,f(1)%vx,f(1)%vy,f(1)%vz, &
       f(1)%ev,f(2)%xx,f(2)%yy,f(2)%zz,f(2)%st,f(2)%vx,f(2)%vy, &
       f(2)%vz,f(2)%ev,f(3)%xx,f(3)%yy,f(3)%zz,f(3)%st,f(3)%vx, &
       f(3)%vy,f(3)%vz,f(3)%ev,f(4)%xx,f(4)%yy,f(4)%zz,f(4)%st, &
       f(4)%vx,f(4)%vy,f(4)%vz,f(4)%ev, &
       yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim,evaporative=.false.)
    endif
  end select

  if(levaporation)then
    call accelerator_maxwell_commit_state(mystart,myend, &
     yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
     jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve, &
     counterlpath,ncounterlpath,maxstress,maxstressposx)
  else
    call accelerator_maxwell_commit_state(mystart,myend, &
     yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
     jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,yev, &
     counterlpath,ncounterlpath,maxstress,maxstressposx, &
     evaporative=.false.)
  endif
  call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
   jetxx,jetyy,jetzz)
! Beads reaching the collector are frozen there (remove_jetbead); with
! insertion the topology step does it.
  if(.not.linserting)call accelerator_freeze_at_collector(inpjet,npjet, &
   collector_h,jetxx,jetfr)
  call accelerator_mark_device_state(.true.)
  timesub=timesub+h

 contains

  subroutine first_stage()
! The stage at the step's initial state; without evaporation the unused
! yev stands for the evaporated volume, which is not referenced.
   if(levaporation)then
     call device_stage(timesub,k,1,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
      jetvz,jetve)
   else
     call device_stage(timesub,k,1,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
      jetvz,yev)
   endif
  end subroutine first_stage

 end subroutine device_rk_step

 subroutine device_platen_step(timesub,h,k)

!***********************************************************************
!
!     One Platen step (Heun for the deterministic part, order 1.5 for the
!     stochastic velocity) on the device-resident state: three force
!     evaluations, at the step's initial state and at the two predicted
!     states y and z, then the per-bead update (velocity with the pool
!     noise of the slice platen/platen_ev reserved for this step,
!     evaporation rate, positions) and the single-gang end of the step
!     (final stress derivative and stress, statistics, inserting bead,
!     topology decisions).  Without evaporation jetvl stands for the
!     evaporated volume, which is not referenced.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h

  call ensure_device_step(scheme_platen,k)

  if(levaporation)then
    call device_stage(timesub,k,1,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
     jetvz,jetve,stochastic=.true.,prep=1)
    call accelerator_platen_predict(mystart,myend,h,airdragamp(1), &
     noisediff,evlim,jetms,jetvl,jetve,jetxx,jetyy,jetzz,jetst,jetvx, &
     jetvy,jetvz,f(1)%xx,f(1)%yy,f(1)%zz,f(1)%st,f(1)%vx,f(1)%vy, &
     f(1)%vz,f(1)%ev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
     zxx,zyy,zzz,zst,zvx,zvy,zvz,zev,linserted,jetfr,.true.)
  else
    call device_stage(timesub,k,1,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
     jetvz,yev,stochastic=.true.,prep=1)
    call accelerator_platen_predict(mystart,myend,h,airdragamp(1), &
     noisediff,evlim,jetms,jetvl,jetvl,jetxx,jetyy,jetzz,jetst,jetvx, &
     jetvy,jetvz,f(1)%xx,f(1)%yy,f(1)%zz,f(1)%st,f(1)%vx,f(1)%vy, &
     f(1)%vz,f(1)%ev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
     zxx,zyy,zzz,zst,zvx,zvy,zvz,zev,linserted,jetfr,.false.)
  endif
  call device_stage(timesub,k,2,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   stochastic=.true.,prep=2)
  call device_stage(timesub,k,3,zxx,zyy,zzz,zst,zvx,zvy,zvz,zev, &
   stochastic=.true.,prep=2)

  if(levaporation)then
    call accelerator_platen_update(mystart,myend,h,npjet,mxnpjet, &
     linserted,.true.,airdragamp(1),noisediff,pfreq,liniperturb,evlim, &
     gaussianhistorybase,gaussianhistorywindow,gaussianhistoryvalues, &
     gaussianhistory,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetms,jetvl, &
     jetve,jetfr,f(1)%xx,f(1)%yy,f(1)%zz,f(1)%vx,f(1)%vy,f(1)%vz, &
     f(1)%ev,f(2)%vx,f(2)%vy,f(2)%vz,f(3)%vx,f(3)%vy,f(3)%vz, &
     yxx,yyy,yzz,yev,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity)
    call accelerator_platen_end_step(mystart,myend,h,npjet,mxnpjet, &
     inpjet,linserted,lremove,collector_h,resolution,dresolution, &
     thresolution,ivelocity,istress,imassa,icharge,ivolume,jetxx,jetyy, &
     jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetve,jetfr, &
     f(1)%st,yst,.true.,consistency,findex,yieldstress,cp0,Bev,mev,tev, &
     counterlpath,ncounterlpath,maxstress,maxstressposx,linserting, &
     .not.host_oracle)
  else
    call accelerator_platen_update(mystart,myend,h,npjet,mxnpjet, &
     linserted,.false.,airdragamp(1),noisediff,pfreq,liniperturb,evlim, &
     gaussianhistorybase,gaussianhistorywindow,gaussianhistoryvalues, &
     gaussianhistory,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetms,jetvl, &
     jetvl,jetfr,f(1)%xx,f(1)%yy,f(1)%zz,f(1)%vx,f(1)%vy,f(1)%vz, &
     f(1)%ev,f(2)%vx,f(2)%vy,f(2)%vz,f(3)%vx,f(3)%vy,f(3)%vz, &
     yxx,yyy,yzz,yev,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity)
    call accelerator_platen_end_step(mystart,myend,h,npjet,mxnpjet, &
     inpjet,linserted,lremove,collector_h,resolution,dresolution, &
     thresolution,ivelocity,istress,imassa,icharge,ivolume,jetxx,jetyy, &
     jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetvl,jetfr, &
     f(1)%st,yst,.false.,consistency,findex,yieldstress,cp0,Bev,mev,tev, &
     counterlpath,ncounterlpath,maxstress,maxstressposx,linserting, &
     .not.host_oracle)
  endif
! Beads reaching the collector are frozen there (remove_jetbead); with
! insertion the end of the step does it.
  if(.not.linserting)call accelerator_freeze_at_collector(inpjet,npjet, &
   collector_h,jetxx,jetfr)
  call accelerator_mark_device_state(.true.)
  timesub=timesub+h

 end subroutine device_platen_step
#endif

 end module device_step_mod
