
module integrator_mod

 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

!***********************************************************************
!     
!     JETSPIN module containing integrators data routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2015
!     
!***********************************************************************

 use version_mod,       only : mystart,myend,mxchunk,sum_world_darr, &
                         set_chunk,set_mxchunk,idrank
 use error_mod,         only : error,warning
 use utility_mod,       only : prepare_gaussian_buffer, &
                         gaussian_buffer_value,prepare_gaussian_history, &
                         gaussian_history_value,gaussianhistory, &
                         mark_gaussianhistory_device_mapped, &
                         begin_gaussian_history_step,gaussianhistoryvalues
 use nanojet_mod,       only : doallocate,mxnpjet,npjet,inpjet,systype, &
                         jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetst,jetvl, &
                         compute_posnoinserted,lKVfluid,levaporation,jetve, &
                         evlim,linserting,lengthscale
 use dynamic_refinement_mod, only : driver_dynamic_refinement
 use profiling_mod, only : profiling_start,profiling_stop,prof_eom, &
                         prof_rk_update
#ifdef _OPENACC
 use accelerator_mod, only : accelerator_mark_device_state, &
                         accelerator_device_state_is_current, &
                         accelerator_is_persistent
#endif
 use coulomb_force_mod, only : smooth_charge,restore_charge,coulforce, &
                         compute_coulomelec_driver
 use driver_eom_mod,    only : xpsys,xpsys_pos,xpsys_stress,xpsys_KV_pos_v, &
                         xpsys_KV_st,xpsys_ev,xpsys_ev_maxwell,xpsys_pos_ev, &
                         xpsys_stress_ev
 use device_step_mod,   only : device_step_eligible, &
                         device_step_request_reset,device_step_supported, &
                         scheme_euler,scheme_rk2,scheme_rk4,scheme_platen
#ifdef _OPENACC
 use device_step_mod,   only : device_rk_step,device_platen_step
#endif

 implicit none

 private
 
 integer, public, save :: integrator
 logical, public, save :: lintegrator=.false.
 
 double precision, public, save :: initime = 0.d0
 double precision, public, save :: endtime = 5.d0
 logical, public, save :: lendtime

 
 public :: driver_integrator
 public :: prepare_integrator_random_history
 public :: reset_persistent_integrator
 public :: integration_last_step

contains

 integer function integration_last_step(h)
! The step that ends the main loop: the first whose time dble(nstep)*h
! reaches endtime, to within a relative 1e-12.  Both are divided by the
! time unit when the input is read, so a final time that is a multiple of
! the timestep gives dble(n)*h one ulp below endtime when the division is
! rounded exactly (GFortran, NVFORTRAN 25.5); until 2026-10-07 such runs
! made one step more than NVFORTRAN 24.3 builds, whose -O3 code divides by
! multiplying with the reciprocal.
  implicit none
  double precision, intent(in) :: h
  integration_last_step=max(1,ceiling((endtime/h)*(1.d0-1.d-12)))
 end function integration_last_step

 subroutine reset_persistent_integrator()
! A capacity change: the device step rebuilds its workspace.
  implicit none
  call device_step_request_reset()
 end subroutine reset_persistent_integrator

 subroutine prepare_integrator_random_history(h)
  implicit none
  double precision, intent(in) :: h
  integer :: nsteps,history_last
  logical :: dynamic_history
! Decided once, before the timestep loop, while a single-bead start still
! has npjet=1: the model conditions of the device step
! (device_step_supported), whatever the build and the number of ranks, so
! that CPU, GPU, serial and MPI runs consume the same noise.  The per-step
! gate (device_step_eligible) opens the device step once the jet has
! grown, in a serial OpenACC run.  Until 2026-09-30 the evaporative
! dynamic run tested its full per-step gate here, which a single-bead start
! never passed; until 2026-10-05 only the evaporative dynamic run, the
! fixed 1000-bead geometries and a separate build option drew the pool;
! until 2026-10-06 the dynamic runs needed refinement and tagged beads and
! the fixed geometries exactly 1000 beads (Example 4 drew its noise step
! by step until then).
  if(integrator/=4 .or. .not.device_step_supported(scheme_platen))return
  dynamic_history=linserting
! The steps of the main loop (integration_last_step).  Until 2026-10-06
! this was nint((endtime-initime)/h): with an endtime that is not a
! multiple of h the loop ran one step more, which read the pool of a fixed
! jet from its start again.
  nsteps=integration_last_step(h)
  history_last=npjet
! A dynamic run takes the whole pool, since its bead count can grow beyond
! the initial capacity.  Beads inserted after initialization then consume
! the same sequential pool on CPU and GPU without drawing random numbers
! inside the timestep loop.
  if(dynamic_history)history_last=mxnpjet
  call prepare_gaussian_history(inpjet,history_last,mxnpjet,3,nsteps, &
   full_pool=dynamic_history)
#ifdef _OPENACC
! The pool is a flat sequence whose size is decided once, above, and never
! changes afterwards.
!$acc enter data copyin(gaussianhistory(0:gaussianhistoryvalues-1))
  call mark_gaussianhistory_device_mapped(.true.)
#endif
 end subroutine prepare_integrator_random_history

  
 subroutine driver_integrator(timesub,h,k,dorefinment)
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which integrate the system 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2016
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(inout) :: dorefinment
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
  integer :: i
  logical :: ltestinst,lthisbeadinst
  logical :: lreportedinst
  integer :: nnaninst,firstnaninst,lastnaninst
  
! perform the dynamic refinement if is requested
  call driver_dynamic_refinement(k,dorefinment)
    
  call set_chunk(inpjet,npjet)
  call set_mxchunk(mxnpjet)
  if(lKVfluid)then
    select case(integrator)
      case(1)
        call eulsys_KV(timesub,h,k)
      case(2)
        call rk2sys_KV(timesub,h,k)
      case(3)
        call rk4sys_KV(timesub,h,k)
      case default
        call error(1)
    end select 
  else
    if(levaporation)then
      select case(integrator)
        case(1)
          call eulsys_ev(timesub,h,k)
        case(2)
          call rk2sys_ev(timesub,h,k)
        case(3)
          call rk4sys_ev(timesub,h,k)
        case(4)
          call platen_ev(timesub,h,k)
        case default
          call error(1)
      end select 
    else
      select case(integrator)
        case(1)
          call eulsys(timesub,h,k)
        case(2)
          call rk2sys(timesub,h,k)
        case(3)
          call rk4sys(timesub,h,k)
        case(4)
          call platen(timesub,h,k)
        case default
          call error(1)
      end select 
    endif
  endif

#ifdef _OPENACC
! The persistent accelerator path leaves the integrated state on the device.
! Host/debug fallbacks keep the host arrays authoritative instead.
  call accelerator_mark_device_state(accelerator_is_persistent())
#endif
  
  ltestinst=.false.
  lreportedinst=.false.
  nnaninst=0
  firstnaninst=-1
  lastnaninst=-1
#ifdef _OPENACC
  if(.not.accelerator_device_state_is_current())then
#endif
    do i=inpjet,npjet
      lthisbeadinst=.false.
      if(ieee_is_nan(dcos(jetxx(i))))lthisbeadinst=.true.
      if(ieee_is_nan(dcos(jetyy(i))))lthisbeadinst=.true.
      if(ieee_is_nan(dcos(jetzz(i))))lthisbeadinst=.true.
      if(ieee_is_nan(dcos(jetst(i))))lthisbeadinst=.true.
      if(ieee_is_nan(dcos(jetvx(i))))lthisbeadinst=.true.
      if(ieee_is_nan(dcos(jetvy(i))))lthisbeadinst=.true.
      if(ieee_is_nan(dcos(jetvz(i))))lthisbeadinst=.true.
      if(lthisbeadinst)ltestinst=.true.
! Report the first offending bead in full detail and count how many beads
! are affected in total. The direct all-to-all Coulomb coupling can spread
! a single bad value to every bead within one timestep, so the count
! distinguishes a widespread propagated failure from an isolated one.
      if(lthisbeadinst)then
        nnaninst=nnaninst+1
        if(firstnaninst==-1)firstnaninst=i
        lastnaninst=i
        if(.not.lreportedinst .and. idrank==0)then
          lreportedinst=.true.
          write(6,'(a,i0,7(a,l1))')'Numerical instability detail: bead=',i, &
           ' x_nan=',ieee_is_nan(dcos(jetxx(i))), &
           ' y_nan=',ieee_is_nan(dcos(jetyy(i))), &
           ' z_nan=',ieee_is_nan(dcos(jetzz(i))), &
           ' st_nan=',ieee_is_nan(dcos(jetst(i))), &
           ' vx_nan=',ieee_is_nan(dcos(jetvx(i))), &
           ' vy_nan=',ieee_is_nan(dcos(jetvy(i))), &
           ' vz_nan=',ieee_is_nan(dcos(jetvz(i)))
! Development-only diagnostic: locate the offending bead relative to the
! nozzle (high index, npjet) and collector (low index, inpjet) ends, using
! the last still-finite neighbour bead since the flagged bead's own
! coordinates/stress may already be NaN.
! Two tests: Fortran does not short-circuit .and., so a single one read
! jetxx(inpjet-1) when the first bead failed (until 2026-10-07).
          if(i-1>=inpjet)then
            if(.not.ieee_is_nan(dcos(jetxx(i-1))))then
              write(6,'(a,i0,a,es14.6,a,es14.6)') &
               'Numerical instability last-good neighbor: bead=',i-1, &
               ' x_cm=',jetxx(i-1)*lengthscale,' reference_volume_cm3=',jetvl(i-1)
            endif
          endif
          write(6,'(a,i0,a,es14.6,a,i0,a,es14.6)') &
           'Numerical instability domain extent: collector_side_bead=',inpjet, &
           ' x_cm=',jetxx(inpjet)*lengthscale,' nozzle_side_bead=',npjet, &
           ' x_cm=',jetxx(npjet)*lengthscale
        endif
      endif
    enddo
#ifdef _OPENACC
  endif
#endif

  if(ltestinst)then
    if(idrank==0)write(6,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
     'Numerical instability context: nstep=',k,' inpjet=',inpjet, &
     ' npjet=',npjet,' nan_bead_count=',nnaninst, &
     ' first_nan_bead=',firstnaninst,' last_nan_bead=',lastnaninst
    call warning(67,dble(k))
    call error(14)
  endif
  
  doallocate=.false.
  
  return
  
 end subroutine
 
 subroutine eulsys(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     first order accurate Euler scheme
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2015
!     
!***********************************************************************
  
  implicit none
  
  
  
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
! service arrays
  double precision, allocatable, dimension (:), save ::  fxx
  double precision, allocatable, dimension (:), save ::  fyy
  double precision, allocatable, dimension (:), save ::  fzz
  double precision, allocatable, dimension (:), save ::  fst
  double precision, allocatable, dimension (:), save ::  fvx
  double precision, allocatable, dimension (:), save ::  fvy
  double precision, allocatable, dimension (:), save ::  fvz
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  
  integer :: ipoint,j
  
  logical, save :: lfirstsub=.true.
  
  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_euler))then
    call device_rk_step(scheme_euler,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(fxx)
        deallocate(fst)
        deallocate(fvx)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
      endif
      allocate(fxx(0:mxchunk))
      allocate(fst(0:mxchunk))
      allocate(fvx(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(fxx)
        deallocate(fyy)
        deallocate(fzz)
        deallocate(fst)
        deallocate(fvx)
        deallocate(fvy)
        deallocate(fvz)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
      endif
      allocate(fxx(0:mxchunk))
      allocate(fyy(0:mxchunk))
      allocate(fzz(0:mxchunk))
      allocate(fst(0:mxchunk))
      allocate(fvx(0:mxchunk))
      allocate(fvy(0:mxchunk))
      allocate(fvz(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif

  
  
! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,coulforce,fxx(j),fyy(j),fzz(j),fst(j), &
          fvx(j),fvy(j),fvz(j),timesub,k)
        j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
	    yxx(ipoint) = jetxx(ipoint) + h*fxx(j)
	    yst(ipoint) = jetst(ipoint) + h*fst(j)
	    yvx(ipoint) = jetvx(ipoint) + h*fvx(j)
	    j=j+1
      enddo
      
      call restore_charge()
      
      timesub=timesub+h
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz,timesub)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,coulforce,fxx(j),fyy(j),fzz(j),fst(j), &
          fvx(j),fvy(j),fvz(j),timesub,k)
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      call restore_charge()
      j=0
      call profiling_start(prof_rk_update)
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
	      yxx(ipoint) = jetxx(ipoint) + h*fxx(j)
	      yyy(ipoint) = jetyy(ipoint) + h*fyy(j)
	      yzz(ipoint) = jetzz(ipoint) + h*fzz(j)
	      yst(ipoint) = jetst(ipoint) + h*fst(j)
	      yvx(ipoint) = jetvx(ipoint) + h*fvx(j)
	      yvy(ipoint) = jetvy(ipoint) + h*fvy(j)
	      yvz(ipoint) = jetvz(ipoint) + h*fvz(j)
	      j=j+1
      enddo
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      call profiling_stop(prof_rk_update)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select
  
  
  return
  
  
 end subroutine eulsys 
 
 subroutine rk2sys(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     second order accurate Heun scheme
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  
  
  integer,intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision,intent(in) :: h
  
  integer :: ipoint,j
! service arrays
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz

  logical, save :: lfirstsub=.true.
  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_rk2))then
    call device_rk_step(scheme_rk2,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1yy(0:mxchunk))
      allocate(f1zz(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f1vy(0:mxchunk))
      allocate(f1vz(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2yy(0:mxchunk))
      allocate(f2zz(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f2vy(0:mxchunk))
      allocate(f2vz(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif

  
  
! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k)
        f1xx(j)=fxx
        f1st(j)=fst
        f1vx(j)=fvx
	    yxx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    yst(ipoint) = jetst(ipoint) + h*f1st(j)
	    yvx(ipoint) = jetvx(ipoint) + h*f1vx(j)
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),fyy,fzz,f2st(j), &
         f2vx(j),fvy,fvz,timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
	    yst(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
	    yvx(ipoint) = jetvx(ipoint) + (h/2.d0)*(f1vx(j)+f2vx(j))
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      timesub=timesub+h
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k)
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      j=0
      call profiling_start(prof_rk_update)
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint)=jetxx(ipoint)+h*f1xx(j)
        yyy(ipoint)=jetyy(ipoint)+h*f1yy(j)
        yzz(ipoint)=jetzz(ipoint)+h*f1zz(j)
        yst(ipoint)=jetst(ipoint)+h*f1st(j)
        yvx(ipoint)=jetvx(ipoint)+h*f1vx(j)
        yvy(ipoint)=jetvy(ipoint)+h*f1vy(j)
        yvz(ipoint)=jetvz(ipoint)+h*f1vz(j)
        j=j+1
      enddo
      call profiling_stop(prof_rk_update)
      call restore_charge()
      
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
      
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),timesub+h,k)
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      j=0
      call profiling_start(prof_rk_update)
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint)=jetxx(ipoint)+(h/2.d0)*(f1xx(j)+f2xx(j))
        yyy(ipoint)=jetyy(ipoint)+(h/2.d0)*(f1yy(j)+f2yy(j))
        yzz(ipoint)=jetzz(ipoint)+(h/2.d0)*(f1zz(j)+f2zz(j))
        yst(ipoint)=jetst(ipoint)+(h/2.d0)*(f1st(j)+f2st(j))
        yvx(ipoint)=jetvx(ipoint)+(h/2.d0)*(f1vx(j)+f2vx(j))
        yvy(ipoint)=jetvy(ipoint)+(h/2.d0)*(f1vy(j)+f2vy(j))
        yvz(ipoint)=jetvz(ipoint)+(h/2.d0)*(f1vz(j)+f2vz(j))
        j=j+1
      enddo
      call profiling_stop(prof_rk_update)
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select

  return
      
 end subroutine rk2sys
 
 subroutine rk4sys(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     fourth order accurate Runge-Kutta scheme
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
! service arrays
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  f3xx
  double precision, allocatable, dimension (:), save ::  f3yy
  double precision, allocatable, dimension (:), save ::  f3zz
  double precision, allocatable, dimension (:), save ::  f3st
  double precision, allocatable, dimension (:), save ::  f3vx
  double precision, allocatable, dimension (:), save ::  f3vy
  double precision, allocatable, dimension (:), save ::  f3vz
  double precision, allocatable, dimension (:), save ::  f4xx
  double precision, allocatable, dimension (:), save ::  f4yy
  double precision, allocatable, dimension (:), save ::  f4zz
  double precision, allocatable, dimension (:), save ::  f4st
  double precision, allocatable, dimension (:), save ::  f4vx
  double precision, allocatable, dimension (:), save ::  f4vy
  double precision, allocatable, dimension (:), save ::  f4vz
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  integer :: ipoint,dm,nv,j
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
  logical, save :: lfirstsub=.true.
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz

  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_rk4))then
    call device_rk_step(scheme_rk4,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f3xx)
        deallocate(f3st)
        deallocate(f3vx)
        deallocate(f4xx)
        deallocate(f4st)
        deallocate(f4vx)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f3xx(0:mxchunk))
      allocate(f3st(0:mxchunk))
      allocate(f3vx(0:mxchunk))
      allocate(f4xx(0:mxchunk))
      allocate(f4st(0:mxchunk))
      allocate(f4vx(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(f3xx)
        deallocate(f3yy)
        deallocate(f3zz)
        deallocate(f3st)
        deallocate(f3vx)
        deallocate(f3vy)
        deallocate(f3vz)
        deallocate(f4xx)
        deallocate(f4yy)
        deallocate(f4zz)
        deallocate(f4st)
        deallocate(f4vx)
        deallocate(f4vy)
        deallocate(f4vz)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1yy(0:mxchunk))
      allocate(f1zz(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f1vy(0:mxchunk))
      allocate(f1vz(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2yy(0:mxchunk))
      allocate(f2zz(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f2vy(0:mxchunk))
      allocate(f2vz(0:mxchunk))
      allocate(f3xx(0:mxchunk))
      allocate(f3yy(0:mxchunk))
      allocate(f3zz(0:mxchunk))
      allocate(f3st(0:mxchunk))
      allocate(f3vx(0:mxchunk))
      allocate(f3vy(0:mxchunk))
      allocate(f3vz(0:mxchunk))
      allocate(f4xx(0:mxchunk))
      allocate(f4yy(0:mxchunk))
      allocate(f4zz(0:mxchunk))
      allocate(f4st(0:mxchunk))
      allocate(f4vx(0:mxchunk))
      allocate(f4vy(0:mxchunk))
      allocate(f4vz(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif

  
! select the proper system type
  select case(systype)
    case(1)
!     1°step
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k)
        f1xx(j)=fxx
        f1st(j)=fst
        f1vx(j)=fvx
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     2°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),fyy,fzz,f2st(j), &
         f2vx(j),fvy,fvz,timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(j)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     3°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f3xx(j),fyy,fzz,f3st(j), &
         f3vx(j),fvy,fvz,timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
        yst(ipoint) = jetst(ipoint) + h*f3st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f3vx(j)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     4°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f4xx(j),fyy,fzz,f4st(j), &
         f4vx(j),fvy,fvz,timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/6.d0)*(f1xx(j)+ &
         2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
        yst(ipoint) = jetst(ipoint) + (h/6.d0)*(f1st(j)+ &
         2.d0*(f2st(j)+f3st(j))+f4st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(j)+ &
         2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      timesub=timesub+h
  case default
!     1°step
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl, &
       jetxx,jetyy,jetzz)
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         jetvl,coulforce,f1xx(j),f1yy(j),f1zz(j),f1st(j),f1vx(j), &
         f1vy(j),f1vz(j),timesub,k)
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      j=0
      call profiling_start(prof_rk_update)
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f1yy(j)
        yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f1zz(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f1vy(j)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f1vz(j)
        j=j+1
      enddo
      call profiling_stop(prof_rk_update)
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     2°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),timesub+h/2.d0,k)
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      call profiling_start(prof_rk_update)
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
        yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f2yy(j)
        yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f2zz(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(j)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f2vy(j)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f2vz(j)
        j=j+1
      enddo
      call profiling_stop(prof_rk_update)
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     3°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
         f3vx(j),f3vy(j),f3vz(j),timesub+h/2.d0,k)
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      call profiling_start(prof_rk_update)
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
        yyy(ipoint) = jetyy(ipoint) + h*f3yy(j)
        yzz(ipoint) = jetzz(ipoint) + h*f3zz(j)
        yst(ipoint) = jetst(ipoint) + h*f3st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f3vx(j)
        yvy(ipoint) = jetvy(ipoint) + h*f3vy(j)
        yvz(ipoint) = jetvz(ipoint) + h*f3vz(j)
        j=j+1
      enddo
      call profiling_stop(prof_rk_update)
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     4°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
         f4vx(j),f4vy(j),f4vz(j),timesub+h,k)
        j=j+1
      enddo
      call profiling_stop(prof_eom)
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      call profiling_start(prof_rk_update)
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/6.d0)*(f1xx(j)+ &
         2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
        yyy(ipoint) = jetyy(ipoint) + (h/6.d0)*(f1yy(j)+ &
         2.d0*(f2yy(j)+f3yy(j))+f4yy(j))
        yzz(ipoint) = jetzz(ipoint) + (h/6.d0)*(f1zz(j)+ &
         2.d0*(f2zz(j)+f3zz(j))+f4zz(j))
        yst(ipoint) = jetst(ipoint) + (h/6.d0)*(f1st(j)+ &
         2.d0*(f2st(j)+f3st(j))+f4st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(j)+ &
         2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
        yvy(ipoint) = jetvy(ipoint) + (h/6.d0)*(f1vy(j)+ &
         2.d0*(f2vy(j)+f3vy(j))+f4vy(j))
        yvz(ipoint) = jetvz(ipoint) + (h/6.d0)*(f1vz(j)+ &
         2.d0*(f2vz(j)+f3vz(j))+f4vz(j))
        j=j+1
      enddo
      call profiling_stop(prof_rk_update)
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select
  
  return
  
 end subroutine rk4sys
 
 
 subroutine platen(timesub,h,k)
 
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the stochastic equation
!     of motion by the 1.5° order accurate Platen scheme
!     for the velocity stochastic part and by the second order accurate
!     Heun scheme for deterministic part
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f1stocvx
  double precision, allocatable, dimension (:), save ::  f1stocvy
  double precision, allocatable, dimension (:), save ::  f1stocvz
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  y1xx
  double precision, allocatable, dimension (:), save ::  y1yy
  double precision, allocatable, dimension (:), save ::  y1zz
  double precision, allocatable, dimension (:), save ::  y1st
  double precision, allocatable, dimension (:), save ::  y1vx
  double precision, allocatable, dimension (:), save ::  y1vy
  double precision, allocatable, dimension (:), save ::  y1vz
  double precision, allocatable, dimension (:), save ::  y2xx
  double precision, allocatable, dimension (:), save ::  y2yy
  double precision, allocatable, dimension (:), save ::  y2zz
  double precision, allocatable, dimension (:), save ::  y2st
  double precision, allocatable, dimension (:), save ::  y2vx
  double precision, allocatable, dimension (:), save ::  y2vy
  double precision, allocatable, dimension (:), save ::  y2vz
  double precision, allocatable, dimension (:), save ::  d3xx,d3yy,d3zz
  double precision, allocatable, dimension (:), save ::  d3st,d3vx,d3vy,d3vz
  
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
  integer :: ipoint,dm,nv,j
  double precision :: dsqrh,tsqh,prefactor1,zztang,u1,u2
  double precision, dimension(1:3) :: ww,zz,utang
  
  logical, save :: lfirstsub=.true.
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  double precision ::  fstocvx
  double precision ::  fstocvy
  double precision ::  fstocvz
  
  double precision ::  f3xx
  double precision ::  f3yy
  double precision ::  f3zz
  double precision ::  f3st
  double precision ::  f3vx
  double precision ::  f3vy
  double precision ::  f3vz
  double precision ::  f3stocvx
  double precision ::  f3stocvy
  double precision ::  f3stocvz


! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1stocvx)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(y1xx)
        deallocate(y1st)
        deallocate(y1vx)
        deallocate(y2xx)
        deallocate(y2st)
        deallocate(y2vx)
      endif
      allocate(f1xx(0:mxnpjet))
      allocate(f1st(0:mxnpjet))
      allocate(f1vx(0:mxnpjet))
      allocate(f1stocvx(0:mxnpjet))
      allocate(f2xx(0:mxnpjet))
      allocate(f2st(0:mxnpjet))
      allocate(f2vx(0:mxnpjet))
      allocate(y1xx(0:mxnpjet))
      allocate(y1st(0:mxnpjet))
      allocate(y1vx(0:mxnpjet))
      allocate(y2xx(0:mxnpjet))
      allocate(y2st(0:mxnpjet))
      allocate(y2vx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f1stocvx)
        deallocate(f1stocvy)
        deallocate(f1stocvz)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(y1xx)
        deallocate(y1yy)
        deallocate(y1zz)
        deallocate(y1st)
        deallocate(y1vx)
        deallocate(y1vy)
        deallocate(y1vz)
        deallocate(y2xx)
        deallocate(y2yy)
        deallocate(y2zz)
        deallocate(y2st)
        deallocate(y2vx)
        deallocate(y2vy)
        deallocate(y2vz)
        deallocate(d3xx,d3yy,d3zz,d3st,d3vx,d3vy,d3vz)
      endif
      allocate(f1xx(0:mxnpjet))
      allocate(f1yy(0:mxnpjet))
      allocate(f1zz(0:mxnpjet))
      allocate(f1st(0:mxnpjet))
      allocate(f1vx(0:mxnpjet))
      allocate(f1vy(0:mxnpjet))
      allocate(f1vz(0:mxnpjet))
      allocate(f1stocvx(0:mxnpjet))
      allocate(f1stocvy(0:mxnpjet))
      allocate(f1stocvz(0:mxnpjet))
      allocate(f2xx(0:mxnpjet))
      allocate(f2yy(0:mxnpjet))
      allocate(f2zz(0:mxnpjet))
      allocate(f2st(0:mxnpjet))
      allocate(f2vx(0:mxnpjet))
      allocate(f2vy(0:mxnpjet))
      allocate(f2vz(0:mxnpjet))
      allocate(y1xx(0:mxnpjet))
      allocate(y1yy(0:mxnpjet))
      allocate(y1zz(0:mxnpjet))
      allocate(y1st(0:mxnpjet))
      allocate(y1vx(0:mxnpjet))
      allocate(y1vy(0:mxnpjet))
      allocate(y1vz(0:mxnpjet))
      allocate(y2xx(0:mxnpjet))
      allocate(y2yy(0:mxnpjet))
      allocate(y2zz(0:mxnpjet))
      allocate(y2st(0:mxnpjet))
      allocate(y2vx(0:mxnpjet))
      allocate(y2vy(0:mxnpjet))
      allocate(y2vz(0:mxnpjet))
      allocate(d3xx(0:mxnpjet),d3yy(0:mxnpjet),d3zz(0:mxnpjet))
      allocate(d3st(0:mxnpjet),d3vx(0:mxnpjet),d3vy(0:mxnpjet),d3vz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif


  dsqrh=dsqrt(dabs(h))
  tsqh=dsqrh**3.d0
  prefactor1=0.5d0/dsqrh

  if(.not.allocated(gaussianhistory))then
    if(systype==1)then
      call prepare_gaussian_buffer(inpjet,npjet,mxnpjet,1)
    else
      call prepare_gaussian_buffer(inpjet,npjet,mxnpjet,3)
    endif
  endif

! Reserve this timestep's slice of the sequential Gaussian pool. Must run
! before any branch below reads the history, including the device step
! that returns early, so it sits ahead of the dispatch. A no-op unless the
! pool is in use.  The slice spans the whole jet, not this rank's beads,
! so that every MPI rank reads the values a serial run reads (mystart..myend
! until 2026-10-06, when MPI runs did not use the pool).
  call begin_gaussian_history_step(inpjet,npjet)

#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod), which
! reads the pool slice just reserved; below it, and for the options that
! step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_platen))then
    call device_platen_step(timesub,h,k)
    return
  endif
#endif

! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      y1xx(:)=0.d0
      y1st(:)=0.d0
      y1vx(:)=0.d0
      y2xx(:)=0.d0
      y2st(:)=0.d0
      y2vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         jetvl,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k, &
         fstocvx,fstocvy,fstocvz)
        f1xx(j)=fxx
        f1st(j)=fst
        f1vx(j)=fvx
        f1stocvx(j)=fstocvx
	    y1xx(ipoint) = jetxx(ipoint) + h*fxx
	    y1st(ipoint) = jetst(ipoint) + h*fst
	    y1vx(ipoint) = jetvx(ipoint) + h*fvx + dsqrh*fstocvx
	    y2xx(ipoint) = jetxx(ipoint) + h*fxx
	    y2st(ipoint) = jetst(ipoint) + h*fst
	    y2vx(ipoint) = jetvx(ipoint) + h*fvx - dsqrh*fstocvx
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      call sum_world_darr(y1vx,npjet+1)
      call sum_world_darr(y2xx,npjet+1)
      call sum_world_darr(y2st,npjet+1)
      call sum_world_darr(y2vx,npjet+1)
      
      call smooth_charge(y1xx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y1xx, &
       y1yy,y1zz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,jetvl, &
         coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k,fstocvx, &
         fstocvy,fstocvz)
        f2xx(j)=fxx
        f2st(j)=fst
        f2vx(j)=fvx
	    j=j+1
      enddo
      call restore_charge()
      
      call smooth_charge(y2xx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y2xx, &
       y2yy,y2zz)
      j=0
      y1vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys(ipoint,y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,jetvl, &
         coulforce,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,timesub,k, &
         f3stocvx,f3stocvy,f3stocvz)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,1,1)
          u2=gaussian_history_value(k,ipoint,1,2)
        else
          u1=gaussian_buffer_value(ipoint,1,1)
          u2=gaussian_buffer_value(ipoint,1,2)
        endif
        ww(1)=(dsqrh*u1)
        zz(1)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
          
        y1vx(ipoint) = jetvx(ipoint) + f1stocvx(j)*ww(1) + &
         prefactor1*(f2vx(j)-f3vx)*zz(1) + &
	     0.25d0*h*(f2vx(j)+2.d0*f1vx(j)+f3vx)
	    
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1vx,npjet+1,jetvx)
      
      j=0
      y1xx(:)=0.d0
      y1st(:)=0.d0
      do ipoint=mystart,myend
	    y1xx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    y1st(ipoint) = jetst(ipoint) + h*f1st(j)
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      
      j=0
      do ipoint=mystart,myend
        call xpsys_pos(ipoint,y1xx,y1yy,y1zz,y1st,jetvx,jetvy,jetvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      y1xx(:)=0.d0
      do ipoint=mystart,myend
        y1xx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1,jetxx)
      
      j=0
      do ipoint=mystart,myend
        call xpsys_stress(ipoint,jetxx,jetyy,jetzz,y1st,jetvx,jetvy, &
         jetvz,jetvl,coulforce,f2st(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      y1st(:)=0.d0
      do ipoint=mystart,myend
	    y1st(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
	    j=j+1
      enddo
      call sum_world_darr(y1st,npjet+1,jetst)
      
      timesub=timesub+h
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      y1xx(:)=0.d0
	  y1yy(:)=0.d0
	  y1zz(:)=0.d0
	  y1st(:)=0.d0
	  y1vx(:)=0.d0
	  y1vy(:)=0.d0
	  y1vz(:)=0.d0
	  y2xx(:)=0.d0
	  y2yy(:)=0.d0
	  y2zz(:)=0.d0
	  y2st(:)=0.d0
	  y2vx(:)=0.d0
	  y2vy(:)=0.d0
	  y2vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         jetvl,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k, &
         fstocvx,fstocvy,fstocvz)
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        f1stocvx(j)=fstocvx
        f1stocvy(j)=fstocvy
        f1stocvz(j)=fstocvz
	    y1xx(ipoint) = jetxx(ipoint) + h*fxx
	    y1yy(ipoint) = jetyy(ipoint) + h*fyy
	    y1zz(ipoint) = jetzz(ipoint) + h*fzz
	    y1st(ipoint) = jetst(ipoint) + h*fst
	    y1vx(ipoint) = jetvx(ipoint) + h*fvx + dsqrh*fstocvx
	    y1vy(ipoint) = jetvy(ipoint) + h*fvy + dsqrh*fstocvy
	    y1vz(ipoint) = jetvz(ipoint) + h*fvz + dsqrh*fstocvz
	    y2xx(ipoint) = jetxx(ipoint) + h*fxx
	    y2yy(ipoint) = jetyy(ipoint) + h*fyy
	    y2zz(ipoint) = jetzz(ipoint) + h*fzz
	    y2st(ipoint) = jetst(ipoint) + h*fst
	    y2vx(ipoint) = jetvx(ipoint) + h*fvx - dsqrh*fstocvx
	    y2vy(ipoint) = jetvy(ipoint) + h*fvy - dsqrh*fstocvy
	    y2vz(ipoint) = jetvz(ipoint) + h*fvz - dsqrh*fstocvz
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1yy,npjet+1)
      call sum_world_darr(y1zz,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      call sum_world_darr(y1vx,npjet+1)
      call sum_world_darr(y1vy,npjet+1)
      call sum_world_darr(y1vz,npjet+1)
      call sum_world_darr(y2xx,npjet+1)
      call sum_world_darr(y2yy,npjet+1)
      call sum_world_darr(y2zz,npjet+1)
      call sum_world_darr(y2st,npjet+1)
      call sum_world_darr(y2vx,npjet+1)
      call sum_world_darr(y2vy,npjet+1)
      call sum_world_darr(y2vz,npjet+1)
      
      call smooth_charge(y1xx,y1yy,y1zz)
      call compute_posnoinserted(y1xx,y1yy,y1zz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y1xx, &
       y1yy,y1zz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,jetvl, &
         coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k,fstocvx, &
         fstocvy,fstocvz)
        f2xx(j)=fxx
        f2yy(j)=fyy
        f2zz(j)=fzz
        f2st(j)=fst
        f2vx(j)=fvx
        f2vy(j)=fvy
        f2vz(j)=fvz
	    j=j+1
      enddo
      call restore_charge()
      
      call smooth_charge(y2xx,y2yy,y2zz)
      call compute_posnoinserted(y2xx,y2yy,y2zz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y2xx, &
       y2yy,y2zz)
      j=0
      y1vx(:)=0.d0
      y1vy(:)=0.d0
      y1vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys(ipoint,y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,jetvl, &
         coulforce,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,timesub,k, &
         f3stocvx,f3stocvy,f3stocvz)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,1,1)
          u2=gaussian_history_value(k,ipoint,1,2)
        else
          u1=gaussian_buffer_value(ipoint,1,1)
          u2=gaussian_buffer_value(ipoint,1,2)
        endif
        ww(1)=(dsqrh*u1)
        zz(1)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,2,1)
          u2=gaussian_history_value(k,ipoint,2,2)
        else
          u1=gaussian_buffer_value(ipoint,2,1)
          u2=gaussian_buffer_value(ipoint,2,2)
        endif
        ww(2)=(dsqrh*u1)
        zz(2)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,3,1)
          u2=gaussian_history_value(k,ipoint,3,2)
        else
          u1=gaussian_buffer_value(ipoint,3,1)
          u2=gaussian_buffer_value(ipoint,3,2)
        endif
        ww(3)=(dsqrh*u1)
        zz(3)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
	    
        y1vx(ipoint) = jetvx(ipoint) + f1stocvx(j)*ww(1) + &
         prefactor1*(f2vx(j)-f3vx)*zz(1) + &
	     0.25d0*h*(f2vx(j)+2.d0*f1vx(j)+f3vx)
	    
	    y1vy(ipoint) = jetvy(ipoint) + f1stocvy(j)*ww(2) + &
	     prefactor1*(f2vy(j)-f3vy)*zz(2) + &
	     0.25d0*h*(f2vy(j)+2.d0*f1vy(j)+f3vy)
	     
	    y1vz(ipoint) = jetvz(ipoint) + f1stocvz(j)*ww(3) + &
	     prefactor1*(f2vz(j)-f3vz)*zz(3) + &
	     0.25d0*h*(f2vz(j)+2.d0*f1vz(j)+f3vz)
	    
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1vx,npjet+1,jetvx)
      call sum_world_darr(y1vy,npjet+1,jetvy)
      call sum_world_darr(y1vz,npjet+1,jetvz)
      
      j=0
      y1xx(:)=0.d0
      y1yy(:)=0.d0
      y1zz(:)=0.d0
      y1st(:)=0.d0
      do ipoint=mystart,myend
	    y1xx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    y1yy(ipoint) = jetyy(ipoint) + h*f1yy(j)
	    y1zz(ipoint) = jetzz(ipoint) + h*f1zz(j)
	    y1st(ipoint) = jetst(ipoint) + h*f1st(j)
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1yy,npjet+1)
      call sum_world_darr(y1zz,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      
      call compute_posnoinserted(y1xx,y1yy,y1zz)
      j=0
      do ipoint=mystart,myend
        call xpsys_pos(ipoint,y1xx,y1yy,y1zz,y1st,jetvx,jetvy,jetvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      y1xx(:)=0.d0
      y1yy(:)=0.d0
      y1zz(:)=0.d0
      do ipoint=mystart,myend
        y1xx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        y1yy(ipoint) = jetyy(ipoint) + (h/2.d0)*(f1yy(j)+f2yy(j))
        y1zz(ipoint) = jetzz(ipoint) + (h/2.d0)*(f1zz(j)+f2zz(j))
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1,jetxx)
      call sum_world_darr(y1yy,npjet+1,jetyy)
      call sum_world_darr(y1zz,npjet+1,jetzz)
      
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      j=0
      do ipoint=mystart,myend
        call xpsys_stress(ipoint,jetxx,jetyy,jetzz,y1st,jetvx,jetvy, &
         jetvz,jetvl,coulforce,f2st(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      y1st(:)=0.d0
      do ipoint=mystart,myend
	    y1st(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
	    j=j+1
      enddo
      call sum_world_darr(y1st,npjet+1,jetst)
      
      timesub=timesub+h
  end select
  
  return
      
 end subroutine platen
 
 subroutine eulsys_KV(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     first order accurate Euler scheme with Kelvin–Voigt model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
  
  
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
! service arrays
  double precision, allocatable, dimension (:), save ::  fxx
  double precision, allocatable, dimension (:), save ::  fyy
  double precision, allocatable, dimension (:), save ::  fzz
  double precision, allocatable, dimension (:), save ::  fst
  double precision, allocatable, dimension (:), save ::  fvx
  double precision, allocatable, dimension (:), save ::  fvy
  double precision, allocatable, dimension (:), save ::  fvz
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  
  integer :: ipoint,j
  
  logical, save :: lfirstsub=.true.
  
  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_euler))then
    call device_rk_step(scheme_euler,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(fxx)
        deallocate(fst)
        deallocate(fvx)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
      endif
      allocate(fxx(0:mxchunk))
      allocate(fst(0:mxchunk))
      allocate(fvx(0:mxnpjet))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(fxx)
        deallocate(fyy)
        deallocate(fzz)
        deallocate(fst)
        deallocate(fvx)
        deallocate(fvy)
        deallocate(fvz)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
      endif
      allocate(fxx(0:mxchunk))
      allocate(fyy(0:mxchunk))
      allocate(fzz(0:mxchunk))
      allocate(fst(0:mxchunk))
      allocate(fvx(0:mxnpjet))
      allocate(fvy(0:mxnpjet))
      allocate(fvz(0:mxnpjet))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif
  
  
! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      fvx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,fxx(j),fyy(j),fzz(j), &
         fvx(ipoint),fvy(ipoint),fvz(ipoint),timesub,k)
        j=j+1
      enddo
      call sum_world_darr(fvx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,fvx,fvy,fvz,fst(j), &
         timesub,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*fxx(j)
        yst(ipoint) = jetst(ipoint) + h*fst(j)
        yvx(ipoint) = jetvx(ipoint) + h*fvx(ipoint)
        j=j+1
      enddo
      
      call restore_charge()
      
      timesub=timesub+h
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz,timesub)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      fvx(:)=0.d0
      fvy(:)=0.d0
      fvz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,&
         jetvl,coulforce,fxx(j),fyy(j),fzz(j), &
         fvx(ipoint),fvy(ipoint),fvz(ipoint),timesub,k)
        j=j+1
      enddo
      call sum_world_darr(fvx,npjet+1)
      call sum_world_darr(fvy,npjet+1)
      call sum_world_darr(fvz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,fvx,fvy,fvz,fst(j), &
         timesub,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*fxx(j)
        yyy(ipoint) = jetyy(ipoint) + h*fyy(j)
        yzz(ipoint) = jetzz(ipoint) + h*fzz(j)
        yst(ipoint) = jetst(ipoint) + h*fst(j)
        yvx(ipoint) = jetvx(ipoint) + h*fvx(ipoint)
        yvy(ipoint) = jetvy(ipoint) + h*fvy(ipoint)
        yvz(ipoint) = jetvz(ipoint) + h*fvz(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select
  
  
  return
  
  
 end subroutine eulsys_KV
 
 subroutine rk2sys_KV(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     second order accurate Heun scheme with Kelvin–Voigt model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
  
  
  integer,intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision,intent(in) :: h
  
  integer :: ipoint,j
! service arrays
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  
  logical, save :: lfirstsub=.true.
  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_rk2))then
    call device_rk_step(scheme_rk2,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxnpjet))
      allocate(f2xx(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxnpjet))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1yy(0:mxchunk))
      allocate(f1zz(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxnpjet))
      allocate(f1vy(0:mxnpjet))
      allocate(f1vz(0:mxnpjet))
      allocate(f2xx(0:mxchunk))
      allocate(f2yy(0:mxchunk))
      allocate(f2zz(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxnpjet))
      allocate(f2vy(0:mxnpjet))
      allocate(f2vz(0:mxnpjet))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif
  
  
! select the proper system type
  select case(systype)
    case(1)
!     1°step
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      f1vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1xx(j),f1yy(j),f1zz(j), &
         f1vx(ipoint),f1vy(ipoint),f1vz(ipoint),timesub,k)
        j=j+1
      enddo
      call sum_world_darr(f1vx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1vx,f1vy,f1vz,f1st(j), &
         timesub,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f1xx(j)
        yst(ipoint) = jetst(ipoint) + h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f1vx(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     2°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f2vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,&
         jetvl,coulforce,f1xx(j),f2yy(j),f2zz(j), &
         f2vx(ipoint),f2vy(ipoint),f2vz(ipoint),timesub+h,k)
        j=j+1
      enddo
      call sum_world_darr(f2vx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst, &
         yvx,yvy,yvz,jetvl,coulforce,f2vx,f2vy,f2vz,f2st(j), &
         timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        yst(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/2.d0)* &
         (f1vx(ipoint)+f2vx(ipoint))
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      timesub=timesub+h
    case default
!     1°step
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      f1vx(:)=0.d0
      f1vy(:)=0.d0
      f1vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1xx(j),f1yy(j),f1zz(j), &
         f1vx(ipoint),f1vy(ipoint),f1vz(ipoint),timesub,k)
        j=j+1
      enddo
      call sum_world_darr(f1vx,npjet+1)
      call sum_world_darr(f1vy,npjet+1)
      call sum_world_darr(f1vz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1vx,f1vy,f1vz,f1st(j), &
         timesub,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f1xx(j)
        yyy(ipoint) = jetyy(ipoint) + h*f1yy(j)
        yzz(ipoint) = jetzz(ipoint) + h*f1zz(j)
        yst(ipoint) = jetst(ipoint) + h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f1vx(ipoint)
        yvy(ipoint) = jetvy(ipoint) + h*f1vy(ipoint)
        yvz(ipoint) = jetvz(ipoint) + h*f1vz(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     2°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f2vx(:)=0.d0
      f2vy(:)=0.d0
      f2vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,&
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j), &
         f2vx(ipoint),f2vy(ipoint),f2vz(ipoint),timesub+h,k)
        j=j+1
      enddo
      call sum_world_darr(f2vx,npjet+1)
      call sum_world_darr(f2vy,npjet+1)
      call sum_world_darr(f2vz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst, &
         yvx,yvy,yvz,jetvl,coulforce,f2vx,f2vy,f2vz,f2st(j), &
         timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        yyy(ipoint) = jetyy(ipoint) + (h/2.d0)*(f1yy(j)+f2yy(j))
        yzz(ipoint) = jetzz(ipoint) + (h/2.d0)*(f1zz(j)+f2zz(j))
        yst(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/2.d0)* &
         (f1vx(ipoint)+f2vx(ipoint))
        yvy(ipoint) = jetvy(ipoint) + (h/2.d0)* &
         (f1vy(ipoint)+f2vy(ipoint))
        yvz(ipoint) = jetvz(ipoint) + (h/2.d0)* &
         (f1vz(ipoint)+f2vz(ipoint))
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select

  return
      
 end subroutine rk2sys_KV
 
 subroutine rk4sys_KV(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     fourth order accurate Runge-Kutta scheme with Kelvin–Voigt model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
! service arrays
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  f3xx
  double precision, allocatable, dimension (:), save ::  f3yy
  double precision, allocatable, dimension (:), save ::  f3zz
  double precision, allocatable, dimension (:), save ::  f3st
  double precision, allocatable, dimension (:), save ::  f3vx
  double precision, allocatable, dimension (:), save ::  f3vy
  double precision, allocatable, dimension (:), save ::  f3vz
  double precision, allocatable, dimension (:), save ::  f4xx
  double precision, allocatable, dimension (:), save ::  f4yy
  double precision, allocatable, dimension (:), save ::  f4zz
  double precision, allocatable, dimension (:), save ::  f4st
  double precision, allocatable, dimension (:), save ::  f4vx
  double precision, allocatable, dimension (:), save ::  f4vy
  double precision, allocatable, dimension (:), save ::  f4vz
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  integer :: ipoint,dm,nv,j
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
  logical, save :: lfirstsub=.true.
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_rk4))then
    call device_rk_step(scheme_rk4,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f3xx)
        deallocate(f3st)
        deallocate(f3vx)
        deallocate(f4xx)
        deallocate(f4st)
        deallocate(f4vx)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxnpjet))
      allocate(f2xx(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxnpjet))
      allocate(f3xx(0:mxchunk))
      allocate(f3st(0:mxchunk))
      allocate(f3vx(0:mxnpjet))
      allocate(f4xx(0:mxchunk))
      allocate(f4st(0:mxchunk))
      allocate(f4vx(0:mxnpjet))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(f3xx)
        deallocate(f3yy)
        deallocate(f3zz)
        deallocate(f3st)
        deallocate(f3vx)
        deallocate(f3vy)
        deallocate(f3vz)
        deallocate(f4xx)
        deallocate(f4yy)
        deallocate(f4zz)
        deallocate(f4st)
        deallocate(f4vx)
        deallocate(f4vy)
        deallocate(f4vz)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1yy(0:mxchunk))
      allocate(f1zz(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxnpjet))
      allocate(f1vy(0:mxnpjet))
      allocate(f1vz(0:mxnpjet))
      allocate(f2xx(0:mxchunk))
      allocate(f2yy(0:mxchunk))
      allocate(f2zz(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxnpjet))
      allocate(f2vy(0:mxnpjet))
      allocate(f2vz(0:mxnpjet))
      allocate(f3xx(0:mxchunk))
      allocate(f3yy(0:mxchunk))
      allocate(f3zz(0:mxchunk))
      allocate(f3st(0:mxchunk))
      allocate(f3vx(0:mxnpjet))
      allocate(f3vy(0:mxnpjet))
      allocate(f3vz(0:mxnpjet))
      allocate(f4xx(0:mxchunk))
      allocate(f4yy(0:mxchunk))
      allocate(f4zz(0:mxchunk))
      allocate(f4st(0:mxchunk))
      allocate(f4vx(0:mxnpjet))
      allocate(f4vy(0:mxnpjet))
      allocate(f4vz(0:mxnpjet))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif
  
! select the proper system type
  select case(systype)
    case(1)
!     1°step
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      f1vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1xx(j),f1yy(j),f1zz(j), &
         f1vx(ipoint),f1vy(ipoint),f1vz(ipoint),timesub,k)
        j=j+1
      enddo
      call sum_world_darr(f1vx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1vx,f1vy,f1vz,f1st(j), &
         timesub,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     2°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f2vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j), &
         f2vx(ipoint),f2vy(ipoint),f2vz(ipoint),timesub+h/2.d0,k)
        j=j+1
      enddo
      call sum_world_darr(f2vx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl, &
         coulforce,f2vx,f2vy,f2vz,f2st(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     3°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f3vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f3xx(j),f3yy(j),f3zz(j), &
         f3vx(ipoint),f3vy(ipoint),f3vz(ipoint),timesub+h/2.d0,k)
        j=j+1
      enddo
      call sum_world_darr(f3vx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl, &
         coulforce,f3vx,f3vy,f3vz,f3st(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
        yst(ipoint) = jetst(ipoint) + h*f3st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f3vx(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
!     4°step       
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f4vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f4xx(j),f4yy(j),f4zz(j), &
         f4vx(ipoint),f4vy(ipoint),f4vz(ipoint),timesub+h,k)
        j=j+1
      enddo
      call sum_world_darr(f4vx,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl, &
         coulforce,f4vx,f4vy,f4vz,f4st(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/6.d0)*(f1xx(j)+ &
         2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
        yst(ipoint) = jetst(ipoint) + (h/6.d0)*(f1st(j)+ &
         2.d0*(f2st(j)+f3st(j))+f4st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(ipoint)+ &
         2.d0*(f2vx(ipoint)+f3vx(ipoint))+f4vx(ipoint))
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      timesub=timesub+h
    case default
!     1°step
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz)
      j=0
      f1vx(:)=0.d0
      f1vy(:)=0.d0
      f1vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1xx(j),f1yy(j),f1zz(j), &
         f1vx(ipoint),f1vy(ipoint),f1vz(ipoint),timesub,k)
        j=j+1
      enddo
      call sum_world_darr(f1vx,npjet+1)
      call sum_world_darr(f1vy,npjet+1)
      call sum_world_darr(f1vz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,coulforce,f1vx,f1vy,f1vz,f1st(j), &
         timesub,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f1yy(j)
        yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f1zz(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(ipoint)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f1vy(ipoint)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f1vz(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     2°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f2vx(:)=0.d0
      f2vy(:)=0.d0
      f2vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j), &
         f2vx(ipoint),f2vy(ipoint),f2vz(ipoint),timesub+h/2.d0,k)
        j=j+1
      enddo
      call sum_world_darr(f2vx,npjet+1)
      call sum_world_darr(f2vy,npjet+1)
      call sum_world_darr(f2vz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl, &
         coulforce,f2vx,f2vy,f2vz,f2st(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
        yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f2yy(j)
        yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f2zz(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(ipoint)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f2vy(ipoint)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f2vz(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     3°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f3vx(:)=0.d0
      f3vy(:)=0.d0
      f3vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f3xx(j),f3yy(j),f3zz(j), &
         f3vx(ipoint),f3vy(ipoint),f3vz(ipoint),timesub+h/2.d0,k)
        j=j+1
      enddo
      call sum_world_darr(f3vx,npjet+1)
      call sum_world_darr(f3vy,npjet+1)
      call sum_world_darr(f3vz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl, &
         coulforce,f3vx,f3vy,f3vz,f3st(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
        yyy(ipoint) = jetyy(ipoint) + h*f3yy(j)
        yzz(ipoint) = jetzz(ipoint) + h*f3zz(j)
        yst(ipoint) = jetst(ipoint) + h*f3st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f3vx(ipoint)
        yvy(ipoint) = jetvy(ipoint) + h*f3vy(ipoint)
        yvz(ipoint) = jetvz(ipoint) + h*f3vz(ipoint)
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
!     4°step       
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      f4vx(:)=0.d0
      f4vy(:)=0.d0
      f4vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f4xx(j),f4yy(j),f4zz(j), &
         f4vx(ipoint),f4vy(ipoint),f4vz(ipoint),timesub+h,k)
        j=j+1
      enddo
      call sum_world_darr(f4vx,npjet+1)
      call sum_world_darr(f4vy,npjet+1)
      call sum_world_darr(f4vz,npjet+1)
      j=0
      do ipoint=mystart,myend
        call xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl, &
         coulforce,f4vx,f4vy,f4vz,f4st(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/6.d0)*(f1xx(j)+ &
         2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
        yyy(ipoint) = jetyy(ipoint) + (h/6.d0)*(f1yy(j)+ &
         2.d0*(f2yy(j)+f3yy(j))+f4yy(j))
        yzz(ipoint) = jetzz(ipoint) + (h/6.d0)*(f1zz(j)+ &
         2.d0*(f2zz(j)+f3zz(j))+f4zz(j))
        yst(ipoint) = jetst(ipoint) + (h/6.d0)*(f1st(j)+ &
         2.d0*(f2st(j)+f3st(j))+f4st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(ipoint)+ &
         2.d0*(f2vx(ipoint)+f3vx(ipoint))+f4vx(ipoint))
        yvy(ipoint) = jetvy(ipoint) + (h/6.d0)*(f1vy(ipoint)+ &
         2.d0*(f2vy(ipoint)+f3vy(ipoint))+f4vy(ipoint))
        yvz(ipoint) = jetvz(ipoint) + (h/6.d0)*(f1vz(ipoint)+ &
         2.d0*(f2vz(ipoint)+f3vz(ipoint))+f4vz(ipoint))
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select
  
  return
  
 end subroutine rk4sys_KV

 
 subroutine eulsys_ev(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     first order accurate Euler scheme
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2015
!     
!***********************************************************************
  
  implicit none
  
  
  
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
! service arrays
  double precision, allocatable, dimension (:), save ::  fxx
  double precision, allocatable, dimension (:), save ::  fyy
  double precision, allocatable, dimension (:), save ::  fzz
  double precision, allocatable, dimension (:), save ::  fst
  double precision, allocatable, dimension (:), save ::  fvx
  double precision, allocatable, dimension (:), save ::  fvy
  double precision, allocatable, dimension (:), save ::  fvz
  double precision, allocatable, dimension (:), save ::  fev
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  double precision, allocatable, dimension (:), save ::  yev
  
  integer :: ipoint,j
  
  logical, save :: lfirstsub=.true.

  
  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_euler))then
    call device_rk_step(scheme_euler,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(fxx)
        deallocate(fst)
        deallocate(fvx)
        deallocate(fev)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yev)
      endif
      allocate(fxx(0:mxchunk))
      allocate(fst(0:mxchunk))
      allocate(fvx(0:mxchunk))
      allocate(fev(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yev(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(fxx)
        deallocate(fyy)
        deallocate(fzz)
        deallocate(fst)
        deallocate(fvx)
        deallocate(fvy)
        deallocate(fvz)
        deallocate(fev)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
        deallocate(yev)
      endif
      allocate(fxx(0:mxchunk))
      allocate(fyy(0:mxchunk))
      allocate(fzz(0:mxchunk))
      allocate(fst(0:mxchunk))
      allocate(fvx(0:mxchunk))
      allocate(fvy(0:mxchunk))
      allocate(fvz(0:mxchunk))
      allocate(fev(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
      allocate(yev(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif
  
  
! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         jetvl,jetve,coulforce,fxx(j),fyy(j),fzz(j),fst(j), &
         fvx(j),fvy(j),fvz(j),fev(j),timesub,k)
        j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*fxx(j)
        yst(ipoint) = jetst(ipoint) + h*fst(j)
        yvx(ipoint) = jetvx(ipoint) + h*fvx(j)
        yev(ipoint) = jetve(ipoint) + h*fev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      
      call restore_charge()
      
      timesub=timesub+h
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yev,npjet+1,jetve)
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz,timesub)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,&
         jetvl,jetve,coulforce,fxx(j),fyy(j),fzz(j),fst(j), &
         fvx(j),fvy(j),fvz(j),fev(j),timesub,k)
        j=j+1
      enddo
      call restore_charge()
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*fxx(j)
        yyy(ipoint) = jetyy(ipoint) + h*fyy(j)
        yzz(ipoint) = jetzz(ipoint) + h*fzz(j)
        yst(ipoint) = jetst(ipoint) + h*fst(j)
        yvx(ipoint) = jetvx(ipoint) + h*fvx(j)
        yvy(ipoint) = jetvy(ipoint) + h*fvy(j)
        yvz(ipoint) = jetvz(ipoint) + h*fvz(j)
        yev(ipoint) = jetve(ipoint) + h*fev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      call sum_world_darr(yev,npjet+1,jetve)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select
  
  
  return
  
  
 end subroutine eulsys_ev
 
 subroutine rk2sys_ev(timesub,h,k)
  
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     second order accurate Heun scheme
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  
  
  integer,intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision,intent(in) :: h
  
  integer :: ipoint,j
! service arrays
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f1ev
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  f2ev
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  double precision, allocatable, dimension (:), save ::  yev
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  double precision ::  fev
  
  logical, save :: lfirstsub=.true.

  
#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_rk2))then
    call device_rk_step(scheme_rk2,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1ev)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2ev)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yev)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f1ev(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f2ev(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yev(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f1ev)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(f2ev)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
        deallocate(yev)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1yy(0:mxchunk))
      allocate(f1zz(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f1vy(0:mxchunk))
      allocate(f1vz(0:mxchunk))
      allocate(f1ev(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2yy(0:mxchunk))
      allocate(f2zz(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f2vy(0:mxchunk))
      allocate(f2vz(0:mxchunk))
      allocate(f2ev(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
      allocate(yev(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif
  
  
! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev, &
          timesub,k)
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        f1ev(j)=fev
        yxx(ipoint) = jetxx(ipoint) + h*f1xx(j)
        yst(ipoint) = jetst(ipoint) + h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f1vx(j)
        yev(ipoint) = jetve(ipoint) + h*f1ev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yev,npjet+1)
      
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),f2ev(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        yst(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/2.d0)*(f1vx(j)+f2vx(j))
        yev(ipoint) = jetve(ipoint) + (h/2.d0)*(f1ev(j)+f2ev(j))
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yev,npjet+1,jetve)
      timesub=timesub+h
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev, &
          timesub,k)
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        f1ev(j)=fev
	    yxx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    yyy(ipoint) = jetyy(ipoint) + h*f1yy(j)
	    yzz(ipoint) = jetzz(ipoint) + h*f1zz(j)
	    yst(ipoint) = jetst(ipoint) + h*f1st(j)
	    yvx(ipoint) = jetvx(ipoint) + h*f1vx(j)
	    yvy(ipoint) = jetvy(ipoint) + h*f1vy(j)
	    yvz(ipoint) = jetvz(ipoint) + h*f1vz(j)
        yev(ipoint) = jetve(ipoint) + h*f1ev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call restore_charge()
      
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
      call sum_world_darr(yev,npjet+1)
      
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),f2ev(j),timesub+h,k)
        j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        yyy(ipoint) = jetyy(ipoint) + (h/2.d0)*(f1yy(j)+f2yy(j))
        yzz(ipoint) = jetzz(ipoint) + (h/2.d0)*(f1zz(j)+f2zz(j))
	    yst(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
	    yvx(ipoint) = jetvx(ipoint) + (h/2.d0)*(f1vx(j)+f2vx(j))
	    yvy(ipoint) = jetvy(ipoint) + (h/2.d0)*(f1vy(j)+f2vy(j))
	    yvz(ipoint) = jetvz(ipoint) + (h/2.d0)*(f1vz(j)+f2vz(j))
        yev(ipoint) = jetve(ipoint) + (h/2.d0)*(f1ev(j)+f2ev(j))
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      call sum_world_darr(yev,npjet+1,jetve)
      timesub=timesub+h
      call compute_posnoinserted(jetxx,jetyy,jetzz)
  end select

  return
      
 end subroutine rk2sys_ev

 subroutine rk4sys_ev(timesub,h,k)

!***********************************************************************
!     
!     JETSPIN subroutine for integrating the system by the 
!     fourth order accurate Runge-Kutta scheme
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
! service arrays
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f1ev
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  f2ev
  double precision, allocatable, dimension (:), save ::  f3xx
  double precision, allocatable, dimension (:), save ::  f3yy
  double precision, allocatable, dimension (:), save ::  f3zz
  double precision, allocatable, dimension (:), save ::  f3st
  double precision, allocatable, dimension (:), save ::  f3vx
  double precision, allocatable, dimension (:), save ::  f3vy
  double precision, allocatable, dimension (:), save ::  f3vz
  double precision, allocatable, dimension (:), save ::  f3ev
  double precision, allocatable, dimension (:), save ::  f4xx
  double precision, allocatable, dimension (:), save ::  f4yy
  double precision, allocatable, dimension (:), save ::  f4zz
  double precision, allocatable, dimension (:), save ::  f4st
  double precision, allocatable, dimension (:), save ::  f4vx
  double precision, allocatable, dimension (:), save ::  f4vy
  double precision, allocatable, dimension (:), save ::  f4vz
  double precision, allocatable, dimension (:), save ::  f4ev
  double precision, allocatable, dimension (:), save ::  yxx
  double precision, allocatable, dimension (:), save ::  yyy
  double precision, allocatable, dimension (:), save ::  yzz
  double precision, allocatable, dimension (:), save ::  yst
  double precision, allocatable, dimension (:), save ::  yvx
  double precision, allocatable, dimension (:), save ::  yvy
  double precision, allocatable, dimension (:), save ::  yvz
  double precision, allocatable, dimension (:), save ::  yev
  integer :: ipoint,dm,nv,j
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
 logical, save :: lfirstsub=.true.

  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  double precision ::  fev


#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod); below it,
! and for the options that step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_rk4))then
    call device_rk_step(scheme_rk4,timesub,h,k)
    return
  endif
#endif
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1ev)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2ev)
        deallocate(f3xx)
        deallocate(f3st)
        deallocate(f3vx)
        deallocate(f3ev)
        deallocate(f4xx)
        deallocate(f4st)
        deallocate(f4vx)
        deallocate(f4ev)
        deallocate(yxx)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yev)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f1ev(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f2ev(0:mxchunk))
      allocate(f3xx(0:mxchunk))
      allocate(f3st(0:mxchunk))
      allocate(f3vx(0:mxchunk))
      allocate(f3ev(0:mxchunk))
      allocate(f4xx(0:mxchunk))
      allocate(f4st(0:mxchunk))
      allocate(f4vx(0:mxchunk))
      allocate(f4ev(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yev(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f1ev)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(f2ev)
        deallocate(f3xx)
        deallocate(f3yy)
        deallocate(f3zz)
        deallocate(f3st)
        deallocate(f3vx)
        deallocate(f3vy)
        deallocate(f3vz)
        deallocate(f3ev)
        deallocate(f4xx)
        deallocate(f4yy)
        deallocate(f4zz)
        deallocate(f4st)
        deallocate(f4vx)
        deallocate(f4vy)
        deallocate(f4vz)
        deallocate(f4ev)
        deallocate(yxx)
        deallocate(yyy)
        deallocate(yzz)
        deallocate(yst)
        deallocate(yvx)
        deallocate(yvy)
        deallocate(yvz)
        deallocate(yev)
      endif
      allocate(f1xx(0:mxchunk))
      allocate(f1yy(0:mxchunk))
      allocate(f1zz(0:mxchunk))
      allocate(f1st(0:mxchunk))
      allocate(f1vx(0:mxchunk))
      allocate(f1vy(0:mxchunk))
      allocate(f1vz(0:mxchunk))
      allocate(f1ev(0:mxchunk))
      allocate(f2xx(0:mxchunk))
      allocate(f2yy(0:mxchunk))
      allocate(f2zz(0:mxchunk))
      allocate(f2st(0:mxchunk))
      allocate(f2vx(0:mxchunk))
      allocate(f2vy(0:mxchunk))
      allocate(f2vz(0:mxchunk))
      allocate(f2ev(0:mxchunk))
      allocate(f3xx(0:mxchunk))
      allocate(f3yy(0:mxchunk))
      allocate(f3zz(0:mxchunk))
      allocate(f3st(0:mxchunk))
      allocate(f3vx(0:mxchunk))
      allocate(f3vy(0:mxchunk))
      allocate(f3vz(0:mxchunk))
      allocate(f3ev(0:mxchunk))
      allocate(f4xx(0:mxchunk))
      allocate(f4yy(0:mxchunk))
      allocate(f4zz(0:mxchunk))
      allocate(f4st(0:mxchunk))
      allocate(f4vx(0:mxchunk))
      allocate(f4vy(0:mxchunk))
      allocate(f4vz(0:mxchunk))
      allocate(f4ev(0:mxchunk))
      allocate(yxx(0:mxnpjet))
      allocate(yyy(0:mxnpjet))
      allocate(yzz(0:mxnpjet))
      allocate(yst(0:mxnpjet))
      allocate(yvx(0:mxnpjet))
      allocate(yvy(0:mxnpjet))
      allocate(yvz(0:mxnpjet))
      allocate(yev(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif
! select the proper system type
  select case(systype)
    case(1)
!     1°step
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev_maxwell(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev, &
          timesub,k)
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        f1ev(j)=fev
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
        yev(ipoint) = jetve(ipoint) + 0.5d0*h*f1ev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yev,npjet+1)
!     2°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),f2ev(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(j)
        yev(ipoint) = jetve(ipoint) + 0.5d0*h*f2ev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yev,npjet+1)
!     3°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
         f3vx(j),f3vy(j),f3vz(j),f3ev(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
        yst(ipoint) = jetst(ipoint) + h*f3st(j)
        yvx(ipoint) = jetvx(ipoint) + h*f3vx(j)
        yev(ipoint) = jetve(ipoint) + h*f3ev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yev,npjet+1)
!     4°step
      call smooth_charge(yxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
         f4vx(j),f4vy(j),f4vz(j),f4ev(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        yxx(ipoint) = jetxx(ipoint) + (h/6.d0)*(f1xx(j)+ &
         2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
        yst(ipoint) = jetst(ipoint) + (h/6.d0)*(f1st(j)+ &
         2.d0*(f2st(j)+f3st(j))+f4st(j))
        yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(j)+ &
         2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
        yev(ipoint) = jetve(ipoint) + (h/6.d0)*(f1ev(j)+ &
         2.d0*(f2ev(j)+f3ev(j))+f4ev(j))
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yev,npjet+1,jetve)
      timesub=timesub+h
    case default
!     1°step
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl, &
       jetxx,jetyy,jetzz,jetve)
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,&
         jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev, &
         timesub,k) 
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        f1ev(j)=fev
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f1yy(j)
        yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f1zz(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f1vy(j)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f1vz(j)
        yev(ipoint) = jetve(ipoint) + 0.5d0*h*f1ev(j)
        if((yev(ipoint)/jetvl(ipoint))<evlim)then
          yev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
      call sum_world_darr(yev,npjet+1)
!     2°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),f2ev(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
! Below the 100-bead gate of the device step the stage updates stay on the
! host, as in the CPU build.  Until 2026-10-06 this host path also
! recomputed the Maxwell evaporative stress and the stage updates on the
! device, copying the arrays at every call.
        do ipoint=mystart,myend
          yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
          yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f2yy(j)
          yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f2zz(j)
          yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
          yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(j)
          yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f2vy(j)
          yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f2vz(j)
          yev(ipoint) = jetve(ipoint) + 0.5d0*h*f2ev(j)
          if((yev(ipoint)/jetvl(ipoint))<evlim)then
            yev(ipoint)=jetvl(ipoint)*evlim
          endif
          j=j+1
        enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
      call sum_world_darr(yev,npjet+1)
!     3°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
         f3vx(j),f3vy(j),f3vz(j),f3ev(j),timesub+h/2.d0,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
        do ipoint=mystart,myend
          yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
          yyy(ipoint) = jetyy(ipoint) + h*f3yy(j)
          yzz(ipoint) = jetzz(ipoint) + h*f3zz(j)
          yst(ipoint) = jetst(ipoint) + h*f3st(j)
          yvx(ipoint) = jetvx(ipoint) + h*f3vx(j)
          yvy(ipoint) = jetvy(ipoint) + h*f3vy(j)
          yvz(ipoint) = jetvz(ipoint) + h*f3vz(j)
          yev(ipoint) = jetve(ipoint) + h*f3ev(j)
          if((yev(ipoint)/jetvl(ipoint))<evlim)then
            yev(ipoint)=jetvl(ipoint)*evlim
          endif
          j=j+1
        enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1)
      call sum_world_darr(yyy,npjet+1)
      call sum_world_darr(yzz,npjet+1)
      call sum_world_darr(yst,npjet+1)
      call sum_world_darr(yvx,npjet+1)
      call sum_world_darr(yvy,npjet+1)
      call sum_world_darr(yvz,npjet+1)
      call sum_world_darr(yev,npjet+1)
!     4°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
         f4vx(j),f4vy(j),f4vz(j),f4ev(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
        do ipoint=mystart,myend
          yxx(ipoint) = jetxx(ipoint) + (h/6.d0)*(f1xx(j)+ &
           2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
          yyy(ipoint) = jetyy(ipoint) + (h/6.d0)*(f1yy(j)+ &
           2.d0*(f2yy(j)+f3yy(j))+f4yy(j))
          yzz(ipoint) = jetzz(ipoint) + (h/6.d0)*(f1zz(j)+ &
           2.d0*(f2zz(j)+f3zz(j))+f4zz(j))
          yst(ipoint) = jetst(ipoint) + (h/6.d0)*(f1st(j)+ &
           2.d0*(f2st(j)+f3st(j))+f4st(j))
          yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(j)+ &
           2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
          yvy(ipoint) = jetvy(ipoint) + (h/6.d0)*(f1vy(j)+ &
           2.d0*(f2vy(j)+f3vy(j))+f4vy(j))
          yvz(ipoint) = jetvz(ipoint) + (h/6.d0)*(f1vz(j)+ &
           2.d0*(f2vz(j)+f3vz(j))+f4vz(j))
          yev(ipoint) = jetve(ipoint) + (h/6.d0)*(f1ev(j)+ &
           2.d0*(f2ev(j)+f3ev(j))+f4ev(j))
          if((yev(ipoint)/jetvl(ipoint))<evlim)then
            yev(ipoint)=jetvl(ipoint)*evlim
          endif
          j=j+1
        enddo
      call restore_charge()
      call sum_world_darr(yxx,npjet+1,jetxx)
      call sum_world_darr(yyy,npjet+1,jetyy)
      call sum_world_darr(yzz,npjet+1,jetzz)
      call sum_world_darr(yst,npjet+1,jetst)
      call sum_world_darr(yvx,npjet+1,jetvx)
      call sum_world_darr(yvy,npjet+1,jetvy)
      call sum_world_darr(yvz,npjet+1,jetvz)
      call sum_world_darr(yev,npjet+1,jetve)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      timesub=timesub+h
  end select
  
  return
  
 end subroutine rk4sys_ev
 
 subroutine platen_ev(timesub,h,k)
 
!***********************************************************************
!     
!     JETSPIN subroutine for integrating the stochastic equation
!     of motion by the 1.5° order accurate Platen scheme
!     for the velocity stochastic part and by the second order accurate
!     Heun scheme for deterministic part
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, allocatable, dimension (:), save ::  f1xx
  double precision, allocatable, dimension (:), save ::  f1yy
  double precision, allocatable, dimension (:), save ::  f1zz
  double precision, allocatable, dimension (:), save ::  f1st
  double precision, allocatable, dimension (:), save ::  f1vx
  double precision, allocatable, dimension (:), save ::  f1vy
  double precision, allocatable, dimension (:), save ::  f1vz
  double precision, allocatable, dimension (:), save ::  f1ev
  double precision, allocatable, dimension (:), save ::  f1stocvx
  double precision, allocatable, dimension (:), save ::  f1stocvy
  double precision, allocatable, dimension (:), save ::  f1stocvz
  double precision, allocatable, dimension (:), save ::  f2xx
  double precision, allocatable, dimension (:), save ::  f2yy
  double precision, allocatable, dimension (:), save ::  f2zz
  double precision, allocatable, dimension (:), save ::  f2st
  double precision, allocatable, dimension (:), save ::  f2vx
  double precision, allocatable, dimension (:), save ::  f2vy
  double precision, allocatable, dimension (:), save ::  f2vz
  double precision, allocatable, dimension (:), save ::  f2ev
  double precision, allocatable, dimension (:), save ::  y1xx
  double precision, allocatable, dimension (:), save ::  y1yy
  double precision, allocatable, dimension (:), save ::  y1zz
  double precision, allocatable, dimension (:), save ::  y1st
  double precision, allocatable, dimension (:), save ::  y1vx
  double precision, allocatable, dimension (:), save ::  y1vy
  double precision, allocatable, dimension (:), save ::  y1vz
  double precision, allocatable, dimension (:), save ::  y1ev
  double precision, allocatable, dimension (:), save ::  y2xx
  double precision, allocatable, dimension (:), save ::  y2yy
  double precision, allocatable, dimension (:), save ::  y2zz
  double precision, allocatable, dimension (:), save ::  y2st
  double precision, allocatable, dimension (:), save ::  y2vx
  double precision, allocatable, dimension (:), save ::  y2vy
  double precision, allocatable, dimension (:), save ::  y2vz
  double precision, allocatable, dimension (:), save ::  y2ev
  double precision, allocatable, dimension (:), save ::  d3xx,d3yy,d3zz
  double precision, allocatable, dimension (:), save ::  d3st,d3vx,d3vy,d3vz
  double precision, allocatable, dimension (:), save ::  d3ev
  
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  
  integer :: ipoint,dm,nv,j
  double precision :: dsqrh,tsqh,prefactor1,zztang,u1,u2
  double precision, dimension(1:3) :: ww,zz,utang
  
  logical, save :: lfirstsub=.true.
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  double precision ::  fev
  double precision ::  fstocvx
  double precision ::  fstocvy
  double precision ::  fstocvz
  
  double precision ::  f3xx
  double precision ::  f3yy
  double precision ::  f3zz
  double precision ::  f3st
  double precision ::  f3vx
  double precision ::  f3vy
  double precision ::  f3vz
  double precision ::  f3ev
  double precision ::  f3stocvx
  double precision ::  f3stocvy
  double precision ::  f3stocvz

  
! check and eventually reallocate the service arrays
  if(doallocate)then
    select case(systype)
    case(1)
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1ev)
        deallocate(f1stocvx)
        deallocate(f2xx)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2ev)
        deallocate(y1xx)
        deallocate(y1st)
        deallocate(y1vx)
        deallocate(y1ev)
        deallocate(y2xx)
        deallocate(y2st)
        deallocate(y2vx)
        deallocate(y2ev)
      endif
      allocate(f1xx(0:mxnpjet))
      allocate(f1st(0:mxnpjet))
      allocate(f1vx(0:mxnpjet))
      allocate(f1ev(0:mxnpjet))
      allocate(f1stocvx(0:mxnpjet))
      allocate(f2xx(0:mxnpjet))
      allocate(f2st(0:mxnpjet))
      allocate(f2vx(0:mxnpjet))
      allocate(f2ev(0:mxnpjet))
      allocate(y1xx(0:mxnpjet))
      allocate(y1st(0:mxnpjet))
      allocate(y1vx(0:mxnpjet))
      allocate(y1ev(0:mxnpjet))
      allocate(y2xx(0:mxnpjet))
      allocate(y2st(0:mxnpjet))
      allocate(y2vx(0:mxnpjet))
      allocate(y2ev(0:mxnpjet))
    case default
      if(.not.lfirstsub)then
        deallocate(f1xx)
        deallocate(f1yy)
        deallocate(f1zz)
        deallocate(f1st)
        deallocate(f1vx)
        deallocate(f1vy)
        deallocate(f1vz)
        deallocate(f1ev)
        deallocate(f1stocvx)
        deallocate(f1stocvy)
        deallocate(f1stocvz)
        deallocate(f2xx)
        deallocate(f2yy)
        deallocate(f2zz)
        deallocate(f2st)
        deallocate(f2vx)
        deallocate(f2vy)
        deallocate(f2vz)
        deallocate(f2ev)
        deallocate(y1xx)
        deallocate(y1yy)
        deallocate(y1zz)
        deallocate(y1st)
        deallocate(y1vx)
        deallocate(y1vy)
        deallocate(y1vz)
        deallocate(y1ev)
        deallocate(y2xx)
        deallocate(y2yy)
        deallocate(y2zz)
        deallocate(y2st)
        deallocate(y2vx)
        deallocate(y2vy)
        deallocate(y2vz)
        deallocate(y2ev)
        deallocate(d3xx,d3yy,d3zz,d3st,d3vx,d3vy,d3vz,d3ev)
      endif
      allocate(f1xx(0:mxnpjet))
      allocate(f1yy(0:mxnpjet))
      allocate(f1zz(0:mxnpjet))
      allocate(f1st(0:mxnpjet))
      allocate(f1vx(0:mxnpjet))
      allocate(f1vy(0:mxnpjet))
      allocate(f1vz(0:mxnpjet))
      allocate(f1ev(0:mxnpjet))
      allocate(f1stocvx(0:mxnpjet))
      allocate(f1stocvy(0:mxnpjet))
      allocate(f1stocvz(0:mxnpjet))
      allocate(f2xx(0:mxnpjet))
      allocate(f2yy(0:mxnpjet))
      allocate(f2zz(0:mxnpjet))
      allocate(f2st(0:mxnpjet))
      allocate(f2vx(0:mxnpjet))
      allocate(f2vy(0:mxnpjet))
      allocate(f2vz(0:mxnpjet))
      allocate(f2ev(0:mxnpjet))
      allocate(y1xx(0:mxnpjet))
      allocate(y1yy(0:mxnpjet))
      allocate(y1zz(0:mxnpjet))
      allocate(y1st(0:mxnpjet))
      allocate(y1vx(0:mxnpjet))
      allocate(y1vy(0:mxnpjet))
      allocate(y1vz(0:mxnpjet))
      allocate(y1ev(0:mxnpjet))
      allocate(y2xx(0:mxnpjet))
      allocate(y2yy(0:mxnpjet))
      allocate(y2zz(0:mxnpjet))
      allocate(y2st(0:mxnpjet))
      allocate(y2vx(0:mxnpjet))
      allocate(y2vy(0:mxnpjet))
      allocate(y2vz(0:mxnpjet))
      allocate(y2ev(0:mxnpjet))
      allocate(d3xx(0:mxnpjet),d3yy(0:mxnpjet),d3zz(0:mxnpjet))
      allocate(d3st(0:mxnpjet),d3vx(0:mxnpjet),d3vy(0:mxnpjet))
      allocate(d3vz(0:mxnpjet),d3ev(0:mxnpjet))
    end select
    lfirstsub=.false.
  endif

  
  dsqrh=dsqrt(dabs(h))
  tsqh=dsqrh**3.d0
  prefactor1=0.5d0/dsqrh

  if(.not.allocated(gaussianhistory))then
    if(systype==1)then
      call prepare_gaussian_buffer(inpjet,npjet,mxnpjet,1)
    else
      call prepare_gaussian_buffer(inpjet,npjet,mxnpjet,3)
    endif
  endif

! Reserve this timestep's slice of the sequential Gaussian pool, ahead of
! the dispatch for the same reason as in platen(). A no-op unless the pool
! is in use; the slice spans the whole jet, as in platen().
  call begin_gaussian_history_step(inpjet,npjet)

#ifdef _OPENACC
! Above its gate the run takes the device step (device_step_mod), which
! reads the pool slice just reserved; below it, and for the options that
! step does not cover, the code of the CPU build.
  if(device_step_eligible(scheme_platen))then
    call device_platen_step(timesub,h,k)
    return
  endif
#endif

! select the proper system type
  select case(systype)
    case(1)
      call smooth_charge(jetxx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      y1xx(:)=0.d0
      y1st(:)=0.d0
      y1vx(:)=0.d0
      y2xx(:)=0.d0
      y2st(:)=0.d0
      y2vx(:)=0.d0
      y2ev(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub, &
         k,fstocvx,fstocvy,fstocvz) 
        f1xx(j)=fxx
        f1st(j)=fst
        f1vx(j)=fvx
        f1ev(j)=fev
        f1stocvx(j)=fstocvx
        y1xx(ipoint) = jetxx(ipoint) + h*fxx
        y1st(ipoint) = jetst(ipoint) + h*fst
        y1vx(ipoint) = jetvx(ipoint) + h*fvx + dsqrh*fstocvx
        y1ev(ipoint) = jetve(ipoint) + h*fev
        if((y1ev(ipoint)/jetvl(ipoint))<evlim)then
          y1ev(ipoint)=jetvl(ipoint)*evlim
        endif
        y2xx(ipoint) = jetxx(ipoint) + h*fxx
        y2st(ipoint) = jetst(ipoint) + h*fst
        y2vx(ipoint) = jetvx(ipoint) + h*fvx - dsqrh*fstocvx
        y2ev(ipoint) = jetve(ipoint) + h*fev
        if((y2ev(ipoint)/jetvl(ipoint))<evlim)then
          y2ev(ipoint)=jetvl(ipoint)*evlim
        endif
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      call sum_world_darr(y1vx,npjet+1)
      call sum_world_darr(y1ev,npjet+1)
      call sum_world_darr(y2xx,npjet+1)
      call sum_world_darr(y2st,npjet+1)
      call sum_world_darr(y2vx,npjet+1)
      call sum_world_darr(y2ev,npjet+1)
      
      call smooth_charge(y1xx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y1xx, &
       y1yy,y1zz,y1ev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,jetvl, &
         y1ev,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub,k, &
         fstocvx,fstocvy,fstocvz)
        f2xx(j)=fxx
        f2st(j)=fst
        f2vx(j)=fvx
        f2ev(j)=fev
        j=j+1
      enddo
      call restore_charge()
      
      call smooth_charge(y2xx)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y2xx, &
       y2yy,y2zz,y2ev)
      j=0
      y1vx(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,jetvl, &
         y2ev,coulforce,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
         timesub,k,f3stocvx,f3stocvy,f3stocvz)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,1,1)
          u2=gaussian_history_value(k,ipoint,1,2)
        else
          u1=gaussian_buffer_value(ipoint,1,1)
          u2=gaussian_buffer_value(ipoint,1,2)
        endif
        ww(1)=(dsqrh*u1)
        zz(1)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
          
        y1vx(ipoint) = jetvx(ipoint) + f1stocvx(j)*ww(1) + &
         prefactor1*(f2vx(j)-f3vx)*zz(1) + &
         0.25d0*h*(f2vx(j)+2.d0*f1vx(j)+f3vx)
        
        j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1vx,npjet+1,jetvx)
      
      j=0
      y1xx(:)=0.d0
      y1st(:)=0.d0
      y1ev(:)=0.d0
      do ipoint=mystart,myend
	    y1xx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    y1st(ipoint) = jetst(ipoint) + h*f1st(j)
	    y1ev(ipoint) = jetve(ipoint) + h*f1ev(j)
	    if((y1ev(ipoint)/jetvl(ipoint))<evlim)then
          y1ev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      call sum_world_darr(y1ev,npjet+1)
      
      j=0
      do ipoint=mystart,myend
        call xpsys_pos_ev(ipoint,y1xx,y1yy,y1zz,y1st,jetvx,jetvy,jetvz, &
         jetvl,y1ev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2ev(j), &
         timesub+h,k)
         j=j+1
      enddo
      j=0
      y1xx(:)=0.d0
      y1ev(:)=0.d0
      do ipoint=mystart,myend
        y1xx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        y1ev(ipoint) = jetve(ipoint) + (h/2.d0)*(f1ev(j)+f2ev(j))
        if((y1ev(ipoint)/jetvl(ipoint))<evlim)then
          y1ev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1,jetxx)
      call sum_world_darr(y1ev,npjet+1,jetve)
      
      j=0
      do ipoint=mystart,myend
        call xpsys_stress_ev(ipoint,jetxx,jetyy,jetzz,y1st,jetvx,jetvy, &
         jetvz,jetvl,jetve,coulforce,f2st(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      y1st(:)=0.d0
      do ipoint=mystart,myend
	    y1st(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
	    j=j+1
      enddo
      call sum_world_darr(y1st,npjet+1,jetst)
      
      timesub=timesub+h
    case default
      call smooth_charge(jetxx,jetyy,jetzz)
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx, &
       jetyy,jetzz,jetve)
      j=0
      y1xx(:)=0.d0
	  y1yy(:)=0.d0
	  y1zz(:)=0.d0
	  y1st(:)=0.d0
	  y1vx(:)=0.d0
	  y1vy(:)=0.d0
	  y1vz(:)=0.d0
	  y1ev(:)=0.d0
	  y2xx(:)=0.d0
	  y2yy(:)=0.d0
	  y2zz(:)=0.d0
	  y2st(:)=0.d0
	  y2vx(:)=0.d0
	  y2vy(:)=0.d0
	  y2vz(:)=0.d0
	  y2ev(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub, &
         k,fstocvx,fstocvy,fstocvz)
        f1xx(j)=fxx
        f1yy(j)=fyy
        f1zz(j)=fzz
        f1st(j)=fst
        f1vx(j)=fvx
        f1vy(j)=fvy
        f1vz(j)=fvz
        f1ev(j)=fev
        f1stocvx(j)=fstocvx
        f1stocvy(j)=fstocvy
        f1stocvz(j)=fstocvz
	    y1xx(ipoint) = jetxx(ipoint) + h*fxx
	    y1yy(ipoint) = jetyy(ipoint) + h*fyy
	    y1zz(ipoint) = jetzz(ipoint) + h*fzz
	    y1st(ipoint) = jetst(ipoint) + h*fst
	    y1vx(ipoint) = jetvx(ipoint) + h*fvx + dsqrh*fstocvx
	    y1vy(ipoint) = jetvy(ipoint) + h*fvy + dsqrh*fstocvy
	    y1vz(ipoint) = jetvz(ipoint) + h*fvz + dsqrh*fstocvz
	    y1ev(ipoint) = jetve(ipoint) + h*fev
	    if((y1ev(ipoint)/jetvl(ipoint))<evlim)then
          y1ev(ipoint)=jetvl(ipoint)*evlim
        endif
	    y2xx(ipoint) = jetxx(ipoint) + h*fxx
	    y2yy(ipoint) = jetyy(ipoint) + h*fyy
	    y2zz(ipoint) = jetzz(ipoint) + h*fzz
	    y2st(ipoint) = jetst(ipoint) + h*fst
	    y2vx(ipoint) = jetvx(ipoint) + h*fvx - dsqrh*fstocvx
	    y2vy(ipoint) = jetvy(ipoint) + h*fvy - dsqrh*fstocvy
	    y2vz(ipoint) = jetvz(ipoint) + h*fvz - dsqrh*fstocvz
	    y2ev(ipoint) = jetve(ipoint) + h*fev
	    if((y2ev(ipoint)/jetvl(ipoint))<evlim)then
          y2ev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1yy,npjet+1)
      call sum_world_darr(y1zz,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      call sum_world_darr(y1vx,npjet+1)
      call sum_world_darr(y1vy,npjet+1)
      call sum_world_darr(y1vz,npjet+1)
      call sum_world_darr(y1ev,npjet+1)
      call sum_world_darr(y2xx,npjet+1)
      call sum_world_darr(y2yy,npjet+1)
      call sum_world_darr(y2zz,npjet+1)
      call sum_world_darr(y2st,npjet+1)
      call sum_world_darr(y2vx,npjet+1)
      call sum_world_darr(y2vy,npjet+1)
      call sum_world_darr(y2vz,npjet+1)
      call sum_world_darr(y2ev,npjet+1)
      
      call smooth_charge(y1xx,y1yy,y1zz)
      call compute_posnoinserted(y1xx,y1yy,y1zz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y1xx, &
       y1yy,y1zz,y1ev)
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,jetvl, &
         y1ev,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub,k, &
         fstocvx,fstocvy,fstocvz)
        f2xx(j)=fxx
        f2yy(j)=fyy
        f2zz(j)=fzz
        f2st(j)=fst
        f2vx(j)=fvx
        f2vy(j)=fvy
        f2vz(j)=fvz
        f2ev(j)=fev
	    j=j+1
      enddo
      call restore_charge()
      
      call smooth_charge(y2xx,y2yy,y2zz)
      call compute_posnoinserted(y2xx,y2yy,y2zz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y2xx, &
       y2yy,y2zz,y2ev)
      j=0
      y1vx(:)=0.d0
      y1vy(:)=0.d0
      y1vz(:)=0.d0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,jetvl, &
         y2ev,coulforce,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
         timesub,k,f3stocvx,f3stocvy,f3stocvz)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,1,1)
          u2=gaussian_history_value(k,ipoint,1,2)
        else
          u1=gaussian_buffer_value(ipoint,1,1)
          u2=gaussian_buffer_value(ipoint,1,2)
        endif
        ww(1)=(dsqrh*u1)
        zz(1)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,2,1)
          u2=gaussian_history_value(k,ipoint,2,2)
        else
          u1=gaussian_buffer_value(ipoint,2,1)
          u2=gaussian_buffer_value(ipoint,2,2)
        endif
        ww(2)=(dsqrh*u1)
        zz(2)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
        if(allocated(gaussianhistory))then
          u1=gaussian_history_value(k,ipoint,3,1)
          u2=gaussian_history_value(k,ipoint,3,2)
        else
          u1=gaussian_buffer_value(ipoint,3,1)
          u2=gaussian_buffer_value(ipoint,3,2)
        endif
        ww(3)=(dsqrh*u1)
        zz(3)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
	    
        y1vx(ipoint) = jetvx(ipoint) + f1stocvx(j)*ww(1) + &
         prefactor1*(f2vx(j)-f3vx)*zz(1) + &
	     0.25d0*h*(f2vx(j)+2.d0*f1vx(j)+f3vx)
	    
	    y1vy(ipoint) = jetvy(ipoint) + f1stocvy(j)*ww(2) + &
	     prefactor1*(f2vy(j)-f3vy)*zz(2) + &
	     0.25d0*h*(f2vy(j)+2.d0*f1vy(j)+f3vy)
	     
	    y1vz(ipoint) = jetvz(ipoint) + f1stocvz(j)*ww(3) + &
	     prefactor1*(f2vz(j)-f3vz)*zz(3) + &
	     0.25d0*h*(f2vz(j)+2.d0*f1vz(j)+f3vz)
	    
	    j=j+1
      enddo
      call restore_charge()
      call sum_world_darr(y1vx,npjet+1,jetvx)
      call sum_world_darr(y1vy,npjet+1,jetvy)
      call sum_world_darr(y1vz,npjet+1,jetvz)
      
      j=0
      y1xx(:)=0.d0
      y1yy(:)=0.d0
      y1zz(:)=0.d0
      y1st(:)=0.d0
      y1ev(:)=0.d0
      do ipoint=mystart,myend
	    y1xx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    y1yy(ipoint) = jetyy(ipoint) + h*f1yy(j)
	    y1zz(ipoint) = jetzz(ipoint) + h*f1zz(j)
	    y1st(ipoint) = jetst(ipoint) + h*f1st(j)
	    y1ev(ipoint) = jetve(ipoint) + h*f1ev(j)
	    if((y1ev(ipoint)/jetvl(ipoint))<evlim)then
          y1ev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1)
      call sum_world_darr(y1yy,npjet+1)
      call sum_world_darr(y1zz,npjet+1)
      call sum_world_darr(y1st,npjet+1)
      call sum_world_darr(y1ev,npjet+1)
      
      call compute_posnoinserted(y1xx,y1yy,y1zz)
      j=0
      do ipoint=mystart,myend
        call xpsys_pos_ev(ipoint,y1xx,y1yy,y1zz,y1st,jetvx,jetvy,jetvz,&
         jetvl,y1ev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2ev(j), &
         timesub+h,k)
         j=j+1
      enddo
      j=0
      y1xx(:)=0.d0
      y1yy(:)=0.d0
      y1zz(:)=0.d0
      y1ev(:)=0.d0
      do ipoint=mystart,myend
        y1xx(ipoint) = jetxx(ipoint) + (h/2.d0)*(f1xx(j)+f2xx(j))
        y1yy(ipoint) = jetyy(ipoint) + (h/2.d0)*(f1yy(j)+f2yy(j))
        y1zz(ipoint) = jetzz(ipoint) + (h/2.d0)*(f1zz(j)+f2zz(j))
        y1ev(ipoint) = jetve(ipoint) + (h/2.d0)*(f1ev(j)+f2ev(j))
        if((y1ev(ipoint)/jetvl(ipoint))<evlim)then
          y1ev(ipoint)=jetvl(ipoint)*evlim
        endif
	    j=j+1
      enddo
      call sum_world_darr(y1xx,npjet+1,jetxx)
      call sum_world_darr(y1yy,npjet+1,jetyy)
      call sum_world_darr(y1zz,npjet+1,jetzz)
      call sum_world_darr(y1ev,npjet+1,jetve)
      
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      j=0
      do ipoint=mystart,myend
        call xpsys_stress_ev(ipoint,jetxx,jetyy,jetzz,y1st,jetvx,jetvy,&
         jetvz,jetvl,jetve,coulforce,f2st(j),timesub+h,k)
         j=j+1
      enddo
      j=0
      y1st(:)=0.d0
      do ipoint=mystart,myend
	    y1st(ipoint) = jetst(ipoint) + (h/2.d0)*(f1st(j)+f2st(j))
	    j=j+1
      enddo
      call sum_world_darr(y1st,npjet+1,jetst)
      
      timesub=timesub+h
  end select
  
  return
      
 end subroutine platen_ev
 
 end module integrator_mod
