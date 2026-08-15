
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
                         set_chunk,set_mxchunk,idrank,mxrank
 use error_mod,         only : error,warning
 use utility_mod,       only : wiener_process1,wiener_process2, &
                         prepare_gaussian_buffer,gaussian_buffer_value, &
                         prepare_gaussian_history,gaussian_history_value, &
                         gaussianhistory,gaussianhistorysteps
 use nanojet_mod,       only : doallocate,mxnpjet,npjet,inpjet,systype,&
                         jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetst, &
                         jetms,jetch,jetvl,compute_posnoinserted, &
                         jetpt,lKVfluid,levaporation,jetve,jetce,evlim,jetfr, &
                         cp0,Bev,mev,tev, &
                         evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity, &
                         resolution, &
                         linserted,liniperturb,lairdrag,lflorentz,luppot, &
                         pfreq,consistency,findex,yieldstress,att,fve,gr, &
                         ks,li,v,velext,linserting,lremove,lmultiplestep, &
                         airdragamp,noisediff,noisefric,ldragvel,typemass, &
                         ltrackbeads,ltagbeads,lbreakup
 use dynamic_refinement_mod, only : driver_dynamic_refinement,lrefinement
 use profiling_mod, only : profiling_start,profiling_stop,prof_eom, &
                         prof_rk_update
#ifdef _OPENACC
 use accelerator_mod, only : accelerator_eom3_stage, &
                         accelerator_maxwell_evap_stage, &
                         accelerator_maxwell_rk4_stage_update, &
                         accelerator_evap_rk2_final_update, &
                         accelerator_maxwell_rk4_final_update, &
                         accelerator_maxwell_commit_state, &
                         accelerator_compute_posnoinserted_3d, &
                         accelerator_mark_device_state, &
                         accelerator_device_state_is_current, &
                         accelerator_maxwell_stress_3d, &
                         accelerator_maxwell_evap_stress_3d, &
                         rheology_maxwell, &
                         accelerator_set_persistent, &
                         accelerator_is_persistent, &
                         accelerator_set_topology_enabled, &
                         accelerator_is_topology_enabled, &
                         accelerator_rk4_final_statistics, &
                         accelerator_euler_final_statistics, &
                         accelerator_rk2_final_statistics, &
                         accelerator_platen_predict, &
                         accelerator_platen_evap_predict, &
                         accelerator_platen_velocity, &
                         accelerator_platen_evap_velocity, &
                         accelerator_platen_positions, &
                         accelerator_platen_evap_positions, &
                         accelerator_platen_stress_statistics
 use statistic_mod, only : counterlpath,ncounterlpath,maxstress, &
                         maxstressposx
#endif
 use electric_field_mod, only : nfieldtype
 use coulomb_force_mod, only : smooth_charge,restore_charge, &
                         coulforce,compute_coulomelec_driver, &
                         set_coulomb_accelerator_persistent, &
                         reset_coulomb_accelerator
 use driver_eom_mod,    only : xpsys,xpsys_pos,xpsys_stress, &
                         xpsys_KV_pos_v,xpsys_KV_st,xpsys_ev, &
                         xpsys_ev_maxwell, &
                         xpsys_pos_ev,xpsys_stress_ev
 use support_functions_mod, only : compute_geometry,compute_tangetversor, &
                          compute_curvcenter,compute_curvature, &
                          project_veltangetversor

 implicit none

 private
 
 integer, public, save :: integrator
 logical, public, save :: lintegrator=.false.
 logical, save :: persistent_reset_requested=.false.
 
 double precision, public, save :: initime = 0.d0
 double precision, public, save :: endtime = 5.d0
 logical, public, save :: lendtime

#ifdef _OPENACC
 ! Shared persistent workspace for the serial dynamic Maxwell evaporation
 ! Euler/RK2 paths and the development-oracle RK4 path.
 double precision, allocatable, save :: maxev_fxx(:,:),maxev_fyy(:,:)
 double precision, allocatable, save :: maxev_fzz(:,:),maxev_fst(:,:)
 double precision, allocatable, save :: maxev_fvx(:,:),maxev_fvy(:,:)
 double precision, allocatable, save :: maxev_fvz(:,:),maxev_fev(:,:)
 double precision, allocatable, save :: maxev_yxx(:),maxev_yyy(:),maxev_yzz(:)
 double precision, allocatable, save :: maxev_yst(:),maxev_yvx(:),maxev_yvy(:)
 double precision, allocatable, save :: maxev_yvz(:),maxev_yev(:)
 logical, save :: maxev_workspace=.false.
 logical, save :: maxev_workspace_device_mapped=.false.
 integer, save :: maxev_workspace_mxnpjet=-1
 integer, save :: maxev_workspace_mxchunk=-1
#endif
 
 public :: driver_integrator
 public :: prepare_integrator_random_history
 public :: reset_persistent_integrator

contains

 subroutine reset_persistent_integrator()
  implicit none
  persistent_reset_requested=.true.
 end subroutine reset_persistent_integrator

 subroutine prepare_integrator_random_history(h)
  implicit none
  double precision, intent(in) :: h
  integer :: nsteps,history_last
  logical :: fixed_history,dynamic_history
  if(integrator/=4 .or. systype/=4)return
  fixed_history=npjet==1000 .and. (fixed_accelerator_geometry() .or. &
   fixed_evaporative_platen_eligible())
  dynamic_history=dynamic_evaporative_platen_eligible()
  if(.not.(fixed_history .or. dynamic_history))return
  nsteps=nint((endtime-initime)/h)
  history_last=npjet
! Generate noise for the full reserved capacity in a dynamic run.  Beads
! inserted after initialization then consume the same indexed history on CPU
! and GPU without drawing random numbers inside the timestep loop.
  if(dynamic_history)history_last=mxnpjet
  call prepare_gaussian_history(inpjet,history_last,mxnpjet,3,nsteps)
#ifdef _OPENACC
!$acc enter data copyin(gaussianhistory(0:(mxnpjet+1)*6* &
!$acc& gaussianhistorysteps-1))
#endif
 end subroutine prepare_integrator_random_history

 logical function fixed_accelerator_geometry()
  implicit none
  fixed_accelerator_geometry=npjet.eq.1000 .and. mxrank.eq.1 .and. &
   mystart.eq.0 .and. myend.eq.npjet .and. linserted .and. &
   .not.linserting .and. .not.lmultiplestep .and. .not.levaporation &
   .and. lairdrag .and. .not.lflorentz .and. .not.luppot .and. &
   nfieldtype.eq.0
 end function fixed_accelerator_geometry

 logical function fixed_accelerator_eligible()
  implicit none
  fixed_accelerator_eligible=systype.eq.3 .and. fixed_accelerator_geometry()
 end function fixed_accelerator_eligible

 logical function fixed_evaporative_platen_eligible()
  implicit none
  fixed_evaporative_platen_eligible=systype.eq.4 .and. levaporation .and. &
   .not.lKVfluid .and. npjet.eq.1000 .and. mxrank.eq.1 .and. &
   mystart.eq.0 .and. myend.eq.npjet .and. linserted .and. &
   .not.linserting .and. .not.lremove .and. .not.lmultiplestep .and. &
   lairdrag .and. .not.lflorentz .and. .not.luppot .and. nfieldtype.eq.0
 end function fixed_evaporative_platen_eligible

 logical function dynamic_evaporative_platen_eligible()
  implicit none
  character(len=16) :: disable_persistent
  disable_persistent=''
  call get_environment_variable('JETSPIN_OPENACC_DISABLE_PERSISTENT', &
   disable_persistent)
  if(trim(disable_persistent)=='1')then
    dynamic_evaporative_platen_eligible=.false.
    return
  endif
  dynamic_evaporative_platen_eligible=systype.eq.4 .and. levaporation .and. &
   .not.lKVfluid .and. integrator.eq.4 .and. npjet>=100 .and. &
   mxnpjet>npjet .and. mxrank.eq.1 .and. mystart.eq.inpjet .and. &
   myend.eq.npjet .and. linserting .and. &
   .not.lmultiplestep .and. lairdrag .and. .not.lflorentz .and. &
   .not.luppot .and. nfieldtype.eq.0 .and. .not.ldragvel .and. &
   typemass.eq.0 .and. .not.ltrackbeads .and. ltagbeads .and. &
   .not.lbreakup .and. lrefinement
 end function dynamic_evaporative_platen_eligible

 logical function dynamic_rk4_accelerator_eligible()
  implicit none
  character(len=16) :: disable_persistent
  disable_persistent=''
  call get_environment_variable('JETSPIN_OPENACC_DISABLE_PERSISTENT', &
   disable_persistent)
  if(trim(disable_persistent)=='1')then
    dynamic_rk4_accelerator_eligible=.false.
    return
  endif
  dynamic_rk4_accelerator_eligible=(systype.eq.3 .and. npjet>=1000 .and. &
   mxnpjet>=1280 .and. mxrank.eq.1 .and. mystart.eq.inpjet .and. &
   myend.eq.npjet .and. linserting .and. lremove .and. &
   .not.lmultiplestep .and. &
   .not.levaporation .and. lairdrag .and. .not.lflorentz .and. &
   .not.luppot .and. nfieldtype.eq.0 .and. .not.ldragvel .and. &
   typemass.eq.0 .and. .not.ltrackbeads .and. .not.ltagbeads .and. &
   .not.lbreakup) .or. small_dynamic_test_eligible()
 end function dynamic_rk4_accelerator_eligible

 logical function small_dynamic_test_eligible()
  implicit none
  small_dynamic_test_eligible=systype.eq.3 .and. npjet>=100 .and. &
   mxnpjet>=100 .and. mxrank.eq.1 .and. mystart.eq.inpjet .and. &
   myend.eq.npjet .and. linserting .and. .not.lremove .and. &
   .not.lmultiplestep .and. .not.levaporation .and. lairdrag .and. &
   .not.lflorentz .and. .not.luppot .and. nfieldtype.eq.0 .and. &
   .not.ldragvel .and. typemass.eq.0 .and. .not.ltrackbeads .and. &
   .not.ltagbeads .and. .not.lbreakup
 end function small_dynamic_test_eligible

 logical function evaporative_dynamic_accelerator_eligible()
  implicit none
  evaporative_dynamic_accelerator_eligible=integrator>=1 .and. &
   integrator<=3 .and. systype.eq.3 .and. npjet>=inpjet .and. &
   mxnpjet>=100 .and. mxrank.eq.1 .and. mystart.eq.inpjet .and. &
   myend.eq.npjet .and. linserting .and. lremove .and. levaporation .and. &
   .not.lKVfluid .and. .not.lmultiplestep .and. lairdrag .and. &
   .not.lflorentz .and. .not.luppot .and. nfieldtype.eq.0 .and. &
   .not.ldragvel .and. typemass.eq.0 .and. .not.ltrackbeads .and. &
   .not.ltagbeads .and. .not.lbreakup
 end function evaporative_dynamic_accelerator_eligible

#if defined(_OPENACC) && defined(JETSPIN_DEV_HOST_FORCE_ORACLE)
 logical function non_evap_host_force_oracle(tstage,k,xs,ys,zs,ss, &
   vxs,vys,vzs,fx,fy,fz,fs,fvxout,fvyout,fvzout)
  implicit none
  integer, intent(in) :: k
  double precision, intent(in) :: tstage
  double precision, allocatable, intent(in) :: xs(:),ys(:),zs(:),ss(:)
  double precision, allocatable, intent(in) :: vxs(:),vys(:),vzs(:)
  double precision, intent(out) :: fx(0:),fy(0:),fz(0:),fs(0:)
  double precision, intent(out) :: fvxout(0:),fvyout(0:),fvzout(0:)
  integer :: ipoint,j
  double precision :: fstocx,fstocy,fstocz

  non_evap_host_force_oracle=.false.
  if(systype/=3 .and. systype/=4)return

  ! Coulomb has already refreshed geometry and shared bead data. Download the
  ! remaining stage state, execute the historical CPU EOM, and upload only
  ! the derivatives consumed by the device integration update.
!$acc update self(xs(0:npjet),ys(0:npjet),zs(0:npjet),ss(0:npjet), &
!$acc& vxs(0:npjet),vys(0:npjet),vzs(0:npjet),jetvl(0:npjet), &
!$acc& coulforce(0:npjet,1:3)) if_present
  j=0
  do ipoint=mystart,myend
    if(systype==4)then
      call xpsys(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,coulforce, &
       fx(j),fy(j),fz(j),fs(j),fvxout(j),fvyout(j),fvzout(j),tstage,k, &
       fstocx,fstocy,fstocz)
    else
      call xpsys(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,coulforce, &
       fx(j),fy(j),fz(j),fs(j),fvxout(j),fvyout(j),fvzout(j),tstage,k)
    endif
    j=j+1
  enddo
!$acc update device(fx(0:myend-mystart),fy(0:myend-mystart), &
!$acc& fz(0:myend-mystart),fs(0:myend-mystart), &
!$acc& fvxout(0:myend-mystart),fvyout(0:myend-mystart), &
!$acc& fvzout(0:myend-mystart)) if_present
  non_evap_host_force_oracle=.true.
 end function non_evap_host_force_oracle
#endif
  
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
  logical, save :: persistent_acc=.false.
  logical :: used_acc_eom
  
  
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

#ifdef _OPENACC
  if(.not.persistent_acc .and. fixed_accelerator_eligible())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet), &
!$acc& jetzz(0:mxnpjet),jetst(0:mxnpjet),jetvx(0:mxnpjet), &
!$acc& jetvy(0:mxnpjet),jetvz(0:mxnpjet),jetvl(0:mxnpjet), &
!$acc& jetms(0:mxnpjet),jetch(0:mxnpjet),jetfr(0:mxnpjet))
!$acc enter data create(fxx(0:mxchunk),fyy(0:mxchunk), &
!$acc& fzz(0:mxchunk),fst(0:mxchunk),fvx(0:mxchunk), &
!$acc& fvy(0:mxchunk),fvz(0:mxchunk))
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(.true.)
    persistent_acc=.true.
  endif
#endif
  
  
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
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub,k,jetxx,jetyy,jetzz, &
       jetst,jetvx,jetvy,jetvz,fxx,fyy,fzz,fst,fvx,fvy,fvz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,jetxx, &
       jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetvl,coulforce,jetms, &
       jetch,jetfr,fxx,fyy,fzz,fst,fvx,fvy,fvz,linserted, &
       liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
       consistency,findex,yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
        do ipoint=mystart,myend
          call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
            jetvl,coulforce,fxx(j),fyy(j),fzz(j),fst(j), &
            fvx(j),fvy(j),fvz(j),timesub,k)
          j=j+1
        enddo
      endif
      call profiling_stop(prof_eom)
      call restore_charge()
      j=0
      call profiling_start(prof_rk_update)
      if(persistent_acc)then
#ifdef _OPENACC
        call accelerator_euler_final_statistics(mystart,myend,h, &
         jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         fxx,fyy,fzz,fst,fvx,fvy,fvz,counterlpath,ncounterlpath, &
         maxstress,maxstressposx)
#endif
      else
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
      endif
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
  logical, save :: persistent_acc=.false.
  logical :: used_acc_eom
  
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

#ifdef _OPENACC
  if(.not.persistent_acc .and. fixed_accelerator_eligible())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet), &
!$acc& jetzz(0:mxnpjet),jetst(0:mxnpjet),jetvx(0:mxnpjet), &
!$acc& jetvy(0:mxnpjet),jetvz(0:mxnpjet),jetvl(0:mxnpjet), &
!$acc& jetms(0:mxnpjet),jetch(0:mxnpjet),jetfr(0:mxnpjet))
!$acc enter data create(yxx(0:mxnpjet),yyy(0:mxnpjet), &
!$acc& yzz(0:mxnpjet),yst(0:mxnpjet),yvx(0:mxnpjet), &
!$acc& yvy(0:mxnpjet),yvz(0:mxnpjet))
!$acc enter data create(f1xx(0:mxchunk),f1yy(0:mxchunk), &
!$acc& f1zz(0:mxchunk),f1st(0:mxchunk),f1vx(0:mxchunk), &
!$acc& f1vy(0:mxchunk),f1vz(0:mxchunk),f2xx(0:mxchunk), &
!$acc& f2yy(0:mxchunk),f2zz(0:mxchunk),f2st(0:mxchunk), &
!$acc& f2vx(0:mxchunk),f2vy(0:mxchunk),f2vz(0:mxchunk))
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(.true.)
    persistent_acc=.true.
  endif
#endif
  
  
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
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub,k,jetxx,jetyy,jetzz, &
       jetst,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,jetxx, &
       jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetvl,coulforce,jetms, &
       jetch,jetfr,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,linserted, &
       liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
       consistency,findex,yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
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
      endif
      call profiling_stop(prof_eom)
      j=0
      call profiling_start(prof_rk_update)
      if(persistent_acc)then
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
!$acc& yxx,yyy,yzz,yst,yvx,yvy,yvz) private(j)
        do ipoint=mystart,myend
          j=ipoint-mystart
          yxx(ipoint)=jetxx(ipoint)+h*f1xx(j)
          yyy(ipoint)=jetyy(ipoint)+h*f1yy(j)
          yzz(ipoint)=jetzz(ipoint)+h*f1zz(j)
          yst(ipoint)=jetst(ipoint)+h*f1st(j)
          yvx(ipoint)=jetvx(ipoint)+h*f1vx(j)
          yvy(ipoint)=jetvy(ipoint)+h*f1vy(j)
          yvz(ipoint)=jetvz(ipoint)+h*f1vz(j)
        enddo
!$acc end parallel loop
#endif
      else
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
      endif
      call profiling_stop(prof_rk_update)
      call restore_charge()
      
      if(.not.persistent_acc)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
      endif
      
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub+h,k,yxx,yyy,yzz, &
       yst,yvx,yvy,yvz,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,yxx,yyy, &
       yzz,yst,yvx,yvy,yvz,jetvl,coulforce,jetms,jetch,jetfr, &
       f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,linserted,liniperturb, &
       lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
       yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
        do ipoint=mystart,myend
          call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
           jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
           f2vx(j),f2vy(j),f2vz(j),timesub+h,k)
          j=j+1
        enddo
      endif
      call profiling_stop(prof_eom)
      j=0
      call profiling_start(prof_rk_update)
      if(persistent_acc)then
#ifdef _OPENACC
        call accelerator_rk2_final_statistics(mystart,myend,h, &
         jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
         f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
         counterlpath,ncounterlpath,maxstress,maxstressposx)
#endif
      else
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
      endif
      call profiling_stop(prof_rk_update)
      call restore_charge()
      if(.not.persistent_acc)then
        call sum_world_darr(yxx,npjet+1,jetxx)
        call sum_world_darr(yyy,npjet+1,jetyy)
        call sum_world_darr(yzz,npjet+1,jetzz)
        call sum_world_darr(yst,npjet+1,jetst)
        call sum_world_darr(yvx,npjet+1,jetvx)
        call sum_world_darr(yvy,npjet+1,jetvy)
        call sum_world_darr(yvz,npjet+1,jetvz)
      endif
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
  logical, save :: persistent_acc=.false.
  logical :: used_acc_eom
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz

#ifdef _OPENACC
! A refinement-driven capacity increase occurs before this integrator call,
! unlike nozzle insertion which requests the reset from the main loop. In
! either case, detach the old scratch arrays before reallocating them.
  if((persistent_reset_requested .or. doallocate) .and. persistent_acc)then
!$acc exit data delete(yxx,yyy,yzz,yst,yvx,yvy,yvz, &
!$acc& f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f2xx,f2yy,f2zz,f2st, &
!$acc& f2vx,f2vy,f2vz,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz, &
!$acc& f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz)
    call reset_coulomb_accelerator()
    call accelerator_set_persistent(.false.)
    call set_coulomb_accelerator_persistent(.false.)
    persistent_acc=.false.
  endif
  persistent_reset_requested=.false.
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

#ifdef _OPENACC
! RK4 supports both the fixed 1,000-bead benchmark and the bounded dynamic
! topology benchmark.  The latter reserves enough capacity before mapping.
  if(.not.persistent_acc .and. (fixed_accelerator_eligible() .or. &
   dynamic_rk4_accelerator_eligible()))then
    ! A capacity rebind performed by the topology driver already owns the
    ! new jet mapping.  Recreate only the RK workspace in that case; entering
    ! the jet arrays again would leave a second OpenACC present reference and
    ! make a later reallocation fail with a partially-present mapping.
    if(.not.accelerator_is_topology_enabled())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet), &
!$acc& jetzz(0:mxnpjet),jetst(0:mxnpjet),jetvx(0:mxnpjet), &
!$acc& jetvy(0:mxnpjet),jetvz(0:mxnpjet),jetvl(0:mxnpjet), &
!$acc& jetms(0:mxnpjet),jetch(0:mxnpjet),jetfr(0:mxnpjet))
    endif
!$acc enter data create(yxx(0:mxnpjet),yyy(0:mxnpjet), &
!$acc& yzz(0:mxnpjet),yst(0:mxnpjet),yvx(0:mxnpjet), &
!$acc& yvy(0:mxnpjet),yvz(0:mxnpjet))
!$acc enter data create(f1xx(0:mxchunk),f1yy(0:mxchunk), &
!$acc& f1zz(0:mxchunk),f1st(0:mxchunk),f1vx(0:mxchunk), &
!$acc& f1vy(0:mxchunk),f1vz(0:mxchunk),f2xx(0:mxchunk), &
!$acc& f2yy(0:mxchunk),f2zz(0:mxchunk),f2st(0:mxchunk), &
!$acc& f2vx(0:mxchunk),f2vy(0:mxchunk),f2vz(0:mxchunk))
!$acc enter data create(f3xx(0:mxchunk),f3yy(0:mxchunk), &
!$acc& f3zz(0:mxchunk),f3st(0:mxchunk),f3vx(0:mxchunk), &
!$acc& f3vy(0:mxchunk),f3vz(0:mxchunk),f4xx(0:mxchunk), &
!$acc& f4yy(0:mxchunk),f4zz(0:mxchunk),f4st(0:mxchunk), &
!$acc& f4vx(0:mxchunk),f4vy(0:mxchunk),f4vz(0:mxchunk))
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(.true.)
    persistent_acc=.true.
  endif
#endif
  
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
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub,k,jetxx,jetyy,jetzz, &
       jetst,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,jetxx, &
       jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetvl,coulforce,jetms, &
       jetch,jetfr,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,linserted, &
       liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
       consistency,findex,yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
        do ipoint=mystart,myend
          call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
           jetvl,coulforce,f1xx(j),f1yy(j),f1zz(j),f1st(j),f1vx(j), &
           f1vy(j),f1vz(j),timesub,k)
          j=j+1
        enddo
      endif
      call profiling_stop(prof_eom)
      j=0
      call profiling_start(prof_rk_update)
      if(persistent_acc)then
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
!$acc& yxx,yyy,yzz,yst,yvx,yvy,yvz) private(j)
        do ipoint=mystart,myend
          j=ipoint-mystart
          yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
          yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f1yy(j)
          yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f1zz(j)
          yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
          yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
          yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f1vy(j)
          yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f1vz(j)
        enddo
!$acc end parallel loop
#endif
      else
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
      endif
      call profiling_stop(prof_rk_update)
      call restore_charge()
      if(.not.persistent_acc)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
      endif
!     2°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub+0.5d0*h,k,yxx,yyy,yzz, &
       yst,yvx,yvy,yvz,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,yxx,yyy, &
       yzz,yst,yvx,yvy,yvz,jetvl,coulforce,jetms,jetch,jetfr,f2xx, &
       f2yy,f2zz,f2st,f2vx,f2vy,f2vz,linserted,liniperturb,lairdrag, &
       lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
       att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
        do ipoint=mystart,myend
          call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
           jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
           f2vx(j),f2vy(j),f2vz(j),timesub+h/2.d0,k)
          j=j+1
        enddo
      endif
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
      if(persistent_acc)then
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
!$acc& yxx,yyy,yzz,yst,yvx,yvy,yvz) private(j)
        do ipoint=mystart,myend
          j=ipoint-mystart
          yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f2xx(j)
          yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f2yy(j)
          yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f2zz(j)
          yst(ipoint) = jetst(ipoint) + 0.5d0*h*f2st(j)
          yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(j)
          yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f2vy(j)
          yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f2vz(j)
        enddo
!$acc end parallel loop
#endif
      else
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
      endif
      call profiling_stop(prof_rk_update)
      call restore_charge()
      if(.not.persistent_acc)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
      endif
!     3°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub+0.5d0*h,k,yxx,yyy,yzz, &
       yst,yvx,yvy,yvz,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,yxx,yyy, &
       yzz,yst,yvx,yvy,yvz,jetvl,coulforce,jetms,jetch,jetfr,f3xx, &
       f3yy,f3zz,f3st,f3vx,f3vy,f3vz,linserted,liniperturb,lairdrag, &
       lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
       att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
        do ipoint=mystart,myend
          call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
           jetvl,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
           f3vx(j),f3vy(j),f3vz(j),timesub+h/2.d0,k)
          j=j+1
        enddo
      endif
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
      if(persistent_acc)then
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz, &
!$acc& yxx,yyy,yzz,yst,yvx,yvy,yvz) private(j)
        do ipoint=mystart,myend
          j=ipoint-mystart
          yxx(ipoint) = jetxx(ipoint) + h*f3xx(j)
          yyy(ipoint) = jetyy(ipoint) + h*f3yy(j)
          yzz(ipoint) = jetzz(ipoint) + h*f3zz(j)
          yst(ipoint) = jetst(ipoint) + h*f3st(j)
          yvx(ipoint) = jetvx(ipoint) + h*f3vx(j)
          yvy(ipoint) = jetvy(ipoint) + h*f3vy(j)
          yvz(ipoint) = jetvz(ipoint) + h*f3vz(j)
        enddo
!$acc end parallel loop
#endif
      else
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
      endif
      call profiling_stop(prof_rk_update)
      call restore_charge()
      if(.not.persistent_acc)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
      endif
!     4°step
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      call profiling_start(prof_eom)
      used_acc_eom=.false.
#ifdef _OPENACC
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
      used_acc_eom=non_evap_host_force_oracle(timesub+h,k,yxx,yyy,yzz, &
       yst,yvx,yvy,yvz,f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz)
#else
      if(systype.eq.3) used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,yxx,yyy, &
       yzz,yst,yvx,yvy,yvz,jetvl,coulforce,jetms,jetch,jetfr,f4xx, &
       f4yy,f4zz,f4st,f4vx,f4vy,f4vz,linserted,liniperturb,lairdrag, &
       lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
       att,fve,gr,ks,li,v,velext,.false.,0.d0)
#endif
#endif
      if(.not.used_acc_eom)then
        do ipoint=mystart,myend
          call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
           jetvl,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
           f4vx(j),f4vy(j),f4vz(j),timesub+h,k)
          j=j+1
        enddo
      endif
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
      if(persistent_acc)then
#ifdef _OPENACC
        call accelerator_rk4_final_statistics(mystart,myend,h, &
         jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
         f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
         f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz, &
         f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz, &
         counterlpath,ncounterlpath,maxstress,maxstressposx)
#endif
      else
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
      endif
      call profiling_stop(prof_rk_update)
      call restore_charge()
      if(.not.persistent_acc)then
        call sum_world_darr(yxx,npjet+1,jetxx)
        call sum_world_darr(yyy,npjet+1,jetyy)
        call sum_world_darr(yzz,npjet+1,jetzz)
        call sum_world_darr(yst,npjet+1,jetst)
        call sum_world_darr(yvx,npjet+1,jetvx)
        call sum_world_darr(yvy,npjet+1,jetvy)
        call sum_world_darr(yvz,npjet+1,jetvz)
      endif
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
  logical, save :: persistent_acc=.false.
  logical :: used_acc_eom
  
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

#ifdef _OPENACC
  if(.not.persistent_acc .and. systype==4 .and. &
   fixed_accelerator_geometry() .and. allocated(gaussianhistory))then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet), &
!$acc& jetzz(0:mxnpjet),jetst(0:mxnpjet),jetvx(0:mxnpjet), &
!$acc& jetvy(0:mxnpjet),jetvz(0:mxnpjet),jetvl(0:mxnpjet), &
!$acc& jetms(0:mxnpjet),jetch(0:mxnpjet),jetfr(0:mxnpjet))
!$acc enter data create(f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,d3xx,d3yy,d3zz,d3st, &
!$acc& d3vx,d3vy,d3vz,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz, &
!$acc& y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz)
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(.true.)
    persistent_acc=.true.
  endif
#endif
  
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
#ifdef _OPENACC
      if(persistent_acc)then
        call smooth_charge(jetxx,jetyy,jetzz)
        call compute_coulomelec_driver(k,timesub,coulforce,jetvl,jetxx,jetyy,jetzz)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        used_acc_eom=non_evap_host_force_oracle(timesub,k,jetxx,jetyy,jetzz, &
         jetst,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz)
!$acc update device(f1xx(0:myend-mystart),f1yy(0:myend-mystart), &
!$acc& f1zz(0:myend-mystart),f1st(0:myend-mystart), &
!$acc& f1vx(0:myend-mystart),f1vy(0:myend-mystart),f1vz(0:myend-mystart))
#else
        used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,jetxx,jetyy, &
         jetzz,jetst,jetvx,jetvy,jetvz,jetvl,coulforce,jetms,jetch,jetfr, &
         f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.true.,noisefric)
#endif
        call accelerator_platen_predict(mystart,myend,h,airdragamp(1), &
         noisediff,jetms,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
         f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,y1xx,y1yy,y1zz,y1st, &
         y1vx,y1vy,y1vz,y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz)
        call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y1xx,y1yy,y1zz)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        used_acc_eom=non_evap_host_force_oracle(timesub,k,y1xx,y1yy,y1zz, &
         y1st,y1vx,y1vy,y1vz,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz)
!$acc update device(f2xx(0:myend-mystart),f2yy(0:myend-mystart), &
!$acc& f2zz(0:myend-mystart),f2st(0:myend-mystart), &
!$acc& f2vx(0:myend-mystart),f2vy(0:myend-mystart),f2vz(0:myend-mystart))
#else
        used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,y1xx,y1yy, &
         y1zz,y1st,y1vx,y1vy,y1vz,jetvl,coulforce,jetms,jetch,jetfr, &
         f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.true.,noisefric)
#endif
        call compute_coulomelec_driver(k,timesub,coulforce,jetvl,y2xx,y2yy,y2zz)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        used_acc_eom=non_evap_host_force_oracle(timesub,k,y2xx,y2yy,y2zz, &
         y2st,y2vx,y2vy,y2vz,d3xx,d3yy,d3zz,d3st,d3vx,d3vy,d3vz)
!$acc update device(d3xx(0:myend-mystart),d3yy(0:myend-mystart), &
!$acc& d3zz(0:myend-mystart),d3st(0:myend-mystart), &
!$acc& d3vx(0:myend-mystart),d3vy(0:myend-mystart),d3vz(0:myend-mystart))
#else
        used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,y2xx,y2yy, &
         y2zz,y2st,y2vx,y2vy,y2vz,jetvl,coulforce,jetms,jetch,jetfr, &
         d3xx,d3yy,d3zz,d3st,d3vx,d3vy,d3vz,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.true.,noisefric)
#endif
        call accelerator_platen_velocity(mystart,myend,mxnpjet, &
         gaussianhistorysteps,k,h, &
         airdragamp(1),noisediff,jetms,gaussianhistory,jetvx,jetvy,jetvz, &
         f1vx,f1vy,f1vz,f2vx,f2vy,f2vz,d3vx,d3vy,d3vz)
        call accelerator_platen_positions(mystart,myend,npjet,h,pfreq, &
         liniperturb,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        used_acc_eom=non_evap_host_force_oracle(timesub+h,k,jetxx,jetyy,jetzz, &
         y1st,jetvx,jetvy,jetvz,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz)
!$acc update device(f2xx(0:myend-mystart),f2yy(0:myend-mystart), &
!$acc& f2zz(0:myend-mystart),f2st(0:myend-mystart), &
!$acc& f2vx(0:myend-mystart),f2vy(0:myend-mystart),f2vz(0:myend-mystart))
#else
        used_acc_eom=accelerator_eom3_stage(mystart,myend,npjet,jetxx,jetyy, &
         jetzz,y1st,jetvx,jetvy,jetvz,jetvl,coulforce,jetms,jetch,jetfr, &
         f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.true.,noisefric)
#endif
        call accelerator_platen_stress_statistics(mystart,myend,h,jetxx, &
         jetyy,jetzz,jetst,f1st,f2st,counterlpath,ncounterlpath,maxstress, &
         maxstressposx)
        call restore_charge()
        timesub=timesub+h
        return
      endif
#endif
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
  logical, save :: persistent_acc=.false.
  
  double precision ::  fxx
  double precision ::  fyy
  double precision ::  fzz
  double precision ::  fst
  double precision ::  fvx
  double precision ::  fvy
  double precision ::  fvz
  
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

#ifdef _OPENACC
 subroutine ensure_maxwell_evap_device_workspace()
  implicit none
  logical :: rebuild

  rebuild=.not.maxev_workspace .or. maxev_workspace_mxnpjet<mxnpjet .or. &
   maxev_workspace_mxchunk<mxchunk .or. persistent_reset_requested
  if(.not.rebuild)then
    call accelerator_set_topology_enabled(.true.)
    call accelerator_set_persistent(.true.)
    call set_coulomb_accelerator_persistent(.true.)
    return
  endif

  if(maxev_workspace)then
    if(maxev_workspace_device_mapped)then
!$acc exit data delete(maxev_fxx,maxev_fyy,maxev_fzz,maxev_fst, &
!$acc& maxev_fvx,maxev_fvy,maxev_fvz,maxev_fev,maxev_yxx,maxev_yyy, &
!$acc& maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev)
      maxev_workspace_device_mapped=.false.
    endif
    deallocate(maxev_fxx,maxev_fyy,maxev_fzz,maxev_fst)
    deallocate(maxev_fvx,maxev_fvy,maxev_fvz,maxev_fev)
    deallocate(maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst)
    deallocate(maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev)
    call reset_coulomb_accelerator(coulforce)
  endif

  allocate(maxev_fxx(0:mxchunk,4),maxev_fyy(0:mxchunk,4))
  allocate(maxev_fzz(0:mxchunk,4),maxev_fst(0:mxchunk,4))
  allocate(maxev_fvx(0:mxchunk,4),maxev_fvy(0:mxchunk,4))
  allocate(maxev_fvz(0:mxchunk,4),maxev_fev(0:mxchunk,4))
  allocate(maxev_yxx(0:mxnpjet),maxev_yyy(0:mxnpjet))
  allocate(maxev_yzz(0:mxnpjet),maxev_yst(0:mxnpjet))
  allocate(maxev_yvx(0:mxnpjet),maxev_yvy(0:mxnpjet))
  allocate(maxev_yvz(0:mxnpjet),maxev_yev(0:mxnpjet))

  if(.not.accelerator_is_topology_enabled())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet),jetzz(0:mxnpjet), &
!$acc& jetst(0:mxnpjet),jetvx(0:mxnpjet),jetvy(0:mxnpjet),jetvz(0:mxnpjet), &
!$acc& jetvl(0:mxnpjet),jetve(0:mxnpjet),jetce(0:mxnpjet),jetms(0:mxnpjet), &
!$acc& jetch(0:mxnpjet),jetfr(0:mxnpjet))
  endif
!$acc enter data create(maxev_fxx,maxev_fyy,maxev_fzz,maxev_fst, &
!$acc& maxev_fvx,maxev_fvy,maxev_fvz,maxev_fev,maxev_yxx,maxev_yyy, &
!$acc& maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev)

  maxev_workspace=.true.
  maxev_workspace_device_mapped=.true.
  maxev_workspace_mxnpjet=mxnpjet
  maxev_workspace_mxchunk=mxchunk
  persistent_reset_requested=.false.
  call accelerator_set_topology_enabled(.true.)
  call accelerator_set_persistent(.true.)
  call set_coulomb_accelerator_persistent(.true.)
 end subroutine ensure_maxwell_evap_device_workspace

 subroutine maxwell_evap_device_stage(tstage,k,xs,ys,zs,ss,vxs,vys,vzs,ves, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,stochastic_model)
  implicit none
  integer, intent(in) :: k
  double precision, intent(in) :: tstage
  double precision, allocatable, intent(inout) :: xs(:),ys(:),zs(:)
  double precision, allocatable, intent(inout) :: ss(:),vxs(:),vys(:),vzs(:)
  double precision, allocatable, intent(inout) :: ves(:)
  double precision, intent(inout) :: fxx(0:),fyy(0:),fzz(0:),fst(0:)
  double precision, intent(inout) :: fvx(0:),fvy(0:),fvz(0:),fev(0:)
  logical, intent(in), optional :: stochastic_model
  integer :: ipoint,j,nactive
  logical :: stochastic_stage
  double precision :: fstocx,fstocy,fstocz

  nactive=myend-mystart
  stochastic_stage=.false.
  if(present(stochastic_model))stochastic_stage=stochastic_model
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
  ! Development oracle: evaluate the complete trusted Maxwell/Yarin EOM on
  ! the host and upload only its derivatives. State updates remain on device.
!$acc update self(xs(0:npjet),ys(0:npjet),zs(0:npjet),ss(0:npjet), &
!$acc& vxs(0:npjet),vys(0:npjet),vzs(0:npjet),ves(0:npjet), &
!$acc& jetvl(0:npjet),jetms(0:npjet),jetch(0:npjet),jetfr(0:npjet)) if_present
  call smooth_charge(xs,ys,zs)
  call compute_posnoinserted(xs,ys,zs)
  call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
  j=0
  do ipoint=mystart,myend
    if(stochastic_stage)then
      call xpsys_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves,coulforce, &
       fxx(j),fyy(j),fzz(j),fst(j),fvx(j),fvy(j),fvz(j),fev(j),tstage,k, &
       fstocx,fstocy,fstocz)
    else
      call xpsys_ev_maxwell(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
       coulforce,fxx(j),fyy(j),fzz(j),fst(j),fvx(j),fvy(j),fvz(j),fev(j), &
       tstage,k)
    endif
    j=j+1
  enddo
  call restore_charge()
!$acc update device(fxx(0:nactive),fyy(0:nactive),fzz(0:nactive), &
!$acc& fst(0:nactive),fvx(0:nactive),fvy(0:nactive), &
!$acc& fvz(0:nactive),fev(0:nactive)) if_present
#else
  call smooth_charge(xs,ys,zs)
  call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution,xs,ys,zs)
  call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
  call accelerator_maxwell_evap_stage(mystart,myend,npjet,xs,ys,zs,ss, &
   vxs,vys,vzs,jetvl,ves,coulforce,jetms,jetch,jetfr, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,linserting,linserted,liniperturb, &
   lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
   yieldstress,att,fve,gr,ks,li,v,velext,stochastic_stage,noisefric,evairv, &
   evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev)
  call restore_charge()
#endif
 end subroutine maxwell_evap_device_stage

 subroutine finish_maxwell_evap_device_step(timesub,h)
  implicit none
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h

  call accelerator_maxwell_commit_state(mystart,myend, &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst, &
   maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve, &
   counterlpath,ncounterlpath,maxstress,maxstressposx)
  call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
   jetxx,jetyy,jetzz)
  call accelerator_mark_device_state(.true.)
  timesub=timesub+h
 end subroutine finish_maxwell_evap_device_step

 subroutine eulsys_maxwell_ev_device(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h

  call ensure_maxwell_evap_device_workspace()
  call maxwell_evap_device_stage(timesub,k,jetxx,jetyy,jetzz,jetst, &
   jetvx,jetvy,jetvz,jetve,maxev_fxx(:,1),maxev_fyy(:,1), &
   maxev_fzz(:,1),maxev_fst(:,1),maxev_fvx(:,1),maxev_fvy(:,1), &
   maxev_fvz(:,1),maxev_fev(:,1))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
  ! NVHPC does not reliably resolve an update through the assumed-shape
  ! column aliases used by maxwell_evap_device_stage.  Name the persistent
  ! workspace columns explicitly in the development oracle.
!$acc update device(maxev_fxx(0:myend-mystart,1), &
!$acc& maxev_fyy(0:myend-mystart,1),maxev_fzz(0:myend-mystart,1), &
!$acc& maxev_fst(0:myend-mystart,1),maxev_fvx(0:myend-mystart,1), &
!$acc& maxev_fvy(0:myend-mystart,1),maxev_fvz(0:myend-mystart,1), &
!$acc& maxev_fev(0:myend-mystart,1))
#endif
  call accelerator_maxwell_rk4_stage_update(mystart,myend,h,3, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,1),maxev_fyy(:,1),maxev_fzz(:,1),maxev_fst(:,1), &
   maxev_fvx(:,1),maxev_fvy(:,1),maxev_fvz(:,1),maxev_fev(:,1), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)
  call finish_maxwell_evap_device_step(timesub,h)
 end subroutine eulsys_maxwell_ev_device

 subroutine rk2sys_maxwell_ev_device(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h

  call ensure_maxwell_evap_device_workspace()
  call maxwell_evap_device_stage(timesub,k,jetxx,jetyy,jetzz,jetst, &
   jetvx,jetvy,jetvz,jetve,maxev_fxx(:,1),maxev_fyy(:,1), &
   maxev_fzz(:,1),maxev_fst(:,1),maxev_fvx(:,1),maxev_fvy(:,1), &
   maxev_fvz(:,1),maxev_fev(:,1))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(maxev_fxx(0:myend-mystart,1), &
!$acc& maxev_fyy(0:myend-mystart,1),maxev_fzz(0:myend-mystart,1), &
!$acc& maxev_fst(0:myend-mystart,1),maxev_fvx(0:myend-mystart,1), &
!$acc& maxev_fvy(0:myend-mystart,1),maxev_fvz(0:myend-mystart,1), &
!$acc& maxev_fev(0:myend-mystart,1))
#endif
  call accelerator_maxwell_rk4_stage_update(mystart,myend,h,3, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,1),maxev_fyy(:,1),maxev_fzz(:,1),maxev_fst(:,1), &
   maxev_fvx(:,1),maxev_fvy(:,1),maxev_fvz(:,1),maxev_fev(:,1), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)

  call maxwell_evap_device_stage(timesub+h,k,maxev_yxx,maxev_yyy,maxev_yzz, &
   maxev_yst,maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev, &
   maxev_fxx(:,2),maxev_fyy(:,2),maxev_fzz(:,2),maxev_fst(:,2), &
   maxev_fvx(:,2),maxev_fvy(:,2),maxev_fvz(:,2),maxev_fev(:,2))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(maxev_fxx(0:myend-mystart,2), &
!$acc& maxev_fyy(0:myend-mystart,2),maxev_fzz(0:myend-mystart,2), &
!$acc& maxev_fst(0:myend-mystart,2),maxev_fvx(0:myend-mystart,2), &
!$acc& maxev_fvy(0:myend-mystart,2),maxev_fvz(0:myend-mystart,2), &
!$acc& maxev_fev(0:myend-mystart,2))
#endif
  call accelerator_evap_rk2_final_update(mystart,myend,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,1),maxev_fyy(:,1),maxev_fzz(:,1),maxev_fst(:,1), &
   maxev_fvx(:,1),maxev_fvy(:,1),maxev_fvz(:,1),maxev_fev(:,1), &
   maxev_fxx(:,2),maxev_fyy(:,2),maxev_fzz(:,2),maxev_fst(:,2), &
   maxev_fvx(:,2),maxev_fvy(:,2),maxev_fvz(:,2),maxev_fev(:,2), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)
  call finish_maxwell_evap_device_step(timesub,h)
 end subroutine rk2sys_maxwell_ev_device

 subroutine rk4sys_maxwell_ev_device(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h

  call ensure_maxwell_evap_device_workspace()

  call maxwell_evap_device_stage(timesub,k,jetxx,jetyy,jetzz,jetst, &
   jetvx,jetvy,jetvz,jetve,maxev_fxx(:,1),maxev_fyy(:,1), &
   maxev_fzz(:,1),maxev_fst(:,1),maxev_fvx(:,1),maxev_fvy(:,1), &
   maxev_fvz(:,1),maxev_fev(:,1))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(maxev_fxx(0:myend-mystart,1), &
!$acc& maxev_fyy(0:myend-mystart,1),maxev_fzz(0:myend-mystart,1), &
!$acc& maxev_fst(0:myend-mystart,1),maxev_fvx(0:myend-mystart,1), &
!$acc& maxev_fvy(0:myend-mystart,1),maxev_fvz(0:myend-mystart,1), &
!$acc& maxev_fev(0:myend-mystart,1))
#endif
  call accelerator_maxwell_rk4_stage_update(mystart,myend,h,1, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,1),maxev_fyy(:,1),maxev_fzz(:,1),maxev_fst(:,1), &
   maxev_fvx(:,1),maxev_fvy(:,1),maxev_fvz(:,1),maxev_fev(:,1), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)

  call maxwell_evap_device_stage(timesub+0.5d0*h,k,maxev_yxx,maxev_yyy, &
   maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev, &
   maxev_fxx(:,2),maxev_fyy(:,2),maxev_fzz(:,2),maxev_fst(:,2), &
   maxev_fvx(:,2),maxev_fvy(:,2),maxev_fvz(:,2),maxev_fev(:,2))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(maxev_fxx(0:myend-mystart,2), &
!$acc& maxev_fyy(0:myend-mystart,2),maxev_fzz(0:myend-mystart,2), &
!$acc& maxev_fst(0:myend-mystart,2),maxev_fvx(0:myend-mystart,2), &
!$acc& maxev_fvy(0:myend-mystart,2),maxev_fvz(0:myend-mystart,2), &
!$acc& maxev_fev(0:myend-mystart,2))
#endif
  call accelerator_maxwell_rk4_stage_update(mystart,myend,h,2, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,2),maxev_fyy(:,2),maxev_fzz(:,2),maxev_fst(:,2), &
   maxev_fvx(:,2),maxev_fvy(:,2),maxev_fvz(:,2),maxev_fev(:,2), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)

  call maxwell_evap_device_stage(timesub+0.5d0*h,k,maxev_yxx,maxev_yyy, &
   maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev, &
   maxev_fxx(:,3),maxev_fyy(:,3),maxev_fzz(:,3),maxev_fst(:,3), &
   maxev_fvx(:,3),maxev_fvy(:,3),maxev_fvz(:,3),maxev_fev(:,3))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(maxev_fxx(0:myend-mystart,3), &
!$acc& maxev_fyy(0:myend-mystart,3),maxev_fzz(0:myend-mystart,3), &
!$acc& maxev_fst(0:myend-mystart,3),maxev_fvx(0:myend-mystart,3), &
!$acc& maxev_fvy(0:myend-mystart,3),maxev_fvz(0:myend-mystart,3), &
!$acc& maxev_fev(0:myend-mystart,3))
#endif
  call accelerator_maxwell_rk4_stage_update(mystart,myend,h,3, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,3),maxev_fyy(:,3),maxev_fzz(:,3),maxev_fst(:,3), &
   maxev_fvx(:,3),maxev_fvy(:,3),maxev_fvz(:,3),maxev_fev(:,3), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)

  call maxwell_evap_device_stage(timesub+h,k,maxev_yxx,maxev_yyy,maxev_yzz, &
   maxev_yst,maxev_yvx,maxev_yvy,maxev_yvz,maxev_yev, &
   maxev_fxx(:,4),maxev_fyy(:,4),maxev_fzz(:,4),maxev_fst(:,4), &
   maxev_fvx(:,4),maxev_fvy(:,4),maxev_fvz(:,4),maxev_fev(:,4))
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(maxev_fxx(0:myend-mystart,4), &
!$acc& maxev_fyy(0:myend-mystart,4),maxev_fzz(0:myend-mystart,4), &
!$acc& maxev_fst(0:myend-mystart,4),maxev_fvx(0:myend-mystart,4), &
!$acc& maxev_fvy(0:myend-mystart,4),maxev_fvz(0:myend-mystart,4), &
!$acc& maxev_fev(0:myend-mystart,4))
#endif
  call accelerator_maxwell_rk4_final_update(mystart,myend,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   maxev_fxx(:,1),maxev_fyy(:,1),maxev_fzz(:,1),maxev_fst(:,1), &
   maxev_fvx(:,1),maxev_fvy(:,1),maxev_fvz(:,1),maxev_fev(:,1), &
   maxev_fxx(:,2),maxev_fyy(:,2),maxev_fzz(:,2),maxev_fst(:,2), &
   maxev_fvx(:,2),maxev_fvy(:,2),maxev_fvz(:,2),maxev_fev(:,2), &
   maxev_fxx(:,3),maxev_fyy(:,3),maxev_fzz(:,3),maxev_fst(:,3), &
   maxev_fvx(:,3),maxev_fvy(:,3),maxev_fvz(:,3),maxev_fev(:,3), &
   maxev_fxx(:,4),maxev_fyy(:,4),maxev_fzz(:,4),maxev_fst(:,4), &
   maxev_fvx(:,4),maxev_fvy(:,4),maxev_fvz(:,4),maxev_fev(:,4), &
   maxev_yxx,maxev_yyy,maxev_yzz,maxev_yst,maxev_yvx,maxev_yvy, &
   maxev_yvz,maxev_yev,evlim)
  call finish_maxwell_evap_device_step(timesub,h)
 end subroutine rk4sys_maxwell_ev_device
#endif
 
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
  if(evaporative_dynamic_accelerator_eligible())then
    call eulsys_maxwell_ev_device(timesub,h,k)
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
  if(evaporative_dynamic_accelerator_eligible())then
    call rk2sys_maxwell_ev_device(timesub,h,k)
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
  logical, save :: persistent_acc=.false.
  logical, save :: workspace_device_mapped=.false.
  
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
  logical :: used_acc_maxwell,device_rk4_chain
#ifdef JETSPIN_COMPARE_MAXWELL_STAGE1
  logical, save :: stage1_reported=.false.
  double precision :: stage1_max_abs,stage1_max_rel,stage1_ref,stage1_gpu
#endif
#ifdef JETSPIN_COMPARE_MAXWELL_STAGE4
  logical, save :: stage4_reported=.false.
  double precision :: stage4_max_abs,stage4_max_rel,stage4_ref,stage4_gpu
  double precision :: stage4_comp(8)
#endif
#ifdef JETSPIN_COMPARE_MAXWELL_STAGE2
  double precision :: compare_comp(8)
#endif
#ifdef JETSPIN_COMPARE_MAXWELL_STAGES
  integer :: compare_stage
  integer :: compare_ipoint
  double precision :: compare_abs,compare_ref,compare_gpu
  double precision :: compare_local
  double precision :: geom_down,geom_up,geom_curv,geom_center(3),geom_norm(3),geom_tan(3)
  logical :: geom_straight
  double precision :: diag_cmass,diag_field,diag_axial,diag_drag,diag_veltan,diag_lup
  double precision :: compare_comp(8)
#endif

#ifdef _OPENACC
#if defined(JETSPIN_DEV_HOST_FORCE_ORACLE) || defined(JETSPIN_DEV_HOST_COULOMB_ORACLE)
  if(evaporative_dynamic_accelerator_eligible())then
    call rk4sys_maxwell_ev_device(timesub,h,k)
    return
  endif
#endif
#endif

  device_rk4_chain=.false.
#if defined(_OPENACC) && !defined(JETSPIN_DEV_HOST_FORCE_ORACLE) && !defined(JETSPIN_DEV_HOST_MAXWELL_STATE_UPDATE) && !defined(JETSPIN_COMPARE_MAXWELL_STAGE1) && !defined(JETSPIN_COMPARE_MAXWELL_STAGE2) && !defined(JETSPIN_COMPARE_MAXWELL_STAGE4) && !defined(JETSPIN_COMPARE_MAXWELL_STAGES) && !defined(JETSPIN_TRACE_MAXWELL_BEAD) && !defined(JETSPIN_PRINT_MAXWELL_STAGE2)
  device_rk4_chain=evaporative_dynamic_accelerator_eligible()
#endif

#ifdef _OPENACC
  if(persistent_reset_requested)then
    if(doallocate .and. workspace_device_mapped)then
!$acc exit data delete(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
!$acc& f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
!$acc& f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
!$acc& f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev)
      workspace_device_mapped=.false.
    endif
    call reset_coulomb_accelerator(coulforce)
    persistent_acc=.false.
    persistent_reset_requested=.false.
  endif
  if(.not.persistent_acc .and. evaporative_dynamic_accelerator_eligible() .and. &
   .not.accelerator_is_topology_enabled())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet),jetzz(0:mxnpjet), &
!$acc& jetst(0:mxnpjet),jetvx(0:mxnpjet),jetvy(0:mxnpjet),jetvz(0:mxnpjet), &
!$acc& jetvl(0:mxnpjet),jetve(0:mxnpjet),jetce(0:mxnpjet),jetms(0:mxnpjet), &
!$acc& jetch(0:mxnpjet),jetfr(0:mxnpjet))
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(device_rk4_chain)
    call accelerator_set_topology_enabled(.true.)
    persistent_acc=.true.
  endif

  if(evaporative_dynamic_accelerator_eligible())then
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(device_rk4_chain)
    persistent_acc=.true.
  endif

#endif

! check and eventually reallocate the service arrays
  if(doallocate)then
#ifdef _OPENACC
    if(workspace_device_mapped)then
!$acc exit data delete(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
!$acc& f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
!$acc& f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
!$acc& f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev)
      workspace_device_mapped=.false.
    endif
#endif
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
#ifdef _OPENACC
    if((persistent_acc .or. accelerator_is_topology_enabled()) .and. systype/=1)then
!$acc enter data create(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
!$acc& f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
!$acc& f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
!$acc& f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev)
      workspace_device_mapped=.true.
      persistent_acc=.true.
    endif
#endif
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
#ifdef _OPENACC
      if(device_rk4_chain)then
        call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
         jetxx,jetyy,jetzz)
      else
#endif
      call compute_posnoinserted(jetxx,jetyy,jetzz)
#ifdef _OPENACC
      endif
#endif
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
      used_acc_maxwell=.false.
#ifdef _OPENACC
      if(evaporative_dynamic_accelerator_eligible())then
        call accelerator_maxwell_evap_stage(mystart,myend,npjet,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetvl,jetve,coulforce,jetms,jetch,jetfr, &
         f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev,linserting,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0,evairv, &
         evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev)
        if(.not.device_rk4_chain)then
!$acc update self(f1xx(0:npjet),f1yy(0:npjet),f1zz(0:npjet),f1st(0:npjet), &
!$acc& f1vx(0:npjet),f1vy(0:npjet),f1vz(0:npjet),f1ev(0:npjet)) if_present
        endif
        used_acc_maxwell=.true.
#ifdef JETSPIN_COMPARE_MAXWELL_STAGE1
        stage1_max_abs=0.d0
        stage1_max_rel=0.d0
        do ipoint=mystart,myend
          j=ipoint-mystart
          call xpsys_ev_maxwell(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
           jetvl,jetve,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub,k)
          stage1_gpu=f1xx(j); stage1_ref=fxx
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1yy(j); stage1_ref=fyy
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1zz(j); stage1_ref=fzz
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1st(j); stage1_ref=fst
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1vx(j); stage1_ref=fvx
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1vy(j); stage1_ref=fvy
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1vz(j); stage1_ref=fvz
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
          stage1_gpu=f1ev(j); stage1_ref=fev
          stage1_max_abs=max(stage1_max_abs,abs(stage1_gpu-stage1_ref))
          stage1_max_rel=max(stage1_max_rel,abs(stage1_gpu-stage1_ref)/max(1.d-30,abs(stage1_ref)))
        enddo
        if(.not.stage1_reported)then
          write(*,'(a,1pe14.6,a,1pe14.6)') 'Maxwell stage-1 CPU/GPU derivative check: max_abs=', &
           stage1_max_abs,' max_rel=',stage1_max_rel
          stage1_reported=.true.
        endif
#endif
      endif
#endif
      if(device_rk4_chain)then
#ifdef _OPENACC
        call accelerator_maxwell_rk4_stage_update(mystart,myend,h,1, &
         jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
         f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
         yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
#endif
      elseif(used_acc_maxwell)then
        j=0
        do ipoint=mystart,myend
          yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
          yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f1yy(j)
          yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f1zz(j)
          yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
          yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
          yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f1vy(j)
          yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f1vz(j)
          yev(ipoint) = jetve(ipoint) + 0.5d0*h*f1ev(j)
          if((yev(ipoint)/jetvl(ipoint))<evlim)yev(ipoint)=jetvl(ipoint)*evlim
          j=j+1
        enddo
      else
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
      endif
      call restore_charge()
      if(.not.device_rk4_chain)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
        call sum_world_darr(yev,npjet+1)
      endif
!     2°step
      call smooth_charge(yxx,yyy,yzz)
#ifdef _OPENACC
      if(device_rk4_chain)then
        call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
         yxx,yyy,yzz)
      else
#endif
      call compute_posnoinserted(yxx,yyy,yzz)
#ifdef _OPENACC
      endif
#endif
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
      j=0
      used_acc_maxwell=.false.
#ifdef _OPENACC
      if(evaporative_dynamic_accelerator_eligible())then
        if(.not.device_rk4_chain)then
!$acc update device(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet),yst(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yev(0:npjet),jetch(0:npjet)) if_present
        endif
        call accelerator_maxwell_evap_stage(mystart,myend,npjet,yxx,yyy,yzz,yst, &
         yvx,yvy,yvz,jetvl,yev,coulforce,jetms,jetch,jetfr, &
         f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev,linserting,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0,evairv, &
         evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev)
        if(.not.device_rk4_chain)then
!$acc update self(f2xx(0:npjet),f2yy(0:npjet),f2zz(0:npjet),f2st(0:npjet), &
!$acc& f2vx(0:npjet),f2vy(0:npjet),f2vz(0:npjet),f2ev(0:npjet)) if_present
        endif
#ifdef JETSPIN_COMPARE_MAXWELL_STAGE2
        compare_comp(:)=0.d0
        do ipoint=mystart,myend
          if(ipoint<=mystart .or. ipoint>=npjet)cycle
          j=ipoint-mystart
          call xpsys_ev_maxwell(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl,yev,coulforce, &
           fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub+h/2.d0,k)
          compare_comp(1)=max(compare_comp(1),abs(f2xx(j)-fxx))
          compare_comp(2)=max(compare_comp(2),abs(f2yy(j)-fyy))
          compare_comp(3)=max(compare_comp(3),abs(f2zz(j)-fzz))
          compare_comp(4)=max(compare_comp(4),abs(f2st(j)-fst))
          compare_comp(5)=max(compare_comp(5),abs(f2vx(j)-fvx))
          compare_comp(6)=max(compare_comp(6),abs(f2vy(j)-fvy))
          compare_comp(7)=max(compare_comp(7),abs(f2vz(j)-fvz))
          compare_comp(8)=max(compare_comp(8),abs(f2ev(j)-fev))
        enddo
        write(*,'(a,8(1pe12.4,1x))') 'Stage-2 abs components [xx yy zz st vx vy vz ev]=',compare_comp
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        j=50-mystart
        call xpsys_ev_maxwell(50,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl,yev,coulforce, &
         fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub+h/2.d0,k)
        write(*,'(a,8(1pe14.6,1x))') 'Stage-2 host-fallback bead50 [gpu cpu]=',f2xx(j),fxx,f2yy(j),fyy,f2zz(j),fzz,f2vx(j),fvx
#endif
#endif
#ifdef JETSPIN_TRACE_MAXWELL_BEAD
        if(timesub==0.d0)write(*,'(a,8(1pe14.6,1x))')'TRACE S2 GPU=',f2xx(50),f2yy(50),f2zz(50),f2st(50),f2vx(50),f2vy(50),f2vz(50),f2ev(50)
#endif
#ifdef JETSPIN_PRINT_MAXWELL_STAGE2
        if(timesub==0.d0)write(*,'(a,8(1pe14.6,1x))')'Stage-2 GPU bead50=', &
         f2xx(50),f2yy(50),f2zz(50),f2st(50),f2vx(50),f2vy(50),f2vz(50),f2ev(50)
#endif
        used_acc_maxwell=.true.
      endif
#endif
      if(.not.used_acc_maxwell)then
      j=0
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),f2ev(j),timesub+h/2.d0,k)
         j=j+1
      enddo
#ifdef JETSPIN_TRACE_MAXWELL_BEAD
      if(timesub==0.d0)write(*,'(a,8(1pe14.6,1x))')'TRACE S2 CPU=',f2xx(50),f2yy(50),f2zz(50),f2st(50),f2vx(50),f2vy(50),f2vz(50),f2ev(50)
#endif
      endif
#ifdef JETSPIN_PRINT_MAXWELL_STAGE2
      if(timesub==0.d0 .and. .not.used_acc_maxwell)write(*,'(a,8(1pe14.6,1x))')'Stage-2 CPU bead50=', &
       f2xx(50),f2yy(50),f2zz(50),f2st(50),f2vx(50),f2vy(50),f2vz(50),f2ev(50)
#endif
#if defined(_OPENACC) && !defined(JETSPIN_DISABLE_MAXWELL_EVAP)
      if(.not.used_acc_maxwell)call accelerator_maxwell_evap_stress_3d(mystart,myend,npjet,linserting,linserted, &
       jetfr,f2ev,f2st,yxx,yyy,yzz,yvx,yvy,yvz,yst,jetvl,yev, &
       evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev, &
       consistency,findex,yieldstress)
#endif
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
#if defined(_OPENACC) && !defined(JETSPIN_DEV_HOST_FORCE_ORACLE) && !defined(JETSPIN_DEV_HOST_MAXWELL_STATE_UPDATE)
      call accelerator_maxwell_rk4_stage_update(mystart,myend,h,2,jetxx,jetyy,jetzz,jetst, &
       jetvx,jetvy,jetvz,jetve,jetvl,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
       yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
      if(.not.device_rk4_chain)then
!$acc update self(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet),yst(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yev(0:npjet)) if_present
      endif
#else
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
#endif
      call restore_charge()
      if(.not.device_rk4_chain)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
        call sum_world_darr(yev,npjet+1)
      endif
!     3°step
      call smooth_charge(yxx,yyy,yzz)
#ifdef _OPENACC
      if(device_rk4_chain)then
        call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
         yxx,yyy,yzz)
      else
#endif
      call compute_posnoinserted(yxx,yyy,yzz)
#ifdef _OPENACC
      endif
#endif
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
#ifdef _OPENACC
      if(evaporative_dynamic_accelerator_eligible())then
        if(.not.device_rk4_chain)then
!$acc update device(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet),yst(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yev(0:npjet),jetch(0:npjet)) if_present
        endif
      endif
#endif
      j=0
      used_acc_maxwell=.false.
#ifdef _OPENACC
      if(evaporative_dynamic_accelerator_eligible())then
        call accelerator_maxwell_evap_stage(mystart,myend,npjet,yxx,yyy,yzz,yst, &
         yvx,yvy,yvz,jetvl,yev,coulforce,jetms,jetch,jetfr, &
         f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev,linserting,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0,evairv, &
         evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev)
        if(.not.device_rk4_chain)then
!$acc update self(f3xx(0:npjet),f3yy(0:npjet),f3zz(0:npjet),f3st(0:npjet), &
!$acc& f3vx(0:npjet),f3vy(0:npjet),f3vz(0:npjet),f3ev(0:npjet)) if_present
        endif
#ifdef JETSPIN_TRACE_MAXWELL_BEAD
        if(timesub==0.d0)write(*,'(a,8(1pe14.6,1x))')'TRACE S3 GPU=',f3xx(50),f3yy(50),f3zz(50),f3st(50),f3vx(50),f3vy(50),f3vz(50),f3ev(50)
#endif
#ifdef JETSPIN_COMPARE_MAXWELL_STAGES
        compare_abs=0.d0
        compare_ipoint=-1
        compare_comp(:)=0.d0
        do ipoint=mystart,myend
          j=ipoint-mystart
          call xpsys_ev_maxwell(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl,yev,coulforce, &
           fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub+h/2.d0,k)
          compare_local=max(abs(f3xx(j)-fxx),abs(f3yy(j)-fyy),abs(f3zz(j)-fzz), &
           abs(f3st(j)-fst),abs(f3vx(j)-fvx),abs(f3vy(j)-fvy),abs(f3vz(j)-fvz),abs(f3ev(j)-fev))
          if(compare_local>compare_abs)compare_ipoint=ipoint
          compare_abs=max(compare_abs,compare_local)
          compare_comp(1)=max(compare_comp(1),abs(f3xx(j)-fxx)); compare_comp(2)=max(compare_comp(2),abs(f3yy(j)-fyy))
          compare_comp(3)=max(compare_comp(3),abs(f3zz(j)-fzz)); compare_comp(4)=max(compare_comp(4),abs(f3st(j)-fst))
          compare_comp(5)=max(compare_comp(5),abs(f3vx(j)-fvx)); compare_comp(6)=max(compare_comp(6),abs(f3vy(j)-fvy))
          compare_comp(7)=max(compare_comp(7),abs(f3vz(j)-fvz)); compare_comp(8)=max(compare_comp(8),abs(f3ev(j)-fev))
        enddo
        write(*,'(a,1pe14.6)') 'Maxwell stage-3 max absolute CPU/GPU difference: ',compare_abs
        write(*,'(a,8(1pe12.4,1x))') 'Stage-3 abs components [xx yy zz st vx vy vz ev]=',compare_comp
        if(compare_ipoint>=0)write(*,'(a,i6,6(1pe14.6,1x))') 'Stage-3 max bead/state [i x y z yve yvl ycf-x]=',compare_ipoint,yxx(compare_ipoint),yyy(compare_ipoint),yzz(compare_ipoint),yev(compare_ipoint),jetvl(compare_ipoint),coulforce(compare_ipoint,1)
        if(compare_ipoint>0 .and. compare_ipoint<npjet)then
          call compute_geometry(compare_ipoint,yxx,yyy,yzz,geom_down,geom_up)
          call compute_tangetversor(compare_ipoint,yxx,yyy,yzz,geom_tan,geom_up)
          call compute_curvcenter(compare_ipoint,yxx,yyy,yzz,geom_center,geom_straight)
          call compute_curvature(compare_ipoint,yxx,yyy,yzz,geom_curv,geom_norm,geom_center,geom_straight)
          write(*,'(a,9(1pe14.6,1x),1x,l1)') 'Stage-3 CPU geometry down up tanx tany tanz curv nx ny nz straight=',geom_down,geom_up,geom_tan(1),geom_tan(2),geom_tan(3),geom_curv,geom_norm(1),geom_norm(2),geom_norm(3),geom_straight
          diag_cmass=yev(compare_ipoint)/jetvl(compare_ipoint)
          call project_veltangetversor(compare_ipoint,yvx,yvy,yvz,diag_veltan,geom_tan)
          diag_lup=geom_up
          diag_field=jetch(compare_ipoint)/(jetms(compare_ipoint)*diag_cmass)*v
          diag_axial=(fve/(jetms(compare_ipoint)*diag_cmass))*yev(compare_ipoint)*yst(compare_ipoint)/diag_lup
          diag_drag=(att/(jetms(compare_ipoint)*diag_cmass))*(abs(diag_lup)**0.905d0)*(abs(diag_veltan)**1.19d0)
          write(*,'(a,8(1pe14.6,1x))') 'Stage-3 CPU fvx terms [total gravity field axial drag coulomb]=',fvx,gr,diag_field,diag_axial,diag_drag,coulforce(compare_ipoint,1),diag_cmass,diag_veltan
        endif
#endif
        used_acc_maxwell=.true.
      endif
#endif
      if(.not.used_acc_maxwell)then
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
         f3vx(j),f3vy(j),f3vz(j),f3ev(j),timesub+h/2.d0,k)
         j=j+1
      enddo
#ifdef JETSPIN_TRACE_MAXWELL_BEAD
      if(timesub==0.d0)write(*,'(a,8(1pe14.6,1x))')'TRACE S3 CPU=',f3xx(50),f3yy(50),f3zz(50),f3st(50),f3vx(50),f3vy(50),f3vz(50),f3ev(50)
#endif
      endif
#if defined(_OPENACC) && !defined(JETSPIN_DISABLE_MAXWELL_EVAP)
      if(.not.used_acc_maxwell)call accelerator_maxwell_evap_stress_3d(mystart,myend,npjet,linserting,linserted, &
       jetfr,f3ev,f3st,yxx,yyy,yzz,yvx,yvy,yvz,yst,jetvl,yev, &
       evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev, &
       consistency,findex,yieldstress)
#endif
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
#if defined(_OPENACC) && !defined(JETSPIN_DEV_HOST_FORCE_ORACLE) && !defined(JETSPIN_DEV_HOST_MAXWELL_STATE_UPDATE)
      call accelerator_maxwell_rk4_stage_update(mystart,myend,h,3,jetxx,jetyy,jetzz,jetst, &
       jetvx,jetvy,jetvz,jetve,jetvl,f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
       yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
      if(.not.device_rk4_chain)then
!$acc update self(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet),yst(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yev(0:npjet)) if_present
      endif
#else
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
#endif
      call restore_charge()
      if(.not.device_rk4_chain)then
        call sum_world_darr(yxx,npjet+1)
        call sum_world_darr(yyy,npjet+1)
        call sum_world_darr(yzz,npjet+1)
        call sum_world_darr(yst,npjet+1)
        call sum_world_darr(yvx,npjet+1)
        call sum_world_darr(yvy,npjet+1)
        call sum_world_darr(yvz,npjet+1)
        call sum_world_darr(yev,npjet+1)
      endif
!     4°step
      call smooth_charge(yxx,yyy,yzz)
#ifdef _OPENACC
      if(device_rk4_chain)then
        call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
         yxx,yyy,yzz)
      else
#endif
      call compute_posnoinserted(yxx,yyy,yzz)
#ifdef _OPENACC
      endif
#endif
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz,yev)
#ifdef _OPENACC
      if(evaporative_dynamic_accelerator_eligible())then
        if(.not.device_rk4_chain)then
!$acc update device(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet),yst(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yev(0:npjet),jetch(0:npjet)) if_present
        endif
      endif
#endif
      j=0
      used_acc_maxwell=.false.
#ifdef _OPENACC
      if(evaporative_dynamic_accelerator_eligible())then
        call accelerator_maxwell_evap_stage(mystart,myend,npjet,yxx,yyy,yzz,yst, &
         yvx,yvy,yvz,jetvl,yev,coulforce,jetms,jetch,jetfr, &
         f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev,linserting,linserted,liniperturb, &
         lairdrag,lflorentz,luppot,nfieldtype,pfreq,consistency,findex, &
         yieldstress,att,fve,gr,ks,li,v,velext,.false.,0.d0,evairv, &
         evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev)
        if(.not.device_rk4_chain)then
!$acc update self(f4xx(0:npjet),f4yy(0:npjet),f4zz(0:npjet),f4st(0:npjet), &
!$acc& f4vx(0:npjet),f4vy(0:npjet),f4vz(0:npjet),f4ev(0:npjet)) if_present
        endif
#ifdef JETSPIN_COMPARE_MAXWELL_STAGES
        compare_abs=0.d0
        compare_comp(:)=0.d0
        do ipoint=mystart,myend
          j=ipoint-mystart
          call xpsys_ev_maxwell(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,jetvl,yev,coulforce, &
           fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub+h/2.d0,k)
          compare_abs=max(compare_abs,abs(f4xx(j)-fxx),abs(f4yy(j)-fyy),abs(f4zz(j)-fzz), &
           abs(f4st(j)-fst),abs(f4vx(j)-fvx),abs(f4vy(j)-fvy),abs(f4vz(j)-fvz),abs(f4ev(j)-fev))
          compare_comp(1)=max(compare_comp(1),abs(f4xx(j)-fxx)); compare_comp(2)=max(compare_comp(2),abs(f4yy(j)-fyy))
          compare_comp(3)=max(compare_comp(3),abs(f4zz(j)-fzz)); compare_comp(4)=max(compare_comp(4),abs(f4st(j)-fst))
          compare_comp(5)=max(compare_comp(5),abs(f4vx(j)-fvx)); compare_comp(6)=max(compare_comp(6),abs(f4vy(j)-fvy))
          compare_comp(7)=max(compare_comp(7),abs(f4vz(j)-fvz)); compare_comp(8)=max(compare_comp(8),abs(f4ev(j)-fev))
        enddo
        write(*,'(a,1pe14.6)') 'Maxwell stage-4 max absolute CPU/GPU difference: ',compare_abs
        write(*,'(a,8(1pe12.4,1x))') 'Stage-4 abs components [xx yy zz st vx vy vz ev]=',compare_comp
#endif
        used_acc_maxwell=.true.
#ifdef JETSPIN_COMPARE_MAXWELL_STAGE4
        stage4_max_abs=0.d0
        stage4_max_rel=0.d0
        stage4_comp(:)=0.d0
        do ipoint=mystart,myend
          j=ipoint-mystart
          call xpsys_ev_maxwell(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
           jetvl,yev,coulforce,fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub+h,k)
          stage4_gpu=f4xx(j); stage4_ref=fxx
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(1)=max(stage4_comp(1),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4yy(j); stage4_ref=fyy
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(2)=max(stage4_comp(2),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4zz(j); stage4_ref=fzz
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(3)=max(stage4_comp(3),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4st(j); stage4_ref=fst
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(4)=max(stage4_comp(4),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4vx(j); stage4_ref=fvx
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(5)=max(stage4_comp(5),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4vy(j); stage4_ref=fvy
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(6)=max(stage4_comp(6),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4vz(j); stage4_ref=fvz
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(7)=max(stage4_comp(7),abs(stage4_gpu-stage4_ref))
          stage4_gpu=f4ev(j); stage4_ref=fev
          stage4_max_abs=max(stage4_max_abs,abs(stage4_gpu-stage4_ref))
          stage4_max_rel=max(stage4_max_rel,abs(stage4_gpu-stage4_ref)/max(1.d-30,abs(stage4_ref)))
          stage4_comp(8)=max(stage4_comp(8),abs(stage4_gpu-stage4_ref))
        enddo
        if(.not.stage4_reported)then
          write(*,'(a,1pe14.6,a,1pe14.6)') 'Maxwell stage-4 CPU/GPU derivative check: max_abs=', &
           stage4_max_abs,' max_rel=',stage4_max_rel
          write(*,'(a,8(1pe12.4,1x))') 'Stage-4 abs components [xx yy zz st vx vy vz ev]=',stage4_comp
          stage4_reported=.true.
        endif
#endif
      endif
#endif
      if(.not.used_acc_maxwell)then
      do ipoint=mystart,myend
        call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,yev,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
         f4vx(j),f4vy(j),f4vz(j),f4ev(j),timesub+h,k)
         j=j+1
      enddo
      endif
#if defined(_OPENACC) && !defined(JETSPIN_DISABLE_MAXWELL_EVAP)
      if(.not.used_acc_maxwell)call accelerator_maxwell_evap_stress_3d(mystart,myend,npjet,linserting,linserted, &
       jetfr,f4ev,f4st,yxx,yyy,yzz,yvx,yvy,yvz,yst,jetvl,yev, &
       evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev, &
       consistency,findex,yieldstress)
#endif
      j=0
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
      yev(:)=0.d0
#if defined(_OPENACC) && !defined(JETSPIN_DEV_HOST_FORCE_ORACLE) && !defined(JETSPIN_DEV_HOST_MAXWELL_STATE_UPDATE)
      call accelerator_maxwell_rk4_final_update(mystart,myend,h,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
       f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
       f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev,f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev, &
       yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
      if(.not.device_rk4_chain)then
!$acc update self(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet),yst(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yev(0:npjet)) if_present
      endif
#else
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
#endif
      call restore_charge()
      if(device_rk4_chain)then
#ifdef _OPENACC
        call accelerator_maxwell_commit_state(mystart,myend, &
         yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
         jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve, &
         counterlpath,ncounterlpath,maxstress,maxstressposx)
        call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
         jetxx,jetyy,jetzz)
        call accelerator_mark_device_state(.true.)
#endif
      else
        call sum_world_darr(yxx,npjet+1,jetxx)
        call sum_world_darr(yyy,npjet+1,jetyy)
        call sum_world_darr(yzz,npjet+1,jetzz)
        call sum_world_darr(yst,npjet+1,jetst)
        call sum_world_darr(yvx,npjet+1,jetvx)
        call sum_world_darr(yvy,npjet+1,jetvy)
        call sum_world_darr(yvz,npjet+1,jetvz)
        call sum_world_darr(yev,npjet+1,jetve)
        call compute_posnoinserted(jetxx,jetyy,jetzz)
#ifdef _OPENACC
        call accelerator_mark_device_state(.false.)
#endif
      endif
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
  logical, save :: persistent_acc=.false.
  
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

#ifdef _OPENACC
! Refinement can grow capacity before this call, so doallocate itself is a
! reset request even when the main topology loop did not issue one.
  if((persistent_reset_requested .or. doallocate) .and. persistent_acc)then
!$acc exit data delete(f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev,d3xx,d3yy,d3zz,d3st, &
!$acc& d3vx,d3vy,d3vz,d3ev,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y1ev, &
!$acc& y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,y2ev)
    call reset_coulomb_accelerator(coulforce)
    call accelerator_set_persistent(.false.)
    call set_coulomb_accelerator_persistent(.false.)
    persistent_acc=.false.
  endif
  persistent_reset_requested=.false.
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

#ifdef _OPENACC
  if(.not.persistent_acc .and. (fixed_evaporative_platen_eligible() .or. &
   dynamic_evaporative_platen_eligible()) .and. &
   allocated(gaussianhistory))then
    if(.not.accelerator_is_topology_enabled())then
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet), &
!$acc& jetzz(0:mxnpjet),jetst(0:mxnpjet),jetvx(0:mxnpjet), &
!$acc& jetvy(0:mxnpjet),jetvz(0:mxnpjet),jetvl(0:mxnpjet), &
!$acc& jetve(0:mxnpjet),jetce(0:mxnpjet),jetms(0:mxnpjet), &
!$acc& jetch(0:mxnpjet),jetfr(0:mxnpjet))
    endif
!$acc enter data create(f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev,d3xx,d3yy,d3zz,d3st, &
!$acc& d3vx,d3vy,d3vz,d3ev,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y1ev, &
!$acc& y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,y2ev)
    call set_coulomb_accelerator_persistent(.true.)
    call accelerator_set_persistent(.true.)
    if(dynamic_evaporative_platen_eligible()) &
     call accelerator_set_topology_enabled(.true.)
    persistent_acc=.true.
  endif
#endif
  
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
#ifdef _OPENACC
      if(persistent_acc)then
        call maxwell_evap_device_stage(timesub,k,jetxx,jetyy,jetzz,jetst, &
         jetvx,jetvy,jetvz,jetve,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
         f1ev,.true.)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(f1xx(0:myend-mystart),f1yy(0:myend-mystart), &
!$acc& f1zz(0:myend-mystart),f1st(0:myend-mystart), &
!$acc& f1vx(0:myend-mystart),f1vy(0:myend-mystart), &
!$acc& f1vz(0:myend-mystart),f1ev(0:myend-mystart))
#endif
        call accelerator_platen_evap_predict(mystart,myend,h,airdragamp(1), &
         noisediff,evlim,jetms,jetvl,jetve,jetxx,jetyy,jetzz,jetst,jetvx, &
         jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev,y1xx,y1yy, &
         y1zz,y1st,y1vx,y1vy,y1vz,y1ev,y2xx,y2yy,y2zz,y2st,y2vx,y2vy, &
         y2vz,y2ev)

        call maxwell_evap_device_stage(timesub,k,y1xx,y1yy,y1zz,y1st, &
         y1vx,y1vy,y1vz,y1ev,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
         f2ev,.true.)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(f2xx(0:myend-mystart),f2yy(0:myend-mystart), &
!$acc& f2zz(0:myend-mystart),f2st(0:myend-mystart), &
!$acc& f2vx(0:myend-mystart),f2vy(0:myend-mystart), &
!$acc& f2vz(0:myend-mystart),f2ev(0:myend-mystart))
#endif
        call maxwell_evap_device_stage(timesub,k,y2xx,y2yy,y2zz,y2st, &
         y2vx,y2vy,y2vz,y2ev,d3xx,d3yy,d3zz,d3st,d3vx,d3vy,d3vz, &
         d3ev,.true.)
#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
!$acc update device(d3xx(0:myend-mystart),d3yy(0:myend-mystart), &
!$acc& d3zz(0:myend-mystart),d3st(0:myend-mystart), &
!$acc& d3vx(0:myend-mystart),d3vy(0:myend-mystart), &
!$acc& d3vz(0:myend-mystart),d3ev(0:myend-mystart))
#endif

        call accelerator_platen_evap_velocity(mystart,myend,mxnpjet, &
         gaussianhistorysteps,k,h,airdragamp(1),noisediff,jetms,jetvl, &
         jetve,gaussianhistory,jetvx,jetvy,jetvz,f1vx,f1vy,f1vz,f2vx, &
         f2vy,f2vz,d3vx,d3vy,d3vz)

#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        call maxwell_evap_device_stage(timesub+h,k,y1xx,y1yy,y1zz,y1st, &
         jetvx,jetvy,jetvz,y1ev,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
         f2ev,.true.)
!$acc update device(f2xx(0:myend-mystart),f2yy(0:myend-mystart), &
!$acc& f2zz(0:myend-mystart),f2st(0:myend-mystart), &
!$acc& f2vx(0:myend-mystart),f2vy(0:myend-mystart), &
!$acc& f2vz(0:myend-mystart),f2ev(0:myend-mystart))
#else
        call accelerator_maxwell_evap_stress_3d(mystart,myend,npjet, &
         linserting,linserted,jetfr,f2ev,f2st,y1xx,y1yy,y1zz,jetvx,jetvy, &
         jetvz,y1st,jetvl,y1ev,evairv,evmasscoeff,sqrevsc,evcsvapour, &
         evumidity,cp0,Bev,mev,tev,consistency,findex,yieldstress)
#endif
        call accelerator_platen_evap_positions(mystart,myend,npjet,h,pfreq, &
         liniperturb,evlim,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetvl,jetve, &
         f1xx,f1yy,f1zz,f1ev,f2ev)

#ifdef JETSPIN_DEV_HOST_FORCE_ORACLE
        call maxwell_evap_device_stage(timesub+h,k,jetxx,jetyy,jetzz,y1st, &
         jetvx,jetvy,jetvz,jetve,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
         f2ev,.true.)
!$acc update device(f2xx(0:myend-mystart),f2yy(0:myend-mystart), &
!$acc& f2zz(0:myend-mystart),f2st(0:myend-mystart), &
!$acc& f2vx(0:myend-mystart),f2vy(0:myend-mystart), &
!$acc& f2vz(0:myend-mystart),f2ev(0:myend-mystart))
#else
        call accelerator_maxwell_stress_3d(mystart,myend,npjet,linserted, &
         jetfr,f2st,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,y1st,jetvl,jetve, &
         cp0,Bev,mev,tev,consistency,findex,yieldstress)
#endif
        call accelerator_platen_stress_statistics(mystart,myend,h,jetxx, &
         jetyy,jetzz,jetst,f1st,f2st,counterlpath,ncounterlpath,maxstress, &
         maxstressposx)
        call accelerator_compute_posnoinserted_3d(npjet,linserted,resolution, &
         jetxx,jetyy,jetzz)
        call accelerator_mark_device_state(.true.)
        timesub=timesub+h
        return
      endif
#endif
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
