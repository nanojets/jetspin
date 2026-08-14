module integrator_kv_ev_mod

 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

!***********************************************************************
! Kelvin-Voigt integration with solvent evaporation.
!
! This module extends the historical JETSPIN Kelvin-Voigt integrators
! to the concentration-dependent material properties used by the
! Yarin-Koombhongse-Reneker evaporation model.  The ordinary Kelvin-
! Voigt path in integrator_mod is left unchanged for backward
! compatibility when evaporation is disabled.
!***********************************************************************

 use version_mod, only : mystart,myend,mxchunk,sum_world_darr, &
                         set_chunk,set_mxchunk,idrank,mxrank
 use error_mod, only : error
 use nanojet_mod, only : mxnpjet,npjet,inpjet,systype,jetxx,jetyy, &
                         jetzz,jetst,jetvx,jetvy,jetvz,jetvl,jetve, &
                         jetms,jetch,jetce, &
                         compute_posnoinserted,evlim,linserting,linserted, &
                         jetfr,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity, &
                         cp0,Bev,mev,tev
 use integrator_mod, only : integrator
 use dynamic_refinement_mod, only : driver_dynamic_refinement
 use coulomb_force_mod, only : smooth_charge,restore_charge,coulforce, &
                               compute_coulomelec_driver
 use eom_ev_mod, only : eom1_KV_pos_v_ev,eom1_KV_st_ev, &
                        eom3_KV_pos_v_ev,eom3_KV_st_ev
#ifdef _OPENACC
 use accelerator_mod, only : accelerator_enabled, &
                             accelerator_kv_evap_stress_3d, &
                             accelerator_set_topology_enabled
#endif

 implicit none
 private
 public :: driver_integrator_KV_ev

 double precision, allocatable, save :: fxx(:,:),fyy(:,:),fzz(:,:)
 double precision, allocatable, save :: fst(:,:),fev(:,:)
 double precision, allocatable, save :: f1vx(:),f1vy(:),f1vz(:)
 double precision, allocatable, save :: f2vx(:),f2vy(:),f2vz(:)
 double precision, allocatable, save :: f3vx(:),f3vy(:),f3vz(:)
 double precision, allocatable, save :: f4vx(:),f4vy(:),f4vz(:)
 double precision, allocatable, save :: yxx(:),yyy(:),yzz(:),yst(:)
 double precision, allocatable, save :: yvx(:),yvy(:),yvz(:),yev(:)
 logical, save :: lworkspace=.false.
 logical, save :: workspace_device_mapped=.false.
 logical, save :: lannounced=.false.
 integer, save :: workspace_mxnpjet=-1
 integer, save :: workspace_mxchunk=-1

 contains

 subroutine driver_integrator_KV_ev(timesub,h,k,dorefinment)

  implicit none
  logical, intent(inout) :: dorefinment
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: i
  logical :: ltestinst

  call driver_dynamic_refinement(k,dorefinment)
  call set_chunk(inpjet,npjet)
  call set_mxchunk(mxnpjet)
  call ensure_workspace()

  if((.not.lannounced).and.idrank==0)then
    write(6,'(/,a,/)') &
     'Kelvin-Voigt evaporation integrator active (Yarin concentration laws)'
    lannounced=.true.
  endif

  select case(integrator)
    case(1)
      call eulsys_KV_ev(timesub,h,k)
    case(2)
      call rk2sys_KV_ev(timesub,h,k)
    case(3)
      call rk4sys_KV_ev(timesub,h,k)
    case default
      call error(1)
  end select

  ltestinst=.false.
  do i=inpjet,npjet
    if(ieee_is_nan(dcos(jetxx(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetyy(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetzz(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetst(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetvx(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetvy(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetvz(i))))ltestinst=.true.
    if(ieee_is_nan(dcos(jetve(i))))ltestinst=.true.
  enddo
  if(ltestinst)call error(14)

  return
 end subroutine driver_integrator_KV_ev

 subroutine ensure_workspace()
  implicit none

  if(lworkspace)then
    if(workspace_mxnpjet>=mxnpjet .and. workspace_mxchunk>=mxchunk)return

#ifdef _OPENACC
  if(workspace_device_mapped)then
!$acc exit data delete(fxx,fyy,fzz,fst,fev,f1vx,f1vy,f1vz,f2vx,f2vy,f2vz, &
!$acc& f3vx,f3vy,f3vz,f4vx,f4vy,f4vz,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
    workspace_device_mapped=.false.
  endif
#endif

    deallocate(fxx,fyy,fzz,fst,fev)
    deallocate(f1vx,f1vy,f1vz,f2vx,f2vy,f2vz)
    deallocate(f3vx,f3vy,f3vz,f4vx,f4vy,f4vz)
    deallocate(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
  endif

  allocate(fxx(0:mxchunk,4),fyy(0:mxchunk,4),fzz(0:mxchunk,4))
  allocate(fst(0:mxchunk,4),fev(0:mxchunk,4))
  allocate(f1vx(0:mxnpjet),f1vy(0:mxnpjet),f1vz(0:mxnpjet))
  allocate(f2vx(0:mxnpjet),f2vy(0:mxnpjet),f2vz(0:mxnpjet))
  allocate(f3vx(0:mxnpjet),f3vy(0:mxnpjet),f3vz(0:mxnpjet))
  allocate(f4vx(0:mxnpjet),f4vy(0:mxnpjet),f4vz(0:mxnpjet))
  allocate(yxx(0:mxnpjet),yyy(0:mxnpjet),yzz(0:mxnpjet))
  allocate(yst(0:mxnpjet),yvx(0:mxnpjet),yvy(0:mxnpjet))
  allocate(yvz(0:mxnpjet),yev(0:mxnpjet))
#ifdef _OPENACC
!$acc enter data copyin(jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms, &
!$acc& jetch,jetvl,jetve,jetce,jetfr)
!$acc enter data create(fxx,fyy,fzz,fst,fev,f1vx,f1vy,f1vz,f2vx,f2vy,f2vz, &
!$acc& f3vx,f3vy,f3vz,f4vx,f4vy,f4vz,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev)
  workspace_device_mapped=.true.
  call accelerator_set_topology_enabled(.true.)
#endif
  workspace_mxnpjet=mxnpjet
  workspace_mxchunk=mxchunk
  lworkspace=.true.

  return
 end subroutine ensure_workspace

 subroutine eval_stage(tstage,k,xs,ys,zs,ss,vxs,vys,vzs,ves, &
                       dx,dy,dz,ds,dve,ax,ay,az)
  implicit none
  integer, intent(in) :: k
  double precision, intent(in) :: tstage
  double precision, allocatable, dimension(:), intent(inout) :: xs,ys,zs
  double precision, allocatable, dimension(:), intent(in) :: ss,vxs,vys,vzs
  double precision, allocatable, dimension(:), intent(in) :: ves
  double precision, dimension(0:), intent(inout) :: dx,dy,dz,ds,dve
  double precision, allocatable, dimension(:), intent(inout) :: ax,ay,az
  integer :: ipoint,j

  dx(:)=0.d0
  dy(:)=0.d0
  dz(:)=0.d0
  ds(:)=0.d0
  dve(:)=0.d0
  ax(:)=0.d0
  ay(:)=0.d0
  az(:)=0.d0

  select case(systype)
    case(1)
      call smooth_charge(xs)
      call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
    case(3)
      call smooth_charge(xs,ys,zs)
      call compute_posnoinserted(xs,ys,zs)
      call compute_coulomelec_driver(k,tstage,coulforce,jetvl,xs,ys,zs,ves)
    case default
      call error(2)
  end select

  j=0
  do ipoint=mystart,myend
    select case(systype)
      case(1)
        call eom1_KV_pos_v_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs, &
         jetvl,ves,coulforce,dx(j),dy(j),dz(j),ax(ipoint),ay(ipoint), &
         az(ipoint),dve(j),tstage,k)
      case(3)
        call eom3_KV_pos_v_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs, &
         jetvl,ves,coulforce,dx(j),dy(j),dz(j),ax(ipoint),ay(ipoint), &
         az(ipoint),dve(j),tstage,k)
    end select
    j=j+1
  enddo

  call sum_world_darr(ax,npjet+1)
  call sum_world_darr(ay,npjet+1)
  call sum_world_darr(az,npjet+1)

  j=0
  do ipoint=mystart,myend
    select case(systype)
      case(1)
        call eom1_KV_st_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,ax,ay,az,dve(j),ds(j),tstage,k)
      case(3)
#ifdef _OPENACC
        if(accelerator_enabled .and. mxrank==1)then
          ds(j)=0.d0
        else
          call eom3_KV_st_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
           coulforce,ax,ay,az,dve(j),ds(j),tstage,k)
        endif
#else
        call eom3_KV_st_ev(ipoint,xs,ys,zs,ss,vxs,vys,vzs,jetvl,ves, &
         coulforce,ax,ay,az,dve(j),ds(j),tstage,k)
#endif
    end select
    j=j+1
  enddo

#ifdef _OPENACC
  if(accelerator_enabled .and. systype==3 .and. mxrank==1)then
    call accelerator_kv_evap_stress_3d(mystart,myend,npjet,linserting,linserted, &
     jetfr,dve,ds,xs,ys,zs,vxs,vys,vzs,ax,ay,az,ss,jetvl,ves,evairv, &
     evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev,evlim)
#ifdef _OPENACC
    !$acc update self(dve(0:myend-mystart),ds(0:myend-mystart))
#endif
  endif
#endif

  call restore_charge()

  return
 end subroutine eval_stage

 subroutine clamp_ev_volume(ipoint,value)
  implicit none
  integer, intent(in) :: ipoint
  double precision, intent(inout) :: value
  if((value/jetvl(ipoint))<evlim)value=jetvl(ipoint)*evlim
  return
 end subroutine clamp_ev_volume

 subroutine gather_state()
  implicit none
  call sum_world_darr(yxx,npjet+1)
  call sum_world_darr(yyy,npjet+1)
  call sum_world_darr(yzz,npjet+1)
  call sum_world_darr(yst,npjet+1)
  call sum_world_darr(yvx,npjet+1)
  call sum_world_darr(yvy,npjet+1)
  call sum_world_darr(yvz,npjet+1)
  call sum_world_darr(yev,npjet+1)
  return
 end subroutine gather_state

 subroutine commit_state()
  implicit none
  call sum_world_darr(yxx,npjet+1,jetxx)
  call sum_world_darr(yyy,npjet+1,jetyy)
  call sum_world_darr(yzz,npjet+1,jetzz)
  call sum_world_darr(yst,npjet+1,jetst)
  call sum_world_darr(yvx,npjet+1,jetvx)
  call sum_world_darr(yvy,npjet+1,jetvy)
  call sum_world_darr(yvz,npjet+1,jetvz)
  call sum_world_darr(yev,npjet+1,jetve)
  return
 end subroutine commit_state

 subroutine zero_state()
  implicit none
  yxx(:)=0.d0
  yyy(:)=0.d0
  yzz(:)=0.d0
  yst(:)=0.d0
  yvx(:)=0.d0
  yvy(:)=0.d0
  yvz(:)=0.d0
  yev(:)=0.d0
  return
 end subroutine zero_state

 subroutine eulsys_KV_ev(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: ipoint,j

  call eval_stage(timesub,k,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetve,fxx(:,1),fyy(:,1),fzz(:,1),fst(:,1),fev(:,1), &
   f1vx,f1vy,f1vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+h*fxx(j,1)
    yyy(ipoint)=jetyy(ipoint)+h*fyy(j,1)
    yzz(ipoint)=jetzz(ipoint)+h*fzz(j,1)
    yst(ipoint)=jetst(ipoint)+h*fst(j,1)
    yvx(ipoint)=jetvx(ipoint)+h*f1vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+h*f1vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+h*f1vz(ipoint)
    yev(ipoint)=jetve(ipoint)+h*fev(j,1)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call commit_state()
  timesub=timesub+h
  if(systype==3)call compute_posnoinserted(jetxx,jetyy,jetzz)
  return
 end subroutine eulsys_KV_ev

 subroutine rk2sys_KV_ev(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: ipoint,j

  call eval_stage(timesub,k,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetve,fxx(:,1),fyy(:,1),fzz(:,1),fst(:,1),fev(:,1), &
   f1vx,f1vy,f1vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+h*fxx(j,1)
    yyy(ipoint)=jetyy(ipoint)+h*fyy(j,1)
    yzz(ipoint)=jetzz(ipoint)+h*fzz(j,1)
    yst(ipoint)=jetst(ipoint)+h*fst(j,1)
    yvx(ipoint)=jetvx(ipoint)+h*f1vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+h*f1vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+h*f1vz(ipoint)
    yev(ipoint)=jetve(ipoint)+h*fev(j,1)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,2),fyy(:,2),fzz(:,2),fst(:,2),fev(:,2),f2vx,f2vy,f2vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+0.5d0*h*(fxx(j,1)+fxx(j,2))
    yyy(ipoint)=jetyy(ipoint)+0.5d0*h*(fyy(j,1)+fyy(j,2))
    yzz(ipoint)=jetzz(ipoint)+0.5d0*h*(fzz(j,1)+fzz(j,2))
    yst(ipoint)=jetst(ipoint)+0.5d0*h*(fst(j,1)+fst(j,2))
    yvx(ipoint)=jetvx(ipoint)+0.5d0*h*(f1vx(ipoint)+f2vx(ipoint))
    yvy(ipoint)=jetvy(ipoint)+0.5d0*h*(f1vy(ipoint)+f2vy(ipoint))
    yvz(ipoint)=jetvz(ipoint)+0.5d0*h*(f1vz(ipoint)+f2vz(ipoint))
    yev(ipoint)=jetve(ipoint)+0.5d0*h*(fev(j,1)+fev(j,2))
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call commit_state()
  timesub=timesub+h
  if(systype==3)call compute_posnoinserted(jetxx,jetyy,jetzz)
  return
 end subroutine rk2sys_KV_ev

 subroutine rk4sys_KV_ev(timesub,h,k)
  implicit none
  integer, intent(in) :: k
  double precision, intent(inout) :: timesub
  double precision, intent(in) :: h
  integer :: ipoint,j

  call eval_stage(timesub,k,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetve,fxx(:,1),fyy(:,1),fzz(:,1),fst(:,1),fev(:,1), &
   f1vx,f1vy,f1vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+0.5d0*h*fxx(j,1)
    yyy(ipoint)=jetyy(ipoint)+0.5d0*h*fyy(j,1)
    yzz(ipoint)=jetzz(ipoint)+0.5d0*h*fzz(j,1)
    yst(ipoint)=jetst(ipoint)+0.5d0*h*fst(j,1)
    yvx(ipoint)=jetvx(ipoint)+0.5d0*h*f1vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+0.5d0*h*f1vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+0.5d0*h*f1vz(ipoint)
    yev(ipoint)=jetve(ipoint)+0.5d0*h*fev(j,1)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+0.5d0*h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,2),fyy(:,2),fzz(:,2),fst(:,2),fev(:,2),f2vx,f2vy,f2vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+0.5d0*h*fxx(j,2)
    yyy(ipoint)=jetyy(ipoint)+0.5d0*h*fyy(j,2)
    yzz(ipoint)=jetzz(ipoint)+0.5d0*h*fzz(j,2)
    yst(ipoint)=jetst(ipoint)+0.5d0*h*fst(j,2)
    yvx(ipoint)=jetvx(ipoint)+0.5d0*h*f2vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+0.5d0*h*f2vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+0.5d0*h*f2vz(ipoint)
    yev(ipoint)=jetve(ipoint)+0.5d0*h*fev(j,2)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+0.5d0*h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,3),fyy(:,3),fzz(:,3),fst(:,3),fev(:,3),f3vx,f3vy,f3vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+h*fxx(j,3)
    yyy(ipoint)=jetyy(ipoint)+h*fyy(j,3)
    yzz(ipoint)=jetzz(ipoint)+h*fzz(j,3)
    yst(ipoint)=jetst(ipoint)+h*fst(j,3)
    yvx(ipoint)=jetvx(ipoint)+h*f3vx(ipoint)
    yvy(ipoint)=jetvy(ipoint)+h*f3vy(ipoint)
    yvz(ipoint)=jetvz(ipoint)+h*f3vz(ipoint)
    yev(ipoint)=jetve(ipoint)+h*fev(j,3)
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call gather_state()

  call eval_stage(timesub+h,k,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   fxx(:,4),fyy(:,4),fzz(:,4),fst(:,4),fev(:,4),f4vx,f4vy,f4vz)

  call zero_state()
  j=0
  do ipoint=mystart,myend
    yxx(ipoint)=jetxx(ipoint)+(h/6.d0)*(fxx(j,1)+ &
     2.d0*(fxx(j,2)+fxx(j,3))+fxx(j,4))
    yyy(ipoint)=jetyy(ipoint)+(h/6.d0)*(fyy(j,1)+ &
     2.d0*(fyy(j,2)+fyy(j,3))+fyy(j,4))
    yzz(ipoint)=jetzz(ipoint)+(h/6.d0)*(fzz(j,1)+ &
     2.d0*(fzz(j,2)+fzz(j,3))+fzz(j,4))
    yst(ipoint)=jetst(ipoint)+(h/6.d0)*(fst(j,1)+ &
     2.d0*(fst(j,2)+fst(j,3))+fst(j,4))
    yvx(ipoint)=jetvx(ipoint)+(h/6.d0)*(f1vx(ipoint)+ &
     2.d0*(f2vx(ipoint)+f3vx(ipoint))+f4vx(ipoint))
    yvy(ipoint)=jetvy(ipoint)+(h/6.d0)*(f1vy(ipoint)+ &
     2.d0*(f2vy(ipoint)+f3vy(ipoint))+f4vy(ipoint))
    yvz(ipoint)=jetvz(ipoint)+(h/6.d0)*(f1vz(ipoint)+ &
     2.d0*(f2vz(ipoint)+f3vz(ipoint))+f4vz(ipoint))
    yev(ipoint)=jetve(ipoint)+(h/6.d0)*(fev(j,1)+ &
     2.d0*(fev(j,2)+fev(j,3))+fev(j,4))
    call clamp_ev_volume(ipoint,yev(ipoint))
    j=j+1
  enddo
  call commit_state()
  timesub=timesub+h
  if(systype==3)call compute_posnoinserted(jetxx,jetyy,jetzz)
  return
 end subroutine rk4sys_KV_ev

end module integrator_kv_ev_mod
