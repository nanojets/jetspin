module accelerator_mod

#ifdef _OPENACC
 use openacc, only : acc_async_sync
#endif
 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.true.
 logical, save :: accelerator_persistent=.false.
 logical, save :: accelerator_topology_enabled=.false.
 logical, save :: accelerator_device_state_authoritative=.false.
 logical, save :: accelerator_statistics_mapped=.false.
 logical, save :: accelerator_eom_env_checked=.false.
 logical, save :: accelerator_eom_disabled=.false.
 double precision, save :: statistics_step_max=-huge(0.d0)
 integer, save :: statistics_step_index=-1
 integer, save :: accelerator_last_host_sync_step=-huge(0)
 double precision, save :: accelerator_smoothed_charge=0.d0
! Outcome of one timestep's device topology checks (npjet, linserted as 0/1,
! bead added, capacity exhausted, bead removed), read back in one transfer
! by accelerator_topology_check.
 integer, save :: accelerator_topology_flags(5)=0
! Set when accelerator_platen_end_step has already made this step's
! topology decisions and statistics update on the device: the following
! accelerator_topology_check only reads the flags back, and
! accelerator_store_statistics does nothing.
 logical, save :: accelerator_topology_decided=.false.
 logical, save :: accelerator_statistics_stored=.false.
! Results of one refinement threshold scan (last over-threshold segment,
! path length, nozzle correction), read back in one transfer.
 double precision, save :: accelerator_refinement_scan(3)=0.d0
! Queue of the persistent dynamic evaporative Platen step: its kernels are
! enqueued without waiting, and the host waits once per step for the
! topology record (accelerator_topology_check) and before every other
! transfer (accelerator_wait).  -1 (acc_async_sync) is synchronous, the
! setting of every other path.
#ifdef _OPENACC
 integer, save, public :: accelerator_queue=acc_async_sync
#else
 integer, save, public :: accelerator_queue=-1
#endif
#ifdef _OPENACC
!$acc declare create(accelerator_smoothed_charge)
!$acc declare create(accelerator_topology_flags)
!$acc declare create(accelerator_refinement_scan)
#endif

 public :: accelerator_prepare
 public :: accelerator_set_async
 public :: accelerator_wait
 public :: accelerator_eom3_stage
 public :: accelerator_maxwell_evap_stage
 public :: accelerator_kv_stage
 public :: accelerator_maxwell_evap_force_correction
 public :: accelerator_maxwell_rk4_stage_update
 public :: accelerator_maxwell_rk4_final_update
 public :: accelerator_maxwell_commit_state
 public :: accelerator_evap_rk4_stage_update
 public :: accelerator_evap_rk2_final_update
 public :: accelerator_evap_rk4_final_update
 public :: accelerator_evap_commit_state
 public :: accelerator_compute_posnoinserted_3d
 public :: accelerator_freeze_at_collector
 public :: accelerator_contact_bead
 public :: accelerator_smooth_charge_3d
 public :: accelerator_restore_charge
 public :: accelerator_coulomb_evap_3d
 public :: accelerator_evaporation_force_3d
 public :: accelerator_set_persistent
 public :: accelerator_set_topology_enabled
 public :: accelerator_is_topology_enabled
 public :: accelerator_rebind_topology
 public :: accelerator_rebind_evaporation
 public :: accelerator_update_device_topology_state
 public :: accelerator_update_device_evaporation_state
 public :: accelerator_refinement_candidate
 public :: accelerator_is_persistent
 public :: accelerator_statistics_on_device
 public :: accelerator_update_host_state
 public :: accelerator_update_host_capacity_state
 public :: accelerator_update_host_evaporation_state
 public :: accelerator_release_jet_capacity
 public :: accelerator_release_evaporation_capacity
 public :: accelerator_update_host_point
 public :: accelerator_update_host_evaporation_point
 public :: accelerator_finish_remove_bead
 public :: accelerator_update_host_removed_evaporation
 public :: accelerator_update_device_removed
 public :: accelerator_update_device_removed_evaporation
 public :: accelerator_topology_check
 public :: accelerator_update_device_added_evaporation
 public :: accelerator_host_state_is_current
 public :: accelerator_mark_device_state
 public :: accelerator_device_state_is_current
 public :: accelerator_store_statistics
 public :: accelerator_update_host_statistics
 public :: accelerator_update_device_statistics
 public :: accelerator_platen_predict
 public :: accelerator_platen_stage_prep
 public :: accelerator_platen_end_step
 public :: accelerator_platen_update

 ! The RK state algebra is common to both evaporation rheologies.  Preserve
 ! the historical specific procedure names while exposing neutral interfaces
 ! for the Kelvin--Voigt integrator.
 interface accelerator_evap_rk4_stage_update
   module procedure accelerator_maxwell_rk4_stage_update
 end interface
 interface accelerator_evap_rk4_final_update
   module procedure accelerator_maxwell_rk4_final_update
 end interface
 interface accelerator_evap_commit_state
   module procedure accelerator_maxwell_commit_state
 end interface

contains

 subroutine accelerator_set_async(enabled)
! Put the persistent step on one asynchronous queue, or back on the
! synchronous one.  JETSPIN_OPENACC_SYNC=1 keeps it synchronous.
  implicit none
  logical, intent(in) :: enabled
  character(len=16) :: env
  call accelerator_wait()
#ifdef _OPENACC
  accelerator_queue=acc_async_sync
  if(.not.enabled)return
  env=''
  call get_environment_variable('JETSPIN_OPENACC_SYNC',env)
  if(trim(env)=='1')return
  accelerator_queue=1
#endif
 end subroutine accelerator_set_async

 subroutine accelerator_wait()
! Complete the queued work before the host reads device data or a
! synchronous construct touches it.  A no-op on the synchronous queue.
  implicit none
#ifdef _OPENACC
  if(accelerator_queue/=acc_async_sync)then
!$acc wait(accelerator_queue)
  endif
#endif
 end subroutine accelerator_wait

 subroutine accelerator_maxwell_evap_force_correction(firstpoint,lastpoint,npjet, &
   linserting,yxx,yyy,yzz,yvx,yvy,yvz,yvl,yve,jetms,jetfr,att,li,fvx,fvy,fvz)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: linserting,jetfr(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yvx(0:),yvy(0:),yvz(0:)
  double precision, intent(in) :: yvl(0:),yve(0:),jetms(0:),att,li
  double precision, intent(inout) :: fvx(0:),fvy(0:),fvz(0:)
  integer :: ipoint,j
  double precision :: dx,dy,dz,lup,ldown,tux,tuy,tuz,tdx,tdy,tdz
  double precision :: v1x,v1y,v1z,v2x,v2y,v2z,l1,l2,dotp,nbx,nby,nbz,lnb
  double precision :: b,c,t,scale1,scale2,rcx,rcy,rcz,radius,curvature
  double precision :: factor3,factor4,factor5,veltangent,cmass,delta
  logical :: straight,use_evap
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx,yyy,yzz,yvx,yvy,yvz,yvl,yve,jetms,jetfr) &
!$acc& present_or_copy(fvx,fvy,fvz) private(j,dx,dy,dz,lup,ldown,tux,tuy,tuz,tdx,tdy,tdz, &
!$acc& v1x,v1y,v1z,v2x,v2y,v2z,l1,l2,dotp,nbx,nby,nbz,lnb,b,c,t,scale1,scale2,rcx,rcy,rcz, &
!$acc& radius,curvature,factor3,factor4,factor5,veltangent,cmass,delta,straight)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    if(jetfr(ipoint) .or. ipoint>=npjet)cycle
    if(ipoint==npjet-1 .and. .not.linserting)cycle
    if(yvl(ipoint)<=0.d0 .or. yve(ipoint)<=0.d0)cycle
    cmass=yve(ipoint)/yvl(ipoint)
    delta=1.d0/cmass-1.d0
    if(delta==0.d0)cycle
    dx=yxx(ipoint)-yxx(ipoint+1); dy=yyy(ipoint)-yyy(ipoint+1); dz=yzz(ipoint)-yzz(ipoint+1)
    lup=dsqrt(dx*dx+dy*dy+dz*dz)
    if(lup<=0.d0)cycle
    tux=dx/lup; tuy=dy/lup; tuz=dz/lup
    veltangent=(yvx(ipoint))*tux+yvy(ipoint)*tuy+yvz(ipoint)*tuz
    factor4=(att/jetms(ipoint))*(dabs(lup)**0.905d0)*(dabs(veltangent)**1.19d0)
    fvx(j)=fvx(j)-delta*factor4*tux
    fvy(j)=fvy(j)-delta*factor4*tuy
    fvz(j)=fvz(j)-delta*factor4*tuz
    if(ipoint==firstpoint)cycle
    dx=yxx(ipoint-1)-yxx(ipoint); dy=yyy(ipoint-1)-yyy(ipoint); dz=yzz(ipoint-1)-yzz(ipoint)
    ldown=dsqrt(dx*dx+dy*dy+dz*dz)
    if(ldown<=0.d0)cycle
    tdx=dx/ldown; tdy=dy/ldown; tdz=dz/ldown
    v1x=-tux; v1y=-tuy; v1z=-tuz
    v2x=dx/ldown; v2y=dy/ldown; v2z=dz/ldown
    l1=1.d0; l2=1.d0; dotp=v2x*v1x+v2y*v1y+v2z*v1z
    nbx=v2x-dotp*v1x; nby=v2y-dotp*v1y; nbz=v2z-dotp*v1z
    lnb=dsqrt(nbx*nbx+nby*nby+nbz*nbz)
    curvature=0.d0; rcx=0.d0; rcy=0.d0; rcz=0.d0
    if(lnb>0.d0)then
      nbx=nbx/lnb; nby=nby/lnb; nbz=nbz/lnb
      b=dx*v1x+dy*v1y+dz*v1z; c=dx*nbx+dy*nby+dz*nbz
      if(c/=0.d0)then
        t=0.5d0*(l1-b)/c; scale1=b/2.d0+c*t; scale2=c/2.d0-b*t
        rcx=scale1*v1x+scale2*nbx; rcy=scale1*v1y+scale2*nby; rcz=scale1*v1z+scale2*nbz
        radius=dsqrt(rcx*rcx+rcy*rcy+rcz*rcz)
        if(radius>0.d0)then
          curvature=1.d0/radius; rcx=rcx/radius; rcy=rcy/radius; rcz=rcz/radius
        endif
      endif
    endif
    factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(lup))+ &
     (dsqrt(yve(ipoint-1))/dsqrt(ldown)))**2.d0
    factor5=factor3*lup*curvature*(veltangent**2.d0)
    factor4=(li/jetms(ipoint))*factor5
    fvx(j)=fvx(j)-delta*factor4*rcx
    fvy(j)=fvy(j)-delta*factor4*rcy
    fvz(j)=fvz(j)-delta*factor4*rcz
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_evap_force_correction

 subroutine accelerator_maxwell_evap_stage(firstpoint,lastpoint,npjet, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf,jetms,jetch,jetfr, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,linserting,linserted,liniperturb,lairdrag, &
   lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
   att,fveparam,gr,ks,li,vfield,velext,stochastic_model,noisefric, &
   evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,nfieldtype
  logical, intent(in) :: linserting,linserted,liniperturb,lairdrag,lflorentz,luppot
  logical, intent(in) :: stochastic_model
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yvl(0:),yve(0:)
  double precision, intent(in) :: ycf(0:,1:),jetms(0:),jetch(0:)
  logical, intent(in) :: jetfr(0:)
  double precision, intent(out) :: fxx(0:),fyy(0:),fzz(0:),fst(0:)
  double precision, intent(out) :: fvx(0:),fvy(0:),fvz(0:),fve(0:)
  double precision, intent(in) :: pfreq,consistency,findex,yieldstress
  double precision, intent(in) :: att,fveparam,gr,ks,li,vfield,velext,noisefric
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity
  double precision, intent(in) :: cp0,Bev,mev,tev
 logical :: ok

  ! The force stage is the already validated serial 3-D accelerator path.
  ! The combined Maxwell kernel then replaces its Newtonian stress and adds
  ! the evaporation derivative: the stage of every Maxwell evaporative
  ! device step (RK and Platen, device_stage).  Once beads have been removed, eom4_ev gives the lead
  ! bead the surface-tension and lift terms computed with the last collected
  ! bead; collector_curvature reproduces them (missing until 2026-09-30).
  ! Since 2026-10-05 the evaporation rate and the Maxwell stress are
  ! computed in the same kernel (fev_evap), not in a second one.
  ok=accelerator_eom3_stage(firstpoint,lastpoint,npjet,yxx,yyy,yzz,yst, &
   yvx,yvy,yvz,yvl,ycf,jetms,jetch,jetfr,fxx,fyy,fzz,fst,fvx,fvy,fvz, &
   linserted,liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
   consistency,findex,yieldstress,att,fveparam,gr,ks,li,vfield,velext, &
   stochastic_model,noisefric,yve,collector_curvature=.true., &
   fev_evap=fve,ev_airv=evairv,ev_masscoeff=evmasscoeff,ev_sqrevsc=sqrevsc, &
   ev_csvapour=evcsvapour,ev_umidity=evumidity,ev_cp0=cp0,ev_bev=Bev, &
   ev_mev=mev,ev_tev=tev)
  if(.not.ok)return
 end subroutine accelerator_maxwell_evap_stage

 subroutine accelerator_kv_stage(firstpoint,lastpoint,npjet, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf,jetms,jetch,jetfr, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,linserting,linserted,liniperturb,lairdrag, &
   lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
   att,fveparam,gr,ks,li,vfield,velext,evairv,evmasscoeff,sqrevsc, &
   evcsvapour,evumidity,cp0,Bev,mev,tev,evlim,evaporative)
! Kelvin-Voigt stage, with evaporation (the default) or without
! (evaporative=.false.: yve and fve not referenced, as in eom3_KV_pos_v and
! eom3_KV_st).  Until 2026-10-06 accelerator_kv_evap_stage, evaporative only.
  implicit none
  logical, intent(in), optional :: evaporative
  integer, intent(in) :: firstpoint,lastpoint,npjet,nfieldtype
  logical, intent(in) :: linserting,linserted,liniperturb,lairdrag
  logical, intent(in) :: lflorentz,luppot,jetfr(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yvl(0:),yve(0:)
  double precision, intent(in) :: ycf(0:,1:),jetms(0:),jetch(0:)
  double precision, intent(out) :: fxx(0:),fyy(0:),fzz(0:),fst(0:)
  double precision, intent(out) :: fvx(0:),fvy(0:),fvz(0:),fve(0:)
  double precision, intent(in) :: pfreq,consistency,findex,yieldstress
  double precision, intent(in) :: att,fveparam,gr,ks,li,vfield,velext
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc
  double precision, intent(in) :: evcsvapour,evumidity,cp0,Bev,mev,tev,evlim
  logical :: ok,evap

  evap=.true.
  if(present(evaporative))evap=evaporative
  ! Reproduce the established CPU eom3_KV_pos_v(_ev) semantics.  Those
  ! routines do not add air drag or lift, even when airdrag is enabled in
  ! the input.
  if(evap)then
    ok=accelerator_eom3_stage(firstpoint,lastpoint,npjet,yxx,yyy,yzz,yst, &
     yvx,yvy,yvz,yvl,ycf,jetms,jetch,jetfr,fxx,fyy,fzz,fst,fvx,fvy,fvz, &
     linserted,liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
     consistency,findex,yieldstress,att,fveparam,gr,ks,li,vfield,velext, &
     .false.,0.d0,yve,apply_airdrag=.false.,collector_curvature=.true.)
  else
    ok=accelerator_eom3_stage(firstpoint,lastpoint,npjet,yxx,yyy,yzz,yst, &
     yvx,yvy,yvz,yvl,ycf,jetms,jetch,jetfr,fxx,fyy,fzz,fst,fvx,fvy,fvz, &
     linserted,liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
     consistency,findex,yieldstress,att,fveparam,gr,ks,li,vfield,velext, &
     .false.,0.d0,apply_airdrag=.false.,collector_curvature=.true.)
  endif
  if(.not.ok)return

  call accelerator_kv_stress_3d(firstpoint,lastpoint,npjet, &
   linserting,linserted,jetfr,fve,fst,yxx,yyy,yzz,yvx,yvy,yvz, &
   fvx,fvy,fvz,yst,yvl,yve,evairv,evmasscoeff,sqrevsc,evcsvapour, &
   evumidity,cp0,Bev,mev,tev,evlim,acceleration_chunked=.true., &
   evaporative=evap)
 end subroutine accelerator_kv_stage

 subroutine accelerator_maxwell_rk4_stage_update(firstpoint,lastpoint,h,stage, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim, &
   evaporative)
! evaporative=.false. (device_step_mod, 2026-10-06): no evaporated volume;
! jetve, fev and yev are then not referenced.
  implicit none
  logical, intent(in), optional :: evaporative
  integer, intent(in) :: firstpoint,lastpoint,stage
  double precision, intent(in) :: h,evlim
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:),jetve(0:),jetvl(0:)
  double precision, intent(in) :: fxx(0:),fyy(0:),fzz(0:),fst(0:)
  double precision, intent(in) :: fvx(0:),fvy(0:),fvz(0:),fev(0:)
  double precision, intent(inout) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(inout) :: yvx(0:),yvy(0:),yvz(0:),yev(0:)
  integer :: ipoint,j
  double precision :: scale,ve
  logical :: evap
  evap=.true.
  if(present(evaporative))evap=evaporative
  if(stage==1 .or. stage==2)then
    scale=0.5d0*h
  else
    scale=h
  endif
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
!$acc& jetve,jetvl) present_or_copyout(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev) &
!$acc& present_or_copyin(fxx,fyy,fzz,fst,fvx,fvy,fvz,fev) &
!$acc& private(j,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    yxx(ipoint)=jetxx(ipoint)+scale*fxx(j)
    yyy(ipoint)=jetyy(ipoint)+scale*fyy(j)
    yzz(ipoint)=jetzz(ipoint)+scale*fzz(j)
    yst(ipoint)=jetst(ipoint)+scale*fst(j)
    yvx(ipoint)=jetvx(ipoint)+scale*fvx(j)
    yvy(ipoint)=jetvy(ipoint)+scale*fvy(j)
    yvz(ipoint)=jetvz(ipoint)+scale*fvz(j)
    if(evap)then
      ve=jetve(ipoint)+scale*fev(j)
      if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
      yev(ipoint)=ve
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_rk4_stage_update

 subroutine accelerator_evap_rk2_final_update(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
   f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim,evaporative)
  implicit none
  logical, intent(in), optional :: evaporative
  integer, intent(in) :: firstpoint,lastpoint
  double precision, intent(in) :: h,evlim
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: jetve(0:),jetvl(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1st(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:),f1ev(0:)
  double precision, intent(in) :: f2xx(0:),f2yy(0:),f2zz(0:),f2st(0:)
  double precision, intent(in) :: f2vx(0:),f2vy(0:),f2vz(0:),f2ev(0:)
  double precision, intent(inout) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(inout) :: yvx(0:),yvy(0:),yvz(0:),yev(0:)
  integer :: ipoint,j
  double precision :: scale,ve
  logical :: evap
  evap=.true.
  if(present(evaporative))evap=evaporative
  scale=0.5d0*h
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,jetve,jetvl,f1xx,f1yy,f1zz,f1st,f1vx,f1vy, &
!$acc& f1vz,f1ev,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev) &
!$acc& present_or_copyout(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev) private(j,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    yxx(ipoint)=jetxx(ipoint)+scale*(f1xx(j)+f2xx(j))
    yyy(ipoint)=jetyy(ipoint)+scale*(f1yy(j)+f2yy(j))
    yzz(ipoint)=jetzz(ipoint)+scale*(f1zz(j)+f2zz(j))
    yst(ipoint)=jetst(ipoint)+scale*(f1st(j)+f2st(j))
    yvx(ipoint)=jetvx(ipoint)+scale*(f1vx(j)+f2vx(j))
    yvy(ipoint)=jetvy(ipoint)+scale*(f1vy(j)+f2vy(j))
    yvz(ipoint)=jetvz(ipoint)+scale*(f1vz(j)+f2vz(j))
    if(evap)then
      ve=jetve(ipoint)+scale*(f1ev(j)+f2ev(j))
      if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
      yev(ipoint)=ve
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_evap_rk2_final_update

 subroutine accelerator_maxwell_rk4_final_update(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
   f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
   f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev, &
   f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim,evaporative)
  implicit none
  logical, intent(in), optional :: evaporative
  integer, intent(in) :: firstpoint,lastpoint
  double precision, intent(in) :: h,evlim
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:),jetve(0:),jetvl(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1st(0:),f1vx(0:),f1vy(0:),f1vz(0:),f1ev(0:)
  double precision, intent(in) :: f2xx(0:),f2yy(0:),f2zz(0:),f2st(0:),f2vx(0:),f2vy(0:),f2vz(0:),f2ev(0:)
  double precision, intent(in) :: f3xx(0:),f3yy(0:),f3zz(0:),f3st(0:),f3vx(0:),f3vy(0:),f3vz(0:),f3ev(0:)
  double precision, intent(in) :: f4xx(0:),f4yy(0:),f4zz(0:),f4st(0:),f4vx(0:),f4vy(0:),f4vz(0:),f4ev(0:)
  double precision, intent(inout) :: yxx(0:),yyy(0:),yzz(0:),yst(0:),yvx(0:),yvy(0:),yvz(0:),yev(0:)
  integer :: ipoint,j
  double precision :: scale,ve
  logical :: evap
  evap=.true.
  if(present(evaporative))evap=evaporative
  scale=h/6.d0
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl) &
!$acc& present_or_copyout(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev) &
!$acc& present_or_copyin(f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev,f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
!$acc& f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz,f3ev,f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz,f4ev) private(j,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    yxx(ipoint)=jetxx(ipoint)+scale*(f1xx(j)+2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
    yyy(ipoint)=jetyy(ipoint)+scale*(f1yy(j)+2.d0*(f2yy(j)+f3yy(j))+f4yy(j))
    yzz(ipoint)=jetzz(ipoint)+scale*(f1zz(j)+2.d0*(f2zz(j)+f3zz(j))+f4zz(j))
    yst(ipoint)=jetst(ipoint)+scale*(f1st(j)+2.d0*(f2st(j)+f3st(j))+f4st(j))
    yvx(ipoint)=jetvx(ipoint)+scale*(f1vx(j)+2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
    yvy(ipoint)=jetvy(ipoint)+scale*(f1vy(j)+2.d0*(f2vy(j)+f3vy(j))+f4vy(j))
    yvz(ipoint)=jetvz(ipoint)+scale*(f1vz(j)+2.d0*(f2vz(j)+f3vz(j))+f4vz(j))
    if(evap)then
      ve=jetve(ipoint)+scale*(f1ev(j)+2.d0*(f2ev(j)+f3ev(j))+f4ev(j))
      if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
      yev(ipoint)=ve
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_rk4_final_update

 subroutine accelerator_maxwell_commit_state(firstpoint,lastpoint, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve, &
   counterlpath,ncounterlpath,maxstress,maxstressposx,evaporative)
  implicit none
  logical, intent(in), optional :: evaporative
  integer, intent(in) :: firstpoint,lastpoint
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yev(0:)
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:),jetve(0:)
  integer, intent(inout) :: ncounterlpath
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint
  double precision :: dx,dy,dz
  logical :: evap
  evap=.true.
  if(present(evaporative))evap=evaporative
  call accelerator_map_statistics(counterlpath,ncounterlpath,maxstress, &
   maxstressposx)
#ifdef _OPENACC
!$acc parallel loop gang vector present(yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
!$acc& jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve, &
!$acc& counterlpath,statistics_step_max) private(dx,dy,dz) &
!$acc& reduction(+:counterlpath) reduction(max:statistics_step_max)
#endif
  do ipoint=firstpoint,lastpoint
    if(ipoint<lastpoint)then
      dx=yxx(ipoint)-yxx(ipoint+1)
      dy=yyy(ipoint)-yyy(ipoint+1)
      dz=yzz(ipoint)-yzz(ipoint+1)
      counterlpath=counterlpath+dsqrt(dx*dx+dy*dy+dz*dz)
    endif
    statistics_step_max=max(statistics_step_max,yst(ipoint))
    jetxx(ipoint)=yxx(ipoint); jetyy(ipoint)=yyy(ipoint); jetzz(ipoint)=yzz(ipoint)
    jetst(ipoint)=yst(ipoint); jetvx(ipoint)=yvx(ipoint); jetvy(ipoint)=yvy(ipoint)
    jetvz(ipoint)=yvz(ipoint)
    if(evap)jetve(ipoint)=yev(ipoint)
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_commit_state

 subroutine accelerator_compute_posnoinserted_3d(npjet,linserted,resolution,yxx,yyy,yzz)
  implicit none
  integer, intent(in) :: npjet
  logical, intent(in) :: linserted
  double precision, intent(in) :: resolution
  double precision, intent(inout) :: yxx(0:),yyy(0:),yzz(0:)
  if(linserted)return
#ifdef _OPENACC
!$acc serial async(accelerator_queue) present(yxx,yyy,yzz)
#endif
  call device_place_inserting_bead(npjet,resolution,yxx,yyy,yzz)
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_compute_posnoinserted_3d

 integer function accelerator_contact_bead(inpjet,npjet,jetfr)
! The last bead of the leading run of beads frozen at the collector (the
! contact bead between the deposited fiber and the free jet), or inpjet:
! the refinement of a run without removal starts there.  One reduction and
! one integer back to the host.
  implicit none
  integer, intent(in) :: inpjet,npjet
  logical, intent(in) :: jetfr(0:)
  integer :: ipoint,firstfree
  firstfree=npjet
#ifdef _OPENACC
  call accelerator_wait()
!$acc parallel loop reduction(min:firstfree) present(jetfr)
#endif
  do ipoint=inpjet,npjet-1
    if(.not.jetfr(ipoint))firstfree=min(firstfree,ipoint)
  enddo
  accelerator_contact_bead=max(inpjet,firstfree-1)
 end function accelerator_contact_bead

 subroutine accelerator_freeze_at_collector(firstpoint,lastpoint,h,jetxx,jetfr)
! remove_jetbead's freezing for the device step of a run without insertion
! (a fixed bead set, which removes nothing on the device): a bead that
! reaches the collector is held there and leaves the Coulomb sums.  Runs
! with insertion freeze in the topology step.
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  double precision, intent(in) :: h
  double precision, intent(inout) :: jetxx(0:)
  logical, intent(inout) :: jetfr(0:)
  integer :: ipoint
#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) present(jetxx,jetfr)
#endif
  do ipoint=firstpoint,lastpoint
    if(jetxx(ipoint)>=h)then
      jetfr(ipoint)=.true.
      jetxx(ipoint)=h
    endif
  enddo
 end subroutine accelerator_freeze_at_collector

 subroutine device_place_inserting_bead(npjet,resolution,yxx,yyy,yzz)
! compute_posnoinserted on one state: the blocked nozzle bead npjet-1 is put
! at resolution from the nozzle bead, towards bead npjet-2.
#ifdef _OPENACC
!$acc routine seq
#endif
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(in) :: resolution
  double precision, intent(inout) :: yxx(0:),yyy(0:),yzz(0:)
  double precision :: dx,dy,dz,distance,scale
  dx=yxx(npjet-2)-yxx(npjet)
  dy=yyy(npjet-2)-yyy(npjet)
  dz=yzz(npjet-2)-yzz(npjet)
  distance=dsqrt(dx*dx+dy*dy+dz*dz)
  if(distance>0.d0)then
    scale=resolution/distance
    yxx(npjet-1)=yxx(npjet)+scale*dx
    yyy(npjet-1)=yyy(npjet)+scale*dy
    yzz(npjet-1)=yzz(npjet)+scale*dz
  endif
 end subroutine device_place_inserting_bead

 subroutine accelerator_smooth_charge_3d(npjet,linserted,thresolution,dresolution, &
   yxx,yyy,yzz,jetch)
  implicit none
  integer, intent(in) :: npjet
  logical, intent(in) :: linserted
  double precision, intent(in) :: thresolution,dresolution
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(inout) :: jetch(0:)
  if(linserted)return
#ifdef _OPENACC
!$acc serial async(accelerator_queue) present(yxx,yyy,yzz,jetch)
#endif
  call device_smooth_charge(npjet,thresolution,dresolution,yxx,yyy,yzz,jetch)
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_smooth_charge_3d

 subroutine device_smooth_charge(npjet,thresolution,dresolution,yxx,yyy,yzz, &
   jetch)
! smooth_charge on one state: keep the charge of the blocked nozzle bead in
! accelerator_smoothed_charge and scale it by the cut-off of its distance.
#ifdef _OPENACC
!$acc routine seq
#endif
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(in) :: thresolution,dresolution
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(inout) :: jetch(0:)
  double precision :: dx,dy,dz,distance,factor
  dx=yxx(npjet-2)-yxx(npjet)
  dy=yyy(npjet-2)-yyy(npjet)
  dz=yzz(npjet-2)-yzz(npjet)
  distance=dsqrt(dx*dx+dy*dy+dz*dz)
  if(distance<thresolution)then
    factor=0.d0
  elseif(distance>=dresolution)then
    factor=1.d0
  else
    factor=-0.5d0*dcos((distance-thresolution)* &
     3.14159265358979323846d0/(dresolution-thresolution))+0.5d0
  endif
  accelerator_smoothed_charge=jetch(npjet-1)
  jetch(npjet-1)=factor*accelerator_smoothed_charge
 end subroutine device_smooth_charge

 subroutine accelerator_platen_stage_prep(npjet,restore,thresolution, &
   dresolution,resolution,yxx,yyy,yzz,jetch)
! One serial kernel before a force evaluation of the fused persistent Platen
! step (accelerator_platen_end_step): restore the charge smoothed for
! the previous evaluation when restore is set, smooth it for this one, and
! place the inserting bead; three kernels before 2026-10-05.  The kernels in
! between do not read jetch.  Only for a blocked nozzle bead.
  implicit none
  integer, intent(in) :: npjet
  logical, intent(in) :: restore
  double precision, intent(in) :: thresolution,dresolution,resolution
  double precision, intent(inout) :: yxx(0:),yyy(0:),yzz(0:),jetch(0:)
#ifdef _OPENACC
!$acc serial async(accelerator_queue) present(yxx,yyy,yzz,jetch)
#endif
  if(restore)jetch(npjet-1)=accelerator_smoothed_charge
  call device_smooth_charge(npjet,thresolution,dresolution,yxx,yyy,yzz,jetch)
  call device_place_inserting_bead(npjet,resolution,yxx,yyy,yzz)
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_platen_stage_prep

 subroutine accelerator_restore_charge(npjet,linserted,jetch)
  implicit none
  integer, intent(in) :: npjet
  logical, intent(in) :: linserted
  double precision, intent(inout) :: jetch(0:)
  if(linserted)return
#ifdef _OPENACC
!$acc serial async(accelerator_queue) present(jetch)
#endif
  jetch(npjet-1)=accelerator_smoothed_charge
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_restore_charge

 subroutine accelerator_evaporation_force_3d(rheology_mode,firstpoint,lastpoint,npjet, &
   linserting,linserted,jetfr,fev,yxx,yyy,yzz,yvx,yvy,yvz,yve,evairv, &
   evmasscoeff,sqrevsc,evcsvapour,evumidity)
  implicit none
  ! This kernel is intentionally rheology-independent. Maxwell and
  ! Kelvin--Voigt stages must call the same evaporation-rate implementation;
  ! only their constitutive stress update belongs in the rheology branch.
  integer, intent(in) :: rheology_mode,firstpoint,lastpoint,npjet
  logical, intent(in) :: linserting,linserted,jetfr(0:)
  double precision, intent(out) :: fev(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yve(0:)
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc
  double precision, intent(in) :: evcsvapour,evumidity
  integer :: ipoint,j
  double precision :: dx,dy,dz,beadlen,vnorm,re
#ifdef _OPENACC
!$acc parallel loop gang vector copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yvx(0:npjet),yvy(0:npjet),yvz(0:npjet), &
!$acc& yve(0:npjet),jetfr(0:npjet)) copyout(fev(0:lastpoint-firstpoint)) &
!$acc& private(j,dx,dy,dz,beadlen,vnorm,re)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    fev(j)=0.d0
    if(jetfr(ipoint))cycle
    if(ipoint>=npjet)cycle
    if(ipoint==npjet-1 .and. .not.linserted)cycle
    if(ipoint==npjet-2 .and. .not.linserted)then
      dx=yxx(npjet)-yxx(ipoint)
      dy=yyy(npjet)-yyy(ipoint)
      dz=yzz(npjet)-yzz(ipoint)
    else
      dx=yxx(ipoint+1)-yxx(ipoint)
      dy=yyy(ipoint+1)-yyy(ipoint)
      dz=yzz(ipoint+1)-yzz(ipoint)
    endif
    beadlen=dsqrt(dx*dx+dy*dy+dz*dz)
    vnorm=dsqrt(yvx(ipoint)*yvx(ipoint)+yvy(ipoint)*yvy(ipoint)+ &
     yvz(ipoint)*yvz(ipoint))
    if(beadlen>0.d0 .and. evairv>0.d0 .and. yve(ipoint)>0.d0)then
      re=(2.d0*dsqrt(yve(ipoint)/(3.14159265358979323846d0*beadlen))* &
       vnorm)/evairv
      fev(j)=-evmasscoeff*0.495d0*(re**(1.d0/3.d0))*sqrevsc* &
       (evcsvapour-evumidity)*3.14159265358979323846d0*beadlen
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_evaporation_force_3d

 subroutine device_maxwell_evap_stress_point(ipoint,j,npjet,linserted,jetfr, &
   fev,fstval,yxx,yyy,yzz,yvx,yvy,yvz,yst,yvl,yve,evairv,evmasscoeff,sqrevsc, &
   evcsvapour,evumidity,cp0,Bev,mev,tev,consistency,findex,yieldstress)
! Evaporation rate fev(j) and Maxwell stress derivative of bead ipoint (Yarin
! law), the body of accelerator_maxwell_evap_stress_3d; accelerator_eom3_stage
! evaluates it inside its own loop when asked to (fev_evap).
#ifdef _OPENACC
!$acc routine seq
#endif
  implicit none
  integer, intent(in) :: ipoint,j,npjet
  logical, intent(in) :: linserted,jetfr(0:)
  double precision, intent(inout) :: fev(0:)
  double precision, intent(out) :: fstval
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yvx(0:),yvy(0:),yvz(0:)
  double precision, intent(in) :: yst(0:),yvl(0:),yve(0:)
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity
  double precision, intent(in) :: cp0,Bev,mev,tev,consistency,findex,yieldstress
  double precision :: dx,dy,dz,beadlen,vnorm,re,beadvel,cp,ratmu,rattao
  fev(j)=0.d0
  fstval=0.d0
  if(jetfr(ipoint) .or. ipoint>=npjet)return
  if(ipoint==npjet-1 .and. .not.linserted)return
  if(ipoint==npjet-2 .and. .not.linserted)then
    dx=yxx(npjet)-yxx(ipoint); dy=yyy(npjet)-yyy(ipoint); dz=yzz(npjet)-yzz(ipoint)
  else
    dx=yxx(ipoint+1)-yxx(ipoint); dy=yyy(ipoint+1)-yyy(ipoint); dz=yzz(ipoint+1)-yzz(ipoint)
  endif
  beadlen=dsqrt(dx*dx+dy*dy+dz*dz)
  vnorm=dsqrt(yvx(ipoint)*yvx(ipoint)+yvy(ipoint)*yvy(ipoint)+yvz(ipoint)*yvz(ipoint))
  if(beadlen>0.d0 .and. evairv>0.d0 .and. yve(ipoint)>0.d0)then
    re=(2.d0*dsqrt(yve(ipoint)/(3.14159265358979323846d0*beadlen))*vnorm)/evairv
    fev(j)=-evmasscoeff*0.495d0*(re**(1.d0/3.d0))*sqrevsc* &
     (evcsvapour-evumidity)*3.14159265358979323846d0*beadlen
  endif
  if(yve(ipoint)<=0.d0 .or. yvl(ipoint)<=0.d0)return
  if(ipoint==npjet-2 .and. .not.linserted)then
    dx=yxx(ipoint)-yxx(npjet); dy=yyy(ipoint)-yyy(npjet); dz=yzz(ipoint)-yzz(npjet)
    beadvel=((yvx(ipoint)-yvx(npjet))*dx+(yvy(ipoint)-yvy(npjet))*dy+ &
     (yvz(ipoint)-yvz(npjet))*dz)/beadlen
  else
    dx=yxx(ipoint)-yxx(ipoint+1); dy=yyy(ipoint)-yyy(ipoint+1); dz=yzz(ipoint)-yzz(ipoint+1)
    beadvel=((yvx(ipoint)-yvx(ipoint+1))*dx+(yvy(ipoint)-yvy(ipoint+1))*dy+ &
     (yvz(ipoint)-yvz(ipoint+1))*dz)/beadlen
  endif
  cp=cp0*yvl(ipoint)/yve(ipoint)
  rattao=(cp/cp0)**tev
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  fstval=(1.d0/rattao)*(yieldstress+consistency*ratmu* &
   (beadvel/beadlen)**findex-yst(ipoint))
 end subroutine device_maxwell_evap_stress_point

 subroutine accelerator_kv_stress_3d(firstpoint,lastpoint,npjet, &
   linserting,linserted,jetfr,fev,fst,yxx,yyy,yzz,yvx,yvy,yvz,yax,yay,yaz, &
   yst,yvl,yve,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity, &
   cp0,Bev,mev,tev,evlim,acceleration_chunked,evaporative)
! Kelvin-Voigt stress derivative from the stage accelerations: the
! concentration-dependent law of eom3_KV_st_ev with the evaporation rate
! (evaporative, the default), or the plain law of eom3_KV_st (strain rate
! plus its time derivative), with fev zero and yve not referenced.  Until
! 2026-10-06 the evaporative kernel was accelerator_kv_evap_stress_3d, and
! accelerator_kv_stress_3d was an unused copy of it.
  implicit none
  logical, intent(in), optional :: evaporative
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: linserting,linserted,jetfr(0:)
  double precision, intent(out) :: fev(0:),fst(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yvx(0:),yvy(0:),yvz(0:)
  double precision, intent(in) :: yax(0:),yay(0:),yaz(0:),yst(0:),yvl(0:),yve(0:)
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity
  double precision, intent(in) :: cp0,Bev,mev,tev,evlim
  logical, intent(in), optional :: acceleration_chunked
  integer :: ipoint,j,ia,ianext
  double precision :: dx,dy,dz,beadlen,vnorm,re,beadvel,beadacc
  double precision :: cp,ratmu,ratg,dcpdt,dratmu,dratg,strain,strainrate,strainacc
  logical :: chunked_acceleration,evap
  chunked_acceleration=.false.
  if(present(acceleration_chunked))chunked_acceleration=acceleration_chunked
  evap=.true.
  if(present(evaporative))evap=evaporative
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yax,yay,yaz, &
!$acc& yst(0:npjet),yvl(0:npjet),yve(0:npjet),jetfr(0:npjet)) &
!$acc& present_or_copyout(fev(0:lastpoint-firstpoint),fst(0:lastpoint-firstpoint)) &
!$acc& private(j,ia,ianext,dx,dy,dz,beadlen,vnorm,re,beadvel,beadacc,cp,ratmu,ratg, &
!$acc& dcpdt,dratmu,dratg,strain,strainrate,strainacc)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    ia=ipoint
    ianext=ipoint+1
    if(chunked_acceleration)then
      ia=j
      ianext=j+1
    endif
    fev(j)=0.d0
    fst(j)=0.d0
    if(jetfr(ipoint) .or. ipoint>=npjet)cycle
    if(ipoint==npjet-1 .and. .not.linserted)cycle
    if(ipoint==npjet-2 .and. .not.linserted)then
      dx=yxx(npjet)-yxx(ipoint); dy=yyy(npjet)-yyy(ipoint); dz=yzz(npjet)-yzz(ipoint)
    else
      dx=yxx(ipoint+1)-yxx(ipoint); dy=yyy(ipoint+1)-yyy(ipoint); dz=yzz(ipoint+1)-yzz(ipoint)
    endif
    beadlen=dsqrt(dx*dx+dy*dy+dz*dz)
    if(evap)then
      vnorm=dsqrt(yvx(ipoint)*yvx(ipoint)+yvy(ipoint)*yvy(ipoint)+yvz(ipoint)*yvz(ipoint))
      if(beadlen>0.d0 .and. evairv>0.d0 .and. yve(ipoint)>0.d0)then
        re=(2.d0*dsqrt(yve(ipoint)/(3.14159265358979323846d0*beadlen))*vnorm)/evairv
        fev(j)=-evmasscoeff*0.495d0*(re**(1.d0/3.d0))*sqrevsc* &
         (evcsvapour-evumidity)*3.14159265358979323846d0*beadlen
      endif
      if(beadlen<=0.d0 .or. yve(ipoint)<=0.d0 .or. yvl(ipoint)<=0.d0)cycle
    else
      if(beadlen<=0.d0)cycle
    endif
    if(ipoint==npjet-2 .and. .not.linserted)then
      dx=yxx(ipoint)-yxx(npjet); dy=yyy(ipoint)-yyy(npjet); dz=yzz(ipoint)-yzz(npjet)
      beadvel=(yvx(ipoint)-yvx(npjet))*dx+(yvy(ipoint)-yvy(npjet))*dy+ &
       (yvz(ipoint)-yvz(npjet))*dz
      if(chunked_acceleration)ianext=npjet-firstpoint
      beadacc=(yax(ia)-yax(ianext))*dx+(yay(ia)-yay(ianext))*dy+ &
       (yaz(ia)-yaz(ianext))*dz
    else
      dx=yxx(ipoint)-yxx(ipoint+1); dy=yyy(ipoint)-yyy(ipoint+1); dz=yzz(ipoint)-yzz(ipoint+1)
      beadvel=(yvx(ipoint)-yvx(ipoint+1))*dx+(yvy(ipoint)-yvy(ipoint+1))*dy+ &
       (yvz(ipoint)-yvz(ipoint+1))*dz
      beadacc=(yax(ia)-yax(ianext))*dx+(yay(ia)-yay(ianext))*dy+ &
       (yaz(ia)-yaz(ianext))*dz
    endif
    strainrate=beadvel/(beadlen*beadlen)
    strainacc=beadacc/(beadlen*beadlen)
    if(.not.evap)then
      fst(j)=strainrate+strainacc
      cycle
    endif
    cp=cp0*yvl(ipoint)/yve(ipoint)
    ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
    ratg=ratmu/(cp/cp0)**tev
    dcpdt=0.d0
    if((yve(ipoint)/yvl(ipoint))>evlim*(1.d0+1.d-12))dcpdt=-cp*fev(j)/yve(ipoint)
    dratmu=ratmu*dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))*dcpdt
    dratg=ratg*(dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))-tev/cp)*dcpdt
    strain=(yst(ipoint)-ratmu*strainrate)/ratg
    fst(j)=ratg*strainrate+ratmu*strainacc+dratg*strain+dratmu*strainrate
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_kv_stress_3d

 subroutine accelerator_coulomb_evap_3d(npjet,inpjet,ycf,yxx,yyy,yzz, &
   yvl,yve,jetms,jetch,jetfr,coulcrossec,q,lmirror,h,ldcutoff,dcutoff, &
   linserting,linserted,icrossec)
  implicit none
  integer, intent(in) :: npjet,inpjet
  double precision, intent(inout) :: ycf(0:,1:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yvl(0:),yve(0:)
  double precision, intent(in) :: jetms(0:),jetch(0:)
  double precision, intent(inout) :: coulcrossec(0:)
  logical, intent(in) :: jetfr(0:),lmirror,ldcutoff,linserting,linserted
  double precision, intent(in) :: q,h,dcutoff,icrossec
  integer :: ipoint,jpoint
  double precision :: dx,dy,dz,norm,cmass1,qt,coef,distance,rmass
  double precision :: sumx,sumy,sumz
  integer :: ihigh
#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yve(0:npjet)) present_or_copyout(coulcrossec(0:npjet)) &
!$acc& private(dx,dy,dz,distance)
#endif
  do ipoint=inpjet,npjet
    if((linserting .and. .not.linserted .and. ipoint>=npjet-2) .or. &
       (linserting .and. linserted .and. ipoint>=npjet-1) .or. &
       (.not.linserting .and. ipoint>=npjet))then
      coulcrossec(ipoint)=icrossec
    else
      dx=yxx(ipoint)-yxx(ipoint+1)
      dy=yyy(ipoint)-yyy(ipoint+1)
      dz=yzz(ipoint)-yzz(ipoint+1)
      distance=dsqrt(dx*dx+dy*dy+dz*dz)
      coulcrossec(ipoint)=dsqrt(yve(ipoint)/ &
       (distance*3.14159265358979323846d0))
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
! One gang per target bead; its vector lanes share the source beads and
! combine their contributions with a reduction.  The former layout ran one
! thread per target over all sources sequentially, which left the device
! almost idle for jets of a few hundred beads.  The summation order differs
! from both the host pair loop and the former kernel only at roundoff.
#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) gang present_or_copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yvl(0:npjet),yve(0:npjet),jetms(0:npjet), &
!$acc& jetch(0:npjet),jetfr(0:npjet),coulcrossec(0:npjet)) &
!$acc& present_or_copyout(ycf(0:npjet,1:3)) private(cmass1,qt,rmass, &
!$acc& sumx,sumy,sumz)
#endif
  do ipoint=inpjet,npjet
    sumx=0.d0
    sumy=0.d0
    sumz=0.d0
    if(.not.jetfr(ipoint))then
      cmass1=yve(ipoint)/yvl(ipoint)
      qt=jetch(ipoint)*q
      rmass=jetms(ipoint)*cmass1
#ifdef _OPENACC
!$acc loop vector reduction(+:sumx,sumy,sumz) private(ihigh,dx,dy,dz,norm,coef)
#endif
      do jpoint=inpjet,npjet
        if(jpoint/=ipoint .and. .not.jetfr(jpoint))then
          dx=yxx(ipoint)-yxx(jpoint)
          dy=yyy(ipoint)-yyy(jpoint)
          dz=yzz(ipoint)-yzz(jpoint)
          norm=dsqrt(dx*dx+dy*dy+dz*dz)
          ! Match the host Maxwell 3D evaporation path: the ordinary
          ! bead-bead interaction does not apply the cutoff. The mirror
          ! interaction below retains the host cutoff behavior.
          if(norm>1.d-30)then
            ! The host pair loop stores the cross-section of the
            ! higher-index bead for both directions.  Use the same
            ! symmetric lookup here.
            ihigh=max(ipoint,jpoint)
            coef=(jetch(jpoint)*qt)/((norm+coulcrossec(ihigh))**2.d0)
            sumx=sumx+coef/rmass*dx/norm
            sumy=sumy+coef/rmass*dy/norm
            sumz=sumz+coef/rmass*dz/norm
          endif
        endif
      enddo
      if(lmirror)then
#ifdef _OPENACC
!$acc loop vector reduction(+:sumx,sumy,sumz) private(dx,dy,dz,norm,coef)
#endif
        do jpoint=inpjet,npjet
          if(.not.jetfr(jpoint))then
            dx=yxx(ipoint)-(dabs(yxx(jpoint)-h)+h)
            dy=yyy(ipoint)-yyy(jpoint)
            dz=yzz(ipoint)-yzz(jpoint)
            norm=dsqrt(dx*dx+dy*dy+dz*dz)
            if(norm>1.d-30 .and. .not.(ldcutoff .and. norm>dcutoff))then
              coef=(jetch(jpoint)*qt)/((norm+coulcrossec(jpoint))**2.d0)
              sumx=sumx-coef/rmass*dx/norm
              sumy=sumy-coef/rmass*dy/norm
              sumz=sumz-coef/rmass*dz/norm
            endif
          endif
        enddo
      endif
    endif
    ycf(ipoint,1)=sumx
    ycf(ipoint,2)=sumy
    ycf(ipoint,3)=sumz
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_coulomb_evap_3d

 subroutine accelerator_release_jet_capacity(mxnpjet,jetxx,jetyy,jetzz,jetst,jetvx, &
   jetvy,jetvz,jetms,jetch,jetvl,jetfr)
  implicit none
  integer, intent(in) :: mxnpjet
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
#ifdef _OPENACC
  call accelerator_wait()
!$acc exit data delete(jetxx(0:mxnpjet),jetyy(0:mxnpjet),jetzz(0:mxnpjet), &
!$acc& jetst(0:mxnpjet),jetvx(0:mxnpjet),jetvy(0:mxnpjet),jetvz(0:mxnpjet), &
!$acc& jetms(0:mxnpjet),jetch(0:mxnpjet),jetvl(0:mxnpjet),jetfr(0:mxnpjet))
#endif
  accelerator_persistent=.false.
  accelerator_topology_enabled=.false.
  accelerator_device_state_authoritative=.false.
 end subroutine accelerator_release_jet_capacity

 subroutine accelerator_release_evaporation_capacity(mxnpjet,jetve,jetce)
  implicit none
  integer, intent(in) :: mxnpjet
  double precision, intent(inout) :: jetve(0:),jetce(0:)
#ifdef _OPENACC
  call accelerator_wait()
!$acc exit data delete(jetve(0:mxnpjet),jetce(0:mxnpjet))
#endif
 end subroutine accelerator_release_evaporation_capacity

 subroutine accelerator_rebind_topology(mxnpjet,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
   jetms,jetch,jetvl,jetfr)
  implicit none
  integer, intent(in) :: mxnpjet
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(in) :: jetfr(0:)
#ifdef _OPENACC
  call accelerator_wait()
!$acc enter data copyin(jetxx(0:mxnpjet),jetyy(0:mxnpjet),jetzz(0:mxnpjet), &
!$acc& jetst(0:mxnpjet),jetvx(0:mxnpjet),jetvy(0:mxnpjet),jetvz(0:mxnpjet), &
!$acc& jetms(0:mxnpjet),jetch(0:mxnpjet),jetvl(0:mxnpjet),jetfr(0:mxnpjet))
#endif
  accelerator_topology_enabled=.true.
 end subroutine accelerator_rebind_topology

 subroutine accelerator_rebind_evaporation(mxnpjet,jetve,jetce)
  implicit none
  integer, intent(in) :: mxnpjet
  double precision, intent(in) :: jetve(0:),jetce(0:)
#ifdef _OPENACC
  call accelerator_wait()
!$acc enter data copyin(jetve(0:mxnpjet),jetce(0:mxnpjet))
#endif
 end subroutine accelerator_rebind_evaporation

 subroutine accelerator_update_device_topology_state(npjet,jetxx,jetyy,jetzz,jetst, &
   jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:),jetms(0:)
  double precision, intent(in) :: jetch(0:),jetvl(0:)
  logical, intent(in) :: jetfr(0:)
  if(.not.accelerator_topology_enabled)return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update device(jetxx(0:npjet),jetyy(0:npjet),jetzz(0:npjet),jetst(0:npjet), &
!$acc& jetvx(0:npjet),jetvy(0:npjet),jetvz(0:npjet),jetms(0:npjet), &
!$acc& jetch(0:npjet),jetvl(0:npjet),jetfr(0:npjet))
#endif
 end subroutine accelerator_update_device_topology_state

 subroutine accelerator_update_device_evaporation_state(npjet,jetve,jetce)
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(in) :: jetve(0:),jetce(0:)
  if(.not.accelerator_topology_enabled)return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update device(jetve(0:npjet),jetce(0:npjet))
#endif
 end subroutine accelerator_update_device_evaporation_state

 subroutine accelerator_refinement_candidate(inpjet,npjet,systype, &
   linserting,linserted,threshold,jetxx,jetyy,jetzz,candidate, &
   refinement_start,path_length,nozzle_correction)
  implicit none
  integer, intent(in) :: inpjet,npjet,systype
  logical, intent(in) :: linserting,linserted
  double precision, intent(in) :: threshold
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:)
  logical, intent(out) :: candidate
  integer, intent(out) :: refinement_start
  double precision, intent(out) :: path_length,nozzle_correction
  integer :: ipoint,lastsegment,scan_start
  double precision :: dx,dy,dz,distance,tdx,tdy,tdz,tail_distance
  double precision :: scan_length,scan_correction

  if(linserting)then
    if(linserted)then
      lastsegment=npjet-2
    else
      lastsegment=npjet-3
    endif
  else
    lastsegment=npjet-1
  endif

  refinement_start=-1
  path_length=0.d0
  nozzle_correction=0.d0
  if(npjet-1>=inpjet)then
#ifdef _OPENACC
! Only the three scan results cross from device to host on an ordinary
! refinement check, in one transfer.  The complete jet state remains resident
! until this scan reports that the historical CPU Akima path may actually
! produce a denser mesh.  One gang reduces the scan and stores the results
! on the device; until 2026-10-01 the three reduction scalars were each
! copied in and out, six transfers per scan.  Only the order of the
! path-length sum differs, and that sum only selects whether to download.
  call accelerator_wait()
!$acc parallel num_gangs(1) vector_length(128) present(jetxx,jetyy,jetzz, &
!$acc& accelerator_refinement_scan) private(scan_start,scan_length,scan_correction)
#endif
    scan_start=-1
    scan_length=0.d0
    scan_correction=0.d0
#ifdef _OPENACC
!$acc loop vector private(dx,dy,dz,distance,tdx,tdy,tdz,tail_distance) &
!$acc& reduction(max:scan_start,scan_correction) reduction(+:scan_length)
#endif
    do ipoint=inpjet,npjet-1
      dx=jetxx(ipoint)-jetxx(ipoint+1)
      if(systype==1)then
        distance=dabs(dx)
      else
        dy=jetyy(ipoint)-jetyy(ipoint+1)
        dz=jetzz(ipoint)-jetzz(ipoint+1)
        distance=dsqrt(dx*dx+dy*dy+dz*dz)
      endif
      scan_length=scan_length+distance
      if(ipoint<=lastsegment .and. distance>threshold) &
       scan_start=max(scan_start,ipoint)

! Match the correction used by check_dynamic_refinement_akima to exclude the
! blocked insertion segment(s) from the requested spline size.
      tail_distance=0.d0
      if(linserting)then
        if(linserted .and. ipoint==npjet-1)then
          tail_distance=distance
        elseif((.not.linserted) .and. ipoint==npjet-2)then
          tdx=jetxx(npjet-2)-jetxx(npjet)
          if(systype==1)then
            tail_distance=dabs(tdx)
          else
            tdy=jetyy(npjet-2)-jetyy(npjet)
            tdz=jetzz(npjet-2)-jetzz(npjet)
            tail_distance=dsqrt(tdx*tdx+tdy*tdy+tdz*tdz)
          endif
        endif
      endif
      scan_correction=max(scan_correction,tail_distance)
    enddo
    accelerator_refinement_scan(1)=dble(scan_start)
    accelerator_refinement_scan(2)=scan_length
    accelerator_refinement_scan(3)=scan_correction
#ifdef _OPENACC
!$acc end parallel
!$acc update self(accelerator_refinement_scan)
#endif
    refinement_start=nint(accelerator_refinement_scan(1))
    path_length=accelerator_refinement_scan(2)
    nozzle_correction=accelerator_refinement_scan(3)
  endif

  candidate=refinement_start>=inpjet
  return
 end subroutine accelerator_refinement_candidate

 subroutine accelerator_platen_predict(firstpoint,lastpoint,h,airamp, &
   noisediff,evlim,jetms,jetvl,jetve,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
   y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y1ev, &
   y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,y2ev,linserted,jetfr,evaporative)
! The two predicted states of a Platen step, with or without evaporation
! (evaporative: jetve, f1ev, y1ev and y2ev are not referenced without it,
! and the mass correction is 1).  Two kernels, accelerator_platen_predict
! and accelerator_platen_evap_predict, until 2026-10-06.
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  logical, intent(in) :: linserted,jetfr(0:),evaporative
  double precision, intent(in) :: h,airamp,noisediff,evlim
  double precision, intent(in) :: jetms(0:),jetvl(0:),jetve(0:)
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1st(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:),f1ev(0:)
  double precision, intent(out) :: y1xx(0:),y1yy(0:),y1zz(0:),y1st(0:)
  double precision, intent(out) :: y1vx(0:),y1vy(0:),y1vz(0:),y1ev(0:)
  double precision, intent(out) :: y2xx(0:),y2yy(0:),y2zz(0:),y2st(0:)
  double precision, intent(out) :: y2vx(0:),y2vy(0:),y2vz(0:),y2ev(0:)
  integer :: ipoint,j
  double precision :: dsqrh,stoc,cmass,ve
  dsqrh=dsqrt(dabs(h))
#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) gang vector present(jetms,jetvl,jetve,jetxx,jetyy, &
!$acc& jetzz,jetst,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy, &
!$acc& f1vz,f1ev,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y1ev,y2xx,y2yy, &
!$acc& y2zz,y2st,y2vx,y2vy,y2vz,y2ev,jetfr) private(j,stoc,cmass,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    cmass=1.d0
    if(evaporative)cmass=jetve(ipoint)/jetvl(ipoint)
    stoc=dsqrt(2.d0*(airamp/(jetms(ipoint)*cmass)+noisediff))
! eom4(_ev) gives no stochastic force to the nozzle, to a collected (frozen)
! bead, or to the bead still being inserted at the nozzle.
    if(ipoint==lastpoint .or. jetfr(ipoint) .or. &
     (ipoint==lastpoint-1 .and. .not.linserted))stoc=0.d0
    y1xx(ipoint)=jetxx(ipoint)+h*f1xx(j)
    y1yy(ipoint)=jetyy(ipoint)+h*f1yy(j)
    y1zz(ipoint)=jetzz(ipoint)+h*f1zz(j)
    y1st(ipoint)=jetst(ipoint)+h*f1st(j)
    y1vx(ipoint)=jetvx(ipoint)+h*f1vx(j)+dsqrh*stoc
    y1vy(ipoint)=jetvy(ipoint)+h*f1vy(j)+dsqrh*stoc
    y1vz(ipoint)=jetvz(ipoint)+h*f1vz(j)+dsqrh*stoc
    y2xx(ipoint)=y1xx(ipoint)
    y2yy(ipoint)=y1yy(ipoint)
    y2zz(ipoint)=y1zz(ipoint)
    y2st(ipoint)=y1st(ipoint)
    y2vx(ipoint)=jetvx(ipoint)+h*f1vx(j)-dsqrh*stoc
    y2vy(ipoint)=jetvy(ipoint)+h*f1vy(j)-dsqrh*stoc
    y2vz(ipoint)=jetvz(ipoint)+h*f1vz(j)-dsqrh*stoc
    if(evaporative)then
      ve=jetve(ipoint)+h*f1ev(j)
      if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
      y1ev(ipoint)=ve
      y2ev(ipoint)=ve
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_predict

 subroutine accelerator_platen_end_step(firstpoint,lastpoint,dt, &
   npjet,mxnpjet,inpjet,linserted,lremove,h,resolution,dresolution, &
   thresolution,ivelocity,istress,imassa,icharge,ivolume,jetxx,jetyy,jetzz, &
   jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetve,jetfr,f1st,yst, &
   evaporative,consistency,findex,yieldstress,cp0,Bev,mev,tev, &
   counterlpath,ncounterlpath,maxstress,maxstressposx,linserting,restore)
! The end of every Platen device step (device_platen_step), evaporative or
! not, in one single-gang kernel, which a few hundred beads keep busy:
! restore the smoothed nozzle charge (restore: not in the oracle builds,
! which smooth on the host), update the stress with the path-length and
! maximum statistics, and, with insertion (linserting), place the inserting
! bead, decide the topology and freeze beads at the collector
! (accelerator_topology_check); then store the step's statistics over the
! active range the topology leaves (accelerator_store_statistics).  A fixed
! bead set makes no topology decision here; device_platen_step freezes its
! beads at the collector (accelerator_freeze_at_collector).  Until
! 2026-10-05 these were eight or nine kernels.  accelerator_topology_check
! then only reads the flags back.  Only the order of the path-length sum
! changes.
! Since 2026-10-06 the stress loop also evaluates the stress derivative at
! the end of the step from the new positions and velocities and the
! predicted stress yst: the Maxwell law of the evaporating jet (as
! xpsys_stress_ev; one kernel before) or, without
! evaporation, the arithmetic of accelerator_eom3_stage (a fourth EOM stage
! before, which computed all the derivatives for this one).  Without
! evaporation jetve is not referenced.
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,mxnpjet,inpjet
  logical, intent(in) :: linserted,lremove
  integer, intent(inout) :: ncounterlpath
  double precision, intent(in) :: dt,h,resolution,dresolution,thresolution
  double precision, intent(in) :: ivelocity,istress,imassa,icharge,ivolume
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  double precision, intent(in) :: jetve(0:)
  logical, intent(inout) :: jetfr(0:)
  double precision, intent(in) :: f1st(0:),yst(0:)
  logical, intent(in) :: evaporative,linserting,restore
  double precision, intent(in) :: consistency,findex,yieldstress
  double precision, intent(in) :: cp0,Bev,mev,tev
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint,j,first,last,lastfreeze,idx
  double precision :: newst,dx,dy,dz,lpsum,stmax
  double precision :: fst2,dxu,dyu,dzu,lup,tux,tuy,tuz,beadvel
  double precision :: beadlen,cp,ratmu,rattao
  call accelerator_map_statistics(counterlpath,ncounterlpath,maxstress,maxstressposx)
  lastfreeze=min(npjet+1,mxnpjet)
#ifdef _OPENACC
!$acc parallel num_gangs(1) vector_length(256) async(accelerator_queue) &
!$acc& present(jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl, &
!$acc& jetve,jetfr,f1st,yst,counterlpath,ncounterlpath,maxstress,maxstressposx, &
!$acc& statistics_step_max,statistics_step_index,accelerator_topology_flags) &
!$acc& private(lpsum,stmax,first,last,idx)
#endif
! restore: the step's force evaluations smoothed the nozzle charge on the
! device (accelerator_platen_stage_prep), which saved the charge to restore;
! the oracle builds smooth and restore it on the host instead.
  if(restore .and. .not.linserted)jetch(npjet-1)=accelerator_smoothed_charge
  lpsum=0.d0
  stmax=-huge(0.d0)
#ifdef _OPENACC
!$acc loop vector reduction(+:lpsum) reduction(max:stmax) &
!$acc& private(j,newst,dx,dy,dz,fst2,dxu,dyu,dzu,lup,tux,tuy,tuz,beadvel, &
!$acc& beadlen,cp,ratmu,rattao)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
! Zero for a frozen bead, the blocked nozzle bead and the last bead; the
! bead before a blocked one is joined to the last.
    fst2=0.d0
    if(.not.jetfr(ipoint) .and. ipoint<npjet .and. &
     .not.(ipoint==npjet-1 .and. .not.linserted))then
      if(ipoint==npjet-2 .and. .not.linserted)then
        dxu=jetxx(ipoint)-jetxx(npjet)
        dyu=jetyy(ipoint)-jetyy(npjet)
        dzu=jetzz(ipoint)-jetzz(npjet)
      else
        dxu=jetxx(ipoint)-jetxx(ipoint+1)
        dyu=jetyy(ipoint)-jetyy(ipoint+1)
        dzu=jetzz(ipoint)-jetzz(ipoint+1)
      endif
      if(evaporative)then
! The Maxwell law of the evaporating jet, as in xpsys_stress_ev.
        beadlen=dsqrt(dxu*dxu+dyu*dyu+dzu*dzu)
        if(beadlen>0.d0 .and. jetve(ipoint)>0.d0 .and. jetvl(ipoint)>0.d0)then
          if(ipoint==npjet-2 .and. .not.linserted)then
            beadvel=((jetvx(ipoint)-jetvx(npjet))*dxu+ &
             (jetvy(ipoint)-jetvy(npjet))*dyu+ &
             (jetvz(ipoint)-jetvz(npjet))*dzu)/beadlen
          else
            beadvel=((jetvx(ipoint)-jetvx(ipoint+1))*dxu+ &
             (jetvy(ipoint)-jetvy(ipoint+1))*dyu+ &
             (jetvz(ipoint)-jetvz(ipoint+1))*dzu)/beadlen
          endif
          cp=cp0*jetvl(ipoint)/jetve(ipoint)
          rattao=(cp/cp0)**tev
          ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
          fst2=(1.d0/rattao)*(yieldstress+consistency*ratmu* &
           (beadvel/beadlen)**findex-yst(ipoint))
        endif
      else
! accelerator_eom3_stage.
        lup=dsqrt(dxu*dxu+dyu*dyu+dzu*dzu)
        tux=dxu/lup
        tuy=dyu/lup
        tuz=dzu/lup
        if(ipoint==npjet-2 .and. .not.linserted)then
          beadvel=(jetvx(ipoint)-jetvx(npjet))*tux+ &
           (jetvy(ipoint)-jetvy(npjet))*tuy+ &
           (jetvz(ipoint)-jetvz(npjet))*tuz
        else
          beadvel=(jetvx(ipoint)-jetvx(ipoint+1))*tux+ &
           (jetvy(ipoint)-jetvy(ipoint+1))*tuy+ &
           (jetvz(ipoint)-jetvz(ipoint+1))*tuz
        endif
        fst2=yieldstress+consistency*(beadvel/lup)**findex-yst(ipoint)
      endif
    endif
    newst=jetst(ipoint)+0.5d0*dt*(f1st(j)+fst2)
    if(ipoint<lastpoint)then
      dx=jetxx(ipoint)-jetxx(ipoint+1)
      dy=jetyy(ipoint)-jetyy(ipoint+1)
      dz=jetzz(ipoint)-jetzz(ipoint+1)
      lpsum=lpsum+dsqrt(dx*dx+dy*dy+dz*dz)
    endif
    jetst(ipoint)=newst
    stmax=max(stmax,newst)
  enddo
  counterlpath=counterlpath+lpsum
  statistics_step_max=max(statistics_step_max,stmax)
! Without insertion (a fixed bead set, Tests 12 and 20; since 2026-10-06
! through this kernel too) there is no topology to decide and the
! statistics cover inpjet..npjet.
  first=inpjet
  last=npjet
  if(linserting)then
    if(.not.linserted)call device_place_inserting_bead(npjet,resolution, &
     jetxx,jetyy,jetzz)
    call device_topology_decide(npjet,mxnpjet,inpjet,linserted,lremove,h, &
     resolution,dresolution,thresolution,ivelocity,istress,imassa,icharge, &
     ivolume,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
! Beads at the collector are frozen with or without removal (remove_jetbead).
#ifdef _OPENACC
!$acc loop vector
#endif
    do ipoint=inpjet,lastfreeze
      if(accelerator_topology_flags(4)==0 .and. &
       ipoint<=accelerator_topology_flags(1))then
        if(jetxx(ipoint)>=h)then
          jetfr(ipoint)=.true.
          jetxx(ipoint)=h
        endif
      endif
    enddo
! accelerator_store_statistics over the range left by the topology step.
    first=inpjet+accelerator_topology_flags(5)
    last=accelerator_topology_flags(1)
  endif
  idx=statistics_step_index
#ifdef _OPENACC
!$acc loop vector reduction(max:idx)
#endif
  do ipoint=first,last
    if(jetst(ipoint)==statistics_step_max)idx=max(idx,ipoint)
  enddo
  statistics_step_index=idx
  ncounterlpath=ncounterlpath+1
  if(statistics_step_max>=maxstress)then
    maxstress=statistics_step_max
    maxstressposx=jetxx(statistics_step_index)
  endif
  statistics_step_max=-huge(0.d0)
  statistics_step_index=-1
#ifdef _OPENACC
!$acc end parallel
#endif
  accelerator_topology_decided=linserting
  accelerator_statistics_stored=.true.
 end subroutine accelerator_platen_end_step

 subroutine accelerator_platen_update(firstpoint,lastpoint,dt,npjet, &
   mxnpjet,linserted,evaporative,airamp,noisediff,pfreq,liniperturb,evlim, &
   historybase,historywindow,historyvalues,gaussianhistory,jetxx,jetyy, &
   jetzz,jetvx,jetvy,jetvz,jetms,jetvl,jetve,jetfr,f1xx,f1yy,f1zz,f1vx, &
   f1vy,f1vz,f1ev,f2vx,f2vy,f2vz,f3vx,f3vy,f3vz,y1xx,y1yy,y1zz,y1ev,evairv, &
   evmasscoeff,sqrevsc,evcsvapour,evumidity)
! The per-bead updates after the last force evaluation of a fused persistent
! dynamic Platen step in one kernel (2026-10-06): velocity
! (accelerator_platen[_evap]_velocity), evaporation rate at the predicted
! positions with the new velocity (accelerator_maxwell_evap_stress_3d), and
! positions and evaporated volume (accelerator_platen[_evap]_positions).
! Each bead reads only its own new velocity and the predicted state, so one
! loop does all three; the Maxwell stress at the new state, which reads the
! neighbours, is evaluated by accelerator_platen_end_step.  Three kernels
! with evaporation (plus the stress kernel), two without, until then.
! Explicit-shape dummies: with one descriptor per assumed-shape array, NVHPC
! 24.3 and 25.5 reject the larger single kernel this replaces (invalid NVVM
! IR) and copy descriptor temporaries for arrays passed on to device
! routines on the asynchronous queue after the routine has returned.
! Without evaporation, jetve, f1ev, y1xx, y1yy, y1zz and y1ev are not
! referenced.
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,mxnpjet
  integer, intent(in) :: historybase,historywindow,historyvalues
  logical, intent(in) :: linserted,evaporative,liniperturb
  double precision, intent(in) :: dt,airamp,noisediff,pfreq,evlim
  double precision, intent(in) :: gaussianhistory(0:historyvalues-1)
  double precision, intent(inout) :: jetxx(0:mxnpjet),jetyy(0:mxnpjet)
  double precision, intent(inout) :: jetzz(0:mxnpjet),jetvx(0:mxnpjet)
  double precision, intent(inout) :: jetvy(0:mxnpjet),jetvz(0:mxnpjet)
  double precision, intent(in) :: jetms(0:mxnpjet),jetvl(0:mxnpjet)
  double precision, intent(inout) :: jetve(0:mxnpjet)
  logical, intent(in) :: jetfr(0:mxnpjet)
  double precision, intent(in) :: f1xx(0:lastpoint-firstpoint)
  double precision, intent(in) :: f1yy(0:lastpoint-firstpoint)
  double precision, intent(in) :: f1zz(0:lastpoint-firstpoint)
  double precision, intent(in) :: f1vx(0:lastpoint-firstpoint)
  double precision, intent(in) :: f1vy(0:lastpoint-firstpoint)
  double precision, intent(in) :: f1vz(0:lastpoint-firstpoint)
  double precision, intent(in) :: f1ev(0:lastpoint-firstpoint)
  double precision, intent(in) :: f2vx(0:lastpoint-firstpoint)
  double precision, intent(in) :: f2vy(0:lastpoint-firstpoint)
  double precision, intent(in) :: f2vz(0:lastpoint-firstpoint)
  double precision, intent(in) :: f3vx(0:lastpoint-firstpoint)
  double precision, intent(in) :: f3vy(0:lastpoint-firstpoint)
  double precision, intent(in) :: f3vz(0:lastpoint-firstpoint)
  double precision, intent(in) :: y1xx(0:npjet),y1yy(0:npjet),y1zz(0:npjet)
  double precision, intent(in) :: y1ev(0:npjet)
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc,evcsvapour
  double precision, intent(in) :: evumidity
  integer :: ipoint,j,component,index1,index2,hbase,hstride,hvalues,hfirst
  double precision :: dx,dy,dz,beadlen,vnorm,re,fev
  double precision :: dsqrh,tsqh,prefactor,stoc,cmass,u1,u2,ww,zz
  double precision :: f2x,f2y,f2z,y1y,y1z,ve
  dsqrh=dsqrt(dabs(dt)); tsqh=dsqrh**3.d0; prefactor=0.5d0/dsqrh
  hbase=historybase; hstride=historywindow; hvalues=historyvalues
  hfirst=firstpoint
#ifdef _OPENACC
!$acc parallel loop gang vector async(accelerator_queue) &
!$acc& present(gaussianhistory,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetms, &
!$acc& jetvl,jetve,jetfr,f1xx,f1yy,f1zz,f1vx,f1vy,f1vz,f1ev,f2vx,f2vy, &
!$acc& f2vz,f3vx,f3vy,f3vz,y1xx,y1yy,y1zz,y1ev) &
!$acc& private(j,component,index1,index2,stoc,cmass,u1,u2,ww,zz,dx,dy,dz, &
!$acc& beadlen,vnorm,re,fev,f2x,f2y,f2z,y1y,y1z,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
! Velocity.
    if(evaporative)then
      cmass=jetve(ipoint)/jetvl(ipoint)
      stoc=dsqrt(2.d0*(airamp/(jetms(ipoint)*cmass)+noisediff))
    else
      stoc=dsqrt(2.d0*(airamp/jetms(ipoint)+noisediff))
    endif
    if(ipoint==lastpoint .or. jetfr(ipoint) .or. &
     (ipoint==lastpoint-1 .and. .not.linserted))stoc=0.d0
    component=1
    index1=mod(hbase+(ipoint-hfirst)+hstride*(component-1),hvalues)
    index2=mod(hbase+(ipoint-hfirst)+hstride*(component-1+3),hvalues)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvx(ipoint)=jetvx(ipoint)+stoc*ww+prefactor*(f2vx(j)-f3vx(j))*zz+ &
     0.25d0*dt*(f2vx(j)+2.d0*f1vx(j)+f3vx(j))
    component=2
    index1=mod(hbase+(ipoint-hfirst)+hstride*(component-1),hvalues)
    index2=mod(hbase+(ipoint-hfirst)+hstride*(component-1+3),hvalues)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvy(ipoint)=jetvy(ipoint)+stoc*ww+prefactor*(f2vy(j)-f3vy(j))*zz+ &
     0.25d0*dt*(f2vy(j)+2.d0*f1vy(j)+f3vy(j))
    component=3
    index1=mod(hbase+(ipoint-hfirst)+hstride*(component-1),hvalues)
    index2=mod(hbase+(ipoint-hfirst)+hstride*(component-1+3),hvalues)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvz(ipoint)=jetvz(ipoint)+stoc*ww+prefactor*(f2vz(j)-f3vz(j))*zz+ &
     0.25d0*dt*(f2vz(j)+2.d0*f1vz(j)+f3vz(j))
! Evaporation rate at the predicted positions with the new velocity
! (device_maxwell_evap_stress_point).
    fev=0.d0
    if(evaporative .and. .not.(jetfr(ipoint) .or. ipoint>=npjet) .and. &
     .not.(ipoint==npjet-1 .and. .not.linserted))then
      if(ipoint==npjet-2 .and. .not.linserted)then
        dx=y1xx(npjet)-y1xx(ipoint); dy=y1yy(npjet)-y1yy(ipoint)
        dz=y1zz(npjet)-y1zz(ipoint)
      else
        dx=y1xx(ipoint+1)-y1xx(ipoint); dy=y1yy(ipoint+1)-y1yy(ipoint)
        dz=y1zz(ipoint+1)-y1zz(ipoint)
      endif
      beadlen=dsqrt(dx*dx+dy*dy+dz*dz)
      vnorm=dsqrt(jetvx(ipoint)*jetvx(ipoint)+jetvy(ipoint)*jetvy(ipoint)+ &
       jetvz(ipoint)*jetvz(ipoint))
      if(beadlen>0.d0 .and. evairv>0.d0 .and. y1ev(ipoint)>0.d0)then
        re=(2.d0*dsqrt(y1ev(ipoint)/(3.14159265358979323846d0*beadlen))* &
         vnorm)/evairv
        fev=-evmasscoeff*0.495d0*(re**(1.d0/3.d0))*sqrevsc* &
         (evcsvapour-evumidity)*3.14159265358979323846d0*beadlen
      endif
    endif
! Positions and evaporated volume.
    f2x=jetvx(ipoint); f2y=jetvy(ipoint); f2z=jetvz(ipoint)
    if(jetfr(ipoint) .or. (ipoint==npjet-1 .and. .not.linserted))then
      f2x=0.d0; f2y=0.d0; f2z=0.d0
    endif
    if(ipoint==npjet)then
      f2x=0.d0
      if(liniperturb)then
        y1y=jetyy(ipoint)+dt*f1yy(j)
        y1z=jetzz(ipoint)+dt*f1zz(j)
        f2y=-pfreq*y1z; f2z=pfreq*y1y
      else
        f2y=0.d0; f2z=0.d0
      endif
    endif
    jetxx(ipoint)=jetxx(ipoint)+0.5d0*dt*(f1xx(j)+f2x)
    jetyy(ipoint)=jetyy(ipoint)+0.5d0*dt*(f1yy(j)+f2y)
    jetzz(ipoint)=jetzz(ipoint)+0.5d0*dt*(f1zz(j)+f2z)
    if(evaporative)then
      ve=jetve(ipoint)+0.5d0*dt*(f1ev(j)+fev)
      if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
      jetve(ipoint)=ve
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_update

 subroutine accelerator_map_statistics(counterlpath,ncounterlpath, &
   maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: ncounterlpath
  double precision, intent(in) :: counterlpath,maxstress,maxstressposx
#ifdef _OPENACC
  if(.not.accelerator_statistics_mapped)then
!$acc enter data copyin(counterlpath,ncounterlpath,maxstress, &
!$acc& maxstressposx,statistics_step_max,statistics_step_index)
    accelerator_statistics_mapped=.true.
  endif
#endif
 end subroutine accelerator_map_statistics

 subroutine accelerator_set_persistent(enabled)
  implicit none
  logical, intent(in) :: enabled
  if(.not.enabled)call accelerator_set_async(.false.)
  if(.not.enabled)accelerator_topology_decided=.false.
  if(.not.enabled)accelerator_statistics_stored=.false.
  accelerator_persistent=enabled
  if(.not.enabled)accelerator_device_state_authoritative=.false.
  return
 end subroutine accelerator_set_persistent

 subroutine accelerator_set_topology_enabled(enabled)
  implicit none
  logical, intent(in) :: enabled
  accelerator_topology_enabled=enabled
  if(.not.enabled)accelerator_device_state_authoritative=.false.
 end subroutine accelerator_set_topology_enabled

 subroutine accelerator_mark_device_state(authoritative)
  implicit none
  logical, intent(in) :: authoritative
  accelerator_device_state_authoritative=authoritative
  if(authoritative)accelerator_last_host_sync_step=-huge(0)
 end subroutine accelerator_mark_device_state

 logical function accelerator_device_state_is_current()
  implicit none
  accelerator_device_state_is_current=accelerator_device_state_authoritative
 end function accelerator_device_state_is_current

 logical function accelerator_is_topology_enabled()
  implicit none
  accelerator_is_topology_enabled=accelerator_topology_enabled
 end function accelerator_is_topology_enabled

 logical function accelerator_is_persistent()
  implicit none
  accelerator_is_persistent=accelerator_persistent
 end function accelerator_is_persistent

 logical function accelerator_statistics_on_device()
! The step statistics live on the device once a device step has mapped
! them, also at a capacity change: there the jet arrays are released and
! rebound (persistent off, topology on) after the step has already added
! its path length on the device.  Until 2026-10-06 that step took the host
! branch of statistic_driver, so the device count missed it and the
! printed path length of that interval was too large by one step in
! 20,000 (RK device paths; the Platen end-of-step kernel counts it itself).
  implicit none
  accelerator_statistics_on_device=accelerator_statistics_mapped .and. &
   (accelerator_persistent .or. accelerator_topology_enabled)
 end function accelerator_statistics_on_device

 logical function accelerator_host_state_is_current(nstep)
  implicit none
  integer, intent(in) :: nstep
  accelerator_host_state_is_current=accelerator_persistent .and. &
   nstep==accelerator_last_host_sync_step
 end function accelerator_host_state_is_current

 subroutine accelerator_update_host_point(ipoint,jetxx,jetyy,jetzz, &
   jetst,jetvx,jetvy,jetvz)
  implicit none
  integer, intent(in) :: ipoint
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  if(.not.(accelerator_persistent .or. accelerator_topology_enabled))return
#ifdef _OPENACC
! The collector radius also needs the adjacent active bead.
  call accelerator_wait()
!$acc update self(jetxx(ipoint:ipoint+1),jetyy(ipoint:ipoint+1), &
!$acc& jetzz(ipoint:ipoint+1),jetst(ipoint),jetvx(ipoint), &
!$acc& jetvy(ipoint),jetvz(ipoint))
#endif
  return
 end subroutine accelerator_update_host_point

 subroutine accelerator_update_host_evaporation_point(ipoint,jetve)
  implicit none
  integer, intent(in) :: ipoint
  double precision, intent(inout) :: jetve(0:)
  if(.not.(accelerator_persistent .or. accelerator_topology_enabled))return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update self(jetve(ipoint))
#endif
 end subroutine accelerator_update_host_evaporation_point

 subroutine accelerator_update_host_state(npjet,jetxx,jetyy,jetzz, &
   jetst,jetvx,jetvy,jetvz,nstep)
  implicit none
  integer, intent(in) :: npjet
  integer, intent(in), optional :: nstep
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  if(.not.accelerator_persistent)return
  if(present(nstep))then
    if(nstep==accelerator_last_host_sync_step)return
  endif
#ifdef _OPENACC
  call accelerator_wait()
!$acc update self(jetxx(0:npjet),jetyy(0:npjet),jetzz(0:npjet), &
!$acc& jetst(0:npjet),jetvx(0:npjet),jetvy(0:npjet),jetvz(0:npjet))
#endif
  if(present(nstep))accelerator_last_host_sync_step=nstep
  return
 end subroutine accelerator_update_host_state

 subroutine accelerator_update_host_capacity_state(npjet,jetxx,jetyy,jetzz, &
   jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr,nstep)
  implicit none
  integer, intent(in) :: npjet
  integer, intent(in), optional :: nstep
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
  if(.not.(accelerator_persistent .or. accelerator_topology_enabled))return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update self(jetxx(0:npjet),jetyy(0:npjet),jetzz(0:npjet), &
!$acc& jetst(0:npjet),jetvx(0:npjet),jetvy(0:npjet),jetvz(0:npjet), &
!$acc& jetms(0:npjet),jetch(0:npjet),jetvl(0:npjet),jetfr(0:npjet))
#endif
  if(present(nstep))accelerator_last_host_sync_step=nstep
 end subroutine accelerator_update_host_capacity_state

 subroutine accelerator_update_host_evaporation_state(npjet,jetve,jetce)
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(inout) :: jetve(0:),jetce(0:)
  if(.not.(accelerator_persistent .or. accelerator_topology_enabled))return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update self(jetve(0:npjet),jetce(0:npjet))
#endif
 end subroutine accelerator_update_host_evaporation_state

 subroutine accelerator_topology_check(npjet,mxnpjet,inpjet,linserted, &
   lremove,h,resolution,dresolution,thresolution,ivelocity,istress,imassa, &
   icharge,ivolume,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch, &
   jetvl,jetfr,ladd,lresize,lrem)
! One timestep of device topology decisions: release of the blocked nozzle
! bead or insertion of a new one (as add_jetbead) and, with removal enabled,
! freezing at the collector and the removal test (as remove_jetbead).  The
! outcome is written to accelerator_topology_flags and returns to the host
! in a single transfer; until 2026-10-01 every decision scalar was copied
! separately (two uploads and five downloads per step).  Bead records cross
! only when an insertion actually happens, or in
! accelerator_finish_remove_bead when a removal does.
  implicit none
  integer, intent(inout) :: npjet
  integer, intent(in) :: mxnpjet,inpjet
  logical, intent(inout) :: linserted
  logical, intent(in) :: lremove
  logical, intent(out) :: ladd,lresize,lrem
  double precision, intent(in) :: h,resolution,dresolution,thresolution
  double precision, intent(in) :: ivelocity,istress,imassa,icharge,ivolume
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
  integer :: ipoint,lastpoint

! accelerator_platen_end_step may already have decided this step on
! the device; then only the readback below remains.
  if(accelerator_topology_decided)then
    accelerator_topology_decided=.false.
  else
#ifdef _OPENACC
!$acc serial async(accelerator_queue) present(jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms, &
!$acc& jetch,jetvl,jetfr,accelerator_topology_flags)
#endif
    call device_topology_decide(npjet,mxnpjet,inpjet,linserted,lremove,h, &
     resolution,dresolution,thresolution,ivelocity,istress,imassa,icharge, &
     ivolume,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
#ifdef _OPENACC
!$acc end serial
#endif
! Beads at the collector are frozen with or without removal
! (remove_jetbead).  The upper bound covers a bead added above; the device
! npjet decides.
    lastpoint=min(npjet+1,mxnpjet)
#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) present(jetxx,jetfr,accelerator_topology_flags)
#endif
    do ipoint=inpjet,lastpoint
      if(accelerator_topology_flags(4)==0 .and. &
       ipoint<=accelerator_topology_flags(1))then
        if(jetxx(ipoint)>=h)then
          jetfr(ipoint)=.true.
          jetxx(ipoint)=h
        endif
      endif
    enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
  endif
#ifdef _OPENACC
!$acc update self(accelerator_topology_flags) async(accelerator_queue)
  call accelerator_wait()
#endif
  npjet=accelerator_topology_flags(1)
  linserted=accelerator_topology_flags(2)==1
  ladd=accelerator_topology_flags(3)==1
  lresize=accelerator_topology_flags(4)==1
  lrem=accelerator_topology_flags(5)==1
#ifdef _OPENACC
  if(ladd)then
! The host statistics need only the injected mass and charge.  The complete
! bead state remains device-resident until an output/checkpoint or resize.
    if(accelerator_device_state_authoritative)then
!$acc update self(jetms(npjet-1),jetch(npjet-1))
    else
! The development host-force path continues on the host after the topology
! event, so it needs the complete new-bead state.
!$acc update self(jetxx(npjet-1:npjet),jetyy(npjet-1:npjet), &
!$acc& jetzz(npjet-1:npjet),jetst(npjet-1:npjet),jetvx(npjet-1:npjet), &
!$acc& jetvy(npjet-1:npjet),jetvz(npjet-1:npjet),jetms(npjet-1:npjet), &
!$acc& jetch(npjet-1:npjet),jetvl(npjet-1:npjet),jetfr(npjet-1:npjet))
    endif
  endif
#endif
 end subroutine accelerator_topology_check

 subroutine device_topology_decide(npjet,mxnpjet,inpjet,linserted,lremove,h, &
   resolution,dresolution,thresolution,ivelocity,istress,imassa,icharge, &
   ivolume,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
! The serial part of one step's topology decisions (release or insertion at
! the nozzle, capacity exhaustion, removal test) on the device, written to
! accelerator_topology_flags; the freezing loop follows in the caller.
#ifdef _OPENACC
!$acc routine seq
#endif
  implicit none
  integer, intent(in) :: npjet,mxnpjet,inpjet
  logical, intent(in) :: linserted,lremove
  double precision, intent(in) :: h,resolution,dresolution,thresolution
  double precision, intent(in) :: ivelocity,istress,imassa,icharge,ivolume
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
  integer :: np,newadd,newresize,newremove
  logical :: ins
  double precision :: dx,dy,dz,distance,scale
  np=npjet
  ins=linserted
  newadd=0
  newresize=0
  newremove=0
  if(.not.ins)then
    dx=jetxx(np-2)-jetxx(np)
    dy=jetyy(np-2)-jetyy(np)
    dz=jetzz(np-2)-jetzz(np)
    distance=dsqrt(dx*dx+dy*dy+dz*dz)
    if(distance>=dresolution)then
      ins=.true.
      jetst(np-1)=0.d0
      jetvx(np-1)=ivelocity
      jetvy(np-1)=0.d0
      jetvz(np-1)=0.d0
    endif
  else
    dx=jetxx(np-1)-jetxx(np)
    dy=jetyy(np-1)-jetyy(np)
    dz=jetzz(np-1)-jetzz(np)
    distance=dsqrt(dx*dx+dy*dy+dz*dz)
    if(distance>=thresolution .and. np>=mxnpjet)then
      newresize=1
    elseif(distance>=thresolution)then
      np=np+1
      jetfr(np)=jetfr(np-1)
      jetxx(np)=jetxx(np-1)
      jetyy(np)=jetyy(np-1)
      jetzz(np)=jetzz(np-1)
      jetst(np)=jetst(np-1)
      jetvx(np)=jetvx(np-1)
      jetvy(np)=jetvy(np-1)
      jetvz(np)=jetvz(np-1)
      jetms(np)=jetms(np-1)
      jetch(np)=jetch(np-1)
      jetvl(np)=jetvl(np-1)
      jetfr(np-1)=.false.
      jetst(np-1)=istress
      jetvx(np-1)=ivelocity
      jetvy(np-1)=0.d0
      jetvz(np-1)=0.d0
      jetms(np-1)=imassa*ivolume
      jetch(np-1)=icharge*ivolume
      jetvl(np-1)=ivolume
      dx=jetxx(np-2)-jetxx(np)
      dy=jetyy(np-2)-jetyy(np)
      dz=jetzz(np-2)-jetzz(np)
      distance=dsqrt(dx*dx+dy*dy+dz*dz)
      scale=resolution/distance
      jetxx(np-1)=jetxx(np)+scale*dx
      jetyy(np-1)=jetyy(np)+scale*dy
      jetzz(np-1)=jetzz(np)+scale*dz
      newadd=1
      ins=.false.
    endif
  endif
! A step that must first grow the capacity stops here: the host grows the
! arrays and decides the step again (main), insertion and removal included.
! Clamping a bead at the collector keeps it at x>=h, so the test can precede
! the freezing loop below.
  if(lremove .and. newresize==0)then
    if(jetxx(inpjet)>=h .and. jetxx(inpjet+1)>=h)newremove=1
  endif
  accelerator_topology_flags(1)=np
  accelerator_topology_flags(2)=merge(1,0,ins)
  accelerator_topology_flags(3)=newadd
  accelerator_topology_flags(4)=newresize
  accelerator_topology_flags(5)=newremove
 end subroutine device_topology_decide

 subroutine accelerator_update_device_added_evaporation(npjet,ivolume,jetve,jetce)
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(in) :: ivolume
  double precision, intent(inout) :: jetve(0:),jetce(0:)
#ifdef _OPENACC
  call accelerator_wait()
!$acc serial present(jetve,jetce)
#endif
  jetve(npjet)=jetve(npjet-1)
  jetce(npjet)=jetce(npjet-1)
  jetve(npjet-1)=ivolume
  jetce(npjet-1)=0.d0
#ifdef _OPENACC
!$acc end serial
  if(.not.accelerator_device_state_authoritative)then
!$acc update self(jetve(npjet-1:npjet),jetce(npjet-1:npjet))
  endif
#endif
 end subroutine accelerator_update_device_added_evaporation

 subroutine accelerator_finish_remove_bead(inpjet,lrem,jetxx,jetyy,jetzz, &
   jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
! Completes a removal decided by accelerator_topology_check: advances the
! active lower bound and returns the records the removal output needs.
  implicit none
  integer, intent(inout) :: inpjet
  logical, intent(in) :: lrem
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
  if(.not.lrem)return
  inpjet=inpjet+1
#ifdef _OPENACC
! Removal observables need both the removed bead and its active neighbour.
  call accelerator_wait()
!$acc update self(jetxx(inpjet-1:inpjet),jetyy(inpjet-1:inpjet), &
!$acc& jetzz(inpjet-1:inpjet), &
!$acc& jetst(inpjet-1),jetvx(inpjet-1),jetvy(inpjet-1),jetvz(inpjet-1), &
!$acc& jetms(inpjet-1),jetch(inpjet-1),jetvl(inpjet-1),jetfr(inpjet-1))
#endif
 end subroutine accelerator_finish_remove_bead

 subroutine accelerator_update_host_removed_evaporation(ipoint,jetve,jetce)
  implicit none
  integer, intent(in) :: ipoint
  double precision, intent(inout) :: jetve(0:),jetce(0:)
#ifdef _OPENACC
  call accelerator_wait()
!$acc update self(jetve(ipoint),jetce(ipoint))
#endif
 end subroutine accelerator_update_host_removed_evaporation

 subroutine accelerator_update_device_removed(first,last,jetst,jetvx,jetvy, &
   jetvz,jetms,jetch,jetvl)
  implicit none
  integer, intent(in) :: first,last
  double precision, intent(in) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: jetms(0:),jetch(0:),jetvl(0:)
  if(.not.(accelerator_persistent .or. accelerator_topology_enabled) .or. last<first)return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update device(jetst(first:last),jetvx(first:last),jetvy(first:last), &
!$acc& jetvz(first:last),jetms(first:last),jetch(first:last),jetvl(first:last))
#endif
 end subroutine accelerator_update_device_removed

 subroutine accelerator_update_device_removed_evaporation(first,last,jetve,jetce)
  implicit none
  integer, intent(in) :: first,last
  double precision, intent(in) :: jetve(0:),jetce(0:)
  if(.not.(accelerator_persistent .or. accelerator_topology_enabled) .or. last<first)return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update device(jetve(first:last),jetce(first:last))
#endif
 end subroutine accelerator_update_device_removed_evaporation

 subroutine accelerator_store_statistics(inpjet,npjet,jetxx,jetyy, &
   jetzz,jetst,counterlpath,ncounterlpath,maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: inpjet,npjet
  integer, intent(inout) :: ncounterlpath
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint

  if(accelerator_statistics_stored)then
! Already done on the device by accelerator_platen_end_step.
    accelerator_statistics_stored=.false.
    return
  endif
  if(.not.accelerator_statistics_on_device())return
#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) gang vector present(jetst,statistics_step_max, &
!$acc& statistics_step_index) reduction(max:statistics_step_index)
#endif
  do ipoint=inpjet,npjet
    if(jetst(ipoint)==statistics_step_max)then
      statistics_step_index=max(statistics_step_index,ipoint)
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
!$acc serial async(accelerator_queue) present(jetxx,counterlpath,ncounterlpath,maxstress, &
!$acc& maxstressposx,statistics_step_max,statistics_step_index)
#endif
  ncounterlpath=ncounterlpath+1
  if(statistics_step_max>=maxstress)then
    maxstress=statistics_step_max
    maxstressposx=jetxx(statistics_step_index)
  endif
  statistics_step_max=-huge(0.d0)
  statistics_step_index=-1
#ifdef _OPENACC
!$acc end serial
#endif
  return
 end subroutine accelerator_store_statistics

 subroutine accelerator_update_host_statistics(counterlpath, &
   ncounterlpath,maxstress,maxstressposx)
  implicit none
  integer, intent(inout) :: ncounterlpath
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  if(.not.accelerator_statistics_mapped)return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update self(counterlpath,ncounterlpath,maxstress,maxstressposx)
#endif
  return
 end subroutine accelerator_update_host_statistics

 subroutine accelerator_update_device_statistics(counterlpath, &
   ncounterlpath,maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: ncounterlpath
  double precision, intent(in) :: counterlpath,maxstress,maxstressposx
  if(.not.accelerator_statistics_mapped)return
#ifdef _OPENACC
  call accelerator_wait()
!$acc update device(counterlpath,ncounterlpath,maxstress,maxstressposx)
#endif
  return
 end subroutine accelerator_update_device_statistics

 subroutine accelerator_prepare()
  implicit none
#ifdef _OPENACC
!$acc init
#endif
  return
 end subroutine accelerator_prepare

 logical function accelerator_eom3_stage(firstpoint,lastpoint,npjet, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf,jetms,jetch,jetfr, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,linserted,liniperturb,lairdrag, &
   lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
   att,fve,gr,ks,li,vfield,velext,stochastic_model,noisefric,yve, &
   apply_airdrag,collector_curvature,fev_evap,ev_airv,ev_masscoeff, &
   ev_sqrevsc,ev_csvapour,ev_umidity,ev_cp0,ev_bev,ev_mev,ev_tev)
! With fev_evap (and the ev_ parameters of the Yarin law) the loop also
! computes the evaporation rate and replaces the Newtonian stress
! derivative with the Maxwell one (device_maxwell_evap_stress_point), the
! work of a second kernel, accelerator_maxwell_evap_stress_3d, until
! 2026-10-05.
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,nfieldtype
  logical, intent(in) :: linserted,liniperturb,lairdrag,lflorentz,luppot
  logical, intent(in) :: stochastic_model
  double precision, intent(in) :: pfreq,consistency,findex,yieldstress
  double precision, intent(in) :: att,fve,gr,ks,li,vfield,velext,noisefric
  double precision, intent(in), optional :: yve(0:)
  logical, intent(in), optional :: apply_airdrag,collector_curvature
  double precision, intent(inout), optional :: fev_evap(0:)
  double precision, intent(in), optional :: ev_airv,ev_masscoeff,ev_sqrevsc
  double precision, intent(in), optional :: ev_csvapour,ev_umidity,ev_cp0
  double precision, intent(in), optional :: ev_bev,ev_mev,ev_tev
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yvl(0:)
  double precision, intent(in) :: ycf(0:,1:),jetms(0:),jetch(0:)
  logical, intent(in) :: jetfr(0:)
  double precision, intent(out) :: fxx(0:),fyy(0:),fzz(0:),fst(0:)
  double precision, intent(out) :: fvx(0:),fvy(0:),fvz(0:)

  integer :: ipoint,j,nout
  double precision :: dxu,dyu,dzu,dxd,dyd,dzd,lup,ldown
  double precision :: tux,tuy,tuz,tdx,tdy,tdz,beadvel
  double precision :: v1x,v1y,v1z,v2x,v2y,v2z,l1,l2,dotp
  double precision :: nbx,nby,nbz,lnb,b,c,t,scale1,scale2
  double precision :: ccx,ccy,ccz,rcx,rcy,rcz,radius,curvature
  double precision :: factor1,factor2,factor3,factor4,factor5
  double precision :: fvet,kst,attt,lit,veltangent,cmass,fvolume,fvolume_prev
  logical :: straight,use_evap,use_airdrag,use_collector_curvature
  logical :: fuse_evap
  double precision :: e_airv,e_masscoeff,e_sqrevsc,e_csvapour,e_umidity
  double precision :: e_cp0,e_bev,e_mev,e_tev,fst_evap

  accelerator_eom3_stage=.false.
  if(.not.accelerator_eom_env_checked)then
    block
      character(len=16) :: env
      env=''
      call get_environment_variable('JETSPIN_OPENACC_DISABLE_EOM',env)
      accelerator_eom_disabled=trim(env)=='1'
    end block
    accelerator_eom_env_checked=.true.
  endif
  if(accelerator_eom_disabled)return
! Without air drag (lairdrag false) the loop leaves out the drag and lift
! terms (use_airdrag), as eom3 does; until 2026-10-06 it returned .false.
  if(lflorentz .or. luppot)return
  if(nfieldtype/=0 .or. lastpoint/=npjet)return
  nout=lastpoint-firstpoint
  ! Evaluate OPTIONAL presence on the host and pass a plain scalar into the
  ! device kernel; PRESENT() itself is not reliable inside OpenACC regions.
  use_evap=present(yve)
  use_airdrag=lairdrag
  if(present(apply_airdrag))use_airdrag=apply_airdrag
  use_collector_curvature=.false.
  if(present(collector_curvature))use_collector_curvature=collector_curvature
  fuse_evap=present(fev_evap) .and. present(yve)
  e_airv=0.d0; e_masscoeff=0.d0; e_sqrevsc=0.d0; e_csvapour=0.d0
  e_umidity=0.d0; e_cp0=0.d0; e_bev=0.d0; e_mev=0.d0; e_tev=0.d0
  if(fuse_evap)then
    e_airv=ev_airv; e_masscoeff=ev_masscoeff; e_sqrevsc=ev_sqrevsc
    e_csvapour=ev_csvapour; e_umidity=ev_umidity; e_cp0=ev_cp0
    e_bev=ev_bev; e_mev=ev_mev; e_tev=ev_tev
  endif

#ifdef _OPENACC
!$acc parallel loop async(accelerator_queue) gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yst(0:npjet),yvx(0:npjet),yvy(0:npjet), &
!$acc& yvz(0:npjet)) present(yvl,jetms,jetch,jetfr) &
#if defined(JETSPIN_DEV_HOST_COULOMB_ORACLE) || defined(JETSPIN_DEV_HOST_FORCE_ORACLE)
!$acc& present(ycf) &
#else
!$acc& present_or_copyin(ycf) &
#endif
!$acc& present(yve) present(fev_evap) &
!$acc& present_or_copyout(fxx,fyy,fzz,fst,fvx,fvy,fvz) &
!$acc& private(j,dxu,dyu,dzu,dxd,dyd,dzd,lup,ldown,tux,tuy,tuz, &
!$acc& tdx,tdy,tdz,beadvel,v1x,v1y,v1z,v2x,v2y,v2z,l1,l2,dotp, &
!$acc& nbx,nby,nbz,lnb,b,c,t,scale1,scale2,ccx,ccy,ccz,rcx,rcy, &
!$acc& rcz,radius,curvature,factor1,factor2,factor3,factor4,factor5, &
!$acc& fvet,kst,attt,lit,veltangent,cmass,fvolume,fvolume_prev,straight, &
!$acc& fst_evap)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    fxx(j)=0.d0
    fyy(j)=0.d0
    fzz(j)=0.d0
    fst(j)=0.d0
    fvx(j)=0.d0
    fvy(j)=0.d0
    fvz(j)=0.d0
! Every bead skipped below before its stress is set (frozen, blocked nozzle
! bead, nozzle) also has a zero Maxwell stress derivative.
    if(fuse_evap)call device_maxwell_evap_stress_point(ipoint,j,npjet, &
     linserted,jetfr,fev_evap,fst_evap,yxx,yyy,yzz,yvx,yvy,yvz,yst,yvl,yve, &
     e_airv,e_masscoeff,e_sqrevsc,e_csvapour,e_umidity,e_cp0,e_bev,e_mev, &
     e_tev,consistency,findex,yieldstress)

    if(jetfr(ipoint))cycle
    if(ipoint==npjet-1 .and. .not.linserted)cycle
    if(ipoint==npjet)then
      if(liniperturb)then
        fyy(j)=-pfreq*yzz(ipoint)
        fzz(j)= pfreq*yyy(ipoint)
        fvy(j)=-(pfreq**2.d0)*yyy(ipoint)
        fvz(j)=-(pfreq**2.d0)*yzz(ipoint)
      endif
      cycle
    endif

    if(ipoint==npjet-2 .and. .not.linserted)then
      dxu=yxx(ipoint)-yxx(npjet)
      dyu=yyy(ipoint)-yyy(npjet)
      dzu=yzz(ipoint)-yzz(npjet)
    else
      dxu=yxx(ipoint)-yxx(ipoint+1)
      dyu=yyy(ipoint)-yyy(ipoint+1)
      dzu=yzz(ipoint)-yzz(ipoint+1)
    endif
    lup=dsqrt(dxu*dxu+dyu*dyu+dzu*dzu)
    tux=dxu/lup
    tuy=dyu/lup
    tuz=dzu/lup
    if(ipoint==npjet-2 .and. .not.linserted)then
      beadvel=(yvx(ipoint)-yvx(npjet))*tux+ &
       (yvy(ipoint)-yvy(npjet))*tuy+ &
       (yvz(ipoint)-yvz(npjet))*tuz
    else
      beadvel=(yvx(ipoint)-yvx(ipoint+1))*tux+ &
       (yvy(ipoint)-yvy(ipoint+1))*tuy+ &
       (yvz(ipoint)-yvz(ipoint+1))*tuz
    endif
    cmass=1.d0
    fvolume=yvl(ipoint)
    fvolume_prev=fvolume
    if(ipoint>firstpoint .or. &
     (ipoint==firstpoint .and. firstpoint>0 .and. use_collector_curvature)) &
     fvolume_prev=yvl(ipoint-1)
    if(use_evap)then
      cmass=yve(ipoint)/yvl(ipoint)
      fvolume=yve(ipoint)
      if(ipoint>firstpoint .or. &
       (ipoint==firstpoint .and. firstpoint>0 .and. use_collector_curvature)) &
       fvolume_prev=yve(ipoint-1)
    endif
    fvet=fve/(jetms(ipoint)*cmass)
    if(stochastic_model .and. yst(ipoint)<=0.d0)then
      factor1=0.d0
    else
      factor1=fvet*fvolume*(yst(ipoint)/lup)
    endif

    fxx(j)=yvx(ipoint)
    fyy(j)=yvy(ipoint)
    fzz(j)=yvz(ipoint)
    if(fuse_evap)then
      fst(j)=fst_evap
    else
      fst(j)=yieldstress+consistency*(beadvel/lup)**findex-yst(ipoint)
    endif
    ! The electric acceleration is divided by the actual (evaporated) bead
    ! mass in the CPU Maxwell evaporation equation.  Keep the same scaling
    ! on device; omitting cmass creates an O(1/cmass) velocity derivative
    ! error as soon as evaporation is enabled.
    fvx(j)=gr+(jetch(ipoint)/(jetms(ipoint)*cmass))*vfield- &
     factor1*tux+ycf(ipoint,1)
    fvy(j)=-factor1*tuy+ycf(ipoint,2)
    fvz(j)=-factor1*tuz+ycf(ipoint,3)

    veltangent=0.d0
    if(use_airdrag)then
      veltangent=(yvx(ipoint)-velext)*tux+yvy(ipoint)*tuy+ &
       yvz(ipoint)*tuz
      attt=att/(jetms(ipoint)*cmass)
      factor4=attt*(dabs(lup)**0.905d0)*(dabs(veltangent)**1.19d0)
      fvx(j)=fvx(j)-factor4*tux
      fvy(j)=fvy(j)-factor4*tuy
      fvz(j)=fvz(j)-factor4*tuz
    endif
    if(stochastic_model)then
      fvx(j)=fvx(j)-noisefric*yvx(ipoint)
      fvy(j)=fvy(j)-noisefric*yvy(ipoint)
      fvz(j)=fvz(j)-noisefric*yvz(ipoint)
    endif

    if(ipoint==firstpoint .and. &
     (firstpoint==0 .or. .not.use_collector_curvature))cycle

    dxd=yxx(ipoint-1)-yxx(ipoint)
    dyd=yyy(ipoint-1)-yyy(ipoint)
    dzd=yzz(ipoint-1)-yzz(ipoint)
    ldown=dsqrt(dxd*dxd+dyd*dyd+dzd*dzd)
    tdx=dxd/ldown
    tdy=dyd/ldown
    tdz=dzd/ldown
    factor2=0.d0
    if(ipoint>firstpoint)then
      if(stochastic_model .and. yst(ipoint-1)<=0.d0)then
        factor2=0.d0
      else
        if(use_evap)then
          factor2=fvet*yve(ipoint-1)*(yst(ipoint-1)/ldown)
        else
          factor2=fvet*yvl(ipoint-1)*(yst(ipoint-1)/ldown)
        endif
      endif
    endif

! Local three-point curvature calculation. Each iteration reads only the
! current bead and its two neighbours.
    v1x=yxx(ipoint+1)-yxx(ipoint)
    v1y=yyy(ipoint+1)-yyy(ipoint)
    v1z=yzz(ipoint+1)-yzz(ipoint)
    v2x=yxx(ipoint-1)-yxx(ipoint)
    v2y=yyy(ipoint-1)-yyy(ipoint)
    v2z=yzz(ipoint-1)-yzz(ipoint)
    l1=dsqrt(v1x*v1x+v1y*v1y+v1z*v1z)
    l2=dsqrt(v2x*v2x+v2y*v2y+v2z*v2z)
    straight=(l1==0.d0 .or. l2==0.d0)
    curvature=0.d0
    rcx=0.d0
    rcy=0.d0
    rcz=0.d0
    if(.not.straight)then
      v1x=v1x/l1
      v1y=v1y/l1
      v1z=v1z/l1
      v2x=v2x/l2
      v2y=v2y/l2
      v2z=v2z/l2
      dotp=v2x*v1x+v2y*v1y+v2z*v1z
      nbx=v2x-dotp*v1x
      nby=v2y-dotp*v1y
      nbz=v2z-dotp*v1z
      lnb=dsqrt(nbx*nbx+nby*nby+nbz*nbz)
      straight=(lnb==0.d0)
      if(.not.straight)then
        nbx=nbx/lnb
        nby=nby/lnb
        nbz=nbz/lnb
        b=(yxx(ipoint-1)-yxx(ipoint))*v1x+ &
         (yyy(ipoint-1)-yyy(ipoint))*v1y+ &
         (yzz(ipoint-1)-yzz(ipoint))*v1z
        c=(yxx(ipoint-1)-yxx(ipoint))*nbx+ &
         (yyy(ipoint-1)-yyy(ipoint))*nby+ &
         (yzz(ipoint-1)-yzz(ipoint))*nbz
        straight=(c==0.d0)
        if(.not.straight)then
          t=0.5d0*(l1-b)/c
          scale1=b/2.d0+c*t
          scale2=c/2.d0-b*t
          rcx=scale1*v1x+scale2*nbx
          rcy=scale1*v1y+scale2*nby
          rcz=scale1*v1z+scale2*nbz
          radius=dsqrt(rcx*rcx+rcy*rcy+rcz*rcz)
          curvature=1.d0/radius
          rcx=rcx/radius
          rcy=rcy/radius
          rcz=rcz/radius
        endif
      endif
    endif

    kst=ks/(jetms(ipoint)*cmass)
    factor3=0.25d0*((dsqrt(fvolume)/dsqrt(lup))+ &
     (dsqrt(fvolume_prev)/dsqrt(ldown)))**2.d0
    fvx(j)=fvx(j)+factor2*tdx+kst*curvature*factor3*rcx
    fvy(j)=fvy(j)+factor2*tdy+kst*curvature*factor3*rcy
    fvz(j)=fvz(j)+factor2*tdz+kst*curvature*factor3*rcz
    if(use_airdrag)then
      lit=li/(jetms(ipoint)*cmass)
      factor5=factor3*lup*curvature*(veltangent**2.d0)
      fvx(j)=fvx(j)-lit*factor5*rcx
      fvy(j)=fvy(j)-lit*factor5*rcy
      fvz(j)=fvz(j)-lit*factor5*rcz
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif

  accelerator_eom3_stage=.true.
 end function accelerator_eom3_stage

end module accelerator_mod
