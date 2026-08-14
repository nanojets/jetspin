module accelerator_mod

 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.true.
 integer, parameter, public :: rheology_maxwell=1
 integer, parameter, public :: rheology_kelvin_voigt=2
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
#ifdef _OPENACC
!$acc declare create(accelerator_smoothed_charge)
#endif

 public :: accelerator_prepare
 public :: accelerator_eom3_stage
 public :: accelerator_maxwell_evap_stage
 public :: accelerator_kv_evap_stage
 public :: accelerator_maxwell_evap_force_correction
 public :: accelerator_maxwell_rk4_stage_update
 public :: accelerator_maxwell_rk4_final_update
 public :: accelerator_maxwell_commit_state
 public :: accelerator_evap_rk4_stage_update
 public :: accelerator_evap_rk2_final_update
 public :: accelerator_evap_rk4_final_update
 public :: accelerator_evap_commit_state
 public :: accelerator_compute_posnoinserted_3d
 public :: accelerator_smooth_charge_3d
 public :: accelerator_restore_charge
 public :: accelerator_coulomb_evap_3d
 public :: accelerator_coulomb_evap_compare
 public :: accelerator_evaporation_force_3d
 public :: accelerator_maxwell_stress_3d
 public :: accelerator_maxwell_evap_stress_3d
 public :: accelerator_kv_stress_3d
 public :: accelerator_kv_evap_stress_3d
 public :: accelerator_evaporation_geometry_3d
 public :: accelerator_set_persistent
 public :: accelerator_set_topology_enabled
 public :: accelerator_is_topology_enabled
 public :: accelerator_rebind_topology
 public :: accelerator_rebind_evaporation
 public :: accelerator_update_device_topology_state
 public :: accelerator_update_device_evaporation_state
 public :: accelerator_is_persistent
 public :: accelerator_update_host_state
 public :: accelerator_update_host_capacity_state
 public :: accelerator_update_host_evaporation_state
 public :: accelerator_release_jet_capacity
 public :: accelerator_release_evaporation_capacity
 public :: accelerator_update_host_point
 public :: accelerator_update_host_evaporation_point
 public :: accelerator_remove_bead
 public :: accelerator_update_host_removed_evaporation
 public :: accelerator_update_device_removed
 public :: accelerator_update_device_removed_evaporation
 public :: accelerator_add_bead
 public :: accelerator_update_device_added_evaporation
 public :: accelerator_host_state_is_current
 public :: accelerator_mark_device_state
 public :: accelerator_device_state_is_current
 public :: accelerator_store_statistics
 public :: accelerator_update_host_statistics
 public :: accelerator_update_device_statistics
 public :: accelerator_rk4_final_statistics
 public :: accelerator_euler_final_statistics
 public :: accelerator_rk2_final_statistics
 public :: accelerator_platen_predict
 public :: accelerator_platen_evap_predict
 public :: accelerator_platen_velocity
 public :: accelerator_platen_evap_velocity
 public :: accelerator_platen_positions
 public :: accelerator_platen_evap_positions
 public :: accelerator_platen_stress_statistics

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
  ! the evaporation derivative.  This is intentionally restricted to the
  ! compatible Test-16 physics until the complete eom4_ev force model is
  ! device-resident.
  ok=accelerator_eom3_stage(firstpoint,lastpoint,npjet,yxx,yyy,yzz,yst, &
   yvx,yvy,yvz,yvl,ycf,jetms,jetch,jetfr,fxx,fyy,fzz,fst,fvx,fvy,fvz, &
   linserted,liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
   consistency,findex,yieldstress,att,fveparam,gr,ks,li,vfield,velext, &
   stochastic_model,noisefric,yve)
  if(.not.ok)return
  call accelerator_maxwell_evap_stress_3d(firstpoint,lastpoint,npjet, &
   linserting=linserting,linserted=linserted,jetfr=jetfr,fev=fve,fst=fst, &
   yxx=yxx,yyy=yyy,yzz=yzz,yvx=yvx,yvy=yvy,yvz=yvz,yst=yst,yvl=yvl,yve=yve, &
   evairv=evairv,evmasscoeff=evmasscoeff,sqrevsc=sqrevsc, &
   evcsvapour=evcsvapour,evumidity=evumidity,cp0=cp0,Bev=Bev,mev=mev, &
   tev=tev,consistency=consistency,findex=findex,yieldstress=yieldstress)
 end subroutine accelerator_maxwell_evap_stage

 subroutine accelerator_kv_evap_stage(firstpoint,lastpoint,npjet, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf,jetms,jetch,jetfr, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,linserting,linserted,liniperturb,lairdrag, &
   lflorentz,luppot,nfieldtype,pfreq,consistency,findex,yieldstress, &
   att,fveparam,gr,ks,li,vfield,velext,evairv,evmasscoeff,sqrevsc, &
   evcsvapour,evumidity,cp0,Bev,mev,tev,evlim)
  implicit none
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
  logical :: ok

  ! Reproduce the established CPU eom3_KV_pos_v_ev semantics.  That routine
  ! does not add air drag or lift, even when airdrag is enabled in the input.
  ok=accelerator_eom3_stage(firstpoint,lastpoint,npjet,yxx,yyy,yzz,yst, &
   yvx,yvy,yvz,yvl,ycf,jetms,jetch,jetfr,fxx,fyy,fzz,fst,fvx,fvy,fvz, &
   linserted,liniperturb,lairdrag,lflorentz,luppot,nfieldtype,pfreq, &
   consistency,findex,yieldstress,att,fveparam,gr,ks,li,vfield,velext, &
   .false.,0.d0,yve,apply_airdrag=.false.,collector_curvature=.true.)
  if(.not.ok)return

  call accelerator_kv_evap_stress_3d(firstpoint,lastpoint,npjet, &
   linserting,linserted,jetfr,fve,fst,yxx,yyy,yzz,yvx,yvy,yvz, &
   fvx,fvy,fvz,yst,yvl,yve,evairv,evmasscoeff,sqrevsc,evcsvapour, &
   evumidity,cp0,Bev,mev,tev,evlim,acceleration_chunked=.true.)
 end subroutine accelerator_kv_evap_stage

 subroutine accelerator_maxwell_rk4_stage_update(firstpoint,lastpoint,h,stage, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
  implicit none
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
    ve=jetve(ipoint)+scale*fev(j)
    if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
    yev(ipoint)=ve
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_rk4_stage_update

 subroutine accelerator_evap_rk2_final_update(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve,jetvl, &
   f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
   f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f2ev, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
  implicit none
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
    ve=jetve(ipoint)+scale*(f1ev(j)+f2ev(j))
    if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
    yev(ipoint)=ve
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
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yev,evlim)
  implicit none
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
    ve=jetve(ipoint)+scale*(f1ev(j)+2.d0*(f2ev(j)+f3ev(j))+f4ev(j))
    if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
    yev(ipoint)=ve
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_rk4_final_update

 subroutine accelerator_maxwell_commit_state(firstpoint,lastpoint, &
   yxx,yyy,yzz,yst,yvx,yvy,yvz,yev, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetve, &
   counterlpath,ncounterlpath,maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yst(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yev(0:)
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:),jetve(0:)
  integer, intent(inout) :: ncounterlpath
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint
  double precision :: dx,dy,dz
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
    jetvz(ipoint)=yvz(ipoint); jetve(ipoint)=yev(ipoint)
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
  double precision :: dx,dy,dz,distance,scale
  if(linserted)return
#ifdef _OPENACC
!$acc serial present(yxx,yyy,yzz) private(dx,dy,dz,distance,scale)
#endif
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
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_compute_posnoinserted_3d

 subroutine accelerator_smooth_charge_3d(npjet,linserted,thresolution,dresolution, &
   yxx,yyy,yzz,jetch)
  implicit none
  integer, intent(in) :: npjet
  logical, intent(in) :: linserted
  double precision, intent(in) :: thresolution,dresolution
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(inout) :: jetch(0:)
  double precision :: dx,dy,dz,distance,factor
  if(linserted)return
#ifdef _OPENACC
!$acc serial present(yxx,yyy,yzz,jetch) private(dx,dy,dz,distance,factor)
#endif
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
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_smooth_charge_3d

 subroutine accelerator_restore_charge(npjet,linserted,jetch)
  implicit none
  integer, intent(in) :: npjet
  logical, intent(in) :: linserted
  double precision, intent(inout) :: jetch(0:)
  if(linserted)return
#ifdef _OPENACC
!$acc serial present(jetch)
#endif
  jetch(npjet-1)=accelerator_smoothed_charge
#ifdef _OPENACC
!$acc end serial
#endif
 end subroutine accelerator_restore_charge

 subroutine accelerator_evaporation_geometry_3d(firstpoint,lastpoint,npjet, &
   linserted,jetfr,yxx,yyy,yzz,yvx,yvy,yvz,yve,evairv,beadlen,reynolds)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: linserted,jetfr(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yve(0:),evairv
  double precision, intent(out) :: beadlen(0:),reynolds(0:)
  integer :: ipoint,j
  double precision :: dx,dy,dz,vnorm
#ifdef _OPENACC
!$acc parallel loop gang vector copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yvx(0:npjet),yvy(0:npjet),yvz(0:npjet), &
!$acc& yve(0:npjet),jetfr(0:npjet)) copyout(beadlen(0:lastpoint-firstpoint), &
!$acc& reynolds(0:lastpoint-firstpoint)) private(j,dx,dy,dz,vnorm)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    beadlen(j)=0.d0
    reynolds(j)=0.d0
    if(jetfr(ipoint) .or. ipoint>=npjet)cycle
    if(ipoint==npjet-1 .and. .not.linserted)cycle
    if(ipoint==npjet-2 .and. .not.linserted)then
      dx=yxx(ipoint)-yxx(npjet)
      dy=yyy(ipoint)-yyy(npjet)
      dz=yzz(ipoint)-yzz(npjet)
    else
      dx=yxx(ipoint)-yxx(ipoint+1)
      dy=yyy(ipoint)-yyy(ipoint+1)
      dz=yzz(ipoint)-yzz(ipoint+1)
    endif
    beadlen(j)=dsqrt(dx*dx+dy*dy+dz*dz)
    vnorm=dsqrt(yvx(ipoint)*yvx(ipoint)+yvy(ipoint)*yvy(ipoint)+ &
     yvz(ipoint)*yvz(ipoint))
    if(beadlen(j)>0.d0 .and. evairv>0.d0 .and. yve(ipoint)>0.d0)then
      reynolds(j)=(2.d0*dsqrt(yve(ipoint)/ &
       (3.14159265358979323846d0*beadlen(j)))*vnorm)/evairv
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_evaporation_geometry_3d

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

 subroutine accelerator_maxwell_stress_3d(firstpoint,lastpoint,npjet,linserted, &
   jetfr,fst,yxx,yyy,yzz,yvx,yvy,yvz,yst,yvl,yve,cp0,Bev,mev,tev, &
   consistency,findex,yieldstress)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: linserted,jetfr(0:)
  double precision, intent(out) :: fst(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:),yst(0:)
  double precision, intent(in) :: yvl(0:),yve(0:)
  double precision, intent(in) :: cp0,Bev,mev,tev,consistency,findex,yieldstress
  integer :: ipoint,j
  double precision :: dx,dy,dz,beadlen,beadvel,cp,ratmu,rattao
#ifdef _OPENACC
!$acc parallel loop gang vector copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yst(0:npjet), &
!$acc& yvl(0:npjet),yve(0:npjet),jetfr(0:npjet)) &
!$acc& copyout(fst(0:lastpoint-firstpoint)) private(j,dx,dy,dz,beadlen, &
!$acc& beadvel,cp,ratmu,rattao)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    fst(j)=0.d0
    if(jetfr(ipoint) .or. ipoint>=npjet)cycle
    if(ipoint==npjet-1 .and. .not.linserted)cycle
    if(ipoint==npjet-2 .and. .not.linserted)then
      dx=yxx(ipoint)-yxx(npjet)
      dy=yyy(ipoint)-yyy(npjet)
      dz=yzz(ipoint)-yzz(npjet)
    else
      dx=yxx(ipoint)-yxx(ipoint+1)
      dy=yyy(ipoint)-yyy(ipoint+1)
      dz=yzz(ipoint)-yzz(ipoint+1)
    endif
    beadlen=dsqrt(dx*dx+dy*dy+dz*dz)
    if(beadlen<=0.d0 .or. yve(ipoint)<=0.d0 .or. yvl(ipoint)<=0.d0)cycle
    if(ipoint==npjet-2 .and. .not.linserted)then
      beadvel=((yvx(ipoint)-yvx(npjet))*dx+ &
       (yvy(ipoint)-yvy(npjet))*dy+(yvz(ipoint)-yvz(npjet))*dz)/beadlen
    else
      beadvel=((yvx(ipoint)-yvx(ipoint+1))*dx+ &
       (yvy(ipoint)-yvy(ipoint+1))*dy+(yvz(ipoint)-yvz(ipoint+1))*dz)/beadlen
    endif
    cp=cp0*yvl(ipoint)/yve(ipoint)
    rattao=(cp/cp0)**tev
    ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
    fst(j)=(1.d0/rattao)*(yieldstress+consistency*ratmu* &
     (beadvel/beadlen)**findex-yst(ipoint))
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_stress_3d

 subroutine accelerator_maxwell_evap_stress_3d(firstpoint,lastpoint,npjet, &
   linserting,linserted,jetfr,fev,fst,yxx,yyy,yzz,yvx,yvy,yvz,yst,yvl,yve, &
   evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity,cp0,Bev,mev,tev, &
   consistency,findex,yieldstress)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: linserting,linserted,jetfr(0:)
  double precision, intent(out) :: fev(0:),fst(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:),yvx(0:),yvy(0:),yvz(0:)
  double precision, intent(in) :: yst(0:),yvl(0:),yve(0:)
  double precision, intent(in) :: evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity
  double precision, intent(in) :: cp0,Bev,mev,tev,consistency,findex,yieldstress
  integer :: ipoint,j
  double precision :: dx,dy,dz,beadlen,vnorm,re,beadvel,cp,ratmu,rattao
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yst(0:npjet),yvl(0:npjet), &
!$acc& yve(0:npjet),jetfr(0:npjet)) present_or_copyout(fev(0:lastpoint-firstpoint), &
!$acc& fst(0:lastpoint-firstpoint)) private(j,dx,dy,dz,beadlen,vnorm,re,beadvel, &
!$acc& cp,ratmu,rattao)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
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
    vnorm=dsqrt(yvx(ipoint)*yvx(ipoint)+yvy(ipoint)*yvy(ipoint)+yvz(ipoint)*yvz(ipoint))
    if(beadlen>0.d0 .and. evairv>0.d0 .and. yve(ipoint)>0.d0)then
      re=(2.d0*dsqrt(yve(ipoint)/(3.14159265358979323846d0*beadlen))*vnorm)/evairv
      fev(j)=-evmasscoeff*0.495d0*(re**(1.d0/3.d0))*sqrevsc* &
       (evcsvapour-evumidity)*3.14159265358979323846d0*beadlen
    endif
    if(yve(ipoint)<=0.d0 .or. yvl(ipoint)<=0.d0)cycle
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
    fst(j)=(1.d0/rattao)*(yieldstress+consistency*ratmu* &
     (beadvel/beadlen)**findex-yst(ipoint))
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_maxwell_evap_stress_3d

 subroutine accelerator_kv_stress_3d(firstpoint,lastpoint,npjet,linserted,jetfr, &
   fst,yxx,yyy,yzz,yvx,yvy,yvz,yax,yay,yaz,yst,yvl,yve,fevlocal, &
   cp0,Bev,mev,tev,evlim)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: linserted,jetfr(0:)
  double precision, intent(out) :: fst(0:)
  double precision, intent(in) :: yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(in) :: yvx(0:),yvy(0:),yvz(0:)
  double precision, intent(in) :: yax(0:),yay(0:),yaz(0:),yst(0:)
  double precision, intent(in) :: yvl(0:),yve(0:),fevlocal(0:)
  double precision, intent(in) :: cp0,Bev,mev,tev,evlim
  integer :: ipoint,j
  double precision :: dx,dy,dz,beadlen,beadvel,beadacc
  double precision :: cp,ratmu,ratg,dcpdt,dratmu,dratg,strain,strainrate,strainacc
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet),yzz(0:npjet), &
!$acc& yvx(0:npjet),yvy(0:npjet),yvz(0:npjet),yax(0:npjet),yay(0:npjet), &
!$acc& yaz(0:npjet),yst(0:npjet),yvl(0:npjet),yve(0:npjet), &
!$acc& fevlocal(0:lastpoint-firstpoint),jetfr(0:npjet)) &
!$acc& copyout(fst(0:lastpoint-firstpoint)) private(j,dx,dy,dz,beadlen,beadvel, &
!$acc& beadacc,cp,ratmu,ratg,dcpdt,dratmu,dratg,strain,strainrate,strainacc)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    fst(j)=0.d0
    if(jetfr(ipoint) .or. ipoint>=npjet)cycle
    if(ipoint==npjet-1 .and. .not.linserted)cycle
    if(ipoint==npjet-2 .and. .not.linserted)then
      dx=yxx(ipoint)-yxx(npjet); dy=yyy(ipoint)-yyy(npjet); dz=yzz(ipoint)-yzz(npjet)
      beadvel=((yvx(ipoint)-yvx(npjet))*dx+(yvy(ipoint)-yvy(npjet))*dy+ &
       (yvz(ipoint)-yvz(npjet))*dz)
      beadacc=((yax(ipoint)-yax(npjet))*dx+(yay(ipoint)-yay(npjet))*dy+ &
       (yaz(ipoint)-yaz(npjet))*dz)
    else
      dx=yxx(ipoint)-yxx(ipoint+1); dy=yyy(ipoint)-yyy(ipoint+1); dz=yzz(ipoint)-yzz(ipoint+1)
      beadvel=((yvx(ipoint)-yvx(ipoint+1))*dx+(yvy(ipoint)-yvy(ipoint+1))*dy+ &
       (yvz(ipoint)-yvz(ipoint+1))*dz)
      beadacc=((yax(ipoint)-yax(ipoint+1))*dx+(yay(ipoint)-yay(ipoint+1))*dy+ &
       (yaz(ipoint)-yaz(ipoint+1))*dz)
    endif
    beadlen=dsqrt(dx*dx+dy*dy+dz*dz)
    if(beadlen<=0.d0 .or. yve(ipoint)<=0.d0 .or. yvl(ipoint)<=0.d0)cycle
    strainrate=beadvel/(beadlen*beadlen)
    strainacc=beadacc/(beadlen*beadlen)
    cp=cp0*yvl(ipoint)/yve(ipoint)
    ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
    ratg=ratmu/(cp/cp0)**tev
    dcpdt=0.d0
    if((yve(ipoint)/yvl(ipoint))>evlim*(1.d0+1.d-12))dcpdt=-cp*fevlocal(j)/yve(ipoint)
    dratmu=ratmu*dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))*dcpdt
    dratg=ratg*(dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))-tev/cp)*dcpdt
    strain=(yst(ipoint)-ratmu*strainrate)/ratg
    fst(j)=ratg*strainrate+ratmu*strainacc+dratg*strain+dratmu*strainrate
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_kv_stress_3d

 subroutine accelerator_kv_evap_stress_3d(firstpoint,lastpoint,npjet, &
   linserting,linserted,jetfr,fev,fst,yxx,yyy,yzz,yvx,yvy,yvz,yax,yay,yaz, &
   yst,yvl,yve,evairv,evmasscoeff,sqrevsc,evcsvapour,evumidity, &
   cp0,Bev,mev,tev,evlim,acceleration_chunked)
  implicit none
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
  logical :: chunked_acceleration
  chunked_acceleration=.false.
  if(present(acceleration_chunked))chunked_acceleration=acceleration_chunked
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
    vnorm=dsqrt(yvx(ipoint)*yvx(ipoint)+yvy(ipoint)*yvy(ipoint)+yvz(ipoint)*yvz(ipoint))
    if(beadlen>0.d0 .and. evairv>0.d0 .and. yve(ipoint)>0.d0)then
      re=(2.d0*dsqrt(yve(ipoint)/(3.14159265358979323846d0*beadlen))*vnorm)/evairv
      fev(j)=-evmasscoeff*0.495d0*(re**(1.d0/3.d0))*sqrevsc* &
       (evcsvapour-evumidity)*3.14159265358979323846d0*beadlen
    endif
    if(beadlen<=0.d0 .or. yve(ipoint)<=0.d0 .or. yvl(ipoint)<=0.d0)cycle
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
 end subroutine accelerator_kv_evap_stress_3d

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
  double precision :: dx,dy,dz,norm,cmass1,cmass2,qt,coef,distance
  integer :: ihigh
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet), &
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
#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yvl(0:npjet),yve(0:npjet),jetms(0:npjet), &
!$acc& jetch(0:npjet),jetfr(0:npjet),coulcrossec(0:npjet)) &
!$acc& present_or_copyout(ycf(0:npjet,1:3)) private(jpoint,ihigh,dx,dy,dz,norm, &
!$acc& cmass1,cmass2,qt,coef)
#endif
  do ipoint=inpjet,npjet
    ycf(ipoint,1:3)=0.d0
    if(jetfr(ipoint))cycle
    cmass1=yve(ipoint)/yvl(ipoint)
    qt=jetch(ipoint)*q
    do jpoint=inpjet,npjet
      if(jpoint==ipoint .or. jetfr(jpoint))cycle
      dx=yxx(ipoint)-yxx(jpoint)
      dy=yyy(ipoint)-yyy(jpoint)
      dz=yzz(ipoint)-yzz(jpoint)
      norm=dsqrt(dx*dx+dy*dy+dz*dz)
      ! Match the host Maxwell 3D evaporation path: the ordinary bead-bead
      ! interaction does not apply the cutoff. The mirror interaction below
      ! retains the host cutoff behavior.
      if(norm>1.d-30)then
        cmass2=yve(jpoint)/yvl(jpoint)
        ! The host pair loop stores the cross-section of the higher-index
        ! bead for both directions.  Use the same symmetric lookup here.
        ihigh=max(ipoint,jpoint)
        coef=(jetch(jpoint)*qt)/((norm+coulcrossec(ihigh))**2.d0)
        ycf(ipoint,1)=ycf(ipoint,1)+coef/(jetms(ipoint)*cmass1)*dx/norm
        ycf(ipoint,2)=ycf(ipoint,2)+coef/(jetms(ipoint)*cmass1)*dy/norm
        ycf(ipoint,3)=ycf(ipoint,3)+coef/(jetms(ipoint)*cmass1)*dz/norm
      endif
    enddo
    if(lmirror)then
      do jpoint=inpjet,npjet
        if(jetfr(jpoint))cycle
        dx=yxx(ipoint)-(dabs(yxx(jpoint)-h)+h)
        dy=yyy(ipoint)-yyy(jpoint)
        dz=yzz(ipoint)-yzz(jpoint)
        norm=dsqrt(dx*dx+dy*dy+dz*dz)
        if(ldcutoff .and. norm>dcutoff)cycle
        if(norm>1.d-30)then
          cmass2=yve(jpoint)/yvl(jpoint)
          coef=(jetch(jpoint)*qt)/((norm+coulcrossec(jpoint))**2.d0)
          ycf(ipoint,1)=ycf(ipoint,1)-coef/(jetms(ipoint)*cmass1)*dx/norm
          ycf(ipoint,2)=ycf(ipoint,2)-coef/(jetms(ipoint)*cmass1)*dy/norm
          ycf(ipoint,3)=ycf(ipoint,3)-coef/(jetms(ipoint)*cmass1)*dz/norm
        endif
      enddo
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_coulomb_evap_3d

 subroutine accelerator_coulomb_evap_compare(npjet,inpjet,ycf,yxx,yyy,yzz, &
   yvl,yve,jetms,jetch,jetfr,coulcrossec,q,lmirror,h,ldcutoff,dcutoff)
  implicit none
  integer, intent(in) :: npjet,inpjet
  double precision, intent(in) :: ycf(0:,1:),yxx(0:),yyy(0:),yzz(0:)
  double precision, intent(in) :: yvl(0:),yve(0:),jetms(0:),jetch(0:)
  double precision, intent(in) :: coulcrossec(0:),q,h,dcutoff
  logical, intent(in) :: jetfr(0:),lmirror,ldcutoff
  double precision, allocatable :: ref(:,:)
  double precision :: dx,dy,dz,norm,cmass1,cmass2,qt,coef,xmirror
  double precision :: maxdiff
  integer :: ipoint,jpoint,imax,jmax,k
  allocate(ref(0:npjet,1:3))
  ref=0.d0
  do ipoint=inpjet,npjet
    if(jetfr(ipoint))cycle
    cmass1=yve(ipoint)/yvl(ipoint)
    qt=jetch(ipoint)*q
    do jpoint=ipoint+1,npjet
      if(jetfr(jpoint))cycle
      dx=yxx(ipoint)-yxx(jpoint)
      dy=yyy(ipoint)-yyy(jpoint)
      dz=yzz(ipoint)-yzz(jpoint)
      norm=dsqrt(dx*dx+dy*dy+dz*dz)
      if(norm>1.d-30)then
        cmass2=yve(jpoint)/yvl(jpoint)
        coef=(jetch(jpoint)*qt)/((norm+coulcrossec(jpoint))**2.d0)
        ref(ipoint,1:3)=ref(ipoint,1:3)+coef/(jetms(ipoint)*cmass1)*(/dx,dy,dz/)/norm
        ref(jpoint,1:3)=ref(jpoint,1:3)-coef/(jetms(jpoint)*cmass2)*(/dx,dy,dz/)/norm
      endif
    enddo
    if(lmirror)then
      do jpoint=inpjet,npjet
        if(jetfr(jpoint))cycle
        xmirror=dabs(yxx(jpoint)-h)+h
        dx=yxx(ipoint)-xmirror
        dy=yyy(ipoint)-yyy(jpoint)
        dz=yzz(ipoint)-yzz(jpoint)
        norm=dsqrt(dx*dx+dy*dy+dz*dz)
        if(ldcutoff .and. norm>dcutoff)cycle
        if(norm>1.d-30)then
          coef=(jetch(jpoint)*qt)/((norm+coulcrossec(jpoint))**2.d0)
          ref(ipoint,1:3)=ref(ipoint,1:3)-coef/(jetms(ipoint)*cmass1)*(/dx,dy,dz/)/norm
        endif
      enddo
    endif
  enddo
  maxdiff=0.d0; imax=-1; jmax=-1
  do ipoint=inpjet,npjet
    do k=1,3
      if(dabs(ycf(ipoint,k)-ref(ipoint,k))>maxdiff)then
        maxdiff=dabs(ycf(ipoint,k)-ref(ipoint,k)); imax=ipoint; jmax=k
      endif
    enddo
  enddo
  write(*,'(A,ES14.6,A,I0,A,I0)') 'COULOMB_DIAGNOSTIC maxdiff=',maxdiff, &
   ' bead=',imax,' component=',jmax
  if(maxdiff>1.d-10 .and. imax>=inpjet)then
    write(*,'(A,2ES24.16,A,4ES24.16)') 'COULOMB_DIAGNOSTIC gpu/ref=', &
     ycf(imax,jmax),ref(imax,jmax),' state x/yve/vefrac/cross=', &
     yxx(imax),yve(imax),yve(imax)/yvl(imax),coulcrossec(imax)
  endif
  deallocate(ref)
 end subroutine accelerator_coulomb_evap_compare

 subroutine accelerator_release_jet_capacity(mxnpjet,jetxx,jetyy,jetzz,jetst,jetvx, &
   jetvy,jetvz,jetms,jetch,jetvl,jetfr)
  implicit none
  integer, intent(in) :: mxnpjet
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
#ifdef _OPENACC
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
!$acc update device(jetve(0:npjet),jetce(0:npjet))
#endif
 end subroutine accelerator_update_device_evaporation_state

 subroutine accelerator_platen_predict(firstpoint,lastpoint,h,airamp,noisediff, &
   jetms,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
   f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
   y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz, &
   y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  double precision, intent(in) :: h,airamp,noisediff,jetms(0:)
  double precision, intent(in) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1st(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:)
  double precision, intent(out) :: y1xx(0:),y1yy(0:),y1zz(0:),y1st(0:)
  double precision, intent(out) :: y1vx(0:),y1vy(0:),y1vz(0:)
  double precision, intent(out) :: y2xx(0:),y2yy(0:),y2zz(0:),y2st(0:)
  double precision, intent(out) :: y2vx(0:),y2vy(0:),y2vz(0:)
  integer :: ipoint,j
  double precision :: dsqrh,stoc
  dsqrh=dsqrt(dabs(h))
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetms,jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
!$acc& y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y2xx,y2yy,y2zz,y2st, &
!$acc& y2vx,y2vy,y2vz) private(j,stoc)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    stoc=dsqrt(2.d0*(airamp/jetms(ipoint)+noisediff))
    if(ipoint==lastpoint)stoc=0.d0
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
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_predict

 subroutine accelerator_platen_evap_predict(firstpoint,lastpoint,h,airamp, &
   noisediff,evlim,jetms,jetvl,jetve,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz,f1ev, &
   y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y1ev, &
   y2xx,y2yy,y2zz,y2st,y2vx,y2vy,y2vz,y2ev)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
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
!$acc parallel loop gang vector present(jetms,jetvl,jetve,jetxx,jetyy, &
!$acc& jetzz,jetst,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy, &
!$acc& f1vz,f1ev,y1xx,y1yy,y1zz,y1st,y1vx,y1vy,y1vz,y1ev,y2xx,y2yy, &
!$acc& y2zz,y2st,y2vx,y2vy,y2vz,y2ev) private(j,stoc,cmass,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    cmass=jetve(ipoint)/jetvl(ipoint)
    stoc=dsqrt(2.d0*(airamp/(jetms(ipoint)*cmass)+noisediff))
    if(ipoint==lastpoint)stoc=0.d0
    y1xx(ipoint)=jetxx(ipoint)+h*f1xx(j)
    y1yy(ipoint)=jetyy(ipoint)+h*f1yy(j)
    y1zz(ipoint)=jetzz(ipoint)+h*f1zz(j)
    y1st(ipoint)=jetst(ipoint)+h*f1st(j)
    y1vx(ipoint)=jetvx(ipoint)+h*f1vx(j)+dsqrh*stoc
    y1vy(ipoint)=jetvy(ipoint)+h*f1vy(j)+dsqrh*stoc
    y1vz(ipoint)=jetvz(ipoint)+h*f1vz(j)+dsqrh*stoc
    ve=jetve(ipoint)+h*f1ev(j)
    if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
    y1ev(ipoint)=ve
    y2xx(ipoint)=y1xx(ipoint)
    y2yy(ipoint)=y1yy(ipoint)
    y2zz(ipoint)=y1zz(ipoint)
    y2st(ipoint)=y1st(ipoint)
    y2vx(ipoint)=jetvx(ipoint)+h*f1vx(j)-dsqrh*stoc
    y2vy(ipoint)=jetvy(ipoint)+h*f1vy(j)-dsqrh*stoc
    y2vz(ipoint)=jetvz(ipoint)+h*f1vz(j)-dsqrh*stoc
    y2ev(ipoint)=ve
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_evap_predict

 subroutine accelerator_platen_velocity(firstpoint,lastpoint,mxnpjet, &
   historysteps,k,h, &
   airamp,noisediff,jetms,gaussianhistory,jetvx,jetvy,jetvz, &
   f1vx,f1vy,f1vz,f2vx,f2vy,f2vz,f3vx,f3vy,f3vz)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,mxnpjet,historysteps,k
  double precision, intent(in) :: h,airamp,noisediff,jetms(0:)
  double precision, intent(in) :: gaussianhistory(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:)
  double precision, intent(in) :: f2vx(0:),f2vy(0:),f2vz(0:)
  double precision, intent(in) :: f3vx(0:),f3vy(0:),f3vz(0:)
  integer :: ipoint,j,component,nperstep,index1,index2,cycle_step
  double precision :: dsqrh,tsqh,prefactor,stoc,u1,u2,ww,zz
  dsqrh=dsqrt(dabs(h)); tsqh=dsqrh**3.d0; prefactor=0.5d0/dsqrh
  nperstep=(mxnpjet+1)*6
  cycle_step=mod(k-1,historysteps)+1
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetms,gaussianhistory,jetvx,jetvy, &
!$acc& jetvz,f1vx,f1vy,f1vz,f2vx,f2vy,f2vz,f3vx,f3vy,f3vz) &
!$acc& private(j,component,index1,index2,stoc,u1,u2,ww,zz)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    stoc=dsqrt(2.d0*(airamp/jetms(ipoint)+noisediff))
    if(ipoint==lastpoint)stoc=0.d0
    component=1
    index1=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1)
    index2=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1+3)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvx(ipoint)=jetvx(ipoint)+stoc*ww+prefactor*(f2vx(j)-f3vx(j))*zz+ &
     0.25d0*h*(f2vx(j)+2.d0*f1vx(j)+f3vx(j))
    component=2
    index1=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1)
    index2=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1+3)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvy(ipoint)=jetvy(ipoint)+stoc*ww+prefactor*(f2vy(j)-f3vy(j))*zz+ &
     0.25d0*h*(f2vy(j)+2.d0*f1vy(j)+f3vy(j))
    component=3
    index1=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1)
    index2=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1+3)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvz(ipoint)=jetvz(ipoint)+stoc*ww+prefactor*(f2vz(j)-f3vz(j))*zz+ &
     0.25d0*h*(f2vz(j)+2.d0*f1vz(j)+f3vz(j))
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_velocity

 subroutine accelerator_platen_evap_velocity(firstpoint,lastpoint,mxnpjet, &
   historysteps,k,h,airamp,noisediff,jetms,jetvl,jetve,gaussianhistory, &
   jetvx,jetvy,jetvz,f1vx,f1vy,f1vz,f2vx,f2vy,f2vz,f3vx,f3vy,f3vz)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,mxnpjet,historysteps,k
  double precision, intent(in) :: h,airamp,noisediff
  double precision, intent(in) :: jetms(0:),jetvl(0:),jetve(0:)
  double precision, intent(in) :: gaussianhistory(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:)
  double precision, intent(in) :: f2vx(0:),f2vy(0:),f2vz(0:)
  double precision, intent(in) :: f3vx(0:),f3vy(0:),f3vz(0:)
  integer :: ipoint,j,component,nperstep,index1,index2,cycle_step
  double precision :: dsqrh,tsqh,prefactor,stoc,cmass,u1,u2,ww,zz
  dsqrh=dsqrt(dabs(h)); tsqh=dsqrh**3.d0; prefactor=0.5d0/dsqrh
  nperstep=(mxnpjet+1)*6
  cycle_step=mod(k-1,historysteps)+1
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetms,jetvl,jetve,gaussianhistory, &
!$acc& jetvx,jetvy,jetvz,f1vx,f1vy,f1vz,f2vx,f2vy,f2vz,f3vx,f3vy,f3vz) &
!$acc& private(j,component,index1,index2,stoc,cmass,u1,u2,ww,zz)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    cmass=jetve(ipoint)/jetvl(ipoint)
    stoc=dsqrt(2.d0*(airamp/(jetms(ipoint)*cmass)+noisediff))
    if(ipoint==lastpoint)stoc=0.d0
    component=1
    index1=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1)
    index2=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1+3)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvx(ipoint)=jetvx(ipoint)+stoc*ww+prefactor*(f2vx(j)-f3vx(j))*zz+ &
     0.25d0*h*(f2vx(j)+2.d0*f1vx(j)+f3vx(j))
    component=2
    index1=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1)
    index2=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1+3)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvy(ipoint)=jetvy(ipoint)+stoc*ww+prefactor*(f2vy(j)-f3vy(j))*zz+ &
     0.25d0*h*(f2vy(j)+2.d0*f1vy(j)+f3vy(j))
    component=3
    index1=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1)
    index2=(cycle_step-1)*nperstep+ipoint+(mxnpjet+1)*(component-1+3)
    u1=gaussianhistory(index1); u2=gaussianhistory(index2)
    ww=dsqrh*u1; zz=0.5d0*tsqh*(u1+u2/dsqrt(3.d0))
    jetvz(ipoint)=jetvz(ipoint)+stoc*ww+prefactor*(f2vz(j)-f3vz(j))*zz+ &
     0.25d0*h*(f2vz(j)+2.d0*f1vz(j)+f3vz(j))
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_evap_velocity

 subroutine accelerator_platen_positions(firstpoint,lastpoint,npjet,h,pfreq, &
   liniperturb,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,f1xx,f1yy,f1zz)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: liniperturb
  double precision, intent(in) :: h,pfreq,jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:)
  integer :: ipoint,j
  double precision :: f2x,f2y,f2z,y1y,y1z
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetvx,jetvy, &
!$acc& jetvz,f1xx,f1yy,f1zz) private(j,f2x,f2y,f2z,y1y,y1z)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    f2x=jetvx(ipoint); f2y=jetvy(ipoint); f2z=jetvz(ipoint)
    if(ipoint==npjet)then
      f2x=0.d0
      if(liniperturb)then
        y1y=jetyy(ipoint)+h*f1yy(j)
        y1z=jetzz(ipoint)+h*f1zz(j)
        f2y=-pfreq*y1z; f2z=pfreq*y1y
      else
        f2y=0.d0; f2z=0.d0
      endif
    endif
    jetxx(ipoint)=jetxx(ipoint)+0.5d0*h*(f1xx(j)+f2x)
    jetyy(ipoint)=jetyy(ipoint)+0.5d0*h*(f1yy(j)+f2y)
    jetzz(ipoint)=jetzz(ipoint)+0.5d0*h*(f1zz(j)+f2z)
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_positions

 subroutine accelerator_platen_evap_positions(firstpoint,lastpoint,npjet,h, &
   pfreq,liniperturb,evlim,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetvl, &
   jetve,f1xx,f1yy,f1zz,f1ev,f2ev)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet
  logical, intent(in) :: liniperturb
  double precision, intent(in) :: h,pfreq,evlim
  double precision, intent(in) :: jetvx(0:),jetvy(0:),jetvz(0:),jetvl(0:)
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetve(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1ev(0:)
  double precision, intent(in) :: f2ev(0:)
  integer :: ipoint,j
  double precision :: f2x,f2y,f2z,y1y,y1z,ve
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetvx,jetvy, &
!$acc& jetvz,jetvl,jetve,f1xx,f1yy,f1zz,f1ev,f2ev) &
!$acc& private(j,f2x,f2y,f2z,y1y,y1z,ve)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    f2x=jetvx(ipoint); f2y=jetvy(ipoint); f2z=jetvz(ipoint)
    if(ipoint==npjet)then
      f2x=0.d0
      if(liniperturb)then
        y1y=jetyy(ipoint)+h*f1yy(j)
        y1z=jetzz(ipoint)+h*f1zz(j)
        f2y=-pfreq*y1z; f2z=pfreq*y1y
      else
        f2y=0.d0; f2z=0.d0
      endif
    endif
    jetxx(ipoint)=jetxx(ipoint)+0.5d0*h*(f1xx(j)+f2x)
    jetyy(ipoint)=jetyy(ipoint)+0.5d0*h*(f1yy(j)+f2y)
    jetzz(ipoint)=jetzz(ipoint)+0.5d0*h*(f1zz(j)+f2z)
    ve=jetve(ipoint)+0.5d0*h*(f1ev(j)+f2ev(j))
    if(ve/jetvl(ipoint)<evlim)ve=jetvl(ipoint)*evlim
    jetve(ipoint)=ve
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_platen_evap_positions

 subroutine accelerator_platen_stress_statistics(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,f1st,f2st,counterlpath,ncounterlpath, &
   maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  integer, intent(inout) :: ncounterlpath
  double precision, intent(in) :: h,jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:)
  double precision, intent(in) :: f1st(0:),f2st(0:)
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint,j
  double precision :: newst,dx,dy,dz
  call accelerator_map_statistics(counterlpath,ncounterlpath,maxstress,maxstressposx)
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst,f1st,f2st, &
!$acc& counterlpath,statistics_step_max) private(j,newst,dx,dy,dz) &
!$acc& reduction(+:counterlpath) reduction(max:statistics_step_max)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    newst=jetst(ipoint)+0.5d0*h*(f1st(j)+f2st(j))
    if(ipoint<lastpoint)then
      dx=jetxx(ipoint)-jetxx(ipoint+1)
      dy=jetyy(ipoint)-jetyy(ipoint+1)
      dz=jetzz(ipoint)-jetzz(ipoint+1)
      counterlpath=counterlpath+dsqrt(dx*dx+dy*dy+dz*dz)
    endif
    jetst(ipoint)=newst
    statistics_step_max=max(statistics_step_max,newst)
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
  call accelerator_store_statistics(firstpoint,lastpoint,jetxx,jetyy, &
   jetzz,jetst,counterlpath,ncounterlpath,maxstress,maxstressposx)
 end subroutine accelerator_platen_stress_statistics

 subroutine accelerator_euler_final_statistics(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,counterlpath,ncounterlpath, &
   maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  integer, intent(inout) :: ncounterlpath
  double precision, intent(in) :: h
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: fxx(0:),fyy(0:),fzz(0:),fst(0:)
  double precision, intent(in) :: fvx(0:),fvy(0:),fvz(0:)
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint,j,jn
  double precision :: newxx,newyy,newzz,newst,nextxx,nextyy,nextzz,dx,dy,dz

  if(.not.(accelerator_persistent .or. accelerator_topology_enabled))return
  call accelerator_map_statistics(counterlpath,ncounterlpath,maxstress, &
   maxstressposx)
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,fxx,fyy,fzz,fst,fvx,fvy,fvz, &
!$acc& counterlpath,statistics_step_max) &
!$acc& private(j,jn,newxx,newyy,newzz,newst,nextxx,nextyy,nextzz, &
!$acc& dx,dy,dz) reduction(+:counterlpath) &
!$acc& reduction(max:statistics_step_max)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    newxx=jetxx(ipoint)+h*fxx(j)
    newyy=jetyy(ipoint)+h*fyy(j)
    newzz=jetzz(ipoint)+h*fzz(j)
    newst=jetst(ipoint)+h*fst(j)
    if(ipoint<lastpoint)then
      jn=j+1
      nextxx=jetxx(ipoint+1)+h*fxx(jn)
      nextyy=jetyy(ipoint+1)+h*fyy(jn)
      nextzz=jetzz(ipoint+1)+h*fzz(jn)
      dx=newxx-nextxx
      dy=newyy-nextyy
      dz=newzz-nextzz
      counterlpath=counterlpath+dsqrt(dx*dx+dy*dy+dz*dz)
    endif
    statistics_step_max=max(statistics_step_max,newst)
    jetxx(ipoint)=newxx
    jetyy(ipoint)=newyy
    jetzz(ipoint)=newzz
    jetst(ipoint)=newst
    jetvx(ipoint)=jetvx(ipoint)+h*fvx(j)
    jetvy(ipoint)=jetvy(ipoint)+h*fvy(j)
    jetvz(ipoint)=jetvz(ipoint)+h*fvz(j)
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_euler_final_statistics

 subroutine accelerator_rk2_final_statistics(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
   f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
   f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
   counterlpath,ncounterlpath,maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  integer, intent(inout) :: ncounterlpath
  double precision, intent(in) :: h
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1st(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:)
  double precision, intent(in) :: f2xx(0:),f2yy(0:),f2zz(0:),f2st(0:)
  double precision, intent(in) :: f2vx(0:),f2vy(0:),f2vz(0:)
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint,j,jn
  double precision :: scale,newxx,newyy,newzz,newst,nextxx,nextyy,nextzz
  double precision :: dx,dy,dz

  if(.not.accelerator_persistent)return
  call accelerator_map_statistics(counterlpath,ncounterlpath,maxstress, &
   maxstressposx)
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,counterlpath, &
!$acc& statistics_step_max) private(j,jn,scale,newxx,newyy,newzz, &
!$acc& newst,nextxx,nextyy,nextzz,dx,dy,dz) &
!$acc& reduction(+:counterlpath) reduction(max:statistics_step_max)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    scale=h/2.d0
    newxx=jetxx(ipoint)+scale*(f1xx(j)+f2xx(j))
    newyy=jetyy(ipoint)+scale*(f1yy(j)+f2yy(j))
    newzz=jetzz(ipoint)+scale*(f1zz(j)+f2zz(j))
    newst=jetst(ipoint)+scale*(f1st(j)+f2st(j))
    if(ipoint<lastpoint)then
      jn=j+1
      nextxx=jetxx(ipoint+1)+scale*(f1xx(jn)+f2xx(jn))
      nextyy=jetyy(ipoint+1)+scale*(f1yy(jn)+f2yy(jn))
      nextzz=jetzz(ipoint+1)+scale*(f1zz(jn)+f2zz(jn))
      dx=newxx-nextxx
      dy=newyy-nextyy
      dz=newzz-nextzz
      counterlpath=counterlpath+dsqrt(dx*dx+dy*dy+dz*dz)
    endif
    statistics_step_max=max(statistics_step_max,newst)
    jetxx(ipoint)=newxx
    jetyy(ipoint)=newyy
    jetzz(ipoint)=newzz
    jetst(ipoint)=newst
    jetvx(ipoint)=jetvx(ipoint)+scale*(f1vx(j)+f2vx(j))
    jetvy(ipoint)=jetvy(ipoint)+scale*(f1vy(j)+f2vy(j))
    jetvz(ipoint)=jetvz(ipoint)+scale*(f1vz(j)+f2vz(j))
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
 end subroutine accelerator_rk2_final_statistics

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

 subroutine accelerator_rk4_final_statistics(firstpoint,lastpoint,h, &
   jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
   f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
   f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz, &
   f3xx,f3yy,f3zz,f3st,f3vx,f3vy,f3vz, &
   f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz, &
   counterlpath,ncounterlpath,maxstress,maxstressposx)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint
  integer, intent(inout) :: ncounterlpath
  double precision, intent(in) :: h
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(in) :: f1xx(0:),f1yy(0:),f1zz(0:),f1st(0:)
  double precision, intent(in) :: f1vx(0:),f1vy(0:),f1vz(0:)
  double precision, intent(in) :: f2xx(0:),f2yy(0:),f2zz(0:),f2st(0:)
  double precision, intent(in) :: f2vx(0:),f2vy(0:),f2vz(0:)
  double precision, intent(in) :: f3xx(0:),f3yy(0:),f3zz(0:),f3st(0:)
  double precision, intent(in) :: f3vx(0:),f3vy(0:),f3vz(0:)
  double precision, intent(in) :: f4xx(0:),f4yy(0:),f4zz(0:),f4st(0:)
  double precision, intent(in) :: f4vx(0:),f4vy(0:),f4vz(0:)
  double precision, intent(inout) :: counterlpath,maxstress,maxstressposx
  integer :: ipoint,j,jn
  double precision :: scale,newxx,newyy,newzz,newst,nextxx,nextyy,nextzz
  double precision :: dx,dy,dz

  if(.not.accelerator_persistent)return
  call accelerator_map_statistics(counterlpath,ncounterlpath,maxstress, &
   maxstressposx)
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetxx,jetyy,jetzz,jetst, &
!$acc& jetvx,jetvy,jetvz,f1xx,f1yy,f1zz,f1st,f1vx,f1vy,f1vz, &
!$acc& f2xx,f2yy,f2zz,f2st,f2vx,f2vy,f2vz,f3xx,f3yy,f3zz, &
!$acc& f3st,f3vx,f3vy,f3vz,f4xx,f4yy,f4zz,f4st,f4vx,f4vy,f4vz, &
!$acc& counterlpath,statistics_step_max) &
!$acc& private(j,jn,scale,newxx,newyy,newzz,newst,nextxx,nextyy, &
!$acc& nextzz,dx,dy,dz) reduction(+:counterlpath) &
!$acc& reduction(max:statistics_step_max)
#endif
  do ipoint=firstpoint,lastpoint
    j=ipoint-firstpoint
    scale=h/6.d0
    newxx=jetxx(ipoint)+scale*(f1xx(j)+2.d0*(f2xx(j)+f3xx(j))+f4xx(j))
    newyy=jetyy(ipoint)+scale*(f1yy(j)+2.d0*(f2yy(j)+f3yy(j))+f4yy(j))
    newzz=jetzz(ipoint)+scale*(f1zz(j)+2.d0*(f2zz(j)+f3zz(j))+f4zz(j))
    newst=jetst(ipoint)+scale*(f1st(j)+2.d0*(f2st(j)+f3st(j))+f4st(j))
    if(ipoint<lastpoint)then
      jn=j+1
      nextxx=jetxx(ipoint+1)+scale*(f1xx(jn)+ &
       2.d0*(f2xx(jn)+f3xx(jn))+f4xx(jn))
      nextyy=jetyy(ipoint+1)+scale*(f1yy(jn)+ &
       2.d0*(f2yy(jn)+f3yy(jn))+f4yy(jn))
      nextzz=jetzz(ipoint+1)+scale*(f1zz(jn)+ &
       2.d0*(f2zz(jn)+f3zz(jn))+f4zz(jn))
      dx=newxx-nextxx
      dy=newyy-nextyy
      dz=newzz-nextzz
      counterlpath=counterlpath+dsqrt(dx*dx+dy*dy+dz*dz)
    endif
    statistics_step_max=max(statistics_step_max,newst)
    jetxx(ipoint)=newxx
    jetyy(ipoint)=newyy
    jetzz(ipoint)=newzz
    jetst(ipoint)=newst
    jetvx(ipoint)=jetvx(ipoint)+scale*(f1vx(j)+ &
     2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
    jetvy(ipoint)=jetvy(ipoint)+scale*(f1vy(j)+ &
     2.d0*(f2vy(j)+f3vy(j))+f4vy(j))
    jetvz(ipoint)=jetvz(ipoint)+scale*(f1vz(j)+ &
     2.d0*(f2vz(j)+f3vz(j))+f4vz(j))
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif
  return
 end subroutine accelerator_rk4_final_statistics

 subroutine accelerator_set_persistent(enabled)
  implicit none
  logical, intent(in) :: enabled
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
!$acc update self(jetve(0:npjet),jetce(0:npjet))
#endif
 end subroutine accelerator_update_host_evaporation_state

 subroutine accelerator_add_bead(npjet,mxnpjet,linserted,ladd,lresize, &
   resolution,dresolution,thresolution,ivelocity,istress,imassa,icharge, &
   ivolume,jetxx,jetyy,jetzz,jetst,jetvx,jetvy, &
   jetvz,jetms,jetch,jetvl,jetfr)
  implicit none
  integer, intent(inout) :: npjet
  integer, intent(in) :: mxnpjet
  logical, intent(inout) :: linserted,ladd,lresize
  double precision, intent(in) :: resolution,dresolution,thresolution
  double precision, intent(in) :: ivelocity,istress,imassa,icharge,ivolume
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
  double precision :: dx,dy,dz,distance,scale

#ifdef _OPENACC
!$acc serial present(jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,jetms, &
!$acc& jetch,jetvl,jetfr) copy(npjet,linserted) copyout(ladd,lresize) &
!$acc& private(dx,dy,dz,distance,scale)
#endif
  ladd=.false.
  lresize=.false.
  if(.not.linserted)then
    dx=jetxx(npjet-2)-jetxx(npjet)
    dy=jetyy(npjet-2)-jetyy(npjet)
    dz=jetzz(npjet-2)-jetzz(npjet)
    distance=dsqrt(dx*dx+dy*dy+dz*dz)
    if(distance>=dresolution)then
      linserted=.true.
      jetst(npjet-1)=0.d0
      jetvx(npjet-1)=ivelocity
      jetvy(npjet-1)=0.d0
      jetvz(npjet-1)=0.d0
    endif
  else
    dx=jetxx(npjet-1)-jetxx(npjet)
    dy=jetyy(npjet-1)-jetyy(npjet)
    dz=jetzz(npjet-1)-jetzz(npjet)
    distance=dsqrt(dx*dx+dy*dy+dz*dz)
    if(distance>=thresolution .and. npjet>=mxnpjet)then
      lresize=.true.
    elseif(distance>=thresolution)then
      npjet=npjet+1
      jetfr(npjet)=jetfr(npjet-1)
      jetxx(npjet)=jetxx(npjet-1)
      jetyy(npjet)=jetyy(npjet-1)
      jetzz(npjet)=jetzz(npjet-1)
      jetst(npjet)=jetst(npjet-1)
      jetvx(npjet)=jetvx(npjet-1)
      jetvy(npjet)=jetvy(npjet-1)
      jetvz(npjet)=jetvz(npjet-1)
      jetms(npjet)=jetms(npjet-1)
      jetch(npjet)=jetch(npjet-1)
      jetvl(npjet)=jetvl(npjet-1)
      jetfr(npjet-1)=.false.
      jetst(npjet-1)=istress
      jetvx(npjet-1)=ivelocity
      jetvy(npjet-1)=0.d0
      jetvz(npjet-1)=0.d0
      jetms(npjet-1)=imassa*ivolume
      jetch(npjet-1)=icharge*ivolume
      jetvl(npjet-1)=ivolume
      dx=jetxx(npjet-2)-jetxx(npjet)
      dy=jetyy(npjet-2)-jetyy(npjet)
      dz=jetzz(npjet-2)-jetzz(npjet)
      distance=dsqrt(dx*dx+dy*dy+dz*dz)
      scale=resolution/distance
      jetxx(npjet-1)=jetxx(npjet)+scale*dx
      jetyy(npjet-1)=jetyy(npjet)+scale*dy
      jetzz(npjet-1)=jetzz(npjet)+scale*dz
      ladd=.true.
      linserted=.false.
    endif
  endif
#ifdef _OPENACC
!$acc end serial
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
 end subroutine accelerator_add_bead

 subroutine accelerator_update_device_added_evaporation(npjet,ivolume,jetve,jetce)
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(in) :: ivolume
  double precision, intent(inout) :: jetve(0:),jetce(0:)
#ifdef _OPENACC
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

 subroutine accelerator_remove_bead(inpjet,npjet,h,jetxx,jetyy,jetzz,jetst, &
   jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr,nremoved,lrem)
  implicit none
  integer, intent(inout) :: inpjet
  integer, intent(in) :: npjet
  integer, intent(out) :: nremoved
  double precision, intent(in) :: h
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:),jetst(0:)
  double precision, intent(inout) :: jetvx(0:),jetvy(0:),jetvz(0:)
  double precision, intent(inout) :: jetms(0:),jetch(0:),jetvl(0:)
  logical, intent(inout) :: jetfr(0:)
  logical, intent(out) :: lrem
  integer :: ipoint,remove_one
#ifdef _OPENACC
!$acc parallel loop present(jetxx,jetfr)
#endif
  do ipoint=inpjet,npjet
    if(jetxx(ipoint)>=h)then
      jetfr(ipoint)=.true.
      jetxx(ipoint)=h
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
!$acc serial present(jetxx) copyout(remove_one)
#endif
  remove_one=0
  if(jetxx(inpjet)>=h .and. jetxx(inpjet+1)>=h)remove_one=1
#ifdef _OPENACC
!$acc end serial
#endif
  nremoved=remove_one
  lrem=remove_one==1
  if(lrem)inpjet=inpjet+1
#ifdef _OPENACC
  if(lrem)then
! Removal observables need both the removed bead and its active neighbour.
!$acc update self(jetxx(inpjet-1:inpjet),jetyy(inpjet-1:inpjet), &
!$acc& jetzz(inpjet-1:inpjet), &
!$acc& jetst(inpjet-1),jetvx(inpjet-1),jetvy(inpjet-1),jetvz(inpjet-1), &
!$acc& jetms(inpjet-1),jetch(inpjet-1),jetvl(inpjet-1),jetfr(inpjet-1))
  endif
#endif
 end subroutine accelerator_remove_bead

 subroutine accelerator_update_host_removed_evaporation(ipoint,jetve,jetce)
  implicit none
  integer, intent(in) :: ipoint
  double precision, intent(inout) :: jetve(0:),jetce(0:)
#ifdef _OPENACC
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

  if(.not.accelerator_persistent)return
#ifdef _OPENACC
!$acc parallel loop gang vector present(jetst,statistics_step_max, &
!$acc& statistics_step_index) reduction(max:statistics_step_index)
#endif
  do ipoint=inpjet,npjet
    if(jetst(ipoint)==statistics_step_max)then
      statistics_step_index=max(statistics_step_index,ipoint)
    endif
  enddo
#ifdef _OPENACC
!$acc end parallel loop
!$acc serial present(jetxx,counterlpath,ncounterlpath,maxstress, &
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
   apply_airdrag,collector_curvature)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,nfieldtype
  logical, intent(in) :: linserted,liniperturb,lairdrag,lflorentz,luppot
  logical, intent(in) :: stochastic_model
  double precision, intent(in) :: pfreq,consistency,findex,yieldstress
  double precision, intent(in) :: att,fve,gr,ks,li,vfield,velext,noisefric
  double precision, intent(in), optional :: yve(0:)
  logical, intent(in), optional :: apply_airdrag,collector_curvature
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
  if(.not.lairdrag .or. lflorentz .or. luppot)return
  if(nfieldtype/=0 .or. lastpoint/=npjet)return
  nout=lastpoint-firstpoint
  ! Evaluate OPTIONAL presence on the host and pass a plain scalar into the
  ! device kernel; PRESENT() itself is not reliable inside OpenACC regions.
  use_evap=present(yve)
  use_airdrag=lairdrag
  if(present(apply_airdrag))use_airdrag=apply_airdrag
  use_collector_curvature=.false.
  if(present(collector_curvature))use_collector_curvature=collector_curvature

#ifdef _OPENACC
!$acc parallel loop gang vector present_or_copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yst(0:npjet),yvx(0:npjet),yvy(0:npjet), &
!$acc& yvz(0:npjet)) present(yvl,jetms,jetch,jetfr) &
#if defined(JETSPIN_DEV_HOST_COULOMB_ORACLE) || defined(JETSPIN_DEV_HOST_FORCE_ORACLE)
!$acc& present(ycf) &
#else
!$acc& present_or_copyin(ycf) &
#endif
!$acc& present(yve) &
!$acc& present_or_copyout(fxx,fyy,fzz,fst,fvx,fvy,fvz) &
!$acc& private(j,dxu,dyu,dzu,dxd,dyd,dzd,lup,ldown,tux,tuy,tuz, &
!$acc& tdx,tdy,tdz,beadvel,v1x,v1y,v1z,v2x,v2y,v2z,l1,l2,dotp, &
!$acc& nbx,nby,nbz,lnb,b,c,t,scale1,scale2,ccx,ccy,ccz,rcx,rcy, &
!$acc& rcz,radius,curvature,factor1,factor2,factor3,factor4,factor5, &
!$acc& fvet,kst,attt,lit,veltangent,cmass,fvolume,fvolume_prev,straight)
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
    fst(j)=yieldstress+consistency*(beadvel/lup)**findex-yst(ipoint)
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
