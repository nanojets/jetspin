module accelerator_mod

 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.true.
 logical, save :: accelerator_persistent=.false.
 logical, save :: accelerator_statistics_mapped=.false.
 double precision, save :: statistics_step_max=-huge(0.d0)
 integer, save :: statistics_step_index=-1

 public :: accelerator_prepare
 public :: accelerator_eom3_stage
 public :: accelerator_set_persistent
 public :: accelerator_is_persistent
 public :: accelerator_update_host_state
 public :: accelerator_store_statistics
 public :: accelerator_update_host_statistics
 public :: accelerator_update_device_statistics
 public :: accelerator_rk4_final_statistics

contains

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
#ifdef _OPENACC
  if(.not.accelerator_statistics_mapped)then
!$acc enter data copyin(counterlpath,ncounterlpath,maxstress, &
!$acc& maxstressposx,statistics_step_max,statistics_step_index)
    accelerator_statistics_mapped=.true.
  endif
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
  return
 end subroutine accelerator_set_persistent

 logical function accelerator_is_persistent()
  implicit none
  accelerator_is_persistent=accelerator_persistent
 end function accelerator_is_persistent

 subroutine accelerator_update_host_state(npjet,jetxx,jetyy,jetzz, &
   jetst,jetvx,jetvy,jetvz)
  implicit none
  integer, intent(in) :: npjet
  double precision, intent(inout) :: jetxx(0:),jetyy(0:),jetzz(0:)
  double precision, intent(inout) :: jetst(0:),jetvx(0:),jetvy(0:),jetvz(0:)
  if(.not.accelerator_persistent)return
#ifdef _OPENACC
!$acc update self(jetxx(0:npjet),jetyy(0:npjet),jetzz(0:npjet), &
!$acc& jetst(0:npjet),jetvx(0:npjet),jetvy(0:npjet),jetvz(0:npjet))
#endif
  return
 end subroutine accelerator_update_host_state

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
   att,fve,gr,ks,li,vfield,velext)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,nfieldtype
  logical, intent(in) :: linserted,liniperturb,lairdrag,lflorentz,luppot
  double precision, intent(in) :: pfreq,consistency,findex,yieldstress
  double precision, intent(in) :: att,fve,gr,ks,li,vfield,velext
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
  double precision :: fvet,kst,attt,lit,veltangent
  logical :: straight

  accelerator_eom3_stage=.false.
  if(.not.linserted .or. .not.lairdrag .or. lflorentz .or. luppot)return
  if(nfieldtype/=0 .or. firstpoint/=0 .or. lastpoint/=npjet)return
  nout=lastpoint-firstpoint

#ifdef _OPENACC
!$acc parallel loop gang vector copyin(yxx(0:npjet),yyy(0:npjet), &
!$acc& yzz(0:npjet),yst(0:npjet),yvx(0:npjet),yvy(0:npjet), &
!$acc& yvz(0:npjet),yvl(0:npjet),ycf(0:npjet,1:3), &
!$acc& jetms(0:npjet),jetch(0:npjet),jetfr(0:npjet)) &
!$acc& copyout(fxx(0:nout),fyy(0:nout),fzz(0:nout),fst(0:nout), &
!$acc& fvx(0:nout),fvy(0:nout),fvz(0:nout)) &
!$acc& private(j,dxu,dyu,dzu,dxd,dyd,dzd,lup,ldown,tux,tuy,tuz, &
!$acc& tdx,tdy,tdz,beadvel,v1x,v1y,v1z,v2x,v2y,v2z,l1,l2,dotp, &
!$acc& nbx,nby,nbz,lnb,b,c,t,scale1,scale2,ccx,ccy,ccz,rcx,rcy, &
!$acc& rcz,radius,curvature,factor1,factor2,factor3,factor4,factor5, &
!$acc& fvet,kst,attt,lit,veltangent,straight)
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
    if(ipoint==npjet)then
      if(liniperturb)then
        fyy(j)=-pfreq*yzz(ipoint)
        fzz(j)= pfreq*yyy(ipoint)
        fvy(j)=-(pfreq**2.d0)*yyy(ipoint)
        fvz(j)=-(pfreq**2.d0)*yzz(ipoint)
      endif
      cycle
    endif

    dxu=yxx(ipoint)-yxx(ipoint+1)
    dyu=yyy(ipoint)-yyy(ipoint+1)
    dzu=yzz(ipoint)-yzz(ipoint+1)
    lup=dsqrt(dxu*dxu+dyu*dyu+dzu*dzu)
    tux=dxu/lup
    tuy=dyu/lup
    tuz=dzu/lup
    beadvel=(yvx(ipoint)-yvx(ipoint+1))*tux+ &
     (yvy(ipoint)-yvy(ipoint+1))*tuy+ &
     (yvz(ipoint)-yvz(ipoint+1))*tuz
    fvet=fve/jetms(ipoint)
    factor1=fvet*yvl(ipoint)*(yst(ipoint)/lup)

    fxx(j)=yvx(ipoint)
    fyy(j)=yvy(ipoint)
    fzz(j)=yvz(ipoint)
    fst(j)=yieldstress+consistency*(beadvel/lup)**findex-yst(ipoint)
    fvx(j)=gr+(jetch(ipoint)/jetms(ipoint))*vfield- &
     factor1*tux+ycf(ipoint,1)
    fvy(j)=-factor1*tuy+ycf(ipoint,2)
    fvz(j)=-factor1*tuz+ycf(ipoint,3)

    veltangent=(yvx(ipoint)-velext)*tux+yvy(ipoint)*tuy+ &
     yvz(ipoint)*tuz
    attt=att/jetms(ipoint)
    factor4=attt*(dabs(lup)**0.905d0)*(dabs(veltangent)**1.19d0)
    fvx(j)=fvx(j)-factor4*tux
    fvy(j)=fvy(j)-factor4*tuy
    fvz(j)=fvz(j)-factor4*tuz

    if(ipoint==firstpoint)cycle

    dxd=yxx(ipoint-1)-yxx(ipoint)
    dyd=yyy(ipoint-1)-yyy(ipoint)
    dzd=yzz(ipoint-1)-yzz(ipoint)
    ldown=dsqrt(dxd*dxd+dyd*dyd+dzd*dzd)
    tdx=dxd/ldown
    tdy=dyd/ldown
    tdz=dzd/ldown
    factor2=fvet*yvl(ipoint-1)*(yst(ipoint-1)/ldown)

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

    kst=ks/jetms(ipoint)
    factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(lup))+ &
     (dsqrt(yvl(ipoint-1))/dsqrt(ldown)))**2.d0
    fvx(j)=fvx(j)+factor2*tdx+kst*curvature*factor3*rcx
    fvy(j)=fvy(j)+factor2*tdy+kst*curvature*factor3*rcy
    fvz(j)=fvz(j)+factor2*tdz+kst*curvature*factor3*rcz
    lit=li/jetms(ipoint)
    factor5=factor3*lup*curvature*(veltangent**2.d0)
    fvx(j)=fvx(j)-lit*factor5*rcx
    fvy(j)=fvy(j)-lit*factor5*rcy
    fvz(j)=fvz(j)-lit*factor5*rcz
  enddo
#ifdef _OPENACC
!$acc end parallel loop
#endif

  accelerator_eom3_stage=.true.
 end function accelerator_eom3_stage

end module accelerator_mod
