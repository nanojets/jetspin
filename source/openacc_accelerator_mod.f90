module accelerator_mod

 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.true.

 public :: accelerator_prepare
 public :: accelerator_eom3_stage

contains

 subroutine accelerator_prepare()
  implicit none
!$acc init
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
!$acc end parallel loop

  accelerator_eom3_stage=.true.
 end function accelerator_eom3_stage

end module accelerator_mod
