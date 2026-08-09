module accelerator_mod

 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.true.
 logical, save :: accelerator_persistent=.false.
 logical, save :: accelerator_statistics_mapped=.false.
 double precision, save :: statistics_step_max=-huge(0.d0)
 integer, save :: statistics_step_index=-1
 integer, save :: accelerator_last_host_sync_step=-huge(0)

 public :: accelerator_prepare
 public :: accelerator_eom3_stage
 public :: accelerator_set_persistent
 public :: accelerator_is_persistent
 public :: accelerator_update_host_state
 public :: accelerator_update_host_point
 public :: accelerator_host_state_is_current
 public :: accelerator_store_statistics
 public :: accelerator_update_host_statistics
 public :: accelerator_update_device_statistics
 public :: accelerator_rk4_final_statistics
 public :: accelerator_euler_final_statistics
 public :: accelerator_rk2_final_statistics
 public :: accelerator_platen_predict
 public :: accelerator_platen_velocity
 public :: accelerator_platen_positions
 public :: accelerator_platen_stress_statistics

contains

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

  if(.not.accelerator_persistent)return
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
  return
 end subroutine accelerator_set_persistent

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
  if(.not.accelerator_persistent)return
#ifdef _OPENACC
!$acc update self(jetxx(ipoint),jetyy(ipoint),jetzz(ipoint), &
!$acc& jetst(ipoint),jetvx(ipoint),jetvy(ipoint),jetvz(ipoint))
#endif
  return
 end subroutine accelerator_update_host_point

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
   att,fve,gr,ks,li,vfield,velext,stochastic_model,noisefric)
  implicit none
  integer, intent(in) :: firstpoint,lastpoint,npjet,nfieldtype
  logical, intent(in) :: linserted,liniperturb,lairdrag,lflorentz,luppot
  logical, intent(in) :: stochastic_model
  double precision, intent(in) :: pfreq,consistency,findex,yieldstress
  double precision, intent(in) :: att,fve,gr,ks,li,vfield,velext,noisefric
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
    if(stochastic_model .and. yst(ipoint)<=0.d0)then
      factor1=0.d0
    else
      factor1=fvet*yvl(ipoint)*(yst(ipoint)/lup)
    endif

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
    if(stochastic_model)then
      fvx(j)=fvx(j)-noisefric*yvx(ipoint)
      fvy(j)=fvy(j)-noisefric*yvy(ipoint)
      fvz(j)=fvz(j)-noisefric*yvz(ipoint)
    endif

    if(ipoint==firstpoint)cycle

    dxd=yxx(ipoint-1)-yxx(ipoint)
    dyd=yyy(ipoint-1)-yyy(ipoint)
    dzd=yzz(ipoint-1)-yzz(ipoint)
    ldown=dsqrt(dxd*dxd+dyd*dyd+dzd*dzd)
    tdx=dxd/ldown
    tdy=dyd/ldown
    tdz=dzd/ldown
    if(stochastic_model .and. yst(ipoint-1)<=0.d0)then
      factor2=0.d0
    else
      factor2=fvet*yvl(ipoint-1)*(yst(ipoint-1)/ldown)
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
