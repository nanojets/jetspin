
 module integrator_mod
 
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
 use utility_mod,       only : gauss,wiener_process1,wiener_process2
 use nanojet_mod,       only : doallocate,mxnpjet,npjet,inpjet,systype,&
                         jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,jetst, &
                         jetms,jetch,jetvl,compute_posnoinserted, &
                         jetpt,lKVfluid
 use dynamic_refinement_mod, only : driver_dynamic_refinement
 use coulomb_force_mod, only : smooth_charge,restore_charge, &
                         coulforce,compute_coulomelec_driver
 use driver_eom_mod,    only : xpsys,xpsys_pos,xpsys_stress, &
                         xpsys_KV_pos_v,xpsys_KV_st

 implicit none

 private
 
 integer, public, save :: integrator
 logical, public, save :: lintegrator=.false.
 
 double precision, public, save :: initime = 0.d0
 double precision, public, save :: endtime = 5.d0
 logical, public, save :: lendtime
 
 public :: driver_integrator

 contains
  
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
  logical :: ltestinst
  
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
  
  ltestinst=.false.
  do i=inpjet,npjet
    if(isnan(dcos(jetxx(i))))ltestinst=.true.
    if(isnan(dcos(jetyy(i))))ltestinst=.true.
    if(isnan(dcos(jetzz(i))))ltestinst=.true.
    if(isnan(dcos(jetst(i))))ltestinst=.true.
    if(isnan(dcos(jetvx(i))))ltestinst=.true.
    if(isnan(dcos(jetvy(i))))ltestinst=.true.
    if(isnan(dcos(jetvz(i))))ltestinst=.true.
  enddo
  
  if(ltestinst)then
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
      do ipoint=mystart,myend
        call xpsys(ipoint,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
          jetvl,coulforce,fxx(j),fyy(j),fzz(j),fst(j), &
          fvx(j),fvy(j),fvz(j),timesub,k)
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
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),timesub+h,k)
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
      yxx(:)=0.d0
      yyy(:)=0.d0
      yzz(:)=0.d0
      yst(:)=0.d0
      yvx(:)=0.d0
      yvy(:)=0.d0
      yvz(:)=0.d0
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
	    yxx(ipoint) = jetxx(ipoint) + h*f1xx(j)
	    yyy(ipoint) = jetyy(ipoint) + h*f1yy(j)
	    yzz(ipoint) = jetzz(ipoint) + h*f1zz(j)
	    yst(ipoint) = jetst(ipoint) + h*f1st(j)
	    yvx(ipoint) = jetvx(ipoint) + h*f1vx(j)
	    yvy(ipoint) = jetvy(ipoint) + h*f1vy(j)
	    yvz(ipoint) = jetvz(ipoint) + h*f1vz(j)
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
      
      call smooth_charge(yxx,yyy,yzz)
      call compute_posnoinserted(yxx,yyy,yzz)
      call compute_coulomelec_driver(k,timesub,coulforce,jetvl,yxx, &
       yyy,yzz)
      j=0
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),timesub+h,k)
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
	    yvx(ipoint) = jetvx(ipoint) + (h/2.d0)*(f1vx(j)+f2vx(j))
	    yvy(ipoint) = jetvy(ipoint) + (h/2.d0)*(f1vy(j)+f2vy(j))
	    yvz(ipoint) = jetvz(ipoint) + (h/2.d0)*(f1vz(j)+f2vz(j))
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
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),timesub+h/2.d0,k)
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
         jetvl,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
         f3vx(j),f3vy(j),f3vz(j),timesub+h/2.d0,k)
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
         jetvl,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
         f4vx(j),f4vy(j),f4vz(j),timesub+h,k)
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
        yxx(ipoint) = jetxx(ipoint) + 0.5d0*h*f1xx(j)
        yyy(ipoint) = jetyy(ipoint) + 0.5d0*h*f1yy(j)
        yzz(ipoint) = jetzz(ipoint) + 0.5d0*h*f1zz(j)
        yst(ipoint) = jetst(ipoint) + 0.5d0*h*f1st(j)
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f1vx(j)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f1vy(j)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f1vz(j)
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
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f2xx(j),f2yy(j),f2zz(j),f2st(j), &
         f2vx(j),f2vy(j),f2vz(j),timesub+h/2.d0,k)
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
        yvx(ipoint) = jetvx(ipoint) + 0.5d0*h*f2vx(j)
        yvy(ipoint) = jetvy(ipoint) + 0.5d0*h*f2vy(j)
        yvz(ipoint) = jetvz(ipoint) + 0.5d0*h*f2vz(j)
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
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f3xx(j),f3yy(j),f3zz(j),f3st(j), &
         f3vx(j),f3vy(j),f3vz(j),timesub+h/2.d0,k)
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
        yvx(ipoint) = jetvx(ipoint) + h*f3vx(j)
        yvy(ipoint) = jetvy(ipoint) + h*f3vy(j)
        yvz(ipoint) = jetvz(ipoint) + h*f3vz(j)
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
      do ipoint=mystart,myend
        call xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz, &
         jetvl,coulforce,f4xx(j),f4yy(j),f4zz(j),f4st(j), &
         f4vx(j),f4vy(j),f4vz(j),timesub+h,k)
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
        yvx(ipoint) = jetvx(ipoint) + (h/6.d0)*(f1vx(j)+ &
         2.d0*(f2vx(j)+f3vx(j))+f4vx(j))
        yvy(ipoint) = jetvy(ipoint) + (h/6.d0)*(f1vy(j)+ &
         2.d0*(f2vy(j)+f3vy(j))+f4vy(j))
        yvz(ipoint) = jetvz(ipoint) + (h/6.d0)*(f1vz(j)+ &
         2.d0*(f2vz(j)+f3vz(j))+f4vz(j))
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
    end select
    lfirstsub=.false.
  endif
  
  dsqrh=dsqrt(dabs(h))
  tsqh=dsqrh**3.d0
  prefactor1=0.5d0/dsqrh
  
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
        u1=gauss()
        u2=gauss()
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
        u1=gauss()
        u2=gauss()
        ww(1)=(dsqrh*u1)
        zz(1)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
        u1=gauss()
        u2=gauss()
        ww(2)=(dsqrh*u1)
        zz(2)=0.5d0*tsqh*(u1+1.d0/(dsqrt(3.d0))*u2)
        u1=gauss()
        u2=gauss()
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
 
 end module integrator_mod

