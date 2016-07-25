 module electric_field_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which manage and compute 
!     the external electric field 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
 use utility_mod,           only : Pi
 use nanojet_mod,           only : jetms,jetch,v,fieldfreq,taoelectr
 
 implicit none
 
 private
 
 integer, public, save :: nfieldtype=0
 
 public :: driver_electric_field
 public :: actual_form_electric_field
 
 contains
 
 subroutine driver_electric_field(ipoint,timesub,yxx,yyy,yzz,vout)

!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute the actual external electric field 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  double precision, intent(in) :: timesub
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  double precision, intent(out), dimension(3) :: vout
  
  double precision :: swave
  
  select case(nfieldtype)
  case(1)
    call electric_field_rectangular_wave(timesub,swave)
    vout(1)=swave*(jetch(ipoint)/jetms(ipoint))*V
    vout(2:3)=0.d0
  case(2)
    call electric_field_rectangular_wave_tao(timesub,swave)
    vout(1)=swave*(jetch(ipoint)/jetms(ipoint))*V
    vout(2:3)=0.d0
  case default
    vout(1)=(jetch(ipoint)/jetms(ipoint))*V
    vout(2:3)=0.d0
  end select
  
  return
  
 end subroutine driver_electric_field
 
 subroutine actual_form_electric_field(timesub,vout)

!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute the actual external electric field 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: timesub
  double precision, intent(out), dimension(3) :: vout
  
  double precision :: swave
  
  select case(nfieldtype)
  case(1)
    call electric_field_rectangular_wave(timesub,swave)
    vout(1)=swave
    vout(2:3)=0.d0
  case(2)
    call electric_field_rectangular_wave_tao(timesub,swave)
    vout(1)=swave
    vout(2:3)=0.d0
  case default
    vout(1)=1.d0
    vout(2:3)=0.d0
  end select
  
  return
  
 end subroutine actual_form_electric_field
 
 subroutine electric_field_rectangular_wave(timesub,swave)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the actual external 
!     electric field which is described by a rectangular wave
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: timesub
  double precision, intent(out) :: swave
  
  double precision :: dtemp
  
  dtemp=dsin(2.d0*pi*fieldfreq*timesub)
  
  if(dtemp>0.d0)then
    swave=1.d0
  elseif(dtemp==0.d0)then
    swave=0.5d0
  else
    swave=0.d0
  endif
  
  return
 
 end subroutine electric_field_rectangular_wave
 
 subroutine electric_field_rectangular_wave_tao(timesub,swave)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the actual external 
!     electric field which is described by a rectangular wave in a RC
!     circuit with a given relaxation time
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: timesub
  double precision, intent(out) :: swave
  
  double precision :: dtemp1,dtemp2,dtemp3,dtemp4
  
  dtemp1=2.d0*pi*fieldfreq*timesub
  dtemp2=dsin(dtemp1)
  dtemp3=dmod(dtemp1,2.d0*pi)
  dtemp4=dtemp3/(2.d0*pi*fieldfreq)
  
  if(dtemp2>0.d0)then
    dtemp4=dtemp3/(2.d0*pi*fieldfreq)
    swave=1.d0-dexp(-dtemp4/taoelectr)
  elseif(dtemp2==0.d0)then
    swave=0.d0
  else
    dtemp4=(dtemp3-pi)/(2.d0*pi*fieldfreq)
    swave=dexp(-dtemp4/taoelectr)
  endif
  
  return
 
 end subroutine electric_field_rectangular_wave_tao
 
 end module electric_field_mod
 
