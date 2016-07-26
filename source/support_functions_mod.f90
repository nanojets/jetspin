 
 module support_functions_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which are called 
!     by different modules (e.g., eom_mod, statistical_mod)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
 use version_mod, only : mystart,myend,sum_world_darr,memyend, &
                   set_chunk,set_mxchunk,or_world_larr
 use utility_mod, only : Pi,modulvec,cross,dot
 use nanojet_mod, only : systype,icrossec,inpjet,npjet,resolution, &
                   velext,linserted,airdragamp,jetms,mxnpjet, &
                   linserting,luppot,kuppot,jetch,jetxx,jetyy,jetzz, &
                   jetpt,jetvl,jetcr,KLor,BLor,lflorentz

 
 implicit none
 
 private
 
 public :: compute_geometry_1d
 public :: compute_geometry_1d_init
 public :: compute_geometry_1d_KV
 public :: beadlength1d
 public :: beadvel1d
 public :: compute_geometry
 public :: compute_geometry_init
 public :: beadlength
 public :: compute_tangetversor
 public :: project_beadveltangetversor
 public :: project_beadacctangetversor
 public :: compute_curvcenter
 public :: compute_curvature
 public :: project_veltangetversor
 public :: compute_stocforce_3d
 public :: compute_crosssec
 public :: compute_length_path
 public :: compute_lorentz_acc
 public :: upwall
 
 contains
 
 subroutine compute_geometry_1d(ipoint,yxx,yst,yvx, &
  beadlendownsub,beadlenupsub,beadvelupsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing geometrical properties 
!     and velocity differences between the i-th bead and the upper
!     and lower bead in the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  
  double precision, intent(out) :: beadlendownsub
  double precision, intent(out) :: beadlenupsub
  double precision, intent(out) :: beadvelupsub
  
  
  beadlenupsub=beadlength1d(ipoint,yxx)
  beadlendownsub=beadlength1d(ipoint-1,yxx)
  beadvelupsub=beadvel1d(ipoint,yvx)
     
  return
  
 end subroutine compute_geometry_1d
 
 subroutine compute_geometry_1d_init(ipoint,yxx,yst,yvx, &
  beadlenupsub,beadvelupsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing geometrical properties 
!     and velocity difference between the first bead and its upper
!     bead in the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  
  double precision, intent(out) :: beadlenupsub
  double precision, intent(out) :: beadvelupsub
  
  
  beadlenupsub=beadlength1d(ipoint,yxx)
  beadvelupsub=beadvel1d(ipoint,yvx)
  
  return
  
 end subroutine compute_geometry_1d_init
 
 subroutine compute_geometry_1d_KV(ipoint,yxx,yst,yvx,yax, &
  beadlenupsub,beadvelupsub,beadaccupsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing geometrical properties 
!     velocity difference and acceleration difference between the first 
!     bead and its upper bead in the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yax
  
  double precision, intent(out) :: beadlenupsub
  double precision, intent(out) :: beadvelupsub
  double precision, intent(out) :: beadaccupsub
  
  
  beadlenupsub=beadlength1d(ipoint,yxx)
  beadvelupsub=beadvel1d(ipoint,yvx)
  beadaccupsub=beadacc1d(ipoint,yax)
  
  return
  
 end subroutine compute_geometry_1d_KV
 
 function beadlength1d(ipoint,yxx)
 
!***********************************************************************
!     
!     JETSPIN function for computing the mutual distance 
!     between the i-th bead and its upper
!     bead in the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  
  double precision :: beadlength1d
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        beadlength1d = yxx(ipoint)-yxx(npjet)
      else
        beadlength1d = yxx(ipoint)-yxx(ipoint+1)
      endif
    else
      beadlength1d = 0.d0
    endif
  else
    beadlength1d = yxx(ipoint)-yxx(ipoint+1)
  endif
     
  return
  
 end function beadlength1d
 
 function beadvel1d(ipoint,yvx)
 
!***********************************************************************
!     
!     JETSPIN function for computing the velocity difference
!     between the i-th bead and its upper
!     bead in the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  
  double precision :: beadvel1d
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        beadvel1d = yvx(ipoint)-yvx(npjet)
      else
        beadvel1d = yvx(ipoint)-yvx(ipoint+1)
      endif
    else
      beadvel1d = 0.d0
    endif
  else
    beadvel1d = yvx(ipoint)-yvx(ipoint+1)
  endif
     
  return
  
 end function beadvel1d
 
 function beadacc1d(ipoint,yax)
 
!***********************************************************************
!     
!     JETSPIN function for computing the acceleration difference
!     between the i-th bead and its upper
!     bead in the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yax
  
  double precision :: beadacc1d
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        beadacc1d = yax(ipoint)-yax(npjet)
      else
        beadacc1d = yax(ipoint)-yax(ipoint+1)
      endif
    else
      beadacc1d = 0.d0
    endif
  else
    beadacc1d = yax(ipoint)-yax(ipoint+1)
  endif
     
  return
  
 end function beadacc1d
 
 subroutine compute_geometry(ipoint,yxx,yyy,yzz, &
  beadlendownsub,beadlenupsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing geometrical properties 
!     and velocity differences between the i-th bead and the upper
!     and lower bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  
  double precision, intent(out) :: beadlendownsub
  double precision, intent(out) :: beadlenupsub
  
  beadlendownsub=beadlength(ipoint-1,yxx,yyy,yzz)
  beadlenupsub=beadlength(ipoint,yxx,yyy,yzz)
     
  return
  
 end subroutine compute_geometry
 
 subroutine compute_geometry_init(ipoint,yxx,yyy,yzz, &
  beadlenupsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing geometrical properties 
!     and velocity differences between the first bead and its upper
!     bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  
  double precision, intent(out) :: beadlenupsub
  
  
  beadlenupsub=beadlength(ipoint,yxx,yyy,yzz)
     
  return
  
 end subroutine compute_geometry_init
 
 function beadlength(ipoint,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN function for computing the mutual distance 
!     between the i-th bead and its upper
!     bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  
  double precision :: beadlength
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        beadlength = dsqrt((yxx(ipoint)-yxx(npjet))**2.d0+ &
         (yyy(ipoint)-yyy(npjet))**2.d0+(yzz(ipoint)-yzz(npjet))**2.d0)
      else
        beadlength = dsqrt((yxx(ipoint)-yxx(ipoint+1))**2.d0+ &
         (yyy(ipoint)-yyy(ipoint+1))**2.d0+(yzz(ipoint)- &
         yzz(ipoint+1))**2.d0)
      endif
    else
      beadlength = 0.d0
    endif
  else
    beadlength = dsqrt((yxx(ipoint)-yxx(ipoint+1))**2.d0+ &
     (yyy(ipoint)-yyy(ipoint+1))**2.d0+(yzz(ipoint)- &
     yzz(ipoint+1))**2.d0)
  endif
     
  return
  
 end function beadlength
 
 subroutine compute_tangetversor(ipoint,yxx,yyy,yzz,utang,beadlenupsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the unit vector pointing
!     the i-th bead and starting by its upper
!     bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) ::  ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  double precision, intent(out), dimension(3) :: utang
  double precision, intent(in) :: beadlenupsub
  
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        utang(1)=(yxx(ipoint)-yxx(npjet))/beadlenupsub
        utang(2)=(yyy(ipoint)-yyy(npjet))/beadlenupsub
        utang(3)=(yzz(ipoint)-yzz(npjet))/beadlenupsub
      else
        utang(1)=(yxx(ipoint)-yxx(ipoint+1))/beadlenupsub
        utang(2)=(yyy(ipoint)-yyy(ipoint+1))/beadlenupsub
        utang(3)=(yzz(ipoint)-yzz(ipoint+1))/beadlenupsub
      endif
    else
      utang(1)=0.d0
      utang(2)=0.d0
      utang(3)=0.d0
    endif
  else
    utang(1)=(yxx(ipoint)-yxx(ipoint+1))/beadlenupsub
    utang(2)=(yyy(ipoint)-yyy(ipoint+1))/beadlenupsub
    utang(3)=(yzz(ipoint)-yzz(ipoint+1))/beadlenupsub
  endif
  
  return
  
 end subroutine compute_tangetversor
 
 subroutine project_beadveltangetversor(ipoint,yvx,yvy,yvz, &
  beadveltang,utang)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the module velocity difference
!     between the i-th bead and its upper
!     bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) ::  ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, intent(out) :: beadveltang
  double precision, intent(in), dimension(3) :: utang
  
  double precision, dimension(3) :: vtang,uvtang
  double precision :: normvtang,dotp
  
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        vtang(1)=yvx(ipoint)-yvx(npjet)
        vtang(2)=yvy(ipoint)-yvy(npjet)
        vtang(3)=yvz(ipoint)-yvz(npjet)
      else
        vtang(1)=yvx(ipoint)-yvx(ipoint+1)
        vtang(2)=yvy(ipoint)-yvy(ipoint+1)
        vtang(3)=yvz(ipoint)-yvz(ipoint+1)
      endif
    else
      vtang(1:3)=0.d0
    endif
  else
    vtang(1)=yvx(ipoint)-yvx(ipoint+1)
    vtang(2)=yvy(ipoint)-yvy(ipoint+1)
    vtang(3)=yvz(ipoint)-yvz(ipoint+1)
  endif
  
  normvtang=dsqrt(vtang(1)**2.d0+vtang(2)**2.d0+vtang(3)**2.d0)
  if(normvtang>0.d0)then
    uvtang(:)=vtang(:)/normvtang
    dotp=dot_product(uvtang,utang)
    beadveltang=normvtang*dotp
  else
    beadveltang=0.d0
  endif
  
  return
  
 end subroutine project_beadveltangetversor
 
 subroutine project_beadacctangetversor(ipoint,yvx,yvy,yvz, &
  beadveltang,utang)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the module acceleration 
!     difference between the i-th bead and its upper
!     bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) ::  ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, intent(out) :: beadveltang
  double precision, intent(in), dimension(3) :: utang
  
  double precision, dimension(3) :: vtang,uvtang
  double precision :: normvtang,dotp
  
  
  if(ipoint==npjet-2)then
    if(ipoint>=0)then
      if(.not.linserted)then
        vtang(1)=yvx(ipoint)-yvx(npjet)
        vtang(2)=yvy(ipoint)-yvy(npjet)
        vtang(3)=yvz(ipoint)-yvz(npjet)
      else
        vtang(1)=yvx(ipoint)-yvx(ipoint+1)
        vtang(2)=yvy(ipoint)-yvy(ipoint+1)
        vtang(3)=yvz(ipoint)-yvz(ipoint+1)
      endif
    else
      vtang(1:3)=0.d0
    endif
  else
    vtang(1)=yvx(ipoint)-yvx(ipoint+1)
    vtang(2)=yvy(ipoint)-yvy(ipoint+1)
    vtang(3)=yvz(ipoint)-yvz(ipoint+1)
  endif
  
  normvtang=dsqrt(vtang(1)**2.d0+vtang(2)**2.d0+vtang(3)**2.d0)
  if(normvtang>0.d0)then
    uvtang(:)=vtang(:)/normvtang
    dotp=dot_product(uvtang,utang)
    beadveltang=normvtang*dotp
  else
    beadveltang=0.d0
  endif
  
  return
  
 end subroutine project_beadacctangetversor
 
 subroutine compute_curvcenter(ipoint,yxx,yyy,yzz,curvcentersub, &
  lstraightsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the center of curvature
!     between three beads: the i-th bead, the upper bead
!     and the lower bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  double precision, dimension (3), intent(inout) :: curvcentersub
  logical, intent(out) :: lstraightsub
  

  double precision, dimension (3) :: p1,p2,p3
  double precision, dimension (3) :: v1,v2
  double precision, dimension (3) :: v1n,v2n
  double precision, dimension (3) :: v1nb,v2nb
  double precision, dimension (2) :: p3_2d
  double precision :: l1,l2,l2nb,dotp,a,b,c,t,scale1,scale2
  
  integer :: i
  
  
! take the three beads coordinates
  p1(1)=yxx(ipoint-1)
  p1(2)=yyy(ipoint-1)
  p1(3)=yzz(ipoint-1)
  p2(1)=yxx(ipoint)
  p2(2)=yyy(ipoint)
  p2(3)=yzz(ipoint)
  if(ipoint+1>npjet)then
    p3(:)=0.d0
  else
    p3(1)=yxx(ipoint+1)
    p3(2)=yyy(ipoint+1)
    p3(3)=yzz(ipoint+1)
  endif
  
! define two vector between the three beads
  v1(1:3)=(p3(1:3)-p2(1:3))
  v2(1:3)=(p1(1:3)-p2(1:3))
  
! compute the respective unit vector
  l1=modulvec(v1)
  l2=modulvec(v2)
  
  if(l1==0.d0 .or. l2==0.d0)then
    lstraightsub=.true.
    curvcentersub(:) = 0.d0
    return
  endif
  
  v1n=v1/l1
  v2n=v2/l2

  v1nb=v1n

! use the Gram–Schmidt process to find an orthogonal basis in 2-D space
  dotp=dot_product(v2n(:),v1nb(:))
  v2nb(:)=v2n(:)-dotp*v1nb(:)

  l2nb=modulvec(v2nb)
  
  if(l2nb==0.d0)then
    lstraightsub=.true.
    curvcentersub(:) = 0.d0
    return
  endif
  
  v2nb(:)=v2nb(:)/l2nb
  

  p3_2d(1) = dot_product(v2,v1nb)
  p3_2d(2) = dot_product(v2,v2nb)
  
  if(p3_2d(2)==0.d0)then
    lstraightsub=.true.
    curvcentersub(:) = 0.d0
    return
  endif
  

! compute the curvature center in this 2-D space
  a = l1
  b = p3_2d(1) 
  c = p3_2d(2)
  t = 0.5d0*(a-b)/c
  scale1 = b/2.d0 + c*t
  scale2 = c/2.d0 - b*t
  
! projection of the curvature center into the 3-D space
  curvcentersub(:) = 0.d0
  do i=1,3
    curvcentersub(i) = p2(i) + scale1*v1nb(i) + scale2*v2nb(i)
  enddo
  
  
! check if the curve center is not infinity 
  if(isnan(curvcentersub(1)).or.isnan(curvcentersub(2)).or. &
   isnan(curvcentersub(3)))then
!   if yes the segment is straight
    lstraightsub=.true.
    curvcentersub(:) = 0.d0
  else
    lstraightsub=.false.
  endif
  
  
  return
  
 end subroutine compute_curvcenter
 
 subroutine compute_curvature(ipoint,yxx,yyy,yzz,curvaturesub, &
  vcurvaturesub,curvcentersub,lstraightsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the curvature of the jet
!     in the segment represented by three beads: the i-th bead,
!     the upper bead and the lower bead in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) ::  ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  double precision, intent(out) :: curvaturesub
  double precision, dimension (3), intent(out) :: vcurvaturesub
  double precision, dimension (3), intent(in) :: curvcentersub
  logical, intent(in) :: lstraightsub
  
  double precision :: radius
  double precision, dimension(3) :: vtemp
  
  if(lstraightsub)then
    vcurvaturesub(1:3)=0.d0
    curvaturesub=0.d0
    return
  endif
  
! compute the radius of curvature
  vtemp(1)=curvcentersub(1)-yxx(ipoint)
  vtemp(2)=curvcentersub(2)-yyy(ipoint)
  vtemp(3)=curvcentersub(3)-yzz(ipoint)
  radius=modulvec(vtemp)
  
! compute the curvature
  curvaturesub=1.d0/radius
  
! compute the unit vector starting  from the i-th nead and pointing the
! curvatre center
  vcurvaturesub(1:3)=vtemp(1:3)/radius
  
  return
  
 end subroutine compute_curvature
 
 subroutine project_veltangetversor(ipoint,yvx,yvy,yvz,veltang,utang)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the module velocity difference
!     projection between the i-th bead and its upper
!     bead along the unit vector jointing this two beads 
!     in the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) ::  ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  
  double precision, intent(out) :: veltang
  double precision, intent(in), dimension(3) :: utang
  
  double precision, dimension(3) :: vtang,uvtang
  double precision :: normvtang,dotp

  
  vtang(1)=yvx(ipoint)-velext
  vtang(2)=yvy(ipoint)
  vtang(3)=yvz(ipoint)
  normvtang=dsqrt(vtang(1)**2.d0+vtang(2)**2.d0+vtang(3)**2.d0)
  if(normvtang==0.d0)then
    veltang=0.d0
  else
    uvtang(:)=vtang(:)/normvtang
    dotp=dot_product(uvtang,utang)
    veltang=normvtang*dotp
  endif
  
  return
  
 end subroutine project_veltangetversor
 
 subroutine compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the stochastic force in the
!     three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) ::  ipoint
  double precision, intent(out) ::  fstocvx
  double precision, intent(out) ::  fstocvy
  double precision, intent(out) ::  fstocvz
  
  double precision :: factorsub,dadt
  
  dadt=airdragamp(1)/jetms(ipoint)
  factorsub=dsqrt(2.d0*dadt)
  
  fstocvx = factorsub
  fstocvy = factorsub
  fstocvz = factorsub  
  
  return
  
 end subroutine compute_stocforce_3d
 
 subroutine compute_crosssec(jetxxs,jetyys,jetzzs,jetvls,jetcrs)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the cross section radius of
!     the nanofiber segment represented by the i-th bead
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), &
   intent(in) ::  jetxxs,jetyys,jetzzs,jetvls
  double precision, allocatable, dimension (:), &
   intent(inout) ::  jetcrs
  
  integer :: ipoint
  double precision :: tempmod0
  
  call set_chunk(inpjet,npjet)
  call set_mxchunk(mxnpjet)
  
  jetcrs(:)=0.d0
    
  if(linserting)then
    if(.not.linserted)then
      select case(systype)
        case(1)
          do ipoint=mystart,myend
            if(ipoint<npjet-2)then
              tempmod0=jetxxs(ipoint)-jetxxs(ipoint+1)
              jetcrs(ipoint) = dsqrt(jetvls(ipoint)/(tempmod0*Pi))
            else
              jetcrs(ipoint) = icrossec
            endif
          enddo
        case default
          do ipoint=mystart,myend
            if(ipoint<npjet-2)then
              tempmod0=dsqrt((jetxxs(ipoint)-jetxxs(ipoint+1))**2.d0+ &
               (jetyys(ipoint)-jetyys(ipoint+1))**2.d0+ &
               (jetzzs(ipoint)-jetzzs(ipoint+1))**2.d0)
              jetcrs(ipoint) = dsqrt(jetvls(ipoint)/(tempmod0*Pi))
            else
              jetcrs(ipoint) = icrossec
            endif
          enddo
      end select
    else
      select case(systype)
        case(1)
          do ipoint=mystart,myend
            if(ipoint<npjet-1)then
              tempmod0=jetxxs(ipoint)-jetxxs(ipoint+1)
              jetcrs(ipoint) = dsqrt(jetvls(ipoint)/(tempmod0*Pi))
            else
              jetcrs(ipoint) = icrossec
            endif
          enddo
        case default
          do ipoint=mystart,myend
            if(ipoint<npjet-1)then
              tempmod0=dsqrt((jetxxs(ipoint)-jetxxs(ipoint+1))**2.d0+ &
               (jetyys(ipoint)-jetyys(ipoint+1))**2.d0+ &
               (jetzzs(ipoint)-jetzzs(ipoint+1))**2.d0)
              jetcrs(ipoint) = dsqrt(jetvls(ipoint)/(tempmod0*Pi))
            else
              jetcrs(ipoint) = icrossec
            endif
          enddo
      end select
    endif
  else
    select case(systype)
      case(1)
        do ipoint=mystart,myend
          if(ipoint<npjet)then
            tempmod0=jetxxs(ipoint)-jetxxs(ipoint+1)
            jetcrs(ipoint) = dsqrt(jetvls(ipoint)/(tempmod0*Pi))
          else
            jetcrs(ipoint) = icrossec
          endif
        enddo
      case default
        do ipoint=mystart,myend
          if(ipoint<npjet)then
            tempmod0=dsqrt((jetxxs(ipoint)-jetxxs(ipoint+1))**2.d0+ &
             (jetyys(ipoint)-jetyys(ipoint+1))**2.d0+ &
             (jetzzs(ipoint)-jetzzs(ipoint+1))**2.d0)
            jetcrs(ipoint) = dsqrt(jetvls(ipoint)/(tempmod0*Pi))
          else
            jetcrs(ipoint) = icrossec
          endif
        enddo
    end select
  endif
  
  call sum_world_darr(jetcrs,npjet+1)
  
  return
  
 end subroutine compute_crosssec
 
 subroutine compute_length_path(yxx,yyy,yzz, &
  ypt,lengthpathsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the path length of the nanofiber
!     between the nozzle and the collector and the curve parameter ypt
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension(:), intent(in) :: yxx
  double precision, allocatable, dimension(:), intent(in) :: yyy
  double precision, allocatable, dimension(:), intent(in) :: yzz
  
  double precision,allocatable,dimension(:),intent(inout) :: ypt
  double precision, intent(out) :: lengthpathsub
  
  integer :: ipoint
  double precision :: tempmod0,tempmod,tempmod1(1)
  
  call set_chunk(inpjet,npjet)
  call set_mxchunk(mxnpjet)
  
  tempmod=0.d0
  ypt(:)=0.d0
  select case(systype)
    case(1)
    do ipoint=mystart+1,memyend
        tempmod0 = yxx(ipoint-1)-yxx(ipoint)
        tempmod = tempmod+tempmod0
        ypt(ipoint) = tempmod
      enddo
    case default
      do ipoint=mystart+1,memyend
        tempmod0 = dsqrt((yxx(ipoint-1)-yxx(ipoint))**2.d0+ &
         (yyy(ipoint-1)-yyy(ipoint))**2.d0+ &
         (yzz(ipoint-1)-yzz(ipoint))**2.d0)
        tempmod = tempmod+tempmod0
        ypt(ipoint) = tempmod
      enddo
  end select
  
  tempmod1(1)=tempmod
  call sum_world_darr(tempmod1,1)
  lengthpathsub=tempmod1(1)
  
  ypt(myend+1:npjet)=tempmod
  
  call sum_world_darr(ypt,npjet+1)
  
  ypt(inpjet:npjet) = ypt(inpjet:npjet)/lengthpathsub
  
  return
  
 end subroutine compute_length_path
 
 
 
 subroutine compute_lorentz_acc(k,vx,vy,vz,aLorx,aLory,aLorz)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the acceleration given by lorentz
!     force
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella, F. Cipolletta
!     last modification April 2016
!     
!***********************************************************************

  implicit none
  integer, intent(in) :: k
  double precision, allocatable, dimension(:), intent(in) :: vx,vy,vz
  double precision, intent(out) :: aLorx,aLory,aLorz                                                             
  
  double precision :: KLorsub
  
  if(.not. lflorentz)then
    aLorx=0.d0
    aLory=0.d0
    aLorz=0.d0
    return
  endif
  
  KLorsub=KLor*jetch(k)/(dsqrt(jetms(k)))
  
  aLorx = KLorsub *(( vy(k)*BLor(3) )-( vz(k)*BLor(2) ))
  aLory = KLorsub *(( vz(k)*BLor(1) )-( vx(k)*BLor(3) ))
  aLorz = KLorsub *(( vx(k)*BLor(2) )-( vy(k)*BLor(1) ))                                       
    
  return
    
 end subroutine compute_lorentz_acc
 
 function upwall(ipoint,yxx)
 
!***********************************************************************
!     
!     JETSPIN function for computing the acceleration given by the 
!     uppest soft wall located behind the nozzle at x equal zero
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  
  double precision :: upwall
  
  if(luppot)then
    if(yxx(ipoint)>=0.d0)then
      upwall=0.d0
    else
      upwall=-1.d0*kuppot*yxx(ipoint)/jetms(ipoint)
    endif
  else
    upwall=0.d0
  endif
  
  return
  
 end function upwall
 
 end module support_functions_mod


