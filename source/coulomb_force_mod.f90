 module coulomb_force_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which manage and compute 
!     the Coulomb force 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 use version_mod,           only : idrank,mxrank,sum_world_darr
 use utility_mod,           only : Pi,modulvec,cross,dot,sig
 use nanojet_mod,           only : jetch,ncutoff,inpjet,npjet,regq,q,&
                             linserted,dresolution,jetms,systype,fcut, &
                             thresolution,lmirror,h,jetfr
 use support_functions_mod, only : beadlength1d,beadlength
 
 implicit none
 
 private
 
 double precision, save :: smoothedcharge
 
 integer, save :: ncoulforce=0
 double precision, allocatable, public, save :: coulforce(:,:)
 
 public :: smooth_charge
 public :: restore_charge
 public :: coulomb_1d
 public :: compute_coulomelec3d
 public :: compute_coulomelec
 
 contains
 
 subroutine allocate_array_coulforce(imiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the array coulforce 
!     which is used to store the Coulomb forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(ncoulforce/=0)then
    if(imiomax>ncoulforce)then
      deallocate(coulforce)
      ncoulforce=imiomax+100
      allocate(coulforce(0:ncoulforce,3))
    endif
  else
    ncoulforce=imiomax+100
    allocate(coulforce(0:ncoulforce,3))
  endif
  
  return
  
 end subroutine allocate_array_coulforce
 
 subroutine smooth_charge(yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for smoothing the charge of the last added
!     bead if the mutual distance between this bead and the nozzle is
!     less than two times the resolution length step (dresolution)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) ::  yxx
  double precision, allocatable, dimension (:), intent(inout), &
    optional ::  yyy
  double precision, allocatable, dimension (:), intent(inout), &
   optional ::  yzz
  
  double precision :: dtemp
 
  if(.not.linserted)then
    select case(systype)
    case(1:2)
      dtemp=fcut(beadlength1d(npjet-2,yxx),thresolution,dresolution)
    case default
      dtemp=fcut(beadlength(npjet-2,yxx,yyy,yzz),thresolution, &
       dresolution)
    end select
    smoothedcharge=jetch(npjet-1)
    jetch(npjet-1)=dtemp*smoothedcharge
  endif
  
 return
 
 end subroutine smooth_charge
 
 subroutine restore_charge()
 
!***********************************************************************
!     
!     JETSPIN subroutine for restoring the smoothed charge of the last 
!     added bead
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
 
  if(.not.linserted)then
    jetch(npjet-1)=smoothedcharge
  endif
  
 return
 
 end subroutine restore_charge
 
 function coulomb_1d(ipoint,yxx)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb force in the
!     one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
 
  integer, intent(in) :: ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  
  double precision :: coulomb_1d
  
  double precision :: totelect,Qt
  integer :: i,ista,ienf
  
  ista=inpjet
  ienf=npjet
  
! compute the Coulomb force
  Qt=jetch(ipoint)*Q/jetms(ipoint)
  totelect=0.d0
  do i=ista,ipoint-1
    totelect=totelect-1.d0*(jetch(i)*Qt)/((dabs(yxx(i)-yxx(ipoint))+ &
     regq)**2.d0)
  enddo
  do i=ipoint+1,ienf
    totelect=totelect+1.d0*(jetch(i)*Qt)/((dabs(yxx(i)-yxx(ipoint))+ &
     regq)**2.d0)
  enddo
  coulomb_1d=totelect
  
  return
  
 end function coulomb_1d
 
 subroutine compute_coulomelec3d(ipoint,yxx,yyy,yzz,coulomelec)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb force in the
!     three dimensional model
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
  double precision, dimension (3), intent(out) ::  coulomelec
  
  double precision, parameter :: onethird=1.d0/(dsqrt(3.d0))
  
  integer :: i,ista,ienf
  double precision :: norm,Qt
  double precision, dimension(3) :: versor,utang
  
  ista=inpjet
  ienf=npjet
  
! compute the Coulomb force
  Qt=jetch(ipoint)*Q/jetms(ipoint)
  coulomelec(:)=0.d0
  do i=ista,ipoint-1
    utang(1)=(yxx(ipoint)-yxx(i))+onethird*regq
    utang(2)=(yyy(ipoint)-yyy(i))+onethird*regq
    utang(3)=(yzz(ipoint)-yzz(i))+onethird*regq
    norm = modulvec(utang)+regq
    versor(1)=utang(1)/norm
    versor(2)=utang(2)/norm
    versor(3)=utang(3)/norm
    coulomelec(1:3)=coulomelec(1:3)+(jetch(i)*Qt)/(norm**2.d0)* &
     versor(1:3)
  enddo
  do i=ipoint+1,ienf
    utang(1)=(yxx(ipoint)-yxx(i))+onethird*regq
    utang(2)=(yyy(ipoint)-yyy(i))+onethird*regq
    utang(3)=(yzz(ipoint)-yzz(i))+onethird*regq
    norm = modulvec(utang)+regq
    versor(1)=utang(1)/norm
    versor(2)=utang(2)/norm
    versor(3)=utang(3)/norm
    coulomelec(1:3)=coulomelec(1:3)+(jetch(i)*Qt)/(norm**2.d0)* &
     versor(1:3)
  enddo
  
     
  return
  
 end subroutine compute_coulomelec3d
 
 subroutine compute_coulomelec(ycf,yxx,yyy,yzz)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb forces of a
!     n-body system by the direct summation method
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  double precision, parameter :: onethird=1.d0/(dsqrt(3.d0))
  
  integer :: ipoint,jpoint,imiomax
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint
  double precision, dimension(3) :: versor,utang
  
  
  
  select case(systype)
    case(1)
    
      imiomax=npjet
      if(ncoulforce/=0)then
        if(imiomax>ncoulforce)then
          deallocate(coulforce)
          ncoulforce=imiomax+100
          allocate(ycf(0:ncoulforce,1))
        endif
      else
        ncoulforce=imiomax+100
        allocate(ycf(0:ncoulforce,1))
      endif
  
      ycf(0:ncoulforce,1:1)=0.d0
      
!     compute the Coulomb forces
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        do jpoint=ipoint+1,npjet
          if(jetfr(jpoint))cycle
          ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
           ((dabs(yxx(jpoint)-yxx(ipoint))+regq)**2.d0)
          ycf(jpoint,1)=ycf(jpoint,1)-1.d0*(jetch(jpoint)*Qt)/ &
           ((dabs(yxx(jpoint)-yxx(ipoint))+regq)**2.d0)
        enddo
        if(lmirror)then
          do jpoint=inpjet,npjet
            if(jetfr(jpoint))cycle
            xjpoint=dabs(yxx(jpoint)-h)+h
            ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
             ((dabs(xjpoint-yxx(ipoint))+regq)**2.d0)
          enddo
        endif
      enddo
      
      imiomax=ncoulforce+1
      call sum_world_darr(ycf,imiomax)
      
    case default
      
      imiomax=npjet
      if(ncoulforce/=0)then
        if(imiomax>ncoulforce)then
          deallocate(coulforce)
          ncoulforce=imiomax+100
          allocate(ycf(0:ncoulforce,3))
        endif
      else
        ncoulforce=imiomax+100
        allocate(ycf(0:ncoulforce,3))
      endif
      
      ycf(0:ncoulforce,1:3)=0.d0
      
!     compute the Coulomb forces
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        do jpoint=ipoint+1,npjet
          if(jetfr(jpoint))cycle
          utang(1)=yxx(ipoint)-yxx(jpoint)
          utang(2)=yyy(ipoint)-yyy(jpoint)
          utang(3)=yzz(ipoint)-yzz(jpoint)
          norm = modulvec(utang)
          if(norm>1.d-30)then
            versor(1)=utang(1)/norm
            versor(2)=utang(2)/norm
            versor(3)=utang(3)/norm
            ycf(ipoint,1:3)=ycf(ipoint,1:3)+ &
             (jetch(jpoint)*Qt)/((norm+regq)**2.d0)*versor(1:3)
            ycf(jpoint,1:3)=ycf(jpoint,1:3)- &
             (jetch(jpoint)*Qt)/((norm+regq)**2.d0)*versor(1:3)
          endif
        enddo
        if(lmirror)then
          do jpoint=inpjet,npjet
            if(jetfr(jpoint))cycle
            xjpoint=dabs(yxx(jpoint)-h)+h
            yjpoint=yyy(jpoint)
            zjpoint=yzz(jpoint)
            utang(1)=yxx(ipoint)-xjpoint
            utang(2)=yyy(ipoint)-yjpoint
            utang(3)=yzz(ipoint)-zjpoint
            norm = modulvec(utang)
            if(norm>1.d-30)then
              versor(1)=utang(1)/norm
              versor(2)=utang(2)/norm
              versor(3)=utang(3)/norm
              ycf(ipoint,1:3)=ycf(ipoint,1:3)- &
               (jetch(jpoint)*Qt)/((norm+regq)**2.d0)*versor(1:3)
            endif
          enddo
        endif
      enddo
      
      imiomax=(ncoulforce+1)*3
      call sum_world_darr(ycf,imiomax)
      
  end select
  
       
  return
  
 end subroutine compute_coulomelec
 
 end module coulomb_force_mod
 
