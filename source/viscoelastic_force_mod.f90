 module viscoelastic_force_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which compute
!     the viscoelastic force 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
 use utility_mod, only : Pi,modulvec,cross,dot
 
 implicit none
 
 private
 
 public :: viscot_3d
 
 contains
 
 function viscot_3d(ipoint,ivar,yst,versor,crosssub)
 
!***********************************************************************
!     
!     JETSPIN function for computing the viscoelastic force in the
!     three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint,ivar
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, dimension(3), intent(in) :: versor
  double precision, intent(in) :: crosssub
  
  double precision :: viscot_3d 
  
  
  viscot_3d=Pi*crosssub*yst(ipoint)*versor(ivar)
    
  return
 
 end function viscot_3d
 
 end module viscoelastic_force_mod
 
