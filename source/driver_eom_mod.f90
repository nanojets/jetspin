 
 module driver_eom_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which manage the call
!     for the subroutines computing the first derivatives of the
!     equations of motion
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
 use error_mod,   only : error
 use utility_mod, only : Pi,modulvec,cross,dot
 use nanojet_mod, only : systype
 use eom_mod,     only : eom1,eom3,eom4,eom4_pos,eom4_stress, &
                          eom1_KV_pos_v,eom1_KV_st, &
                          eom3_KV_pos_v,eom3_KV_st
 use eom_ev_mod,  only : eom1_ev,eom3_ev,eom4_ev,eom4_pos_ev, &
                          eom4_stress_ev, &
                          eom1_KV_pos_v_ev,eom1_KV_st_ev, &
                          eom3_KV_pos_v_ev,eom3_KV_st_ev
                        
 
 implicit none
 
 private
 
 public :: xpsys
 public :: xpsys_pos
 public :: xpsys_stress
 public :: xpsys_KV_pos_v
 public :: xpsys_KV_st
 public :: xpsys_ev
 public :: xpsys_ev_maxwell
 public :: xpsys_pos_ev
 public :: xpsys_stress_ev
 
 contains
 
 subroutine xpsys(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k,fstocvx,fstocvy,fstocvz) 
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute the first derivatives of the system 
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
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fxx
  double precision, intent(inout) ::  fyy
  double precision, intent(inout) ::  fzz
  double precision, intent(inout) ::  fst
  double precision, intent(inout) ::  fvx
  double precision, intent(inout) ::  fvy
  double precision, intent(inout) ::  fvz
  double precision, optional, intent(inout) ::  fstocvx
  double precision, optional, intent(inout) ::  fstocvy
  double precision, optional, intent(inout) ::  fstocvz
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  select case(systype)
    case (1)
      call eom1(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k)
    case (3)
      call eom3(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k)
    case (4)
      call eom4(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k,fstocvx,fstocvy,fstocvz)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys 
 
 subroutine xpsys_pos(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,timesub,k) 
       
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute only the position first derivatives 
!     of the system 
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
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fxx
  double precision, intent(inout) ::  fyy
  double precision, intent(inout) ::  fzz
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  select case(systype)
    case (4)
      call eom4_pos(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,timesub,k)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_pos
 
 subroutine xpsys_stress(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute only the stress first derivative 
!     of the system 
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
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fst
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  select case(systype)
    case (4)
      call eom4_stress(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fst,timesub,k)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_stress
 
 subroutine xpsys_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fvx,fvy,fvz,timesub,k) 
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute the first derivatives of the system
!     whenever the Kelvin–Voigt model is activated
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fxx
  double precision, intent(inout) ::  fyy
  double precision, intent(inout) ::  fzz
  double precision, intent(inout) ::  fvx
  double precision, intent(inout) ::  fvy
  double precision, intent(inout) ::  fvz
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  select case(systype)
    case (1)
      call eom1_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fvx,fvy,fvz,timesub,k)
    case (3)
      call eom3_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fvx,fvy,fvz,timesub,k)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_KV_pos_v
 
 subroutine xpsys_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       yax,yay,yaz,fst,timesub,k) 
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute the first derivatives of the system
!     whenever the Kelvin–Voigt model is activated
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: ipoint
  double precision, allocatable, dimension (:), intent(in) ::  yxx
  double precision, allocatable, dimension (:), intent(in) ::  yyy
  double precision, allocatable, dimension (:), intent(in) ::  yzz
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, allocatable, dimension (:), intent(in) ::  yax
  double precision, allocatable, dimension (:), intent(in) ::  yay
  double precision, allocatable, dimension (:), intent(in) ::  yaz
  double precision, intent(inout) ::  fst
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  select case(systype)
    case (1)
      call eom1_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       yax,yay,yaz,fst,timesub,k)
    case (3)
      call eom3_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       yax,yay,yaz,fst,timesub,k) 
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_KV_st
 
 subroutine xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,timesub,k,fstocvx,fstocvy, &
       fstocvz) 
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute the first derivatives of the system 
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
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:), intent(in) ::  yve
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fxx
  double precision, intent(inout) ::  fyy
  double precision, intent(inout) ::  fzz
  double precision, intent(inout) ::  fst
  double precision, intent(inout) ::  fvx
  double precision, intent(inout) ::  fvy
  double precision, intent(inout) ::  fvz
  double precision, intent(inout) ::  fve
  double precision, optional, intent(inout) ::  fstocvx
  double precision, optional, intent(inout) ::  fstocvy
  double precision, optional, intent(inout) ::  fstocvz
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  select case(systype)
    case (1)
      call eom1_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,timesub,k)
    case (3)
      call eom3_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,timesub,k)
    case (4)
      call eom4_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,timesub,k,fstocvx,fstocvy, &
       fstocvz)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_ev

 subroutine xpsys_ev_maxwell(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,timesub,k)
  implicit none
  integer, intent(in) :: ipoint,k
  double precision, allocatable, intent(in) :: yxx(:),yyy(:),yzz(:),yst(:)
  double precision, allocatable, intent(in) :: yvx(:),yvy(:),yvz(:),yvl(:),yve(:)
  double precision, allocatable, intent(in) :: ycf(:,:)
  double precision, intent(inout) :: fxx,fyy,fzz,fst,fvx,fvy,fvz,fve
  double precision, intent(in) :: timesub
  call xpsys_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
   fxx,fyy,fzz,fst,fvx,fvy,fvz,fve,timesub,k)
 end subroutine xpsys_ev_maxwell

 subroutine xpsys_pos_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve, &
       ycf,fxx,fyy,fzz,fve,timesub,k) 
       
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute only the position first derivatives 
!     of the system 
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
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:), intent(in) ::  yve
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fxx
  double precision, intent(inout) ::  fyy
  double precision, intent(inout) ::  fzz
  double precision, intent(inout) ::  fve
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  select case(systype)
    case (4)
      call eom4_pos_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fve,timesub,k)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_pos_ev
 
 subroutine xpsys_stress_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve, &
       ycf,fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which compute only the stress first derivative 
!     of the system 
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
  double precision, allocatable, dimension (:), intent(in) ::  yst
  double precision, allocatable, dimension (:), intent(in) ::  yvx
  double precision, allocatable, dimension (:), intent(in) ::  yvy
  double precision, allocatable, dimension (:), intent(in) ::  yvz
  double precision, allocatable, dimension (:), intent(in) ::  yvl
  double precision, allocatable, dimension (:), intent(in) ::  yve
  double precision, allocatable, dimension (:,:), intent(in) ::  ycf
  double precision, intent(inout) ::  fst
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  select case(systype)
    case (4)
      call eom4_stress_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve, &
       ycf,fst,timesub,k)
    case default
      call error(2)
  end select
  
  return
  
 end subroutine xpsys_stress_ev

 end module driver_eom_mod
