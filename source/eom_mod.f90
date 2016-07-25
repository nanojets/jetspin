 
 module eom_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which compute 
!     the first derivatives of the system 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
 use version_mod, only : idrank
 use utility_mod,           only : Pi,modulvec,cross,dot
 use nanojet_mod,           only : jetms,jetch,inpjet,npjet, &
                             consistency,findex,yieldstress, &
                             liniperturb,linserted,pfreq,att,fve, &
                             gr,ks,li,lrg,v,jetfr
 use support_functions_mod, only : compute_geometry, &
                             compute_tangetversor, &
                             project_beadveltangetversor, &
                             compute_geometry_init, &
                             compute_curvcenter,compute_curvature,&
                             project_veltangetversor, &
                             compute_stocforce_3d, &
                             compute_geometry_1d, &
                             compute_geometry_1d_init,upwall, &
                             compute_lorentz_acc, &
                             project_beadacctangetversor, &
                             compute_geometry_1d_kv
 use electric_field_mod,    only : driver_electric_field
 
 implicit none
 
 private
 
 public :: eom1
 public :: eom3
 public :: eom4
 public :: eom4_pos
 public :: eom4_stress
 public :: eom1_KV_pos_v
 public :: eom1_KV_st
 public :: eom3_KV_pos_v
 public :: eom3_KV_st
 
 contains
 
 subroutine eom1(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the one dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
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
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  double precision :: beadlendown,beadlenup,beadvelup
  
  double precision :: Vtvec(3),Fvet,coulomelec
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fst=0.d0
    fvx=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    call compute_geometry_1d_init(ipoint,yxx,yst,yvx,beadlenup, &
     beadvelup)
    coulomelec=ycf(ipoint,1)
    call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
    
    
    Fvet=Fve/jetms(ipoint)
    fxx = yvx(ipoint) 
    fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
     yst(ipoint)
    fvx = Gr+Vtvec(1)-Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)+ &
     coulomelec+upwall(ipoint,yxx)
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown, &
         beadlenup,beadvelup)
        coulomelec=ycf(ipoint,1)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
        
        
        Fvet=Fve/jetms(ipoint)
        fxx = yvx(ipoint) 
        fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
         yst(ipoint)
        fvx = Gr+Vtvec(1)-Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)+ &
         Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
          upwall(ipoint,yxx)
      else
        fxx=0.d0
        fst=0.d0
        fvx=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    fxx=0.d0
    fst=0.d0
    fvx=0.d0
    return
  endif
  
! ordinary case
  call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown,beadlenup, &
   beadvelup)
  coulomelec=ycf(ipoint,1)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
  
  
  Fvet=Fve/jetms(ipoint)
  fxx = yvx(ipoint) 
  fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
   yst(ipoint)
  fvx = Gr+Vtvec(1)-Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)+ &
   Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
   upwall(ipoint,yxx)

  return
  
 end subroutine eom1
 
 subroutine eom1_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fvx,fvy,fvz,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the one dimensional model
!     with Kelvin–Voigt model activated
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
  
  double precision :: beadlendown,beadlenup,beadvelup
  
  double precision :: Vtvec(3),Fvet,coulomelec
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fvx=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    call compute_geometry_1d_init(ipoint,yxx,yst,yvx,beadlenup, &
     beadvelup)
    coulomelec=ycf(ipoint,1)
    call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
    
    
    Fvet=Fve/jetms(ipoint)
    fxx = yvx(ipoint) 
    fvx = Gr+Vtvec(1)-Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)+ &
     coulomelec+upwall(ipoint,yxx)
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown, &
         beadlenup,beadvelup)
        coulomelec=ycf(ipoint,1)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
        
        
        Fvet=Fve/jetms(ipoint)
        fxx = yvx(ipoint) 
        fvx = Gr+Vtvec(1)-Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)+ &
         Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
          upwall(ipoint,yxx)
      else
        fxx=0.d0
        fvx=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    fxx=0.d0
    fvx=0.d0
    return
  endif
  
! ordinary case
  call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown,beadlenup, &
   beadvelup)
  coulomelec=ycf(ipoint,1)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
  
  
  Fvet=Fve/jetms(ipoint)
  fxx = yvx(ipoint) 
  fvx = Gr+Vtvec(1)-Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)+ &
   Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
   upwall(ipoint,yxx)

  return
  
 end subroutine eom1_KV_pos_v
 
 subroutine eom1_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       yax,yay,yaz,fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the one dimensional model 
!     with Kelvin–Voigt model activated
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
  
  double precision :: beadlenup,beadvelup,beadaccup
  
  double precision :: Vtvec(3),Fvet,coulomelec
  
  
! special cases
  
  if(jetfr(ipoint))then
    fst=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    call compute_geometry_1d_KV(ipoint,yxx,yst,yvx,yax,beadlenup, &
     beadvelup,beadaccup)
    fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_1d_KV(ipoint,yxx,yst,yvx,yax, &
         beadlenup,beadvelup,beadaccup)
        fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)
      else
        fst=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    fst=0.d0
    return
  endif
  
! ordinary case
  call compute_geometry_1d_KV(ipoint,yxx,yst,yvx,yax,beadlenup, &
   beadvelup,beadaccup)
  fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)

  return
  
 end subroutine eom1_KV_st
 
 subroutine eom3(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
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
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision :: curvature
  double precision, dimension(3) :: tangentversorup,tangentversordown
  double precision, dimension(3) :: curvcenter,vcurvature,coulomelec
  double precision :: aLorx,aLory,aLorz
  logical :: lstraight
  
  double precision :: Vtvec(3),Fvet,Kst,factor1,factor2,factor3
  
  integer,save :: ij
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fst=0.d0
    fvx=0.d0
    fvy=0.d0
    fvz=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
      
      Fvet=Fve/jetms(ipoint)
      
      factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
      
      fxx = yvx(ipoint) 
      fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
       yst(ipoint)
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+coulomelec(3)+aLorz
      
    else
      call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
      call compute_curvature(ipoint,yxx,yyy,yzz,curvature,vcurvature, &
       curvcenter,lstraight)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
      
      Fvet=Fve/jetms(ipoint)
      Kst=Ks/jetms(ipoint)
      
      factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
      factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
       (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
      
      fxx = yvx(ipoint) 
      fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
       yst(ipoint)
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
       Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+Kst*curvature*factor3* &
       vcurvature(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+Kst*curvature*factor3* &
       vcurvature(3)+coulomelec(3)+aLorz
        
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
        call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
         beadlenup)
        call compute_tangetversor(ipoint-1,yxx,yyy,yzz, &
         tangentversordown,beadlendown)
        call project_beadveltangetversor(ipoint,yvx,yvy,yvz, &
         beadvelup,tangentversorup)
        call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
        call compute_curvature(ipoint,yxx,yyy,yzz,curvature, &
         vcurvature,curvcenter,lstraight)
        coulomelec(1:3)=ycf(ipoint,1:3)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
        call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
        
        Fvet=Fve/jetms(ipoint)
        Kst=Ks/jetms(ipoint)
      
        factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
        factor2=Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)
        factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
         (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
        
        fxx = yvx(ipoint) 
        fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
         yst(ipoint)
        fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
         factor2*tangentversordown(1)+ &
         Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
         upwall(ipoint,yxx)+aLorx
        
        fyy = yvy(ipoint) 
        fvy =Vtvec(2)-factor1*tangentversorup(2)+ &
         factor2*tangentversordown(2)+ &
         Kst*curvature*factor3*vcurvature(2)+coulomelec(2)+aLory
        
        fzz = yvz(ipoint) 
        fvz =Vtvec(3)-factor1*tangentversorup(3)+ &
         factor2*tangentversordown(3)+ &
         Kst*curvature*factor3*vcurvature(3)+coulomelec(3)+aLorz
        
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fst=0.d0
        fvx=0.d0
        fvy=0.d0
        fvz=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fxx=0.d0
      fyy = -1.d0*pfreq*yzz(ipoint) 
      fzz = pfreq*yyy(ipoint)
      fst=0.d0
      fvx=0.d0
      fvy = -1.d0*pfreq**2.d0*yyy(ipoint)
      fvz = -1.d0*pfreq**2.d0*yzz(ipoint) 
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
      fst=0.d0
      fvx=0.d0
      fvy=0.d0
      fvz=0.d0
    endif
    return
  endif
   
! ordinary case
  call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
  call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
   beadlenup)
  call compute_tangetversor(ipoint-1,yxx,yyy,yzz,tangentversordown, &
   beadlendown)
  call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
   tangentversorup)
  call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
  call compute_curvature(ipoint,yxx,yyy,yzz,curvature,vcurvature, &
   curvcenter,lstraight)
  coulomelec(1:3)=ycf(ipoint,1:3)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
  call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
  
  Fvet=Fve/jetms(ipoint)
  Kst=Ks/jetms(ipoint)
  
  factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
  factor2=Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)
  factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
   (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
  
  fxx = yvx(ipoint) 
  fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
   yst(ipoint)
  fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
   factor2*tangentversordown(1)+ &
   Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
   upwall(ipoint,yxx)+aLorx
  
  fyy = yvy(ipoint) 
  fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
   factor2*tangentversordown(2)+ &
   Kst*curvature*factor3*vcurvature(2)+coulomelec(2)+aLory
  
  fzz = yvz(ipoint) 
  fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
   factor2*tangentversordown(3)+ &
   Kst*curvature*factor3*vcurvature(3)+coulomelec(3)+aLorz
  
  
  return
  
 end subroutine eom3 
 
 subroutine eom3_KV_pos_v(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fvx,fvy,fvz,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional model
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
  
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision :: curvature
  double precision, dimension(3) :: tangentversorup,tangentversordown
  double precision, dimension(3) :: curvcenter,vcurvature,coulomelec
  double precision :: aLorx,aLory,aLorz
  logical :: lstraight
  
  double precision :: Vtvec(3),Fvet,Kst,factor1,factor2,factor3
  
  integer,save :: ij
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fvx=0.d0
    fvy=0.d0
    fvz=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
      
      Fvet=Fve/jetms(ipoint)
      
      factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
      
      fxx = yvx(ipoint) 
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+coulomelec(3)+aLorz
      
    else
      call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
      call compute_curvature(ipoint,yxx,yyy,yzz,curvature,vcurvature, &
       curvcenter,lstraight)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
      
      Fvet=Fve/jetms(ipoint)
      Kst=Ks/jetms(ipoint)
      
      factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
      factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
       (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
      
      fxx = yvx(ipoint) 
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
       Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+Kst*curvature*factor3* &
       vcurvature(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+Kst*curvature*factor3* &
       vcurvature(3)+coulomelec(3)+aLorz
        
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
        call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
         beadlenup)
        call compute_tangetversor(ipoint-1,yxx,yyy,yzz, &
         tangentversordown,beadlendown)
        call project_beadveltangetversor(ipoint,yvx,yvy,yvz, &
         beadvelup,tangentversorup)
        call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
        call compute_curvature(ipoint,yxx,yyy,yzz,curvature, &
         vcurvature,curvcenter,lstraight)
        coulomelec(1:3)=ycf(ipoint,1:3)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
        call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
        
        Fvet=Fve/jetms(ipoint)
        Kst=Ks/jetms(ipoint)
      
        factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
        factor2=Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)
        factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
         (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
        
        fxx = yvx(ipoint) 
        fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
         factor2*tangentversordown(1)+ &
         Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
         upwall(ipoint,yxx)+aLorx
        
        fyy = yvy(ipoint) 
        fvy =Vtvec(2)-factor1*tangentversorup(2)+ &
         factor2*tangentversordown(2)+ &
         Kst*curvature*factor3*vcurvature(2)+coulomelec(2)+aLory
        
        fzz = yvz(ipoint) 
        fvz =Vtvec(3)-factor1*tangentversorup(3)+ &
         factor2*tangentversordown(3)+ &
         Kst*curvature*factor3*vcurvature(3)+coulomelec(3)+aLorz
        
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fvx=0.d0
        fvy=0.d0
        fvz=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fxx=0.d0
      fyy = -1.d0*pfreq*yzz(ipoint) 
      fzz = pfreq*yyy(ipoint)
      fvx=0.d0
      fvy = -1.d0*pfreq**2.d0*yyy(ipoint)
      fvz = -1.d0*pfreq**2.d0*yzz(ipoint) 
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
      fvx=0.d0
      fvy=0.d0
      fvz=0.d0
    endif
    return
  endif
   
! ordinary case
  call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
  call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
   beadlenup)
  call compute_tangetversor(ipoint-1,yxx,yyy,yzz,tangentversordown, &
   beadlendown)
  call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
   tangentversorup)
  call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
  call compute_curvature(ipoint,yxx,yyy,yzz,curvature,vcurvature, &
   curvcenter,lstraight)
  coulomelec(1:3)=ycf(ipoint,1:3)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
  call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
  
  Fvet=Fve/jetms(ipoint)
  Kst=Ks/jetms(ipoint)
  
  factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
  factor2=Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)
  factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
   (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
  
  fxx = yvx(ipoint) 
  fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
   factor2*tangentversordown(1)+ &
   Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
   upwall(ipoint,yxx)+aLorx
  
  fyy = yvy(ipoint) 
  fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
   factor2*tangentversordown(2)+ &
   Kst*curvature*factor3*vcurvature(2)+coulomelec(2)+aLory
  
  fzz = yvz(ipoint) 
  fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
   factor2*tangentversordown(3)+ &
   Kst*curvature*factor3*vcurvature(3)+coulomelec(3)+aLorz
  
  
  return
  
 end subroutine eom3_KV_pos_v
 
 subroutine eom3_KV_st(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       yax,yay,yaz,fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional model
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
  
  
  double precision :: beadlenup,beadvelup,beadaccup
  double precision :: curvature
  double precision, dimension(3) :: tangentversorup,tangentversordown
  double precision, dimension(3) :: curvcenter,vcurvature,coulomelec
  double precision :: aLorx,aLory,aLorz
  logical :: lstraight
  
  double precision :: Vtvec(3),Fvet,Kst,factor1,factor2,factor3
  
  integer,save :: ij
  
  
! special cases
  
  if(jetfr(ipoint))then
    fst=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call project_beadacctangetversor(ipoint,yax,yay,yaz,beadaccup, &
       tangentversorup)
      
      fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)
      
    else
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call project_beadacctangetversor(ipoint,yax,yay,yaz,beadaccup, &
       tangentversorup)
      
      fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)
        
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
        call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
         beadlenup)
        call project_beadveltangetversor(ipoint,yvx,yvy,yvz, &
         beadvelup,tangentversorup)
        call project_beadacctangetversor(ipoint,yax,yay,yaz,beadaccup, &
         tangentversorup)
        
        fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)
        
      else
        fst=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fst=0.d0
    else
      fst=0.d0
    endif
    return
  endif
   
! ordinary case
  call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
  call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
   beadlenup)
  call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
   tangentversorup)
  call project_beadacctangetversor(ipoint,yax,yay,yaz,beadaccup, &
   tangentversorup)
  
  fst = (beadvelup/beadlenup)+(beadaccup/beadlenup)
  
  
  return
  
 end subroutine eom3_KV_st
 
 subroutine eom4(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,timesub,k,fstocvx,fstocvy,fstocvz) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional stochastic model
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
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
  double precision, optional ::  fstocvx
  double precision, optional ::  fstocvy
  double precision, optional ::  fstocvz
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision :: curvature,veltangent
  double precision, dimension(3) :: tangentversorup,tangentversordown, &
   friction
  double precision, dimension(3) :: curvcenter,vcurvature,coulomelec
  logical :: lstraight
  
  double precision :: Vtvec(3),Fvet,Kst,attt,Lit,factor1,factor2,factor3
  double precision :: factor4,factor5
  
  double precision :: aLorx,aLory,aLorz
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fst=0.d0
    fvx=0.d0
    fvy=0.d0
    fvz=0.d0
    fstocvx=0.d0
    fstocvy=0.d0
    fstocvz=0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
       tangentversorup)
      call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
      
      Fvet=Fve/jetms(ipoint)
      attt=att/jetms(ipoint)
      
      if(yst(ipoint)>0.d0)then
        factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
      else
        factor1=0.d0
      endif
      factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
  
      fxx = yvx(ipoint) 
      fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
       yst(ipoint)
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+coulomelec(1)- &
       factor4*tangentversorup(1)+upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+coulomelec(2)- &
       factor4*tangentversorup(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+coulomelec(3)- &
       factor4*tangentversorup(3)+aLorz
      
    else
      
      call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
      call compute_curvature(ipoint,yxx,yyy,yzz,curvature,vcurvature, &
       curvcenter,lstraight)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
       tangentversorup)
      call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
      
      Fvet=Fve/jetms(ipoint)
      Kst=Ks/jetms(ipoint)
      attt=att/jetms(ipoint)
      Lit=Li/jetms(ipoint)
      
      if(yst(ipoint)>0.d0)then
        factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
      else
        factor1=0.d0
      endif
      factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
       (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
      factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
      factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
      
      fxx = yvx(ipoint) 
      fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
       yst(ipoint)
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
       Kst*curvature*factor3*vcurvature(1)+ &
       coulomelec(1)-factor4*tangentversorup(1)- &
       Lit*factor5*vcurvature(1)+upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
       Kst*curvature*factor3*vcurvature(2)+ &
       coulomelec(2)-factor4*tangentversorup(2)- &
       Lit*factor5*vcurvature(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
       Kst*curvature*factor3*vcurvature(3)+ &
       coulomelec(3)-factor4*tangentversorup(3)- &
       Lit*factor5*vcurvature(3)+aLorz
      
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        
        call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
        call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
         beadlenup)
        call compute_tangetversor(ipoint-1,yxx,yyy,yzz, &
         tangentversordown,beadlendown)
        call project_beadveltangetversor(ipoint,yvx,yvy,yvz, &
         beadvelup,tangentversorup)
        call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
        call compute_curvature(ipoint,yxx,yyy,yzz,curvature, &
         vcurvature,curvcenter,lstraight)
        coulomelec(1:3)=ycf(ipoint,1:3)
        call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
         tangentversorup)
        call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
        call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
        
        Fvet=Fve/jetms(ipoint)
        Kst=Ks/jetms(ipoint)
        attt=att/jetms(ipoint)
        Lit=Li/jetms(ipoint)
        
        if(yst(ipoint)>0.d0)then
          factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
        else
          factor1=0.d0
        endif
        if(yst(ipoint-1)>0.d0)then
          factor2=Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)
        else
          factor2=0.d0
        endif
        factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
         (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
        factor4=attt*(dabs(beadlenup)**0.905d0)* &
         (dabs(veltangent)**1.19d0)
        factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
        
        fxx = yvx(ipoint) 
        fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
         yst(ipoint)
        fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
         factor2*tangentversordown(1)+ &
         Kst*curvature*factor3*vcurvature(1)+coulomelec(1)- &
         factor4*tangentversorup(1)- &
         Lit*factor5*vcurvature(1)+upwall(ipoint,yxx)+aLorx
        
        fyy = yvy(ipoint) 
        fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
         factor2*tangentversordown(2)+ &
         Kst*curvature*factor3*vcurvature(2)+coulomelec(2)- &
         factor4*tangentversorup(2)- &
         Lit*factor5*vcurvature(2)+aLory
        
        fzz = yvz(ipoint) 
        fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
         factor2*tangentversordown(3)+ &
         Kst*curvature*factor3*vcurvature(3)+coulomelec(3)- &
         factor4*tangentversorup(3)- &
         Lit*factor5*vcurvature(3)+aLorz
       
        
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fst=0.d0
        fvx=0.d0
        fvy=0.d0
        fvz=0.d0
        fstocvx=0.d0
        fstocvy=0.d0
        fstocvz=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fxx=0.d0
      fyy = -1.d0*pfreq*yzz(ipoint) 
      fzz = pfreq*yyy(ipoint)
      fst=0.d0
      fvx=0.d0
      fvy = -1.d0*pfreq**2.d0*yyy(ipoint)
      fvz = -1.d0*pfreq**2.d0*yzz(ipoint) 
      fstocvx=0.d0
      fstocvy=0.d0
      fstocvz=0.d0
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
      fst=0.d0
      fvx=0.d0
      fvy=0.d0
      fvz=0.d0
      fstocvx=0.d0
      fstocvy=0.d0
      fstocvz=0.d0
    endif
    return
  endif
  
  
! ordinary case
  call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
  call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
   beadlenup)
  call compute_tangetversor(ipoint-1,yxx,yyy,yzz,tangentversordown, &
   beadlendown)
  call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
   tangentversorup)
  call compute_curvcenter(ipoint,yxx,yyy,yzz,curvcenter,lstraight)       
  call compute_curvature(ipoint,yxx,yyy,yzz,curvature,vcurvature, &
   curvcenter,lstraight)
  coulomelec(1:3)=ycf(ipoint,1:3)
  call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
   tangentversorup)
  call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec)
  call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz)
  
  Fvet=Fve/jetms(ipoint)
  Kst=Ks/jetms(ipoint)
  attt=att/jetms(ipoint)
  Lit=Li/jetms(ipoint)
  
  if(yst(ipoint)>0.d0)then
    factor1=Fvet*yvl(ipoint)*(yst(ipoint)/beadlenup)
  else
    factor1=0.d0
  endif
  if(yst(ipoint-1)>0.d0)then
    factor2=Fvet*yvl(ipoint-1)*(yst(ipoint-1)/beadlendown)
  else
    factor2=0.d0
  endif
  factor3=0.25d0*((dsqrt(yvl(ipoint))/dsqrt(beadlenup))+ &
   (dsqrt(yvl(ipoint-1))/dsqrt(beadlendown)))**2.d0
  factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
  factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
  
  fxx = yvx(ipoint) 
  fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
   yst(ipoint)
  fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
   factor2*tangentversordown(1)+ &
   Kst*curvature*factor3*vcurvature(1)+coulomelec(1)- &
   factor4*tangentversorup(1)-Lit*factor5*vcurvature(1)+ &
   upwall(ipoint,yxx)+aLorx
  
  fyy = yvy(ipoint) 
  fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
   factor2*tangentversordown(2)+ &
   Kst*curvature*factor3*vcurvature(2)+coulomelec(2)- &
   factor4*tangentversorup(2)-Lit*factor5*vcurvature(2)+aLory
  
  fzz = yvz(ipoint) 
  fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
   factor2*tangentversordown(3)+ &
   Kst*curvature*factor3*vcurvature(3)+coulomelec(3)- &
   factor4*tangentversorup(3)-Lit*factor5*vcurvature(3)+aLorz
  
  
  return
  
 end subroutine eom4
 
 subroutine eom4_pos(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fxx,fyy,fzz,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing only the position first 
!     derivatives of the system for the three dimensional stochastic 
!     model
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
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx = 0.d0
    fyy = 0.d0
    fzz = 0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      fxx = yvx(ipoint)
      fyy = yvy(ipoint) 
      fzz = yvz(ipoint) 
    else
      fxx = yvx(ipoint)
      fyy = yvy(ipoint) 
      fzz = yvz(ipoint) 
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        fxx = yvx(ipoint)
        fyy = yvy(ipoint) 
        fzz = yvz(ipoint) 
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fxx=0.d0
      fyy = -1.d0*pfreq*yzz(ipoint) 
      fzz = pfreq*yyy(ipoint)
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
    endif
    return
  endif
  
  
! ordinary case
  fxx = yvx(ipoint)
  fyy = yvy(ipoint) 
  fzz = yvz(ipoint) 
  
  
  return
  
 end subroutine eom4_pos
 
  subroutine eom4_stress(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,ycf, &
       fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing only the stress first 
!     derivative of the system for the three dimensional stochastic 
!     model
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
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision, dimension(3) :: tangentversorup
  
  
! special cases
  
  if(jetfr(ipoint))then
    fst = 0.d0
    return
  endif
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      
      fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
       yst(ipoint)
      
    else
      call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      
      fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
       yst(ipoint)
      
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
        call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
         beadlenup)
       
        call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
         tangentversorup)
       
        fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
         yst(ipoint)
        
      else
        fst=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fst=0.d0
    else
      fst=0.d0
    endif
    return
  endif
  
  
! ordinary case
  call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
  call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
   beadlenup)
  
  call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
   tangentversorup)
  
  fst = yieldstress+consistency*(beadvelup/beadlenup)**findex- &
   yst(ipoint)
  
  
  return
  
 end subroutine eom4_stress

 end module eom_mod


