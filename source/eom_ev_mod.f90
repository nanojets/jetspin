 
 module eom_ev_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which compute 
!     the first derivatives of the system including evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!***********************************************************************
 
 use version_mod, only : idrank,finalize_world
 use utility_mod,           only : Pi,modulvec,cross,dot
 use nanojet_mod,           only : jetms,jetch,inpjet,npjet, &
                             consistency,findex,yieldstress, &
                             liniperturb,linserted,pfreq,att,fve, &
                             gr,ks,li,lrg,v,jetfr,noisefric,cp0,Bev, &
                             mev,evairv,evmasscoeff,sqrevsc,evumidity, &
                             evcsvapour,tev,lengthscale,tao,lairdrag,evlim
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
 
 public :: eom1_ev
 public :: eom3_ev
 public :: eom4_ev
 public :: eom4_pos_ev
 public :: eom4_stress_ev
 public :: eom1_KV_pos_v_ev
 public :: eom1_KV_st_ev
 public :: eom3_KV_pos_v_ev
 public :: eom3_KV_st_ev
 
 contains
 
 subroutine eom1_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the one dimensional model with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, intent(inout) ::  fev
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  double precision :: beadlendown,beadlenup,beadvelup
  
  double precision :: Vtvec(3),Fvet,coulomelec
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,Re
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fst=0.d0
    fvx=0.d0
    fev=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  
  
  if(ipoint==inpjet)then
    call compute_geometry_1d_init(ipoint,yxx,yst,yvx,beadlenup, &
     beadvelup)
    coulomelec=ycf(ipoint,1)
    call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
    
    !compute Reynolds number
    Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*dabs(yvx(ipoint)))/evairv
    Fvet=Fve/(jetms(ipoint)*cmass)
    fxx = yvx(ipoint) 
    fst = (1.d0/rattao)*(yieldstress+ &
     consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
    fvx = Gr+Vtvec(1)-Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)+ &
     coulomelec+upwall(ipoint,yxx)
    
    fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
     sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown, &
         beadlenup,beadvelup)
        coulomelec=ycf(ipoint,1)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec, &
         cmass)
        
        !compute Reynolds number
        Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*dabs(yvx(ipoint)))/ &
         evairv
        Fvet=Fve/(jetms(ipoint)*cmass)
        fxx = yvx(ipoint) 
        fst = (1.d0/rattao)*(yieldstress+ &
         consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
        fvx = Gr+Vtvec(1)-Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)+ &
         Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
         upwall(ipoint,yxx)
        
        fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
         sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      else
        fxx=0.d0
        fst=0.d0
        fvx=0.d0
        fev=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    fxx=0.d0
    fst=0.d0
    fvx=0.d0
    fev=0.d0
    return
  endif
  
! ordinary case
  call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown,beadlenup, &
   beadvelup)
  coulomelec=ycf(ipoint,1)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
  
  !compute Reynolds number
  Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*dabs(yvx(ipoint)))/evairv
  Fvet=Fve/(jetms(ipoint)*cmass)
  fxx = yvx(ipoint) 
  fst = (1.d0/rattao)*(yieldstress+ &
   consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
  fvx = Gr+Vtvec(1)-Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)+ &
   Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
   upwall(ipoint,yxx)
  
  fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
   sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup

  return
  
 end subroutine eom1_ev
 
 subroutine eom1_KV_pos_v_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl, &
       yve,ycf,fxx,fyy,fzz,fvx,fvy,fvz,fev,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the one dimensional model
!     with Kelvin–Voigt model activated with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, intent(inout) ::  fvx
  double precision, intent(inout) ::  fvy
  double precision, intent(inout) ::  fvz
  double precision, intent(inout) ::  fev
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  double precision :: beadlendown,beadlenup,beadvelup
  
  double precision :: Vtvec(3),Fvet,coulomelec
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fvx=0.d0
    fev=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    call compute_geometry_1d_init(ipoint,yxx,yst,yvx,beadlenup, &
     beadvelup)
    coulomelec=ycf(ipoint,1)
    call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
    
    Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*dabs(yvx(ipoint)))/ &
     evairv
    Fvet=Fve/(jetms(ipoint)*cmass)
    fxx = yvx(ipoint) 
    fvx = Gr+Vtvec(1)-Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)+ &
     coulomelec+upwall(ipoint,yxx)
    
    fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
     sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown, &
         beadlenup,beadvelup)
        coulomelec=ycf(ipoint,1)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec, &
         cmass)
        
        Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*dabs(yvx(ipoint)))/ &
         evairv
        Fvet=Fve/(jetms(ipoint)*cmass)
        fxx = yvx(ipoint) 
        fvx = Gr+Vtvec(1)-Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)+ &
         Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
         upwall(ipoint,yxx)
        
        fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
         sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup
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
    fev=0.d0
    return
  endif
  
! ordinary case
  call compute_geometry_1d(ipoint,yxx,yst,yvx,beadlendown,beadlenup, &
   beadvelup)
  coulomelec=ycf(ipoint,1)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
  
  Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*dabs(yvx(ipoint)))/evairv
  Fvet=Fve/(jetms(ipoint)*cmass)
  fxx = yvx(ipoint) 
  fvx = Gr+Vtvec(1)-Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)+ &
   Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)+coulomelec+ &
   upwall(ipoint,yxx)
  
  fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
   sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup

  return
  
 end subroutine eom1_KV_pos_v_ev
 
 subroutine eom1_KV_st_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve, &
       ycf,yax,yay,yaz,fevlocal,fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the one dimensional model 
!     with Kelvin–Voigt model activated with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, allocatable, dimension (:), intent(in) ::  yax
  double precision, allocatable, dimension (:), intent(in) ::  yay
  double precision, allocatable, dimension (:), intent(in) ::  yaz
  double precision, intent(in) :: fevlocal
  double precision, intent(inout) ::  fst
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  double precision :: beadlenup,beadvelup,beadaccup
  
  double precision :: Vtvec(3),Fvet,coulomelec
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re
  
! special cases
  
  if(jetfr(ipoint))then
    fst=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    call compute_geometry_1d_KV(ipoint,yxx,yst,yvx,yax,beadlenup, &
     beadvelup,beadaccup)
    call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry_1d_KV(ipoint,yxx,yst,yvx,yax, &
         beadlenup,beadvelup,beadaccup)
        call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)
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
  call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)

  return
  
 end subroutine eom1_KV_st_ev
 
 subroutine eom3_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional model with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, intent(inout) ::  fev
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision :: curvature,veltangent
  double precision, dimension(3) :: tangentversorup,tangentversordown
  double precision, dimension(3) :: curvcenter,vcurvature,coulomelec
  double precision :: aLorx,aLory,aLorz
  logical :: lstraight
  
  double precision :: Vtvec(3),Fvet,Kst,attt,Lit,factor1,factor2,factor3
  double precision :: factor4,factor5
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re,vnorm
  
  integer,save :: ij=0
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fst=0.d0
    fvx=0.d0
    fvy=0.d0
    fvz=0.d0
    fev=0.d0
    return
  endif
  
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
       cmass)
      
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      Fvet=Fve/(jetms(ipoint)*cmass)
      
      factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
      
      fxx = yvx(ipoint) 
      fst = (1.d0/rattao)*(yieldstress+ &
       consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+coulomelec(3)+aLorz
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup
       
      if(lairdrag)then
        call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
         tangentversorup)
        attt=att/(jetms(ipoint)*cmass)
        factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
        fvx = fvx-factor4*tangentversorup(1)
        fvy = fvy-factor4*tangentversorup(2)
        fvz = fvz-factor4*tangentversorup(3)
      endif
      
      
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
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
       cmass)
      
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      Fvet=Fve/(jetms(ipoint)*cmass)
      Kst=Ks/(jetms(ipoint)*cmass)
      
      factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
      factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
       (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
      
      fxx = yvx(ipoint) 
      fst = (1.d0/rattao)*(yieldstress+ &
       consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
       Kst*curvature*factor3*vcurvature(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+Kst*curvature*factor3* &
       vcurvature(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+Kst*curvature*factor3* &
       vcurvature(3)+coulomelec(3)+aLorz
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      
      if(lairdrag)then
        call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
         tangentversorup)
        attt=att/(jetms(ipoint)*cmass)
        Lit=Li/(jetms(ipoint)*cmass)
        factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
        factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
        fvx = fvx-factor4*tangentversorup(1)-Lit*factor5*vcurvature(1)
        fvy = fvy-factor4*tangentversorup(2)-Lit*factor5*vcurvature(2)
        fvz = fvz-factor4*tangentversorup(3)-Lit*factor5*vcurvature(3)
      endif
      
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
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec, &
         cmass)
        call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
         cmass)
        
        !compute Reynolds number
        vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+ &
         yvz(ipoint)**2.d0)
        Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
        Fvet=Fve/(jetms(ipoint)*cmass)
        Kst=Ks/(jetms(ipoint)*cmass)
      
        factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
        factor2=Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)
        factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
         (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
        
        fxx = yvx(ipoint) 
        fst = (1.d0/rattao)*(yieldstress+ &
         consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
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
        
        fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
         sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
         
        if(lairdrag)then
          call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
           tangentversorup)
          attt=att/(jetms(ipoint)*cmass)
          Lit=Li/(jetms(ipoint)*cmass)
          factor4=attt*(dabs(beadlenup)**0.905d0)* &
           (dabs(veltangent)**1.19d0)
          factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
          fvx = fvx-factor4*tangentversorup(1)-Lit*factor5*vcurvature(1)
          fvy = fvy-factor4*tangentversorup(2)-Lit*factor5*vcurvature(2)
          fvz = fvz-factor4*tangentversorup(3)-Lit*factor5*vcurvature(3)
        endif
         
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fst=0.d0
        fvx=0.d0
        fvy=0.d0
        fvz=0.d0
        fev=0.d0
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
      fev=0.d0
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
      fst=0.d0
      fvx=0.d0
      fvy=0.d0
      fvz=0.d0
      fev=0.d0
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
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
  call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
   cmass)
  
  !compute Reynolds number
  vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
  Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
  Fvet=Fve/(jetms(ipoint)*cmass)
  Kst=Ks/(jetms(ipoint)*cmass)
  
  factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
  factor2=Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)
  factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
   (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
  
  fxx = yvx(ipoint) 
  fst = (1.d0/rattao)*(yieldstress+ &
   consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
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
  
  fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
   sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup
  
  if(lairdrag)then
    call project_veltangetversor(ipoint,yvx,yvy,yvz,veltangent, &
     tangentversorup)
    attt=att/(jetms(ipoint)*cmass)
    Lit=Li/(jetms(ipoint)*cmass)
    factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
    factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
    fvx = fvx-factor4*tangentversorup(1)-Lit*factor5*vcurvature(1)
    fvy = fvy-factor4*tangentversorup(2)-Lit*factor5*vcurvature(2)
    fvz = fvz-factor4*tangentversorup(3)-Lit*factor5*vcurvature(3)
  endif
  
  return
  
 end subroutine eom3_ev
 
 subroutine eom3_KV_pos_v_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl, &
       yve,ycf,fxx,fyy,fzz,fvx,fvy,fvz,fev,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional model with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, intent(inout) ::  fvx
  double precision, intent(inout) ::  fvy
  double precision, intent(inout) ::  fvz
  double precision, intent(inout) ::  fev
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision :: curvature
  double precision, dimension(3) :: tangentversorup,tangentversordown
  double precision, dimension(3) :: curvcenter,vcurvature,coulomelec
  double precision :: aLorx,aLory,aLorz
  logical :: lstraight
  
  double precision :: Vtvec(3),Fvet,Kst,factor1,factor2,factor3
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re,vnorm
  
  integer,save :: ij=0
  
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fvx=0.d0
    fvy=0.d0
    fvz=0.d0
    fev=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      coulomelec(1:3)=ycf(ipoint,1:3)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
       cmass)
      
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      Fvet=Fve/(jetms(ipoint)*cmass)
      
      factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
      
      fxx = yvx(ipoint) 
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+coulomelec(1)+ &
       upwall(ipoint,yxx)+aLorx
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+coulomelec(2)+aLory
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+coulomelec(3)+aLorz
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      
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
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
       cmass)
      
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      Fvet=Fve/(jetms(ipoint)*cmass)
      Kst=Ks/(jetms(ipoint)*cmass)
      
      factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
      factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
       (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
      
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
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
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
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec, &
         cmass)
        call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
         cmass)
        
        !compute Reynolds number
        vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+ &
         yvz(ipoint)**2.d0)
        Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
        Fvet=Fve/(jetms(ipoint)*cmass)
        Kst=Ks/(jetms(ipoint)*cmass)
      
        factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
        factor2=Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)
        factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
         (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
        
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
        
        fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
         sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fvx=0.d0
        fvy=0.d0
        fvz=0.d0
        fev=0.d0
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
      fev=0.d0
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
      fvx=0.d0
      fvy=0.d0
      fvz=0.d0
      fev=0.d0
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
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
  call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
   cmass)
  
  !compute Reynolds number
  vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
  Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
  Fvet=Fve/(jetms(ipoint)*cmass)
  Kst=Ks/(jetms(ipoint)*cmass)
  
  factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
  factor2=Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)
  factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
   (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
  
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
  
  fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
   sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
  
  return
  
 end subroutine eom3_KV_pos_v_ev
 
 subroutine eom3_KV_st_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       yax,yay,yaz,fevlocal,fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional model with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, allocatable, dimension (:), intent(in) ::  yax
  double precision, allocatable, dimension (:), intent(in) ::  yay
  double precision, allocatable, dimension (:), intent(in) ::  yaz
  double precision, intent(in) :: fevlocal
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
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re,vnorm
  
  integer,save :: ij=0
  
  
! special cases
  
  if(jetfr(ipoint))then
    fst=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call project_beadacctangetversor(ipoint,yax,yay,yaz,beadaccup, &
       tangentversorup)
      
      call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)
      
    else
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      call project_beadacctangetversor(ipoint,yax,yay,yaz,beadaccup, &
       tangentversorup)
      
      call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)
        
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
        
        call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)
        
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
  
  call kv_ev_stress_rate(cp,ratmu,ratg,yve(ipoint),yvl(ipoint), &
     fevlocal,yst(ipoint),beadvelup/beadlenup, &
     beadaccup/beadlenup,fst)
  
  
  return
  
 end subroutine eom3_KV_st_ev
 
 subroutine eom4_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &
       fxx,fyy,fzz,fst,fvx,fvy,fvz,fev,timesub,k,fstocvx,fstocvy, &
       fstocvz) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives 
!     of the system for the three dimensional stochastic model 
!     with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, intent(inout) ::  fev
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
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re,vnorm
  
  integer,save :: ij=0
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fst=0.d0
    fvx=0.d0
    fvy=0.d0
    fvz=0.d0
    fev=0.d0
    fstocvx=0.d0
    fstocvy=0.d0
    fstocvz=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
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
      call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz,cmass)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
       cmass)
      
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      Fvet=Fve/(jetms(ipoint)*cmass)
      attt=att/(jetms(ipoint)*cmass)
      
      if(yst(ipoint)>0.d0)then
        factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
      else
        factor1=0.d0
      endif
      factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
  
      fxx = yvx(ipoint) 
      fst = (1.d0/rattao)*(yieldstress+ &
       consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+coulomelec(1)- &
       factor4*tangentversorup(1)+upwall(ipoint,yxx)+aLorx- &
       noisefric*yvx(ipoint) 
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+coulomelec(2)- &
       factor4*tangentversorup(2)+aLory-noisefric*yvy(ipoint) 
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+coulomelec(3)- &
       factor4*tangentversorup(3)+aLorz-noisefric*yvz(ipoint) 
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      
       
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
      call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz,cmass)
      call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
      call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
       cmass)
      
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      Fvet=Fve/(jetms(ipoint)*cmass)
      Kst=Ks/(jetms(ipoint)*cmass)
      attt=att/(jetms(ipoint)*cmass)
      Lit=Li/(jetms(ipoint)*cmass)
      
      if(yst(ipoint)>0.d0)then
        factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
      else
        factor1=0.d0
      endif
      factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
       (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
      factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
      factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
      
      fxx = yvx(ipoint) 
      fst = (1.d0/rattao)*(yieldstress+ &
       consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
      fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
       Kst*curvature*factor3*vcurvature(1)+ &
       coulomelec(1)-factor4*tangentversorup(1)- &
       Lit*factor5*vcurvature(1)+upwall(ipoint,yxx)+aLorx- &
       noisefric*yvx(ipoint) 
      
      fyy = yvy(ipoint) 
      fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
       Kst*curvature*factor3*vcurvature(2)+ &
       coulomelec(2)-factor4*tangentversorup(2)- &
       Lit*factor5*vcurvature(2)+aLory-noisefric*yvy(ipoint)
      
      fzz = yvz(ipoint) 
      fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
       Kst*curvature*factor3*vcurvature(3)+ &
       coulomelec(3)-factor4*tangentversorup(3)- &
       Lit*factor5*vcurvature(3)+aLorz-noisefric*yvz(ipoint)
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
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
        call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz,cmass)
        call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec, &
         cmass)
        call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz, &
         cmass)
        
        !compute Reynolds number
        vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+ &
         yvz(ipoint)**2.d0)
        Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
        Fvet=Fve/(jetms(ipoint)*cmass)
        Kst=Ks/(jetms(ipoint)*cmass)
        attt=att/(jetms(ipoint)*cmass)
        Lit=Li/(jetms(ipoint)*cmass)
        
        if(yst(ipoint)>0.d0)then
          factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
        else
          factor1=0.d0
        endif
        if(yst(ipoint-1)>0.d0)then
          factor2=Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)
        else
          factor2=0.d0
        endif
        factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
         (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
        factor4=attt*(dabs(beadlenup)**0.905d0)* &
         (dabs(veltangent)**1.19d0)
        factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
        
        fxx = yvx(ipoint) 
        fst = (1.d0/rattao)*(yieldstress+ &
         consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
        fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
         factor2*tangentversordown(1)+ &
         Kst*curvature*factor3*vcurvature(1)+coulomelec(1)- &
         factor4*tangentversorup(1)- &
         Lit*factor5*vcurvature(1)+upwall(ipoint,yxx)+aLorx- &
         noisefric*yvx(ipoint) 
        
        fyy = yvy(ipoint) 
        fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
         factor2*tangentversordown(2)+ &
         Kst*curvature*factor3*vcurvature(2)+coulomelec(2)- &
         factor4*tangentversorup(2)- &
         Lit*factor5*vcurvature(2)+aLory-noisefric*yvy(ipoint) 
        
        fzz = yvz(ipoint) 
        fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
         factor2*tangentversordown(3)+ &
         Kst*curvature*factor3*vcurvature(3)+coulomelec(3)- &
         factor4*tangentversorup(3)- &
         Lit*factor5*vcurvature(3)+aLorz-noisefric*yvz(ipoint) 
        
        fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
         sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fst=0.d0
        fvx=0.d0
        fvy=0.d0
        fvz=0.d0
        fev=0.d0
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
      fev=0.d0
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
      fev=0.d0
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
  call compute_stocforce_3d(ipoint,fstocvx,fstocvy,fstocvz,cmass)
  call driver_electric_field(ipoint,timesub,yxx,yyy,yzz,Vtvec,cmass)
  call compute_lorentz_acc(ipoint,yvx,yvy,yvz,aLorx,aLory,aLorz,cmass)
  
  !compute Reynolds number
  vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
  Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
  Fvet=Fve/(jetms(ipoint)*cmass)
  Kst=Ks/(jetms(ipoint)*cmass)
  attt=att/(jetms(ipoint)*cmass)
  Lit=Li/(jetms(ipoint)*cmass)
  
  if(yst(ipoint)>0.d0)then
    factor1=Fvet*yve(ipoint)*(yst(ipoint)/beadlenup)
  else
    factor1=0.d0
  endif
  if(yst(ipoint-1)>0.d0)then
    factor2=Fvet*yve(ipoint-1)*(yst(ipoint-1)/beadlendown)
  else
    factor2=0.d0
  endif
  factor3=0.25d0*((dsqrt(yve(ipoint))/dsqrt(beadlenup))+ &
   (dsqrt(yve(ipoint-1))/dsqrt(beadlendown)))**2.d0
  factor4=attt*(dabs(beadlenup)**0.905d0)*(dabs(veltangent)**1.19d0)
  factor5=factor3*beadlenup*curvature*(veltangent**2.d0)
  
  fxx = yvx(ipoint) 
  fst = (1.d0/rattao)*(yieldstress+ &
   consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
  fvx = Gr+Vtvec(1)-factor1*tangentversorup(1)+ &
   factor2*tangentversordown(1)+ &
   Kst*curvature*factor3*vcurvature(1)+coulomelec(1)- &
   factor4*tangentversorup(1)-Lit*factor5*vcurvature(1)+ &
   upwall(ipoint,yxx)+aLorx-noisefric*yvx(ipoint) 
  
  fyy = yvy(ipoint) 
  fvy = Vtvec(2)-factor1*tangentversorup(2)+ &
   factor2*tangentversordown(2)+ &
   Kst*curvature*factor3*vcurvature(2)+coulomelec(2)- &
   factor4*tangentversorup(2)-Lit*factor5*vcurvature(2)+aLory- &
   noisefric*yvy(ipoint) 
  
  fzz = yvz(ipoint) 
  fvz = Vtvec(3)-factor1*tangentversorup(3)+ &
   factor2*tangentversordown(3)+ &
   Kst*curvature*factor3*vcurvature(3)+coulomelec(3)- &
   factor4*tangentversorup(3)-Lit*factor5*vcurvature(3)+aLorz- &
   noisefric*yvz(ipoint) 
  
  fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
   sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
  
  return
  
 end subroutine eom4_ev
 
 subroutine eom4_pos_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf,&
       fxx,fyy,fzz,fev,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing only the position first 
!     derivatives of the system for the three dimensional stochastic 
!     model with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  double precision, intent(inout) ::  fev
  double precision, intent(in) :: timesub
  integer, intent(in) :: k
  
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re,vnorm
  double precision :: beadlendown,beadlenup
  
! special cases
  
  if(jetfr(ipoint))then
    fxx=0.d0
    fyy=0.d0
    fzz=0.d0
    fev=0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      fxx = yvx(ipoint)
      fyy = yvy(ipoint) 
      fzz = yvz(ipoint)
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
    else
      call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
      !compute Reynolds number
      vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
      Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
      fxx = yvx(ipoint)
      fyy = yvy(ipoint) 
      fzz = yvz(ipoint) 
      
      fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
       sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
    endif
    return
  endif
  
  if(ipoint==npjet-1)then
    if(ipoint>0)then
      if(linserted)then
        call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
        !compute Reynolds number
        vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+ &
         yvz(ipoint)**2.d0)
        Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
        fxx = yvx(ipoint)
        fyy = yvy(ipoint) 
        fzz = yvz(ipoint) 
        
        fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
         sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
      else
        fxx=0.d0
        fyy=0.d0
        fzz=0.d0
        fev=0.d0
      endif
    endif
    return
  endif
  
  if(ipoint==npjet)then
    if(liniperturb)then
      fxx=0.d0
      fyy = -1.d0*pfreq*yzz(ipoint) 
      fzz = pfreq*yyy(ipoint)
      fev=0.d0
    else
      fxx=0.d0
      fyy=0.d0
      fzz=0.d0
      fev=0.d0
    endif
    return
  endif
  
  
! ordinary case
  call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
!compute Reynolds number
  vnorm=dsqrt(yvx(ipoint)**2.d0+yvy(ipoint)**2.d0+yvz(ipoint)**2.d0)
  Re=(2.d0*dsqrt(yve(ipoint)/(Pi*beadlenup))*vnorm)/evairv
  fxx = yvx(ipoint)
  fyy = yvy(ipoint) 
  fzz = yvz(ipoint) 
  
  fev = -evmasscoeff*0.495d0*(Re**(1.d0/3.d0))* &
   sqrevsc*(evcsvapour-evumidity)*Pi*beadlenup 
  
  
  return
  
 end subroutine eom4_pos_ev
 
  subroutine eom4_stress_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl, &
       yve,ycf,fst,timesub,k) 
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing only the stress first 
!     derivative of the system for the three dimensional stochastic 
!     model with evaporation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
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
  
  double precision :: beadlendown,beadlenup,beadvelup
  double precision, dimension(3) :: tangentversorup
  double precision :: newtao,ratmu,rattao,cmass,cp,cs,ratg,Re,vnorm
  
  
! special cases
  
  if(jetfr(ipoint))then
    fst = 0.d0
    return
  endif
  
  !mass fraction of actual polymer
  cp=cp0*yvl(ipoint)/yve(ipoint)
  !mass fraction of actual solvent
  cs=1.d0-cp
  !ratio between corrected for evaporation tao and old tao
  rattao=(cp/cp0)**tev
  !ratio between corrected for evaporation mu and old mu
  ratmu=10.d0**(Bev*((cp**mev)-(cp0**mev)))
  !correction factor for the evaporated mass
  cmass=yve(ipoint)/yvl(ipoint)
  !ratio between corrected for evaporation G and old G
  ratg=ratmu/rattao
  
  if(ipoint==inpjet)then
    if(ipoint==0)then
      call compute_geometry_init(ipoint,yxx,yyy,yzz,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      
      fst = (1.d0/rattao)*(yieldstress+ &
       consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
      
    else
      call compute_geometry(ipoint,yxx,yyy,yzz,beadlendown,beadlenup)
      call compute_tangetversor(ipoint,yxx,yyy,yzz,tangentversorup, &
       beadlenup)
      call project_beadveltangetversor(ipoint,yvx,yvy,yvz,beadvelup, &
       tangentversorup)
      
      fst = (1.d0/rattao)*(yieldstress+ &
       consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
      
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
       
        fst = (1.d0/rattao)*(yieldstress+ &
         consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
        
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
  
  fst = (1.d0/rattao)*(yieldstress+ &
   consistency*ratmu*(beadvelup/beadlenup)**findex-yst(ipoint))
  
  
  return
  
 end subroutine eom4_stress_ev


 subroutine kv_ev_stress_rate(cp,ratmu,ratg,yve,yvl,fevlocal,stress, &
                              strainrate,strainacc,fst)

!***********************************************************************
! Product-rule Kelvin-Voigt stress rate with concentration-dependent
! viscosity and elastic modulus.  All quantities are in the standard
! JETSPIN nondimensionalization.  The historical JETSPIN Kelvin-Voigt
! kinematics is retained: strainrate=(1/l) dl/dt and the acceleration
! contribution is represented by strainacc=(1/l) dv_parallel/dt.
!***********************************************************************

  implicit none
  double precision, intent(in) :: cp,ratmu,ratg,yve,yvl,fevlocal
  double precision, intent(in) :: stress,strainrate,strainacc
  double precision, intent(out) :: fst
  double precision :: dcpdt,dratmu,dratg,strain

  dcpdt=0.d0
  if(yve>0.d0 .and. yvl>0.d0)then
    if((yve/yvl)>evlim*(1.d0+1.d-12))then
      dcpdt=-cp*fevlocal/yve
    endif
  endif

  dratmu=ratmu*dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))*dcpdt
  dratg=ratg*(dlog(10.d0)*Bev*mev*(cp**(mev-1.d0))-tev/cp)*dcpdt

  strain=(stress-ratmu*strainrate)/ratg
  fst=ratg*strainrate+ratmu*strainacc+dratg*strain+dratmu*strainrate

  return
 end subroutine kv_ev_stress_rate

 end module eom_ev_mod


