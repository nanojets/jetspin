 
 module dynamic_refinement_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which are dealing
!     the dynamic refinement of the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!*********************************************************************** 
 
 use version_mod, only : mystart,myend,memyend,idrank,mxrank, &
                   sum_world_darr,or_world_larr,bcast_world_darr, &
                   sum_world_iarr
 use error_mod,   only : error,warning
 use utility_mod, only : Pi,buffservice,allocate_array_buffservice
 use nanojet_mod, only : resolution,inpjet,npjet,jetpt,jetxx,jetyy, &
                  jetzz,jetvx,jetvy,jetvz,jetst,jetms,jetch,jetvl, &
                  jetcr,linserted,linserting,systype,lengthpath, &
                  doallocate,doreorder,mxnpjet,incnpjet,lengthscale, &
                  ivolume,jetbd,typemass,massratio,imassa,jetfr,h, &
                  lenthresholdbead,ltagbeads
 use fit_mod,     only : jetptc,allocate_arrayspline,create_spline, &
                   allocate_array_jetptc,driver_fit_spline, &
                   looking_indexes_2,looking_indexes_4, &
                   cubic_interpolation,fit_akima,allocate_arrayakima, &
                   jetbdc
 use support_functions_mod, only : compute_length_path,compute_crosssec
 
 implicit none
 
 private
 
 integer, save :: nfitspline
 integer, save :: nfitlinear
 integer, save :: nfitting
 integer, save :: oldlowerbound
 integer, save :: oldlowerbuff
 integer, save :: newlowerbound
 integer, save :: newlowerbuff
 integer, save :: newupperbound
 
 logical, public, save :: lrefinement=.false.
 logical, public, save :: lrefinementthreshold=.false.
 logical, public, save :: lrefinementevery=.false.
 logical, public, save :: lrefinementstart=.false.
 
 logical, public, save :: lrefbeadstart=.false.
 
 integer, public, save :: irefbeadstart=0
 integer, public, save :: irefinementstart=0
 integer, public, save :: irefinementevery=0
 integer, public, save :: irefinementdone=0
 
 integer, public, save :: refinementcons=0
 
 double precision, public, save :: refinementthreshold=0.d0
 
 double precision, public, save :: refbeadstartfit=0.d0
 
 double precision, allocatable, save :: buff(:)
 integer, save :: nbuff=0
 
 logical, save :: mydoallocate
 
 integer, save :: nmyindex=0
 integer, allocatable, save :: myindex(:,:)
 
 integer, save :: nkeepvol=0
 double precision, dimension(:), allocatable, save :: keepvol
 
 integer, save :: inpjetbackup,npjetbackup
 integer, save :: nbackup=0
 integer, save :: njetlthr=0
 double precision, save :: timebackup
 logical, dimension(:), allocatable, save :: jetlthr
 double precision, dimension(:), allocatable, save :: jetxxbak, &
  jetyybak,jetzzbak,jetvxbak,jetvybak,jetvzbak,jetstbak,jetmsbak, &
  jetchbak,jetvlbak
 
 
 public :: driver_dynamic_refinement
 public :: set_refinement_threshold
 
 contains
 
 subroutine driver_dynamic_refinement(k,dorefinment)

!***********************************************************************
!     
!     JETSPIN subroutine for driving the dynamic refinement
!     of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************

implicit none
  
  integer, intent(in) :: k
  logical, intent(inout) :: dorefinment
  
  if(.not.lrefinement)return
  if(k<irefinementstart)return
  
! check if the dynamic refinement is necessary
  call check_dynamic_refinement_akima(k,dorefinment)
  
  if(.not.dorefinment)return
  
! apply the dynamic refinement if requested
  call dynamic_refinement_akima(k,dorefinment)
  
  return
  
 end subroutine driver_dynamic_refinement
 
 subroutine set_refinement_threshold()
 
!***********************************************************************
!     
!     JETSPIN subroutine for setting the threshold of the dynamics 
!     refinement of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  if(lrefinement)then
    if(.not.lrefinementthreshold)then
      refinementthreshold=10.d0*resolution
      refbeadstartfit=10.d0*resolution
      lenthresholdbead=refinementthreshold
      lrefbeadstart=.true.
      lrefinementthreshold=.true.
      call warning(66,refinementthreshold)
    else
      if(refinementthreshold<(2.d0*resolution))call error(16)
    endif
  endif
  
  return
  
 end subroutine set_refinement_threshold
 
 subroutine check_dynamic_refinement_akima(nstep,dorefinment)
 
!***********************************************************************
!     
!     JETSPIN subroutine for checking the dynamic refinement
!     of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  logical, intent(out) :: dorefinment
 
  integer, save :: icounter=0
  double precision :: tempmod0
  logical :: ldorefinement
  
  dorefinment=.false.
  
  if(.not.lrefinement)return
  if(linserting .and. (npjet-inpjet+1)<10)return
  
  icounter=icounter+1
  if(icounter<irefinementevery)return
  
  call test_length_threshold_akima(nstep,jetxx,jetyy,jetzz, &
   refinementthreshold,ldorefinement,jetpt,lengthpath)
  
  if(.not.ldorefinement)return
  
  if(linserting)then
    if(.not.linserted)then
      select case(systype)
      case(1)
        tempmod0 = jetxx(npjet-2)-jetxx(npjet)
      case default
        tempmod0 = dsqrt((jetxx(npjet-2)-jetxx(npjet))**2.d0+ &
         (jetyy(npjet-2)-jetyy(npjet))**2.d0+ &
         (jetzz(npjet-2)-jetzz(npjet))**2.d0)
      end select
      nfitspline=nint((lengthpath-tempmod0)/(refinementthreshold))
      nfitting=nfitspline+3
      if(nfitting<=nint(dble(npjet-inpjet+1)+1.d0))return
    else
      select case(systype)
      case(1)
        tempmod0 = jetxx(npjet-1)-jetxx(npjet)
      case default
        tempmod0 = dsqrt((jetxx(npjet-1)-jetxx(npjet))**2.d0+ &
         (jetyy(npjet-1)-jetyy(npjet))**2.d0+ &
         (jetzz(npjet-1)-jetzz(npjet))**2.d0)
      end select
      nfitspline=nint((lengthpath-tempmod0)/(refinementthreshold))
      nfitting=nfitspline+2
      if(nfitting<=nint(dble(npjet-inpjet+1)+1.d0))return
    endif
  else
    nfitspline=nint((lengthpath)/(refinementthreshold))
    nfitting=nfitspline
    if(nfitting<=nint(dble(npjet-inpjet+1)+1.d0))return
  endif
  
  dorefinment=.true.
  icounter=0
  irefinementdone=irefinementdone+1
  
  return
  
 end subroutine check_dynamic_refinement_akima
 
 subroutine dynamic_refinement_akima(nstep,dorefinment)
 
!***********************************************************************
!     
!     JETSPIN subroutine for performing the dynamic refinement
!     of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  logical, intent(in) :: dorefinment
  
  integer :: i,j,jptinit,jptend,ipoint,totjptend,nmassbd,irefbeadstop
  integer :: totjptend0,istart
  integer, parameter :: strategytype=0
  integer, allocatable :: massbd(:,:),massbdpoint(:)
  double precision, allocatable :: massbddist(:),massbdgap(:)
  double precision :: tempmod0,tempmod1,tempmod2,voltot,newvoltot
  double precision :: mastot,newmastot
  double precision :: chrtot,newchrtot
  
  
  integer, save :: icounter=0
  
  if(.not.dorefinment)return
  
  mydoallocate=.false.
      
  irefbeadstop=-1
  do ipoint=inpjet,irefbeadstart
    if(jetxx(ipoint)>=h)irefbeadstop=max(irefbeadstop,ipoint)
  enddo
    
  if(irefbeadstop==-1)then
    istart=inpjet
  else
    istart=irefbeadstop
  endif
  
  nmassbd=0
  do ipoint=istart,irefbeadstart
    if(jetbd(ipoint))nmassbd=nmassbd+1
  enddo
  
  allocate(massbd(0:nmassbd,2),massbdpoint(0:nmassbd))
  allocate(massbddist(0:nmassbd),massbdgap(0:nmassbd))
    
  i=0
  massbd(0,1)=istart
  do ipoint=massbd(0,1),irefbeadstart
    if(jetbd(ipoint))then
      i=i+1
      massbd(i-1,2)=ipoint
      massbd(i,1)=ipoint
    endif
  enddo
  massbd(i,2)=irefbeadstart
  
  nfitspline=0
  do i=0,nmassbd
    tempmod0=(jetpt(massbd(i,2))-jetpt(massbd(i,1)))
    massbddist(i)=tempmod0
    massbdpoint(i)=ceiling((tempmod0*lengthpath)/(refinementthreshold))
    massbdgap(i)=(tempmod0/dble(massbdpoint(i)))
    nfitspline=massbdpoint(i)+nfitspline
  enddo
  
  if(irefbeadstop<inpjet)then 
    nfitspline=nfitspline+(irefbeadstop-inpjet+1)+2
    nfitting=nfitspline+(npjet-irefbeadstart)
  else
    if(irefbeadstop==inpjet)then
      nfitspline=nfitspline+(irefbeadstop-inpjet+1)
      nfitting=nfitspline+(npjet-irefbeadstart)
    else
      nfitspline=nfitspline+(irefbeadstop-inpjet+1)
      nfitting=nfitspline+(npjet-irefbeadstart)
    endif
  endif
  
  if(nfitting<=(npjet-inpjet+1))return
  
  call allocate_arrayakima()
  
  call convert_to_density(jetms,jetch,jetvl)
  
  voltot=0.d0
  mastot=0.d0
  chrtot=0.d0
  do i=inpjet,irefbeadstart
    voltot=voltot+jetvl(i)
    mastot=mastot+jetms(i)
    chrtot=chrtot+jetch(i)
  enddo
  call allocate_array_keepvol(npjet-irefbeadstart)
  i=0
  do ipoint=irefbeadstart+1,npjet
    i=i+1
    keepvol(i)=jetvl(ipoint)
  enddo
  
  call compute_crosssec(jetxx,jetyy,jetzz,jetvl,jetcr)
  
  call define_akima_bounds()
  call allocate_array_jetptc()
  j=0
  jptinit=newlowerbuff+1
  jptend=jptinit+nfitspline
  
  jetbdc(:)=.false.
  
  j=jptinit
  jetptc(j)=0.d0
  if(jetbd(inpjet))jetbdc(j)=.true.
  do ipoint=inpjet+1,irefbeadstop
    j=j+1
    jetptc(j)=jetpt(ipoint)
    if(jetbd(ipoint))jetbdc(j)=.true.
  enddo
  do i=0,nmassbd
    tempmod0=massbdgap(i)
    do ipoint=1,massbdpoint(i)
      j=j+1
      jetptc(j)=jetptc(j-1)+tempmod0
      if(ipoint==massbdpoint(i))then
        if(i/=nmassbd)then
          jetbdc(j)=.true.
          jetptc(j)=jetpt(massbd(i,2))
        else
          jetptc(j)=jetpt(irefbeadstart)
        endif         
      endif
    enddo
  enddo
  totjptend0=j
  do ipoint=irefbeadstart+1,npjet
    j=j+1
    jetptc(j) = jetpt(ipoint)
    if(jetbd(ipoint))then
      jetbdc(j)=.true.     
    endif
  enddo
  totjptend=j
    
  
  if(mydoallocate)then
    deallocate(jetbd)
    allocate(jetbd(0:mxnpjet))
  endif
  jetbd(0:mxnpjet)=.false.
  do ipoint=jptinit,totjptend
    jetbd(ipoint)=jetbdc(ipoint)
  enddo
  
  deallocate(massbd,massbdpoint)
  deallocate(massbddist,massbdgap)
  
   
  call fit_jet_akima(jptinit,totjptend)
  
    
  do i=jptinit,jptend
    tempmod0=lengthpath*(jetptc(i+1)-jetptc(i))
    jetvl(i)=tempmod0*Pi*(jetcr(i))**2.d0
  enddo
  
  i=0
  do ipoint=jptend+1,totjptend
    i=i+1
    jetvl(ipoint) = keepvol(i)
  enddo
  
  newvoltot=0.d0
  do i=jptinit,jptend
    newvoltot=newvoltot+jetvl(i)
  enddo
  
  jetvl(jptinit:jptend)=jetvl(jptinit:jptend)*voltot/newvoltot
  
  newmastot=0.d0
  newchrtot=0.d0
  do i=jptinit,jptend
    newmastot=newmastot+jetms(i)
    newchrtot=newchrtot+jetch(i)
  enddo
  
  if(refinementcons==0)then
    jetms(jptinit:jptend)=jetms(jptinit:jptend)*mastot/newmastot
    jetch(jptinit:jptend)=jetch(jptinit:jptend)*chrtot/newchrtot
  endif
  
  call convert_from_density(jetms,jetch,jetvl)
  
  doallocate=(doallocate .or. mydoallocate)
  
  return
  
 end subroutine dynamic_refinement_akima
 
 subroutine test_length_threshold_akima(nstep,yxx,yyy,yzz,dthreshold, &
  lthreshold,ypt,lengthpathsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for testing if the mutual distance between
!     any two nanofiber beads is beyond a given threshold
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  double precision, allocatable, dimension(:), intent(in) :: yxx
  double precision, allocatable, dimension(:), intent(in) :: yyy
  double precision, allocatable, dimension(:), intent(in) :: yzz
  double precision, intent(in) :: dthreshold
  
  logical, intent(inout) :: lthreshold
  double precision, allocatable, dimension(:), intent(inout) :: ypt
  double precision, intent(inout) :: lengthpathsub
  
  integer :: ipoint
  double precision :: tempmod0,tempmod,tempmod1(1)
  logical, dimension(1) :: logicaltemp
  
  call compute_length_path(yxx,yyy,yzz, &
   ypt,lengthpathsub)
  
  logicaltemp(1)=.false.
  
  irefbeadstart=0
  if(linserting)then
    if(.not.linserted)then
      do ipoint=inpjet,npjet-3
        tempmod0 = (ypt(ipoint+1)-ypt(ipoint))*lengthpathsub
        if(tempmod0 > refbeadstartfit )then
          logicaltemp(1)=.true.
          irefbeadstart=ipoint
        endif
      enddo
    else
      do ipoint=inpjet,npjet-2
        tempmod0 = (ypt(ipoint+1)-ypt(ipoint))*lengthpathsub
        if(tempmod0 > refbeadstartfit )then
          logicaltemp(1)=.true.
          irefbeadstart=ipoint
        endif
      enddo 
    endif
  else
    do ipoint=inpjet,npjet-1
      tempmod0 = (ypt(ipoint+1)-ypt(ipoint))*lengthpathsub
      if(tempmod0 > refbeadstartfit )then
        logicaltemp(1)=.true.
        irefbeadstart=ipoint
      endif
    enddo
  endif
  
  lthreshold=logicaltemp(1)
  
  return
  
 end subroutine test_length_threshold_akima
 
 subroutine define_akima_bounds() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining the new bounds of arrays which 
!     describe the nanojet if the spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer :: ncutoffsub
  
  mydoallocate=.false.
  ncutoffsub=1
  
  if(inpjet==0)then
    oldlowerbound=0
    oldlowerbuff=inpjet-1
    
    newlowerbound=0
    newlowerbuff=oldlowerbuff-oldlowerbound
    
    mydoallocate=.false.
    newupperbound=nfitting
    if(newupperbound>mxnpjet)then
      mxnpjet=newupperbound+incnpjet
      mydoallocate=.true.
    endif
    
  elseif(inpjet>0)then
    oldlowerbound=max(0,inpjet-ncutoffsub)
    oldlowerbuff=inpjet-1
    
    newlowerbound=0
    newlowerbuff=oldlowerbuff-oldlowerbound
    
    newupperbound=nfitting+newlowerbuff+1
    if(newupperbound>mxnpjet)then
      mxnpjet=newupperbound+incnpjet
      mydoallocate=.true.
    endif
    
  endif
  
  return
  
 end subroutine define_akima_bounds
 
 subroutine fit_jet_akima(jptinit,jptend) 
 
!***********************************************************************
!     
!     JETSPIN subroutine for managing the cubic spline interpolation
!     if any two nanofiber beads are beyond a given threshold
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: jptinit,jptend
  
  
  integer :: i,ipoint
  double precision :: tempmod0
  
  call allocate_array_buffservice(newupperbound)
  
  doreorder=.true.
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetxx(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetxx,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetxx)
    allocate(jetxx(0:mxnpjet))
  endif
  jetxx(:)=0.d0
  jetxx(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetyy(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetyy,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetyy)
    allocate(jetyy(0:mxnpjet))
  endif
  jetyy(:)=0.d0
  jetyy(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetzz(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetzz,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetzz)
    allocate(jetzz(0:mxnpjet))
  endif
  jetzz(:)=0.d0
  jetzz(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetvx(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetvx,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetvx)
    allocate(jetvx(0:mxnpjet))
  endif
  jetvx(:)=0.d0
  jetvx(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
    
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetvy(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetvy,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetvy)
    allocate(jetvy(0:mxnpjet))
  endif
  jetvy(:)=0.d0
  jetvy(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetvz(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetvz,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetvz)
    allocate(jetvz(0:mxnpjet))
  endif
  jetvz(:)=0.d0
  jetvz(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetst(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetst,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetst)
    allocate(jetst(0:mxnpjet))
  endif
  jetst(:)=0.d0
  jetst(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetms(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetms,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetms)
    allocate(jetms(0:mxnpjet))
  endif
  jetms(:)=0.d0
  jetms(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetch(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetch,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetch)
    allocate(jetch(0:mxnpjet))
  endif
  jetch(:)=0.d0
  jetch(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetcr(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetcr,jetptc,buffservice)
  if(mydoallocate)then
    deallocate(jetcr)
    allocate(jetcr(0:mxnpjet))
  endif
  jetcr(:)=0.d0
  jetcr(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  if(mydoallocate)then
    deallocate(jetvl)
    allocate(jetvl(0:mxnpjet))
  endif
  jetvl(:)=0.d0
  
  if(mydoallocate)then
    deallocate(jetpt)
    allocate(jetpt(0:mxnpjet))
  endif
  jetpt(:)=0.d0
  
  
  inpjet=jptinit
  npjet=jptend
  
  do ipoint=inpjet,npjet
    jetpt(ipoint)=jetptc(ipoint)
  enddo
    
  if(mydoallocate)then
    deallocate(jetfr)
    allocate(jetfr(0:mxnpjet))
  endif
  jetfr(:)=.false.
  do ipoint=inpjet,npjet
    if(jetxx(ipoint)>=h)then
      jetfr(ipoint)=.true.
      jetxx(ipoint)=h
    endif
  enddo
  jetxx(npjet)=0.d0
  
 ! call compute_length_path(jetxx,jetyy,jetzz,jetpt,lengthpath)
  
  return
  
 end subroutine fit_jet_akima
 
 subroutine convert_to_density(jms,jch,jvl)

!***********************************************************************
!     
!     JETSPIN subroutine for converting the mass and charge data
!     to their respective densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) :: jms,jch
  double precision, allocatable, dimension (:), intent(inout) :: jvl
  
  integer :: ipoint
  double precision :: extramass,extravol
  
  
  if(typemass==3  .or. ltagbeads)then
    extramass=imassa*(massratio-1.d0)*ivolume
    extravol=(massratio-1.d0)*ivolume
    do ipoint=inpjet,npjet
      if(jetbd(ipoint))then
        jvl(ipoint)=jvl(ipoint)-extravol
        jms(ipoint)=(jms(ipoint)-extramass)/jvl(ipoint)
        jch(ipoint)=jch(ipoint)/jvl(ipoint)
      else
        jms(ipoint)=jms(ipoint)/jvl(ipoint)
        jch(ipoint)=jch(ipoint)/jvl(ipoint)
      endif
    enddo
  else
    jms(:)=jms(:)/jvl(:)
    jch(:)=jch(:)/jvl(:)
  endif
  
  
  return
  
 end subroutine convert_to_density
 
 subroutine convert_from_density(jms,jch,jvl)

!***********************************************************************
!     
!     JETSPIN subroutine for converting the mass and charge data
!     from their respective densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) :: jms,jch
  double precision, allocatable, dimension (:), intent(inout) :: jvl
  
  integer :: ipoint
  double precision :: extramass,extravol
  
  if(typemass==3  .or. ltagbeads)then
    extramass=imassa*(massratio-1.d0)*ivolume
    extravol=(massratio-1.d0)*ivolume
    do ipoint=inpjet,npjet
      if(jetbd(ipoint))then
        jms(ipoint)=jms(ipoint)*jvl(ipoint)+extramass
        jch(ipoint)=jch(ipoint)*jvl(ipoint)
        jvl(ipoint)=jvl(ipoint)+extravol
      else
        jms(ipoint)=jms(ipoint)*jvl(ipoint)
        jch(ipoint)=jch(ipoint)*jvl(ipoint)
      endif
    enddo
  else
    jms(:)=jms(:)*jvl(:)
    jch(:)=jch(:)*jvl(:)
  endif
  
  return
  
 end subroutine convert_from_density
 
 subroutine fit_jetpt_spline() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the array jetpt which 
!     parametrizes the nanojet if the spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  if(mydoallocate)then
    deallocate(jetpt)
    allocate(jetpt(0:mxnpjet))
  endif
  jetpt(0:npjet)=jetptc(0:npjet)
  jetpt(npjet+1:mxnpjet)=0.d0
  
  
  return
  
 end subroutine fit_jetpt_spline
 
 subroutine store_backup(timesub)

!***********************************************************************
!     
!     JETSPIN subroutine for storing all the data on the service 
!     backup arrays if it is requested
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(inout) :: timesub
  
  inpjetbackup=inpjet
  npjetbackup=npjet
  timebackup=timesub
  
  call allocate_backup(timesub)
  
  if(idrank==0)then
  
  jetxxbak(0:npjet)=jetxx(0:npjet)
  jetyybak(0:npjet)=jetyy(0:npjet)
  jetzzbak(0:npjet)=jetzz(0:npjet)
  jetvxbak(0:npjet)=jetvx(0:npjet)
  jetvybak(0:npjet)=jetvy(0:npjet)
  jetvzbak(0:npjet)=jetvz(0:npjet)
  jetstbak(0:npjet)=jetst(0:npjet)
  jetmsbak(0:npjet)=jetms(0:npjet)
  jetchbak(0:npjet)=jetch(0:npjet)
  jetvlbak(0:npjet)=jetvl(0:npjet)
  
  endif
  
  return
  
 end subroutine store_backup
 
 subroutine restore_backup(timesub)

!***********************************************************************
!     
!     JETSPIN subroutine for restoring all the data from the service 
!     backup arrays if it is requested
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(inout) :: timesub
  
  inpjet=inpjetbackup
  npjet=npjetbackup
  timesub=timebackup
  
  if(idrank==0)then
  
  jetxx(0:npjet)=jetxxbak(0:npjet)
  jetyy(0:npjet)=jetyybak(0:npjet)
  jetzz(0:npjet)=jetzzbak(0:npjet)
  jetvx(0:npjet)=jetvxbak(0:npjet)
  jetvy(0:npjet)=jetvybak(0:npjet)
  jetvz(0:npjet)=jetvzbak(0:npjet)
  jetst(0:npjet)=jetstbak(0:npjet)
  jetms(0:npjet)=jetmsbak(0:npjet)
  jetch(0:npjet)=jetchbak(0:npjet)
  jetvl(0:npjet)=jetvlbak(0:npjet)
  
  
  endif
  
  call bcast_world_darr(jetxx,npjet+1)
  call bcast_world_darr(jetyy,npjet+1)
  call bcast_world_darr(jetzz,npjet+1)
  call bcast_world_darr(jetvx,npjet+1)
  call bcast_world_darr(jetvy,npjet+1)
  call bcast_world_darr(jetvz,npjet+1)
  call bcast_world_darr(jetst,npjet+1)
  call bcast_world_darr(jetms,npjet+1)
  call bcast_world_darr(jetch,npjet+1)
  call bcast_world_darr(jetvl,npjet+1)
  
  return
  
 end subroutine restore_backup
 
 subroutine allocate_backup(timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the service backup arrays
!     which are used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(inout) :: timesub
  
  logical, save :: mioallocfirst=.true.
  
  if(idrank==0)then
    if(mxnpjet>nbackup)then
      nbackup=mxnpjet
      if(.not.mioallocfirst)then
        deallocate(jetxxbak,jetyybak,jetzzbak)
        deallocate(jetvxbak,jetvybak,jetvzbak)
        deallocate(jetstbak,jetmsbak,jetchbak)
        deallocate(jetvlbak)
      endif
      allocate(jetxxbak(0:mxnpjet),jetyybak(0:mxnpjet), &
       jetzzbak(0:mxnpjet))
      allocate(jetvxbak(0:mxnpjet),jetvybak(0:mxnpjet), &
       jetvzbak(0:mxnpjet))
      allocate(jetstbak(0:mxnpjet),jetmsbak(0:mxnpjet), &
       jetchbak(0:mxnpjet))
      allocate(jetvlbak(0:mxnpjet))
    endif
  endif
  
  mioallocfirst=.false.
  
  return
  
 end subroutine allocate_backup
 
 subroutine allocate_jetlthr()
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the service backup arrays
!     which are used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, save :: mioallocfirst=.true.
  
  if(mxnpjet>njetlthr)then
    njetlthr=mxnpjet
    if(.not.mioallocfirst)then
      deallocate(jetlthr)
    endif
    allocate(jetlthr(0:mxnpjet))
  endif
  
  mioallocfirst=.false.
  
  return
  
 end subroutine allocate_jetlthr
 
 subroutine allocate_array_myindex(jmiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array myindex 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jmiomax
  
  if(nmyindex/=0)then
    if(jmiomax>nmyindex)then
      deallocate(myindex)
      nmyindex=jmiomax+10
      allocate(myindex(nmyindex,6))
    endif
  else
    nmyindex=jmiomax+10
    allocate(myindex(nmyindex,6))
  endif
  
  return
  
 end subroutine allocate_array_myindex
 
 subroutine allocate_array_keepvol(jmiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array myindex 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jmiomax
  
  if(nkeepvol/=0)then
    if(jmiomax>nkeepvol)then
      deallocate(keepvol)
      nkeepvol=jmiomax+10
      allocate(keepvol(nkeepvol))
    endif
  else
    nkeepvol=jmiomax+10
    allocate(keepvol(nkeepvol))
  endif
  
  return
  
 end subroutine allocate_array_keepvol
 
 end module dynamic_refinement_mod
 
