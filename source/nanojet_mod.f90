
 module nanojet_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which are defining
!     and dealing the quantities describing the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************

 use version_mod, only : idrank,set_mxchunk,bcast_world_d, &
                   bcast_world_l
 use error_mod,   only : error,warning
 use utility_mod, only : Pi,modulvec,cross,dot,gauss,ibuffservice, &
                   allocate_array_ibuffservice,buffservice, &
                   allocate_array_buffservice,lbuffservice, &
                   allocate_array_lbuffservice

 implicit none

 private
 
 integer, public, save :: nvariable
 integer, public, save :: ndimension
 integer, public, save :: inpjet
 integer, public, save :: npjet
 integer, public, save :: mxnpjet
 integer, public, parameter :: incnpjet=100
 integer, public, save :: systype
 integer, public, save :: units
 integer, public, save :: insertmode
 integer, public, save :: myseed=1
 integer, public, save :: naddtrack=0
 integer, public, save :: nremtrack=0
 integer, public, save :: typemass=0
 integer, public, save :: typedragvel=0
 integer, public, save :: nmulstep=0
 integer, save, public :: nmulstepdone=0
 integer, save, public :: nmultisteperror=0
 logical, public, save :: ldragvel=.false.
 logical, public, save :: doallocate
 logical, public, save :: doreorder
 logical, public, save :: lbreakup=.false.
 logical, public, save :: lfieldtype=.false.
 logical, public, save :: lfieldfreq=.false.
 logical, public, save :: lneighlistdo=.true.
 logical, public, save :: lmultisteperror=.false.
 logical, public, save :: lmaxdispl=.false.
 logical, public, save :: ltaoelectr=.false.
 logical, public, save :: lreadrest=.false.
 logical, public, save :: lfirstmass=.true.
 logical, public, save :: lflorentz=.false.
 logical, public, save :: lmagneticfield=.false.
 double precision, public, save :: multisteperror=0.d0
 double precision, public, save :: q=0.d0
 double precision, public, save :: sqrq=0.d0
 double precision, public, save :: v=0.d0
 double precision, public, save :: fve=0.d0
 double precision, public, save :: Hg=0.d0
 double precision, public, save :: Lrg=0.d0
 double precision, public, save :: att=0.d0
 double precision, public, save :: ks=0.d0
 double precision, public, save :: Li=0.d0
 double precision, public, save :: Gr=0.d0
 double precision, public, save :: lengthscale
 double precision, public, save :: imassa!=2.25487812711407d-4 !g
 double precision, public, save :: lencorrmassa=0.0d0
 double precision, public, save :: icharge!=9.37212742146083d-4 !g^0.5 cm^1.5 sec^-1
 double precision, public, save :: mu!=12211.32d0 !g cm^-1 sec^-1
 double precision, public, save :: G!=122132.d0 ! g cm^-1 sec^-2
 double precision, public, save :: surfacet!=700.d0 ! g sec^-2
 double precision, public, save :: icrossec!=1.5d-4 ! cm
 double precision, public, save :: V0!=30.69980124d0 ! cm^0.5 g^0.5 sec^-1
 double precision, public, save :: h!=2.d0 ! cm
 double precision, public, save :: ilength!=3.19d-3 ! cm
 double precision, public, save :: ivelocity!=0.d0 ! cm s^-1
 double precision, public, save :: istress!=0.d0 ! g cm^-1 s^-2
 double precision, public, save :: yieldstress!=0.d0 ! g cm^-1 s^-2
 double precision, public, save :: findex 
 double precision, public, save :: consistency
 double precision, public, save :: resolution!=3.9875d-5 !cm
 double precision, public, save :: tao ! s
 double precision, public, save :: hresolution
 double precision, public, save :: thresolution
 double precision, public, save :: dresolution
 double precision, public, save :: pfreq ! s^-1
 double precision, public, save :: pampl ! cm
 double precision, public, save :: airdragamp(1:3)=0.d0 ! cm^2 s^-2
 double precision, public, save :: gamma(1:3)=0.d0 ! s^-1
 double precision, public, save :: aird=0.d0 
 double precision, public, save :: airv=0.d0 
 double precision, public, save :: velext=0.d0 ! cm s^-1
 double precision, public, save :: current=0.d0 !g^0.5 cm^1.5 sec^-2
 double precision, public, save :: timedeposition=0.d0 !sec
 double precision, public, save :: ultimatestrength=0.d0
 double precision, public, save :: meancharge
 double precision, public, save :: meanmass
 double precision, public, save :: regq
 double precision, public, save :: tanfriction
 double precision, public, save :: massscale
 double precision, public, save :: chargescale
 double precision, public, save :: tstep=0.d0
 double precision, public, save :: xyzrescale=1.d0
 double precision, public, save :: pdbrescale=1.d0
 double precision, public, save :: insvx=0.d0
 double precision, public, save :: insvy=0.d0
 double precision, public, save :: insvz=0.d0
 double precision, public, save :: dcutoff=0.d0
 double precision, public, save :: ivolume=0.d0
 double precision, public, save :: fvere=0.d0
 double precision, public, save :: ksre=0.d0
 double precision, public, save :: lire=0.d0
 double precision, public, save :: lengthpath=0.d0
 double precision, public, save :: kuppot=0.d0
 double precision, public, save :: lenprobmassa=0.d0
 double precision, public, save :: massratio=1.d0
 double precision, public, save :: lenthresholdbead=0.d0
 double precision, public, save :: fieldfreq=0.d0
 double precision, public, save :: maxdispl=0.d0
 double precision, public, save :: taoelectr=0.d0
 double precision, public, save, dimension(3) :: BLor=(/0.d0,0.d0,0.d0/)
 double precision, public, save :: KLor=0.d0
 double precision, public, save :: oldgaussn,corr,radcorr
 
 logical, allocatable, public, save :: jetbd(:)
 logical, allocatable, public, save :: jetbr(:)
 logical, allocatable, public, save :: jetfr(:)
 logical, allocatable, public, save :: jetfm(:)
 integer, allocatable, public, save :: jetlb(:)
 double precision, allocatable, public, save :: jetpt(:)
 double precision, allocatable, public, save :: jetxx(:)
 double precision, allocatable, public, save :: jetyy(:)
 double precision, allocatable, public, save :: jetzz(:)
 double precision, allocatable, public, save :: jetst(:)
 double precision, allocatable, public, save :: jetvx(:)
 double precision, allocatable, public, save :: jetvy(:)
 double precision, allocatable, public, save :: jetvz(:)
 double precision, allocatable, public, save :: jetms(:) !g
 double precision, allocatable, public, save :: jetch(:) !g^0.5 cm^1.5 sec^-1
 double precision, allocatable, public, save :: jetcr(:) ! cm
 double precision, allocatable, public, save :: jetvl(:) ! cm^3
 
 logical, private, parameter :: listresscompensate=.false.
 logical, public, parameter :: ldevelopers=.false.
 logical, public, save :: lnthreads=.false.
 logical, public, save :: lsystype=.false.
 logical, public, save :: lresolution=.false.
 logical, public, save :: lilengt=.false.
 logical, public, save :: limass=.false.
 logical, public, save :: licharge=.false.
 logical, public, save :: licrossec=.false.
 logical, public, save :: listress=.false.
 logical, public, save :: livelocity=.false.
 logical, public, save :: lmu=.false.
 logical, public, save :: lG=.false.
 logical, public, save :: lh=.false.
 logical, public, save :: lV0=.false.
 logical, public, save :: lsurfacet=.false.
 logical, public, save :: lnpjet=.false.
 logical, public, save :: lunits=.false.
 logical, public, save :: lyieldstress=.false.
 logical, public, save :: linserting=.false.
 logical, public, save :: lremove=.false.
 logical, public, save :: linserted=.true. 
 logical, public, save :: liniperturb=.false.
 logical, public, save :: lmultiplestep=.false.
 logical, public, save :: ldcutoff=.false.
 logical, public, save :: lpfreq=.false.
 logical, public, save :: lpampl=.false.
 logical, public, save :: lairdrag=.false.
 logical, public, save :: lairdragamp=.false.
 logical, public, save :: laird=.false.
 logical, public, save :: lairv=.false.
 logical, public, save :: lairvel=.false.
 logical, public, save :: llift=.false.
 logical, public, save :: lmyseed=.false.
 logical, public, save :: lHBfluid=.false.
 logical, public, save :: lKVfluid=.false.
 logical, public, save :: lconsistency=.false.
 logical, public, save :: lfindex=.false.
 logical, public, save :: ltstep=.false.
 logical, public, save :: lgravity=.false.
 logical, public, save :: lmassavariable=.false.
 logical, public, save :: llencorrmassa=.false.
 logical, public, save :: lxyzrescale=.false.
 logical, public, save :: lpdbrescale=.false.
 logical, public, save :: lpdbtagbeads=.false.
 logical, public, save :: ltrackbeads=.false.
 logical, public, save :: lreordertrack=.false.
 logical, public, save :: luppot=.false.
 logical, public, save :: lmirror=.false.
 logical, public, save :: ltagbeads=.false.
 logical, public, save :: lultimatestrength=.false.
 
 
 public :: set_resolution_length
 public :: allocate_jet
 public :: deallocate_jet
 public :: set_initial_jet
 public :: add_jetbead
 public :: reallocate_jet
 public :: fcut
 public :: remove_jetbead
 public :: erase_jetbead
 public :: compute_posnoinserted
 
 
 contains
 
 subroutine set_resolution_length() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for setting the length step in which 
!     the nanofiber is discretised
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  inpjet=0
  
! first strategy to set the length step (resolution):
! if the the parameter resolution is defined in input file the
! length step is simply equal to resolution and the number of beads 
! is computed by a simple division ilength/resolution
! ilength= initial length of the nanofiber
  if(lresolution)then
    npjet=nint(ilength/resolution)
    call warning(21,dble(npjet))
    if(mod(ilength,resolution)/=0.d0)then
      call warning(1,mod(ilength,resolution))
      call warning(2)
      resolution=ilength/dble(npjet)
      call warning(3,resolution)
    endif
  endif
  
! second strategy to set the length step:
! if the the parameter resolution is not defined in input file but the
! the number of beads (npjet) in which the jet is initially discretised 
! is given in input file, the length step is computed by the division
! ilength/npjet
  if(lnpjet)then
    resolution=ilength/dble(npjet)
    call warning(21,dble(npjet))
    call warning(23,resolution)
  endif
  
! if it was not set in input file, set the cutoff of the primary shell 
! for the Coulomb forces which are explicitly computed by direct sum
  if(lmultiplestep)then
    if(.not.ldcutoff)then
       dcutoff=h
       call warning(30,dcutoff)
    endif
  endif
  
! the length step has to be grater than the cross section radius
! of jet otherwise the model loses validity
  if(resolution<icrossec)then
    call warning(62)
    call error(8)
  endif
  
  return
  
 end subroutine set_resolution_length
 
 subroutine allocate_jet(miomax) 
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating arrays which 
!     describe the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in), optional :: miomax
  
  
  
  
  select case(systype)
    case(1:2)
      nvariable = 3
      ndimension = 1
    case default
      nvariable = 3
      ndimension = 3
  end select
  
  if(present(miomax))then
    if(.not. miomax)then
      mxnpjet=incnpjet
    endif
  else
    mxnpjet=incnpjet
  endif
  mxnpjet=max(mxnpjet,npjet)
  
  call set_mxchunk(mxnpjet)
  allocate(jetfr(0:mxnpjet))
  if(lmultiplestep)allocate(jetfm(0:mxnpjet))
  allocate(jetpt(0:mxnpjet))
  allocate(jetxx(0:mxnpjet))
  allocate(jetyy(0:mxnpjet))
  allocate(jetzz(0:mxnpjet))
  allocate(jetst(0:mxnpjet))
  allocate(jetvx(0:mxnpjet))
  allocate(jetvy(0:mxnpjet))
  allocate(jetvz(0:mxnpjet))
  allocate(jetms(0:mxnpjet))
  allocate(jetch(0:mxnpjet))
  allocate(jetcr(0:mxnpjet))
  allocate(jetvl(0:mxnpjet))
  if(ltrackbeads .and. idrank==0)allocate(jetlb(0:mxnpjet))
  if(typemass==3 .or. ltagbeads)allocate(jetbd(0:mxnpjet))
  if(lbreakup)allocate(jetbr(0:mxnpjet))
  
  doallocate=.true.
  
  return
  
 end subroutine allocate_jet
 
 subroutine deallocate_jet() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for deallocating arrays which 
!     describe the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  deallocate(jetfr)
  if(lmultiplestep)deallocate(jetfm)
  deallocate(jetpt)
  deallocate(jetxx)
  deallocate(jetyy)
  deallocate(jetzz)
  deallocate(jetst)
  deallocate(jetvx)
  deallocate(jetvy)
  deallocate(jetvz)
  deallocate(jetms)
  deallocate(jetch)
  deallocate(jetcr)
  deallocate(jetvl)
  if(ltrackbeads .and. idrank==0)deallocate(jetlb)
  if(typemass==3 .or. ltagbeads)deallocate(jetbd)
  if(lbreakup)deallocate(jetbr)
  
  return
  
 end subroutine deallocate_jet
 
 subroutine reallocate_jet() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for reallocating arrays which 
!     describe the nanojet if the system size is changed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  integer :: oldinpjet,oldnpjet,oldmxnpjet
  integer :: newinpjet,newnpjet
  integer :: oldinit,oldend
  integer :: newinit,newend,ncutoffsub
  
  
  ncutoffsub=1
  doallocate=.false.
  doreorder=.true.
  
  if(inpjet==0)then
    oldinit=0
    oldend=npjet-1
    oldinpjet=inpjet
    oldnpjet=npjet-1
    newinit=0
    newend=oldend-oldinit
    
    newinpjet=oldinpjet-oldinit
    newnpjet=oldnpjet-oldinit
    
    if((oldend-newend)<incnpjet)then
      oldmxnpjet=mxnpjet
      mxnpjet=oldmxnpjet+incnpjet
      call set_mxchunk(mxnpjet)
      doallocate=.true.
    endif
  elseif(inpjet>0)then
    oldinit=max(0,inpjet-ncutoffsub)
    oldend=npjet-1
    oldinpjet=inpjet
    oldnpjet=npjet-1
    newinit=0
    newend=oldend-oldinit
    
    newinpjet=oldinpjet-oldinit
    newnpjet=oldnpjet-oldinit
    
    if((oldend-newend)<incnpjet)then
      oldmxnpjet=mxnpjet
      mxnpjet=oldmxnpjet+incnpjet
      call set_mxchunk(mxnpjet)
      doallocate=.true.
    endif
  endif
  
  call allocate_array_lbuffservice(newend)
  
  lbuffservice(newinit:newend)=jetfr(oldinit:oldend)
  if(doallocate)then
    deallocate(jetfr)
    allocate(jetfr(0:mxnpjet))
  endif
  jetfr(:)=.false.
  jetfr(newinit:newend)=lbuffservice(newinit:newend)
  
  if(doallocate .and. lmultiplestep)then
    deallocate(jetfm)
    allocate(jetfm(0:mxnpjet))
  endif
  
  call allocate_array_buffservice(newend)
  
  buffservice(newinit:newend)=jetpt(oldinit:oldend)
  if(doallocate)then
    deallocate(jetpt)
    allocate(jetpt(0:mxnpjet))
  endif
  jetpt(:)=0.d0
  jetpt(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetxx(oldinit:oldend)
  if(doallocate)then
    deallocate(jetxx)
    allocate(jetxx(0:mxnpjet))
  endif
  jetxx(:)=0.d0
  jetxx(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetyy(oldinit:oldend)
  if(doallocate)then
    deallocate(jetyy)
    allocate(jetyy(0:mxnpjet))
  endif
  jetyy(:)=0.d0
  jetyy(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetzz(oldinit:oldend)
  if(doallocate)then
    deallocate(jetzz)
    allocate(jetzz(0:mxnpjet))
  endif
  jetzz(:)=0.d0
  jetzz(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetst(oldinit:oldend)
  if(doallocate)then
    deallocate(jetst)
    allocate(jetst(0:mxnpjet))
  endif
  jetst(:)=0.d0
  jetst(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetvx(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvx)
    allocate(jetvx(0:mxnpjet))
  endif
  jetvx(:)=0.d0
  jetvx(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetvy(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvy)
    allocate(jetvy(0:mxnpjet))
  endif
  jetvy(:)=0.d0
  jetvy(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetvz(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvz)
    allocate(jetvz(0:mxnpjet))
  endif
  jetvz(:)=0.d0
  jetvz(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetms(oldinit:oldend)
  if(doallocate)then
    deallocate(jetms)
    allocate(jetms(0:mxnpjet))
  endif
  jetms(:)=0.d0
  jetms(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetch(oldinit:oldend)
  if(doallocate)then
    deallocate(jetch)
    allocate(jetch(0:mxnpjet))
  endif
  jetch(:)=0.d0
  jetch(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetcr(oldinit:oldend)
  if(doallocate)then
    deallocate(jetcr)
    allocate(jetcr(0:mxnpjet))
  endif
  jetcr(:)=0.d0
  jetcr(newinit:newend)=buffservice(newinit:newend)
  
  buffservice(newinit:newend)=jetvl(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvl)
    allocate(jetvl(0:mxnpjet))
  endif
  jetvl(:)=0.d0
  jetvl(newinit:newend)=buffservice(newinit:newend)
  
  if(typemass==3 .or. ltagbeads)then
    lbuffservice(newinit:newend)=jetbd(oldinit:oldend)
    if(doallocate)then
      deallocate(jetbd)
      allocate(jetbd(0:mxnpjet))
    endif
    jetbd(:)=.false.
    jetbd(newinit:newend)=lbuffservice(newinit:newend)
  endif
  
  if(lbreakup)then
    lbuffservice(newinit:newend)=jetbr(oldinit:oldend)
    if(doallocate)then
      deallocate(jetbr)
      allocate(jetbr(0:mxnpjet))
    endif
    jetbr(:)=.false.
    jetbr(newinit:newend)=lbuffservice(newinit:newend)
  endif
  
  if(ltrackbeads .and. idrank==0)then
    call allocate_array_ibuffservice(newend)
    ibuffservice(newinit:newend)=jetlb(oldinit:oldend)
    if(doallocate)then
      deallocate(jetlb)
      allocate(jetlb(0:mxnpjet))
    endif
    jetlb(:)=-1
    jetlb(newinit:newend)=ibuffservice(newinit:newend)
  endif
  
  inpjet=newinpjet
  npjet=newnpjet+1
  
  return
  
 end subroutine reallocate_jet
 
 subroutine set_initial_jet(initime,endtime,refinementthreshold, &
  refbeadstartfit)
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining the initial quantities and
!     dimensional scales of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(inout) :: initime,endtime, &
   refinementthreshold,refbeadstartfit
  integer :: i
  double precision :: di
  
  ivolume=Pi*(icrossec**2.d0)*resolution
  meanmass=imassa*ivolume
  massscale=meanmass
  
  meancharge=icharge*ivolume 
  chargescale=meancharge
  regq=0.d0 !  icrossec  !0.d0 !meancharge/1.d10
  
  if(.not.lV0)then
    call warning(17)
    V0=0.d0
    call warning(65,V0)
  endif
  
  if(.not.lsurfacet.and.systype/=1)then
    call warning(18)
    surfacet=0.d0
    call warning(64,surfacet)
  endif
  
  if(laird.and.lairv)then
    tanfriction=0.65d0*Pi*aird*(2.d0/airv)**(-0.81d0)* &
    (resolution**0.095d0)*(icrossec**0.19d0)
    airdragamp(1:3)=airdragamp(1:3)*tanfriction/meanmass
  else
    tanfriction=0.d0
    airdragamp(1:3)=0.d0
  endif
  
  lengthscale=dsqrt((chargescale**2.d0)/(Pi*(icrossec**2.d0)*G))
  
  call set_dimensionless_patameters(initime,endtime)
  
  hresolution=resolution/2.d0
  thresolution=resolution/2.d0*3.d0
  dresolution=resolution*2.d0
  
  if(.not.livelocity)then
    call warning(13)
    ivelocity=0.d0
    call warning(26,ivelocity)
  endif
  
  if(.not.listress)then
    call warning(12)
    if(listresscompensate)then
      di=0.d0
      do i=1,npjet
        di=di+(1.d0/(dble(i)**2.d0))
      enddo
      istress=meancharge**2.d0/(Pi*(icrossec*resolution)**2.d0)*di
    else
      istress=0.d0
    endif
    call warning(24,istress)
  endif
  
  if(.not.lyieldstress)then
    call warning(27)
    if(listresscompensate)then
      yieldstress=istress
    else
      yieldstress=0.d0
    endif
    call warning(28,yieldstress)
  endif
  
  jetfr(:)=.false.
  
  jetxx(:)=0.d0
  jetyy(:)=0.d0
  jetzz(:)=0.d0
  
  
  do i=0,npjet
    jetxx(i) = ilength-(resolution*dble(i))
  enddo
  
  if(typemass==3 .or. ltagbeads)then
    jetbd(:)=.false.
  endif
  
  if(lbreakup)then
    jetbr(:)=.false.
  endif
  
  if(ltrackbeads .and. idrank==0)then
    jetlb(:)=-1
    do i=0,npjet
      jetlb(i) = trackbeads_number()
    enddo
  endif
  
  jetst(:)=0.d0
  
  jetvx(:)=0.d0
  jetvy(:)=0.d0
  jetvz(:)=0.d0
  
  jetms(:)=0.d0
  jetch(:)=0.d0
  jetvl(:)=0.d0
  
  do i=0,npjet
    jetst(i)=istress
    jetvx(i)=ivelocity
    jetms(i)=imassa*ivolume
    jetch(i)=icharge*ivolume
    jetvl(i)=ivolume
  enddo
  
  call conv_dimensionless_unit_jet(initime,endtime, &
   refinementthreshold,refbeadstartfit)
  call add_initial_perturbation()
  
  return
  
 end subroutine set_initial_jet
 
 subroutine set_dimensionless_patameters(initime,endtime)
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining the
!     dimensionless parameters of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(inout) :: initime,endtime
  
  tao=(mu/G)
  
  q=(chargescale*mu)**2.d0/((lengthscale**3.d0)*massscale*(G**2.d0))
  v=(chargescale*V0*(mu**2.d0))/(h*lengthscale*massscale*(G**2.d0))
  fve=lengthscale*(mu**2.d0)/(massscale*G)
  fvere=Pi*(icrossec*mu)**2.d0/(massscale*G*lengthscale)
  Hg=h/lengthscale
  Lrg=resolution/lengthscale
  if(lflorentz)then
    KLor=(chargescale)/(dsqrt(massscale)*dsqrt(lengthscale)* &
     29979245800.d0)
  endif
  if(systype/=1)then
    ks=surfacet*(mu**2.d0)/(massscale*(G**2.d0))
    ksre=surfacet*Pi*(mu**2.d0)*(icrossec**2.d0)/(massscale*(G**2.d0)* &
     (lengthscale**2.d0))
  endif
  if(lairdrag)then
    Li=aird*(lengthscale**3.d0)/massscale
    Lire=aird*Pi*(icrossec**2.d0)/massscale
    att=tanfriction*tao*(lengthscale**0.905d0)* &
     ((lengthscale/tao)**0.19d0)/massscale
  endif
  if(lgravity)then
    Gr=980.665d0*(tao**2.d0)/(lengthscale)
  endif
  
  sqrq=dsqrt(q)
  
  return
  
 end subroutine set_dimensionless_patameters
 
 subroutine conv_dimensionless_unit_jet(initime,endtime, &
  refinementthreshold,refbeadstartfit)
 
!***********************************************************************
!     
!     JETSPIN subroutine for converting quatities in
!     dimensionless units
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(inout) :: initime,endtime
  double precision, intent(inout) :: refinementthreshold, &
   refbeadstartfit
  
  
  tstep=tstep/tao
  initime=initime/tao
  endtime=endtime/tao
  jetxx(:)=jetxx(:)/lengthscale
  jetyy(:)=jetyy(:)/lengthscale
  jetzz(:)=jetzz(:)/lengthscale
  jetst(:)=jetst(:)/G
  jetvx(:)=jetvx(:)*(tao/lengthscale)
  jetvy(:)=jetvy(:)*(tao/lengthscale)
  jetvz(:)=jetvz(:)*(tao/lengthscale)
  jetms(:)=jetms(:)/massscale
  jetch(:)=jetch(:)/chargescale
  jetvl(:)=jetvl(:)/(lengthscale**3.d0)
  regq=regq/lengthscale
  resolution=resolution/lengthscale
  hresolution=hresolution/lengthscale
  thresolution=thresolution/lengthscale
  dresolution=dresolution/lengthscale
  refinementthreshold=refinementthreshold/lengthscale
  refbeadstartfit=refbeadstartfit/lengthscale
  h=h/lengthscale
  ilength=ilength/lengthscale
  icrossec=icrossec/lengthscale
  istress=istress/G
  yieldstress=yieldstress/G
  ivelocity=ivelocity*(tao/lengthscale)
  velext=velext*(tao/lengthscale)
  consistency=consistency/(mu*((tao)**(findex-1.d0)))
  meanmass=meanmass/massscale
  meancharge=meancharge/chargescale
  imassa=imassa/massscale*(lengthscale**3.d0)
  icharge=icharge/chargescale*(lengthscale**3.d0)
  BLor(1:3)=BLor(1:3)*dsqrt(lengthscale)*tao/dsqrt(massscale)
  if(lultimatestrength)then
    ultimatestrength=ultimatestrength/G
  endif
  if(ltagbeads)then
    lenthresholdbead=lenthresholdbead/lengthscale
  endif
  if(lmassavariable)then
    lencorrmassa=lencorrmassa/lengthscale
  endif
  ivolume=ivolume/(lengthscale**3.d0)
  pampl=pampl/lengthscale
  if(liniperturb)then
    pfreq=pfreq*tao
  endif
  if(lairdrag)then
    airdragamp(1:3)=airdragamp(1:3)*(tao**3.d0)/(lengthscale**2.d0)
  endif
  if(.not.lxyzrescale)then
    xyzrescale=1.d0
  else
    xyzrescale=xyzrescale*lengthscale
  endif
  if(.not.lpdbrescale)then
    pdbrescale=1.d0
  else
    pdbrescale=pdbrescale*lengthscale
  endif
  if(luppot)then
    kuppot=kuppot*(tao**2.d0)/(massscale)
  endif
  if(lfieldfreq)then
    fieldfreq=fieldfreq*tao
  endif
  if(ltaoelectr)then
    taoelectr=taoelectr/tao
  endif
  if(lmultiplestep)then
    dcutoff=dcutoff/lengthscale
    maxdispl=maxdispl/lengthscale
  endif
  if(typemass==3)then
    lenprobmassa=lenprobmassa/lengthscale
    massratio=massratio/massscale*(lengthscale**3.d0)
  endif
    
  return
  
 end subroutine conv_dimensionless_unit_jet
 
 subroutine deconv_dimensionless_unit_jet()
 
!***********************************************************************
!     
!     JETSPIN subroutine for deconverting quatities from
!     dimensionless units to cgs units
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  jetxx(:)=jetxx(:)*ilength
  jetyy(:)=jetyy(:)*ilength
  jetzz(:)=jetzz(:)*ilength
  jetst(:)=jetst(:)*G
  jetvx(:)=jetvx(:)/(tao/ilength)
  jetvy(:)=jetvy(:)/(tao/ilength)
  jetvz(:)=jetvz(:)/(tao/ilength)
  
  
  return
  
 end subroutine deconv_dimensionless_unit_jet
 
 subroutine add_jetbead(nstepsub,timesub,ladd)
 
!***********************************************************************
!     
!     JETSPIN subroutine for adding a new bead at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  double precision, intent(in) :: timesub
  logical, intent(inout) :: ladd
  
  double precision :: tempmod,tempmod0,actualimass,actualicharge
  
  if(.not.linserting)return
  timedeposition=timedeposition+tstep
  if(.not.linserted)then
!   if a new bead was inserted but blocked check if its time evolution
!   should be started
    select case(systype)
    case(1:2)
      tempmod=jetxx(npjet-2)-jetxx(npjet)
    case default
      tempmod=dsqrt((jetxx(npjet-2)-jetxx(npjet))**2.d0+ &
       (jetyy(npjet-2)-jetyy(npjet))**2.d0+ &
       (jetzz(npjet-2)-jetzz(npjet))**2.d0)
    end select
    if(tempmod>=dresolution)then
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      linserted=.true.
      jetst(npjet-1)=0.d0!jetst(npjet-2)
      if(ldragvel)then
        if(typedragvel==1)then
          tempmod=dsqrt(((jetvx(npjet-2)-ivelocity)/2.d0)**2.d0+ &
           (jetvy(npjet-2)/2.d0)**2.d0+(jetvz(npjet-2)/2.d0)**2.d0)
          jetvx(npjet-1)=tempmod+ivelocity
          jetvy(npjet-1)=0.d0
          jetvz(npjet-1)=0.d0
        elseif(typedragvel==2)then
          jetvx(npjet-1)=(jetvx(npjet-2)-ivelocity)/2.d0+ivelocity
          jetvy(npjet-1)=0.d0
          jetvz(npjet-1)=0.d0
        elseif(typedragvel==3)then
          tempmod=dsqrt((jetvy(npjet-2)/2.d0)**2.d0+ &
           (jetvz(npjet-2)/2.d0)**2.d0)/ &
           (dsqrt(((jetvx(npjet-2)-ivelocity)/2.d0)**2.d0+ &
           (jetvy(npjet-2)/2.d0)**2.d0+(jetvz(npjet-2)/2.d0)**2.d0))
          jetvx(npjet-1)=(jetvx(npjet-2)-ivelocity)/2.d0+ivelocity
          jetvy(npjet-1)=(1.d0-tempmod)*jetvy(npjet-2)/2.d0
          jetvz(npjet-1)=(1.d0-tempmod)*jetvz(npjet-2)/2.d0
        else
          jetvx(npjet-1)=(jetvx(npjet-2)-ivelocity)/2.d0+ivelocity
          jetvy(npjet-1)=jetvy(npjet-2)/2.d0
          jetvz(npjet-1)=jetvz(npjet-2)/2.d0
        endif
      else
        jetvx(npjet-1)=ivelocity
        jetvy(npjet-1)=0.d0
        jetvz(npjet-1)=0.d0
      endif
      insvx=jetvx(npjet-1)
      insvy=jetvy(npjet-1)
      insvz=jetvz(npjet-1)
    endif
  else
!   check if a new bead should be inserted
    select case(systype)
    case(1:2)
      tempmod0=jetxx(npjet-1)-jetxx(npjet)
    case default
      tempmod0=dsqrt((jetxx(npjet-1)-jetxx(npjet))**2.d0+ &
       (jetyy(npjet-1)-jetyy(npjet))**2.d0+ &
       (jetzz(npjet-1)-jetzz(npjet))**2.d0)
    end select
    if(tempmod0>=thresolution)then
      npjet=npjet+1
      if(npjet>mxnpjet)call reallocate_jet()
      jetfr(npjet)=jetfr(npjet-1)
      jetxx(npjet)=jetxx(npjet-1)
      jetyy(npjet)=jetyy(npjet-1)
      jetzz(npjet)=jetzz(npjet-1)
      jetst(npjet)=jetst(npjet-1)
      jetvx(npjet)=jetvx(npjet-1)
      jetvy(npjet)=jetvy(npjet-1)
      jetvz(npjet)=jetvz(npjet-1)
      jetms(npjet)=jetms(npjet-1)
      jetch(npjet)=jetch(npjet-1)
      jetvl(npjet)=jetvl(npjet-1)
      if(typemass==3 .or. ltagbeads)jetbd(npjet)=jetbd(npjet-1)
      if(lbreakup)jetbr(npjet)=jetbr(npjet-1)
      if(ltrackbeads .and. idrank==0)jetlb(npjet)=jetlb(npjet-1)
      jetfr(npjet-1)=.false.
      jetxx(npjet-1)=resolution
      jetyy(npjet-1)=0.d0
      jetzz(npjet-1)=0.d0
      jetst(npjet-1)=istress
      jetvx(npjet-1)=ivelocity
      jetvy(npjet-1)=0.d0
      jetvz(npjet-1)=0.d0
      call tag_beads()
      call extract_actual_massdensity(actualimass,timesub)
      jetms(npjet-1)=actualimass*ivolume 
      jetch(npjet-1)=icharge*ivolume
      jetvl(npjet-1)=ivolume
      if(ltrackbeads .and. idrank==0)jetlb(npjet-1)=trackbeads_number()
      ladd=.true.
      naddtrack=naddtrack+1
      timedeposition=0.d0
      if(lmultiplestep)lneighlistdo=.true.
!     the new bead is blocked until a given condition is not satisfied
      linserted=.false.
      call compute_posnoinserted(jetxx,jetyy,jetzz)
    endif
  endif
  
  return
 
 end subroutine add_jetbead
 
 subroutine add_initial_perturbation()
 
!***********************************************************************
!     
!     JETSPIN subroutine for adding the initial position of the nozzle
!     on the y-z plane
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision :: drand,theta
  
  if(.not.liniperturb)return
  
  select case(systype)
  case(1:2)
    return
  case(3:4)
    if(idrank==0)then
      call random_number(drand)
    endif
    call bcast_world_d(drand)
    theta=2.d0*Pi*drand
    jetyy(npjet)=pampl*dcos(theta)
    jetzz(npjet)=pampl*dsin(theta)
  end select
  
  return
  
 end subroutine add_initial_perturbation
 
 pure function fcut(r,inner_cut,outer_cut)
 
!***********************************************************************
!     
!     JETSPIN function for fading an observable (r) within a given 
!     interval
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
      
  implicit none
  
  double precision,intent(in) :: r,inner_cut,outer_cut
  double precision :: fcut
  
  if ( r < inner_cut ) then
    fcut = 0.d0
  elseif ( r >= inner_cut ) then
    if ( r >= outer_cut ) then
      fcut = 1.d0
    else
      fcut = -0.5d0*cos((r-inner_cut)*Pi/(outer_cut-inner_cut))+0.5d0
    endif
  end if
      
  return
  
 end function fcut
 
 subroutine remove_jetbead(nstepsub,nremoved,timesub,lrem,lremdat)
 
!***********************************************************************
!     
!     JETSPIN subroutine for removing the last bead at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  integer, intent(out) :: nremoved
  double precision, intent(in) :: timesub
  logical, intent(inout) :: lrem,lremdat
  
  integer :: ipoint,newinpjet
  integer, parameter :: strategysub=0
  
  if(.not.lremove)return
  
  
  if(strategysub==1)then
    lrem=.false.
    do ipoint=inpjet,npjet
      if(jetxx(ipoint)>=h)then
        newinpjet=ipoint
        lrem=.true.
        if(lmultiplestep)lneighlistdo=.true.
      endif
    enddo
    
    if(lrem)then
      nremoved=newinpjet-inpjet+1
      nremtrack=nremtrack+nremoved
      inpjet=newinpjet+1
      lremdat=.true.
    else
      lremdat=.false.
    endif
    
  else
    
    do ipoint=inpjet,npjet
      if(jetxx(ipoint)>=h)then
        jetfr(ipoint)=.true.
        jetxx(ipoint)=h
      endif
    enddo
    
    lrem=.false.
    ipoint=inpjet
    if(jetxx(ipoint)>=h)then
      if(jetxx(ipoint+1)>=h)then
        newinpjet=ipoint
        lrem=.true.
      endif
    endif
    
    if(lrem)then
      nremoved=newinpjet-inpjet+1
      nremtrack=nremtrack+nremoved
      inpjet=newinpjet+1
      lremdat=.true.
    else
      lremdat=.false.
    endif
  
  endif
    
    
  return
  
 end subroutine remove_jetbead
 
 subroutine erase_jetbead(nstepsub,timesub,ladd,lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for erasing the data of the last bead 
!     at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstepsub,nremovedsub
  double precision, intent(in) :: timesub
  logical, intent(inout) :: ladd,lrem
  
  integer :: ipoint
  
  if(lrem)then
    do ipoint=inpjet-nremovedsub,inpjet-1
      jetvl(ipoint)=0.d0
      jetch(ipoint)=0.d0
      jetms(ipoint)=0.d0
      jetst(ipoint)=0.d0
      jetvx(ipoint)=0.d0
      jetvy(ipoint)=0.d0
      jetvz(ipoint)=0.d0
      if(typemass==3 .or. ltagbeads)jetbd(ipoint)=.false.
    enddo
  endif
  
  ladd=.false.
  lrem=.false.
  
  return
  
 end subroutine erase_jetbead
 
 subroutine compute_posnoinserted(yxx,yyy,yzz,timesub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for updating the position of the last 
!     inserted bead at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) ::  yxx
  double precision, allocatable, dimension (:), intent(inout) ::  yyy
  double precision, allocatable, dimension (:), intent(inout) ::  yzz
  double precision,optional, intent(in) :: timesub
  
  double precision, dimension (3) :: v1
  double precision :: l1,sol1  
  
  if(linserted)return
  
! the position is computed by a linear interpolation
  v1(1)=(yxx(npjet-2)-yxx(npjet))
  v1(2)=(yyy(npjet-2)-yyy(npjet))
  v1(3)=(yzz(npjet-2)-yzz(npjet))
  l1=modulvec(v1)
  sol1=resolution/l1
  yxx(npjet-1)=yxx(npjet)+sol1*v1(1)
  yyy(npjet-1)=yyy(npjet)+sol1*v1(2)
  yzz(npjet-1)=yzz(npjet)+sol1*v1(3)
  
  return
  
 end subroutine compute_posnoinserted
 
 subroutine tag_beads()
 
!***********************************************************************
!     
!     JETSPIN subroutine for tagging the beads if the dynamic
!     refinement procedure is activated
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, save :: lfirst=.true.
  logical :: ljetbd
  double precision, save :: olddistancesub
  
  if(.not. ltagbeads)return
  
  if(lfirst)then
    lfirst=.false.
    olddistancesub=0.d0
    if(olddistancesub<lenthresholdbead)then
      ljetbd=.false.
    else
      ljetbd=.true.
    endif
  else
    olddistancesub=olddistancesub+resolution
    if(olddistancesub<lenthresholdbead)then
      ljetbd=.false.
    else
      ljetbd=.true.
      olddistancesub=0.d0
    endif
  endif
  
  jetbd(npjet-1)=ljetbd
  
  return
  
 end subroutine tag_beads
  
 subroutine extract_actual_massdensity(actualmass,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for extracting randomically the mass of 
!     the last inserted bead at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2016
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(out) :: actualmass
  double precision, intent(in), optional :: timesub
  
  logical, save :: lfirst=.true.
  
  double precision :: actualgaussn,newgaussn,volumenp
  
  actualmass=0.d0
  if(lfirst .and. lfirstmass)then
    lfirst=.false.
    lfirstmass=.false.
    if(lmassavariable .and. ldevelopers)then
      if(idrank==0)then
        if(ltagbeads)then
          if(jetbd(npjet-1))then
            volumenp=(4.d0/3.d0)*Pi*lenprobmassa**3.d0
            actualmass=(imassa*(ivolume-volumenp)+massratio*volumenp)/ &
             ivolume
          else
            actualmass=imassa
          endif
        endif
      endif
      call bcast_world_d(actualmass)
    else
      actualmass=imassa
    endif
  else
    if(lmassavariable .and. ldevelopers)then
      if(idrank==0)then
        if(ltagbeads)then
          if(jetbd(npjet-1))then
            volumenp=(4.d0/3.d0)*Pi*lenprobmassa**3.d0
            actualmass=(imassa*(ivolume-volumenp)+massratio*volumenp)/ &
             ivolume
          else
            actualmass=imassa
          endif
        endif
      endif
      call bcast_world_d(actualmass)
    else
      actualmass=imassa
    endif
  endif
  
  return
  
 end subroutine extract_actual_massdensity
 
 function trackbeads_number()
 
  implicit none
  
  integer :: trackbeads_number
  integer, parameter :: numlimit=huge(1)
  integer, save :: lastnumber=0
  
  lastnumber=lastnumber+1
  if(lastnumber>numlimit)then
    lastnumber=1
    lreordertrack=.true.
  endif
  
  trackbeads_number=lastnumber
  
  return
  
 end function trackbeads_number
 
 end module nanojet_mod


