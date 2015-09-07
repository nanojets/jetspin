
 module nanojet_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which are defining
!     and dealing the quantities describing the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************

 use version_mod, only : idrank,set_mxchunk,bcast_world_d
 use error_mod,   only : error,warning
 use utility_mod, only : Pi,modulvec,cross,dot,gauss

 implicit none

 private
 
 integer, public, save :: nvariable
 integer, public, save :: ndimension
 integer, public, save :: inpjet
 integer, public, save :: npjet
 integer, public, save :: mxnpjet
 integer, private, parameter :: incnpjet=100
 integer, public, save :: systype
 integer, public, save :: units
 integer, public, save :: ncutoff
 integer, public, save :: insertmode
 integer, public, save :: myseed=1
 logical, public, save :: doallocate
 logical, public, save :: doreorder
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
 double precision, public, save :: imassa!=2.25487812711407d-4 !grammi
 double precision, public, save :: imassadev=0.d0 !grammi^0.5 cm^1.5 sec^-1
 double precision, public, save :: lencorrmassa=0.0d0
 double precision, public, save :: icharge!=9.37212742146083d-4 !grammi^0.5 cm^1.5 sec^-1
 double precision, public, save :: ichargedev=0.d0 !grammi^0.5 cm^1.5 sec^-1
 double precision, public, save :: lencorrcharge=0.d0
 double precision, public, save :: mu!=12211.32d0 !grammi cm^-1 sec^-1
 double precision, public, save :: G!=122132.d0 ! grammi cm^-1 sec^-2
 double precision, public, save :: surfacet!=700.d0 ! grammi sec^-2
 double precision, public, save :: icrossec!=1.5d-4 ! cm
 double precision, public, save :: V0!=30.69980124d0 ! cm^0.5 grammi^0.5 sec^-1
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
 double precision, public, save :: current=0.d0 !grammi^0.5 cm^1.5 sec^-2
 double precision, public, save :: timedeposition=0.d0 !sec
 double precision, public, save :: meancharge
 double precision, public, save :: meanmass
 double precision, public, save :: regq
 double precision, public, save :: tanfriction
 double precision, public, save :: massscale
 double precision, public, save :: chargescale
 double precision, public, save :: tstep=0.d0
 double precision, public, save :: xyzrescale=1.d0
 double precision, public, save :: insvx=0.d0
 double precision, public, save :: insvy=0.d0
 double precision, public, save :: insvz=0.d0
 double precision, public, save :: dcutoff=0.d0
 
 double precision, private, save :: vdumbbell
 
 double precision, allocatable, public, save :: jetxx(:)
 double precision, allocatable, public, save :: jetyy(:)
 double precision, allocatable, public, save :: jetzz(:)
 double precision, allocatable, public, save :: jetst(:)
 double precision, allocatable, public, save :: jetvx(:)
 double precision, allocatable, public, save :: jetvy(:)
 double precision, allocatable, public, save :: jetvz(:)
 double precision, allocatable, public, save :: jetms(:) !grammi
 double precision, allocatable, public, save :: jetch(:) !grammi^0.5 cm^1.5 sec^-1
 double precision, allocatable, public, save :: jetcr(:)! cm
 
 logical, private, parameter :: listresscompensate=.false.
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
 logical, public, save :: lncutoff=.false.
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
 logical, public, save :: lconsistency=.false.
 logical, public, save :: lfindex=.false.
 logical, public, save :: ltstep=.false.
 logical, public, save :: lgravity=.false.
 logical, public, save :: lmassavariable=.false.
 logical, public, save :: lchargevariable=.false.
 logical, public, save :: limassadev=.false.
 logical, public, save :: llencorrmassa=.false.
 logical, public, save :: llencorrcharge=.false.
 logical, public, save :: lichargedev=.false.
 logical, public, save :: lxyzrescale=.false.
 
 
 public :: set_resolution_length
 public :: allocate_jet
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
! the number of beads (npjet) in which the jet is initially discretised is
! given in input file, the length step is computed by the division
! ilength/npjet
  if(lnpjet)then
    resolution=ilength/dble(npjet)
    call warning(21,dble(npjet))
    call warning(23,resolution)
  endif
  
! set the cutoff for the Coulomb force if it was not set in input file
  if(.not.lncutoff)then
     dcutoff=h
     ncutoff=ceiling(dcutoff/resolution)
     call warning(30,dcutoff)
  else
    ncutoff=ceiling(dcutoff/resolution)
  endif
  
! the length step has to be grater than the cross section radius
! of jet otherwise the model loses validity
  if(resolution<icrossec)then
    call warning(62)
    call error(8)
  endif
  
  return
  
 end subroutine set_resolution_length
 
 subroutine allocate_jet() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating arrays which 
!     describe the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  
  select case(systype)
    case(1:2)
      nvariable = 3
      ndimension = 1
    case default
      nvariable = 3
      ndimension = 3
  end select
  
  mxnpjet=incnpjet
  mxnpjet=max(mxnpjet,npjet)
  call set_mxchunk(mxnpjet)
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
  
  
  doallocate=.true.
  
  
  return
  
 end subroutine allocate_jet
 
 subroutine reallocate_jet() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for reallocating arrays which 
!     describe the nanojet if the system size is changed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer :: oldinpjet,oldnpjet,oldmxnpjet
  integer :: newinpjet,newnpjet
  integer :: oldinit,oldend
  integer :: newinit,newend,ncutoffsub
  
  double precision, allocatable :: buff(:)
  
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
  
  allocate(buff(newinit:newend))
  
  buff(newinit:newend)=jetxx(oldinit:oldend)
  if(doallocate)then
    deallocate(jetxx)
    allocate(jetxx(0:mxnpjet))
  endif
  jetxx(:)=0.d0
  jetxx(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetyy(oldinit:oldend)
  if(doallocate)then
    deallocate(jetyy)
    allocate(jetyy(0:mxnpjet))
  endif
  jetyy(:)=0.d0
  jetyy(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetzz(oldinit:oldend)
  if(doallocate)then
    deallocate(jetzz)
    allocate(jetzz(0:mxnpjet))
  endif
  jetzz(:)=0.d0
  jetzz(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetst(oldinit:oldend)
  if(doallocate)then
    deallocate(jetst)
    allocate(jetst(0:mxnpjet))
  endif
  jetst(:)=0.d0
  jetst(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetvx(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvx)
    allocate(jetvx(0:mxnpjet))
  endif
  jetvx(:)=0.d0
  jetvx(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetvy(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvy)
    allocate(jetvy(0:mxnpjet))
  endif
  jetvy(:)=0.d0
  jetvy(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetvz(oldinit:oldend)
  if(doallocate)then
    deallocate(jetvz)
    allocate(jetvz(0:mxnpjet))
  endif
  jetvz(:)=0.d0
  jetvz(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetms(oldinit:oldend)
  if(doallocate)then
    deallocate(jetms)
    allocate(jetms(0:mxnpjet))
  endif
  jetms(:)=0.d0
  jetms(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetch(oldinit:oldend)
  if(doallocate)then
    deallocate(jetch)
    allocate(jetch(0:mxnpjet))
  endif
  jetch(:)=0.d0
  jetch(newinit:newend)=buff(newinit:newend)
  
  buff(newinit:newend)=jetcr(oldinit:oldend)
  if(doallocate)then
    deallocate(jetcr)
    allocate(jetcr(0:mxnpjet))
  endif
  jetcr(:)=0.d0
  jetcr(newinit:newend)=buff(newinit:newend)
  
  deallocate(buff)
  
  inpjet=newinpjet
  npjet=newnpjet+1
  
  return
  
 end subroutine reallocate_jet
 
 subroutine set_initial_jet(tstep,initime,endtime)
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining the initial quantities and
!     dimensional scales of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(inout) :: tstep,initime,endtime
  integer :: i
  double precision :: di
  
  vdumbbell=Pi*(icrossec**2.d0)*resolution
  jetms(:)=imassa*vdumbbell
  meanmass=imassa*vdumbbell
  massscale=meanmass
  
  jetch(:)=icharge*vdumbbell
  meancharge=icharge*vdumbbell 
  chargescale=meancharge
  regq=0.d0!meancharge/1.d10
  
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
  
  call set_dimensionless_patameters(tstep,initime,endtime)
  
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
  
  
  
  jetxx(:)=0.d0
  jetyy(:)=0.d0
  jetzz(:)=0.d0
  
  
  do i=0,npjet
    jetxx(i) = ilength-(resolution*dble(i))
  enddo
  
  jetst(:)=istress
  
  jetvx(:)=ivelocity
  jetvy(:)=0.d0
  jetvz(:)=0.d0
  
  
  call conv_dimensionless_unit_jet(tstep,initime,endtime)
  call add_initial_perturbation()
  
  return
  
 end subroutine set_initial_jet
 
 subroutine set_dimensionless_patameters(tstep,initime,endtime)
 
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
  
  double precision, intent(inout) :: tstep,initime,endtime
  
  tao=(mu/G)
  
  q=(chargescale*mu)**2.d0/((lengthscale**3.d0)*massscale*(G**2.d0))
  v=(chargescale*V0*(mu**2.d0))/(h*lengthscale*massscale*(G**2.d0))
  fve=Pi*(icrossec*mu)**2.d0/(massscale*G*lengthscale)
  Hg=h/lengthscale
  Lrg=resolution/lengthscale
  if(systype/=1)then
    ks=surfacet*Pi*(mu**2.d0)*(icrossec**2.d0)/(massscale*(G**2.d0)* &
     (lengthscale**2.d0))
  endif
  if(lairdrag)then
    Li=aird*Pi*(icrossec**2.d0)/massscale
    att=tanfriction*tao*(lengthscale**0.905d0)*((lengthscale/tao)**0.19d0)/massscale
  endif
  if(lgravity)then
    Gr=980.665d0*(tao**2.d0)/(lengthscale)
  endif
  
  sqrq=dsqrt(q)
  
  return
  
 end subroutine set_dimensionless_patameters
 
 subroutine conv_dimensionless_unit_jet(tstep,initime,endtime)
 
!***********************************************************************
!     
!     JETSPIN subroutine for converting quatities in
!     dimensionless units
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(inout) :: tstep,initime,endtime
  
  
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
  regq=regq/lengthscale
  resolution=resolution/lengthscale
  hresolution=hresolution/lengthscale
  thresolution=thresolution/lengthscale
  dresolution=dresolution/lengthscale
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
  if(lmassavariable)then
    lencorrmassa=lencorrmassa/lengthscale
    imassadev=imassadev/massscale*(lengthscale**3.d0)
  endif
  if(lchargevariable)then
    lencorrcharge=lencorrcharge/lengthscale
    ichargedev=ichargedev/chargescale*(lengthscale**3.d0)
  endif  
  vdumbbell=vdumbbell/(lengthscale**3.d0)
  pampl=pampl/lengthscale
  if(liniperturb)then
    pfreq=pfreq*tao
  endif
  if(lairdrag)then
    airdragamp(1:3)=airdragamp(1:3)*(tao**3.d0)/(lengthscale**2.d0)
  endif
  if(lxyzrescale)then
    xyzrescale=1.d0
  else
    xyzrescale=lengthscale
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
 
 subroutine add_jetbead(nstepsub,timesub,tstepsub,ladd)
 
!***********************************************************************
!     
!     JETSPIN subroutine for adding a new bead at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  double precision, intent(in) :: timesub,tstepsub
  logical, intent(inout) :: ladd
  
  double precision :: tempmod,tempmod0,actualimass,actualicharge
  
  if(.not.linserting)return
  timedeposition=timedeposition+tstepsub
  if(.not.linserted)then
!   if a new bead was inserted but blocked check if its time evolution
!   should be started
    select case(systype)
    case(1:2)
      tempmod=jetxx(npjet-2)-jetxx(npjet)
    case default
      tempmod=dsqrt((jetxx(npjet-2)-jetxx(npjet))**2.d0+ &
       (jetyy(npjet-2)-jetyy(npjet))**2.d0+(jetzz(npjet-2)-jetzz(npjet))**2.d0)
    end select
    if(tempmod>=dresolution)then
      call compute_posnoinserted(jetxx,jetyy,jetzz)
      linserted=.true.
      jetst(npjet-1)=0.d0!jetst(npjet-2)
      jetvx(npjet-1)=(jetvx(npjet-2)-ivelocity)/2.d0+ivelocity
      jetvy(npjet-1)=jetvy(npjet-2)/2.d0
      jetvz(npjet-1)=jetvz(npjet-2)/2.d0
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
       (jetyy(npjet-1)-jetyy(npjet))**2.d0+(jetzz(npjet-1)-jetzz(npjet))**2.d0)
    end select
    if(tempmod0>=thresolution)then
      npjet=npjet+1
      if(npjet>mxnpjet)call reallocate_jet()
      jetxx(npjet)=jetxx(npjet-1)
      jetyy(npjet)=jetyy(npjet-1)
      jetzz(npjet)=jetzz(npjet-1)
      jetst(npjet)=jetst(npjet-1)
      jetvx(npjet)=jetvx(npjet-1)
      jetvy(npjet)=jetvy(npjet-1)
      jetvz(npjet)=jetvz(npjet-1)
      jetms(npjet)=jetms(npjet-1)
      jetch(npjet)=jetch(npjet-1)
      jetxx(npjet-1)=resolution
      jetyy(npjet-1)=0.d0
      jetzz(npjet-1)=0.d0
      jetst(npjet-1)=istress
      jetvx(npjet-1)=ivelocity
      jetvy(npjet-1)=0.d0
      jetvz(npjet-1)=0.d0
      call extract_actual_massdensity(actualimass)
      call extract_actual_chargedensity(actualicharge)
      jetms(npjet-1)=actualimass*vdumbbell 
      jetch(npjet-1)=actualicharge*vdumbbell
      ladd=.true.
      timedeposition=0.d0
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
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision :: drand,theta
  
  if(.not.liniperturb)return
  
  select case(systype)
  case(1:2)
    return
  case(3:4)
    call random_number(drand)
    theta=2.d0*Pi*drand
    theta=0.d0
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
 
 subroutine remove_jetbead(nstepsub,timesub,lrem)
 
!***********************************************************************
!     
!     JETSPIN subroutine for removing the last bead at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstepsub
  double precision, intent(in) :: timesub
  logical, intent(inout) :: lrem
  
  if(.not.lremove)return
  if(jetxx(inpjet)>h)then
    lrem=.true.
    inpjet=inpjet+1
  endif
  
  return
  
 end subroutine remove_jetbead
 
 subroutine erase_jetbead(nstepsub,timesub,ladd,lrem)
 
!***********************************************************************
!     
!     JETSPIN subroutine for erasing the data of the last bead 
!     at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstepsub
  double precision, intent(in) :: timesub
  logical, intent(inout) :: ladd,lrem
  
  if(lrem)then
    jetch(inpjet-1)=0.d0
    jetms(inpjet-1)=0.d0
    jetst(inpjet-1)=0.d0
    jetvx(inpjet-1)=0.d0
    jetvy(inpjet-1)=0.d0
    jetvz(inpjet-1)=0.d0
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
 
 subroutine extract_actual_massdensity(actualmass)
 
!***********************************************************************
!     
!     JETSPIN subroutine for extracting randomically the mass of 
!     the last inserted bead at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(out) :: actualmass
  
  logical, save :: lfirst=.true.
  double precision, save :: oldgaussn,corr,radcorr
  double precision :: actualgaussn,newgaussn
  
  if(lfirst)then
    lfirst=.false.
    if(lmassavariable)then
      if(lencorrmassa==0.d0)then
        corr=0.d0
        radcorr=1.d0
      else
!       we consider a spatial correlation
        corr=dexp(-1.d0*resolution/lencorrmassa)
        radcorr=dsqrt(1.d0-corr**2.d0)
      endif
      if(idrank==0)then
        actualmass=0.d0
!       negative mass are not permitted
        do while(actualmass<=0.d0)
          newgaussn=gauss()
          actualgaussn=newgaussn
          actualmass=imassa+imassadev*actualgaussn
        enddo
        oldgaussn=newgaussn
      endif
      call bcast_world_d(actualmass)
    else
      actualmass=imassa
    endif
  else
    if(lmassavariable)then
      if(idrank==0)then
        actualmass=0.d0
!       negative mass are not permitted  
        do while(actualmass<=0.d0)
          newgaussn=gauss()
          actualgaussn=corr*oldgaussn+radcorr*newgaussn
          actualmass=imassa+imassadev*actualgaussn
        enddo
        oldgaussn=newgaussn
      endif
      call bcast_world_d(actualmass)
    else
      actualmass=imassa
    endif
  endif
  
  return
  
 end subroutine extract_actual_massdensity
 
 subroutine extract_actual_chargedensity(actualcharge)
 
!***********************************************************************
!     
!     JETSPIN subroutine for extracting randomically the charge of 
!     the last inserted bead at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(out) :: actualcharge
  
  logical, save :: lfirst=.true.
  double precision, save :: oldgaussn,corr,radcorr
  double precision :: actualgaussn,newgaussn
  
  if(lfirst)then
    lfirst=.false.
    if(lmassavariable)then
      if(lencorrcharge==0.d0)then
        corr=0.d0
        radcorr=1.d0
      else
!       we consider a spatial correlation
        corr=dexp(-1.d0*resolution/lencorrcharge)
        radcorr=dsqrt(1.d0-corr**2.d0)
      endif
      if(idrank==0)then
        actualcharge=0.d0
!       negative charge are not permitted 
        do while(actualcharge<=0.d0)
          newgaussn=gauss()
          actualgaussn=newgaussn
          actualcharge=icharge+ichargedev*actualgaussn
        enddo
        oldgaussn=newgaussn
      endif
      call bcast_world_d(actualcharge)
    else
      actualcharge=icharge
    endif
  else
    if(lmassavariable)then
      if(idrank==0)then
        actualcharge=0.d0
!       negative charge are not permitted        
        do while(actualcharge<=0.d0)
          newgaussn=gauss()
          actualgaussn=corr*oldgaussn+radcorr*newgaussn
          actualcharge=icharge+ichargedev*actualgaussn
        enddo
        oldgaussn=newgaussn
      endif
      call bcast_world_d(actualcharge)
    else
      actualcharge=icharge
    endif
  endif
  
  return
  
 end subroutine extract_actual_chargedensity
 
 end module nanojet_mod


