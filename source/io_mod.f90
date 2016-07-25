
 module io_mod
 
!***********************************************************************
!     
!     JETSPIN module for input/output routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
 use version_mod
 use parse_module
 use error_mod
 use utility_mod,           only : write_fmtnumb,pi
 use nanojet_mod,           only : airdragamp,doreorder,tao,aird,airv,&
                             chargescale,consistency,findex,g,h,&
                             icharge,icrossec,ilength,&
                             imassa,inpjet,istress,ivelocity,&
                             lencorrmassa,lengthscale,&
                             massscale,mu,jetbr,npjet,pampl,pfreq,&
                             resolution,surfacet,systype,tstep,units,&
                             V0,velext,yieldstress,dcutoff,lairdrag,&
                             lairvel,&
                             lconsistency,lfindex,lgravity,lHBfluid,&
                             lilengt,liniperturb,linserting,listress,&
                             livelocity,lmassavariable,lmyseed,&
                             ldcutoff,lnpjet,lremove,lresolution,&
                             lxyzrescale,lyieldstress,myseed,&
                             xyzrescale,laird,lairdragamp,lairv,lg,lh,&
                             licharge,licrossec,limass,&
                             llencorrmassa,lKVfluid,&
                             lmu,lpampl,lpfreq,lsurfacet,lsystype,&
                             ltstep,lunits,lv0,att,fve,gr,hg,ks,li,lrg,&
                             q,v,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,&
                             jetst,jetms,jetch,jetvl,jetlb,fvere,ksre, &
                             lire,ltrackbeads,lreordertrack, &
                             naddtrack,nremtrack,linserted, &
                             allocate_jet,deallocate_jet,mxnpjet, &
                             luppot,kuppot,lmirror, &
                             lenprobmassa,massratio,ldragvel, &
                             typedragvel,jetbd,ltagbeads,ldevelopers, &
                             lenthresholdbead,lbreakup,lfieldtype, &
                             fieldfreq,lfieldfreq,ultimatestrength,& 
                             lultimatestrength,lpdbrescale,pdbrescale, &
                             nmulstep,lmultiplestep,lmultisteperror, &
                             lmaxdispl,maxdispl,taoelectr,ltaoelectr, &
                             lpdbtagbeads,nmulstepdone,nmultisteperror,&
                             multisteperror,lreadrest,oldgaussn,corr, &
                             radcorr,lfirstmass,lneighlistdo,jetfr, &
                             jetfm,BLor,lflorentz,KLor,lmagneticfield
 use dynamic_refinement_mod, only : lrefinement,lrefinementthreshold,&
                             refinementthreshold,lrefinementevery, &
                             irefinementevery,lrefinementstart, &
                             irefinementstart,lrefbeadstart, &
                             refbeadstartfit,irefinementdone, &
                             llenthresholdbead
 use electric_field_mod,    only : nfieldtype
 use integrator_mod,        only : integrator,endtime,lendtime,&
                             lintegrator
 use statistic_mod,         only : nmaxstatdata,reprinttime,statdata,&
                             compute_statistic,addedmass,removedmass, &
                             addedcharge,removedcharge,countericurr, &
                             counterecurr,meanicurr,meanecurr, &
                             counterimass,counteremass,meanimass, &
                             meanemass,counterevel,counterevelrel, &
                             meanevel,meanevelrel,counterelen, &
                             counterecross,meanelen,meanecross, &
                             meancputime,counterivel,meanivel, &
                             counterlpath,meanlpath,ncounterevel, &
                             ncounterevelrel,ncountergeom,reprinttime, &
                             ncounterivel,ncounterlpath

 implicit none

 private
 
 integer, parameter :: maxlen=150
 
 integer, public, save :: nprintlist
 integer, public, save, allocatable, dimension(:) :: printlist
 integer, public, save :: nprintlist2
 integer, public, save, allocatable, dimension(:) :: printlist2
 logical, public, save :: lprintlist=.false.
 logical, public, save :: lprintlist2=.false.
 logical, public, save :: lprintxyz=.false.
 logical, public, save :: lprintxyzsing=.false.
 logical, public, save :: lprintpdbsing=.false.
 logical, public, save :: lprintdat=.false.
 logical, public, save :: lprintstatdat=.false.
 logical, save :: lnmulstep=.false.
 double precision, save :: printtime
 logical, public, save :: lprinttime=.false.
 logical, public, save :: lsprintdat=.false.
 logical, public, save :: lprintdatrem=.false.
 integer, public, save :: iprinttime=0
 integer, public, save :: eprinttime=0
 integer, public, save :: iprintxyz=0
 integer, public, save :: iprintxyzsing=0
 integer, public, save :: iprintpdbsing=0
 integer, public, save :: iprintdat=0
 integer, public, save :: sprintdat=1
 integer, public, save :: maxnumxyz=100
 double precision, save :: printxyz
 double precision, save :: printxyzsing
 double precision, save :: printpdbsing
 double precision, save :: printdat
 double precision, save :: refinementevery
 double precision, save :: refinementstart=0.d0
 double precision, public, allocatable, save :: xprint(:),xprint2(:)
 character(len=11),save :: namefile
 
 character(len=20) , public, allocatable :: printarg(:)
 character(len=20) , public, allocatable :: printarg2(:)
 
 logical :: lrestartreset=.false.
 logical :: ldragvelfound=.false.
 
 public :: print_logo
 public :: allocate_print,outprint_driver,read_input
 public :: print_internal_units
 public :: print_adimensional_parameters
 public :: print_adimensional_groups
 public :: print_legend_observables
 public :: finish_print
 public :: open_dat_file
 public :: write_dat_parameter
 public :: write_dat_frame
 public :: close_dat_file
 public :: open_datrem_file
 public :: write_datrem_frame
 public :: close_datrem_file
 public :: open_xyz_file
 public :: write_xyz_frame
 public :: close_xyz_file
 public :: write_xyz_singlefile
 public :: write_pdb_singlefile
 public :: write_restart_file
 public :: read_restart_file
 
 contains
 
 subroutine allocate_print()
  
  implicit none
  
  namefile='statout.dat'
  
  allocate(xprint(1:nprintlist))
  allocate(xprint2(1:nprintlist2))
  
  return
  
 end subroutine allocate_print
 
 subroutine print_logo(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the logo
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  character(len=*),parameter :: of='(a)'
  character(len=50) :: stringversion
  
  if(idrank/=0)return
  
  call print_version(stringversion)
  
  write(iu,*)
  write(iu,*)
  write(iu,of)"*******************************************************************************"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                 O U T P U T                                 *"
  write(iu,of)"*                                    from                                     *"
  write(iu,of)"*            __   _______ .___________.    _______..______    __  .__   __.   *"
  write(iu,of)"*          |  | |   ____||           |   /       ||   _  \  |  | |  \ |  |    *"
  write(iu,of)"*          |  | |  |__   `---|  |----`  |   (----`|  |_)  | |  | |   \|  |    *"
  write(iu,of)"*    .--.  |  | |   __|      |  |        \   \    |   ___/  |  | |  . `  |    *"
  write(iu,of)"*    |  `--'  | |  |____     |  |    .----)   |   |  |      |  | |  |\   |    *"
  write(iu,of)"*     \______/  |_______|    |__|    |_______/    | _|      |__| |__| \__|    *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*******************************************************************************"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*         =========================================================           *"
  write(iu,of)"*                        Program for Electrospinning                          *"
  write(iu,of)"*                         Simulations of Nanofibers                           *"
  write(iu,of)"*         =========================================================           *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Version 1.21 (July 2016)                                                 *"
  if(ldevelopers) &
  write(iu,of)"*    Compiled in developer mode                                               *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Written by:                                                              *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Marco Lauricella         IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    with contributions from:                                                 *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Ivan Coluzza             University of Vienna              Austria       *"
  write(iu,of)"*    Dario Pisignano          University of Salento             Italy         *"
  write(iu,of)"*    Giuseppe Pontrelli       IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*    Sauro Succi              IAC-CNR, Rome                     Italy         *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    This is an experimental code. The authors accept no responsibility       *"
  write(iu,of)"*    for the performance of the code or for the correctness of the results.   *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    The code is licensed under Open Software License v. 3.0 (OSL-3.0).       *"
  write(iu,of)"*    The full text of the licence can be found on the website:                *"
  write(iu,of)"*    http://opensource.org/licenses/OSL-3.0                                   *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    A brief explanation of this license is available on the website:         *"
  write(iu,of)"*    http://rosenlaw.com/OSL3.0-explained.htm                                 *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    The software development process has received funding from the           *"
  write(iu,of)"*    European Research Council under the European Union's Seventh Framework   *"
  write(iu,of)"*    Programme (FP/2007-2013)/ERC Grant Agreement n. 306357 (NANO-JETS).      *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    If results obtained with this code are published, an                     *"
  write(iu,of)"*    appropriate citation would be:                                           *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*    Marco Lauricella, Giuseppe Pontrelli, Ivan Coluzza,                      *"
  write(iu,of)"*    Dario Pisignano, Sauro Succi,                                            *"
  write(iu,of)"*    JETSPIN: A specific-purpose open-source software for                     *"
  write(iu,of)"*    electrospinning simulations of nanofibers,                               *"
  write(iu,of)"*    Computer Physics Communications, 197 (2015), pp. 227-238.                *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*******************************************************************************"
  write(iu,of)"*                                                                             *"
  write(iu,'(a5,a50,a24)')"*    ",stringversion,"                       *"
  write(iu,of)"*                                                                             *"
  write(iu,of)"*******************************************************************************"
  write(iu,of)"                                                                               "
  
  return
  
 end subroutine print_logo
 
 subroutine print_adimensional_parameters(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the dimensionless parameters
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  
  character(len=10) :: r_char
  
  if(idrank/=0)return
  
  write(iu,*)
  write(iu,'(a)')"Adimensional parameters:"
  write(iu,'(a,g20.10,a,g20.10,a)')"timestep                   = ", &
   tstep," = ",tstep*tao," s"
  write(iu,'(a,g20.10,a,g20.10,a)')"final time of integration  = ", &
   endtime," = ",endtime*tao," s"
  write(iu,'(a,g20.10,a,g20.10,a)')"resolution                 = ", &
   resolution," = ",resolution*lengthscale," cm"
  write(iu,'(a,g20.10,a,g20.10,a)')"initial length             = ", &
   ilength," = ",ilength*lengthscale," cm"
  if(listress)then
  write(iu,'(a,g20.10,a,g20.10,a)')"nozzle stress              = ", &
   istress," = ",istress*G," g cm^-1 s^-2"
  endif
  if(livelocity)then
  write(iu,'(a,g20.10,a,g20.10,a)')"nozzle velocity            = ", &
   ivelocity," = ",ivelocity*(tao/lengthscale)," cm s^-1"
  endif
  if(lyieldstress)then
  write(iu,'(a,g20.10,a,g20.10,a)')"yield stress               = ", &
   yieldstress," = ",yieldstress*G," g cm^-1 s^-2"
  endif
  write(iu,'(a,g20.10,a,g20.10,a)')"density mass               = ", &
   imassa," = ",imassa*massscale*(lengthscale**3.d0)," g cm^-3"
  write(iu,'(a,g20.10,a,g20.10,a)')"density charge             = ", &
   icharge," = ",icharge*chargescale*(lengthscale**3.d0)," statC cm^-3"
  write(iu,'(a,g20.10,a,g20.10,a)')"nozzle cross section       = ", &
   icrossec," = ",icrossec*lengthscale," cm"
  if(lairvel)then
  write(iu,'(a,g20.10,a,g20.10,a)')"airdrag airvelocity        = ", &
   velext," = ",velext*(lengthscale/tao)," cm s^-1"
  endif
  if(lairdragamp)then
  write(iu,'(a,g20.10,a,g20.10,a)')"airdrag amplitude          = ", &
   airdragamp(1)," = ",airdragamp(1)*(lengthscale**2.d0)/(tao**3.d0),&
    " cm^2 s^-3"
  endif
  if(lconsistency)then
  write (r_char,'(f10.2)')findex-2.d0
  write(iu,'(a,g20.10,a,g20.10,2a)')"hbfluid consistency        = ", &
   consistency," = ",consistency*(mu*((tao)**(findex-1.d0))),&
    " g cm^-1 s^",trim(adjustl(r_char))
  endif
  if(lmassavariable)then
  write(iu,'(a,g20.10,a,g20.10,a)')"variable mass distance NP  = ", &
   lencorrmassa," = ",lencorrmassa*lengthscale," cm"
  endif
  if(lrefinementthreshold)then
  write(iu,'(a,g20.10,a,g20.10,a)')"dynamic refinement threshold= ", &
   refinementthreshold," = ",refinementthreshold*lengthscale," cm"
  endif
  if(lmagneticfield)then
  write(iu,'(a,g20.10,a,g20.10,a)')"magnetic field along x    = ", &
   BLor(1)," = ",BLor(1)*dsqrt(massscale)/(dsqrt(lengthscale)*tao), &
   " cm^-0.5 g^0.5 s^-1"
  write(iu,'(a,g20.10,a,g20.10,a)')"magnetic field along y    = ", &
   BLor(2)," = ",BLor(2)*dsqrt(massscale)/(dsqrt(lengthscale)*tao), &
   " cm^-0.5 g^0.5 s^-1"
  write(iu,'(a,g20.10,a,g20.10,a)')"magnetic field along z    = ", &
   BLor(3)," = ",BLor(3)*dsqrt(massscale)/(dsqrt(lengthscale)*tao), &
   " cm^-0.5 g^0.5 s^-1"
  endif
  
  write(iu,*)
  
  return
  
 end subroutine print_adimensional_parameters
  
 subroutine print_adimensional_groups(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the dimensionless groups
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  
  if(idrank/=0)return
  
  write(iu,*)
  write(iu,'(a)')"Adimensional groups:"
  write(iu,'(a,g20.10)')"Q     = ",q
  write(iu,'(a,g20.10)')"Phi   = ",v
  write(iu,'(a,g20.10,a,g20.10,a)')"Fve   = ",fve," = ",fvere, &
   " in Reneker model"
  write(iu,'(a,g20.10)')"H     = ",Hg
  write(iu,'(a,g20.10)')"Lstep = ",Lrg
  if(lflorentz)then
    write(iu,'(a,g20.10)')"KLor = ",KLor
  endif
  if(lgravity)then
    write(iu,'(a,g20.10)')"Fg    = ",Gr
  endif
  if(systype/=1)then
  write(iu,'(a,g20.10,a,g20.10,a)')"A     = ",ks," = ",ksre, &
   " in Reneker model"
  endif
  if(liniperturb)then
  write(iu,'(a,g20.10)')"Ks    = ",pfreq
  endif
  if(lairdrag)then
  write(iu,'(a,g20.10)')"Gamma = ",att
  write(iu,'(a,g20.10,a,g20.10,a)')"Lambda= ",Li," = ",lire, &
   " in Reneker model"
  write(iu,'(a,g20.10)')"Theta = ",airdragamp(1)
  endif
  write(iu,*)
  
  return
  
 end subroutine print_adimensional_groups
 
 subroutine print_legend_observables(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the legend of the observables
!     requested in input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  
  integer :: i
  
  if(idrank/=0)return
  
  write(iu,*)
  write(iu,'(a)')"Printed observables as specified in input file:"
  do i=1,nprintlist
    write(iu,'(a)')legendobs(i,printlist)
  enddo
  write(iu,*)
  write(iu,*)  
  write(iu,'(a)')'Start simulation'
  
  return
  
 end subroutine print_legend_observables
 
 function legendobs(iarg,printcodsub)
 
!***********************************************************************
!     
!     JETSPIN function for returning the legend which is associated to
!     the integer contained in the printcod array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iarg
  integer, intent(in), allocatable, dimension(:) :: printcodsub
  
  character(len=48) :: legendobs
  
  legendobs=repeat(' ',48)
  if(printcodsub(iarg)==1)then
    legendobs='t    =  unscaled time                            '
  elseif(printcodsub(iarg)==24)then
    legendobs='ts   =  scaled time                              '
  elseif(printcodsub(iarg)==2)then
    legendobs='x    =  x unscaled coordinates                   '
  elseif(printcodsub(iarg)==3)then
    legendobs='y    =  y unscaled coordinates                   '                 
  elseif(printcodsub(iarg)==4)then
    legendobs='z    =  z unscaled coordinates                   '
  elseif(printcodsub(iarg)==25)then
    legendobs='xs   =  x scaled coordinates                     '
  elseif(printcodsub(iarg)==26)then
    legendobs='ys   =  y scaled coordinates                     '
  elseif(printcodsub(iarg)==27)then
    legendobs='zs   =  z scaled coordinates                     '
  elseif(printcodsub(iarg)==5)then
    legendobs='st   =  unscaled stress                          '
  elseif(printcodsub(iarg)==28)then
    legendobs='sts  =  scaled stress                            '
  elseif(printcodsub(iarg)==6)then
    legendobs='vx   =  x unscaled velocities                    '
  elseif(printcodsub(iarg)==7)then
    legendobs='vy   =  y unscaled velocities                    '
  elseif(printcodsub(iarg)==8)then
    legendobs='vz   =  z unscaled velocities                    '
  elseif(printcodsub(iarg)==29)then
    legendobs='vxs  =  x  scaled velocities                     '
  elseif(printcodsub(iarg)==30)then
    legendobs='vys  =  y  scaled velocities                     '
  elseif(printcodsub(iarg)==31)then
    legendobs='vzs  =  z  scaled velocities                     '
  elseif(printcodsub(iarg)==9)then
    legendobs='yz   =  unscaled normal distance from the x axis '
  elseif(printcodsub(iarg)==32)then
    legendobs='yzs  =  scaled normal distance from the x axis   '
  elseif(printcodsub(iarg)==34)then
    legendobs='mass =  last inserted mass at the nozzle         '
  elseif(printcodsub(iarg)==35)then
    legendobs='q    =  last inserted charge at the nozzle       '
  elseif(printcodsub(iarg)==22)then
    legendobs='cpu  =  time for every print interval            '
  elseif(printcodsub(iarg)==23)then
    legendobs='cpur =  remaining time to the end                '
  elseif(printcodsub(iarg)==33)then
    legendobs='cpue =  elapsed time                             '
  elseif(printcodsub(iarg)==10)then
    legendobs='n    =  bead number of jet discretization        '
  elseif(printcodsub(iarg)==11)then
    legendobs='f    =  index of the first bead                  '
  elseif(printcodsub(iarg)==12)then
    legendobs='l    =  index of the last bead                   '
  elseif(printcodsub(iarg)==13)then
    legendobs='curn =  current at the nozzle                    '
  elseif(printcodsub(iarg)==14)then
    legendobs='curc =  current at the collector                 '
  elseif(printcodsub(iarg)==17)then
    legendobs='vc   =  velocity modulus of jet at the collector '
  elseif(printcodsub(iarg)==18)then
    legendobs='svc  =  strain velocity at the collector         '
  elseif(printcodsub(iarg)==15)then
    legendobs='mfn  =  mass flux at the nozzle                  '
  elseif(printcodsub(iarg)==16)then
    legendobs='mfc  =  mass flux at the collector               '
  elseif(printcodsub(iarg)==20)then
    legendobs='rc   =  radius of jet  at the collector          '
  elseif(printcodsub(iarg)==21)then
    legendobs='rrr  =  radius reduction ratio of jet            '
  elseif(printcodsub(iarg)==36)then
    legendobs='vn   =  velocity modulus of jet at the nozzle    '
  elseif(printcodsub(iarg)==37)then
    legendobs='lp   =  length path of jet                       '
  elseif(printcodsub(iarg)==38)then
    legendobs='rlp  =  length path of jet / collector distance  '
  elseif(printcodsub(iarg)==39)then
    legendobs='nref =  number of performed dynamic refinements  '
  elseif(printcodsub(iarg)==19)then
    legendobs='mlc  =  mutual length segment at the collector   '
  elseif(printcodsub(iarg)==40)then
    legendobs='angl =  instantaneus angular aperture            '
  elseif(printcodsub(iarg)==41)then
    legendobs='mxst =  maximum stress along the jet             '
  elseif(printcodsub(iarg)==42)then
    legendobs='mxsx =  x position of the jet maximum stress     '
  elseif(printcodsub(iarg)==43)then
    legendobs='nms  =  number of multiple step procedure done   '
  elseif(printcodsub(iarg)==44)then
    legendobs='erms =  maximum error of multiple step approach  '
  elseif(printcodsub(iarg)==45)then
    legendobs='v    =  value of the external electric potential '
  endif
  legendobs=adjustl(legendobs)
  
  return
  
 end function legendobs
 
 subroutine print_internal_units(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the internat units
!     which are used in JETSPIN
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: iu
  
  if(idrank/=0)return
  
  write(iu,*)
  write(iu,'(a)')"Internal units:"
  write(iu,'(a,g20.10,a)')"Time unit     = ",tao,' s'
  write(iu,'(a,g20.10,a)')"Length unit   = ",lengthscale,' cm'
  write(iu,'(a,g20.10,a)')"Stress unit   = ",G,' g cm^-1 s^-2'
  write(iu,'(a,g20.10,a)')"Velocity unit = ",lengthscale/tao,' cm s^-1'
  write(iu,'(a,g20.10,a)')"Mass unit     = ",massscale,' g'
  write(iu,'(a,g20.10,a)')"Charge unit   = ",chargescale,' statC'
  write(iu,*)
  
  return
  
 end subroutine print_internal_units
 
 subroutine read_input(inputunit,inputname)
 
!***********************************************************************
!     
!     JETSPIN subroutine for reading the input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: inputunit
  character(len=*), intent(in) :: inputname
  
  character(len=maxlen) :: redstring,directive
  logical :: safe,lredo,ltest,lprintlisterror,lprintlisterror2
  logical :: ltestread,lexists,lfoundprint
  integer :: inumchar,i,nwords,iline,itest
  character(len=maxlen) ,allocatable :: outwords(:),outwords2(:)
  
! initialize parameters  
  systype=0      
  integrator=0   
  nprintlist=0
  nprintlist2=0
  nfieldtype=0
  ilength=0.d-0
  istress=0.d0 
  ivelocity=0.d0 
  icrossec=0.d0 
  resolution=0.d0
  imassa=0.d0 
  icharge=0.d0 
  mu=0.0d0 
  G=0.d0 
  h=0.d0 
  V0=0.d0 
  surfacet=0.d0 
  tstep=0.d0
  pfreq=0.d0
  pampl=0.d0
  myseed=1
  sprintdat=1
  aird=0.d0
  airv=1.d0
  velext=0.d0
  findex=1.d0
  consistency=1.d0 
  refinementthreshold=0.d0
  xyzrescale=1.d0
  pdbrescale=0.d0
  refinementstart=0.d0
  kuppot=0.d0
  fieldfreq=0.d0
  ultimatestrength=0.d0
  
  lsystype=.false.
  lintegrator=.false.
  lprintlist=.false.
  lprintlist2=.false.
  lresolution=.false.
  lilengt=.false.
  limass=.false.
  licharge=.false.
  licrossec=.false.
  listress=.false.
  livelocity=.false.
  lmu=.false.
  lG=.false.
  lh=.false.
  lV0=.false.
  lsurfacet=.false.
  ltstep=.false.
  lendtime=.false.
  lprinttime=.false.
  lnpjet=.false.
  lunits=.false.
  lyieldstress=.false.
  linserting=.false.
  liniperturb=.false.
  lpfreq=.false.
  lpampl=.false.
  ldcutoff=.false.
  lairdrag=.false.
  lairdragamp=.false.
  lremove=.false.
  lmyseed=.false.
  lprintxyz=.false.
  lprintxyzsing=.false.
  lprintpdbsing=.false.
  lxyzrescale=.false.
  lpdbrescale=.false.
  lpdbtagbeads=.false.
  lprintdat=.false.
  laird=.false.
  lairv=.false.
  lairvel=.false.
  lHBfluid=.false.
  lKVfluid=.false.
  lgravity=.false.
  lprintstatdat=.false.
  
  lprintlisterror=.false.
  lprintlisterror2=.false.
  lfoundprint=.false.
  lsprintdat=.false.
  lprintdatrem=.false.
  
  ltestread=.false.
  
  lrefinement=.false.
  lrefinementthreshold=.false.
  
  lrefbeadstart=.false.
  luppot=.false.
  lmirror=.false.
  ltagbeads=.false.
  
  ldragvelfound=.false.
  
  lbreakup=.false.
  lfieldtype=.false.
  lfieldfreq=.false.
  lultimatestrength=.false.
  lnmulstep=.false.
  lmultisteperror=.false.
  
  llenthresholdbead=.false.
  lflorentz=.false.
  lmagneticfield=.false.
  
! note the parameters are read only by the zero node
  if(idrank==0)then
    
!   check if the input file exist
    inquire(file=inputname,exist=lexists)
    if(.not.lexists)call error(10)
    
!   open the inout file
    open(unit=inputunit,file=inputname,status='old',action='read', &
    iostat=itest)
    if(itest/=0)call error(11)
    
!   initialize the lredo condition
    lredo=.true.
    
!   counter the line which are read in the input file
    iline=0
    
!   read the input file as long as the finish directive is read
    do while(lredo)
      call getline(safe,inputunit,maxlen,redstring)
      if(.not.safe)call error(5)
      iline=iline+1
      call strip(redstring,maxlen)
      call lowcase(redstring,maxlen)
      call copystring(redstring,directive,maxlen)
      
      if(redstring(1:1)=='#')then
        cycle
      elseif(redstring(1:1)=='!')then
        cycle
      elseif(redstring(1:1)==' ')then
        cycle
      elseif(findstring('primary',directive,inumchar,maxlen))then
        if(findstring('cutoff',directive,inumchar,maxlen))then
          ldcutoff=.true.
          dcutoff=dblstr(directive,maxlen,inumchar)
        endif
      elseif(findstring('maximum',directive,inumchar,maxlen))then
        if(findstring('displ',directive,inumchar,maxlen))then
          lmaxdispl=.true.
          maxdispl=dblstr(directive,maxlen,inumchar)
        endif
      elseif(findstring('multiple',directive,inumchar,maxlen))then
        if(findstring('step',directive,inumchar,maxlen))then
          if(findstring('yes',directive,inumchar,maxlen))then
            lmultiplestep=.true.
          elseif(findstring('every',directive,inumchar,maxlen))then
            lnmulstep=.true.
            nmulstep=intstr(directive,maxlen,inumchar)
          endif
        endif
      elseif(findstring('system',directive,inumchar,maxlen))then
        systype=intstr(directive,maxlen,inumchar)
        lsystype=.true.
      elseif(findstring('integrator',directive,inumchar,maxlen))then
        integrator=intstr(directive,maxlen,inumchar)
        lintegrator=.true.
      elseif(findstring('units',directive,inumchar,maxlen))then
        lunits=.true.
        if(findstring('cgs',directive,inumchar,maxlen))then
          units=1
        elseif(findstring('real',directive,inumchar,maxlen))then
          units=2
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('seed',directive,inumchar,maxlen))then
        myseed=intstr(directive,maxlen,inumchar)
        lmyseed=.true.
      elseif(findstring('resolution',directive,inumchar,maxlen))then
        resolution=dblstr(directive,maxlen,inumchar)
        lresolution=.true.
      elseif(findstring('points',directive,inumchar,maxlen))then
        npjet=intstr(directive,maxlen,inumchar)
        lnpjet=.true.
      elseif(findstring('inserting',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          linserting=.true.
        elseif(findstring('no',directive,inumchar,maxlen))then
          linserting=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('removing',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          lremove=.true.
        elseif(findstring('no',directive,inumchar,maxlen))then
          lremove=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('lorentz',directive,inumchar,maxlen).and. &
       ldevelopers)then
        if(findstring('yes',directive,inumchar,maxlen))then
          lflorentz=.true.
        else
          lflorentz=.false.
        endif
      elseif(findstring('magnetic',directive,inumchar,maxlen).and. &
       ldevelopers)then
        if(findstring('field',directive,inumchar,maxlen))then
          lmagneticfield=.true.
          BLor(1)=dblstr(directive,maxlen,inumchar)
          BLor(2)=dblstr(directive,maxlen,inumchar)
          BLor(3)=dblstr(directive,maxlen,inumchar)
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('gravity',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          lgravity=.true.
        elseif(findstring('no',directive,inumchar,maxlen))then
          lgravity=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('density',directive,inumchar,maxlen))then
        if(findstring('mass',directive,inumchar,maxlen))then
          imassa=dblstr(directive,maxlen,inumchar)
          limass=.true.
        elseif(findstring('charge',directive,inumchar,maxlen))then
          icharge=dblstr(directive,maxlen,inumchar)
          licharge=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('initial',directive,inumchar,maxlen))then
        if(findstring('length',directive,inumchar,maxlen))then
          ilength=dblstr(directive,maxlen,inumchar)
          lilengt=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('nozzle',directive,inumchar,maxlen))then
        if(findstring('cross',directive,inumchar,maxlen))then
          icrossec=dblstr(directive,maxlen,inumchar)
          licrossec=.true.
        elseif(findstring('stress',directive,inumchar,maxlen))then
          istress=dblstr(directive,maxlen,inumchar)
          listress=.true.
        elseif(findstring('velocity',directive,inumchar,maxlen))then
          ivelocity=dblstr(directive,maxlen,inumchar)
          livelocity=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('perturb',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          liniperturb=.true.
        elseif(findstring('freq',directive,inumchar,maxlen))then
          lpfreq=.true.
          pfreq=dblstr(directive,maxlen,inumchar)
        elseif(findstring('ampl',directive,inumchar,maxlen))then
          lpampl=.true.
          pampl=dblstr(directive,maxlen,inumchar)
        elseif(findstring('no',directive,inumchar,maxlen))then
          liniperturb=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('airdrag',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          lairdrag=.true.
        elseif(findstring('amplitude',directive,inumchar,maxlen))then
          lairdragamp=.true.
          airdragamp(1)=dblstr(directive,maxlen,inumchar)
          airdragamp(2)=airdragamp(1)
          airdragamp(3)=airdragamp(2)
        elseif(findstring('airdensity',directive,inumchar,maxlen))then
          laird=.true.
          aird=dblstr(directive,maxlen,inumchar)
        elseif(findstring('airviscosity',directive,inumchar,maxlen))then
          lairv=.true.
          airv=dblstr(directive,maxlen,inumchar)
        elseif(findstring('airvelocity',directive,inumchar,maxlen))then
          lairvel=.true.
          velext=dblstr(directive,maxlen,inumchar)
        elseif(findstring('no',directive,inumchar,maxlen))then
          lairdrag=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('viscosity',directive,inumchar,maxlen))then
        mu=dblstr(directive,maxlen,inumchar)
        lmu=.true.
      elseif(findstring('elastic',directive,inumchar,maxlen))then
        if(findstring('modulus',directive,inumchar,maxlen))then
          G=dblstr(directive,maxlen,inumchar)
          lG=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('yield',directive,inumchar,maxlen))then
        if(findstring('stress',directive,inumchar,maxlen))then
          yieldstress=dblstr(directive,maxlen,inumchar)
          lyieldstress=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('collector',directive,inumchar,maxlen))then
        if(findstring('distance',directive,inumchar,maxlen))then
          h=dblstr(directive,maxlen,inumchar)
          lh=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('external',directive,inumchar,maxlen))then
        if(findstring('potential',directive,inumchar,maxlen))then  
          if(findstring('type',directive,inumchar,maxlen))then
            nfieldtype=intstr(directive,maxlen,inumchar)
            lfieldtype=.true.
          elseif(findstring('freq',directive,inumchar,maxlen))then
            fieldfreq=dblstr(directive,maxlen,inumchar)
            lfieldfreq=.true.
          elseif(findstring('time',directive,inumchar,maxlen))then
            taoelectr=dblstr(directive,maxlen,inumchar)
            ltaoelectr=.true.
          else
            V0=dblstr(directive,maxlen,inumchar)
            lV0=.true.
          endif
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('surface',directive,inumchar,maxlen))then
        if(findstring('tension',directive,inumchar,maxlen))then
          lsurfacet=.true.
          surfacet=dblstr(directive,maxlen,inumchar)
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('timestep',directive,inumchar,maxlen))then
        tstep=dblstr(directive,maxlen,inumchar)
        ltstep=.true.
      elseif(findstring('final',directive,inumchar,maxlen))then
        if(findstring('time',directive,inumchar,maxlen))then
          endtime=dblstr(directive,maxlen,inumchar)
          lendtime=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('printstat',directive,inumchar,maxlen))then
        if(findstring('list',directive,inumchar,maxlen))then
          lprintlist2=.true.
          lprintlisterror2=.false.
          call findwords(nwords,outwords,directive,maxlen)
          nprintlist2=max(0,nwords-2)
          directive=outwords(1)
          if(.not.findstring('printstat',directive,inumchar,maxlen))then
            lprintlisterror2=.true.
          endif
          directive=outwords(2)
          if(.not.findstring('list',directive,inumchar,maxlen))then
            lprintlisterror2=.true.
          endif
          if(nprintlist2>0)allocate(printlist2(nprintlist2))
          if(allocated(outwords2))deallocate(outwords2)
          if(nprintlist2>0)allocate(outwords2(nprintlist2))
          do i=1,nprintlist2
            outwords2(i)=outwords(i+2)
            call identify_argument(i,outwords2,printlist2,maxlen, &
             lfoundprint)
            lprintlisterror2=(lprintlisterror2 .or. (.not.lfoundprint))
          enddo
        elseif(findstring('binary',directive,inumchar,maxlen))then
          lprintstatdat=.true.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('print',directive,inumchar,maxlen))then
        if(findstring('time',directive,inumchar,maxlen))then
          printtime=dblstr(directive,maxlen,inumchar)
          lprinttime=.true.
        elseif(findstring('xyz',directive,inumchar,maxlen))then
          if(findstring('maxnum',directive,inumchar,maxlen))then
            maxnumxyz=intstr(directive,maxlen,inumchar)
          elseif(findstring('frame',directive,inumchar,maxlen))then
            if(findstring('rescalexyz',directive,inumchar,maxlen))then
              call warning(61,dble(iline),redstring)
              ltestread=.true.
            endif
            printxyzsing=dblstr(directive,maxlen,inumchar)
            lprintxyzsing=.true.
          elseif(findstring('rescalexyz',directive,inumchar,maxlen))then
            if(findstring('frame',directive,inumchar,maxlen))then
              call warning(61,dble(iline),redstring)
              ltestread=.true.
            endif
            lxyzrescale=.true.
            xyzrescale=dblstr(directive,maxlen,inumchar)
          else
            printxyz=dblstr(directive,maxlen,inumchar)
            lprintxyz=.true.
          endif
        elseif(findstring('pdb',directive,inumchar,maxlen))then
          if(findstring('frame',directive,inumchar,maxlen))then
            if(findstring('rescalepdb',directive,inumchar,maxlen))then
              call warning(61,dble(iline),redstring)
              ltestread=.true.
            endif
            printpdbsing=dblstr(directive,maxlen,inumchar)
            lprintpdbsing=.true.
          elseif(findstring('rescalepdb',directive,inumchar,maxlen))then
            if(findstring('frame',directive,inumchar,maxlen))then
              call warning(61,dble(iline),redstring)
              ltestread=.true.
            endif
            lpdbrescale=.true.
            pdbrescale=dblstr(directive,maxlen,inumchar)
          elseif(findstring('tagbead',directive,inumchar,maxlen) .and. &
           ldevelopers)then
            lpdbtagbeads=.true.
          else
            call warning(61,dble(iline),redstring)
            ltestread=.true.
          endif
        elseif(findstring('binary',directive,inumchar,maxlen))then
          if(findstring('removed',directive,inumchar,maxlen))then
            lprintdatrem=.true.
          else
            lprintdat=.true.
            printdat=dblstr(directive,maxlen,inumchar)
            if(findstring('style',directive,inumchar,maxlen))then
              lsprintdat=.true.
              sprintdat=intstr(directive,maxlen,inumchar)
            endif
          endif
        elseif(findstring('list',directive,inumchar,maxlen))then
          lprintlist=.true.
          lprintlisterror=.false.
          call findwords(nwords,outwords,directive,maxlen)
          nprintlist=max(0,nwords-2)
          directive=outwords(1)
          if(.not.findstring('print',directive,inumchar,maxlen))then
            lprintlisterror=.true.
          endif
          directive=outwords(2)
          if(.not.findstring('list',directive,inumchar,maxlen))then
            lprintlisterror=.true.
          endif
          if(nprintlist>0)allocate(printlist(nprintlist))
          if(allocated(outwords2))deallocate(outwords2)
          if(nprintlist>0)allocate(outwords2(nprintlist))
          do i=1,nprintlist
            outwords2(i)=outwords(i+2)
            call identify_argument(i,outwords2,printlist,maxlen, &
             lfoundprint)
            lprintlisterror=(lprintlisterror .or. (.not.lfoundprint))
          enddo
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('hbfluid',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          lHBfluid=.true.
        elseif(findstring('consistency',directive,inumchar,maxlen))then
          lconsistency=.true.
          consistency=dblstr(directive,maxlen,inumchar)
        elseif(findstring('index',directive,inumchar,maxlen))then
          lfindex=.true.
          findex=dblstr(directive,maxlen,inumchar)
        elseif(findstring('no',directive,inumchar,maxlen))then
          lHBfluid=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('kvfluid',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          lKVfluid=.true.
        elseif(findstring('no',directive,inumchar,maxlen))then
          lKVfluid=.false.
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('variable',directive,inumchar,maxlen) .and. &
       ldevelopers)then
        if(findstring('mass',directive,inumchar,maxlen))then
          if(findstring('yes',directive,inumchar,maxlen))then
            lmassavariable=.true.
            ltagbeads=.true.
          elseif(findstring('npdist',directive,inumchar,maxlen) &
           )then
            llencorrmassa=.true.
            lencorrmassa=dblstr(directive,maxlen,inumchar)
          elseif(findstring('no',directive,inumchar,maxlen))then
            lmassavariable=.false.
          elseif(findstring('npradius',directive,inumchar,maxlen))then
            lenprobmassa=dblstr(directive,maxlen,inumchar)
          elseif(findstring('npmass',directive,inumchar,maxlen))then
            massratio=dblstr(directive,maxlen,inumchar)
          else
            call warning(61,dble(iline),redstring)
          endif
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('dynam',directive,inumchar,maxlen))then
        if(findstring('refin',directive,inumchar,maxlen))then
          if(findstring('yes',directive,inumchar,maxlen))then
            lrefinement=.true.
            ltagbeads=.true.
          elseif(findstring('threshold',directive,inumchar,maxlen))then
            lrefinementthreshold=.true.
            refinementthreshold=dblstr(directive,maxlen,inumchar)
            lrefbeadstart=.true.
            refbeadstartfit=refinementthreshold
          elseif(findstring('anchor',directive,inumchar,maxlen))then
            llenthresholdbead=.true.
            lenthresholdbead=dblstr(directive,maxlen,inumchar)
          elseif(findstring('every',directive,inumchar,maxlen))then
            lrefinementevery=.true.
            refinementevery=dblstr(directive,maxlen,inumchar)
          elseif(findstring('start',directive,inumchar,maxlen))then
            lrefinementstart=.true.
            refinementstart=dblstr(directive,maxlen,inumchar)
          elseif(findstring('no',directive,inumchar,maxlen))then
            lrefinement=.false.
          else
            call warning(61,dble(iline),redstring)
            ltestread=.true.
          endif
        endif
      elseif(findstring('restart',directive,inumchar,maxlen))then
        if(findstring('yes',directive,inumchar,maxlen))then
          lreadrest=.true.
        elseif(findstring('reset',directive,inumchar,maxlen))then
          lrestartreset=.true.
        elseif(findstring('no',directive,inumchar,maxlen))then
          lrestartreset=.false.
        endif
      elseif(findstring('wall',directive,inumchar,maxlen) .and. &
       ldevelopers)then
        if(findstring('yes',directive,inumchar,maxlen))then
          luppot=.true.
          kuppot=dblstr(directive,maxlen,inumchar)
        endif
      elseif(findstring('dragvel',directive,inumchar,maxlen))then
        ldragvelfound=.true.
        if(findstring('yes',directive,inumchar,maxlen))then
          ldragvel=.true.
          if(findstring('type',directive,inumchar,maxlen) .and. &
           ldevelopers)then
            typedragvel=intstr(directive,maxlen,inumchar)
          endif
        elseif(findstring('no',directive,inumchar,maxlen))then
          ldragvel=.false.
        endif
      elseif(findstring('mirror',directive,inumchar,maxlen) .and. &
       ldevelopers)then
        if(findstring('yes',directive,inumchar,maxlen))then
          lmirror=.true.
        endif
      elseif(findstring('breakup',directive,inumchar,maxlen) .and. &
       ldevelopers)then
        if(findstring('yes',directive,inumchar,maxlen))then
          lbreakup=.true.
        endif
      elseif(findstring('ultimate',directive,inumchar,maxlen) .and. &
       ldevelopers)then
        if(findstring('strength',directive,inumchar,maxlen))then
          lultimatestrength=.true.
          ultimatestrength=dblstr(directive,maxlen,inumchar)
        else
          call warning(61,dble(iline),redstring)
          ltestread=.true.
        endif
      elseif(findstring('finish',directive,inumchar,maxlen))then
        lredo=.false.
      else
        call warning(37,1.d0,redstring)
        call error(6)
      endif
    enddo
    
    close(inputunit)
    
  endif
  
! send the read parameters to all the nodes
  call bcast_world_l(ldcutoff)
  call bcast_world_d(dcutoff)
  call bcast_world_i(systype)
  call bcast_world_l(lsystype)
  call bcast_world_i(integrator)
  call bcast_world_l(lintegrator)
  call bcast_world_l(lunits)
  call bcast_world_i(units)
  call bcast_world_i(nprintlist)
  call bcast_world_l(lprintlist)
  call bcast_world_l(lprintlisterror)
  call bcast_world_i(nprintlist2)
  call bcast_world_l(lprintlist2)
  call bcast_world_l(lprintlisterror2)
  call bcast_world_i(myseed)
  call bcast_world_l(lmyseed)
  call bcast_world_d(resolution)
  call bcast_world_l(lresolution)
  call bcast_world_i(npjet)
  call bcast_world_l(lnpjet)
  call bcast_world_l(linserting)
  call bcast_world_l(lremove)
  call bcast_world_l(lgravity)
  call bcast_world_d(ilength)
  call bcast_world_l(lilengt)
  call bcast_world_d(imassa)
  call bcast_world_l(limass)
  call bcast_world_d(icharge)
  call bcast_world_l(licharge)
  call bcast_world_d(icrossec)
  call bcast_world_l(licrossec)
  call bcast_world_d(istress)
  call bcast_world_l(listress)
  call bcast_world_d(ivelocity)
  call bcast_world_l(livelocity)
  call bcast_world_l(liniperturb)
  call bcast_world_l(lpfreq)
  call bcast_world_d(pfreq)
  call bcast_world_l(lpampl)
  call bcast_world_d(pampl)
  call bcast_world_l(lairdrag)
  call bcast_world_l(lairdragamp)
  call bcast_world_darr(airdragamp,3)
  call bcast_world_d(mu)
  call bcast_world_l(lmu)
  call bcast_world_d(G)
  call bcast_world_l(lG)
  call bcast_world_d(yieldstress)
  call bcast_world_l(lyieldstress)
  call bcast_world_d(h)
  call bcast_world_l(lh)
  call bcast_world_d(V0)
  call bcast_world_l(lV0)
  call bcast_world_l(lsurfacet)
  call bcast_world_d(surfacet)
  call bcast_world_d(tstep)
  call bcast_world_l(ltstep)
  call bcast_world_d(endtime)
  call bcast_world_l(lendtime)
  call bcast_world_d(printtime)
  call bcast_world_l(lprinttime)
  call bcast_world_i(maxnumxyz)
  call bcast_world_d(printxyz)
  call bcast_world_l(lprintxyz)
  call bcast_world_l(lprintxyzsing)
  call bcast_world_d(xyzrescale)
  call bcast_world_l(lxyzrescale)
  call bcast_world_l(lprintpdbsing)
  call bcast_world_d(pdbrescale)
  call bcast_world_l(lpdbrescale)
  call bcast_world_l(lpdbtagbeads)
  call bcast_world_l(lprintdat)
  call bcast_world_d(printdat)
  call bcast_world_i(sprintdat)
  call bcast_world_l(lsprintdat)
  call bcast_world_l(laird)
  call bcast_world_l(lairv)
  call bcast_world_l(lairvel)
  call bcast_world_d(aird)
  call bcast_world_d(airv)
  call bcast_world_d(velext)
  call bcast_world_l(lHBfluid)
  call bcast_world_l(lKVfluid)
  call bcast_world_l(lconsistency)
  call bcast_world_l(lfindex)
  call bcast_world_d(consistency)
  call bcast_world_d(findex)
  call bcast_world_l(lmassavariable)
  call bcast_world_d(lencorrmassa)
  call bcast_world_l(llencorrmassa)
  call bcast_world_l(lprintstatdat)
  call bcast_world_l(ltestread)
  call bcast_world_l(lrefinement)
  call bcast_world_l(ltagbeads)
  call bcast_world_l(lrefinementthreshold)
  call bcast_world_d(refinementthreshold)
  call bcast_world_l(lrefinementevery)
  call bcast_world_d(refinementevery)
  call bcast_world_l(lrefinementstart)
  call bcast_world_d(refinementstart)
  call bcast_world_l(lprintdatrem)
  call bcast_world_l(lrefbeadstart)
  call bcast_world_d(refbeadstartfit)
  call bcast_world_l(lreadrest)
  call bcast_world_l(lrestartreset)
  call bcast_world_l(luppot)
  call bcast_world_d(kuppot)
  call bcast_world_l(lmirror)
  call bcast_world_d(lenprobmassa)
  call bcast_world_d(massratio)
  call bcast_world_l(ldragvel)
  call bcast_world_i(typedragvel)
  call bcast_world_d(lenthresholdbead)
  call bcast_world_l(ldragvelfound)
  call bcast_world_l(lbreakup)
  call bcast_world_l(lfieldtype)
  call bcast_world_l(lfieldfreq)
  call bcast_world_i(nfieldtype)
  call bcast_world_d(fieldfreq)
  call bcast_world_l(lultimatestrength)
  call bcast_world_d(ultimatestrength)
  call bcast_world_l(lmultiplestep)
  call bcast_world_i(nmulstep)
  call bcast_world_l(lnmulstep)
  call bcast_world_l(lmultisteperror)
  call bcast_world_d(maxdispl)
  call bcast_world_l(lmaxdispl)
  call bcast_world_d(taoelectr)
  call bcast_world_l(ltaoelectr)
  call bcast_world_l(llenthresholdbead)
  call bcast_world_l(lflorentz)
  call bcast_world_l(lmagneticfield)
  call bcast_world_darr(BLor,3)
  
  if(lprintlist)then
    if(idrank/=0)allocate(printlist(nprintlist))
    call bcast_world_iarr(printlist,nprintlist)
    call label_argument(nprintlist,printlist,printarg)
  endif
  
  if(lprintlist2)then
    if(idrank/=0)allocate(printlist2(nprintlist2))
    call bcast_world_iarr(printlist2,nprintlist2)
    call label_argument(nprintlist2,printlist2,printarg2)
  endif
  
  ltest=ltestread
  
! check if an error is present in input file
  if(.not.lsystype)then
    call warning(4)
    ltest=.true.
  endif
  if(.not.lintegrator)then
    call warning(5)
    ltest=.true.
  endif
  if(.not.lresolution.and.(.not.lnpjet))then
    call warning(7)
    ltest=.true.
  endif
  if(.not.lilengt)then
    call warning(8)
    ltest=.true.
  endif
  if(.not.limass)then
    call warning(9)
    ltest=.true.
  endif
  if(.not.licharge)then
    call warning(10)
    ltest=.true.
  endif
  if(.not.licrossec)then
    call warning(11)
    ltest=.true.
  endif
  if(.not.lmu)then
    call warning(14)
    ltest=.true.
  endif
  if(.not.lG)then
    call warning(15)
    ltest=.true.
  endif
  if(.not.lh)then
    call warning(16)
    ltest=.true.
  endif
  if(.not.ltstep)then
    call warning(19)
    ltest=.true.
  endif
  if(.not.lendtime)then
    call warning(20)
    ltest=.true.
  endif
  if(lnpjet.and.lresolution)then
    call warning(22)
    ltest=.true.
  endif
  if(liniperturb)then
    if(.not.lpfreq)then
      call warning(31)
      ltest=.true.
    endif
    if(.not.lpampl)then
      call warning(32)
      ltest=.true.
    endif
  endif
  if(lairdrag)then
    if(.not.lairdragamp)then
      call warning(33)
      ltest=.true.
    endif
    if(integrator<4)then
      call warning(34,dble(integrator))
      ltest=.true.
    endif
    if(systype<4)then
      call warning(35,dble(systype))
      ltest=.true.
    endif
    if(.not.laird)then
      call warning(38)
      ltest=.true.
    endif
    if(.not.lairv)then
      call warning(39)
      ltest=.true.
    endif
  endif
  if(lHBfluid)then
    if(.not.lconsistency)then
      call warning(40)
      ltest=.true.
    endif
    if(.not.lfindex)then
      call warning(41)
      ltest=.true.
    endif
  else
    consistency=1.d0
    findex=1.d0
  endif
  if(ilength<icrossec)then
    call warning(42,ilength)
    call warning(43,icrossec)
    call warning(44)
    ltest=.true.
  endif
  
  if(pampl>ilength)then
    call warning(45,pampl)
    call warning(46,ilength)
    call warning(47)
    ltest=.true.
  endif
  
  if(lprintlisterror)then
    call warning(48)
    call warning(60)
    ltest=.true.
  endif
  
  if(lprintlisterror2)then
    call warning(58)
    call warning(60)
    ltest=.true.
  endif
  
  if(lmassavariable .and. ldevelopers)then
    if(.not.llencorrmassa)then
      call warning(51)
      ltest=.true.
    else
      if(lencorrmassa<0.d0)then
        call warning(55,lencorrmassa)
        ltest=.true.
      endif
    endif
  endif
  
  if(lsprintdat)then
    if(sprintdat>7 .or. sprintdat<1)then
      call warning(63,dble(sprintdat))
      ltest=.true.
    endif
    if(sprintdat==4)then
      ltrackbeads=.true.
      if(lrefinement)then
        call warning(67,dble(sprintdat))
        ltest=.true.
      endif
    endif
  endif
  
  if(lrefinement)then
    if(.not. lrefinementevery)then
      call warning(82)
      ltest=.true.
    endif
    if(.not. lrefinementthreshold)then
      call warning(83)
      ltest=.true.
    endif
  endif
  
  if(lKVfluid .and. lHBfluid)then
    call warning(87)
    ltest=.true.
  endif
  
  if(lKVfluid)then
    if(systype/=1 .and. systype/=3)then
      call warning(89,dble(systype))
      ltest=.true.
    endif
  endif
  
  if(ldevelopers)then
    if(lmassavariable .and. lrefinement)then
      if(llenthresholdbead)then
        lenthresholdbead=lencorrmassa
        call warning(84,lenthresholdbead)
      endif
    endif
    if(lmassavariable .and. (.not. ltagbeads))then
      call warning(70)
      ltest=.true.
    endif
  endif
  
  if(lfieldtype .and. nfieldtype/=0)then
    if(.not. lfieldfreq)then
      call warning(71)
      ltest=.true.
    endif
  endif
  
  if(lfieldtype)then
    if(nfieldtype<0 .or. nfieldtype>2)then
      call warning(73,dble(nfieldtype))
      ltest=.true.
    endif
    if(nfieldtype==2)then
      if(.not. ltaoelectr)then
        call warning(78)
        ltest=.true.
      endif
    endif
  endif
  
  if(lbreakup)then
    if(.not. lultimatestrength)then
      call warning(72)
      ltest=.true.
    endif
  endif
  
  if(lmultiplestep)then
    if(.not. ldcutoff)then
      call warning(74)
      ltest=.true.
    else
      if(nmulstep<=2)then
        call warning(75,dble(nmulstep))
        ltest=.true.
      endif
    endif
    if(.not. lnmulstep)then
      call warning(76)
      ltest=.true.
    endif
  endif
  
  if(sprintdat==6)then
    if(.not.(ltagbeads))then
      call warning(77,dble(sprintdat))
      ltest=.true.
    endif
  endif
  
  if(lpdbtagbeads)then
    if(.not. ltagbeads)then
      call warning(79)
      ltest=.true.
    endif
  endif
  
  if(lflorentz)then
    if(.not.lmagneticfield)then
      call warning(80)
      ltest=.true.
    endif
  endif
  
  if(lG)then
    if(G<1.d-4)then
      call warning(88,G) 
      ltest=.true.
    endif
  endif
  
  if(.not.lunits)units=1
  
  
  if(ltest)call error(7)
  
! convert the unit of measurement if they are not provided in cgs
! (only for developers)
  call conv_print_unit_jet()
  
  if(.not.lprintlist)then
    call warning(6)
  endif
  
! check air drag variables (only for developers)  
  if(lairdrag .and. (.not.lairvel))then
    call warning(49,velext)
  endif
  
  if((.not.lprinttime).and.ltstep)then
    lprinttime=.true.
    printtime=tstep*1000.d0
    call warning(59,printtime)
  endif
  
  if(lprintxyz)iprintxyz=nint(printxyz/tstep)
  if(lprintxyzsing)iprintxyzsing=nint(printxyzsing/tstep)
  if(lprintpdbsing)iprintpdbsing=nint(printpdbsing/tstep)
  if(lprintdat)iprintdat=nint(printdat/tstep)
  if(lrefinement)then
    if(lrefinementevery)then
      irefinementevery=nint(refinementevery/tstep)
    else
      refinementevery=1.d4*tstep
      irefinementevery=nint(refinementevery/tstep)
      call warning(68,refinementevery)
    endif
    if(lrefinementstart)then
      irefinementstart=nint(refinementstart/tstep)
    endif
  endif
  iprinttime=nint(printtime/tstep)
  eprinttime=nint(endtime/tstep)
  reprinttime=nint(dble(eprinttime)/dble(iprinttime))
  
  return
  
 end subroutine read_input
 
 subroutine identify_argument(iarg,arg,printcodsub,lenstring,lfound)
 
!***********************************************************************
!     
!     JETSPIN subroutine for identifying the symbolic string of the 
!     output observables
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer ,intent(in) :: iarg
  character(len=lenstring) ,allocatable, intent(in) :: arg(:)
  integer, intent(in) :: lenstring
  logical, intent(out) :: lfound
  integer, allocatable, intent(inout) :: printcodsub(:)
  
  integer :: inumchar
  character(len=lenstring) :: temps
  
  lfound=.false.
  temps=arg(iarg)
  
! for each symbolic string we associate an integer defined in 
!  printcodsub

  if(findstring('angl',temps,inumchar,lenstring))then
    printcodsub(iarg)=40
    lfound=.true.
  elseif(findstring('curc',temps,inumchar,lenstring))then
    printcodsub(iarg)=14
    lfound=.true.
  elseif(findstring('cpur',temps,inumchar,lenstring))then
    printcodsub(iarg)=23
    lfound=.true.
  elseif(findstring('cpue',temps,inumchar,lenstring))then
    printcodsub(iarg)=33
    lfound=.true.
  elseif(findstring('curn',temps,inumchar,lenstring))then
    printcodsub(iarg)=13
    lfound=.true.
  elseif(findstring('erms',temps,inumchar,lenstring))then
    printcodsub(iarg)=44
    lfound=.true.
    lmultisteperror=.true.
  elseif(findstring('mass',temps,inumchar,lenstring))then
    printcodsub(iarg)=34
    lfound=.true.
  elseif(findstring('mxsx',temps,inumchar,lenstring))then
    printcodsub(iarg)=42
    lfound=.true.
  elseif(findstring('mxst',temps,inumchar,lenstring))then
    printcodsub(iarg)=41
    lfound=.true.
  elseif(findstring('nref',temps,inumchar,lenstring))then
    printcodsub(iarg)=39
    lfound=.true.
  elseif(findstring('cpu',temps,inumchar,lenstring))then
    printcodsub(iarg)=22
    lfound=.true.
  elseif(findstring('mfc',temps,inumchar,lenstring))then
    printcodsub(iarg)=16
    lfound=.true.
  elseif(findstring('mfn',temps,inumchar,lenstring))then
    printcodsub(iarg)=15
    lfound=.true.
  elseif(findstring('nms',temps,inumchar,lenstring))then
    printcodsub(iarg)=43
    lfound=.true.
  elseif(findstring('rlp',temps,inumchar,lenstring))then
    printcodsub(iarg)=38
    lfound=.true.
  elseif(findstring('rrr',temps,inumchar,lenstring))then
    printcodsub(iarg)=21
    lfound=.true.
  elseif(findstring('sts',temps,inumchar,lenstring))then
    printcodsub(iarg)=28
    lfound=.true.
  elseif(findstring('svc',temps,inumchar,lenstring))then
    printcodsub(iarg)=18
    lfound=.true.
  elseif(findstring('vxs',temps,inumchar,lenstring))then
    printcodsub(iarg)=29
    lfound=.true.
  elseif(findstring('vys',temps,inumchar,lenstring))then
    printcodsub(iarg)=30
    lfound=.true.
  elseif(findstring('vzs',temps,inumchar,lenstring))then
    printcodsub(iarg)=31
    lfound=.true.
  elseif(findstring('yzs',temps,inumchar,lenstring))then
    printcodsub(iarg)=32
    lfound=.true.
  elseif(findstring('lp',temps,inumchar,lenstring))then
    printcodsub(iarg)=37
    lfound=.true.
  elseif(findstring('rc',temps,inumchar,lenstring))then
    printcodsub(iarg)=20
    lfound=.true.
  elseif(findstring('st',temps,inumchar,lenstring))then
    printcodsub(iarg)=5
    lfound=.true.
  elseif(findstring('ts',temps,inumchar,lenstring))then
    printcodsub(iarg)=24
    lfound=.true.
  elseif(findstring('vc',temps,inumchar,lenstring))then
    printcodsub(iarg)=17
    lfound=.true.
  elseif(findstring('vn',temps,inumchar,lenstring))then
    printcodsub(iarg)=36
    lfound=.true.
  elseif(findstring('vx',temps,inumchar,lenstring))then
    printcodsub(iarg)=6
    lfound=.true.
  elseif(findstring('vy',temps,inumchar,lenstring))then
    printcodsub(iarg)=7
    lfound=.true.
  elseif(findstring('vz',temps,inumchar,lenstring))then
    printcodsub(iarg)=8
    lfound=.true.
  elseif(findstring('xs',temps,inumchar,lenstring))then
    printcodsub(iarg)=25
    lfound=.true.
  elseif(findstring('ys',temps,inumchar,lenstring))then
    printcodsub(iarg)=26
    lfound=.true.
  elseif(findstring('zs',temps,inumchar,lenstring))then
    printcodsub(iarg)=27
    lfound=.true.
  elseif(findstring('yz',temps,inumchar,lenstring))then
    printcodsub(iarg)=9
    lfound=.true.
  elseif(findstring('f',temps,inumchar,lenstring))then
    printcodsub(iarg)=11
    lfound=.true.
  elseif(findstring('l',temps,inumchar,lenstring))then
    printcodsub(iarg)=12
    lfound=.true.
  elseif(findstring('n',temps,inumchar,lenstring))then
    printcodsub(iarg)=10
    lfound=.true.
  elseif(findstring('q',temps,inumchar,lenstring))then
    printcodsub(iarg)=35
    lfound=.true.
  elseif(findstring('t',temps,inumchar,lenstring))then
    printcodsub(iarg)=1
    lfound=.true.
  elseif(findstring('v',temps,inumchar,lenstring))then
    printcodsub(iarg)=45
    lfound=.true.
  elseif(findstring('x',temps,inumchar,lenstring))then
    printcodsub(iarg)=2
    lfound=.true.
  elseif(findstring('y',temps,inumchar,lenstring))then
    printcodsub(iarg)=3
    lfound=.true.
  elseif(findstring('z',temps,inumchar,lenstring))then
    printcodsub(iarg)=4
    lfound=.true.
  else
    lfound=.false.
  endif
  
  
  return
  
 end subroutine identify_argument
 
 subroutine label_argument(narg,printcodsub,printlisub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for associating a description string 
!     for each output observables requested in input file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer ,intent(in) :: narg
  integer, allocatable, intent(in) :: printcodsub(:)
  character(len=20), allocatable, intent(out) :: printlisub(:)
  
  integer :: iarg
  
  if(allocated(printlisub))deallocate(printlisub)
  allocate(printlisub(narg))
  
  do iarg=1,narg
  printlisub(iarg)=repeat(' ',20)
  if(printcodsub(iarg)==1)then
    printlisub(iarg)='t (s)'
  elseif(printcodsub(iarg)==24)then
    printlisub(iarg)='ts'
  elseif(printcodsub(iarg)==2)then
    printlisub(iarg)='x (cm)'
  elseif(printcodsub(iarg)==3)then
    printlisub(iarg)='y (cm)'
  elseif(printcodsub(iarg)==4)then
    printlisub(iarg)='z (cm)'
  elseif(printcodsub(iarg)==25)then
    printlisub(iarg)='xs'
  elseif(printcodsub(iarg)==26)then
    printlisub(iarg)='ys'
  elseif(printcodsub(iarg)==27)then
    printlisub(iarg)='zs'
  elseif(printcodsub(iarg)==5)then
    printlisub(iarg)='st (g cm^-1 s^-2)'
  elseif(printcodsub(iarg)==28)then
    printlisub(iarg)='sts'
  elseif(printcodsub(iarg)==6)then
    printlisub(iarg)='vx (cm s^-1)'
  elseif(printcodsub(iarg)==7)then
    printlisub(iarg)='vy (cm s^-1)'
  elseif(printcodsub(iarg)==8)then
    printlisub(iarg)='vz (cm s^-1)'
  elseif(printcodsub(iarg)==29)then
    printlisub(iarg)='vxs'
  elseif(printcodsub(iarg)==30)then
    printlisub(iarg)='vys'
  elseif(printcodsub(iarg)==31)then
    printlisub(iarg)='vzs'
  elseif(printcodsub(iarg)==9)then
    printlisub(iarg)='yz (cm)'
  elseif(printcodsub(iarg)==32)then
    printlisub(iarg)='yzs'
  elseif(printcodsub(iarg)==34)then
    printlisub(iarg)='mass (g)'
  elseif(printcodsub(iarg)==35)then
    printlisub(iarg)='q (statC)'
  elseif(printcodsub(iarg)==22)then
    printlisub(iarg)='cpu (s)'
  elseif(printcodsub(iarg)==23)then
    printlisub(iarg)='cpur (s)'
  elseif(printcodsub(iarg)==33)then
    printlisub(iarg)='cpue (s)'
  elseif(printcodsub(iarg)==10)then
    printlisub(iarg)='n'
  elseif(printcodsub(iarg)==11)then
    printlisub(iarg)='f'
  elseif(printcodsub(iarg)==12)then
    printlisub(iarg)='l'
  elseif(printcodsub(iarg)==13)then
    printlisub(iarg)='curn (statC s^-1)'
  elseif(printcodsub(iarg)==14)then
    printlisub(iarg)='curc (statC s^-1)'
  elseif(printcodsub(iarg)==36)then
    printlisub(iarg)='vn  (cm s^-1)'
  elseif(printcodsub(iarg)==37)then
    printlisub(iarg)='lp  (cm)'
  elseif(printcodsub(iarg)==38)then
    printlisub(iarg)='rlp'
  elseif(printcodsub(iarg)==39)then
    printlisub(iarg)='nref'
  elseif(printcodsub(iarg)==40)then
    printlisub(iarg)='angl (°)'
  elseif(printcodsub(iarg)==17)then
    printlisub(iarg)='vc  (cm s^-1)'
  elseif(printcodsub(iarg)==18)then
    printlisub(iarg)='svc (s^-1)'
  elseif(printcodsub(iarg)==15)then
    printlisub(iarg)='mfn (g s^-1)'
  elseif(printcodsub(iarg)==16)then
    printlisub(iarg)='mfc (g s^-1)'
  elseif(printcodsub(iarg)==20)then
    printlisub(iarg)='rc (cm)'
  elseif(printcodsub(iarg)==21)then
    printlisub(iarg)='rrr'
  elseif(printcodsub(iarg)==41)then
    printlisub(iarg)='mxst (g cm^-1 s^-2)'
  elseif(printcodsub(iarg)==42)then
    printlisub(iarg)='mxsx (cm)'
  elseif(printcodsub(iarg)==43)then
    printlisub(iarg)='nms'
  elseif(printcodsub(iarg)==44)then
    printlisub(iarg)='erms (dyne)'
  elseif(printcodsub(iarg)==45)then
    printlisub(iarg)='v (statV)'
  endif
  printlisub(iarg)=adjustr(printlisub(iarg))
  enddo
  
  return
  
 end subroutine label_argument
 
 subroutine conv_print_unit_jet()
 
!***********************************************************************
!     
!     JETSPIN subroutine for converting and printing on terminal
!     the input parameters
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
  
  implicit none
  
  character(len=10) :: r_char
  character(len=11),parameter :: form1='(a30,a3,i8)'
  character(len=4),parameter :: form2='(3a)'
  character(len=3),parameter :: form3='(a)'
  character(len=17),parameter :: form4='(a30,a3,g20.10,a)'
  character(len=15),parameter :: form5='(a30,a3,g20.10)'
  character(len=12),parameter :: form6='(a30,a3,i10)'
  character(len=14),parameter :: form7='(a30,a3,i10,a)'
  
  character(len=30) :: labelsub
  character(len=3),parameter :: ugualab=' = '
  
  if(idrank==0)then
    write(6,'(/,a,/)')'Start printing input file'
      labelsub='system'
      write(6,form1)labelsub,ugualab,systype
      labelsub='integrator'
      write(6,form1)labelsub,ugualab,integrator
    if(lmyseed)then
      labelsub='seed'
      write(6,form1)labelsub,ugualab,myseed
    endif
  endif
  
  select case(units)
  case (2)
    if(idrank==0)then
      write(6,form3)'units real'
      if(lresolution)then
      labelsub='resolution'
      write(6,form4)labelsub,ugualab,resolution,' cm'
      endif
      if(lnpjet)then
      labelsub='points of jet discretization'
      write(6,form6)labelsub,ugualab,npjet
      endif
      labelsub='timestep of integration'
      write(6,form4)labelsub,ugualab,tstep,' s'
      labelsub='final time of integration'
      write(6,form4)labelsub,ugualab,endtime,' s'
      if(lprinttime)then
      labelsub='print time every'
      write(6,form4)labelsub,ugualab,printtime,' s'
      endif
      if(lprintxyz)then
      labelsub='print xyz file every'
      write(6,form4)labelsub,ugualab,printxyz,' s'
      endif
      if(lprintxyzsing)then
      labelsub='print a single xyz file for each frame every'
      write(6,form4)labelsub,ugualab,printxyzsing,' s'
      endif
      if(lprintpdbsing)then
      labelsub='print a single pdb file for each frame every'
      write(6,form4)labelsub,ugualab,printpdbsing,' s'
      endif
      if(lxyzrescale)then
      labelsub='print xyz file rescaled by'
      write(6,form4)labelsub,ugualab,xyzrescale,' factor'
      endif
      if(lpdbrescale)then
      labelsub='print pdb file rescaled by'
      write(6,form4)labelsub,ugualab,pdbrescale,' factor'
      endif
      if(lpdbtagbeads)then
      write(6,form3)'print pdb with tagged beads yes'
      endif
      if(lprintdat)then
      labelsub='print binary file every'
      write(6,form4)labelsub,ugualab,printdat,' s'
      endif
      if(lprintstatdat)then
      labelsub='printstat binary file every'
      write(6,form4)labelsub,ugualab,printtime,' s'
      endif
      if(lilengt)then
      labelsub='initial length of jet'
      write(6,form4)labelsub,ugualab,ilength,' cm'
      endif
      if(listress)then
      labelsub='nozzle stress of jet'
      write(6,form4)labelsub,ugualab,istress,' g cm^-1 s^-2'
      endif
      if(livelocity)then
      labelsub='nozzle velocity of jet'
      write(6,form4)labelsub,ugualab,ivelocity,' cm s^-1'
      endif
      if(lHBfluid)then
      write(6,form3)"hbfluid yes - Herschel-Bulkley mode"
      endif
      if(lKVfluid)then
      write(6,form3)"kvfluid yes - Kelvin–Voigt mode"
      endif
      if(lconsistency)then
      labelsub='consistency of Herschel-Bulkley jet = '
      write(6,form4)labelsub,ugualab,consistency, &
      ' times the viscosity of jet'
      endif
      if(lfindex)then
      labelsub='index of Herschel-Bulkley jet'
      write(6,form5)labelsub,ugualab,findex
      endif
      if(lyieldstress)then
      labelsub='yield stress of jet'
      write(6,form4)labelsub,ugualab,yieldstress,' g cm^-1 s^-2'
      endif
      labelsub='mass of jet'
      write(6,form4)labelsub,ugualab,imassa,' g cm^-3'
      labelsub='charge of jet'
      write(6,form4)labelsub,ugualab,icharge,' C L^-1'
      labelsub='nozzle cross section of jet'
      write(6,form4)labelsub,ugualab,icrossec,' cm'
      labelsub='viscosity of jet'
      write(6,form4)labelsub,ugualab,mu,' g cm^-1 s^-1'
      labelsub='elastic modulus of jet'
      write(6,form4)labelsub,ugualab,G,' g cm^-1 s^-2'
      labelsub='collector distance'
      write(6,form4)labelsub,ugualab,h,' cm'
      labelsub='external potential'
      write(6,form4)labelsub,ugualab,V0,' V'
      if(lfieldtype)then
      labelsub='external potential type'
      write(6,form6)labelsub,ugualab,nfieldtype
      endif
      if(lfieldfreq)then
      labelsub='external potential frequency'
      write(6,form4)labelsub,ugualab,fieldfreq,' s^-1'
      endif
      if(ltaoelectr)then
      labelsub='external potential time'
      write(6,form4)labelsub,ugualab,taoelectr,' s'
      endif
      labelsub='surface tension of the jet'
      write(6,form4)labelsub,ugualab,surfacet,' g s^-2'
      if(liniperturb)then
      labelsub='nozzle perturbation frequency'
      write(6,form4)labelsub,ugualab,pfreq,' s^-1'
      labelsub='nozzle perturbation amplitude'
      write(6,form4)labelsub,ugualab,pampl,' cm'
      endif
      if(lmultiplestep)then
      write(6,form3)"multiple step yes"
      labelsub='compute all the forces every'
      write(6,form6)labelsub,ugualab,nmulstep
      endif
      if(lmultiplestep .and. ldcutoff)then
      labelsub='cutoff of the primary shell'
      write(6,form4)labelsub,ugualab,dcutoff,' cm'
      endif
      if(lmultiplestep .and. lmaxdispl)then
      labelsub='maximum displacement'
      write(6,form4)labelsub,ugualab,maxdispl,' cm'
      endif
      if(lairdrag)then
      write(6,form3)"airdrag yes - stochastic mode"
      labelsub='airdrag amplitude'
      write(6,form5)labelsub,ugualab,airdragamp(1)
      labelsub='airdrag airdensity'
      write(6,form4)labelsub,ugualab,aird,' g cm^-3'
      labelsub='airdrag airviscosity'
      write(6,form4)labelsub,ugualab,airv,' cm^2 s^-1'
      if(lairvel)then
      labelsub='airdrag airvelocity'
      write(6,form4)labelsub,ugualab,velext,' cm s^-1'
      endif
      endif
      if(lmassavariable .and. ldevelopers)then
      write(6,form3)"variable mass yes"
      labelsub='variable mass distance between NPs'
      write(6,form4)labelsub,ugualab,lencorrmassa,' cm'
      labelsub='npradius'
      write(6,form4)labelsub,ugualab,lenprobmassa,' ONLY DEVELOPERS'
      labelsub='npmass'
      write(6,form4)labelsub,ugualab,massratio,' ONLY DEVELOPERS'
      endif
      if(lrefinement)then
        write(6,form3) &
         "dynamic refinement yes"
      if(lrefinementstart)then
      labelsub='dynamic refinement start'
      write(6,form4)labelsub,ugualab,refinementstart,' s'
      endif
      if(lrefinementthreshold)then
      labelsub='dynamic refinement threshold'
      write(6,form4)labelsub,ugualab,refinementthreshold,' cm'
      endif
      if(llenthresholdbead)then
      labelsub='dynamic refinement anchor bead distance'
      write(6,form4)labelsub,ugualab,lenthresholdbead, &
       ' cm'
      endif
      if(lrefinementevery)then
      labelsub='dynamic refinement every'
      write(6,form4)labelsub,ugualab,refinementevery,' s'
      endif
      endif
      if(lreadrest)then
      write(6,form3)"restart yes"
      if(lrestartreset)then
      write(6,form3)"reset statistical data"
      endif
      endif
      if(luppot)then
      write(6,form3)"wall yes"
      if(kuppot/=0.d0)then
      labelsub='wall constant'
      write(6,form4)labelsub,ugualab,kuppot,' dyne/cm'
      endif
      endif
      if(lmirror)then
      write(6,form3)"mirror yes"
      endif
      if(lbreakup)then
      write(6,form3)"breakup yes"
      endif
      if(lultimatestrength)then
      labelsub='ultimate strength of jet'
      write(6,form4)labelsub,ugualab,ultimatestrength,' g cm^-1 s^-2'
      endif
      if(ldragvelfound)then
        if(ldragvel)then
          write(6,form3)"velocity drag yes"
        else
          write(6,form3)"velocity drag no"
        endif
      else
        write(6,form3)"velocity drag yes by default"
      endif
      if(typedragvel/=0 .and. ldevelopers)then
        labelsub='type of velocity drag'
        write(6,form6)labelsub,ugualab,typedragvel
      endif
      if(lflorentz)then
        write(6,form3)"lorentz yes"
        if(lmagneticfield)then
          labelsub='magnetic field along x'
          write(6,form4)labelsub,ugualab,BLor(1),' cm^-0.5 g^0.5 s^-1'
          labelsub='magnetic field along y'
          write(6,form4)labelsub,ugualab,BLor(2),' cm^-0.5 g^0.5 s^-1'
          labelsub='magnetic field along z'
          write(6,form4)labelsub,ugualab,BLor(3),' cm^-0.5 g^0.5 s^-1'
        endif
      endif
    endif
    icharge=icharge/1000.d0*2997919999.93d0
    V0=V0/299.792458d0 !statV=cm^0.5 g^0.5 s^-1
    consistency=consistency*mu
  case default
    if(idrank==0)then
      write(6,form3)'units cgs'
      if(lresolution)then
      labelsub='resolution'
      write(6,form4)labelsub,ugualab,resolution,' cm'
      endif
      if(lnpjet)then
      labelsub='points of jet discretization'
      write(6,form6)labelsub,ugualab,npjet
      endif
      labelsub='timestep of integration'
      write(6,form4)labelsub,ugualab,tstep,' s'
      labelsub='final time of integration'
      write(6,form4)labelsub,ugualab,endtime,' s'
      if(lprinttime)then
      labelsub='print time every'
      write(6,form4)labelsub,ugualab,printtime,' s'
      endif
      if(lprintxyz)then
      labelsub='print xyz file every'
      write(6,form4)labelsub,ugualab,printxyz,' s'
      endif
      if(lprintxyzsing)then
      labelsub='print a single xyz file for each frame every'
      write(6,form4)labelsub,ugualab,printxyzsing,' s'
      endif
      if(lprintpdbsing)then
      labelsub='print a single pdb file for each frame every'
      write(6,form4)labelsub,ugualab,printpdbsing,' s'
      endif
      if(lxyzrescale)then
      labelsub='print xyz file rescaled by'
      write(6,form4)labelsub,ugualab,xyzrescale,' factor'
      endif
      if(lpdbrescale)then
      labelsub='print pdb file rescaled by'
      write(6,form4)labelsub,ugualab,pdbrescale,' factor'
      endif
      if(lpdbtagbeads)then
      write(6,form3)'print pdb with tagged beads yes'
      endif
      if(lprintdat)then
      labelsub='print binary file every'
      write(6,form4)labelsub,ugualab,printdat,' s'
      endif
      if(lprintstatdat)then
      labelsub='printstat binary file every'
      write(6,form4)labelsub,ugualab,printtime,' s'
      endif
      if(lilengt)then
      labelsub='initial length of jet'
      write(6,form4)labelsub,ugualab,ilength,' cm'
      endif
      if(listress)then
      labelsub='nozzle stress of jet'
      write(6,form4)labelsub,ugualab,istress,' g cm^-1 s^-2'
      endif
      if(livelocity)then
      labelsub='nozzle velocity of jet'
      write(6,form4)labelsub,ugualab,ivelocity,' cm s^-1'
      endif
      if(lHBfluid)then
      write(6,form3)"hbfluid yes - Herschel-Bulkley mode"
      endif
      if(lKVfluid)then
      write(6,form3)"kvfluid yes - Kelvin–Voigt mode"
      endif
      if(lconsistency)then
      labelsub='consistency of Herschel-Bulkley jet = '
      write(6,form4)labelsub,ugualab,consistency, &
      ' times the viscosity of jet'
      endif
      if(lfindex)then
      labelsub='index of Herschel-Bulkley jet'
      write(6,form5)labelsub,ugualab,findex
      endif
      if(lyieldstress)then
      labelsub='yield stress of jet'
      write(6,form4)labelsub,ugualab,yieldstress,' g cm^-1 s^-2'
      endif
      labelsub='density mass of jet'
      write(6,form4)labelsub,ugualab,imassa,' g cm^-3'
      labelsub='density charge of jet'
      write(6,form4)labelsub,ugualab,icharge,' statC cm^-3'
      labelsub='nozzle cross section of jet'
      write(6,form4)labelsub,ugualab,icrossec,' cm'
      labelsub='viscosity of jet'
      write(6,form4)labelsub,ugualab,mu,' g cm^-1 s^-1'
      labelsub='elastic modulus of jet'
      write(6,form4)labelsub,ugualab,G,' g cm^-1 s^-2'
      labelsub='collector distance'
      write(6,form4)labelsub,ugualab,h,' cm'
      labelsub='external potential'
      write(6,form4)labelsub,ugualab,V0,' statV'
      if(lfieldtype)then
      labelsub='external potential type'
      write(6,form6)labelsub,ugualab,nfieldtype
      endif
      if(lfieldfreq)then
      labelsub='external potential frequency'
      write(6,form4)labelsub,ugualab,fieldfreq,' s^-1'
      endif
      if(ltaoelectr)then
      labelsub='external potential time'
      write(6,form4)labelsub,ugualab,taoelectr,' s'
      endif
      labelsub='surface tension of the jet'
      write(6,form4)labelsub,ugualab,surfacet,' g s^-2'
      if(liniperturb)then
      labelsub='nozzle perturbation frequency'
      write(6,form4)labelsub,ugualab,pfreq,' s^-1'
      labelsub='nozzle perturbation amplitude'
      write(6,form4)labelsub,ugualab,pampl,' cm'
      endif
      if(lmultiplestep)then
      write(6,form3)"multiple step yes"
      labelsub='compute all the forces every'
      write(6,form6)labelsub,ugualab,nmulstep
      endif
      if(lmultiplestep .and. ldcutoff)then
      labelsub='cutoff of the primary shell'
      write(6,form4)labelsub,ugualab,dcutoff,' cm'
      endif
      if(lmultiplestep .and. lmaxdispl)then
      labelsub='maximum displacement'
      write(6,form4)labelsub,ugualab,maxdispl,' cm'
      endif
      if(lairdrag)then
      write(6,form3)"airdrag yes - stochastic mode"
      labelsub='airdrag amplitude'
      write(6,form5)labelsub,ugualab,airdragamp(1)
      labelsub='airdrag airdensity'
      write(6,form4)labelsub,ugualab,aird,' g cm^-3'
      labelsub='airdrag airviscosity'
      write(6,form4)labelsub,ugualab,airv,' cm^2 s^-1'
      if(lairvel)then
      labelsub='airdrag airvelocity'
      write(6,form4)labelsub,ugualab,velext,' cm s^-1'
      endif
      endif
      if(lmassavariable .and. ldevelopers)then
      write(6,form3)"variable mass yes"
      labelsub='variable mass distance between NPs'
      write(6,form4)labelsub,ugualab,lencorrmassa,' cm'
      labelsub='npradius'
      write(6,form4)labelsub,ugualab,lenprobmassa,' ONLY DEVELOPERS'
      labelsub='npmass'
      write(6,form4)labelsub,ugualab,massratio,' ONLY DEVELOPERS'
      endif
      if(lrefinement)then
        write(6,form3) &
         "dynamic refinement yes"
      if(lrefinementstart)then
      labelsub='dynamic refinement start'
      write(6,form4)labelsub,ugualab,refinementstart,' s'
      endif
      if(lrefinementthreshold)then
      labelsub='dynamic refinement threshold'
      write(6,form4)labelsub,ugualab,refinementthreshold,' cm'
      endif
      if(llenthresholdbead)then
      labelsub='dynamic refinement anchor bead distance'
      write(6,form4)labelsub,ugualab,lenthresholdbead, &
       ' cm'
      endif
      if(lrefinementevery)then
      labelsub='dynamic refinement every'
      write(6,form4)labelsub,ugualab,refinementevery,' s'
      endif
      endif
      if(lreadrest)then
      write(6,form3)"restart yes"
      if(lrestartreset)then
      write(6,form3)"reset statistical data"
      endif
      endif
      if(luppot)then
      write(6,form3)"wall yes"
      if(kuppot/=0.d0)then
      labelsub='wall constant'
      write(6,form4)labelsub,ugualab,kuppot,' dyne/cm'
      endif
      endif
      if(lmirror)then
      write(6,form3)"mirror yes"
      endif
      if(lbreakup)then
      write(6,form3)"breakup yes"
      endif
      if(lultimatestrength)then
      labelsub='ultimate strength of jet'
      write(6,form4)labelsub,ugualab,ultimatestrength,' g cm^-1 s^-2'
      endif
      if(ldragvelfound)then
        if(ldragvel)then
          write(6,form3)"velocity drag yes"
        else
          write(6,form3)"velocity drag no"
        endif
      else
        write(6,form3)"velocity drag yes by default"
      endif
      if(typedragvel/=0 .and. ldevelopers)then
        labelsub='type of velocity drag'
        write(6,form6)labelsub,ugualab,typedragvel
      endif
      if(lflorentz)then
        write(6,form3)"lorentz yes"
        if(lmagneticfield)then
          labelsub='magnetic field along x'
          write(6,form4)labelsub,ugualab,BLor(1),' cm^-0.5 g^0.5 s^-1'
          labelsub='magnetic field along y'
          write(6,form4)labelsub,ugualab,BLor(2),' cm^-0.5 g^0.5 s^-1'
          labelsub='magnetic field along z'
          write(6,form4)labelsub,ugualab,BLor(3),' cm^-0.5 g^0.5 s^-1'
        endif
      endif
    endif
    consistency=consistency*mu
  end select
  
  if(idrank==0)then
    if(lprintlist)then
      if(linserting)then
      write(6,form3)"inserting yes - inserting mode"
      endif
      if(lremove)then
      write(6,form3)"removing yes - removing mode"
      endif
      if(lgravity)then
      write(6,form3)"gravity yes - gravity mode"
      endif
      if(nprintlist>0)then
      write (r_char,'(i10)')nprintlist
      write(6,form2)"print list specified with ", &
       trim(adjustl(r_char))," arguments"
      endif
      if(nprintlist2>0)then
      write (r_char,'(i10)')nprintlist2
      write(6,form2)"printstat list specified with ", &
       trim(adjustl(r_char))," arguments"
      endif
    endif
    write(6,'(/,a,/)')'Finish printing input file'
  endif
  
  return
  
 end subroutine conv_print_unit_jet
 
 subroutine outprint_driver(k,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for managing the output subroutines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: k
  double precision, intent(in) :: timesub
  
  integer :: i,j
  double precision :: tempint
  
  if(.not.lprinttime)return
  if(mod(k,iprinttime)/=0)return
  
  tempint=dble(iprinttime)*tstep
  call compute_statistic(tempint,timesub,k)
  
  if(idrank/=0)return
  
  do i=1,nprintlist
    j=printlist(i)
    xprint(i)=statdata(j)
  enddo

  call outprint_term(k,nprintlist,xprint,printarg)
  
  if(lprintlist2)then
    
    do i=1,nprintlist2
      j=printlist2(i)
      xprint2(i)=statdata(j)
    enddo
    
    call outprint_file(16,namefile,k,nprintlist2,xprint2,printarg2)
    
  endif
  
  call write_stat_dat(lprintstatdat,170,'statdat.dat',k)
  
  return
  
 end subroutine outprint_driver
 
 subroutine outprint_term(k,nprint1,dprint,printargsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing records on the terminal
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in) :: k,nprint1
  double precision,intent(in),dimension(1:nprint1) :: dprint
  character(len=20),intent(in),dimension(1:nprint1) :: printargsub
  
  integer :: i
  
  logical, save :: lfirst=.true.
  integer,save :: icount=0
  
  if(idrank/=0)return
  
  if(lfirst .or. (mod(icount,100)==0))then
    lfirst=.false.    
    write(6,'(a10)',advance="no")'#*********'
    do i=1,nprint1
      write(6,'(a20)',advance="no")' *******************'
    enddo
    write(6,'(1x)')           
    write(6,'(a10)',advance="no")'#    nstep'
    do i=1,nprint1
      write(6,'(a20)',advance="no")printargsub(i)
    enddo
    write(6,'(1x)')
    write(6,'(a10)',advance="no")'#*********'
    do i=1,nprint1
      write(6,'(a20)',advance="no")' *******************'
    enddo
    write(6,'(1x)')
  endif
  
  write(6,'(i10)',advance="no")k
  do i=1,nprint1
    write(6,'(g20.10)',advance="no")dprint(i)
  enddo
  
  write(6,'(1x)')
  call flush(6)
  
  icount=icount+1
  
  return
  
 end subroutine outprint_term
 
 subroutine finish_print(k,timesub,itime,ftime)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the last record on the terminal
!     and close the output file 'statdat.dat'
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: k
  double precision, intent(in) :: timesub,itime,ftime
  
  if(idrank/=0)return
  write(6,'(/,a,/)')'Program closed correctly'
  if(mxrank==1)then
    write(6,'(a,g20.10,a,i2,a,/)')'CPU time = ',(ftime-itime), &
     ' seconds on ',mxrank,' CPU'
  else
    write(6,'(a,g20.10,a,i2,a,/)')'CPU time = ',(ftime-itime), &
     ' seconds on ',mxrank,' CPUs'
  endif
  call flush(6)
  close(16)
  
  return
  
 end subroutine finish_print
 
 subroutine outprint_file(outp,soutp,k,nprint1,dprint,printargsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing records on the output file 
!     'statdat.dat'
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in) :: outp,k,nprint1
  character(len=*),intent(in) :: soutp
  double precision,intent(in),dimension(1:nprint1) :: dprint
  character(len=20),intent(in),dimension(1:nprint1) :: printargsub
  
  logical, save :: lfirst=.true.
  
  integer :: i
  
  if(idrank/=0)return
  
  if(lfirst)then
    lfirst=.false.
    if(lreadrest)then
      open(unit=outp,file=trim(soutp),form='formatted', &
       status='old',action='write',position='append')
    else
      open(unit=outp,file=trim(soutp),form='formatted', &
       status='replace',action='write')
    endif  
    write(outp,'(a10)',advance="no")'#*********'
    do i=1,nprint1
      write(outp,'(a20)',advance="no")' *******************'
    enddo
    write(outp,'(1x)')           
    write(outp,'(a10)',advance="no")'#    nstep'
    do i=1,nprint1
      write(outp,'(a20)',advance="no")printargsub(i)
    enddo
    write(outp,'(1x)')
    write(outp,'(a10)',advance="no")'#*********'
    do i=1,nprint1
      write(outp,'(a20)',advance="no")' *******************'
    enddo
    write(outp,'(1x)')
  endif
  
  
  write(outp,'(i10)',advance="no")k
  do i=1,nprint1
    write(outp,'(g20.10)',advance="no")dprint(i)
  enddo
  
  write(outp,'(1x)')
  
  call flush(outp)
  
  return
  
 end subroutine outprint_file
 
 subroutine open_dat_file(lprintdatsub,fileout,filename)
 
!***********************************************************************
!     
!     JETSPIN subroutine for opening the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lprintdatsub
  integer, intent(in) :: fileout
  character(len=*), intent(in) :: filename
  
  if(idrank/=0)return
  if(.not.lprintdatsub)return
  
  if(lreadrest)then
    open(fileout,file=trim(filename),form='unformatted', &
     status='old',action='write',position='append')
  else
    open(fileout,file=trim(filename),form='unformatted', &
     status='replace',action='write')
  endif
  
  return
  
 end subroutine open_dat_file
 
 subroutine open_datrem_file(lprintdatsub,fileout,filename)
 
!***********************************************************************
!     
!     JETSPIN subroutine for opening the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lprintdatsub
  integer, intent(in) :: fileout
  character(len=*), intent(in) :: filename
  
  if(idrank/=0)return
  
  if(.not.lprintdatrem)return
  
  if(lreadrest)then
    open(fileout,file=trim(filename),form='unformatted', &
     status='old',action='write',position='append')
  else
    open(fileout,file=trim(filename),form='unformatted', &
     status='replace',action='write')
  endif
  
  return
  
 end subroutine open_datrem_file
 
 subroutine write_restart_file(ievery,fileout,filename,k,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for writing the binary restart.dat file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ievery,fileout,k
  character(len=*), intent(in) :: filename
  double precision, intent(in) :: timesub
  
  integer :: mysprintdat
  
  if(idrank/=0)return
  
  if(mod(k,ievery)/=0)return
  
  open(fileout,file=filename,form='unformatted',status='replace', &
   action='write')
   
  call write_dat_parameter(.true.,fileout,timesub)
  call set_sprintdat(mysprintdat)
  call write_dat_frame(.true.,fileout,k,timesub,1, &
  inpjet,npjet,mysprintdat,systype,linserted)
  call close_dat_file(.true.,fileout)
  
  return
  
 end subroutine write_restart_file
 
 subroutine read_restart_file(fileout,filename, &
  k,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for restarting JETSPIN from a binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: fileout
  integer, intent(inout) :: k
  character(len=*), intent(in) :: filename
  double precision, intent(inout) :: timesub
  
  logical :: ltest,lexist
  
  if(.not.lreadrest)return
  
  inquire(file=filename,exist=lexist)
  if(.not.lexist)call error(15)
  
  if(idrank==0)then
    open(fileout,file=filename,form='unformatted',status='old', &
     action='read')
  endif
   
  call read_dat_parameter(fileout,timesub)
  call read_dat_restart(fileout,k,timesub)
  
  k=nint(timesub/tstep)
  
  if(idrank==0)then
    close(fileout)
  endif
  
  call warning(69,timesub*tao)
  
  return
  
 end subroutine read_restart_file
 
 subroutine read_dat_parameter(fileout,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for reading the parameters from 
!     the restart file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  

  integer, intent(in) :: fileout
  double precision, intent(inout) :: timesub
  
  integer :: natms,i,mioind
  integer :: sprintdatsub
  logical :: lstart
  
  double precision :: dtemp(15)
  integer :: itemp(15)
  logical :: ltemp(15)
  
  if(idrank==0)then
    read(fileout)lstart
    read(fileout)natms,sprintdatsub,systype,timesub,lstart
  endif
  
  call bcast_world_i(systype)
  call bcast_world_d(timesub)
  
  
  if(idrank==0)then
    read(fileout)itemp(1),systype,integrator,units
    read(fileout)resolution,ilength, &
     imassa,icharge,icrossec, &
     istress,ivelocity
  endif
   
  call bcast_world_i(systype)
  call bcast_world_i(integrator)
  call bcast_world_i(units)
  call bcast_world_d(resolution)
  call bcast_world_d(ilength)
  call bcast_world_d(imassa)
  call bcast_world_d(icharge)
  call bcast_world_d(icrossec)
  call bcast_world_d(istress)
  call bcast_world_d(ivelocity)
   
  if(idrank==0)then
    read(fileout)pfreq,pampl, &
     airdragamp(1), &
     airdragamp(2), &
     airdragamp(3),aird,airv
    read(fileout)mu,G,yieldstress,h,V0,surfacet, &
     tstep,consistency,findex
    read(fileout)massscale,chargescale,lengthscale,tao
  endif
  
  call bcast_world_d(pfreq)
  call bcast_world_d(pampl)
  call bcast_world_darr(airdragamp,3)
  call bcast_world_d(aird)
  call bcast_world_d(airv)
  call bcast_world_d(mu)
  call bcast_world_d(G)
  call bcast_world_d(yieldstress)
  call bcast_world_d(h)
  call bcast_world_d(V0)
  call bcast_world_d(surfacet)
  call bcast_world_d(tstep)
  call bcast_world_d(consistency)
  call bcast_world_d(findex)
  call bcast_world_d(massscale)
  call bcast_world_d(chargescale)
  call bcast_world_d(lengthscale)
  call bcast_world_d(tao)
  
  resolution=(resolution/lengthscale)
  ilength=(ilength/lengthscale)
  imassa=imassa/(massscale)
  icharge=(icharge/chargescale)
  icrossec=(icrossec/lengthscale)
  istress=(istress/G)
  ivelocity=(ivelocity*tao/lengthscale) 
  yieldstress=(yieldstress/G)
  h=(h/lengthscale)
  tstep=(tstep/tao)
  consistency=(consistency/(mu*((tao)**(findex-1.d0))))
  
  timesub=(timesub/tao)
  
  pfreq=(pfreq*tao)
  pampl=(pampl/lengthscale)
  airdragamp(1)=(airdragamp(1)/(lengthscale**2.d0)*(tao**3.d0))
  airdragamp(2)=(airdragamp(2)/(lengthscale**2.d0)*(tao**3.d0))
  airdragamp(3)=(airdragamp(3)/(lengthscale**2.d0)*(tao**3.d0))
  
  if(idrank==0)then
    read(fileout)dtemp(1), &
     lencorrmassa, &
     dtemp(2), &
     dtemp(3),velext
  endif
  
  call bcast_world_d(lencorrmassa)
  call bcast_world_d(velext)
  
  lencorrmassa=(lencorrmassa/lengthscale)
  velext=(velext*tao/lengthscale)
  
  if(idrank==0)then
    read(fileout)mioind
  endif
  call bcast_world_i(mioind)
  
  if(mioind>=1)then
    if(idrank==0)then
      read(fileout)q,v,fve,fvere,Hg,Lrg,Gr,ks,ksre,att,Li,lire
    endif
    call bcast_world_d(q)
    call bcast_world_d(v)
    call bcast_world_d(fve)
    call bcast_world_d(fvere)
    call bcast_world_d(Hg)
    call bcast_world_d(Lrg)
    call bcast_world_d(Gr)
    call bcast_world_d(Ks)
    call bcast_world_d(Ksre)
    call bcast_world_d(att)
    call bcast_world_d(Li)
    call bcast_world_d(lire)
  endif
  
  
  if(mioind>=2)then
    if(idrank==0)then
      read(fileout)naddtrack,nremtrack,linserting
    endif
    call bcast_world_i(naddtrack)
    call bcast_world_i(nremtrack)
    call bcast_world_l(linserting)
  endif
  
  if(lrestartreset)then
  
    if(mioind>=3)then
      if(idrank==0)then
        read(fileout)(dtemp(i),i=1,12)
      endif
    endif
    if(mioind>=4)then
      if(idrank==0)then
        read(fileout)(dtemp(i),i=1,12)
      endif
    endif
    if(mioind>=5)then
      if(idrank==0)then
        read(fileout)(dtemp(i),i=1,12)
      endif
    endif
    if(mioind>=6)then
      if(idrank==0)then
        read(fileout)(itemp(i),i=1,12)
      endif
    endif
    if(mioind>=7)then
      if(idrank==0)then
        read(fileout)(ltemp(i),i=1,15)
      endif
    endif
    if(mioind>=8)then
      if(idrank==0)then
        read(fileout)(itemp(i),i=1,12)
      endif
    endif
    if(mioind>=9)then
      if(idrank==0)then
        read(fileout)(dtemp(i),i=1,12)
      endif
    endif
  
  else
  
    if(mioind>=3)then
      if(idrank==0)then
        read(fileout)addedcharge,removedcharge,countericurr, &
         counterecurr,meanicurr,meanecurr,addedmass,removedmass, &
         counterimass,counteremass,meanimass,meanemass
      endif
      call bcast_world_d(addedcharge)
      call bcast_world_d(removedcharge)
      call bcast_world_d(countericurr)
      call bcast_world_d(counterecurr)
      call bcast_world_d(meanicurr)
      call bcast_world_d(meanecurr)
      call bcast_world_d(addedmass)
      call bcast_world_d(removedmass)
      call bcast_world_d(counterimass)
      call bcast_world_d(counteremass)
      call bcast_world_d(meanimass)
      call bcast_world_d(meanemass)
    endif
    
    if(mioind>=4)then
      if(idrank==0)then
        read(fileout)counterevel,counterevelrel,meanevel,meanevelrel, &
         counterelen,counterecross,meanelen,meanecross,meancputime, &
         counterivel,meanivel,counterlpath
      endif
      call bcast_world_d(counterevel)
      call bcast_world_d(counterevelrel)
      call bcast_world_d(meanevel)
      call bcast_world_d(meanevelrel)
      call bcast_world_d(counterelen)
      call bcast_world_d(counterecross)
      call bcast_world_d(meanelen)
      call bcast_world_d(meanecross)
      call bcast_world_d(meancputime)
      call bcast_world_d(counterivel)
      call bcast_world_d(meanivel)
      call bcast_world_d(counterlpath)
    endif
    
    if(mioind>=5)then
      if(idrank==0)then
        read(fileout)meanlpath,multisteperror,(dtemp(i),i=1,10)
      endif
      call bcast_world_d(meanlpath)
      call bcast_world_d(multisteperror)
    endif
    
    if(mioind>=6)then
      if(idrank==0)then
          read(fileout)ncounterevel,ncounterevelrel,ncountergeom, &
       reprinttime,ncounterivel,ncounterlpath,irefinementdone, &
         nmulstep,nmulstepdone,nmultisteperror,mxnpjet,itemp(1)
      endif
      call bcast_world_i(ncounterevel)
      call bcast_world_i(ncounterevelrel)
      call bcast_world_i(ncountergeom)
      call bcast_world_i(reprinttime)
      call bcast_world_i(ncounterivel)
      call bcast_world_i(ncounterlpath)
      call bcast_world_i(irefinementdone)
      call bcast_world_i(nmulstep)
      call bcast_world_i(nmulstepdone)
      call bcast_world_i(nmultisteperror)
      call bcast_world_i(mxnpjet)
      if(nmulstep>0)lmultiplestep=.true.
    endif
    
    if(mioind>=7)then
      if(idrank==0)then
        read(fileout)lrefinement,ltagbeads,lrefinementthreshold, &
         lrefbeadstart,llenthresholdbead, &
         lrefinementevery,lrefinementstart,lmassavariable,ltemp(1), &
         llencorrmassa,lmassavariable,lfirstmass
      endif
      call bcast_world_l(lrefinement)
      call bcast_world_l(ltagbeads)
      call bcast_world_l(lrefinementthreshold)
      call bcast_world_l(lrefbeadstart)
      call bcast_world_l(llenthresholdbead)
      call bcast_world_l(lrefinementevery)
      call bcast_world_l(lrefinementstart)
      call bcast_world_l(lmassavariable)
      call bcast_world_l(llencorrmassa)
      call bcast_world_l(lmassavariable)
      call bcast_world_l(lfirstmass)
    endif
    
    if(mioind>=8)then
      if(idrank==0)then
        read(fileout)itemp(11),itemp(12),(itemp(i),i=1,10)
      endif
    endif
    
    if(mioind>=9)then
      if(idrank==0)then
        read(fileout)refinementthreshold,refbeadstartfit,lencorrmassa, &
         lenprobmassa,massratio,dtemp(1),oldgaussn,corr,radcorr, &
         lenthresholdbead,refinementevery,refinementstart
      endif
      call bcast_world_d(refinementthreshold)
      call bcast_world_d(refbeadstartfit)
      call bcast_world_d(lencorrmassa)
      call bcast_world_d(lenprobmassa)
      call bcast_world_d(massratio)
      call bcast_world_d(oldgaussn)
      call bcast_world_d(corr)
      call bcast_world_d(radcorr)
      call bcast_world_d(lenthresholdbead)
      call bcast_world_d(refinementevery)
      call bcast_world_d(refinementstart)
    endif
    
  endif
  
  return
    
 end subroutine read_dat_parameter
 
 subroutine read_dat_restart(fileout,k,timesub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for reading the jet geometry 
!     from the restart file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: fileout
  integer, intent(inout) ::k
  double precision, intent(inout) :: timesub
  
  integer :: natms,i,j,sprintdatsub,systypesub,ipoint
  logical :: lstart
  real(4) :: dtemp(10)
  
  
  if(idrank==0)then
    read(fileout)lstart
    read(fileout)natms,sprintdatsub,systypesub,timesub,linserted
  endif
  call bcast_world_d(timesub)
  call bcast_world_l(linserted)

  timesub=(timesub/tao)
  if(idrank==0)read(fileout)doreorder,inpjet,npjet
  
  call bcast_world_i(inpjet)
  call bcast_world_i(npjet)
  
  doreorder=.false.
  call deallocate_jet()
  call allocate_jet(.true.) 
  lreordertrack=.false.
  naddtrack=0
  nremtrack=0
  if(idrank==0)then
    if(sprintdat==5)then
      select case(systype)
      case (1:2)
        jetxx(0:mxnpjet)=0.d0
        jetst(0:mxnpjet)=0.d0
        jetvx(0:mxnpjet)=0.d0
        jetms(0:mxnpjet)=0.d0
        jetch(0:mxnpjet)=0.d0
        jetvl(0:mxnpjet)=0.d0
        do i=inpjet,npjet
          read(fileout)(dtemp(j),j=1,6)
          jetxx(i)=dble(dtemp(1))
          jetst(i)=dble(dtemp(2))
          jetvx(i)=dble(dtemp(3))
          jetms(i)=dble(dtemp(4))
          jetch(i)=dble(dtemp(5))
          jetvl(i)=dble(dtemp(6))
        end do
      case default
        jetxx(0:mxnpjet)=0.d0
        jetyy(0:mxnpjet)=0.d0
        jetzz(0:mxnpjet)=0.d0
        jetst(0:mxnpjet)=0.d0
        jetvx(0:mxnpjet)=0.d0
        jetvy(0:mxnpjet)=0.d0
        jetvz(0:mxnpjet)=0.d0
        jetms(0:mxnpjet)=0.d0
        jetch(0:mxnpjet)=0.d0
        jetvl(0:mxnpjet)=0.d0
        do i=inpjet,npjet
          read(fileout)(dtemp(j),j=1,10)
          jetxx(i)=dble(dtemp(1))
          jetyy(i)=dble(dtemp(2))
          jetzz(i)=dble(dtemp(3))
          jetst(i)=dble(dtemp(4))
          jetvx(i)=dble(dtemp(5))
          jetvy(i)=dble(dtemp(6))
          jetvz(i)=dble(dtemp(7))
          jetms(i)=dble(dtemp(8))
          jetch(i)=dble(dtemp(9))
          jetvl(i)=dble(dtemp(10))
         end do
      end select
    elseif(sprintdat==6)then
      select case(systype)
      case (1:2)
        jetxx(0:mxnpjet)=0.d0
        jetst(0:mxnpjet)=0.d0
        jetvx(0:mxnpjet)=0.d0
        jetms(0:mxnpjet)=0.d0
        jetch(0:mxnpjet)=0.d0
        jetvl(0:mxnpjet)=0.d0
        do i=inpjet,npjet
          read(fileout)(dtemp(j),j=1,6),jetbd(i)
          jetxx(i)=dble(dtemp(1))
          jetst(i)=dble(dtemp(2))
          jetvx(i)=dble(dtemp(3))
          jetms(i)=dble(dtemp(4))
          jetch(i)=dble(dtemp(5))
          jetvl(i)=dble(dtemp(6))
        end do
      case default
        jetxx(0:mxnpjet)=0.d0
        jetyy(0:mxnpjet)=0.d0
        jetzz(0:mxnpjet)=0.d0
        jetst(0:mxnpjet)=0.d0
        jetvx(0:mxnpjet)=0.d0
        jetvy(0:mxnpjet)=0.d0
        jetvz(0:mxnpjet)=0.d0
        jetms(0:mxnpjet)=0.d0
        jetch(0:mxnpjet)=0.d0
        jetvl(0:mxnpjet)=0.d0
        do i=inpjet,npjet
          read(fileout)(dtemp(j),j=1,10),jetbd(i)
          jetxx(i)=dble(dtemp(1))
          jetyy(i)=dble(dtemp(2))
          jetzz(i)=dble(dtemp(3))
          jetst(i)=dble(dtemp(4))
          jetvx(i)=dble(dtemp(5))
          jetvy(i)=dble(dtemp(6))
          jetvz(i)=dble(dtemp(7))
          jetms(i)=dble(dtemp(8))
          jetch(i)=dble(dtemp(9))
          jetvl(i)=dble(dtemp(10))
         end do
      end select
    endif
  endif
  
  if(sprintdat==5)then
    select case(systype)
    case (1:2)
      call bcast_world_darr(jetxx,mxnpjet+1)
      call bcast_world_darr(jetst,mxnpjet+1)
      call bcast_world_darr(jetvx,mxnpjet+1)
      call bcast_world_darr(jetms,mxnpjet+1)
      call bcast_world_darr(jetch,mxnpjet+1)
      call bcast_world_darr(jetvl,mxnpjet+1)
    case default
      call bcast_world_darr(jetxx,mxnpjet+1)
      call bcast_world_darr(jetyy,mxnpjet+1)
      call bcast_world_darr(jetzz,mxnpjet+1)
      call bcast_world_darr(jetst,mxnpjet+1)
      call bcast_world_darr(jetvx,mxnpjet+1)
      call bcast_world_darr(jetvy,mxnpjet+1)
      call bcast_world_darr(jetvz,mxnpjet+1)
      call bcast_world_darr(jetms,mxnpjet+1)
      call bcast_world_darr(jetch,mxnpjet+1)
      call bcast_world_darr(jetvl,mxnpjet+1)
    end select
  elseif(sprintdat==6)then
    select case(systype)
    case (1:2)
      call bcast_world_darr(jetxx,mxnpjet+1)
      call bcast_world_darr(jetst,mxnpjet+1)
      call bcast_world_darr(jetvx,mxnpjet+1)
      call bcast_world_darr(jetms,mxnpjet+1)
      call bcast_world_darr(jetch,mxnpjet+1)
      call bcast_world_darr(jetvl,mxnpjet+1)
      call bcast_world_larr(jetbd,mxnpjet+1)
    case default
      call bcast_world_darr(jetxx,mxnpjet+1)
      call bcast_world_darr(jetyy,mxnpjet+1)
      call bcast_world_darr(jetzz,mxnpjet+1)
      call bcast_world_darr(jetst,mxnpjet+1)
      call bcast_world_darr(jetvx,mxnpjet+1)
      call bcast_world_darr(jetvy,mxnpjet+1)
      call bcast_world_darr(jetvz,mxnpjet+1)
      call bcast_world_darr(jetms,mxnpjet+1)
      call bcast_world_darr(jetch,mxnpjet+1)
      call bcast_world_darr(jetvl,mxnpjet+1)
      call bcast_world_larr(jetbd,mxnpjet+1)
    end select
  endif
  
  jetfr(:)=.false.
  do ipoint=inpjet,npjet
    if(jetxx(ipoint)>=h)then
      jetfr(ipoint)=.true.
      jetxx(ipoint)=h
    endif
  enddo
  jetxx(npjet)=0.d0
  
  if(lmultiplestep)then
    jetfm(:)=.false.
    lneighlistdo=.true.
  endif
  
  
  return
    
 end subroutine read_dat_restart
 
 subroutine write_dat_parameter(lprintdatsub,fileout,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for writing the parameters on the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lprintdatsub

  integer, intent(in) :: fileout
  double precision, intent(in) :: timesub
  
  integer :: natms,i
  integer, parameter :: sprintdatsub=0
  logical, parameter :: lstart=.true.
  
  integer, parameter :: mioind=9
  
  double precision :: dtemp(15)
  integer :: itemp(15)
  logical :: ltemp(15)
  
  if(idrank/=0)return
  if(.not.lprintdatsub)return

  
  natms=npjet-inpjet+1
  write(fileout)lstart
  write(fileout)natms,sprintdatsub,systype,(timesub*tao),lstart
  
  itemp(:)=0
  write(fileout)itemp(1),systype,integrator,units
  write(fileout)(resolution*lengthscale),(ilength*lengthscale), &
   imassa*(massscale)/(lengthscale**3.d0), &
   icharge*(chargescale)/(lengthscale**3.d0),(icrossec*lengthscale), &
   (istress*G),(ivelocity*lengthscale/tao)
  write(fileout)(pfreq/tao),(pampl*lengthscale), &
   (airdragamp(1)*(lengthscale**2.d0)/(tao**3.d0)), &
   (airdragamp(2)*(lengthscale**2.d0)/(tao**3.d0)), &
   (airdragamp(3)*(lengthscale**2.d0)/(tao**3.d0)),aird,airv
  write(fileout)mu,G,(yieldstress*G),(h*lengthscale),V0,surfacet, &
   (tstep*tao),(consistency*(mu*((tao)**(findex-1.d0)))),findex
  write(fileout)massscale,chargescale,lengthscale,tao
  write(fileout)dtemp(1), &
   (lencorrmassa*lengthscale), &
   dtemp(2), &
   dtemp(3),(velext*lengthscale/tao)
  
  write(fileout)mioind
  write(fileout)q,v,fve,fvere,Hg,Lrg,Gr,ks,ksre,att,Li,lire
  write(fileout)naddtrack,nremtrack,linserting
  write(fileout)addedcharge,removedcharge,countericurr,counterecurr, &
   meanicurr,meanecurr,addedmass,removedmass,counterimass, &
   counteremass,meanimass,meanemass
  write(fileout)counterevel,counterevelrel,meanevel,meanevelrel, &
  counterelen,counterecross,meanelen,meanecross,meancputime, &
  counterivel,meanivel,counterlpath
  dtemp(:)=0.d0
  write(fileout)meanlpath,multisteperror,(dtemp(i),i=1,10)
  itemp(:)=0
  if(lmultiplestep)then
    write(fileout)ncounterevel,ncounterevelrel,ncountergeom, &
     reprinttime,ncounterivel,ncounterlpath,irefinementdone,nmulstep, &
     nmulstepdone,nmultisteperror,mxnpjet,itemp(1)
  else
    write(fileout)ncounterevel,ncounterevelrel,ncountergeom, &
     reprinttime,ncounterivel,ncounterlpath,irefinementdone, &
     (itemp(i),i=1,3),mxnpjet,itemp(4)
  endif
  
  ltemp(:)=.false.
  write(fileout)lrefinement,ltagbeads,lrefinementthreshold, &
         lrefbeadstart,llenthresholdbead, &
         lrefinementevery,lrefinementstart,lmassavariable,lstart, &
         llencorrmassa,lmassavariable,lfirstmass
     
  write(fileout)itemp(11),itemp(12),(itemp(i),i=1,10)
     
  write(fileout)refinementthreshold,refbeadstartfit,lencorrmassa, &
         lenprobmassa,massratio,dtemp(1),oldgaussn,corr,radcorr, &
         lenthresholdbead,refinementevery,refinementstart
  
      
  call flush(fileout)
  
  return
    
 end subroutine write_dat_parameter
 
 subroutine write_dat_frame(lprintdatsub,fileout,k,timesub,iprintdat, &
  inpjet,npjet,sprintdat,systype, &
  linserted)
  
!***********************************************************************
!     
!     JETSPIN subroutine for writing the jet geometry on the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lprintdatsub,linserted
  integer, intent(in) :: fileout,k,iprintdat,inpjet,npjet,sprintdat, &
   systype
  double precision, intent(in) :: timesub
  
  integer :: natms,i
  logical, parameter :: lstart=.true.
  
  if(idrank/=0)return
  if(.not.lprintdatsub)return
  if(mod(k,iprintdat)/=0)return
  
  natms=npjet-inpjet+1
  write(fileout)lstart
  write(fileout)natms,sprintdat,systype,(timesub*tao),linserted
  if(sprintdat==1)then
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4),real(jetvx(i),4)
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
        real(jetzz(i),4), &
        real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
        real(jetvz(i),4)
      end do
    end select
  elseif(sprintdat==2)then
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4), &
        real(jetvx(i),4), &
        real(jetms(i),4),real(jetch(i),4)
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
         real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4)
      end do
    end select
  elseif(sprintdat==3)then
    write(fileout)doreorder,inpjet,npjet
    doreorder=.false.
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4), &
        real(jetvx(i),4), &
        real(jetms(i),4),real(jetch(i),4),i
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
         real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),i
      end do
    end select
  elseif(sprintdat==4)then
    write(fileout)lreordertrack,naddtrack,nremtrack
    lreordertrack=.false.
    naddtrack=0
    nremtrack=0
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4), &
        real(jetvx(i),4), &
        real(jetms(i),4),real(jetch(i),4),jetlb(i)
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
         real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),jetlb(i)
      end do
    end select
  elseif(sprintdat==5)then
    write(fileout)doreorder,inpjet,npjet
    lreordertrack=.false.
    naddtrack=0
    nremtrack=0
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4), &
        real(jetvx(i),4), &
        real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4)
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
         real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4)
      end do
    end select
  elseif(sprintdat==6)then
    write(fileout)doreorder,inpjet,npjet
    lreordertrack=.false.
    naddtrack=0
    nremtrack=0
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4), &
        real(jetvx(i),4), &
        real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4),jetbd(i)
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
         real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4),jetbd(i)
      end do
    end select
  elseif(sprintdat==7)then
    write(fileout)doreorder,inpjet,npjet
    lreordertrack=.false.
    naddtrack=0
    nremtrack=0
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetst(i),4), &
        real(jetvx(i),4), &
        real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4),jetbd(i), &
        jetbr(i)
      end do
    case default
      do i=inpjet,npjet
        write(fileout)real(jetxx(i),4),real(jetyy(i),4), &
         real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4),jetbd(i), &
         jetbr(i)
      end do
    end select
  else
    call error(9)
  endif
  
  call flush(fileout)
  
  return
    
 end subroutine write_dat_frame
 
 subroutine write_datrem_frame(lprintdatsub,fileout,k,timesub, &
  iprintdat,inpjet,npjet,sprintdat,systype,lremdat,nremovedsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for writing the jet beads on the binary file
!     which have hit the collector
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lprintdatsub
  integer, intent(in) :: fileout,k,iprintdat,inpjet,npjet,sprintdat, &
   systype,nremovedsub
  double precision, intent(in) :: timesub
  logical, intent(inout) :: lremdat
  
  integer :: natms,i
  logical, parameter :: lstart=.true.
  double precision :: tempmod0,tempcr
  
  if(idrank/=0)return
  
  if(.not.lprintdatrem)return
  
  if(lremdat)then
    do i=inpjet-nremovedsub,inpjet-1
      if(systype==1)then
        tempmod0 = jetxx(i+1)-jetxx(i)
      else
        tempmod0 = dsqrt((jetxx(i+1)-jetxx(i))**2.d0+ &
         (jetyy(i+1)-jetyy(i))**2.d0+ &
         (jetzz(i+1)-jetzz(i))**2.d0)
      endif
      tempcr=dsqrt((jetvl(i))/(Pi*tempmod0))
      write(fileout)lstart
      if(ltagbeads)then
        write(fileout)real(timesub,4),real(jetxx(i),4), &
         real(jetyy(i),4),real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4), &
         real(tempcr,4),real(tempmod0,4),jetbd(i)
      else
        write(fileout)real(timesub,4),real(jetxx(i),4), &
         real(jetyy(i),4),real(jetzz(i),4), &
         real(jetst(i),4),real(jetvx(i),4),real(jetvy(i),4), &
         real(jetvz(i),4), &
         real(jetms(i),4),real(jetch(i),4),real(jetvl(i),4), &
         real(tempcr,4),real(tempmod0,4)
      endif
    enddo
  endif
  
  lremdat=.false.
  
  return
  
 end subroutine write_datrem_frame
 
 subroutine close_dat_file(lprintdatsub,fileout)
 
!***********************************************************************
!     
!     JETSPIN subroutine for closing the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lprintdatsub
  integer, intent(in) :: fileout
  
  
  if(idrank/=0)return
  if(.not.lprintdatsub)return
  
  close(fileout)
  
  return
  
 end subroutine close_dat_file
 
 subroutine close_datrem_file(lprintdatsub,fileout)
 
!***********************************************************************
!     
!     JETSPIN subroutine for closing the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lprintdatsub
  integer, intent(in) :: fileout
  
  
  if(idrank/=0)return
  
  if(.not.lprintdatrem)return
  
  close(fileout)
  
  return
  
 end subroutine close_datrem_file
 
 subroutine open_xyz_file(lprintxyz,fileout,filename)
 
!***********************************************************************
!     
!     JETSPIN subroutine for opening the XYZ formatted output file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lprintxyz
  integer, intent(in) :: fileout
  character(len=*), intent(in) :: filename
  
  if(idrank/=0)return
  if(.not.lprintxyz)return
  
  if(lreadrest)then
    open(fileout,file=trim(filename),form='formatted', &
     status='old',action='write',position='append')
  else
    open(fileout,file=trim(filename),form='formatted', &
     status='replace',action='write')
  endif
  
  return
  
 end subroutine open_xyz_file
 
 subroutine write_xyz_frame(lprintxyz,fileout,k,timesub,iprintxyz, &
  inpjet,npjet,rescale,systype,linserted,totnatms,lprintnpjet)
  
!***********************************************************************
!     
!     JETSPIN subroutine for printing the data on the XYZ 
!     formatted output file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lprintxyz,linserted
  integer, intent(in) :: fileout,k,iprintxyz,inpjet,npjet,systype, &
   totnatms
  double precision, intent(in) :: timesub,rescale
  logical, intent(in) :: lprintnpjet
  
  integer :: natms,i,j
  character(len=8),parameter :: atname='C       '
  double precision :: mytime
  
  if(idrank/=0)return
  if(.not.lprintxyz)return
  if(mod(k,iprintxyz)/=0)return
  
  mytime=timesub*tao
  
  j=0
  if(linserted)then
    natms=npjet-inpjet+1
    select case(systype)
      case (1:2)
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
           0.d0,0.d0
        end do
        do while(j<totnatms)
          j=j+1
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        end do
      case default
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet-1
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
           jetyy(i)*rescale,jetzz(i)*rescale
        end do
        j=j+1
        if(j>totnatms)return
        if(.not.lprintnpjet)then
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        else
          i=npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
            jetyy(i)*rescale,jetzz(i)*rescale
        endif
        do while(j<totnatms)
          j=j+1
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        end do
    end select
  else
    natms=npjet-inpjet
    select case(systype)
      case (1:2)
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet-2
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
           0.d0,0.d0
        end do
        j=j+1
        if(j>totnatms)return
        i=npjet
        write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
         0.d0,0.d0
        do while(j<totnatms)
          j=j+1
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        end do
      case default
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet-2
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
           jetyy(i)*rescale,jetzz(i)*rescale
        end do
        j=j+1
        if(j>totnatms)return
        if(.not.lprintnpjet)then
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        else
          i=npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*rescale, &
            jetyy(i)*rescale,jetzz(i)*rescale
        endif
        do while(j<totnatms)
          j=j+1
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        end do
    end select
  endif
  
  
  
  call flush(fileout)
  
  return
    
 end subroutine write_xyz_frame
 
 subroutine close_xyz_file(lprintxyz,fileout)
 
!***********************************************************************
!     
!     JETSPIN subroutine for closing the XYZ formatted output file
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, intent(in) :: lprintxyz
  integer, intent(in) :: fileout
  
  if(idrank/=0)return
  if(.not.lprintxyz)return
  close(fileout)
  
  return
  
 end subroutine close_xyz_file
 
 subroutine write_xyz_singlefile(filename,lprintxyzsingsub,fileout,k, &
  timesub,iprintxyzsub,inpjet,npjet,rescale,systype,linserted, &
  lprintnpjet)
  
!***********************************************************************
!     
!     JETSPIN subroutine for printing the data on the XYZ 
!     formatted output file in a single file for frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lprintxyzsingsub,linserted
  integer, intent(in) :: fileout,k,iprintxyzsub,inpjet,npjet,systype
  double precision, intent(in) :: timesub,rescale
  logical, intent(in) :: lprintnpjet
  
  
  integer :: i,j,totnatms
  character(len=8),parameter :: atname='C       '
  
  character(len=150) :: nomefile
  integer :: l
  double precision :: mytime,myrescale
  
  if(idrank/=0)return
  if(.not.lprintxyzsingsub)return
  if(mod(k,iprintxyzsub)/=0)return
  
  mytime=timesub*tao
  myrescale=rescale*lengthscale
  
  l=k/iprintxyzsub
  nomefile=trim(filename)//write_fmtnumb(l)//'.xyz'
  
  open(fileout,file=trim(nomefile),status='replace',action='write')
  
  if(linserted)then
    totnatms=npjet-inpjet+1
    select case(systype)
      case (1:2)
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
           0.d0,0.d0
        end do
      case default
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet-1
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
           jetyy(i)*myrescale,jetzz(i)*myrescale
        end do
        if(.not.lprintnpjet)then
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        else
          i=npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
            jetyy(i)*myrescale,jetzz(i)*myrescale
        endif
    end select
  else
    totnatms=npjet-inpjet
    select case(systype)
      case (1:2)
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet-2
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
           0.d0,0.d0
        end do
        i=npjet
        write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
         0.d0,0.d0
      case default
        write(fileout,'(i10)') totnatms
        write(fileout,'(a,g20.10,i10)')'frame at time',mytime,k
        do i=inpjet,npjet-2
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
           jetyy(i)*myrescale,jetzz(i)*myrescale
        end do
        if(.not.lprintnpjet)then
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        else
          i=npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)*myrescale, &
            jetyy(i)*myrescale,jetzz(i)*myrescale
        endif
    end select
  endif
  
  
  close(fileout)
  
  return
    
 end subroutine write_xyz_singlefile
 
 subroutine write_pdb_singlefile(filename,lprintpdbsingsub,fileout,k, &
  timesub,iprintpdbsub,inpjet,npjet,rescale,systype,linserted, &
  lprintnpjet)
  
!***********************************************************************
!     
!     JETSPIN subroutine for printing the data on the PDB 
!     formatted output file in a single file for frame
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lprintpdbsingsub,linserted
  integer, intent(in) :: fileout,k,iprintpdbsub,inpjet,npjet,systype
  double precision, intent(in) :: timesub,rescale
  logical, intent(in) :: lprintnpjet
  
  
  integer :: i,j,totnatms
  character(len=8),parameter :: atname='C       '
  
  character(len=150) :: nomefile
  integer :: l
  double precision :: mytime,myrescale
  
  if(idrank/=0)return
  if(.not.lprintpdbsingsub)return
  if(mod(k,iprintpdbsub)/=0)return
  
  mytime=timesub*tao
  myrescale=rescale*lengthscale
  
  l=k/iprintpdbsub
  nomefile=trim(filename)//write_fmtnumb(l)//'.pdb'
  
  open(fileout,file=trim(nomefile),status='replace',action='write')
  
  if(lpdbtagbeads)then
    if(linserted)then
      totnatms=npjet-inpjet+1
      select case(systype)
        case (1:2)
          j=0
          do i=inpjet,npjet
            j=j+1
            if(jetbd(i))then
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'N  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
            else
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
            endif
          end do
        case default
          j=0
          do i=inpjet,npjet-1
            j=j+1
            if(jetbd(i))then
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'N  ','ALA',j,jetxx(i)*myrescale, &
               jetyy(i)*myrescale,jetzz(i)*myrescale
            else
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
               jetyy(i)*myrescale,jetzz(i)*myrescale
            endif
          end do
          if(.not.lprintnpjet)then
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,0.d0,0.d0,0.d0
          else
            i=npjet
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
             jetyy(i)*myrescale,jetzz(i)*myrescale
          endif
      end select
    else
      totnatms=npjet-inpjet
      select case(systype)
        case (1:2)
          j=0
          do i=inpjet,npjet-2
            j=j+1
            if(jetbd(i))then
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'N  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
            else
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
            endif
          end do
          i=npjet
          j=j+1
          write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
           'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
        case default
          j=0
          do i=inpjet,npjet-2
            j=j+1
            if(jetbd(i))then
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'N  ','ALA',j,jetxx(i)*myrescale, &
               jetyy(i)*myrescale,jetzz(i)*myrescale
            else
              write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
               'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
               jetyy(i)*myrescale,jetzz(i)*myrescale
            endif
          end do
          if(.not.lprintnpjet)then
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,0.d0,0.d0,0.d0
          else
            i=npjet
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
             jetyy(i)*myrescale,jetzz(i)*myrescale
          endif
      end select
    endif
  else
    if(linserted)then
      totnatms=npjet-inpjet+1
      select case(systype)
        case (1:2)
          j=0
          do i=inpjet,npjet
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
          end do
        case default
          j=0
          do i=inpjet,npjet-1
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
             jetyy(i)*myrescale,jetzz(i)*myrescale
          end do
          if(.not.lprintnpjet)then
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,0.d0,0.d0,0.d0
          else
            i=npjet
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
             jetyy(i)*myrescale,jetzz(i)*myrescale
          endif
      end select
    else
      totnatms=npjet-inpjet
      select case(systype)
        case (1:2)
          j=0
          do i=inpjet,npjet-2
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
          end do
          i=npjet
          j=j+1
          write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
           'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale,0.d0,0.d0
        case default
          j=0
          do i=inpjet,npjet-2
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
             jetyy(i)*myrescale,jetzz(i)*myrescale
          end do
          if(.not.lprintnpjet)then
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,0.d0,0.d0,0.d0
          else
            i=npjet
            j=j+1
            write(fileout,'(a4,2x,i5,2x,a3,1x,a3,1x,i5,4x,3f8.2)') &
             'ATOM',j,'C  ','ALA',j,jetxx(i)*myrescale, &
             jetyy(i)*myrescale,jetzz(i)*myrescale
          endif
      end select
    endif
  endif 
  close(fileout)
  
  call write_psf_singlefile(fileout,filename,l,totnatms)
  
  return
    
 end subroutine write_pdb_singlefile
 
 subroutine write_psf_singlefile(psfunit,psfname,k,numpoints)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing a PSF formatted output file
!     for frame containing all the bonds for any consecutive pair
!     of jet beads
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: psfunit,k
  character(len=*), intent(in) :: psfname
  integer,intent(in) :: numpoints
  
  character(len=150) :: nomefile
  integer :: i,j,ismax,l,m,p,ipoint
  character*1, allocatable, dimension(:) :: sstring
  character(len=10) :: r_char,r2_char
  character(len=46) :: fstringout
  logical :: ladvance
  
  nomefile=trim(psfname)//write_fmtnumb(k)//'.psf'
  
  open(unit=psfunit,file=trim(nomefile),status='replace',action='write')
  
  write(psfunit,'(a,/)')'PSF'
  write(psfunit,'(i8,1x,a,/,/)')1,'!TITLE'
  write(psfunit,'(i8,1x,a)')numpoints,'!NATOM'
  
  if(lpdbtagbeads)then
    l=0
    ipoint=inpjet
    do i=1,numpoints
      if(i>=100000)then
        write (r_char,'(i10)')i
        allocate(sstring(6))
        r2_char=adjustl(r_char)
        do j=1,6
          sstring(j)=r2_char(j:j)
        enddo
        ismax=5
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      elseif(i>=10000)then
        write (r_char,'(i10)')i
        allocate(sstring(5))
        r2_char=adjustl(r_char)
        do j=1,5
          sstring(j)=r2_char(j:j)
        enddo
        ismax=5
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      elseif(i>=1000)then
        write (r_char,'(i10)')i
        allocate(sstring(4))
        r2_char=adjustl(r_char)
        do j=1,4
          sstring(j)=r2_char(j:j)
        enddo
        ismax=4
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      else
        write (r_char,'(i10)')i
        allocate(sstring(3))
        r2_char=adjustl(r_char)
        do j=1,3
          sstring(j)=r2_char(j:j)
        enddo
        ismax=3
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      endif
      if(jetbd(ipoint))then
        write(psfunit,fstringout) &
         i,'PROT',(sstring(j),j=1,ismax),'ALA  N    N',0.0,14.007,0
      else
        write(psfunit,fstringout) &
         i,'PROT',(sstring(j),j=1,ismax),'ALA  C    C',0.0,12.011,0
      endif
      deallocate(sstring)
      l=l+1
      ipoint=ipoint+1
    enddo
  else
    l=0
    do i=1,numpoints
      if(i>=100000)then
        write (r_char,'(i10)')i
        allocate(sstring(6))
        r2_char=adjustl(r_char)
        do j=1,6
          sstring(j)=r2_char(j:j)
        enddo
        ismax=5
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      elseif(i>=10000)then
        write (r_char,'(i10)')i
        allocate(sstring(5))
        r2_char=adjustl(r_char)
        do j=1,5
          sstring(j)=r2_char(j:j)
        enddo
        ismax=5
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      elseif(i>=1000)then
        write (r_char,'(i10)')i
        allocate(sstring(4))
        r2_char=adjustl(r_char)
        do j=1,4
          sstring(j)=r2_char(j:j)
        enddo
        ismax=4
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      else
        write (r_char,'(i10)')i
        allocate(sstring(3))
        r2_char=adjustl(r_char)
        do j=1,3
          sstring(j)=r2_char(j:j)
        enddo
        ismax=3
        write(fstringout,'(a,i2,a)')'(i8,1x,a4,1x,',ismax, &
         'a1,2x,a,6x,f8.6,8x,f7.4,2x,i10)'
      endif
      write(psfunit,fstringout) &
       i,'PROT',(sstring(j),j=1,ismax),'ALA  C    C',0.0,12.011,0
      deallocate(sstring)
      l=l+1
    enddo
  endif
  write(psfunit,*)
  write(psfunit,'(i8,1x,a)')l-1,'!NBOND'
  m=1
  p=1
  ladvance=.true.
  do while(m<=l)
    if(p==8)then
      write(psfunit,'(i8)')m
      p=1
    else
      write(psfunit,'(i8)',advance="no")m
      p=p+1
    endif
    if(ladvance)then
      m=m+1
      ladvance=.false.
    else
      if(m==l)exit
      ladvance=.true.
    endif
  enddo
  
  close(psfunit)
  
  return
  
 end subroutine write_psf_singlefile
 
 subroutine write_stat_dat(lprintstatdat,fileout,filename,k)
 
  implicit none
  
!***********************************************************************
!     
!     JETSPIN subroutine for opening and printing all the output data
!     on the binary file 'statdat.dat'
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  logical, intent(in) :: lprintstatdat
  integer, intent(in) :: fileout,k
  character(len=*), intent(in) :: filename
  
  logical, save :: lfirst=.true.
  real(4) :: rtemp(nmaxstatdata)
  integer :: i
  
  if(idrank/=0)return
  if(.not.lprintstatdat)return
  
  if(lfirst)then
    lfirst=.false.
    if(lreadrest)then
      open(fileout,file=trim(filename),form='unformatted', &
       status='old',action='write',position='append')
    else
      open(fileout,file=trim(filename),form='unformatted', &
       status='replace',action='write')
    endif
  endif
  
  do i=1,nmaxstatdata
    rtemp(i)=real(statdata(i),4)
  enddo
  
  write(fileout)k,nmaxstatdata,(rtemp(i),i=1,nmaxstatdata)
  
  return
  
 end subroutine write_stat_dat
 
 subroutine set_sprintdat(myoutdata)
 
!***********************************************************************
!     
!     JETSPIN subroutine for setting the sprintdat variable
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2016
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(out) :: myoutdata
  
  if(ltagbeads)then
    myoutdata=6
  else
    myoutdata=5
  endif
  
  return
  
 end subroutine set_sprintdat
 
 end module io_mod


