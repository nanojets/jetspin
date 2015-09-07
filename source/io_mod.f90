
 module io_mod
 
!***********************************************************************
!     
!     JETSPIN module for input/output routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
 use version_mod
 use parse_module
 use error_mod
 use nanojet_mod,           only : airdragamp,doreorder,tao,aird,airv,&
                             chargescale,consistency,findex,g,h,&
                             icharge,ichargedev,icrossec,ilength,&
                             imassa,imassadev,inpjet,istress,ivelocity,&
                             lencorrcharge,lencorrmassa,lengthscale,&
                             massscale,mu,ncutoff,npjet,pampl,pfreq,&
                             resolution,surfacet,systype,tstep,units,&
                             V0,velext,yieldstress,dcutoff,lairdrag,&
                             lairvel,lchargevariable,lchargevariable,&
                             lconsistency,lfindex,lgravity,lHBfluid,&
                             lilengt,liniperturb,linserting,listress,&
                             livelocity,lmassavariable,lmyseed,&
                             lncutoff,lnpjet,lremove,lresolution,&
                             lxyzrescale,lyieldstress,myseed,&
                             xyzrescale,laird,lairdragamp,lairv,lg,lh,&
                             licharge,lichargedev,licrossec,limass,&
                             limassadev,llencorrcharge,llencorrmassa,&
                             lmu,lpampl,lpfreq,lsurfacet,lsystype,&
                             ltstep,lunits,lv0,att,fve,gr,hg,ks,li,lrg,&
                             q,v,jetxx,jetyy,jetzz,jetvx,jetvy,jetvz,&
                             jetst,jetms,jetch
 use integrator_mod,        only : integrator,endtime,lendtime,&
                             lintegrator
 use statistic_mod,         only : nmaxstatdata,reprinttime,statdata,&
                             compute_statistic

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
 logical, public, save :: lprintdat=.false.
 logical, public, save :: lprintstatdat=.false.
 double precision, save :: printtime
 logical, public, save :: lprinttime=.false.
 logical, public, save :: lsprintdat=.false.
 integer, public, save :: iprinttime=0
 integer, public, save :: eprinttime=0
 integer, public, save :: iprintxyz=0
 integer, public, save :: iprintdat=0
 integer, public, save :: sprintdat=1
 integer, public, save :: maxnumxyz=100
 double precision, save :: printxyz
 double precision, save :: printdat
 double precision, public, allocatable, save :: xprint(:),xprint2(:)
 character(len=11),save :: namefile
 
 character(len=20) , public, allocatable :: printarg(:)
 character(len=20) , public, allocatable :: printarg2(:)
 
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
 public :: open_xyz_file
 public :: write_xyz_frame
 public :: close_xyz_file

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
!     last modification March 2015
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
  write(iu,of)"*    Version 1.00  (July 2015)                                                *"
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
  write(iu,of)"*    Computer Physics Communications, 2015, doi:10.1016/j.cpc.2015.08.013.    *"
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
!     last modification March 2015
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
  write(iu,'(a,g20.10,a,g20.10,a)')"variable mass deviation    = ", &
   imassadev," = ",imassadev*massscale/(lengthscale**3.d0)," g cm^-3"
  write(iu,'(a,g20.10,a,g20.10,a)')"variable mass correlation  = ", &
   lencorrmassa," = ",lencorrmassa*lengthscale," cm"
  endif
  if(lchargevariable)then
  write(iu,'(a,g20.10,a,g20.10,a)')"variable charge deviation  = ", &
   ichargedev," = ",ichargedev*chargescale/(lengthscale**3.d0), &
    " statC cm^-3"
  write(iu,'(a,g20.10,a,g20.10,a)')"variable charge correlation= ", &
   lencorrcharge," = ",lencorrcharge*lengthscale," cm"
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
  write(iu,'(a,g20.10)')"Q    = ",q
  write(iu,'(a,g20.10)')"V    = ",v
  write(iu,'(a,g20.10)')"Fve  = ",fve
  write(iu,'(a,g20.10)')"H    = ",Hg
  write(iu,'(a,g20.10)')"Lstep= ",Lrg
  if(lgravity)then
    write(iu,'(a,g20.10)')"Fg   = ",Gr
  endif
  if(systype/=1)then
  write(iu,'(a,g20.10)')"A    = ",ks
  endif
  if(liniperturb)then
  write(iu,'(a,g20.10)')"Ks   = ",pfreq
  endif
  if(lairdrag)then
  write(iu,'(a,g20.10)')"At   = ",att
  write(iu,'(a,g20.10)')"Li   = ",Li
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
  elseif(printcodsub(iarg)==19)then
    legendobs='mlc  =  mutual length segment at the collector   '
  endif
  legendobs=adjustl(legendobs)
  
  return
  
 end function legendobs
 
 subroutine print_internal_units(iu)
  
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
!     last modification March 2015
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
  ncutoff=0
  myseed=1
  sprintdat=1
  aird=0.d0
  airv=1.d0
  velext=0.d0
  findex=1.d0
  consistency=1.d0 
  
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
  lncutoff=.false.
  lairdrag=.false.
  lairdragamp=.false.
  lremove=.false.
  lmyseed=.false.
  lprintxyz=.false.
  lxyzrescale=.false.
  lprintdat=.false.
  laird=.false.
  lairv=.false.
  lairvel=.false.
  lHBfluid=.false.
  lgravity=.false.
  lprintstatdat=.false.
  
  lprintlisterror=.false.
  lprintlisterror2=.false.
  lfoundprint=.false.
  lsprintdat=.false.
  
  ltestread=.false.
  
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
      elseif(findstring('cutoff',directive,inumchar,maxlen))then
        lncutoff=.true.
        dcutoff=dblstr(directive,maxlen,inumchar)
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
          V0=dblstr(directive,maxlen,inumchar)
          lV0=.true.
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
          else
            printxyz=dblstr(directive,maxlen,inumchar)
            lprintxyz=.true.
          endif
        elseif(findstring('binary',directive,inumchar,maxlen))then
          lprintdat=.true.
          printdat=dblstr(directive,maxlen,inumchar)
          if(findstring('style',directive,inumchar,maxlen))then
            lsprintdat=.true.
            sprintdat=intstr(directive,maxlen,inumchar)
          endif
        elseif(findstring('rescalexyz',directive,inumchar,maxlen))then
          lxyzrescale=.true.
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
      elseif(findstring('variable',directive,inumchar,maxlen))then
        if(findstring('mass',directive,inumchar,maxlen))then
          if(findstring('yes',directive,inumchar,maxlen))then
            lmassavariable=.true.
          elseif(findstring('deviation',directive,inumchar,maxlen))then
            limassadev=.true.
            imassadev=dblstr(directive,maxlen,inumchar)
          elseif(findstring('correlation',directive,inumchar,maxlen) &
           )then
            llencorrmassa=.true.
            lencorrmassa=dblstr(directive,maxlen,inumchar)
          elseif(findstring('no',directive,inumchar,maxlen))then
            lmassavariable=.false.
          else
            call warning(61,dble(iline),redstring)
          endif
        elseif(findstring('charge',directive,inumchar,maxlen) &
         )then
          if(findstring('yes',directive,inumchar,maxlen))then
            lchargevariable=.true.
          elseif(findstring('deviation',directive,inumchar,maxlen))then
            lichargedev=.true.
            ichargedev=dblstr(directive,maxlen,inumchar)
          elseif(findstring('correlation',directive,inumchar,maxlen) &
           )then
            llencorrcharge=.true.
            lencorrcharge=dblstr(directive,maxlen,inumchar)
          elseif(findstring('no',directive,inumchar,maxlen))then
            lmassavariable=.false.
          else
            call warning(61,dble(iline),redstring)
            ltestread=.true.
          endif
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
  call bcast_world_l(lncutoff)
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
  call bcast_world_d(xyzrescale)
  call bcast_world_l(lxyzrescale)
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
  call bcast_world_l(lconsistency)
  call bcast_world_l(lfindex)
  call bcast_world_d(consistency)
  call bcast_world_d(findex)
  call bcast_world_l(lmassavariable)
  call bcast_world_d(imassadev)
  call bcast_world_d(lencorrmassa)
  call bcast_world_l(lchargevariable)
  call bcast_world_d(ichargedev)
  call bcast_world_d(lencorrcharge)
  call bcast_world_l(limassadev)
  call bcast_world_l(llencorrmassa)
  call bcast_world_l(lichargedev)
  call bcast_world_l(llencorrcharge)
  call bcast_world_l(lprintstatdat)
  call bcast_world_l(ltestread)
  
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
  
  if(lmassavariable)then
    if(.not.limassadev)then
      call warning(50)
      ltest=.true.
    else
      if(imassadev<0.d0)then
        call warning(54,imassadev)
        ltest=.true.
      endif
    endif
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
  
  if(lchargevariable)then
    if(.not.lichargedev)then
      call warning(52)
      ltest=.true.
    else
      if(ichargedev<0.d0)then
        call warning(56,ichargedev)
        ltest=.true.
      endif
    endif
    if(.not.llencorrcharge)then
      call warning(53)
      ltest=.true.
    else
      if(lencorrcharge<0.d0)then
        call warning(57,lencorrcharge)
        ltest=.true.
      endif
    endif
  endif
  if(lsprintdat)then
    if(sprintdat>3 .or. sprintdat<1)then
      call warning(63,dble(sprintdat))
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
  if(lprintdat)iprintdat=nint(printdat/tstep)
  
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
!     last modification March 2015
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
  
  if(findstring('curc',temps,inumchar,lenstring))then
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
  elseif(findstring('mass',temps,inumchar,lenstring))then
    printcodsub(iarg)=34
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
!     last modification March 2015
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
!     last modification March 2015
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
      if(lxyzrescale)then
      labelsub='print xyz file rescaled by'
      write(6,form4)labelsub,ugualab,xyzrescale,' factor'
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
      labelsub='surface tension of the jet'
      write(6,form4)labelsub,ugualab,surfacet,' g s^-2'
      if(liniperturb)then
      labelsub='nozzle perturbation frequency'
      write(6,form4)labelsub,ugualab,pfreq,' s^-1'
      labelsub='nozzle perturbation amplitude'
      write(6,form4)labelsub,ugualab,pampl,' cm'
      endif
      if(lncutoff)then
      labelsub='cutoff points'
      write(6,form4)labelsub,ugualab,dcutoff,' cm'
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
      if(lmassavariable)then
      write(6,form3)"variable mass yes"
      labelsub='variable mass deviation of jet'
      write(6,form4)labelsub,ugualab,imassadev,' g cm^-3'
      labelsub='variable mass correlation of jet'
      write(6,form4)labelsub,ugualab,lencorrmassa,' cm'
      endif
      if(lchargevariable)then
      write(6,form3)"variable charge yes"
      labelsub='variable charge deviation of jet'
      write(6,form4)labelsub,ugualab,ichargedev,' C L^-1'
      labelsub='variable charge correlation of jet'
      write(6,form4)labelsub,ugualab,lencorrcharge,' cm'
      endif
    endif
    icharge=icharge/1000.d0*2997919999.93d0
    ichargedev=ichargedev/1000.d0*2997919999.93d0
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
      if(lxyzrescale)then
      labelsub='print xyz file rescaled by'
      write(6,form4)labelsub,ugualab,xyzrescale,' factor'
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
      labelsub='surface tension of the jet'
      write(6,form4)labelsub,ugualab,surfacet,' g s^-2'
      if(liniperturb)then
      labelsub='nozzle perturbation frequency'
      write(6,form4)labelsub,ugualab,pfreq,' s^-1'
      labelsub='nozzle perturbation amplitude'
      write(6,form4)labelsub,ugualab,pampl,' cm'
      endif
      if(lncutoff)then
      labelsub='cutoff points'
      write(6,form4)labelsub,ugualab,dcutoff,' cm'
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
      if(lmassavariable)then
      write(6,form3)"variable mass yes"
      labelsub='variable mass deviation of jet'
      write(6,form4)labelsub,ugualab,imassadev,' g cm^-3'
      labelsub='variable mass correlation of jet'
      write(6,form4)labelsub,ugualab,lencorrmassa,' cm'
      endif
      if(lchargevariable)then
      write(6,form3)"variable charge yes"
      labelsub='variable charge deviation of jet'
      write(6,form4)labelsub,ugualab,ichargedev,' statC cm^-3'
      labelsub='variable charge correlation of jet'
      write(6,form4)labelsub,ugualab,lencorrcharge,' cm'
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
  call compute_statistic(tempint,timesub)
  
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
    open(unit=outp,file=soutp,status='replace',action='write')    
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
 
 subroutine open_dat_file(lprintdat,fileout,filename)
 
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
  
  logical, intent(in) :: lprintdat
  integer, intent(in) :: fileout
  character(len=*), intent(in) :: filename
  
  if(idrank/=0)return
  if(.not.lprintdat)return
  
  open(fileout,file=filename,form='unformatted',status='replace', &
   action='write')
  
  return
  
 end subroutine open_dat_file
 
 subroutine write_dat_parameter(lprintdat,fileout,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for writing the parameters on the binary file
!     ONLY FOR DEVELOPERS
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical, intent(in) :: lprintdat

  integer, intent(in) :: fileout
  double precision, intent(in) :: timesub
  
  integer :: natms,i
  integer, parameter :: sprintdatsub=0
  logical, parameter :: lstart=.true.
  
  if(idrank/=0)return
  if(.not.lprintdat)return

  
  natms=npjet-inpjet+1
  write(fileout)lstart
  write(fileout)natms,sprintdatsub,systype,(timesub*tao),lstart
  
  write(fileout)ncutoff,systype,integrator,units
  write(fileout)(resolution*lengthscale),(ilength*lengthscale), &
   imassa*(massscale),(icharge*chargescale),(icrossec*lengthscale), &
   (istress*G),(ivelocity*lengthscale/tao)
  write(fileout)(pfreq/tao),(pampl*lengthscale), &
   (airdragamp(1)*(lengthscale**2.d0)/(tao**3.d0)), &
   (airdragamp(2)*(lengthscale**2.d0)/(tao**3.d0)), &
   (airdragamp(3)*(lengthscale**2.d0)/(tao**3.d0)),aird,airv
  write(fileout)mu,G,(yieldstress*G),(h*lengthscale),V0,surfacet, &
   (tstep*tao),(consistency*(mu*((tao)**(findex-1.d0)))),findex
  write(fileout)massscale,chargescale,lengthscale,tao
  write(fileout)(imassadev*massscale/(lengthscale**3.d0)), &
   (lencorrmassa*lengthscale), &
   (ichargedev*chargescale/(lengthscale**3.d0)), &
   (lencorrcharge*lengthscale),(velext*lengthscale/tao)
  
  
  call flush(fileout)
  
  return
    
 end subroutine write_dat_parameter
 
 subroutine write_dat_frame(lprintdat,fileout,k,timesub,iprintdat, &
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
  
  logical, intent(in) :: lprintdat,linserted
  integer, intent(in) :: fileout,k,iprintdat,inpjet,npjet,sprintdat, &
   systype
  double precision, intent(in) :: timesub
  
  integer :: natms,i
  logical, parameter :: lstart=.true.
  
  if(idrank/=0)return
  if(.not.lprintdat)return
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
  else
    call error(9)
  endif
  
  call flush(fileout)
  
  return
    
 end subroutine write_dat_frame
 
 subroutine close_dat_file(lprintdat,fileout)
 
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
  
  logical, intent(in) :: lprintdat
  integer, intent(in) :: fileout
  
  logical, parameter :: lend=.false.
  
  if(idrank/=0)return
  if(.not.lprintdat)return
  write(fileout)lend
  close(fileout)
  
  return
  
 end subroutine close_dat_file
 
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
  
  open(fileout,file=filename,status='replace',action='write')
  
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
  
  if(idrank/=0)return
  if(.not.lprintxyz)return
  if(mod(k,iprintxyz)/=0)return
  j=0
  if(linserted)then
    natms=npjet-inpjet+1
    select case(systype)
      case (1:2)
        write(fileout,*) totnatms
        write(fileout,'(a,g20.10)')'frame at time',timesub
        do i=inpjet,npjet
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
           0.d0,0.d0
        end do
        do while(j<totnatms)
          j=j+1
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        end do
      case default
        write(fileout,*) totnatms
        write(fileout,'(a,g20.10)')'frame at time',timesub
        do i=inpjet,npjet-1
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
           jetyy(i)/rescale,jetzz(i)/rescale
        end do
        j=j+1
        if(j>totnatms)return
        if(.not.lprintnpjet)then
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        else
          i=npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
            jetyy(i)/rescale,jetzz(i)/rescale
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
        write(fileout,*) totnatms
        write(fileout,'(a,g20.10)')'frame at time',timesub
        do i=inpjet,npjet-2
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
           0.d0,0.d0
        end do
        j=j+1
        if(j>totnatms)return
        i=npjet
        write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
         0.d0,0.d0
        do while(j<totnatms)
          j=j+1
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        end do
      case default
        write(fileout,*) totnatms
        write(fileout,'(a,g20.10)')'frame at time',timesub
        do i=inpjet,npjet-2
          j=j+1
          if(j>totnatms)return
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
           jetyy(i)/rescale,jetzz(i)/rescale
        end do
        j=j+1
        if(j>totnatms)return
        if(.not.lprintnpjet)then
          write(fileout,"(a8,3f10.5)") atname,0.d0,0.d0,0.d0
        else
          i=npjet
          write(fileout,"(a8,3f10.5)") atname,jetxx(i)/rescale, &
            jetyy(i)/rescale,jetzz(i)/rescale
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
    open(fileout,file=filename,form='unformatted',status='replace', &
     action='write')
  endif
  
  do i=1,nmaxstatdata
    rtemp(i)=real(statdata(i),4)
  enddo
  
  write(fileout)k,nmaxstatdata,(rtemp(i),i=1,nmaxstatdata)
  
  return
  
 end subroutine write_stat_dat

 end module io_mod


