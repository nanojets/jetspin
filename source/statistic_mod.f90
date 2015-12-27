
 module statistic_mod
 
!***********************************************************************
!     
!     JETSPIN module containing statistical data routines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
 use version_mod,           only : time_world
 use utility_mod,           only : Pi
 use nanojet_mod,           only : icrossec,inpjet,resolution,systype,&
                             insvx,insvy,insvz,jetxx,jetyy,jetzz,jetvx,&
                             jetvy,jetvz,jetst,jetms,jetch,npjet,g,h,&
                             chargescale,jetcr,lengthscale,massscale,&
                             tao,jetpt,jetvl
 use dynamic_refinement_mod,only : irefinementdone
 use support_functions_mod, only : compute_crosssec,compute_length_path
 
 implicit none
 
 private
 
 integer, public, parameter :: nmaxstatdata=40
 
 double precision, public, save, dimension(nmaxstatdata) :: statdata
 
 double precision, public, save :: addedcharge=0.d0
 double precision, public, save :: removedcharge=0.d0
 
 double precision, public, save :: countericurr=0.d0 
 double precision, public, save :: counterecurr=0.d0
 double precision, public, save :: meanicurr=0.d0
 double precision, public, save :: meanecurr=0.d0
 
 double precision, public, save :: addedmass=0.d0
 double precision, public, save :: removedmass=0.d0
 
 double precision, public, save :: counterimass=0.d0
 double precision, public, save :: counteremass=0.d0
 double precision, public, save :: meanimass=0.d0
 double precision, public, save :: meanemass=0.d0
 
 integer, public, save :: ncounterevel=0
 integer, public, save :: ncounterevelrel=0
 double precision, public, save :: counterevel=0.d0
 double precision, public, save :: counterevelrel=0.d0
 double precision, public, save :: meanevel=0.d0
 double precision, public, save :: meanevelrel=0.d0
 
 integer, public, save :: ncountergeom=0
 double precision, public, save :: counterelen=0.d0
 double precision, public, save :: counterecross=0.d0
 double precision, public, save :: meanelen=0.d0
 double precision, public, save :: meanecross=0.d0
 
 double precision, public, save :: meancputime=0.d0
 
 double precision, public, save :: elapsedcputime=0.d0
 
 integer, public, save :: reprinttime=0
 
 integer, public, save :: ncounterivel=0
 double precision, public, save :: counterivel=0.d0
 double precision, public, save :: meanivel=0.d0
 
 integer, public, save :: ncounterlpath=0
 double precision, public, save :: counterlpath=0.d0
 double precision, public, save :: meanlpath=0.d0
 
 public :: statistic_driver
 public :: compute_statistic
 
 contains
 
 subroutine statistic_driver(timesub,tstepsub,nstepsub,nremovedsub, &
  ladd,lrem)
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling calls to
!     subroutines which store statistical data 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in) :: timesub,tstepsub
  integer, intent(in) :: nstepsub,nremovedsub
  logical, intent(in) :: ladd,lrem
  
! update the counter for all the observables
  call store_addedcharge(ladd)
  call store_current_init(ladd)
  call store_addedmass(ladd)
  call store_mass_init(ladd)
  
  call store_removedcharge(lrem,nremovedsub)
  call store_current_end(lrem,nremovedsub)
  call store_removedmass(lrem,nremovedsub)
  call store_mass_end(lrem,nremovedsub)
  call store_vel_end(lrem,nremovedsub)
  call store_velrel_end(lrem,nremovedsub)
  call store_geometry_end(lrem,nremovedsub)
  
  call store_vel_init()
  call store_length_path()
  
  return
  
 end subroutine statistic_driver
 
 subroutine compute_statistic(tempint,timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for controlling subroutine calls to
!     subroutines which compute statistical data 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in) :: tempint,timesub
  
  integer, save :: icount=0
  integer :: ipoint
  
  
! compute all the observables
  call compute_crosssec(jetxx,jetyy,jetzz,jetvl,jetcr)
  call compute_current_init(tempint)
  call compute_mass_init(tempint)
  
  call compute_current_end(tempint)
  call compute_mass_end(tempint)
  call compute_vel_end()
  call compute_velrel_end()
  call compute_geometry_end()
  
  call compute_vel_init()
  call compute_mean_length_path()
  
  call compute_elapsed_cpu_time()
  call compute_cpu_time()
  
! store all the observables in the statdata array to be printed
  ipoint=inpjet
  statdata(1)=timesub*tao
  statdata(2)=jetxx(ipoint)*lengthscale
  statdata(3)=jetyy(ipoint)*lengthscale
  statdata(4)=jetzz(ipoint)*lengthscale
  statdata(5)=jetst(ipoint)*G
  statdata(6)=jetvx(ipoint)*(lengthscale/tao)
  statdata(7)=jetvy(ipoint)*(lengthscale/tao)
  statdata(8)=jetvz(ipoint)*(lengthscale/tao)
  statdata(9)=dsqrt(jetyy(ipoint)**2.d0+jetzz(ipoint)**2.d0)*lengthscale
  statdata(10)=dble(npjet-inpjet)
  statdata(11)=dble(inpjet)
  statdata(12)=dble(npjet)
  statdata(13)=meanicurr*chargescale/tao
  statdata(15)=meanimass*massscale/tao
  if(isnan(meanemass) .or. isnan(meanecross) .or. isnan(meanecurr))then
    statdata(14)=0.d0 
    statdata(16)=0.d0
    statdata(17)=0.d0
    statdata(18)=0.d0
    statdata(19)=0.d0
    statdata(20)=0.d0
    statdata(21)=0.d0
  else
    statdata(14)=meanecurr*chargescale/tao
    statdata(16)=meanemass*massscale/tao
    statdata(17)=meanevel*(lengthscale/tao)
    statdata(18)=meanevelrel*(1.d0/tao)
    statdata(19)=meanelen*lengthscale
    statdata(20)=meanecross*lengthscale
    statdata(21)=meanecross/icrossec
  endif
  statdata(22)=meancputime
  
! estimate the remaining CPU time
  statdata(23)=meancputime*(reprinttime-icount)
  statdata(24)=timesub
  statdata(25)=jetxx(ipoint)
  statdata(26)=jetyy(ipoint)
  statdata(27)=jetzz(ipoint)
  statdata(28)=jetst(ipoint)
  statdata(29)=jetvx(ipoint)
  statdata(30)=jetvy(ipoint)
  statdata(31)=jetvz(ipoint)
  statdata(32)=dsqrt(jetyy(ipoint)**2.d0+jetzz(ipoint)**2.d0)
  statdata(33)=elapsedcputime
  statdata(34)=jetms(npjet-1)*massscale
  statdata(35)=jetch(npjet-1)*chargescale
  statdata(36)=meanivel*(lengthscale/tao)
  statdata(37)=meanlpath*lengthscale
  statdata(38)=meanlpath*lengthscale/h
  statdata(39)=dble(irefinementdone)
  if(jetxx(ipoint)/=0.d0)then
    statdata(40)=2.d0*datan(dsqrt(jetyy(ipoint)**2.d0+ &
     jetzz(ipoint)**2.d0)/dabs(jetxx(ipoint)))*(180.d0/Pi)
  else
    statdata(40)=0.d0
  endif
  
! update the counter
  icount=icount+1
  
  return
  
 end subroutine compute_statistic
 
 subroutine store_addedcharge(ladd)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the total charge added to
!     the system
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical,intent(in) :: ladd
  
  if(ladd)addedcharge=addedcharge+jetch(npjet-1)
  
  return
  
 end subroutine store_addedcharge
 
 subroutine store_removedcharge(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the total charge removed from
!     the system
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  
  if(lrem)then
    do ipoint=inpjet-nremovedsub,inpjet-1
      removedcharge=removedcharge+jetch(ipoint)
    enddo
  endif
  
  return
  
 end subroutine store_removedcharge
 
 subroutine store_current_init(ladd)
 
!***********************************************************************
!     
!     JETSPIN subroutine for counting the current at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: ladd
  
  if(ladd)countericurr=countericurr+jetch(npjet-1)
  
  return
  
 end subroutine store_current_init
 
 subroutine compute_current_init(tempint)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the current at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: tempint
  
  meanicurr=countericurr/tempint
  
  countericurr=0.d0
  
  return
  
 end subroutine compute_current_init
 
 subroutine store_current_end(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for counting the current at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  
  if(lrem)then
    do ipoint=inpjet-nremovedsub,inpjet-1
      counterecurr=counterecurr+jetch(ipoint)
    enddo
  endif
  
  return
  
 end subroutine store_current_end
 
 subroutine compute_current_end(tempint)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the current at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: tempint
  
  meanecurr=counterecurr/tempint
  
  counterecurr=0.d0
  
  return
  
 end subroutine compute_current_end
 
 subroutine store_addedmass(ladd)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the total mass added to
!     the system
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  logical,intent(in) :: ladd
  
  if(ladd)addedmass=addedmass+jetms(npjet-1)
  
  return
  
 end subroutine store_addedmass
 
 subroutine store_removedmass(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the total mass removed from
!     the system
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
  
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  
  if(lrem)then
    do ipoint=inpjet-nremovedsub,inpjet-1
      removedmass=removedmass+jetms(ipoint)
    enddo
  endif
  
  return
  
 end subroutine store_removedmass
 
 subroutine store_mass_init(ladd)
 
!***********************************************************************
!     
!     JETSPIN subroutine for counting the mass flux at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: ladd
  
  if(ladd)counterimass=counterimass+jetms(npjet-1)
  
  return
  
 end subroutine store_mass_init
 
 subroutine compute_mass_init(tempint)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the mass flux at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: tempint
  
  meanimass=counterimass/tempint
  
  counterimass=0.d0
  
  return
  
 end subroutine compute_mass_init
 
 subroutine store_mass_end(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for counting the mass flux at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  
  if(lrem)then
    do ipoint=inpjet-nremovedsub,inpjet-1
      counteremass=counteremass+jetms(ipoint)
    enddo
  endif
  
  return
  
 end subroutine store_mass_end
 
 subroutine compute_mass_end(tempint)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the mass flux at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: tempint
  
  meanemass=counteremass/tempint
  
  counteremass=0.d0
  
  return
  
 end subroutine compute_mass_end
 
 subroutine store_length_path()
 
!***********************************************************************
!     
!     JETSPIN subroutine for storing the actual path length of the
!     nanofiber between the nozzle and the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision :: tempmod0
  
  
  call compute_length_path(jetxx,jetyy,jetzz, &
  jetpt,tempmod0)
  
  ncounterlpath=ncounterlpath+1
  counterlpath=counterlpath+tempmod0
  
  return
  
 end subroutine store_length_path
 
 subroutine compute_mean_length_path()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the mean path length of the
!     nanofiber between the nozzle and the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  if(ncounterlpath/=0)then
    meanlpath=counterlpath/dble(ncounterlpath)
  else
    meanlpath=0.d0
  endif
  
  ncounterlpath=0
  counterlpath=0.d0
  
  return
  
 end subroutine compute_mean_length_path
 
 subroutine store_vel_init()
 
!***********************************************************************
!     
!     JETSPIN subroutine for storing the actual velocity module vector
!     of the nanofiber at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision :: tempmod0
  
  
  select case(systype)
    case(1:2)
      tempmod0=insvx
    case default
      tempmod0=dsqrt(insvx**2.d0+ &
       insvy**2.d0+insvz**2.d0)
  end select
  
  ncounterivel=ncounterivel+1
  counterivel=counterivel+tempmod0
  
  return
  
 end subroutine store_vel_init
 
 subroutine compute_vel_init()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the mean velocity module vector
!     of the nanofiber at the nozzle
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  if(ncounterivel/=0)then
    meanivel=counterivel/dble(ncounterivel)
  else
    meanivel=0.d0
  endif
  
  ncounterivel=0
  counterivel=0.d0
  
  return
  
 end subroutine compute_vel_init
 
 subroutine store_vel_end(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for storing the actual velocity module vector
!     of the nanofiber at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  double precision :: tempmod0
  
  if(.not.lrem)return
  
  do ipoint=inpjet-nremovedsub,inpjet-1
  
    select case(systype)
      case(1:2)
        tempmod0=jetvx(ipoint)
      case default
        tempmod0=dsqrt(jetvx(ipoint)**2.d0+ &
         jetvy(ipoint)**2.d0+jetvz(ipoint)**2.d0)
    end select
    
    ncounterevel=ncounterevel+1
    counterevel=counterevel+tempmod0
    
  enddo
  
  return
  
 end subroutine store_vel_end
 
 subroutine compute_vel_end()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the mean velocity module vector
!     of the nanofiber at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  if(ncounterevel/=0)then
    meanevel=counterevel/dble(ncounterevel)
  else
    meanevel=0.d0
  endif
  
  ncounterevel=0
  counterevel=0.d0
  
  return
  
 end subroutine compute_vel_end
 
 subroutine store_velrel_end(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for storing the actual velocity strain value
!     of the nanofiber at the collector
!     NOTE: the velocity strain is defined as the difference in velocity
!     between the two beads closest to the collector divided by their
!     mutual distance
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  double precision :: tempmod0
  
  if(.not.lrem)return
  
  do ipoint=inpjet-nremovedsub,inpjet-1
  
    select case(systype)
      case(1:2)
        continue
        tempmod0=jetvx(ipoint)-jetvx(ipoint+1)
      case default
        tempmod0=dsqrt((jetvx(ipoint)-jetvx(ipoint+1))**2.d0+ &
         (jetvy(ipoint)-jetvy(ipoint+1))**2.d0+(jetvz(ipoint)- &
          jetvz(ipoint+1))**2.d0)
    end select
  
    ncounterevelrel=ncounterevelrel+1
    counterevelrel=counterevelrel+tempmod0
  
  enddo
  
  return
  
 end subroutine store_velrel_end
 
 subroutine compute_velrel_end()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computine the mean velocity strain value
!     of the nanofiber at the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  if(ncounterevelrel/=0)then
    meanevelrel=counterevelrel/dble(ncounterevelrel)
  else
    meanevelrel=0.d0
  endif
  
  ncounterevelrel=0
  counterevelrel=0.d0
  
  return
  
 end subroutine compute_velrel_end
 
 subroutine store_geometry_end(lrem,nremovedsub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for storing the actual mutual distance
!     between the two beads closest to the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  logical,intent(in) :: lrem
  integer,intent(in) :: nremovedsub
  
  integer :: ipoint
  double precision :: tempmod0,tempmod1
  
  if(.not.lrem)return
  
  do ipoint=inpjet-nremovedsub,inpjet-1
    select case(systype)
      case(1:2)
        tempmod0=jetxx(ipoint)-jetxx(ipoint+1)
      case default
        tempmod0=dsqrt((jetxx(ipoint)-jetxx(ipoint+1))**2.d0+ &
         (jetyy(ipoint)-jetyy(ipoint+1))**2.d0+(jetzz(ipoint)- &
          jetzz(ipoint+1))**2.d0)
      end select
  
    tempmod1=dsqrt((icrossec**2.d0)*resolution/tempmod0)
  
    ncountergeom=ncountergeom+1
    counterelen=counterelen+tempmod0
    counterecross=counterecross+tempmod1
    
  enddo
  
  return
  
 end subroutine store_geometry_end
 
 subroutine compute_geometry_end()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computine the mean mutual distance
!     between the two beads closest to the collector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  if(ncountergeom/=0)then
    meanelen=counterelen/dble(ncountergeom)
    meanecross=counterecross/dble(ncountergeom)
  else
    meanelen=0.d0
    meanecross=0.d0
  endif
  
  ncountergeom=0
  counterelen=0.d0
  counterecross=0.d0
  
  return
  
 end subroutine compute_geometry_end
 
 subroutine compute_cpu_time()
 
!***********************************************************************
!     
!     JETSPIN subroutine for registering the actual CPU time
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, save :: timeold=0.d0
  double precision :: timenew
  
  call time_world(timenew)
  meancputime=timenew-timeold
  timeold=timenew
  
  return
  
 end subroutine compute_cpu_time
 
 subroutine compute_elapsed_cpu_time()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the elapsed CPU time
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, save :: initialtime=0.d0
  logical, save :: lfirst=.true.
  
  double precision :: timenew
  
  if(lfirst)then
    lfirst=.false.
    call time_world(initialtime)
  endif
  
  call time_world(timenew)
  elapsedcputime=timenew-initialtime
  
  return
  
 end subroutine compute_elapsed_cpu_time
  
 end module statistic_mod


