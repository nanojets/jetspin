 
module utility_mod

 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
 
!***********************************************************************
!     
!     JETSPIN module containing generic supporting subroutines called 
!     by different modules 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2016
!     
!***********************************************************************
 
 use version_mod, only : idrank,bcast_world_darr
 
 implicit none
 
 private
 
 integer, public, save :: nlbuffservice=0
 integer, public, save :: nibuffservice=0
 integer, public, save :: nbuffservice=0
 logical, public, allocatable, save :: lbuffservice(:)
 integer, public, allocatable, save :: ibuffservice(:)
 double precision, public, allocatable, save :: buffservice(:)
  
 double precision, public, parameter :: & 
  Pi=3.141592653589793238462643383279502884d0
 double precision, allocatable,save :: wienerlist(:)
 double precision, allocatable,save :: gaussianbuffer(:)
 double precision, public, allocatable, save :: gaussianhistory(:)
 integer,save :: ngaussianbuffer=-1
 integer,save :: ngaussianhistory=-1
 integer, public, save :: gaussianhistorysteps=-1
 ! Size of the pre-generated Gaussian pool, in values (input directive
! "noise pool"). The pool is read cyclically, so the noise sequence repeats
! after maxgaussianhistory/(6*active beads) timesteps; a larger pool
! lengthens that period at the cost of 8 bytes per value on the host and,
! for OpenACC builds, again on the device.
 integer, public, save :: maxgaussianhistory=100000000
 integer, public, parameter :: mingaussianhistory=1000000
 integer, public, parameter :: maxgaussianhistorylimit=2000000000
! Sequential-pool consumption state.
!
! The history is one flat random sequence, allocated once and never
! remapped. Each timestep consumes exactly (active beads)*6 consecutive
! values starting at gaussianhistorybase, and gaussianhistorycursor walks
! forward until the whole sequence is used before wrapping to the start.
! (Until 2026-09-30 the default layout indexed a 4D array (step, bead,
! component, draw) with the bead stride hard-wired to mxnpjet+1, so every
! capacity change forced a host rebuild plus a device remap and shortened
! the covered cycle; the pool was then available only with
! JETSPIN_GPU_DYNAMIC_PLATEN.) gaussianhistoryvalues stays <= 0 while no
! history is prepared, which makes begin_gaussian_history_step a no-op.
 integer, public, save :: gaussianhistoryvalues=-1
 integer, public, save :: gaussianhistorybase=0
 integer, public, save :: gaussianhistorywindow=0
 integer, public, save :: gaussianhistoryfirst=0
 integer, save :: gaussianhistorycursor=0
 double precision,save :: hwiener
 integer,save :: winenernodes
! Tracks, independently of any caller's guess, whether gaussianhistory is
! currently present on the accelerator device. The initial mapping happens
! in integrator_mod.f90 (prepare_integrator_random_history), outside this
! module, so that caller reports it here via mark_gaussianhistory_device_
! mapped once done. resize_gaussian_history then trusts this flag instead
! of a value recomputed once per timestep elsewhere, which can be one step
! stale relative to an integrator activating its persistent path in the
! very same step that a capacity-growing refinement event fires.
 logical, save :: gaussianhistory_device_mapped=.false.

 public :: allocate_array_lbuffservice
 public :: allocate_array_ibuffservice
 public :: allocate_array_buffservice
 public :: init_random_seed,gauss,wiener_process1,wiener_process2,wiener
 public :: prepare_gaussian_buffer,gaussian_buffer_value
 public :: prepare_gaussian_history,resize_gaussian_history
 public :: gaussian_history_value
 public :: mark_gaussianhistory_device_mapped
 public :: begin_gaussian_history_step
 public :: modulvec
 public :: dot
 public :: cross
 public :: sig
 public :: write_fmtnumb
 public :: get_prntime
 
 contains
 
 subroutine allocate_array_lbuffservice(imiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array buff 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nlbuffservice/=0)then
    if(imiomax>nlbuffservice)then
      deallocate(lbuffservice)
      nlbuffservice=imiomax+100
      allocate(lbuffservice(0:nlbuffservice))
    endif
  else
    nlbuffservice=imiomax+100
    allocate(lbuffservice(0:nlbuffservice))
  endif
  
  return
  
 end subroutine allocate_array_lbuffservice
 
 subroutine allocate_array_ibuffservice(imiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array buff 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nibuffservice/=0)then
    if(imiomax>nibuffservice)then
      deallocate(ibuffservice)
      nibuffservice=imiomax+100
      allocate(ibuffservice(0:nibuffservice))
    endif
  else
    nibuffservice=imiomax+100
    allocate(ibuffservice(0:nibuffservice))
  endif
  
  return
  
 end subroutine allocate_array_ibuffservice
 
 subroutine allocate_array_buffservice(imiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array buff 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  if(nbuffservice/=0)then
    if(imiomax>nbuffservice)then
      deallocate(buffservice)
      nbuffservice=imiomax+100
      allocate(buffservice(0:nbuffservice))
    endif
  else
    nbuffservice=imiomax+100
    allocate(buffservice(0:nbuffservice))
  endif
  
  return
  
 end subroutine allocate_array_buffservice
 
 subroutine init_random_seed(myseed)
 
!***********************************************************************
!     
!     JETSPIN subroutine for initialising the random generator
!     by the seed given in input file or by a random seed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in),optional :: myseed
  integer :: i, n, clock
  
  integer, allocatable :: seed(:)
          
  call random_seed(size = n)
  
  allocate(seed(n))
  
  if(present(myseed))then
!   If the seed is given in input
    seed = myseed*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
  else
!   If the seed is not given in input it is generated by the clock
    call system_clock(count=clock)
         
    seed = clock*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
  endif
  
  call random_seed(put = seed)
       
  deallocate(seed)
  
  return
 
 end subroutine init_random_seed
 
 function gauss()
 
!***********************************************************************
!     
!     JETSPIN subroutine for generating random number normally
!     distributed by the Box-Muller transformation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision :: gauss
  double precision :: dtemp1,dtemp2
  logical :: lredo
  
  call random_number(dtemp1)
  call random_number(dtemp2)
  
  lredo=.true.
  
! the number is extract again if it is nan
  do while(lredo)
    lredo=.false.
!   Box-Muller transformation
    gauss=dsqrt(-2.d0*dlog(dtemp1))*dcos(2*pi*dtemp2)
    if(ieee_is_nan(dcos(gauss)))lredo=.true.
  enddo
  
 end function gauss

 subroutine prepare_gaussian_buffer(inpnt,npnt,mxpnt,ndim)

!***********************************************************************
!
!     Generate one rank-independent block of Gaussian random numbers.
!     Rank 0 advances the pseudo-random sequence in global bead order and
!     broadcasts the complete block.  Every rank can consequently address
!     a value by global bead index, component, and draw number without the
!     result depending on the MPI domain decomposition.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: inpnt,npnt,mxpnt,ndim
  integer :: ipoint,icomponent,idraw,nvalues

  if(mxpnt<0)stop "Invalid Gaussian-buffer capacity"
  if(ndim<1 .or. ndim>3)stop "Invalid Gaussian-buffer dimension"

  if((.not.allocated(gaussianbuffer)) .or. mxpnt/=ngaussianbuffer)then
    if(allocated(gaussianbuffer))deallocate(gaussianbuffer)
    allocate(gaussianbuffer(0:(mxpnt+1)*3*2-1))
    ngaussianbuffer=mxpnt
  endif

  gaussianbuffer(:)=0.d0
  if(idrank==0)then
    do ipoint=inpnt,npnt
      do icomponent=1,ndim
        do idraw=1,2
          gaussianbuffer(ipoint+(mxpnt+1)*((icomponent-1)+ &
           3*(idraw-1)))=gauss()
        enddo
      enddo
    enddo
  endif

  nvalues=(mxpnt+1)*3*2
  call bcast_world_darr(gaussianbuffer,nvalues)

  return

 end subroutine prepare_gaussian_buffer

 function gaussian_buffer_value(ipoint,icomponent,idraw)

  implicit none

  integer, intent(in) :: ipoint,icomponent,idraw
  double precision :: gaussian_buffer_value

  if(.not.allocated(gaussianbuffer))stop "Gaussian buffer is not prepared"
  if(ipoint<0 .or. ipoint>ngaussianbuffer)stop "Invalid Gaussian bead index"
  if(icomponent<1 .or. icomponent>3)stop "Invalid Gaussian component"
  if(idraw<1 .or. idraw>2)stop "Invalid Gaussian draw index"

  gaussian_buffer_value=gaussianbuffer(ipoint+(ngaussianbuffer+1)* &
   ((icomponent-1)+3*(idraw-1)))

  return

 end function gaussian_buffer_value

 subroutine prepare_gaussian_history(inpnt,npnt,mxpnt,ndim,nsteps,full_pool)
  implicit none
  integer, intent(in) :: inpnt,npnt,mxpnt,ndim,nsteps
  logical, intent(in), optional :: full_pool
  integer :: nperstep,nvalues,index,nwindow,istep,ipoint,icomponent,idraw
  logical :: whole_pool

  if(mxpnt<0 .or. nsteps<1)stop "Invalid Gaussian-history extent"
  if(ndim<1 .or. ndim>3)stop "Invalid Gaussian-history dimension"
  if(npnt<inpnt)stop "Invalid Gaussian-history bead range"
  whole_pool=.false.
  if(present(full_pool))whole_pool=full_pool
! One flat random sequence, allocated once: capacity growth later in the
! run needs no rebuild here and no device remap (see
! resize_gaussian_history).  Each timestep reads (active beads)*6
! consecutive values (begin_gaussian_history_step).
  if(whole_pool)then
! A run whose bead count can grow (insertion, dynamic refinement) takes the
! whole pool at once, filled in index order: the per-step window changes
! during the run, so no estimate from the initial beads would fit.
    nvalues=maxgaussianhistory
    nperstep=(npnt-inpnt+1)*3*2
    if(nvalues<nperstep)stop "One Gaussian timestep exceeds history limit"
    gaussianhistorysteps=nvalues/nperstep
  else
! A fixed-topology run reads the same window npnt-inpnt+1 at every step.
! Fill each step's slice in the draw order of prepare_gaussian_buffer
! (bead, component, draw), so that a serial run reading the pool uses
! exactly the numbers that an MPI or non-history run draws step by step.
! Divide before multiplying: nsteps*nperstep would overflow a default
! integer for a long run, while maxgaussianhistory/nperstep cannot.
    nperstep=(npnt-inpnt+1)*3*2
    gaussianhistorysteps=min(nsteps,maxgaussianhistory/nperstep)
    if(gaussianhistorysteps<1)stop "One Gaussian timestep exceeds history limit"
    nvalues=gaussianhistorysteps*nperstep
  endif
  gaussianhistoryvalues=nvalues
  gaussianhistorycursor=0
  gaussianhistorybase=0
  gaussianhistorywindow=0
  gaussianhistoryfirst=inpnt
  if(allocated(gaussianhistory))deallocate(gaussianhistory)
  allocate(gaussianhistory(0:nvalues-1))
  gaussianhistory(:)=0.d0
  if(idrank==0)then
    if(whole_pool)then
      do index=0,nvalues-1
        gaussianhistory(index)=gauss()
      enddo
    else
      nwindow=npnt-inpnt+1
      do istep=0,gaussianhistorysteps-1
        do ipoint=inpnt,npnt
          do icomponent=1,ndim
            do idraw=1,2
              index=istep*nperstep+(ipoint-inpnt)+nwindow* &
               ((icomponent-1)+3*(idraw-1))
              gaussianhistory(index)=gauss()
            enddo
          enddo
        enddo
      enddo
    endif
  endif
  call bcast_world_darr(gaussianhistory,nvalues)
  ngaussianhistory=mxpnt
! The step count is exact for a fixed window; for a dynamic run it refers
! to npnt-inpnt+1 active beads, and the period scales inversely with the
! active count (docs/introduction/random-numbers.md).
  if(idrank==0)then
    write(6,'(a,i0,a,i0,a,i0,a)')'Gaussian history pool: values=',nvalues, &
     ' covers ',gaussianhistorysteps,' steps at ',npnt-inpnt+1,' beads'
  endif
 end subroutine prepare_gaussian_history

 subroutine begin_gaussian_history_step(firstpoint,lastpoint)

!***********************************************************************
!
!     JETSPIN subroutine for reserving this timestep's slice of the
!     sequential Gaussian pool and advancing the cursor for the next one.
!     Must be called exactly once per timestep, before any read, on every
!     execution path of the integrator that consumes the history.
!
!     A no-op unless the pool layout is in use (gaussianhistoryvalues<=0
!     otherwise), so call sites need no conditional compilation.
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!
!***********************************************************************

  implicit none

  integer, intent(in) :: firstpoint,lastpoint

  if(gaussianhistoryvalues<=0)return
  gaussianhistorybase=gaussianhistorycursor
  gaussianhistoryfirst=firstpoint
  gaussianhistorywindow=max(1,lastpoint-firstpoint+1)
! Wrap only once the whole stored sequence has been consumed. The slice is
! read with the same modulo, so a slice straddling the end of the pool is
! served correctly rather than by discarding the tail.
  gaussianhistorycursor=mod(gaussianhistorycursor+gaussianhistorywindow*6, &
   gaussianhistoryvalues)

 end subroutine begin_gaussian_history_step

 subroutine mark_gaussianhistory_device_mapped(mapped)
  implicit none
  logical, intent(in) :: mapped
  gaussianhistory_device_mapped=mapped
 end subroutine mark_gaussianhistory_device_mapped

 subroutine resize_gaussian_history(mxpnt,device_mapped)
! The flat pool is independent of the bead count, so a capacity increase
! needs no rebuild, no move_alloc and no device remap. Keeping this entry
! point as a no-op leaves the callers (capacity growth in reallocate_jet
! and in dynamic refinement) unchanged and structurally removes the
! stale-device-mapping failure mode that the former rebuild path had.
  implicit none
  integer, intent(in) :: mxpnt
  logical, intent(in), optional :: device_mapped
  return
 end subroutine resize_gaussian_history

 function gaussian_history_value(istep,ipoint,icomponent,idraw)
  implicit none
  integer, intent(in) :: istep,ipoint,icomponent,idraw
  integer :: index
  double precision :: gaussian_history_value
  if(.not.allocated(gaussianhistory))stop "Gaussian history is not prepared"
  if(istep<1)stop "Invalid Gaussian-history step"
  if(icomponent<1 .or. icomponent>3)stop "Invalid Gaussian component"
  if(idraw<1 .or. idraw>2)stop "Invalid Gaussian draw index"
! Read inside this timestep's reserved slice. The valid bead range is the
! live window set by begin_gaussian_history_step, not a fixed capacity:
! the pool never grows, so bounds must follow the beads, not the array.
! istep is unused here; the position comes from the cursor instead, which
! is what lets the sequence be consumed end to end.
  if(ipoint<gaussianhistoryfirst .or. &
   ipoint>gaussianhistoryfirst+gaussianhistorywindow-1) &
   stop "Invalid Gaussian bead index"
  index=mod(gaussianhistorybase+(ipoint-gaussianhistoryfirst)+ &
   gaussianhistorywindow*((icomponent-1)+3*(idraw-1)),gaussianhistoryvalues)
  gaussian_history_value=gaussianhistory(index)
 end function gaussian_history_value
  
  subroutine wiener_process1(inpnt,npnt,nvar,ndim,h,fwienersub1)
  
!***********************************************************************
!     
!     JETSPIN subroutine for generating and storing a single
!     wiener process
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: inpnt,npnt,nvar,ndim
  double precision, intent(in) :: h
  double precision,intent(inout),allocatable :: fwienersub1(:,:,:)
  double precision :: hw1,psi
  integer :: i,j
  
  fwienersub1(:,:,:)=0.d0
  
  hw1=dsqrt(h)
  
  call flush(6)
  do i=inpnt,npnt
    do j=1,ndim
      psi=gauss()
      fwienersub1(i,3,j)=hw1*psi
    enddo
  enddo
    
  return
  
 end subroutine wiener_process1
 
 subroutine wiener_process2(inpnt,npnt,nvar,ndim,h,fwienersub1, &
  fwienersub2)
 
!***********************************************************************
!     
!     JETSPIN subroutine for generating and storing two different
!     wiener processes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: inpnt,npnt,nvar,ndim
  double precision, intent(in) :: h
  double precision,intent(inout),allocatable :: fwienersub1(:,:,:), &
   fwienersub2(:,:,:)
  double precision :: hw1,hw2,co1,co2,psi,theta
  integer :: i,j
  
  fwienersub1(:,:,:)=0.d0
  fwienersub2(:,:,:)=0.d0
  
  hw1=dsqrt(h)
  hw2=(hw1)**3.d0
  co1=0.5d0
  co2=1.d0/(2.d0*dsqrt(3.d0))
  
  do i=inpnt,npnt
    do j=1,ndim
      psi=gauss()
      theta=gauss()
      fwienersub1(i,3,j)=hw1*psi
      fwienersub2(i,3,j)=hw2*(co1*psi+co2*theta)
    enddo
  enddo
    
  return
  
 end subroutine wiener_process2
 
 function wiener(t)
 
!***********************************************************************
!     
!     JETSPIN function for extracting the index of a wiener process
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision,intent(in) :: t
  integer :: ind
  double precision :: wiener
  
  ind=nint((t-wienerlist(0))/hwiener)
  
  wiener=wienerlist(ind)
  
  
  return
  
 end function wiener
 
 pure function modulvec(a)
 
!***********************************************************************
!     
!     JETSPIN function for computing the module of a vector
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision :: modulvec
  double precision, dimension(3), intent(in) :: a

  modulvec = dsqrt(a(1)**2.d0 + a(2)**2.d0+ a(3)**2.d0)
  
  return
  
 end function modulvec
 
 pure function cross(a,b)
 
!***********************************************************************
!     
!     JETSPIN function for computing the cross product of two vectors
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, dimension(3) :: cross
  double precision, dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
  
  return
  
 end function cross
 
 pure function dot(a,b)
 
!***********************************************************************
!     
!     JETSPIN function for computing the dot product of two vectors
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision :: dot
  double precision, dimension(3), intent(in) :: a, b
  
  dot=dot_product(a,b)
  
  return
  
 end function dot
 
 pure function sig(num)
 
!***********************************************************************
!     
!     JETSPIN function for returning the sign of a floating number
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in):: num
  
  double precision :: sig
  
  if(num>0)then
    sig=1.d0
  elseif(num==0)then
    sig=0.d0
  else
    sig=-1.d0
  endif
  
  return
 
 end function sig
 
 function dimenumb(inum)
 
!***********************************************************************
!     
!     JETSPIN function for returning the number of digits
!     of an integer number
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************

  implicit none

  integer,intent(in) :: inum
  integer :: dimenumb
  integer :: i
  double precision :: tmp

  i=1
  tmp=dble(inum)
  do 
    if(tmp<10.d0)exit
    i=i+1
    tmp=tmp/10.d0
  enddo

  dimenumb=i

  return

 end function dimenumb

 function write_fmtnumb(inum)
 
!***********************************************************************
!     
!     JETSPIN function for returning the string of six characters 
!     with integer digits and leading zeros to the left
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: write_fmtnumb
  integer :: numdigit,irest
  double precision :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=6-numdigit
  write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
  write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum


  return

 end function write_fmtnumb
 
 subroutine get_prntime(hms,timelp,prntim)
  
!***********************************************************************
!     
!     JETSPIN subroutine for casting cpu elapsed time into days, hours,
!     minutes and seconds for printing (input timelp in seconds)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2016
!     
!***********************************************************************
  
  implicit none
  
  character(len=1), intent(out) :: hms
  double precision, intent(in) :: timelp
  double precision, intent(out) :: prntim
  
  if(timelp.ge.8.64d4)then
    hms='d'
    prntim=timelp/8.64d4
  elseif(timelp.ge.3.6d3)then
    hms='h'
    prntim=timelp/3.6d3
  elseif(timelp.ge.6.0d1)then
    hms='m'
    prntim=timelp/6.0d1
  else
    hms='s'
    prntim=timelp
  endif
  
  return
  
 end subroutine get_prntime
 
 end module utility_mod
