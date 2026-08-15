 
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
 integer, public, parameter :: maxgaussianhistory=100000000
 double precision,save :: hwiener
 integer,save :: winenernodes
 
 public :: allocate_array_lbuffservice
 public :: allocate_array_ibuffservice
 public :: allocate_array_buffservice
 public :: init_random_seed,gauss,wiener_process1,wiener_process2,wiener
 public :: prepare_gaussian_buffer,gaussian_buffer_value
 public :: prepare_gaussian_history,resize_gaussian_history
 public :: gaussian_history_value
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

 subroutine prepare_gaussian_history(inpnt,npnt,mxpnt,ndim,nsteps)
  implicit none
  integer, intent(in) :: inpnt,npnt,mxpnt,ndim,nsteps
  integer :: istep,ipoint,icomponent,idraw,nperstep,nvalues,index

  if(mxpnt<0 .or. nsteps<1)stop "Invalid Gaussian-history extent"
  if(ndim<1 .or. ndim>3)stop "Invalid Gaussian-history dimension"
  nperstep=(mxpnt+1)*3*2
  gaussianhistorysteps=min(nsteps,maxgaussianhistory/nperstep)
  if(gaussianhistorysteps<1)stop "One Gaussian timestep exceeds history limit"
  nvalues=nperstep*gaussianhistorysteps
  if(allocated(gaussianhistory))deallocate(gaussianhistory)
  allocate(gaussianhistory(0:nvalues-1))
  gaussianhistory(:)=0.d0
  if(idrank==0)then
    do istep=1,gaussianhistorysteps
      do ipoint=inpnt,npnt
        do icomponent=1,ndim
          do idraw=1,2
            index=(istep-1)*nperstep+ipoint+(mxpnt+1)* &
             ((icomponent-1)+3*(idraw-1))
            gaussianhistory(index)=gauss()
          enddo
        enddo
      enddo
    enddo
  endif
  call bcast_world_darr(gaussianhistory,nvalues)
  ngaussianhistory=mxpnt
 end subroutine prepare_gaussian_history

 subroutine resize_gaussian_history(mxpnt,device_mapped)
  implicit none
  integer, intent(in) :: mxpnt
  logical, intent(in), optional :: device_mapped
  integer :: oldcapacity,oldsteps,newsteps
  integer :: oldnperstep,newnperstep,oldnvalues,newnvalues
  integer :: istep,ipoint,icomponent,idraw,oldindex,newindex
  logical :: mapped
  double precision, allocatable :: resizedhistory(:)

  if(.not.allocated(gaussianhistory))return
  if(mxpnt<=ngaussianhistory)return

  mapped=.false.
  if(present(device_mapped))mapped=device_mapped
  oldcapacity=ngaussianhistory
  oldsteps=gaussianhistorysteps
  oldnperstep=(oldcapacity+1)*3*2
  newnperstep=(mxpnt+1)*3*2
  oldnvalues=oldnperstep*oldsteps
  newsteps=min(oldsteps,maxgaussianhistory/newnperstep)
  if(newsteps<1)stop "One resized Gaussian timestep exceeds history limit"
  newnvalues=newnperstep*newsteps
  allocate(resizedhistory(0:newnvalues-1))
  resizedhistory(:)=0.d0

! Preserve every existing bead/component/draw value for the retained cycle.
! The storage stride changes with capacity, so this must be a semantic copy
! rather than a contiguous prefix copy.
  do istep=1,newsteps
    do ipoint=0,oldcapacity
      do icomponent=1,3
        do idraw=1,2
          oldindex=(istep-1)*oldnperstep+ipoint+(oldcapacity+1)* &
           ((icomponent-1)+3*(idraw-1))
          newindex=(istep-1)*newnperstep+ipoint+(mxpnt+1)* &
           ((icomponent-1)+3*(idraw-1))
          resizedhistory(newindex)=gaussianhistory(oldindex)
        enddo
      enddo
    enddo
  enddo

! New capacity slots receive a single rank-independent extension of the
! initial Gaussian history. No random extraction is introduced in the time
! integration loop.
  if(idrank==0)then
    do istep=1,newsteps
      do ipoint=oldcapacity+1,mxpnt
        do icomponent=1,3
          do idraw=1,2
            newindex=(istep-1)*newnperstep+ipoint+(mxpnt+1)* &
             ((icomponent-1)+3*(idraw-1))
            resizedhistory(newindex)=gauss()
          enddo
        enddo
      enddo
    enddo
  endif
  call bcast_world_darr(resizedhistory,newnvalues)

#ifdef _OPENACC
! The old allocation must be detached before move_alloc changes its host
! address. Re-enter the resized history once; ordinary timesteps remain free
! of random-history transfers.
!$acc exit data delete(gaussianhistory(0:oldnvalues-1)) if(mapped)
#endif
  call move_alloc(resizedhistory,gaussianhistory)
  ngaussianhistory=mxpnt
  gaussianhistorysteps=newsteps
#ifdef _OPENACC
!$acc enter data copyin(gaussianhistory(0:newnvalues-1)) if(mapped)
#endif

  if(idrank==0)then
    write(6,'(a,i0,a,i0,a,i0,a,i0)') &
     'Gaussian history capacity: old=',oldcapacity,' new=',mxpnt, &
     ' retained_steps=',newsteps,' values=',newnvalues
  endif
 end subroutine resize_gaussian_history

 function gaussian_history_value(istep,ipoint,icomponent,idraw)
  implicit none
  integer, intent(in) :: istep,ipoint,icomponent,idraw
  integer :: nperstep,index,cycle_step
  double precision :: gaussian_history_value
  if(.not.allocated(gaussianhistory))stop "Gaussian history is not prepared"
  if(istep<1)stop "Invalid Gaussian-history step"
  if(ipoint<0 .or. ipoint>ngaussianhistory)stop "Invalid Gaussian bead index"
  if(icomponent<1 .or. icomponent>3)stop "Invalid Gaussian component"
  if(idraw<1 .or. idraw>2)stop "Invalid Gaussian draw index"
  nperstep=(ngaussianhistory+1)*3*2
  cycle_step=mod(istep-1,gaussianhistorysteps)+1
  index=(cycle_step-1)*nperstep+ipoint+(ngaussianhistory+1)* &
   ((icomponent-1)+3*(idraw-1))
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
