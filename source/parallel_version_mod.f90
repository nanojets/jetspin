
 module version_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which manage all the 
!     communication tasks for the code parallel version and assign
!     a subset of beads to a specific node
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
 implicit none
 
 private
 
 integer, public, save :: idrank=0
 integer, public, save :: mxrank=1
 integer, public, save :: nchunkmin=10
 
 integer, public, save :: paralleltype
 integer, public, save :: mxchunk
 integer, public, save :: mystart
 integer, public, save :: myend
 integer, public, save, allocatable, dimension(:) :: inidom,enddom,minidom,menddom
 
 public :: print_version
 public :: alloc_domain
 public :: set_mxchunk
 public :: set_chunk
 public :: get_rank_world
 public :: get_size_world
 public :: get_sync_world
 public :: init_world
 public :: finalize_world
 public :: abort_world
 public :: time_world
 public :: bcast_world_i
 public :: bcast_world_l
 public :: bcast_world_d
 public :: bcast_world_iarr
 public :: bcast_world_larr
 public :: bcast_world_darr
 public :: sum_world_iarr
 public :: sum_world_darr
 public :: min_world_iarr
 public :: min_world_darr
 public :: max_world_iarr
 public :: max_world_darr
 public :: and_world_larr
 public :: or_world_larr
 
  integer, save, allocatable, dimension(:) :: nrank
 integer, allocatable, dimension(:), save :: ibuffer
 double precision, allocatable, dimension(:), save :: dbuffer
 logical, allocatable, dimension(:), save :: lbuffer
 integer, save :: nibuffer,ndbuffer,nlbuffer
 logical, save :: aibuffer=.false.
 logical, save :: adbuffer=.false.
 logical, save :: albuffer=.false.
 integer, parameter :: nincrement=100
 integer, parameter :: minchunk=10
 
 contains
 
 subroutine allocate_ibuffer(narr)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the buffer array of integer
!     type
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(aibuffer)then
    if(narr>nibuffer)then
      deallocate(ibuffer)
      nibuffer=narr+nincrement
      allocate(ibuffer(nibuffer))
    endif
  else
    nibuffer=narr+nincrement
    allocate(ibuffer(nibuffer))
    aibuffer=.true.
  endif
  
  return
  
 end subroutine allocate_ibuffer
 
 subroutine allocate_dbuffer(narr)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the buffer array of double
!     precision type
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(adbuffer)then
    if(narr>ndbuffer)then
      deallocate(dbuffer)
      ndbuffer=narr+nincrement
      allocate(dbuffer(ndbuffer))
    endif
  else
    ndbuffer=narr+nincrement
    allocate(dbuffer(ndbuffer))
    adbuffer=.true.
  endif
  
  return
  
 end subroutine allocate_dbuffer
 
 subroutine allocate_lbuffer(narr)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the buffer array of logical
!     type
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: narr
  
  if(albuffer)then
    if(narr>nlbuffer)then
      deallocate(lbuffer)
      nlbuffer=narr+nincrement
      allocate(lbuffer(nlbuffer))
    endif
  else
    nlbuffer=narr+nincrement
    allocate(lbuffer(nlbuffer))
    albuffer=.true.
  endif
  
  return
  
 end subroutine allocate_lbuffer
 
 subroutine print_version(iu)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing the code version which is
!     currently in use
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=50), intent(out) :: iu
  
  paralleltype=2
  

  if(mxrank==1)then
    write(iu,'(a40,i4,a6)')'The code is running in parallel mode on ',mxrank, &
       ' CPUs '
  else
    write(iu,'(a40,i4,a6)')'The code is running in parallel mode on ',mxrank, &
       ' CPUs '
  endif
  
  return
  
 end subroutine print_version
 
 subroutine alloc_domain()
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating service bookkeeping arrays
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  allocate(inidom(0:mxrank-1))
  allocate(enddom(0:mxrank-1))
  allocate(minidom(0:mxrank-1))
  allocate(menddom(0:mxrank-1))
  allocate(nrank(0:mxrank-1))
  
  return
  
 end subroutine alloc_domain
 
 subroutine set_mxchunk(mxnpjet)
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining the number of beads which are 
!     dealt from a specific node
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: mxnpjet
  
  integer :: tchunk
  
  tchunk=ceiling(dble(mxnpjet)/dble(mxrank))+1
  
  mxchunk=max(nchunkmin,tchunk)
  
  return
  
 end subroutine set_mxchunk
 
 subroutine set_chunk(inpjet,npjet)
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining and assigning a subset of beads
!     to a specific node
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: inpjet
  integer, intent(in) :: npjet
  
  integer :: nattemp,nchunktemp,tchunk,remain,i
  
  nattemp=npjet-inpjet
 
  nchunktemp=nchunkmin*mxrank
  tchunk=floor(dble(nattemp+1)/dble(mxrank))
  
  if(tchunk<nchunkmin)then
    i=0
    nrank(:)=0
    remain=nattemp+1
    do while (remain>0)
      if(nrank(i)<nchunkmin)then
        nrank(i)=nrank(i)+1
        remain=remain-1
      else
        i=i+1
        if(i>=mxrank)exit
      endif
    enddo
  else
    nrank(:)=tchunk
    remain=nattemp+1-(tchunk*mxrank)
    do i=0,mxrank-1
      remain=remain-1
      if(remain<0)exit
      nrank(i)=nrank(i)+1
    enddo
  endif
  
  inidom(:)=0
  enddom(:)=-1
  
  minidom(:)=0
  menddom(:)=-1
  
  inidom(0)=inpjet
  enddom(0)=nrank(0)-1+inidom(0)
 
  do i=1,mxrank-1
    if(nrank(i)==0)exit
    inidom(i)=enddom(i-1)+1
    enddom(i)=inidom(i)+nrank(i)-1
  enddo
  
  minidom(0)=inpjet
  if(enddom(0)+1>npjet)then
    menddom(0)=npjet
  else
    menddom(0)=enddom(0)+1
    do i=1,mxrank-1
      minidom(i)=inidom(i)-1
      if(enddom(i)+1>npjet)then
        menddom(i)=npjet
        exit
      else
        menddom(i)=enddom(i)+1
      endif
    enddo
  endif
  
  mystart=inidom(idrank)
  myend=enddom(idrank)
  
  
  return
  
 end subroutine set_chunk
 
 subroutine get_rank_world()
 
!***********************************************************************
!     
!     JETSPIN subroutine to determine identity of processing node 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  include 'mpif.h'
  
  integer ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,idrank,ier)
  
  return
  
 end subroutine get_rank_world
 
 subroutine get_size_world()
 
!***********************************************************************
!     
!     JETSPIN subroutine to determine the number of processing nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  include 'mpif.h'
  
  integer ier
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mxrank,ier)
  
  return
  
 end subroutine get_size_world
 
 subroutine get_sync_world()
 
!***********************************************************************
!     
!     JETSPIN subroutine to synchronize the processing nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  include 'mpif.h'
  
  integer ier
  
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine get_sync_world
 
 subroutine init_world()
 
!***********************************************************************
!     
!     JETSPIN subroutine to initialize the MPI work
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************

  implicit none
  
  include 'mpif.h'
  
  integer ier
  
  call MPI_INIT(ier)
  
  return
 
 end subroutine init_world
 
 subroutine finalize_world()
 
!***********************************************************************
!     
!     JETSPIN subroutine to finalize the MPI work
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  integer ier
  
  call MPI_FINALIZE(ier)
  
  return
  
 end subroutine finalize_world
 
 subroutine abort_world()
 
!***********************************************************************
!     
!     JETSPIN subroutine to abort the MPI work
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  include 'mpif.h'
  
  integer ier
 
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  
  return
  
 end subroutine abort_world
 
 subroutine time_world(timecpu)
 
!***********************************************************************
!     
!     JETSPIN subroutine to use the MPI CPU clock
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  include 'mpif.h'
  
  double precision, intent(out) :: timecpu
 
  timecpu = MPI_WTIME()
  
  return
  
 end subroutine time_world
 
 subroutine bcast_world_i(argument)
 
!***********************************************************************
!     
!     JETSPIN subroutine to broadcast an integer number to all other 
!     nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  integer, intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_i
 
 subroutine bcast_world_l(argument)
 
!***********************************************************************
!     
!     JETSPIN subroutine to broadcast a logical variable to all other 
!     nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  logical, intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_l
 
 subroutine bcast_world_d(argument)
 
!***********************************************************************
!     
!     JETSPIN subroutine to broadcast an double precision number to all 
!     other nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  double precision, intent(inout) :: argument
  
  integer ier
  
  call MPI_BCAST(argument,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_d
 
 subroutine bcast_world_iarr(argument,narr)
 
!***********************************************************************
!     
!     JETSPIN subroutine to broadcast an integer array to all other 
!     nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_iarr
 
 subroutine bcast_world_larr(argument,narr)
 
!***********************************************************************
!     
!     JETSPIN subroutine to broadcast a logical array to all other 
!     nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_larr
 
 subroutine bcast_world_darr(argument,narr)
  
!***********************************************************************
!     
!     JETSPIN subroutine to broadcast an double precision array to all 
!     other nodes
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  
  integer ier
  
  call MPI_BCAST(argument,narr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  
  return
  
 end subroutine bcast_world_darr
 
 subroutine sum_world_iarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global summation subroutine for a integer array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_INTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
  return
  
 end subroutine sum_world_iarr
 
 subroutine sum_world_darr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global summation subroutine for a double precision array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  double precision, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_DOUBLE_PRECISION, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_dbuffer(narr)
    
    call MPI_ALLREDUCE(argument,dbuffer,narr,MPI_DOUBLE_PRECISION, &
      MPI_SUM,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=dbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine sum_world_darr
 
 subroutine min_world_iarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global minimum subroutine for an integer array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_INTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
      MPI_MIN,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=ibuffer(1:narr)
    
  endif
  
  return
  
 end subroutine min_world_iarr
 
 subroutine min_world_darr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global minimum subroutine for a double precision array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  double precision, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
  
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_DOUBLE_PRECISION, &
      MPI_MIN,MPI_COMM_WORLD,ier)
  
  else
  
    call allocate_dbuffer(narr)
  
    call MPI_ALLREDUCE(argument,dbuffer,narr,MPI_DOUBLE_PRECISION, &
      MPI_MIN,MPI_COMM_WORLD,ier)
    
    argument(1:narr)=dbuffer(1:narr)
  
  endif
  
  return
  
 end subroutine min_world_darr
 
 subroutine max_world_iarr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global maximum subroutine for an integer array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_INTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_ibuffer(narr)
    
    call MPI_ALLREDUCE(argument,ibuffer,narr,MPI_INTEGER, &
      MPI_MAX,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=ibuffer(1:narr)
  
  endif
  
  return
  
 end subroutine max_world_iarr
 
 subroutine max_world_darr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global maximum subroutine for a double precision array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  double precision, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_DOUBLE_PRECISION, &
      MPI_MAX,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_dbuffer(narr)
    
    call MPI_ALLREDUCE(argument,dbuffer,narr,MPI_DOUBLE_PRECISION, &
      MPI_MAX,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=dbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine max_world_darr
 
 subroutine and_world_larr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global 'logical and' subroutine for a logical array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_LOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MPI_LOGICAL, &
      MPI_LAND,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=lbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine and_world_larr
 
 subroutine or_world_larr(argument,narr,buffersub)
 
!***********************************************************************
!     
!     JETSPIN global 'logical or' subroutine for a logical array
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  include 'mpif.h'
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))then
    
    call MPI_ALLREDUCE(argument,buffersub,narr,MPI_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
    
  else
    
    call allocate_lbuffer(narr)
    
    call MPI_ALLREDUCE(argument,lbuffer,narr,MPI_LOGICAL, &
      MPI_LOR,MPI_COMM_WORLD,ier)
      
    argument(1:narr)=lbuffer(1:narr)
    
  endif
  
  return
  
 end subroutine or_world_larr
  
 end module version_mod
