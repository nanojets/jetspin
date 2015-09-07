
 module version_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which replace all the 
!     communication tasks providing a serial version of the code
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
 
 contains
 
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

  paralleltype=0
  
  write(iu,'(a50)')'The code is running in serial mode                '
  
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
  
  mxchunk=mxnpjet
  
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
  
  inidom(:)=inpjet
  minidom(:)=inpjet
  
  enddom(:)=npjet
  menddom(:)=npjet
  
  mystart=inpjet
  myend=npjet
  
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
  
  idrank=0
  
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
  
  mxrank=1
  
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
  
  integer ier
  
  return
  
 end subroutine abort_world
 
 subroutine time_world(timecpu)
 
!***********************************************************************
!     
!     JETSPIN subroutine to use the CPU clock
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(out) :: timecpu
  
  call cpu_time(timecpu)
  
  return
  
 end subroutine time_world
 
 subroutine bcast_world_i(buffer)
 
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
  
  integer, intent(inout) :: buffer
  
  return
  
 end subroutine bcast_world_i
 
 subroutine bcast_world_l(buffer)
 
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
  
  logical, intent(inout) :: buffer
  
  return
  
 end subroutine bcast_world_l
 
 subroutine bcast_world_d(buffer)
 
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
  
  double precision, intent(inout) :: buffer
  
  return
  
 end subroutine bcast_world_d
 
 subroutine bcast_world_iarr(buffer,narr)
 
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
  
  integer, intent(inout), dimension(narr) :: buffer
  integer, intent(in) :: narr
  
  return
  
 end subroutine bcast_world_iarr
 
 subroutine bcast_world_larr(buffer,narr)
 
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
  
  logical, intent(inout), dimension(narr) :: buffer
  integer, intent(in) :: narr
  
  return
  
 end subroutine bcast_world_larr
 
 subroutine bcast_world_darr(buffer,narr)
 
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
  
  double precision, intent(inout), dimension(narr) :: buffer
  integer, intent(in) :: narr
  
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
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  double precision, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  double precision, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  integer, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  integer, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  double precision, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  double precision, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
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
  
  logical, intent(inout), dimension(narr) :: argument
  integer, intent(in) :: narr
  logical, intent(inout), dimension(narr), optional :: buffersub
  
  integer ier
  
  if(present(buffersub))buffersub(1:narr)=argument(1:narr)
  
  return
  
 end subroutine or_world_larr
 
 end module version_mod
