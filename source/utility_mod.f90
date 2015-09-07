 
 module utility_mod
 
!***********************************************************************
!     
!     JETSPIN module containing generic supporting subroutines called 
!     by different modules 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
 use version_mod, only : idrank
 
 implicit none
 
 private
  
 double precision, public, parameter :: & 
  Pi=3.141592653589793238462643383279502884d0
 double precision, allocatable,save :: wienerlist(:)
 double precision,save :: hwiener
 integer,save :: winenernodes
 
 public :: init_random_seed,gauss,wiener_process1,wiener_process2,wiener
 public :: modulvec
 public :: dot
 public :: cross
 public :: sig
 
 contains
 
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
    if(isnan(dcos(gauss)))lredo=.true.
  enddo
  
  end function gauss
  
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
 
 end module utility_mod


