
 module fit_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which are defining
!     and dealing interpolation fit of the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!*********************************************************************** 
 
 use error_mod,   only : error
 use version_mod ,only : idrank,mxrank,sum_world_darr,mystart,myend
 use nanojet_mod, only : mxnpjet,npjet,inpjet,typemass,ltagbeads, &
                          lbreakup
  
 
 implicit none
 
 private
 
 logical, allocatable, dimension(:), save, public :: jetbdc
 logical, allocatable, dimension(:), save, public :: jetbrc
 
 double precision, allocatable, dimension(:), save, public :: jetptc
 double precision, allocatable, dimension(:), save, public :: jetxxc
 double precision, allocatable, dimension(:), save, public :: jetyyc
 double precision, allocatable, dimension(:), save, public :: jetzzc
 double precision, allocatable, dimension(:), save, public :: jetvxc
 double precision, allocatable, dimension(:), save, public :: jetvyc
 double precision, allocatable, dimension(:), save, public :: jetvzc
 double precision, allocatable, dimension(:), save, public :: jetstc
 double precision, allocatable, dimension(:), save, public :: jetmsc
 double precision, allocatable, dimension(:), save, public :: jetchc
 double precision, allocatable, dimension(:), save, public :: jetcrc
 double precision, allocatable, dimension(:), save :: jetak1
 double precision, allocatable, dimension(:), save :: jetak2
 double precision, allocatable, dimension(:), save :: jetak3
 double precision, allocatable, dimension(:), save :: jetak4
 double precision, allocatable, dimension(:), save :: mak
 double precision, allocatable, dimension(:), save :: tak

 
 logical, save :: larrayspline=.false.
 integer, save :: narrayspline=0
 
 logical, save :: larrayjetptc=.false.
 integer, save :: narrayjetptc=0
 
 logical, save :: larrayakima=.false.
 integer, save :: narrayakima=0
 
 integer, save :: inpjetspline
 integer, save :: npjetspline
 
 public :: create_spline
 public :: fit_spline
 public :: allocate_arrayspline
 public :: allocate_array_jetptc
 public :: driver_fit_spline
 public :: looking_indexes_2
 public :: driver_fit_curve2
 public :: looking_indexes_4
 public :: cubic_interpolation
 public :: fit_akima
 public :: allocate_arrayakima
 public :: findcurve4
 
 contains
 
 
 subroutine allocate_arrayspline()

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating arrays which are used to store
!     the spline coefficients if the spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  
  if(larrayspline)then
    if(mxnpjet>narrayspline)then
      deallocate(jetxxc,jetyyc,jetzzc)
      deallocate(jetvxc,jetvyc,jetvzc)
      deallocate(jetstc,jetmsc,jetchc,jetcrc)
      narrayspline=mxnpjet
      allocate(jetxxc(0:narrayspline), &
       jetyyc(0:narrayspline),jetzzc(0:narrayspline))
      allocate(jetvxc(0:narrayspline), &
       jetvyc(0:narrayspline),jetvzc(0:narrayspline))
      allocate(jetstc(0:narrayspline), &
       jetmsc(0:narrayspline),jetchc(0:narrayspline), &
       jetcrc(0:narrayspline))
    endif
  else
    narrayspline=mxnpjet
    allocate(jetxxc(0:narrayspline), &
     jetyyc(0:narrayspline),jetzzc(0:narrayspline))
    allocate(jetvxc(0:narrayspline), &
     jetvyc(0:narrayspline),jetvzc(0:narrayspline))
    allocate(jetstc(0:narrayspline), &
     jetmsc(0:narrayspline),jetchc(0:narrayspline), &
     jetcrc(0:narrayspline))
    larrayspline=.true.
  endif
  
  return
  
 end subroutine allocate_arrayspline
 
 subroutine allocate_array_jetptc()

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the array jetptc which 
!     parametrizes the nanojet if the spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  integer :: ierr=0
  
  if(larrayjetptc)then
    if(mxnpjet>narrayjetptc)then
      deallocate(jetptc)
      narrayjetptc=mxnpjet
      allocate(jetptc(0:narrayjetptc),stat=ierr)
      if(typemass==3  .or. ltagbeads)then
        deallocate(jetbdc)
        allocate(jetbdc(0:narrayjetptc),stat=ierr)
        jetbdc(0:narrayjetptc)=.false.
        if(lbreakup)then
          deallocate(jetbrc)
          allocate(jetbrc(0:narrayjetptc),stat=ierr)
          jetbrc(0:narrayjetptc)=.false.
        endif
      endif
    endif
  else
    narrayjetptc=mxnpjet
    allocate(jetptc(0:narrayjetptc),stat=ierr)
    larrayjetptc=.true.
    if(typemass==3  .or. ltagbeads)then
      allocate(jetbdc(0:narrayjetptc),stat=ierr)
      jetbdc(0:narrayjetptc)=.false.
      if(lbreakup)then
        allocate(jetbrc(0:narrayjetptc),stat=ierr)
        jetbrc(0:narrayjetptc)=.false.
      endif
    endif
  endif
  
  if(ierr/=0)call error(13)
  
  return
  
 end subroutine allocate_array_jetptc
 
 subroutine create_spline(jpt,jxx,jyy,jzz,jvx,jvy,jvz,jst,jms,jch,jcr)

!***********************************************************************
!     
!     JETSPIN subroutine for preparing the cubic 
!     spline interpolation of the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(in) :: jpt
  double precision, allocatable, dimension (:), intent(in) :: jxx
  double precision, allocatable, dimension (:), intent(in) :: jyy
  double precision, allocatable, dimension (:), intent(in) :: jzz
  double precision, allocatable, dimension (:), intent(in) :: jvx
  double precision, allocatable, dimension (:), intent(in) :: jvy
  double precision, allocatable, dimension (:), intent(in) :: jvz
  double precision, allocatable, dimension (:), intent(in) :: jst
  double precision, allocatable, dimension (:), intent(in) :: jms
  double precision, allocatable, dimension (:), intent(in) :: jch
  double precision, allocatable, dimension (:), intent(in) :: jcr
  
  double precision :: derinit,derend
  integer :: i,j
  
  
  jetxxc(:)=0.d0
  jetyyc(:)=0.d0
  jetzzc(:)=0.d0
  jetvxc(:)=0.d0
  jetvyc(:)=0.d0
  jetvzc(:)=0.d0
  jetstc(:)=0.d0
  jetmsc(:)=0.d0
  jetchc(:)=0.d0
  jetcrc(:)=0.d0
  
  inpjetspline=inpjet
  npjetspline=npjet
  
  do i=idrank+1,10,mxrank
  
    select case(i)
      case(1)
        call define_derivatives_spline(jpt,jxx,derinit,derend)
        call spline(jpt,jxx,derinit,derend,jetxxc)
      case(2)
        call define_derivatives_spline(jpt,jyy,derinit,derend)
        call spline(jpt,jyy,derinit,derend,jetyyc)
      case(3)
        call define_derivatives_spline(jpt,jzz,derinit,derend)
        call spline(jpt,jzz,derinit,derend,jetzzc)
      case(4)
        call define_derivatives_spline(jpt,jvx,derinit,derend)
        call spline(jpt,jvx,derinit,derend,jetvxc)
      case(5)
        call define_derivatives_spline(jpt,jvy,derinit,derend)
        call spline(jpt,jvy,derinit,derend,jetvyc)
      case(6)
        call define_derivatives_spline(jpt,jvz,derinit,derend)
        call spline(jpt,jvz,derinit,derend,jetvzc)
      case(7)
        call define_derivatives_spline(jpt,jst,derinit,derend)
        call spline(jpt,jst,derinit,derend,jetstc)
      case(8)
        call define_derivatives_spline(jpt,jms,derinit,derend)
        call spline(jpt,jms,derinit,derend,jetmsc)
      case(9)
        call define_derivatives_spline(jpt,jch,derinit,derend)
        call spline(jpt,jch,derinit,derend,jetchc)
      case(10)
        call define_derivatives_spline(jpt,jcr,derinit,derend)
        call spline(jpt,jcr,derinit,derend,jetcrc)
      case default
        continue
    end select
    
  enddo
    
  call sum_world_darr(jetxxc,npjetspline+1)
  call sum_world_darr(jetyyc,npjetspline+1)
  call sum_world_darr(jetzzc,npjetspline+1)
  call sum_world_darr(jetvxc,npjetspline+1)
  call sum_world_darr(jetvyc,npjetspline+1)
  call sum_world_darr(jetvzc,npjetspline+1)
  call sum_world_darr(jetstc,npjetspline+1)
  call sum_world_darr(jetmsc,npjetspline+1)
  call sum_world_darr(jetchc,npjetspline+1)
  call sum_world_darr(jetcrc,npjetspline+1)
  
  return
  
 end subroutine create_spline
 
 subroutine driver_fit_spline(i,myinit,myend,jpt,jvr,jpf,jvrf)

!***********************************************************************
!     
!     JETSPIN subroutine for driving the cubic 
!     spline interpolation of the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer,intent(in) :: i,myinit,myend
  double precision, allocatable, dimension(:), intent(in) :: jpt
  double precision, allocatable, dimension(:), intent(in) :: jvr
  double precision, allocatable, dimension(:), intent(in) :: jpf
  double precision, allocatable, dimension(:), intent(inout) :: jvrf
  
  select case(i)
    case(1)
      call fit_spline(myinit,myend,jpt,jvr,jetxxc,jpf,jvrf)
    case(2)
      call fit_spline(myinit,myend,jpt,jvr,jetyyc,jpf,jvrf)
    case(3)
      call fit_spline(myinit,myend,jpt,jvr,jetzzc,jpf,jvrf)
    case(4)
      call fit_spline(myinit,myend,jpt,jvr,jetvxc,jpf,jvrf)
    case(5)
      call fit_spline(myinit,myend,jpt,jvr,jetvyc,jpf,jvrf)
    case(6)
      call fit_spline(myinit,myend,jpt,jvr,jetvzc,jpf,jvrf)
    case(7)
      call fit_spline(myinit,myend,jpt,jvr,jetstc,jpf,jvrf)
    case(8)
      call fit_spline(myinit,myend,jpt,jvr,jetmsc,jpf,jvrf)
    case(9)
      call fit_spline(myinit,myend,jpt,jvr,jetchc,jpf,jvrf)
    case(10)
      call fit_spline(myinit,myend,jpt,jvr,jetcrc,jpf,jvrf)
    case default
      continue
  end select
  
  return
  
 end subroutine driver_fit_spline
 
 subroutine fit_spline(myinit,myend,jpt,jvr,jvrc,jpf,jvrf)

!***********************************************************************
!     
!     JETSPIN subroutine for interpolating the cubic 
!     spline interpolations of the nanojet given the points in jpt
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************

  implicit none
  
  integer,intent(in) :: myinit,myend
  double precision, allocatable, dimension(:), intent(in) :: jpt
  double precision, allocatable, dimension(:), intent(in) :: jvr
  double precision, allocatable, dimension(:), intent(in) :: jvrc
  double precision, allocatable, dimension(:), intent(in) :: jpf
  double precision, allocatable, dimension(:), intent(inout) :: jvrf
  
  integer :: i
  
  logical :: lerror
  
  do i=myinit+idrank,myend,mxrank
    call splint(jpt, jvr, jvrc,jpf(i), jvrf(i), lerror)
    if(lerror)call error(12)
  enddo
  
  call sum_world_darr(jvrf,myend+1)
  
  return
  
 end subroutine fit_spline
 
 subroutine define_derivatives_spline(jpt,jft,dinit,dend)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the first derivatives at the 
!     extremes in order to perform a natural spline interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension(:), intent(in) :: jpt,jft
  double precision, intent(out) :: dinit,dend
  
  integer :: ipoint,i
  double precision, dimension(3) :: xpoints,fpoints,interpc
  
  i=0
  do ipoint=inpjetspline,inpjetspline+2
    i=i+1
    xpoints(i)=jpt(ipoint)
    fpoints(i)=jft(ipoint)
  enddo
  call findcurve3(xpoints,fpoints,interpc)
  dinit=intdercurve3(jpt(inpjetspline),interpc)
  
  
  i=0
  do ipoint=npjetspline-2,npjetspline
    i=i+1
    xpoints(i)=jpt(ipoint)
    fpoints(i)=jft(ipoint)
  enddo
  call findcurve3(xpoints,fpoints,interpc)
  dend=intdercurve3(jpt(npjetspline),interpc)
  
  return
  
 end subroutine define_derivatives_spline
 
 subroutine spline(x,y,yp1,ypn,y2)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the cubic 
!     spline interpolation coefficients
!     (adopted from Numerical Recipes in FORTRAN 77)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  
  double precision, intent(in) :: yp1, ypn
  double precision, allocatable, intent(in), dimension(:) :: x,y 
  double precision, allocatable, intent(inout), dimension(:) :: y2
  
  integer :: i,k
  double precision :: p, qn, sig, un
  integer, save :: nmaxu
  double precision, allocatable, save :: u(:)
  
  if(allocated(u))then
    if(mxnpjet > nmaxu)then
      nmaxu=mxnpjet
      deallocate(u)
      allocate(u(0:mxnpjet))
    endif
  else 
    nmaxu=mxnpjet
    allocate(u(0:mxnpjet))
  endif

  if(yp1 > 0.99d30)then
    y2(inpjetspline)=0.d0
    u(inpjetspline)=0.d0
  else
    y2(inpjetspline)=-0.5d0
    u(inpjetspline)=(3.d0/(x(inpjetspline+1)-x(inpjetspline)))* &
     ((y(inpjetspline+1)-y(inpjetspline))/(x(inpjetspline+1)- &
     x(inpjetspline))-yp1)
  endif

  do i=inpjetspline+1, npjetspline-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/   &
     (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo

  if(ypn > 0.99d30) then
    qn=0.d0
    un=0.d0
  else
    qn=0.5d0
    un=(3.d0/(x(npjetspline)-x(npjetspline-1)))*(ypn-(y(npjetspline)- &
     y(npjetspline-1))/(x(npjetspline)-x(npjetspline-1)))
  endif

  y2(npjetspline)=(un-qn*u(npjetspline-1))/(qn*y2(npjetspline-1)+1.d0)

  do k=npjetspline-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  return
  
 end subroutine spline
 
 subroutine splint(xa,ya,y2a,x,y,lerror,y1,y2)

!***********************************************************************
!     
!     JETSPIN subroutine for interpolating a cubic 
!     spline
!     (adopted from Numerical Recipes in FORTRAN 77)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in) :: x
  double precision, intent(inout) :: y 
  double precision, allocatable, intent(in), dimension(:) :: xa,ya,y2a
  logical, intent(out) :: lerror
  double precision, intent(inout), optional :: y1,y2 
  integer :: k, khi, klo
  double precision :: a, b, h

  klo=inpjetspline
  khi=npjetspline
  do while((khi-klo)>1)
    k=(khi+klo)/2
    if(xa(k)>x) then
      khi=k
    else
      klo=k
    endif
  enddo
  
  h=xa(khi)-xa(klo)
  
  lerror=.false.
  if(h==0.d0)lerror=.true.
  
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))* &
   (h**2.d0)/6.d0
   
  if(present(y1))then
    y1=((3.d0*a**2.d0-1.d0)/6.d0)*(xa(khi)-xa(klo))*y2a(klo)+ &
       ((3.d0*b**2.d0-1.d0)/6.d0)*(xa(khi)-xa(klo))*y2a(khi)
  endif
  
  if(present(y2))then
    y2=a*y2a(klo)+b*y2a(khi)
  endif
  
  return
  
 end subroutine splint
 
 subroutine allocate_arrayakima()

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating arrays which are used to store
!     the spline coefficients if the Akima spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  
  if(larrayakima)then
    if(mxnpjet>narrayakima)then
      deallocate(jetak1,jetak2,jetak3,jetak4,mak,tak)
      narrayakima=mxnpjet
      allocate(jetak1(0:narrayakima),jetak2(0:narrayakima), &
       jetak3(0:narrayakima),jetak4(0:narrayakima))
      allocate(mak(-2:narrayakima+1),tak(0:narrayakima))
    endif
  else
    narrayakima=mxnpjet
    allocate(jetak1(0:narrayakima),jetak2(0:narrayakima), &
     jetak3(0:narrayakima),jetak4(0:narrayakima))
    allocate(mak(-2:narrayakima+1),tak(0:narrayakima))
    larrayakima=.true.
  endif
  
  return
  
 end subroutine allocate_arrayakima
 
 subroutine fit_akima(jptinit,jptend,jpt,jvr,jptc,jvrfit)

!***********************************************************************
!     
!     JETSPIN subroutine for performing the Akima 
!     spline interpolation of the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jptinit,jptend
  double precision, allocatable, dimension (:), intent(in) :: jpt
  double precision, allocatable, dimension (:), intent(in) :: jvr
  double precision, allocatable, dimension (:), intent(in) :: jptc
  double precision, allocatable, dimension (:), intent(inout) :: jvrfit
  
  double precision :: derinit,derend
  integer :: i,j
  
  inpjetspline=inpjet
  npjetspline=npjet
  
  call setup_akima( jpt, jvr, jetak1,jetak2,jetak3,jetak4)
  
  call interp_akima(jptinit,jptend,jpt,jetak1,jetak2,jetak3,jetak4, &
    jptc,jvrfit)
  
  return
  
 end subroutine fit_akima
 
 subroutine setup_akima( xpt, ypt, p0, p1, p2, p3)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the Akima 
!     spline interpolation coefficients
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, dimension(:),allocatable, intent(in) :: xpt, ypt
  double precision, dimension(:), allocatable, intent(inout) :: p0, p1, p2, p3
  
  integer :: i,ndarray
  double precision :: m1, m2, m3, m4, w1, w2
  double precision :: t1, t2, dx
  double precision, parameter :: eps = 1d-30
  
  mak(-2:narrayakima+1)=0.d0
! segment slopes are computed
  do i=mystart,myend
    mak(i)=(ypt(i+1)-ypt(i))/(xpt(i+1)-xpt(i))
  end do
  
  ndarray=narrayakima+4
  call sum_world_darr(mak,ndarray)
    
! segment slopes for the initial and last points
  mak(inpjetspline-1) = 2.d0*mak(inpjetspline) - mak(inpjetspline+1)
  mak(inpjetspline-2) = 2.d0*mak(inpjetspline-1) - mak(inpjetspline)
  mak(npjetspline) = 2.d0*mak(npjetspline-1) - mak(npjetspline-2)
  mak(npjetspline+1) = 2.d0*mak(npjetspline) - mak(npjetspline-1)
  
  tak(0:narrayakima)=0.d0
! slope at knots are computed
  do i=mystart,myend
    m1=mak(i-2)
    m2=mak(i-1)
    m3=mak(i)
    m4=mak(i+1)
    w1=dabs(m4-m3)
    w2=dabs(m2-m1)
    if(w1<eps .and. w2<eps)then
!   the division by zero is avoided
      tak(i)=0.5d0*(m2 + m3)
    else
      tak(i)=(w1*m2+w2*m3)/(w1+w2)
    end if
  end do
  
  ndarray=narrayakima+1
  call sum_world_darr(tak,ndarray)
  
  
  p0(0:narrayakima)=0.d0
  p1(0:narrayakima)=0.d0
  p2(0:narrayakima)=0.d0
  p3(0:narrayakima)=0.d0
! compute the polynomial cofficients
  do i=mystart,myend
    if(i==npjetspline)cycle
    dx=xpt(i+1)-xpt(i)
    t1=tak(i)
    t2=tak(i+1)
    p0(i)=ypt(i)
    p1(i)=t1
    p2(i)=(3.d0*mak(i)-2.d0*t1-t2)/dx
    p3(i)=(t1+t2-2.d0*mak(i))/dx**2.d0
  end do
  
  ndarray=narrayakima+1
  call sum_world_darr(p0,ndarray)
  call sum_world_darr(p1,ndarray)
  call sum_world_darr(p2,ndarray)
  call sum_world_darr(p3,ndarray)
  
  return
  
 end subroutine setup_akima
 
 subroutine interp_akima(iinterp,ninterp,xpt,p0,p1,p2,p3,x,y,dydx)

!***********************************************************************
!     
!     JETSPIN subroutine for interpolating the Akima spline
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
    
  implicit none
  
  integer, intent(in) :: iinterp,ninterp
  double precision, allocatable, dimension(:), intent(in) :: xpt
  double precision, allocatable, dimension(:), intent(in) :: p0,p1,p2,p3
  double precision, allocatable, intent(in) :: x(:)
  
  double precision, allocatable, intent(inout) :: y(:)
  double precision, allocatable, intent(inout), optional :: dydx(:)
  
  integer :: i, j, k
  double precision :: dx
  
! interpolate at each point
  do i = iinterp+idrank,ninterp,mxrank
!   look for the location in array (use end segments if out of bounds)
    if(x(i)<xpt(inpjetspline))then
      j = inpjetspline
    else
!     look for the index
      do j=npjetspline-1,inpjetspline,-1
        if(x(i)>=xpt(j))then
          exit
        end if
      end do
    end if
!   evaluate polynomial (and derivative if requested)
    dx=(x(i)-xpt(j))
    y(i)=p0(j)+p1(j)*dx+p2(j)*dx**2.d0+p3(j)*dx**3.d0
    if(present(dydx))then
      dydx(i)=p1(j)+ 2.d0*p2(j)*dx+3.d0*p3(j)*dx**2.d0
    endif
  end do
  
  call sum_world_darr(y,ninterp+1)
  if(present(dydx))call sum_world_darr(dydx,ninterp+1)
  
  return
 
 end subroutine interp_akima
 
 subroutine looking_indexes_2(ipoint,jpoint,indexwork,nstep,jetptl)

!***********************************************************************
!     
!     JETSPIN subroutine for defining the indexes of arrays which 
!     describe the nanojet before and after the adding of a bead by
!     the dynamic refinement
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer,intent(in) :: ipoint,jpoint,nstep
  
  integer, allocatable, intent(inout) :: indexwork(:,:)
  logical, allocatable, intent(inout) :: jetptl(:)
  
  integer :: i,j,startcopy(2),endcopy(2),befipoint,istart
  logical ::lfirst,lredo
  
  lfirst=.true.
  lredo=.true.
  j=0
  i=inpjet-1
  
  do while(lredo)
    i=i+1
    if(jetptl(i))then
      j=j+1
      if(j==jpoint-1)befipoint=i
      if(j==jpoint)then
        lredo=.false.
      endif
    endif
  enddo
  
  
  if(jpoint==1)then
    if(ipoint==inpjet)then
      startcopy(1)=inpjet
      startcopy(2)=inpjet
      endcopy(1)=inpjet
      endcopy(2)=inpjet
    else
      startcopy(1)=inpjet
      startcopy(2)=inpjet
      endcopy(1)=ipoint
      endcopy(2)=ipoint
    endif
  else
    istart=befipoint+1
    startcopy(1)=istart+(jpoint-1)
    startcopy(2)=istart
    endcopy(1)=ipoint+(jpoint-1)
    endcopy(2)=ipoint
  endif
    
    
  indexwork(jpoint,1)=ipoint
  indexwork(jpoint,2)=startcopy(1)
  indexwork(jpoint,3)=startcopy(2)
  indexwork(jpoint,4)=endcopy(1)
  indexwork(jpoint,5)=endcopy(2)
  indexwork(jpoint,6)=0
  
  return
  
 end subroutine looking_indexes_2
 
 subroutine driver_fit_curve2(jpoint,npoint,indexwork,ypt,yvr,yvrf)

!***********************************************************************
!     
!     JETSPIN subroutine for performing a linear interpolation
!     if any two nanofiber beads are beyond a given threshold
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jpoint,npoint
  integer, allocatable, intent(inout) :: indexwork(:,:)
  double precision, allocatable, dimension(:), intent(in) :: ypt
  double precision, allocatable, dimension(:), intent(in) :: yvr
  double precision, allocatable, dimension(:), intent(inout) :: yvrf
  
  integer :: i,k,newnpjet
  integer, dimension(2) :: startcopy,endcopy
  double precision, dimension(2) :: myx,myfx,coeffout
  double precision :: fitx,fitfx
  
  i=indexwork(jpoint,1)
  k=indexwork(jpoint,2)
  startcopy(1)=indexwork(jpoint,3)
  startcopy(2)=indexwork(jpoint,4)
  endcopy(1)=indexwork(jpoint,5)
  endcopy(2)=indexwork(jpoint,6)
  

  myx(1)=ypt(i)
  myx(2)=ypt(i+1)
  myfx(1)=yvr(i)
  myfx(2)=yvr(i+1)
  fitx=(myx(2)-myx(1))/2.d0+myx(1)
      
  call findcurve2(myx,myfx,coeffout)
  fitfx=intcurve2(fitx,coeffout)
  
  yvrf(startcopy(1):endcopy(1))=yvr(startcopy(2):endcopy(2))

  yvrf(endcopy(1)+1)=fitfx
  
  newnpjet=npjet+npoint
  if(jpoint==npoint)then
    yvrf(endcopy(1)+2:newnpjet)=yvr(endcopy(2)+1:npjet)
  endif
  
  return
  
 end subroutine driver_fit_curve2
 
 subroutine findcurve2(xsub,fxsub,coeffoutsub)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the coefficients of
!     the linear interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in), dimension(2) :: xsub,fxsub
  double precision, intent(out), dimension(2) :: coeffoutsub
  
  double precision :: detnom,detdenom
  double precision, dimension(2,2) :: matnom,matdenom
  
  
  matdenom(1:2,1)=xsub(1:2)
  matdenom(1:2,2)=1.d0
  
  call compute_det2(matdenom,detdenom)
  
  matnom=matdenom
  matnom(1:2,1)=fxsub(1:2)
  
  call compute_det2(matnom,detnom)
  coeffoutsub(1)=detnom/detdenom
  
  matnom=matdenom
  matnom(1:2,2)=fxsub(1:2)
  
  call compute_det2(matnom,detnom)
  coeffoutsub(2)=detnom/detdenom
  
  return
  
 end subroutine findcurve2
 
 subroutine compute_det2(matrix,determinant)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the determinant of
!     a 2x2 matrix
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in), dimension(2,2) :: matrix
  double precision, intent(out) :: determinant
  
  determinant=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
  
  return
  
 end subroutine compute_det2
 
 function intcurve2(xsub,coeffinsub)

!***********************************************************************
!     
!     JETSPIN function for interpolating the data by a linear
!     interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: xsub
  double precision, intent(in), dimension(2) :: coeffinsub
  
  double precision :: intcurve2
  
  intcurve2=coeffinsub(1)*xsub+coeffinsub(2)
  
  return
  
 end function intcurve2
 
 function intdercurve2(xsub,coeffinsub)

!***********************************************************************
!     
!     JETSPIN function for interpolating the first derivative by a 
!     linear interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: xsub
  double precision, intent(in), dimension(2) :: coeffinsub
  
  double precision :: intdercurve2
  
  intdercurve2=coeffinsub(1)
  
  return
  
 end function intdercurve2
 
 subroutine findcurve3(xsub,fxsub,coeffoutsub)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the coefficients of
!     the quadratic interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in), dimension(3) :: xsub,fxsub
  double precision, intent(out), dimension(3) :: coeffoutsub
  
  integer :: i
  double precision :: detnom,detdenom
  double precision, dimension(3,3) :: matnom,matdenom
  

  matdenom(1:3,1)=xsub(1:3)**2.d0
  matdenom(1:3,2)=xsub(1:3)
  matdenom(1:3,3)=1.d0
  
  call compute_det3(matdenom,detdenom)
  
  matnom=matdenom
  matnom(1:3,1)=fxsub(1:3)
  
  call compute_det3(matnom,detnom)
  coeffoutsub(1)=detnom/detdenom
  
  matnom=matdenom
  matnom(1:3,2)=fxsub(1:3)
  
  call compute_det3(matnom,detnom)
  coeffoutsub(2)=detnom/detdenom
  
  matnom=matdenom
  matnom(1:3,3)=fxsub(1:3)
  
  call compute_det3(matnom,detnom)
  coeffoutsub(3)=detnom/detdenom
  
  return
  
 end subroutine findcurve3
 
 subroutine compute_det3(matrix,determinant)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the determinant of
!     a 3x3 matrix
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in), dimension(3,3) :: matrix
  double precision, intent(out) :: determinant
  
  determinant=matrix(1,1)*matrix(2,2)*matrix(3,3)+ &
   matrix(1,2)*matrix(2,3)*matrix(3,1)+ &
   matrix(1,3)*matrix(2,1)*matrix(3,2)- &
   matrix(1,3)*matrix(2,2)*matrix(3,1)- &
   matrix(1,2)*matrix(2,1)*matrix(3,3)- &
   matrix(1,1)*matrix(2,3)*matrix(3,2)
  
  
  return
  
 end subroutine compute_det3
 
 function intcurve3(xsub,coeffinsub)

!***********************************************************************
!     
!     JETSPIN function for interpolating the data by a quadratic
!     interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: xsub
  double precision, intent(in), dimension(3) :: coeffinsub
  
  double precision :: intcurve3
  
  intcurve3=coeffinsub(1)*(xsub**2.d0)+coeffinsub(2)*xsub+coeffinsub(3)
  
  return
  
 end function intcurve3
 
 function intdercurve3(xsub,coeffinsub)

!***********************************************************************
!     
!     JETSPIN function for interpolating the first derivative by a 
!     quadratic interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: xsub
  double precision, intent(in), dimension(3) :: coeffinsub
  
  double precision :: intdercurve3
  
  intdercurve3=coeffinsub(1)*(2.d0*xsub)+coeffinsub(2)
  
  return
  
 end function intdercurve3
 
 subroutine looking_indexes_4(ipoint,jpoint,indexwork,nstep,jetptl)

!***********************************************************************
!     
!     JETSPIN subroutine for defining the indexes of arrays which 
!     describe the nanojet before and after the adding of a bead by
!     the dynamic refinement
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer,intent(in) :: ipoint,jpoint,nstep
  integer, allocatable, intent(inout) :: indexwork(:,:)
  logical, allocatable, intent(inout) :: jetptl(:)
  
  integer :: i,j,k,startcopy(2),endcopy(2),befipoint,istart
  logical ::lfirst,lredo
  
  lfirst=.true.
  lredo=.true.
  j=0
  i=inpjet-1
  
  do while(lredo)
    i=i+1
    if(jetptl(i))then
      j=j+1
      if(j==jpoint-1)befipoint=i
      if(j==jpoint)then
        if(i-1>=inpjet .and. i+2<=npjet)then
          k=2
        else
          if(i-1<inpjet)then
            k=1
          else
            k=3
          endif
        endif
        lredo=.false.
      endif
    endif
  enddo
  
  
  if(jpoint==1)then
    if(ipoint==inpjet)then
      startcopy(1)=inpjet
      startcopy(2)=inpjet
      endcopy(1)=inpjet
      endcopy(2)=inpjet
    else
      startcopy(1)=inpjet
      startcopy(2)=inpjet
      endcopy(1)=ipoint
      endcopy(2)=ipoint
    endif
  else
    istart=befipoint+1
    startcopy(1)=istart+(jpoint-1)
    startcopy(2)=istart
    endcopy(1)=ipoint+(jpoint-1)
    endcopy(2)=ipoint
  endif
    
    
  indexwork(jpoint,1)=ipoint
  indexwork(jpoint,2)=startcopy(1)
  indexwork(jpoint,3)=startcopy(2)
  indexwork(jpoint,4)=endcopy(1)
  indexwork(jpoint,5)=endcopy(2)
  indexwork(jpoint,6)=k
  
  return
  
 end subroutine looking_indexes_4
 
 subroutine findcurve4(xsub,fxsub,coeffoutsub)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the coefficients of
!     the cubic interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in), dimension(4) :: xsub,fxsub
  double precision, intent(out), dimension(4) :: coeffoutsub
  
  double precision :: detnom,detdenom
  double precision, dimension(4,4) :: matnom,matdenom
  
  matdenom(1:4,1)=xsub(1:4)**3.d0
  matdenom(1:4,2)=xsub(1:4)**2.d0
  matdenom(1:4,3)=xsub(1:4)
  matdenom(1:4,4)=1.d0
  
  call compute_det4(matdenom,detdenom)
  
  matnom=matdenom
  matnom(1:4,1)=fxsub(1:4)
  
  call compute_det4(matnom,detnom)
  coeffoutsub(1)=detnom/detdenom
  
  matnom=matdenom
  matnom(1:4,2)=fxsub(1:4)
  
  call compute_det4(matnom,detnom)
  coeffoutsub(2)=detnom/detdenom
  
  matnom=matdenom
  matnom(1:4,3)=fxsub(1:4)
  
  call compute_det4(matnom,detnom)
  coeffoutsub(3)=detnom/detdenom
  
  matnom=matdenom
  matnom(1:4,4)=fxsub(1:4)
  
  call compute_det4(matnom,detnom)
  coeffoutsub(4)=detnom/detdenom
  
  return
  
 end subroutine findcurve4
 
 subroutine compute_det4(matrix,determinant)

!***********************************************************************
!     
!     JETSPIN subroutine for computing the determinant of
!     a 4x4 matrix
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in), dimension(4,4) :: matrix
  double precision, intent(out) :: determinant

  determinant =  matrix(1,1)*(matrix(2,2)*(matrix(3,3)*matrix(4,4)- &
   matrix(3,4)*matrix(4,3))+matrix(2,3)*(matrix(3,4)*matrix(4,2)- &
   matrix(3,2)*matrix(4,4))+matrix(2,4)*(matrix(3,2)*matrix(4,3)- &
   matrix(3,3)*matrix(4,2)))-matrix(1,2)*(matrix(2,1)*(matrix(3,3)* &
   matrix(4,4)-matrix(3,4)*matrix(4,3))+matrix(2,3)*(matrix(3,4)* &
   matrix(4,1)-matrix(3,1)*matrix(4,4))+matrix(2,4)*(matrix(3,1)* &
   matrix(4,3)-matrix(3,3)*matrix(4,1)))+matrix(1,3)*(matrix(2,1)* &
   (matrix(3,2)*matrix(4,4)-matrix(3,4)*matrix(4,2))+matrix(2,2)* &
   (matrix(3,4)*matrix(4,1)-matrix(3,1)*matrix(4,4))+matrix(2,4)* &
   (matrix(3,1)*matrix(4,2)-matrix(3,2)*matrix(4,1)))-matrix(1,4)* &
   (matrix(2,1)*(matrix(3,2)*matrix(4,3)-matrix(3,3)*matrix(4,2))+ &
   matrix(2,2)*(matrix(3,3)*matrix(4,1)-matrix(3,1)*matrix(4,3))+ &
   matrix(2,3)*(matrix(3,1)*matrix(4,2)-matrix(3,2)*matrix(4,1)))
             
  return
  
 end subroutine compute_det4
 
 function intcurve4(xsub,coeffinsub)

!***********************************************************************
!     
!     JETSPIN function for interpolating the data by a cubic
!     interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: xsub
  double precision, intent(in), dimension(4) :: coeffinsub
  
  double precision :: intcurve4
  
  intcurve4=coeffinsub(1)*(xsub**3.d0)+ &
  coeffinsub(2)*(xsub**2.d0)+coeffinsub(3)*xsub+coeffinsub(4)
  
  return
  
 end function intcurve4
 
 function intdercurve4(xsub,coeffinsub)

!***********************************************************************
!     
!     JETSPIN function for interpolating the first derivative by a 
!     cubic interpolation
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: xsub
  double precision, intent(in), dimension(4) :: coeffinsub
  
  double precision :: intdercurve4
  
  intdercurve4=coeffinsub(1)*3.d0*(xsub**2.d0)+ &
   coeffinsub(2)*2.d0*xsub+coeffinsub(3)
  
  return
  
 end function intdercurve4
 
 subroutine cubic_interpolation(i,k,ypt,yvr,fitfx)

!***********************************************************************
!     
!     JETSPIN subroutine for performing a cubic interpolation
!     if any two nanofiber beads are beyond a given threshold
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: i,k
  double precision, allocatable, dimension(:), intent(in) :: ypt
  double precision, allocatable, dimension(:), intent(in) :: yvr
  double precision, intent(out) :: fitfx
  
  double precision, dimension(4) :: myx,myfx,coeffout
  double precision :: fitx
  
  
  select case (k)
    case(1)
      myx(1)=ypt(i)
      myx(2)=ypt(i+1)
      myx(3)=ypt(i+2)
      myx(4)=ypt(i+3)
      myfx(1)=yvr(i)
      myfx(2)=yvr(i+1)
      myfx(3)=yvr(i+2)
      myfx(4)=yvr(i+3)
      fitx=(myx(2)-myx(1))/2.d0+myx(1)
    case (2)
      myx(1)=ypt(i-1)
      myx(2)=ypt(i)
      myx(3)=ypt(i+1)
      myx(4)=ypt(i+2)
      myfx(1)=yvr(i-1)
      myfx(2)=yvr(i)
      myfx(3)=yvr(i+1)
      myfx(4)=yvr(i+2)
      fitx=(myx(3)-myx(2))/2.d0+myx(2)
    case (3)
      myx(1)=ypt(i-2)
      myx(2)=ypt(i-1)
      myx(3)=ypt(i)
      myx(4)=ypt(i+1)
      myfx(1)=yvr(i-2)
      myfx(2)=yvr(i-1)
      myfx(3)=yvr(i)
      myfx(4)=yvr(i+1)
      fitx=(myx(4)-myx(3))/2.d0+myx(3)
  end select
      
  call findcurve4(myx,myfx,coeffout)
  fitfx=intcurve4(fitx,coeffout)
  
  return
  
 end subroutine cubic_interpolation
 
 end module fit_mod
