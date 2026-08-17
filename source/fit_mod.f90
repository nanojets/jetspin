
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
 logical, save :: lakima_accelerator_coordinates=.false.
 
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
 public :: begin_akima_accelerator_data
 public :: end_akima_accelerator_data
 
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
 
 subroutine fit_akima(jptinit,jptend,jpt,jvr,jptc,jvrfit, &
  use_accelerator,field_name)

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
  logical, intent(in), optional :: use_accelerator
  character(len=*), intent(in), optional :: field_name

  logical :: run_accelerator
  logical :: lmonotonefield
#ifdef _OPENACC
  double precision :: coefficient_max_abs,coefficient_max_rel
#ifdef JETSPIN_COMPARE_AKIMA
  character(len=16) :: diagnostic_name
  double precision, allocatable :: host_reference(:)
  double precision :: value_max_abs,value_max_rel
#endif
#endif
  
  inpjetspline=inpjet
  npjetspline=npjet

! Never physically meaningful to overshoot outside the local data range
! for these fields (cross-section/evaporation area, mass/charge density);
! position/velocity/stress keep the unmodified Akima behaviour already
! validated by the existing test suite.
  lmonotonefield=.false.
  if(present(field_name))then
    lmonotonefield=trim(field_name)=='radius_area' .or. &
     trim(field_name)=='evap_radius_area' .or. &
     trim(field_name)=='mass_density' .or. &
     trim(field_name)=='charge_density'
  endif

  run_accelerator=.false.
#ifdef _OPENACC
  if(present(use_accelerator))run_accelerator=use_accelerator
  run_accelerator=run_accelerator .and. mxrank==1
#ifdef JETSPIN_DEV_HOST_AKIMA
  run_accelerator=.false.
#endif
#endif

#ifdef _OPENACC
  if(run_accelerator)then
#ifdef JETSPIN_COMPARE_AKIMA
    allocate(host_reference(lbound(jvrfit,1):ubound(jvrfit,1)))
    host_reference=jvrfit
    call setup_akima(jpt,jvr,jetak1,jetak2,jetak3,jetak4,lmonotonefield)
    call interp_akima(jptinit,jptend,jpt,jetak1,jetak2,jetak3,jetak4, &
     jptc,host_reference)
#endif
    call fit_akima_accelerator(jptinit,jptend,jpt,jvr,jptc,jvrfit, &
     coefficient_max_abs,coefficient_max_rel,lmonotonefield)
#ifdef JETSPIN_COMPARE_AKIMA
    value_max_abs=maxval(dabs(jvrfit(jptinit:jptend)- &
     host_reference(jptinit:jptend)))
    value_max_rel=value_max_abs/max( &
     maxval(dabs(host_reference(jptinit:jptend))),1.d-30)
    diagnostic_name='unnamed'
    if(present(field_name))diagnostic_name=field_name
    if(idrank==0)write(6,'(a,a,4(a,es12.4))') &
     'Akima device comparison: field=',trim(diagnostic_name), &
     ' coefficient_max_abs=',coefficient_max_abs, &
     ' coefficient_max_rel=',coefficient_max_rel, &
     ' value_max_abs=',value_max_abs,' value_max_rel=',value_max_rel
    deallocate(host_reference)
#endif
  else
#endif
    call setup_akima(jpt,jvr,jetak1,jetak2,jetak3,jetak4,lmonotonefield)
    call interp_akima(jptinit,jptend,jpt,jetak1,jetak2,jetak3,jetak4, &
     jptc,jvrfit)
#ifdef _OPENACC
  endif
#endif

! Development-only diagnostic: expose the endpoint segment slopes, the
! extrapolated tangent actually used at the lower boundary knot, and the
! raw (pre-dabs) fitted value there and at its global minimum, to check
! whether the classic Akima endpoint-extrapolation formula is producing an
! overshoot/undershoot for this field somewhere in the fitted segment
! (inpjetspline==jptinit is the collector-side endpoint here).
  if(idrank==0 .and. present(field_name))then
    if(trim(field_name)=='radius_area' .or. trim(field_name)=='mass_density' &
     .or. trim(field_name)=='charge_density' .or. trim(field_name)=='stress' &
     .or. trim(field_name)=='vx' .or. trim(field_name)=='vy' &
     .or. trim(field_name)=='vz')then
      write(6,'(a,a,6(a,es14.6),a,i0)')'Akima endpoint diagnostic: field=', &
       trim(field_name), &
       ' m_first_segment=',mak(inpjetspline), &
       ' m_second_segment=',mak(inpjetspline+1), &
       ' tangent_at_endpoint=',tak(inpjetspline), &
       ' raw_fit_at_endpoint=',jvrfit(jptinit), &
       ' raw_fit_min=',minval(jvrfit(jptinit:jptend)), &
       ' raw_fit_max=',maxval(jvrfit(jptinit:jptend)), &
       ' raw_fit_min_at=',jptinit-1+minloc(jvrfit(jptinit:jptend),1)
    endif
  endif

  return

 end subroutine fit_akima

 subroutine begin_akima_accelerator_data(xpt,x,enable)

  implicit none

  double precision, allocatable, intent(in) :: xpt(:),x(:)
  logical, intent(in) :: enable

  lakima_accelerator_coordinates=.false.
#ifdef _OPENACC
  if(enable)then
! Normalized source coordinates and target knots are shared by all eleven
! interpolated fields. Map them once per accepted event instead of once per
! field.
!$acc enter data copyin(xpt,x)
    lakima_accelerator_coordinates=.true.
  endif
#endif

  return

 end subroutine begin_akima_accelerator_data

 subroutine end_akima_accelerator_data(xpt,x)

  implicit none

  double precision, allocatable, intent(in) :: xpt(:),x(:)

#ifdef _OPENACC
  if(lakima_accelerator_coordinates)then
!$acc exit data delete(xpt,x)
  endif
#endif
  lakima_accelerator_coordinates=.false.

  return

 end subroutine end_akima_accelerator_data

#ifdef _OPENACC
 subroutine fit_akima_accelerator(iinterp,ninterp,xpt,ypt,x,y, &
  coefficient_max_abs,coefficient_max_rel,lmonotone)

!***********************************************************************
!
!     OpenACC Akima coefficient construction and interpolation.
!     Segment slopes, knot tangents, polynomial coefficients, and target
!     points are independent. Only the four endpoint-slope extrapolations,
!     and the optional monotonicity limiter below, use one serial device
!     thread.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: iinterp,ninterp
  double precision, allocatable, intent(in) :: xpt(:),ypt(:),x(:)
  double precision, allocatable, intent(inout) :: y(:)
  double precision, intent(out) :: coefficient_max_abs
  double precision, intent(out) :: coefficient_max_rel
  logical, intent(in), optional :: lmonotone

  integer :: i,j,ipass
  double precision :: m1,m2,m3,m4,w1,w2,t1,t2,dx
  double precision :: coefficient_scale
  double precision :: alpha,beta,tau
  double precision, parameter :: eps=1.d-30
  double precision, allocatable :: slopes(:),tangents(:)
  double precision, allocatable :: p0(:),p1(:),p2(:),p3(:)
  logical :: uselimiter

  allocate(slopes(inpjetspline-2:npjetspline+1))
  allocate(tangents(inpjetspline:npjetspline))
  allocate(p0(inpjetspline:npjetspline-1))
  allocate(p1(inpjetspline:npjetspline-1))
  allocate(p2(inpjetspline:npjetspline-1))
  allocate(p3(inpjetspline:npjetspline-1))

#ifdef JETSPIN_COMPARE_AKIMA
!$acc data create(xpt,ypt,x,slopes,tangents,y) &
!$acc& copyout(p0,p1,p2,p3)
#else
!$acc data create(xpt,ypt,x,slopes,tangents,p0,p1,p2,p3,y)
#endif
! The source arrays can already belong to the persistent jet mapping. Their
! host values include the event-time density/cross-section preparation, so
! explicitly refresh the active ranges before constructing coefficients.
  if(.not.lakima_accelerator_coordinates)then
!$acc update device(xpt(inpjetspline:npjetspline),x(iinterp:ninterp))
  endif
!$acc update device(ypt(inpjetspline:npjetspline))

!$acc parallel loop gang vector present(xpt,ypt,slopes)
  do i=inpjetspline,npjetspline-1
    slopes(i)=(ypt(i+1)-ypt(i))/(xpt(i+1)-xpt(i))
  enddo
!$acc end parallel loop

!$acc serial present(slopes)
  slopes(inpjetspline-1)=2.d0*slopes(inpjetspline)- &
   slopes(inpjetspline+1)
  slopes(inpjetspline-2)=2.d0*slopes(inpjetspline-1)- &
   slopes(inpjetspline)
  slopes(npjetspline)=2.d0*slopes(npjetspline-1)- &
   slopes(npjetspline-2)
  slopes(npjetspline+1)=2.d0*slopes(npjetspline)- &
   slopes(npjetspline-1)
!$acc end serial

!$acc parallel loop gang vector present(slopes,tangents)
  do i=inpjetspline,npjetspline
    m1=slopes(i-2)
    m2=slopes(i-1)
    m3=slopes(i)
    m4=slopes(i+1)
    w1=dabs(m4-m3)
    w2=dabs(m2-m1)
    if(w1<eps .and. w2<eps)then
      tangents(i)=0.5d0*(m2+m3)
    else
      tangents(i)=(w1*m2+w2*m3)/(w1+w2)
    endif
  enddo
!$acc end parallel loop

! Device counterpart of limit_akima_tangents_monotone (setup_akima, host):
! the classic Fritsch-Carlson sufficient condition for a monotone cubic
! Hermite interpolant, restricting the tangents just computed above. Same
! algorithm, applied to slopes/tangents instead of the host's module-level
! mak/tak. A single serial device thread, like the endpoint extrapolation
! just above: this only ever runs once per accepted refinement event, on
! a segment of at most a few hundred knots, so there is no performance
! reason to parallelize an inherently sequential fixed-point iteration.
  uselimiter=.false.
  if(present(lmonotone))uselimiter=lmonotone
  if(uselimiter)then
!$acc serial present(slopes,tangents)
    do i=inpjetspline,npjetspline
      if(slopes(i-1)*slopes(i)<=0.d0)tangents(i)=0.d0
    enddo
    do ipass=1,3
      do i=inpjetspline,npjetspline-1
        if(dabs(slopes(i))<eps)then
          tangents(i)=0.d0
          tangents(i+1)=0.d0
          cycle
        endif
        alpha=tangents(i)/slopes(i)
        beta=tangents(i+1)/slopes(i)
        if(alpha<0.d0)tangents(i)=0.d0
        if(beta<0.d0)tangents(i+1)=0.d0
        alpha=tangents(i)/slopes(i)
        beta=tangents(i+1)/slopes(i)
        if(alpha**2.d0+beta**2.d0>9.d0)then
          tau=3.d0/dsqrt(alpha**2.d0+beta**2.d0)
          tangents(i)=tau*alpha*slopes(i)
          tangents(i+1)=tau*beta*slopes(i)
        endif
      enddo
    enddo
!$acc end serial
  endif

!$acc parallel loop gang vector present(xpt,ypt,slopes,tangents,p0,p1,p2,p3)
  do i=inpjetspline,npjetspline-1
    dx=xpt(i+1)-xpt(i)
    t1=tangents(i)
    t2=tangents(i+1)
    p0(i)=ypt(i)
    p1(i)=t1
    p2(i)=(3.d0*slopes(i)-2.d0*t1-t2)/dx
    p3(i)=(t1+t2-2.d0*slopes(i))/dx**2.d0
  enddo
!$acc end parallel loop

! Each target knot owns one thread. The descending interval search is local
! to that knot and exactly follows the historical host selection rule.
!$acc parallel loop gang vector private(j,dx) present(xpt,x,y,p0,p1,p2,p3)
  do i=iinterp,ninterp
    if(x(i)<xpt(inpjetspline))then
      j=inpjetspline
    else
      do j=npjetspline-1,inpjetspline,-1
        if(x(i)>=xpt(j))exit
      enddo
    endif
    dx=x(i)-xpt(j)
    y(i)=p0(j)+p1(j)*dx+p2(j)*dx**2.d0+p3(j)*dx**3.d0
  enddo
!$acc end parallel loop
! Only the new interpolated values are required on the host. The prefix in
! the shared service buffer was never overwritten on the device.
!$acc update self(y(iinterp:ninterp))
!$acc end data

  coefficient_max_abs=0.d0
  coefficient_max_rel=0.d0
#ifdef JETSPIN_COMPARE_AKIMA
  coefficient_max_abs=max( &
   maxval(dabs(p0-jetak1(inpjetspline:npjetspline-1))), &
   maxval(dabs(p1-jetak2(inpjetspline:npjetspline-1))), &
   maxval(dabs(p2-jetak3(inpjetspline:npjetspline-1))), &
   maxval(dabs(p3-jetak4(inpjetspline:npjetspline-1))))
  coefficient_scale=max( &
   maxval(dabs(jetak1(inpjetspline:npjetspline-1))), &
   maxval(dabs(jetak2(inpjetspline:npjetspline-1))), &
   maxval(dabs(jetak3(inpjetspline:npjetspline-1))), &
   maxval(dabs(jetak4(inpjetspline:npjetspline-1))))
  coefficient_max_rel=coefficient_max_abs/max(coefficient_scale,1.d-30)
#endif

  deallocate(slopes,tangents,p0,p1,p2,p3)

  return

 end subroutine fit_akima_accelerator
#endif
 
 subroutine setup_akima( xpt, ypt, p0, p1, p2, p3, lmonotone)

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
  logical, intent(in), optional :: lmonotone

  integer :: i,ndarray
  double precision :: m1, m2, m3, m4, w1, w2
  double precision :: t1, t2, dx
  double precision, parameter :: eps = 1d-30
  logical :: uselimiter
  
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

! Classic Akima gives no guarantee that the resulting cubic stays within
! the range spanned by the two knots of any segment it interpolates --
! it can overshoot/undershoot well beyond it even when the source data
! is itself smooth, and that overshoot compounds across repeated fits
! (an accepted dynamic-refinement event re-fits from the previous
! event's already-fitted state). Limiting the knot tangents to satisfy
! the Fritsch-Carlson sufficient condition for a monotone cubic Hermite
! removes that oscillation at its source, for every fit, regardless of
! how the data has already evolved. This is opted into only for fields
! where an overshoot is never physically meaningful (bead cross-section
! and evaporation area); position/velocity/stress keep the unmodified
! Akima behaviour already validated by the existing test suite.
  uselimiter=.false.
  if(present(lmonotone))uselimiter=lmonotone
  if(uselimiter)call limit_akima_tangents_monotone()

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

 subroutine limit_akima_tangents_monotone()

!***********************************************************************
!
!     JETSPIN subroutine for limiting the Akima knot tangents (tak) in
!     place to satisfy the Fritsch-Carlson sufficient condition for a
!     monotone cubic Hermite interpolant on every segment where the
!     underlying source data (mak, the segment secant slopes) is itself
!     monotone. A knot where the two adjacent secant slopes disagree in
!     sign is a local extremum of the source data; any nonzero tangent
!     there would force the cubic to overshoot past it, so the tangent
!     is zeroed. On every remaining segment the two endpoint tangents
!     are then rescaled, if needed, to stay within the classic
!     circle-of-radius-3 bound relative to the segment's own secant
!     slope. Operates directly on the module-level mak/tak arrays over
!     the current inpjetspline:npjetspline source range; iterated a few
!     passes since a tangent shared between two segments can only
!     shrink further on a later pass, never grow back, so this converges
!     quickly to a fixed point.
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2026
!
!***********************************************************************

  implicit none

  double precision, parameter :: eps=1.d-30
  double precision :: alpha,beta,tau
  integer :: i,ipass

  do i=inpjetspline,npjetspline
    if(mak(i-1)*mak(i)<=0.d0)tak(i)=0.d0
  enddo

  do ipass=1,3
    do i=inpjetspline,npjetspline-1
      if(dabs(mak(i))<eps)then
        tak(i)=0.d0
        tak(i+1)=0.d0
        cycle
      endif
      alpha=tak(i)/mak(i)
      beta=tak(i+1)/mak(i)
      if(alpha<0.d0)tak(i)=0.d0
      if(beta<0.d0)tak(i+1)=0.d0
      alpha=tak(i)/mak(i)
      beta=tak(i+1)/mak(i)
      if(alpha**2.d0+beta**2.d0>9.d0)then
        tau=3.d0/dsqrt(alpha**2.d0+beta**2.d0)
        tak(i)=tau*alpha*mak(i)
        tak(i+1)=tau*beta*mak(i)
      endif
    enddo
  enddo

  return

 end subroutine limit_akima_tangents_monotone
 
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
