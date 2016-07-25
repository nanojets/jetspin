
 module breaking_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which manage the 
!     simulation of the nanofiber breaking up
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
 use version_mod, only : mystart,myend,or_world_larr
 use utility_mod, only : Pi,sig
 use nanojet_mod, only : inpjet,npjet,jetms,jetms,jetxx,jetyy,jetzz, &
                   jetpt,jetvl,jetcr,lbreakup,jetbr,jetbd,typemass, &
                   ltagbeads,jetvx,jetvy,jetvz,jetch,jetst
 use fit_mod,     only : findcurve4,jetbdc,jetbrc,jetptc,jetxxc, &
                   jetyyc,jetzzc,jetvxc,jetvyc,jetvzc,jetstc,jetmsc, &
                   jetchc,jetcrc
 use support_functions_mod, only : compute_length_path,compute_crosssec
 
 implicit none
 
 private
 
 integer, public, save :: condition_breakup=0
 
 public :: ckeck_breakup
 public :: clean_breakup
 
 contains
 
 subroutine ckeck_breakup(timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for checking if a breaking up is happening
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(in) :: timesub
  
  integer :: istart,iend,ipoint
  double precision :: lengthpathsub
  
  if(.not. lbreakup)return
  
  if(mystart==inpjet)then
    istart=mystart+1
  else
    istart=mystart
  endif
  
  if(myend==npjet)then
    iend=myend-2
  else
    iend=myend
  endif
  
  call compute_length_path(jetxx,jetyy,jetzz,jetpt,lengthpathsub)
  
  call compute_crosssec(jetxx,jetyy,jetzz,jetvl,jetcr)
  
  select case(condition_breakup)
  case(1)
    do ipoint=istart,iend
      call condition_breakup_1(ipoint,jetbr(ipoint))
    enddo
  case default
    do ipoint=istart,iend
      call condition_breakup_1(ipoint,jetbr(ipoint))
    enddo
  end select
  
  call or_world_larr(jetbr,npjet+1)
  
  if(typemass==3 .or. ltagbeads)then
    do ipoint=inpjet,npjet
      if(jetbr(ipoint))then
        jetbd(ipoint)=.true.
        jetbd(ipoint+1)=.true.
      endif
    enddo
  endif
  
  return
  
 end subroutine ckeck_breakup
 
 subroutine condition_breakup_1(ipoint,lchecksub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for checking if a breaking up condition is
!     verified. 
!     CONDITION 1: The condition is verified if a cubic interpolating
!     curve is zero.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: ipoint
  logical, intent(inout) :: lchecksub
  
  double precision, dimension(4) :: myx,myfx,coeffout
  double precision :: a,b,c,q,r,theta,ar,br,x1,x2,x3
  logical :: lcheck
  
  myx(1)=jetpt(ipoint-1)
  myx(2)=jetpt(ipoint)
  myx(3)=jetpt(ipoint+1)
  myx(4)=jetpt(ipoint+2)
  myfx(1)=jetcr(ipoint-1)
  myfx(2)=jetcr(ipoint)
  myfx(3)=jetcr(ipoint+1)
  myfx(4)=jetcr(ipoint+2)
  call findcurve4(myx,myfx,coeffout)
  a=coeffout(2)/coeffout(1)
  b=coeffout(3)/coeffout(1)
  c=coeffout(4)/coeffout(1)
  q=(a**2.d0-3.d0*b)/9.d0
  r=(2.d0*a**3.d0-9.d0*a*b+27.d0*c)/54.d0
  lcheck=.false.
  if(r**2.d0<q**3.d0)then
    theta=dacos(r/(dsqrt(q**3.d0)))
    x1=-2.d0*dsqrt(q)*dcos(theta/3.d0)-a/3.d0
    x2=-2.d0*dsqrt(q)*dcos((theta+2.d0*pi)/3.d0)-a/3.d0
    x3=-2.d0*dsqrt(q)*dcos((theta-2.d0*pi)/3.d0)-a/3.d0
    if(x1>myx(1) .and. x1<myx(4))lcheck=.true.
    if(x2>myx(1) .and. x2<myx(4))lcheck=.true.
    if(x3>myx(1) .and. x3<myx(4))lcheck=.true.
  else
    ar=-1.d0*sig(r)*(dabs(r)+dsqrt(r**2.d0-q**3.d0))**(1.d0/3.d0)
    if(ar==0.d0)then
      br=0.d0
    else
      br=q/ar
    endif
    x1=(ar+br)-a/3.d0
    if(x1>myx(1) .and. x1<myx(4))lcheck=.true.
  endif
  
  lchecksub = ( lchecksub .or. lcheck )
  
  return
  
 end subroutine condition_breakup_1
 
 subroutine clean_breakup(jptinit,jptend,totjptend)
  
!***********************************************************************
!     
!     JETSPIN subroutine for removing extra beads after the dynamic
!     refunement procedure: an extra bead is the bead inserted in
!     a point where the filament is break up
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jptinit,jptend,totjptend
  
  integer :: ipoint,j
  logical :: ltest
  
  
  jetbdc(inpjet:npjet)=jetbd(inpjet:npjet)
  jetbrc(inpjet:npjet)=jetbr(inpjet:npjet)
  jetptc(inpjet:npjet)=jetvl(inpjet:npjet)
  jetxxc(inpjet:npjet)=jetxx(inpjet:npjet)
  jetyyc(inpjet:npjet)=jetyy(inpjet:npjet)
  jetzzc(inpjet:npjet)=jetzz(inpjet:npjet)
  jetvxc(inpjet:npjet)=jetvx(inpjet:npjet)
  jetvyc(inpjet:npjet)=jetvy(inpjet:npjet)
  jetvzc(inpjet:npjet)=jetvz(inpjet:npjet)
  jetstc(inpjet:npjet)=jetst(inpjet:npjet)
  jetmsc(inpjet:npjet)=jetms(inpjet:npjet)
  jetchc(inpjet:npjet)=jetch(inpjet:npjet)
  jetcrc(inpjet:npjet)=jetcr(inpjet:npjet)
  
  jetbd(inpjet:npjet)=.false.
  jetbr(inpjet:npjet)=.false.
  jetvl(inpjet:npjet)=0.d0
  jetxx(inpjet:npjet)=0.d0
  jetyy(inpjet:npjet)=0.d0
  jetzz(inpjet:npjet)=0.d0
  jetvx(inpjet:npjet)=0.d0
  jetvy(inpjet:npjet)=0.d0
  jetvz(inpjet:npjet)=0.d0
  jetst(inpjet:npjet)=0.d0
  jetms(inpjet:npjet)=0.d0
  jetch(inpjet:npjet)=0.d0
  jetcr(inpjet:npjet)=0.d0
  
  ltest=.false.
  j=jptinit-1
  do ipoint=jptinit,jptend
    if(jetbr(ipoint))then
      ltest=.true.
      j=j+1
      jetbd(j)=jetbdc(ipoint)
      jetbr(j)=jetbrc(ipoint)
      jetvl(j)=jetptc(ipoint)
      jetxx(j)=jetxxc(ipoint)
      jetyy(j)=jetyyc(ipoint)
      jetzz(j)=jetzzc(ipoint)
      jetvx(j)=jetvxc(ipoint)
      jetvy(j)=jetvyc(ipoint)
      jetvz(j)=jetvzc(ipoint)
      jetst(j)=jetstc(ipoint)
      jetms(j)=jetmsc(ipoint)
      jetch(j)=jetchc(ipoint)
      jetcr(j)=jetcrc(ipoint)
    else
      if(jetbd(ipoint))ltest=.false.
      if(.not. ltest)then
        j=j+1
        jetbd(j)=jetbdc(ipoint)
        jetbr(j)=jetbrc(ipoint)
        jetvl(j)=jetptc(ipoint)
        jetxx(j)=jetxxc(ipoint)
        jetyy(j)=jetyyc(ipoint)
        jetzz(j)=jetzzc(ipoint)
        jetvx(j)=jetvxc(ipoint)
        jetvy(j)=jetvyc(ipoint)
        jetvz(j)=jetvzc(ipoint)
        jetst(j)=jetstc(ipoint)
        jetms(j)=jetmsc(ipoint)
        jetch(j)=jetchc(ipoint)
        jetcr(j)=jetcrc(ipoint)
      endif
    endif
  enddo
  
  do ipoint=jptend+1,totjptend
    j=j+1
    jetbd(j)=jetbdc(ipoint)
    jetbr(j)=jetbrc(ipoint)
    jetvl(j)=jetptc(ipoint)
    jetxx(j)=jetxxc(ipoint)
    jetyy(j)=jetyyc(ipoint)
    jetzz(j)=jetzzc(ipoint)
    jetvx(j)=jetvxc(ipoint)
    jetvy(j)=jetvyc(ipoint)
    jetvz(j)=jetvzc(ipoint)
    jetst(j)=jetstc(ipoint)
    jetms(j)=jetmsc(ipoint)
    jetch(j)=jetchc(ipoint)
    jetcr(j)=jetcrc(ipoint)
  enddo
  
  inpjet=jptinit
  npjet=j
  
  return
  
 end subroutine clean_breakup
 
 end module breaking_mod
 
