 module coulomb_force_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which manage and compute 
!     the Coulomb force 
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification April 2016
!     
!***********************************************************************
 use version_mod,           only : idrank,mxrank,sum_world_darr, &
                             max_world_darr,sum_world_iarr
 use error_mod
 use utility_mod,           only : Pi,modulvec,cross,dot,sig
 use nanojet_mod,           only : jetch,jetfr,h,inpjet,npjet,q,&
                             linserted,dresolution,jetms,systype,fcut, &
                             thresolution,lmirror,mxnpjet,nmulstep, &
                             lmultiplestep,dcutoff,lneighlistdo,tstep, &
                             lmultisteperror,multisteperror,lremove, &
                             jetfm,nmulstepdone,nmultisteperror, &
                             ldcutoff,incnpjet,maxdispl,lmaxdispl
 use support_functions_mod, only : beadlength1d,beadlength, &
                             compute_crosssec
 
 implicit none
 
 private
 
 logical, save :: lcomputfder=.false.
 logical, save :: lmscomputed=.false.
 
 integer, save :: ncoulforce=0
 integer, save :: maxneighlist=50
 integer, save :: maxlistchunk=0
 
 integer, save :: msinpjet=0
 integer, save :: msnpjet=0
 
 double precision, save :: smoothedcharge
 double precision, save :: timemultistep
 
 double precision, allocatable, save :: coulservicearr(:,:)
 double precision, allocatable, save :: coulforcems(:,:)
 double precision, allocatable, save :: oldcoulforcems(:,:)
 double precision, allocatable, save :: dvcoulforcems(:,:)
 double precision, allocatable, save :: coulcrossec(:)
 
 integer, save, allocatable, dimension(:,:) :: neighlist
 integer, save, allocatable, dimension(:) :: neighlentry
 integer, save, allocatable, dimension(:) :: neighlentrymirr
 
 double precision, allocatable, public, save :: coulforce(:,:)
 
 integer, save :: ncoulcrossec=0
 
 public :: smooth_charge
 public :: restore_charge
 public :: compute_coulomelec_driver
 
 contains
 
 subroutine allocate_coulcrossec(imiomax)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating coulcrossec array necessary for computing
!     storaging the cross section
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification April 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax
  
  
  if(ncoulcrossec/=0)then
    if(imiomax>ncoulcrossec)then
      deallocate(coulcrossec)
      ncoulcrossec=imiomax
      allocate(coulcrossec(0:ncoulcrossec))
    endif
  else
    ncoulcrossec=imiomax
    allocate(coulcrossec(0:ncoulcrossec))
  endif
  
  return
  
 end subroutine allocate_coulcrossec
 
 subroutine smooth_charge(yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for smoothing the charge of the last added
!     bead if the mutual distance between this bead and the nozzle is
!     less than two times the resolution length step (dresolution)
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) ::  yxx
  double precision, allocatable, dimension (:), intent(inout), &
    optional ::  yyy
  double precision, allocatable, dimension (:), intent(inout), &
   optional ::  yzz
  
  double precision :: dtemp
 
  if(.not.linserted)then
    select case(systype)
    case(1:2)
      dtemp=fcut(beadlength1d(npjet-2,yxx),thresolution,dresolution)
    case default
      dtemp=fcut(beadlength(npjet-2,yxx,yyy,yzz),thresolution, &
       dresolution)
    end select
    smoothedcharge=jetch(npjet-1)
    jetch(npjet-1)=dtemp*smoothedcharge
  endif
  
 return
 
 end subroutine smooth_charge
 
 subroutine restore_charge()
 
!***********************************************************************
!     
!     JETSPIN subroutine for restoring the smoothed charge of the last 
!     added bead
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
 
  implicit none
 
  if(.not.linserted)then
    jetch(npjet-1)=smoothedcharge
  endif
  
 return
 
 end subroutine restore_charge
 
 subroutine compute_coulomelec_driver(nstep,timesub,ycf,yvl,yxx,yyy, &
   yzz)
  
!***********************************************************************
!     
!     JETSPIN subroutine for driving the Coulomb forces subroutines
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification April 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) ::  nstep
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(in) ::  yvl(:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in) ::  yyy(:)
  double precision, allocatable, intent(in) ::  yzz(:)
  
  call allocate_coulcrossec(npjet)
  call compute_crosssec(yxx,yyy,yzz,yvl,coulcrossec)
  select case(systype)
  case(1)
    if(lmultiplestep)then
      call compute_coulomelec_multistep(nstep,timesub,ycf,yxx)
    else
      call compute_coulomelec(nstep,timesub,ycf,yxx)
    endif
  case default
    if(lmultiplestep)then
      call compute_coulomelec_multistep(nstep,timesub,ycf,yxx,yyy,yzz)
    else
      call compute_coulomelec(nstep,timesub,ycf,yxx,yyy,yzz)
    endif
  end select
  
  return
  
 end subroutine compute_coulomelec_driver
  
 subroutine compute_coulomelec(nstep,timesub,ycf,yxx,yyy,yzz)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb forces of a
!     n-body system by the direct summation method
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification july 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) ::  nstep
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  double precision, parameter :: onethird=1.d0/(dsqrt(3.d0))
  
  integer :: ipoint,jpoint,imiomax
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint,dtemp
  double precision, dimension(3) :: versor,utang
  
  
  select case(systype)
    case(1)
    
      imiomax=mxnpjet
      if(ncoulforce/=0)then
        if(imiomax>ncoulforce)then
          deallocate(ycf)
          ncoulforce=imiomax
          allocate(ycf(0:ncoulforce,1))
        endif
      else
        ncoulforce=imiomax
        allocate(ycf(0:ncoulforce,1))
      endif
  
      ycf(0:ncoulforce,1:1)=0.d0
      
!     compute the Coulomb forces
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        do jpoint=ipoint+1,npjet
          if(jetfr(jpoint))cycle
          if(ldcutoff)then
            dtemp=dabs(yxx(jpoint)-yxx(ipoint))
            if(dtemp>dcutoff)cycle
          endif
          ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
           ((dabs(yxx(jpoint)-yxx(ipoint))+coulcrossec(jpoint))**2.d0)
          ycf(jpoint,1)=ycf(jpoint,1)-1.d0*(jetch(jpoint)*Qt)/ &
           ((dabs(yxx(jpoint)-yxx(ipoint))+coulcrossec(jpoint))**2.d0)
        enddo
        if(lmirror)then
          do jpoint=inpjet,npjet
            if(jetfr(jpoint))cycle
            xjpoint=dabs(yxx(jpoint)-h)+h
    !        if(ldcutoff)then
    !          dtemp=dabs(xjpoint-yxx(ipoint))
    !          if(dtemp>dcutoff)cycle
    !        endif
            ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
             ((dabs(xjpoint-yxx(ipoint))+coulcrossec(jpoint))**2.d0)
          enddo
        endif
      enddo
      
      imiomax=ncoulforce+1
      call sum_world_darr(ycf,imiomax)
      
    case default
      
      imiomax=mxnpjet
      if(ncoulforce/=0)then
        if(imiomax>ncoulforce)then
          deallocate(ycf)
          ncoulforce=imiomax
          allocate(ycf(0:ncoulforce,3))
        endif
      else
        ncoulforce=imiomax
        allocate(ycf(0:ncoulforce,3))
      endif
      
      ycf(0:ncoulforce,1:3)=0.d0
      
!     compute the Coulomb forces
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        do jpoint=ipoint+1,npjet
          if(jetfr(jpoint))cycle
          utang(1)=yxx(ipoint)-yxx(jpoint)
          utang(2)=yyy(ipoint)-yyy(jpoint)
          utang(3)=yzz(ipoint)-yzz(jpoint)
          norm = modulvec(utang)
   !       if(ldcutoff)then
   !         if(norm>dcutoff)cycle
   !       endif
          if(norm>1.d-30)then
            versor(1)=utang(1)/norm
            versor(2)=utang(2)/norm
            versor(3)=utang(3)/norm
            ycf(ipoint,1:3)=ycf(ipoint,1:3)+ &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
             versor(1:3)
            ycf(jpoint,1:3)=ycf(jpoint,1:3)- &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
             versor(1:3)
          endif
        enddo
        if(lmirror)then
          do jpoint=inpjet,npjet
            if(jetfr(jpoint))cycle
            xjpoint=dabs(yxx(jpoint)-h)+h
            yjpoint=yyy(jpoint)
            zjpoint=yzz(jpoint)
            utang(1)=yxx(ipoint)-xjpoint
            utang(2)=yyy(ipoint)-yjpoint
            utang(3)=yzz(ipoint)-zjpoint
            norm = modulvec(utang)
            if(ldcutoff)then
              if(norm>dcutoff)cycle
            endif
            if(norm>1.d-30)then
              versor(1)=utang(1)/norm
              versor(2)=utang(2)/norm
              versor(3)=utang(3)/norm
              ycf(ipoint,1:3)=ycf(ipoint,1:3)- &
               (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
               versor(1:3)
            endif
          enddo
        endif
      enddo
      
      imiomax=(ncoulforce+1)*3
      call sum_world_darr(ycf,imiomax)
      
  end select
       
  return
  
 end subroutine compute_coulomelec
 
 subroutine compute_coulomelec_multistep(nstep,timesub,ycf,yxx,yyy,yzz)
  
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb forces of a
!     n-body system by the multiple step approach
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) ::  nstep
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  double precision, parameter :: onethird=1.d0/(dsqrt(3.d0))
  
  integer :: ipoint,jpoint,imiomax,isub,jsub
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint,meanerrms,dt
  double precision, dimension(3) :: versor,utang
  
  logical :: lneighlistdosub,ldodirectsum
  double precision, save :: myoldtime=0.d0
  
  integer :: inpjetmio,npjetmio
  double precision :: daje(3),daje2(3)
  logical, save :: lfirst=.true.
  logical, save :: lcheckerr=.true.
  
  
  lneighlistdosub=.false.
  meanerrms=0.d0
  
  
  if(lremove)then
    do ipoint=inpjet,npjet
      if(jetfr(ipoint))then
        if(.not. jetfm(ipoint))lneighlistdosub=.true.
      endif
    enddo
  endif
  
  select case(systype)
    case(1)
      
      call allocate_coularrays(mxnpjet,ycf,1,lneighlistdosub)
      ldodirectsum=(mod(nstep,nmulstep)==0)
      if(ldodirectsum)lcheckerr=.true.
      lneighlistdo=(lneighlistdo .or. lneighlistdosub .or. ldodirectsum)
      lneighlistdo=(lneighlistdo .and. (.not. lcomputfder))
      call list_test(lneighlistdo,timesub,yxx)
      if(lneighlistdo)then
        lneighlistdo=.false.
        lcomputfder=.true.
        ycf(:,:)=0.d0
        call compute_neighlist_and_coulomelec(ycf,oldcoulforcems, &
         timesub,yxx)
        myoldtime=timesub
        nmulstepdone=nmulstepdone+1
        lmscomputed=.false.
      elseif(lcomputfder)then
        dt=timesub-myoldtime
        ycf(:,:)=0.d0
        call compute_coulomelec_and_external(nstep,ycf, &
         coulforcems,timesub,yxx)
        if(dt/=0.d0)then
          
          lcomputfder=.false.
!         compute first derivative of force
          isub=0
          do ipoint=inpjet+idrank,npjet,mxrank
            if(jetfr(ipoint))cycle
            isub=isub+1
            dvcoulforcems(isub,1)=(coulforcems(isub,1)- &
             oldcoulforcems(isub,1))/dt
          enddo
          
!         compute time corresponding to the derivative
          timemultistep=(timesub+myoldtime)/2.d0
          
!         compute relative forces
          isub=0
          do ipoint=inpjet+idrank,npjet,mxrank
            if(jetfr(ipoint))cycle
            isub=isub+1
            coulforcems(isub,1)=(coulforcems(isub,1)+ &
             oldcoulforcems(isub,1))/2.d0
          enddo
          
          lmscomputed=.true.
          
        endif 
      elseif(lmscomputed)then
        
        if(mod(nstep+1,nmulstep)==0)then
          if(lcheckerr)then
            lcheckerr=.false.
            call compute_multistep_error(nstep,lmultisteperror, &
             meanerrms,timesub,ycf,oldcoulforcems,yxx)
            multisteperror=meanerrms+multisteperror
            nmultisteperror=nmultisteperror+1
          else
            ycf(:,:)=0.d0
            call compute_inner_coulomelec(timesub,ycf,yxx)
            call compute_outer_coulomelec(timesub,ycf,yxx)
            imiomax=(ncoulforce+1)
            call sum_world_darr(ycf,imiomax)
          endif
        else
          ycf(:,:)=0.d0
          call compute_inner_coulomelec(timesub,ycf,yxx)
          call compute_outer_coulomelec(timesub,ycf,yxx)
          imiomax=(ncoulforce+1)
          call sum_world_darr(ycf,imiomax) 
        endif
      else
        call error(17)
      endif
      
    case default
      
      call allocate_coularrays(mxnpjet,ycf,3,lneighlistdosub)
      ldodirectsum=(mod(nstep,nmulstep)==0)
      if(ldodirectsum)lcheckerr=.true.
      lneighlistdo=(lneighlistdo .or. lneighlistdosub .or. ldodirectsum)
      lneighlistdo=(lneighlistdo .and. (.not. lcomputfder))
      call list_test(lneighlistdo,timesub,yxx,yyy,yzz)
      if(lneighlistdo)then
        
        lneighlistdo=.false.
        lcomputfder=.true.
        ycf(:,:)=0.d0
        call compute_neighlist_and_coulomelec(ycf,oldcoulforcems, &
         timesub,yxx,yyy,yzz)
        myoldtime=timesub
        nmulstepdone=nmulstepdone+1
        lmscomputed=.false.
      elseif(lcomputfder)then
        
        dt=timesub-myoldtime
        ycf(:,:)=0.d0
        call compute_coulomelec_and_external(nstep,ycf, &
         coulforcems,timesub,yxx,yyy,yzz)
         
        if(dt/=0.d0)then
          lcomputfder=.false.
!         compute first derivative of force
          isub=0
          do ipoint=inpjet+idrank,npjet,mxrank
            if(jetfr(ipoint))cycle
            isub=isub+1
            dvcoulforcems(isub,1:3)=(coulforcems(isub,1:3)- &
             oldcoulforcems(isub,1:3))/dt
          enddo
          
!         compute time corresponding to the derivative
          timemultistep=(timesub+myoldtime)/2.d0
          
!         compute relative forces
          isub=0
          do ipoint=inpjet+idrank,npjet,mxrank
            if(jetfr(ipoint))cycle
            isub=isub+1
            coulforcems(isub,1:3)=(coulforcems(isub,1:3)+ &
             oldcoulforcems(isub,1:3))/2.d0
          enddo
          
          lmscomputed=.true.
          
        endif
         
      elseif(lmscomputed)then
        
        if(mod(nstep+1,nmulstep)==0)then
          if(lcheckerr)then
            lcheckerr=.false.
            call compute_multistep_error(nstep,lmultisteperror, &
             meanerrms,timesub,ycf,oldcoulforcems,yxx,yyy,yzz)
            multisteperror=meanerrms+multisteperror
            nmultisteperror=nmultisteperror+1
          else
            ycf(:,:)=0.d0
            call compute_inner_coulomelec(timesub,ycf,yxx,yyy,yzz)
            call compute_outer_coulomelec(timesub,ycf,yxx,yyy,yzz)
            imiomax=(ncoulforce+1)*3
            call sum_world_darr(ycf,imiomax)
          endif
        else
          ycf(:,:)=0.d0
          call compute_inner_coulomelec(timesub,ycf,yxx,yyy,yzz)
          call compute_outer_coulomelec(timesub,ycf,yxx,yyy,yzz)
          imiomax=(ncoulforce+1)*3
          call sum_world_darr(ycf,imiomax)
        endif
      else
        call error(17)
      endif
        
  end select
  
  return
  
 end subroutine compute_coulomelec_multistep
 
 subroutine allocate_coularrays(imiomax,ycf,indarr,lneighlistdosub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating arrays necessary for computing
!     the Coulomb forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: imiomax,indarr
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  logical, intent(inout), optional :: lneighlistdosub
  
  integer :: mychunk
  
  
  if(ncoulforce/=0)then
    if(imiomax>ncoulforce)then
      deallocate(ycf)
      ncoulforce=imiomax
      allocate(ycf(0:ncoulforce,indarr))
      if(lmultiplestep)then
        if(present(lneighlistdosub))lneighlistdosub=.true.
        deallocate(coulforcems)
        deallocate(oldcoulforcems)
        deallocate(dvcoulforcems)
        deallocate(neighlentry)
        deallocate(coulservicearr)
        mychunk=ceiling(dble(imiomax)/dble(mxrank))
        allocate(coulforcems(mychunk,indarr))
        allocate(oldcoulforcems(mychunk,indarr))
        allocate(dvcoulforcems(mychunk,indarr))
        allocate(neighlentry(mychunk))
        allocate(coulservicearr(0:ncoulforce,indarr))
        if(lmirror)then
          deallocate(neighlentrymirr)
          deallocate(neighlist)
          allocate(neighlentrymirr(mychunk))
          allocate(neighlist(2*mychunk,maxneighlist))
        else
          deallocate(neighlist)
          allocate(neighlist(mychunk,maxneighlist))
        endif
      endif
    endif
  else
    ncoulforce=imiomax
    allocate(ycf(0:ncoulforce,indarr))
    if(lmultiplestep)then
      if(present(lneighlistdosub))lneighlistdosub=.true.
      mychunk=ceiling(dble(imiomax)/dble(mxrank))
      allocate(coulforcems(mychunk,indarr))
      allocate(oldcoulforcems(mychunk,indarr))
      allocate(dvcoulforcems(mychunk,indarr))
      allocate(neighlentry(mychunk))
      allocate(coulservicearr(0:ncoulforce,indarr))
      if(lmirror)then
        allocate(neighlentrymirr(mychunk))
        allocate(neighlist(2*mychunk,maxneighlist))
      else
        allocate(neighlist(mychunk,maxneighlist))
      endif
    endif
  endif
  
  return
  
 end subroutine allocate_coularrays
 
 subroutine compute_neighlist_and_coulomelec(ycf,outerycf, &
  timesub,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating arrays necessary for computing
!     the Coulomb forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(inout) :: outerycf(:,:)
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  integer :: isub,imiomax,ipoint,jpoint
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint,dtemp(1),dt,deno
  double precision, dimension(3) :: versor,utang
  
  ycf(:,:)=0.d0
  outerycf(:,:)=0.d0
  coulservicearr(:,:)=0.d0
  jetfm(:)=.false.
  jetfm(inpjet:npjet)=jetfr(inpjet:npjet)
  msinpjet=inpjet
  msnpjet=npjet
  neighlentry(:)=0
  neighlist(:,:)=0
  if(lmirror)neighlentrymirr(:)=0
  
  
  select case(systype)
    case(1)
!     compute the Coulomb forces
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        do jpoint=ipoint+1,npjet
          if(jetfr(jpoint))cycle
          norm=dabs(yxx(jpoint)-yxx(ipoint))
          ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
           ((norm+coulcrossec(jpoint))**2.d0)
          ycf(jpoint,1)=ycf(jpoint,1)-1.d0*(jetch(jpoint)*Qt)/ &
           ((norm+coulcrossec(jpoint))**2.d0)
          if(norm<=dcutoff)then
            neighlentry(isub)=neighlentry(isub)+1
            if(neighlentry(isub)>maxneighlist)then
              call reallocate_neighlist()
            endif
            neighlist(isub,neighlentry(isub))=jpoint
          else
            coulservicearr(ipoint,1)=coulservicearr(ipoint,1)+1.d0* &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)
            coulservicearr(jpoint,1)=coulservicearr(jpoint,1)-1.d0* &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)
          endif
        enddo
        if(lmirror)then
          neighlentrymirr(isub)=neighlentry(isub)
          do jpoint=inpjet,npjet
            if(jetfr(jpoint))cycle
            xjpoint=dabs(yxx(jpoint)-h)+h
            norm=dabs(xjpoint-yxx(ipoint))
            ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
             ((norm+coulcrossec(jpoint))**2.d0)
            if(norm<=dcutoff)then
              neighlentrymirr(isub)=neighlentrymirr(isub)+1
              if(neighlentrymirr(isub)>maxneighlist)then
                call reallocate_neighlist()
              endif
              neighlist(isub,neighlentrymirr(isub))=jpoint
            else
              coulservicearr(ipoint,1)=coulservicearr(ipoint,1)+1.d0* &
               (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)
            endif
          enddo
        endif
      enddo
      
      imiomax=ncoulforce+1
      call sum_world_darr(ycf,imiomax)
      call sum_world_darr(coulservicearr,imiomax)
      
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        outerycf(isub,1)=coulservicearr(ipoint,1)
      enddo
    
      
    case default
    
!     compute the Coulomb forces
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        do jpoint=ipoint+1,npjet
          if(jetfr(jpoint))cycle
          utang(1)=yxx(ipoint)-yxx(jpoint)
          utang(2)=yyy(ipoint)-yyy(jpoint)
          utang(3)=yzz(ipoint)-yzz(jpoint)
          norm = modulvec(utang)
          if(norm>1.d-30)then
            versor(1)=utang(1)/norm
            versor(2)=utang(2)/norm
            versor(3)=utang(3)/norm
            ycf(ipoint,1:3)=ycf(ipoint,1:3)+ &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
             versor(1:3)
            ycf(jpoint,1:3)=ycf(jpoint,1:3)- &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
             versor(1:3)
            if(norm<=dcutoff)then
              neighlentry(isub)=neighlentry(isub)+1
              if(neighlentry(isub)>maxneighlist)then
                call reallocate_neighlist()
              endif
              neighlist(isub,neighlentry(isub))=jpoint
            else
              coulservicearr(ipoint,1:3)=coulservicearr(ipoint,1:3)+ &
               (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
               versor(1:3)
              coulservicearr(jpoint,1:3)=coulservicearr(jpoint,1:3)- &
               (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
               versor(1:3)
            endif
          endif
        enddo
        if(lmirror)then
          neighlentrymirr(isub)=neighlentry(isub)
          do jpoint=inpjet,npjet
            if(jetfr(jpoint))cycle
            xjpoint=dabs(yxx(jpoint)-h)+h
            yjpoint=yyy(jpoint)
            zjpoint=yzz(jpoint)
            utang(1)=yxx(ipoint)-xjpoint
            utang(2)=yyy(ipoint)-yjpoint
            utang(3)=yzz(ipoint)-zjpoint
            norm = modulvec(utang)
            if(norm>1.d-30)then
              versor(1)=utang(1)/norm
              versor(2)=utang(2)/norm
              versor(3)=utang(3)/norm
              ycf(ipoint,1:3)=ycf(ipoint,1:3)- &
               (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
               versor(1:3)
              if(norm<=dcutoff)then
                neighlentrymirr(isub)=neighlentrymirr(isub)+1
                if(neighlentrymirr(isub)>maxneighlist)then
                  call reallocate_neighlist()
                endif
                neighlist(isub,neighlentrymirr(isub))=jpoint
              else
                coulservicearr(ipoint,1:3)=coulservicearr(ipoint,1:3)- &
                 (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)*&
                 versor(1:3)
              endif
            endif
          enddo
        endif
      enddo
      
      imiomax=(ncoulforce+1)*3
      call sum_world_darr(ycf,imiomax)
      call sum_world_darr(coulservicearr,imiomax)
      
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        outerycf(isub,1:3)=coulservicearr(ipoint,1:3)
      enddo
      
  end select
  
  
  return
  
 end subroutine compute_neighlist_and_coulomelec
 
 subroutine compute_coulomelec_and_external(nstep,ycf,outerycf, &
  timesub,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating arrays necessary for computing
!     the Coulomb forces
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(inout) :: outerycf(:,:)
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  integer :: isub,imiomax,ipoint,jpoint
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint,dtemp(1),dt,deno
  double precision, dimension(3) :: versor,utang
  
  ycf(:,:)=0.d0
  outerycf(:,:)=0.d0
  coulservicearr(:,:)=0.d0
  
  
  select case(systype)
    case(1)
    
!     compute the total Coulomb forces
      call compute_coulomelec(nstep,timesub,ycf,yxx)
      
!     compute the inner Coulomb forces
      call compute_inner_coulomelec(timesub,coulservicearr,yxx)
        
      imiomax=(ncoulforce+1)
      call sum_world_darr(coulservicearr,imiomax)
      
!     compute the outer Coulomb forces as difference
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        outerycf(isub,1)=ycf(ipoint,1)-coulservicearr(ipoint,1)
      enddo
      
    case default
    
!     compute the total Coulomb forces
      call compute_coulomelec(nstep,timesub,ycf,yxx,yyy,yzz)
      
!     compute the inner Coulomb forces
      call compute_inner_coulomelec(timesub,coulservicearr,yxx,yyy,yzz)
        
      imiomax=(ncoulforce+1)*3
      call sum_world_darr(coulservicearr,imiomax)
      
!     compute the outer Coulomb forces as difference
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        outerycf(isub,1:3)=ycf(ipoint,1:3)-coulservicearr(ipoint,1:3)
      enddo
      
  end select
  
  
  return
  
 end subroutine compute_coulomelec_and_external
 
 subroutine compute_outer_coulomelec(timesub,ycf,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb forces belonging to
!     the outer shell by the multiple step approach
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************

  implicit none
 
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  integer :: ipoint,isub
  double precision :: dt
 
  select case(systype)
    case(1)
      dt=timesub-timemultistep
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        ycf(ipoint,1)=ycf(ipoint,1)+coulforcems(isub,1)+ &
         dvcoulforcems(isub,1)*dt
      enddo
    case default
      dt=timesub-timemultistep
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        ycf(ipoint,1:3)=ycf(ipoint,1:3)+coulforcems(isub,1:3)+ &
         dvcoulforcems(isub,1:3)*dt
      enddo
  end select
      
  return
  
 end subroutine compute_outer_coulomelec
 
 subroutine compute_inner_coulomelec(timesub,ycf,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb forces belonging to
!     the inner shell by the multiple step approach
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
  
  implicit none
  
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  
  integer :: ipoint,isub,mymax,mymaxmirr,jsub,jpoint
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint
  double precision, dimension(3) :: versor,utang
  
  select case(systype)
    case(1)
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        mymax=neighlentry(isub)
        do jsub=1,mymax
          jpoint=neighlist(isub,jsub)
          norm=dabs(yxx(jpoint)-yxx(ipoint))
          ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
           ((norm+coulcrossec(jpoint))**2.d0)
          ycf(jpoint,1)=ycf(jpoint,1)-1.d0*(jetch(jpoint)*Qt)/ &
           ((norm+coulcrossec(jpoint))**2.d0)
        enddo
        if(lmirror)then
          mymaxmirr=neighlentrymirr(isub)
          do jsub=mymax+1,mymaxmirr
            jpoint=neighlist(isub,jsub)
            xjpoint=dabs(yxx(jpoint)-h)+h
            norm=dabs(xjpoint-yxx(ipoint))
            ycf(ipoint,1)=ycf(ipoint,1)+1.d0*(jetch(jpoint)*Qt)/ &
             ((norm+coulcrossec(jpoint))**2.d0)
          enddo
        endif
      enddo     
    case default
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        Qt=jetch(ipoint)*Q/jetms(ipoint)
        mymax=neighlentry(isub)
        do jsub=1,mymax
          jpoint=neighlist(isub,jsub)
          utang(1)=yxx(ipoint)-yxx(jpoint)
          utang(2)=yyy(ipoint)-yyy(jpoint)
          utang(3)=yzz(ipoint)-yzz(jpoint)
          norm = modulvec(utang)
          if(norm>1.d-30)then
            versor(1)=utang(1)/norm
            versor(2)=utang(2)/norm
            versor(3)=utang(3)/norm
            ycf(ipoint,1:3)=ycf(ipoint,1:3)+ &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
             versor(1:3)
            ycf(jpoint,1:3)=ycf(jpoint,1:3)- &
             (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
             versor(1:3)
          endif
        enddo
        if(lmirror)then
          mymaxmirr=neighlentrymirr(isub)
          do jsub=mymax+1,mymaxmirr
            jpoint=neighlist(isub,jsub)
            xjpoint=dabs(yxx(jpoint)-h)+h
            yjpoint=yyy(jpoint)
            zjpoint=yzz(jpoint)
            utang(1)=yxx(ipoint)-xjpoint
            utang(2)=yyy(ipoint)-yjpoint
            utang(3)=yzz(ipoint)-zjpoint
            norm = modulvec(utang)
            if(norm>1.d-30)then
              versor(1)=utang(1)/norm
              versor(2)=utang(2)/norm
              versor(3)=utang(3)/norm
              ycf(ipoint,1:3)=ycf(ipoint,1:3)- &
               (jetch(jpoint)*Qt)/((norm+coulcrossec(jpoint))**2.d0)* &
               versor(1:3)
            endif
          enddo
        endif
      enddo
  end select
  
  return
      
 end subroutine compute_inner_coulomelec
 
 subroutine compute_multistep_error(nstep,lmeanerrmssub, &
  meanerrmssub,timesub,ycf,outerycf,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the error introduced by the
!     multiple step approach
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  logical, intent(in) :: lmeanerrmssub
  double precision, intent(out) :: meanerrmssub
  double precision, intent(in) ::  timesub
  double precision, allocatable, intent(inout) ::  ycf(:,:)
  double precision, allocatable, intent(inout) :: outerycf(:,:)
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
  
  double precision, save :: meanerrforce
  
  
  double precision :: norma,normb,dtemp(1),myerr
  
  integer :: imiomax,ipoint,isub,mymax,mymaxmirr,jsub,jpoint
  double precision :: norm,Qt,xjpoint,yjpoint,zjpoint,dt
  double precision, dimension(3) :: versor,utang
  
  integer, save :: icounter=0
  
  meanerrmssub=0.d0
  
  if(.not.(lmeanerrmssub .and. lmscomputed))return
  if(msinpjet/=inpjet)return
  if(msnpjet/=npjet)return
  
  ycf(:,:)=0.d0
  outerycf(:,:)=0.d0
  coulservicearr(:,:)=0.d0
  
  
  select case(systype)
    case(1)
    
!     compute the total Coulomb forces
      call compute_coulomelec(nstep,timesub,ycf,yxx)
      
!     compute the inner Coulomb forces
      call compute_inner_coulomelec(timesub,coulservicearr,yxx)
        
      imiomax=(ncoulforce+1)
      call sum_world_darr(coulservicearr,imiomax)
      
!     compute the outer Coulomb forces as difference
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        outerycf(isub,1)=ycf(ipoint,1)-coulservicearr(ipoint,1)
      enddo
      
      coulservicearr(:,:)=0.d0
      
      call compute_outer_coulomelec(timesub,coulservicearr,yxx)
      
      
!     compute the error
      dtemp(1)=0.d0
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        dtemp(1)=dtemp(1)+ &
         dabs(coulservicearr(ipoint,1)-outerycf(isub,1))
      enddo
      call sum_world_darr(dtemp(1),1)
      meanerrmssub=dtemp(1)/dble(npjet-inpjet+1)
      
    case default
    
!     compute the total Coulomb forces
      call compute_coulomelec(nstep,timesub,ycf,yxx,yyy,yzz)
      
!     compute the inner Coulomb forces
      call compute_inner_coulomelec(timesub,coulservicearr,yxx,yyy,yzz)
        
      imiomax=(ncoulforce+1)*3
      call sum_world_darr(coulservicearr,imiomax)
      
!     compute the outer Coulomb forces as difference
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        outerycf(isub,1:3)=ycf(ipoint,1:3)-coulservicearr(ipoint,1:3)
      enddo
      
      coulservicearr(:,:)=0.d0
      
      call compute_outer_coulomelec(timesub,coulservicearr,yxx,yyy,yzz)
      
      
!     compute the error
      dtemp(1)=0.d0
      isub=0
      do ipoint=inpjet+idrank,npjet,mxrank
        if(jetfr(ipoint))cycle
        isub=isub+1
        utang(1)=coulservicearr(ipoint,1)-outerycf(isub,1)
        utang(2)=coulservicearr(ipoint,2)-outerycf(isub,2)
        utang(3)=coulservicearr(ipoint,3)-outerycf(isub,3)
        norm = modulvec(utang)
        dtemp(1)=dtemp(1) + norm
      enddo
      call sum_world_darr(dtemp(1),1)
      meanerrmssub=dtemp(1)/dble(npjet-inpjet+1)
    
  end select
      
      
  return
  
 end subroutine compute_multistep_error
 
 subroutine reallocate_neighlist()
 
!***********************************************************************
!     
!     JETSPIN subroutine for computing the Coulomb forces belonging to
!     the inner shell by the multiple step approach
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
  
  integer, save :: nservicelist1=0
  integer, save :: nservicelist2=0
  integer, save, allocatable, dimension(:,:) :: servicelist
 
  integer :: mychunk,newmaxneighlist
  integer :: newnservicelist1,newnservicelist2
  logical :: ldoreallocate
  
  call error(18)
  
  mychunk=ceiling(dble(mxnpjet)/dble(mxrank))
  newmaxneighlist=maxneighlist+incnpjet
  
  if(lmirror)then
    newnservicelist1=2*mychunk
    allocate(servicelist(newnservicelist1,maxneighlist))
    nservicelist1=newnservicelist1
    nservicelist2=newmaxneighlist
    servicelist(:,:)=0
    servicelist(1:nservicelist1,1:maxneighlist)= &
     neighlist(1:nservicelist1,1:maxneighlist)
    deallocate(neighlist)
    allocate(neighlist(nservicelist1,newmaxneighlist))
    neighlist(:,:)=0
    neighlist(1:nservicelist1,1:maxneighlist)= &
     servicelist(1:nservicelist1,1:maxneighlist)
  else
    newnservicelist1=mychunk
    allocate(servicelist(newnservicelist1,maxneighlist))
    nservicelist1=newnservicelist1
    nservicelist2=newmaxneighlist
    servicelist(:,:)=0
    servicelist(1:nservicelist1,1:maxneighlist)= &
     neighlist(1:nservicelist1,1:maxneighlist)
    deallocate(neighlist)
    allocate(neighlist(nservicelist1,newmaxneighlist))
    neighlist(:,:)=0
    neighlist(1:nservicelist1,1:maxneighlist)= &
     servicelist(1:nservicelist1,1:maxneighlist)
  endif
  
  maxneighlist=newmaxneighlist
  
  deallocate(servicelist)
  
  return
  
 end subroutine reallocate_neighlist
 
 subroutine list_test(newlist,timesub,yxx,yyy,yzz)
 
!***********************************************************************
!     
!     JETSPIN subroutine to test for updating of neighbour list
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification January 2016
!     
!***********************************************************************
 
  implicit none
      
  logical, intent(inout) :: newlist
  double precision, intent(in) :: timesub
  double precision, allocatable, intent(in) ::  yxx(:)
  double precision, allocatable, intent(in), optional ::  yyy(:)
  double precision, allocatable, intent(in), optional ::  yzz(:)
      
  logical, save :: newjob=.true.
  integer, save :: mychunkold=0
  integer :: ipoint,isub,jsub,moved,ibuff(1),mychunk
  double precision :: rmax,dr
      
  double precision, allocatable, save :: xold(:),yold(:),zold(:)
  double precision, allocatable, save :: xdif(:),ydif(:),zdif(:)
  
  if(.not.lmaxdispl)return
  
  mychunk=ceiling(dble(mxnpjet)/dble(mxrank))
  if(mychunk>mychunkold)newlist=.true.
  
  select case(systype)
    case(1)
      
      if(newjob .or. newlist)then
        
        if(newjob)then
          mychunkold=mychunk
          allocate (xold(mychunk))
          allocate (xdif(mychunk))
        else
          if(mychunk>mychunkold)then
            mychunkold=mychunk
            deallocate(xold)
            deallocate(xdif)
            allocate (xold(mychunk))
            allocate (xdif(mychunk))
          endif
        endif
        
        isub=0
        do ipoint=inpjet+idrank,npjet,mxrank
          isub=isub+1
          xold(isub)=yxx(ipoint)
        enddo
        newjob=.false.
        newlist=.true.
         
      else
        
!       calculate atomic shifts
        isub=0
        do ipoint=inpjet+idrank,npjet,mxrank
          isub=isub+1
          xdif(isub)=yxx(ipoint)-xold(isub)
        enddo
        
!       maximum displacement 
        rmax=(maxdispl/2.d0)**2.d0
        
!       test atomic displacements
        moved=0
        do jsub=1,isub
          dr=(xdif(jsub)**2.d0)
          if(dr>rmax)moved=moved+1
        enddo
        
        ibuff(1)=moved
        call sum_world_iarr(ibuff,1)
        moved=ibuff(1)
        
!       test for new list
        newlist=(moved>=2)
        
!       update stored positions
        if(newlist)then
          isub=0
          do ipoint=inpjet+idrank,npjet,mxrank
            isub=isub+1
            xold(isub)=yxx(ipoint)
          enddo
        endif
      endif
      
    case default
       
      if(newjob .or. newlist)then
        
        if(newjob)then
          mychunkold=mychunk
          allocate (xold(mychunk),yold(mychunk),zold(mychunk))
          allocate (xdif(mychunk),ydif(mychunk),zdif(mychunk))
        else
          if(mychunk>mychunkold)then
            mychunkold=mychunk
            deallocate(xold,yold,zold)
            deallocate(xdif,ydif,zdif)
            allocate (xold(mychunk),yold(mychunk),zold(mychunk))
            allocate (xdif(mychunk),ydif(mychunk),zdif(mychunk))
          endif
        endif
        
        isub=0
        do ipoint=inpjet+idrank,npjet,mxrank
          isub=isub+1
          xold(isub)=yxx(ipoint)
          yold(isub)=yyy(ipoint)
          zold(isub)=yzz(ipoint)
        enddo
        newjob=.false.
        newlist=.true.
         
      else
        
!       calculate atomic shifts
        isub=0
        do ipoint=inpjet+idrank,npjet,mxrank
          isub=isub+1
          xdif(isub)=yxx(ipoint)-xold(isub)
          ydif(isub)=yyy(ipoint)-yold(isub)
          zdif(isub)=yzz(ipoint)-zold(isub)
        enddo
        
!       maximum displacement 
        rmax=(maxdispl/2.d0)**2.d0
        
!       test atomic displacements
        moved=0
        do jsub=1,isub
          dr=(xdif(jsub)**2.d0+ydif(jsub)**2.d0+zdif(jsub)**2.d0)
          if(dr>rmax)moved=moved+1
        enddo
        
        ibuff(1)=moved
        call sum_world_iarr(ibuff,1)
        moved=ibuff(1)
        
!       test for new list
        newlist=(moved>=2)
        
!       update stored positions
        if(newlist)then
          isub=0
          do ipoint=inpjet+idrank,npjet,mxrank
            isub=isub+1
            xold(isub)=yxx(ipoint)
            yold(isub)=yyy(ipoint)
            zold(isub)=yzz(ipoint)
          enddo
        endif
      endif
  end select
  
       
  return
   
 end subroutine list_test
 
 end module coulomb_force_mod
 
