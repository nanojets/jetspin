 
 module dynamic_refinement_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which are dealing
!     the dynamic refinement of the nanojet
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!*********************************************************************** 
 
 use version_mod, only : mystart,myend,memyend,idrank,mxrank, &
                   sum_world_darr,or_world_larr,bcast_world_darr, &
                   sum_world_iarr
 use error_mod,   only : error,warning
 use utility_mod, only : Pi,buffservice,allocate_array_buffservice, &
                   resize_gaussian_history
 use nanojet_mod, only : resolution,inpjet,npjet,jetpt,jetxx,jetyy, &
                  jetzz,jetvx,jetvy,jetvz,jetst,jetms,jetch,jetvl, &
                  jetcr,jetve,linserted,linserting,systype,lengthpath, &
                  doallocate,doreorder,mxnpjet,incnpjet,lengthscale, &
                  ivolume,jetbd,massratio,imassa,jetfr,h, &
                  lenthresholdbead,ltagbeads,lbreakup,jetbr, &
                  lmultiplestep,lneighlistdo,jetfm,lmassavariable, &
                  lenprobmassa,icharge,jetce,levaporation
 use fit_mod,     only : jetptc,allocate_arrayspline,create_spline, &
                   allocate_array_jetptc,driver_fit_spline, &
                   looking_indexes_2,looking_indexes_4, &
                   cubic_interpolation,fit_akima,allocate_arrayakima, &
                   jetbdc,jetbrc,begin_akima_accelerator_data, &
                   end_akima_accelerator_data
 use support_functions_mod, only : compute_length_path, &
                   compute_crosssec
 use breaking_mod, only : clean_breakup
#ifdef _OPENACC
 use accelerator_mod, only : accelerator_device_state_is_current, &
                   accelerator_refinement_candidate, &
                   accelerator_update_host_capacity_state, &
                   accelerator_update_host_evaporation_state, &
                   accelerator_release_jet_capacity, &
                   accelerator_release_evaporation_capacity, &
                   accelerator_is_topology_enabled, &
                   accelerator_rebind_topology, &
                   accelerator_rebind_evaporation, &
                   accelerator_update_device_topology_state, &
                   accelerator_update_device_evaporation_state, &
                   accelerator_mark_device_state
#endif
 
 implicit none
 
 private
 
 integer, save :: nfitspline
 integer, save :: nfitlinear
 integer, save :: nfitting
 integer, save :: oldlowerbound
 integer, save :: oldlowerbuff
 integer, save :: newlowerbound
 integer, save :: newlowerbuff
 integer, save :: newupperbound
 
 logical, public, save :: lrefinement=.false.
 logical, public, save :: lrefinementthreshold=.false.
 logical, public, save :: lrefinementevery=.false.
 logical, public, save :: lrefinementstart=.false.
 logical, public, save :: llenthresholdbead=.false.
 
 logical, public, save :: lrefbeadstart=.false.
 
 integer, public, save :: irefbeadstart=0
 integer, public, save :: irefinementstart=0
 integer, public, save :: irefinementevery=0
 integer, public, save :: irefinementdone=0
 
 double precision, public, save :: refinementthreshold=0.d0
 
 double precision, public, save :: refbeadstartfit=0.d0
 
 double precision, allocatable, save :: buff(:)
 integer, save :: nbuff=0
 
 logical, save :: mydoallocate
 
 integer, save :: nmyindex=0
 integer, allocatable, save :: myindex(:,:)
 
 integer, save :: nkeepvol=0
 double precision, dimension(:), allocatable, save :: keepvol
 double precision, dimension(:), allocatable, save :: keepvolev
 
 integer, save :: inpjetbackup,npjetbackup
 integer, save :: nbackup=0
 integer, save :: njetlthr=0
 double precision, save :: timebackup
 logical, dimension(:), allocatable, save :: jetlthr
 double precision, dimension(:), allocatable, save :: jetxxbak, &
  jetyybak,jetzzbak,jetvxbak,jetvybak,jetvzbak,jetstbak,jetmsbak, &
  jetchbak,jetvlbak,jetvebak
 
 
 public :: driver_dynamic_refinement
 public :: set_refinement_threshold
 
 contains
 
 subroutine driver_dynamic_refinement(k,dorefinment)

!***********************************************************************
!     
!     JETSPIN subroutine for driving the dynamic refinement
!     of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************

implicit none
  
  integer, intent(in) :: k
  logical, intent(inout) :: dorefinment
  logical :: device_refinement

  device_refinement=.false.
#ifdef _OPENACC
  device_refinement=accelerator_device_state_is_current()
#endif
  
  if(.not.lrefinement)return
  if(k<irefinementstart)return
  
! check if the dynamic refinement is necessary
  call check_dynamic_refinement_akima(k,dorefinment)
  
  if(.not.dorefinment)return
  
! apply the dynamic refinement if requested
  call dynamic_refinement_akima(k,dorefinment,device_refinement)

#ifdef _OPENACC
! Only an accepted event reaches this point. Target-mesh construction and
! conservation leave the completed state on the host after the device Akima
! kernels; upload it once and resume persistent device integration.
  if(device_refinement)then
    if(accelerator_is_topology_enabled())then
      call accelerator_update_device_topology_state(npjet,jetxx,jetyy,jetzz, &
       jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
      if(levaporation) &
       call accelerator_update_device_evaporation_state(npjet,jetve,jetce)
    else
! Capacity growth replaced the host allocations after their old mappings
! were released. Bind the new addresses and copy the completed Akima
! state once before persistent integration resumes.
      call accelerator_rebind_topology(mxnpjet,jetxx,jetyy,jetzz,jetst, &
       jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
      if(levaporation) &
       call accelerator_rebind_evaporation(mxnpjet,jetve,jetce)
    endif
    call accelerator_mark_device_state(.true.)
  endif
#endif
  
  return
  
 end subroutine driver_dynamic_refinement
 
 subroutine set_refinement_threshold()
 
!***********************************************************************
!     
!     JETSPIN subroutine for setting the threshold of the dynamics 
!     refinement of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  
  if(lrefinement)then
    if(.not.lrefinementthreshold)then
      refinementthreshold=10.d0*resolution
      refbeadstartfit=10.d0*resolution
      lrefbeadstart=.true.
      lrefinementthreshold=.true.
      call warning(66,refinementthreshold)
    else
      if(refinementthreshold<(2.d0*resolution))call error(16)
      if(refinementthreshold<5.d0*resolution)then
        refinementthreshold=5.d0*resolution
        refbeadstartfit=refinementthreshold
        call warning(86,refinementthreshold)
      endif
    endif
    if(.not. llenthresholdbead)then
      llenthresholdbead=.true.
      lenthresholdbead=5.d0*resolution
      call warning(81,lenthresholdbead)
    else
      if(lenthresholdbead<resolution)then
        lenthresholdbead=5.d0*resolution
        call warning(85,lenthresholdbead)
      endif
    endif
  endif
  
  return
  
 end subroutine set_refinement_threshold
 
 subroutine check_dynamic_refinement_akima(nstep,dorefinment)
 
!***********************************************************************
!     
!     JETSPIN subroutine for checking the dynamic refinement
!     of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  logical, intent(out) :: dorefinment
 
  integer, save :: icounter=0
  integer :: device_refinement_start,device_nfitting
  double precision :: tempmod0,device_lengthpath,device_nozzle_correction
  logical :: ldorefinement,device_candidate
  
  dorefinment=.false.
  
  if(.not.lrefinement)return
  if(linserting .and. (npjet-inpjet+1)<10)return
  
  icounter=icounter+1
  if(icounter<irefinementevery)return

#ifdef _OPENACC
  if(accelerator_device_state_is_current())then
    call accelerator_refinement_candidate(inpjet,npjet,systype,linserting, &
     linserted,refbeadstartfit,jetxx,jetyy,jetzz,device_candidate, &
     device_refinement_start,device_lengthpath,device_nozzle_correction)
    if(.not.device_candidate)return

! Apply the inexpensive part of the historical CPU acceptance test using
! only device reductions.  A reduction-order difference can move a marginal
! event by one timestep, but should not trigger repeated full-state downloads.
    if(linserting)then
      device_nfitting=nint((device_lengthpath-device_nozzle_correction)/ &
       refinementthreshold)
      if(linserted)then
        device_nfitting=device_nfitting+2
      else
        device_nfitting=device_nfitting+3
      endif
    else
      device_nfitting=nint(device_lengthpath/refinementthreshold)
    endif
    if(device_nfitting<=npjet-inpjet+2)return

! The threshold scan found a possible event. Synchronize the complete state
! once for host target-mesh preparation; coefficient construction and
! interpolation can then execute through the device path.
    call accelerator_update_host_capacity_state(npjet,jetxx,jetyy,jetzz, &
     jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr,nstep)
    if(levaporation) &
     call accelerator_update_host_evaporation_state(npjet,jetve,jetce)
  endif
#endif

  call test_length_threshold_akima(nstep,jetxx,jetyy,jetzz, &
   refinementthreshold,ldorefinement,jetpt,lengthpath)
  
  if(.not.ldorefinement)return
  
  if(linserting)then
    if(.not.linserted)then
      select case(systype)
      case(1)
        tempmod0 = jetxx(npjet-2)-jetxx(npjet)
      case default
        tempmod0 = dsqrt((jetxx(npjet-2)-jetxx(npjet))**2.d0+ &
         (jetyy(npjet-2)-jetyy(npjet))**2.d0+ &
         (jetzz(npjet-2)-jetzz(npjet))**2.d0)
      end select
      nfitspline=nint((lengthpath-tempmod0)/(refinementthreshold))
      nfitting=nfitspline+3
      if(nfitting<=nint(dble(npjet-inpjet+1)+1.d0))return
    else
      select case(systype)
      case(1)
        tempmod0 = jetxx(npjet-1)-jetxx(npjet)
      case default
        tempmod0 = dsqrt((jetxx(npjet-1)-jetxx(npjet))**2.d0+ &
         (jetyy(npjet-1)-jetyy(npjet))**2.d0+ &
         (jetzz(npjet-1)-jetzz(npjet))**2.d0)
      end select
      nfitspline=nint((lengthpath-tempmod0)/(refinementthreshold))
      nfitting=nfitspline+2
      if(nfitting<=nint(dble(npjet-inpjet+1)+1.d0))return
    endif
  else
    nfitspline=nint((lengthpath)/(refinementthreshold))
    nfitting=nfitspline
    if(nfitting<=nint(dble(npjet-inpjet+1)+1.d0))return
  endif
  
  dorefinment=.true.
  icounter=0
  irefinementdone=irefinementdone+1
  
  return
  
 end subroutine check_dynamic_refinement_akima
 
 subroutine dynamic_refinement_akima(nstep,dorefinment,device_refinement)
 
!***********************************************************************
!     
!     JETSPIN subroutine for performing the dynamic refinement
!     of the nanofiber
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: nstep
  logical, intent(in) :: dorefinment
  logical, intent(in) :: device_refinement
  
  integer :: i,j,jptinit,jptend,ipoint,totjptend,nmassbd,irefbeadstop
  integer :: totjptend0,istart
  integer :: oldactive,newactive,oldanchorcount,newanchorcount
  integer :: oldcapacity
  integer, parameter :: strategytype=0
  integer, allocatable :: massbd(:,:),massbdpoint(:)
  double precision, allocatable :: massbddist(:),massbdgap(:)
  double precision, allocatable :: anchorpos(:,:),anchorvel(:,:)
  double precision, allocatable :: anchorstress(:),anchorradius(:,:)
  double precision :: tempmod0,tempmod1,tempmod2,voltot,newvoltot
  double precision :: voltotev,newvoltotev
  double precision :: anchorerror,volerror,voleverror
  double precision :: anchorvelerror,anchorstresserror
  double precision :: anchorradiuserror,anchorevradiuserror
  double precision :: minsegmentlength,minradius,minvolume
  integer :: minvolumeat
  double precision :: oldmass,newmass,oldcharge,newcharge
  double precision :: masserror,chargeerror
  logical :: orderedmesh,device_capacity_rebind
  logical :: execute_device_reconstruction

  integer, save :: icounter=0
  
  if(.not.dorefinment)return
  
  mydoallocate=.false.
      
  irefbeadstop=-1
  do ipoint=inpjet,irefbeadstart
    if(jetxx(ipoint)>=h)irefbeadstop=max(irefbeadstop,ipoint)
  enddo
    
  if(irefbeadstop==-1)then
    istart=inpjet
  else
    istart=irefbeadstop
  endif
  
  nmassbd=0
  do ipoint=istart,irefbeadstart
    if(jetbd(ipoint))nmassbd=nmassbd+1
  enddo
  
  allocate(massbd(0:nmassbd,2),massbdpoint(0:nmassbd))
  allocate(massbddist(0:nmassbd),massbdgap(0:nmassbd))
    
  i=0
  massbd(0,1)=istart
  do ipoint=massbd(0,1),irefbeadstart
    if(jetbd(ipoint))then
      i=i+1
      massbd(i-1,2)=ipoint
      massbd(i,1)=ipoint
    endif
  enddo
  massbd(i,2)=irefbeadstart
  
  nfitspline=0
  do i=0,nmassbd
    tempmod0=(jetpt(massbd(i,2))-jetpt(massbd(i,1)))
    massbddist(i)=tempmod0
    massbdpoint(i)=ceiling((tempmod0*lengthpath)/(refinementthreshold))
    massbdgap(i)=(tempmod0/dble(massbdpoint(i)))
    nfitspline=massbdpoint(i)+nfitspline
  enddo
  
  if(irefbeadstop<inpjet)then 
    nfitspline=nfitspline+(irefbeadstop-inpjet+1)+2
    nfitting=nfitspline+(npjet-irefbeadstart)
  else
    if(irefbeadstop==inpjet)then
      nfitspline=nfitspline+(irefbeadstop-inpjet+1)
      nfitting=nfitspline+(npjet-irefbeadstart)
    else
      nfitspline=nfitspline+(irefbeadstop-inpjet+1)
      nfitting=nfitspline+(npjet-irefbeadstart)
    endif
  endif
  
  if(nfitting<=(npjet-inpjet+1))return

  oldactive=npjet-inpjet
! The interpolation algorithm needs only anchors inside the fitted segment,
! hence nmassbd above.  The invariant check must instead follow every active
! anchor: anchors in the copied prefix and untouched tail are just as
! important during repeated refinement events.
  oldanchorcount=count(jetbd(inpjet:npjet))
  if(oldanchorcount>0)then
    allocate(anchorpos(3,oldanchorcount))
    allocate(anchorvel(3,oldanchorcount))
    allocate(anchorstress(oldanchorcount),anchorradius(2,oldanchorcount))
    j=0
    do ipoint=inpjet,npjet
      if(.not.jetbd(ipoint))cycle
      j=j+1
      anchorpos(1,j)=jetxx(ipoint)
      anchorpos(2,j)=jetyy(ipoint)
      anchorpos(3,j)=jetzz(ipoint)
      anchorvel(1,j)=jetvx(ipoint)
      anchorvel(2,j)=jetvy(ipoint)
      anchorvel(3,j)=jetvz(ipoint)
      anchorstress(j)=jetst(ipoint)
      anchorradius(1,j)=0.d0
      anchorradius(2,j)=0.d0
    enddo
  endif
  
  call allocate_arrayakima()

  oldmass=sum(jetms(inpjet:irefbeadstart))
  oldcharge=sum(jetch(inpjet:irefbeadstart))
  
  call convert_to_density(jetms,jetch,jetvl)
  
  voltot=0.d0
  do i=inpjet,irefbeadstart
    voltot=voltot+jetvl(i)
  enddo
  call allocate_array_keepvol(npjet-irefbeadstart)
  i=0
  do ipoint=irefbeadstart+1,npjet
    i=i+1
    keepvol(i)=jetvl(ipoint)
  enddo
  
  call compute_crosssec(jetxx,jetyy,jetzz,jetvl,jetcr)
  j=0
  do ipoint=inpjet,npjet
    if(.not.jetbd(ipoint))cycle
    j=j+1
    anchorradius(1,j)=jetcr(ipoint)
  enddo
  
  if(levaporation)then
  
    voltotev=0.d0
    do i=inpjet,irefbeadstart
      voltotev=voltotev+jetve(i)
    enddo
    
    i=0
    do ipoint=irefbeadstart+1,npjet
      i=i+1
      keepvolev(i)=jetve(ipoint)
    enddo
  
    call compute_crosssec(jetxx,jetyy,jetzz,jetve,jetce)
    j=0
    do ipoint=inpjet,npjet
      if(.not.jetbd(ipoint))cycle
      j=j+1
      anchorradius(2,j)=jetce(ipoint)
    enddo
  
  endif
  
  oldcapacity=mxnpjet
  call define_akima_bounds()
  if(mydoallocate .and. idrank==0)then
    write(6,'(a,i0,a,i0)')'Dynamic refinement capacity: old=', &
     oldcapacity,' new=',mxnpjet
  endif
  device_capacity_rebind=.false.
#ifdef _OPENACC
  device_capacity_rebind=mydoallocate .and. &
   accelerator_device_state_is_current()
  if(device_capacity_rebind)then
! The active state was downloaded by the acceptance check. Detach the old
! capacity before fit_jet_akima deallocates any persistently mapped array.
    if(levaporation) &
     call accelerator_release_evaporation_capacity(oldcapacity,jetve,jetce)
    call accelerator_release_jet_capacity(oldcapacity,jetxx,jetyy,jetzz, &
     jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr)
    if(idrank==0)write(6,'(a,i0,a,i0)') &
     'OpenACC refinement capacity rebind: old=',oldcapacity, &
     ' new=',mxnpjet
  endif
#endif
! Gaussian values belonging to existing bead indices are preserved while
! any new capacity slots are generated once on the host. In OpenACC builds
! the resized history is detached and rebound inside this routine.
  if(mydoallocate) &
   call resize_gaussian_history(mxnpjet,device_capacity_rebind)
  call allocate_array_jetptc()
  j=0
  jptinit=newlowerbuff+1
  jptend=jptinit+nfitspline
  
  jetbdc(:)=.false.
  if(lbreakup)jetbrc(:)=.false.
  
  j=jptinit
  jetptc(j)=0.d0
  jetbdc(j)=jetbd(inpjet)
  if(lbreakup)jetbrc(j)=jetbr(inpjet)
  
  do ipoint=inpjet+1,irefbeadstop
    j=j+1
    jetptc(j)=jetpt(ipoint)
    jetbdc(j)=jetbd(ipoint)
    if(lbreakup)jetbrc(j)=jetbr(ipoint)
  enddo
  do i=0,nmassbd
    tempmod0=massbdgap(i)
    do ipoint=1,massbdpoint(i)
      j=j+1
      jetptc(j)=jetptc(j-1)+tempmod0
      if(ipoint==massbdpoint(i))then
        if(i/=nmassbd)then
          jetbdc(j)=.true.
          jetptc(j)=jetpt(massbd(i,2))
          if(lbreakup)jetbrc(j)=jetbr(massbd(i,2))
        else
          jetptc(j)=jetpt(irefbeadstart)
          if(lbreakup)jetbrc(j)=jetbr(irefbeadstart)
        endif         
      endif
    enddo
  enddo
  totjptend0=j
  do ipoint=irefbeadstart+1,npjet
    j=j+1
    jetptc(j) = jetpt(ipoint)
    jetbdc(j)=jetbd(ipoint)
    if(lbreakup)jetbrc(j)=jetbr(ipoint) 
  enddo
  totjptend=j
    
  
  if(mydoallocate)then
    deallocate(jetbd)
    allocate(jetbd(0:mxnpjet))
  endif
  jetbd(0:mxnpjet)=.false.
  do ipoint=jptinit,totjptend
    jetbd(ipoint)=jetbdc(ipoint)
  enddo
  
  if(lbreakup)then
    if(mydoallocate)then
      deallocate(jetbr)
      allocate(jetbr(0:mxnpjet))
    endif
    jetbr(0:mxnpjet)=.false.
    do ipoint=jptinit,totjptend
      jetbr(ipoint)=jetbrc(ipoint)
    enddo
  endif
  
  deallocate(massbd,massbdpoint)
  deallocate(massbddist,massbdgap)
  
  call fit_jet_akima(jptinit,totjptend,device_refinement)

! The bead volume/evaporation-volume reconstruction from the freshly
! interpolated cross-section radius, their reference-volume conservation
! rescale, and the density-to-quantity conversion of mass and charge are
! all independent per-bead updates or single reductions over the new
! mesh, unlike the data-dependent mass-boundary/anchor bookkeeping above.
! Route them to the same validated device path used by Akima whenever it
! is active; the developer-gated variable-mass branch of
! convert_from_density is intentionally excluded and always falls back
! to the host.
  execute_device_reconstruction=.false.
#ifdef _OPENACC
  execute_device_reconstruction=device_refinement .and. mxrank==1 &
   .and. .not.lmassavariable
#ifdef JETSPIN_DEV_HOST_AKIMA
  execute_device_reconstruction=.false.
#endif
#endif

#ifdef _OPENACC
  if(execute_device_reconstruction)then
#ifdef JETSPIN_COMPARE_REFINEMENT_ASSEMBLY
! This diagnostic evaluates the host reference from the same
! pre-reconstruction inputs and leaves the device result as the final,
! authoritative state; it replaces (not supplements) the plain device
! call below so the device kernel still runs exactly once per event.
    call compare_refinement_assembly_device(jptinit,jptend,totjptend, &
     voltot,voltotev,levaporation)
#else
    call accelerator_reconstruct_refinement_state(jptinit,jptend, &
     totjptend,keepvol,keepvolev,voltot,voltotev,levaporation)
#endif
  else
#endif
    call reconstruct_refinement_state_host(jptinit,jptend,totjptend, &
     voltot,voltotev,levaporation)
#ifdef _OPENACC
  endif
#endif

! Check the remeshing invariants while both the old anchor coordinates and
! the completed target state are available.  These diagnostics are emitted
! only at an accepted refinement event, not at every timestep.
  newactive=npjet-inpjet
  newanchorcount=0
  anchorerror=0.d0
  anchorvelerror=0.d0
  anchorstresserror=0.d0
  anchorradiuserror=0.d0
  anchorevradiuserror=0.d0
  do ipoint=jptinit,totjptend
    if(jetbd(ipoint))then
      newanchorcount=newanchorcount+1
      if(newanchorcount<=oldanchorcount)then
        anchorerror=max(anchorerror, &
         dsqrt((jetxx(ipoint)-anchorpos(1,newanchorcount))**2.d0+ &
               (jetyy(ipoint)-anchorpos(2,newanchorcount))**2.d0+ &
               (jetzz(ipoint)-anchorpos(3,newanchorcount))**2.d0))
        anchorvelerror=max(anchorvelerror, &
         dsqrt((jetvx(ipoint)-anchorvel(1,newanchorcount))**2.d0+ &
               (jetvy(ipoint)-anchorvel(2,newanchorcount))**2.d0+ &
               (jetvz(ipoint)-anchorvel(3,newanchorcount))**2.d0))
        anchorstresserror=max(anchorstresserror, &
         dabs(jetst(ipoint)-anchorstress(newanchorcount)))
        anchorradiuserror=max(anchorradiuserror, &
         dabs(jetcr(ipoint)-anchorradius(1,newanchorcount)))
        if(levaporation)then
          anchorevradiuserror=max(anchorevradiuserror, &
           dabs(jetce(ipoint)-anchorradius(2,newanchorcount)))
        endif
      endif
    endif
  enddo

  newvoltot=sum(jetvl(jptinit:jptend))
  volerror=dabs(newvoltot-voltot)/max(dabs(voltot),tiny(1.d0))
  newmass=sum(jetms(jptinit:jptend))
  newcharge=sum(jetch(jptinit:jptend))
  masserror=dabs(newmass-oldmass)/max(dabs(oldmass),tiny(1.d0))
  chargeerror=dabs(newcharge-oldcharge)/max(dabs(oldcharge),tiny(1.d0))
  voleverror=0.d0
  if(levaporation)then
    newvoltotev=sum(jetve(jptinit:jptend))
    voleverror=dabs(newvoltotev-voltotev)/ &
     max(dabs(voltotev),tiny(1.d0))
  endif

  orderedmesh=.true.
  minsegmentlength=huge(1.d0)
  do ipoint=inpjet,npjet-1
    if(jetpt(ipoint+1)<=jetpt(ipoint))orderedmesh=.false.
    minsegmentlength=min(minsegmentlength,jetpt(ipoint+1)-jetpt(ipoint))
  enddo
  minradius=minval(jetcr(inpjet:npjet))
  minvolume=minval(jetvl(inpjet:npjet))
  minvolumeat=inpjet-1+minloc(jetvl(inpjet:npjet),1)

  if(idrank==0)then
    write(6,'(a,i0,4(a,i0))')'Dynamic refinement event: step=',nstep, &
     ' active_before=',oldactive,' active_after=',newactive, &
     ' anchors_before=',oldanchorcount,' anchors_after=',newanchorcount
    write(6,'(a,3(a,es12.4),2(a,i0))')'Dynamic refinement geometry check:', &
     ' min_segment_length_cm=',minsegmentlength*lengthscale, &
     ' min_radius_cm=',minradius*lengthscale, &
     ' min_reference_volume_cm3=',minvolume, &
     ' min_volume_bead=',minvolumeat,' jptinit=',jptinit
    write(6,'(a,3(a,i0))')'Dynamic refinement geometry range:', &
     ' jptend=',jptend,' totjptend=',totjptend,' npjet=',npjet
    write(6,'(a,3(a,es12.4),a,l1)')'Dynamic refinement invariants:', &
     ' anchor_position_max_displacement_cm=',anchorerror*lengthscale, &
     ' reference_volume_relative_difference=',volerror, &
     ' evaporation_volume_relative_difference=',voleverror, &
     ' ordered_path=',orderedmesh
    write(6,'(a,4(a,es12.4))')'Dynamic refinement anchor fields:', &
     ' velocity_max_difference_internal=',anchorvelerror, &
     ' stress_max_difference_internal=',anchorstresserror, &
     ' radius_max_difference_cm=',anchorradiuserror*lengthscale, &
     ' evaporation_radius_max_difference_cm=', &
     anchorevradiuserror*lengthscale
    write(6,'(a,2(a,es12.4))')'Dynamic refinement conserved amounts:', &
     ' mass_relative_difference=',masserror, &
     ' charge_relative_difference=',chargeerror
  endif

  if(allocated(anchorpos))deallocate(anchorpos)
  if(allocated(anchorvel))deallocate(anchorvel)
  if(allocated(anchorstress))deallocate(anchorstress)
  if(allocated(anchorradius))deallocate(anchorradius)
  
  doallocate=(doallocate .or. mydoallocate)
  if(lmultiplestep)lneighlistdo=.true.
  if(lbreakup)call clean_breakup(jptinit,jptend,totjptend)
  
  return
  
 end subroutine dynamic_refinement_akima
 
 subroutine test_length_threshold_akima(nstep,yxx,yyy,yzz,dthreshold, &
  lthreshold,ypt,lengthpathsub)
  
!***********************************************************************
!     
!     JETSPIN subroutine for testing if the mutual distance between
!     any two nanofiber beads is beyond a given threshold
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification December 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: nstep
  double precision, allocatable, dimension(:), intent(in) :: yxx
  double precision, allocatable, dimension(:), intent(in) :: yyy
  double precision, allocatable, dimension(:), intent(in) :: yzz
  double precision, intent(in) :: dthreshold
  
  logical, intent(inout) :: lthreshold
  double precision, allocatable, dimension(:), intent(inout) :: ypt
  double precision, intent(inout) :: lengthpathsub
  
  integer :: ipoint
  double precision :: tempmod0,tempmod,tempmod1(1)
  logical, dimension(1) :: logicaltemp
  
  call compute_length_path(yxx,yyy,yzz, &
   ypt,lengthpathsub)
  
  logicaltemp(1)=.false.
  
  irefbeadstart=0
  if(linserting)then
    if(.not.linserted)then
      do ipoint=inpjet,npjet-3
        tempmod0 = (ypt(ipoint+1)-ypt(ipoint))*lengthpathsub
        if(tempmod0 > refbeadstartfit )then
          logicaltemp(1)=.true.
          irefbeadstart=ipoint
        endif
      enddo
    else
      do ipoint=inpjet,npjet-2
        tempmod0 = (ypt(ipoint+1)-ypt(ipoint))*lengthpathsub
        if(tempmod0 > refbeadstartfit )then
          logicaltemp(1)=.true.
          irefbeadstart=ipoint
        endif
      enddo 
    endif
  else
    do ipoint=inpjet,npjet-1
      tempmod0 = (ypt(ipoint+1)-ypt(ipoint))*lengthpathsub
      if(tempmod0 > refbeadstartfit )then
        logicaltemp(1)=.true.
        irefbeadstart=ipoint
      endif
    enddo
  endif
  
  lthreshold=logicaltemp(1)
  
  return
  
 end subroutine test_length_threshold_akima
 
 subroutine define_akima_bounds() 
 
!***********************************************************************
!     
!     JETSPIN subroutine for defining the new bounds of arrays which 
!     describe the nanojet if the spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  integer :: ncutoffsub,growthincrement,growthiostat
  character(len=32) :: growthenv
  
  mydoallocate=.false.
  ncutoffsub=1
  growthincrement=incnpjet
! Production refinement continues to grow by incnpjet entries.  Test 22 can
! reduce only this event-time increment to exercise several release/rebind
! cycles in a short run without changing the physical input.
  growthenv=''
  call get_environment_variable('JETSPIN_REFINEMENT_GROWTH_INCREMENT', &
   growthenv)
  if(len_trim(growthenv)>0)then
    read(growthenv,*,iostat=growthiostat)growthincrement
    if(growthiostat/=0 .or. growthincrement<1)growthincrement=incnpjet
  endif
  
  if(inpjet==0)then
    oldlowerbound=0
    oldlowerbuff=inpjet-1
    
    newlowerbound=0
    newlowerbuff=oldlowerbuff-oldlowerbound
    
    mydoallocate=.false.
    newupperbound=nfitting
    
  elseif(inpjet>0)then
    oldlowerbound=max(0,inpjet-ncutoffsub)
    oldlowerbuff=inpjet-1
    
    newlowerbound=0
    newlowerbuff=oldlowerbuff-oldlowerbound
    
    newupperbound=nfitting+newlowerbuff+1
    
  endif

  if(newupperbound>mxnpjet)then
    mxnpjet=newupperbound+growthincrement
    mydoallocate=.true.
  endif
  
  return
  
 end subroutine define_akima_bounds
 
 subroutine fit_jet_akima(jptinit,jptend,device_akima)
 
!***********************************************************************
!     
!     JETSPIN subroutine for managing the cubic spline interpolation
!     if any two nanofiber beads are beyond a given threshold
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: jptinit,jptend
  logical, intent(in) :: device_akima
  
  
  integer :: i,ipoint
  double precision :: tempmod0
  logical :: execute_device_akima

  execute_device_akima=device_akima
#ifndef _OPENACC
  execute_device_akima=.false.
#else
#ifdef JETSPIN_DEV_HOST_AKIMA
  execute_device_akima=.false.
#endif
#endif
  
  call allocate_array_buffservice(newupperbound)
  call begin_akima_accelerator_data(jetpt,jetptc,execute_device_akima)
  
  doreorder=.true.
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetxx(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetxx,jetptc,buffservice, &
   execute_device_akima,'x')
  if(mydoallocate)then
    deallocate(jetxx)
    allocate(jetxx(0:mxnpjet))
  endif
  jetxx(:)=0.d0
  jetxx(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetyy(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetyy,jetptc,buffservice, &
   execute_device_akima,'y')
  if(mydoallocate)then
    deallocate(jetyy)
    allocate(jetyy(0:mxnpjet))
  endif
  jetyy(:)=0.d0
  jetyy(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetzz(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetzz,jetptc,buffservice, &
   execute_device_akima,'z')
  if(mydoallocate)then
    deallocate(jetzz)
    allocate(jetzz(0:mxnpjet))
  endif
  jetzz(:)=0.d0
  jetzz(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetvx(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetvx,jetptc,buffservice, &
   execute_device_akima,'vx')
  if(mydoallocate)then
    deallocate(jetvx)
    allocate(jetvx(0:mxnpjet))
  endif
  jetvx(:)=0.d0
  jetvx(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
    
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetvy(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetvy,jetptc,buffservice, &
   execute_device_akima,'vy')
  if(mydoallocate)then
    deallocate(jetvy)
    allocate(jetvy(0:mxnpjet))
  endif
  jetvy(:)=0.d0
  jetvy(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetvz(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetvz,jetptc,buffservice, &
   execute_device_akima,'vz')
  if(mydoallocate)then
    deallocate(jetvz)
    allocate(jetvz(0:mxnpjet))
  endif
  jetvz(:)=0.d0
  jetvz(newlowerbound:newupperbound)= &
   buffservice(newlowerbound:newupperbound)
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetst(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetst,jetptc,buffservice, &
   execute_device_akima,'stress')
  if(mydoallocate)then
    deallocate(jetst)
    allocate(jetst(0:mxnpjet))
  endif
  jetst(:)=0.d0
  jetst(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetms(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetms,jetptc,buffservice, &
   execute_device_akima,'mass_density')
  if(mydoallocate)then
    deallocate(jetms)
    allocate(jetms(0:mxnpjet))
  endif
  jetms(:)=0.d0
  jetms(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetch(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetch,jetptc,buffservice, &
   execute_device_akima,'charge_density')
  if(mydoallocate)then
    deallocate(jetch)
    allocate(jetch(0:mxnpjet))
  endif
  jetch(:)=0.d0
  jetch(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  buffservice(:)=0.d0
  if(idrank==0)then
    buffservice(newlowerbound:newlowerbuff)= &
     jetcr(oldlowerbound:oldlowerbuff)
  endif
  call fit_akima(jptinit,jptend,jetpt,jetcr,jetptc,buffservice, &
   execute_device_akima,'radius')
  if(mydoallocate)then
    deallocate(jetcr)
    allocate(jetcr(0:mxnpjet))
  endif
  jetcr(:)=0.d0
  jetcr(newlowerbound:newupperbound)= &
   dabs(buffservice(newlowerbound:newupperbound))
  
  if(levaporation)then
    buffservice(:)=0.d0
    if(idrank==0)then
      buffservice(newlowerbound:newlowerbuff)= &
       jetce(oldlowerbound:oldlowerbuff)
    endif
    call fit_akima(jptinit,jptend,jetpt,jetce,jetptc,buffservice, &
     execute_device_akima,'evap_radius')
    if(mydoallocate)then
      deallocate(jetce)
      allocate(jetce(0:mxnpjet))
    endif
    jetce(:)=0.d0
    jetce(newlowerbound:newupperbound)= &
     dabs(buffservice(newlowerbound:newupperbound))
  endif
  call end_akima_accelerator_data(jetpt,jetptc)
  
  if(mydoallocate)then
    deallocate(jetvl)
    allocate(jetvl(0:mxnpjet))
  endif
  jetvl(:)=0.d0
  
  if(levaporation)then
    if(mydoallocate)then
      deallocate(jetve)
      allocate(jetve(0:mxnpjet))
    endif
    jetve(:)=0.d0
  endif
  
  if(mydoallocate)then
    deallocate(jetpt)
    allocate(jetpt(0:mxnpjet))
  endif
  jetpt(:)=0.d0
  
  
  inpjet=jptinit
  npjet=jptend
  
  do ipoint=inpjet,npjet
    jetpt(ipoint)=jetptc(ipoint)
  enddo
    
  if(mydoallocate)then
    deallocate(jetfr)
    allocate(jetfr(0:mxnpjet))
    if(lmultiplestep)then
      deallocate(jetfm)
      allocate(jetfm(0:mxnpjet))
    endif
  endif
  jetfr(:)=.false.
  do ipoint=inpjet,npjet
    if(jetxx(ipoint)>=h)then
      jetfr(ipoint)=.true.
      jetxx(ipoint)=h
    endif
  enddo
  jetxx(npjet)=0.d0
  
 ! call compute_length_path(jetxx,jetyy,jetzz,jetpt,lengthpath)
  
  return
  
 end subroutine fit_jet_akima
 
 subroutine convert_to_density(jms,jch,jvl)

!***********************************************************************
!     
!     JETSPIN subroutine for converting the mass and charge data
!     to their respective densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) :: jms,jch
  double precision, allocatable, dimension (:), intent(inout) :: jvl
  
  integer :: ipoint
  double precision :: extramass,extravol
  
  
  if(lmassavariable)then
    extravol=(4.d0/3.d0)*Pi*lenprobmassa**3.d0
    extramass=massratio*extravol
    do ipoint=inpjet,npjet
      if(jetbd(ipoint))then
        jvl(ipoint)=jvl(ipoint)-extravol
        jms(ipoint)=imassa
        jch(ipoint)=icharge
      else
        jms(ipoint)=imassa
        jch(ipoint)=icharge
      endif
    enddo
  else
    jms(:)=imassa
    jch(:)=icharge
  endif
  
  
  return
  
 end subroutine convert_to_density
 
 subroutine convert_from_density(jms,jch,jvl)

!***********************************************************************
!     
!     JETSPIN subroutine for converting the mass and charge data
!     from their respective densities
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, allocatable, dimension (:), intent(inout) :: jms,jch
  double precision, allocatable, dimension (:), intent(inout) :: jvl
  
  integer :: ipoint
  double precision :: extramass,extravol
  
  if(lmassavariable)then
    extravol=(4.d0/3.d0)*Pi*lenprobmassa**3.d0
    extramass=massratio*extravol
    do ipoint=inpjet,npjet
      if(jetbd(ipoint))then
        jms(ipoint)=imassa*jvl(ipoint)+extramass
        jch(ipoint)=icharge*jvl(ipoint)
        jvl(ipoint)=jvl(ipoint)+extravol
      else
        jms(ipoint)=imassa*jvl(ipoint)
        jch(ipoint)=icharge*jvl(ipoint)
      endif
    enddo
  else
    jms(:)=jms(:)*jvl(:)
    jch(:)=jch(:)*jvl(:)
  endif
  
  return
  
 end subroutine convert_from_density

 subroutine reconstruct_refinement_state_host(jptinit,jptend,totjptend, &
  voltotsub,voltotevsub,do_evaporation)

!***********************************************************************
!
!     JETSPIN subroutine for reconstructing, on the host, the bead
!     volume and evaporation-volume of an accepted dynamic-refinement
!     event from the Akima-interpolated cross-section radius, applying
!     their reference-volume conservation rescale, and converting the
!     interpolated mass/charge densities back to per-bead quantities.
!     This is the trusted reference path, used directly when the device
!     path is unavailable or excluded, and as the oracle compared
!     against the OpenACC path under JETSPIN_COMPARE_REFINEMENT_ASSEMBLY.
!
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification August 2026
!
!***********************************************************************

  implicit none

  integer, intent(in) :: jptinit,jptend,totjptend
  double precision, intent(in) :: voltotsub,voltotevsub
  logical, intent(in) :: do_evaporation

  integer :: i,ipoint
  double precision :: tempmod0,newvoltot,newvoltotev

! When the fitted segment reaches the true jet endpoint (no preserved tail
! beyond it, jptend==totjptend), bead totjptend has no following point to
! define a forward segment length. This case never arose while every
! validated refinement test kept an un-refined tail beyond the fitted
! region; it does arise once the active jet is short enough that the whole
! mesh lies inside the refinement zone (e.g. a jet still growing from a
! single nozzle bead). Use the preceding segment length there instead of
! reading past the last active bead; fall back to the elementary
! resolution length in the degenerate single-bead case.
  do i=jptinit,jptend
    if(i==totjptend)then
      if(i>jptinit)then
        tempmod0=lengthpath*(jetptc(i)-jetptc(i-1))
      else
        tempmod0=resolution
      endif
    else
      tempmod0=lengthpath*(jetptc(i+1)-jetptc(i))
    endif
    jetvl(i)=tempmod0*Pi*(jetcr(i))**2.d0
  enddo

  i=0
  do ipoint=jptend+1,totjptend
    i=i+1
    jetvl(ipoint) = keepvol(i)
  enddo

  newvoltot=0.d0
  do i=jptinit,jptend
    newvoltot=newvoltot+jetvl(i)
  enddo

  jetvl(jptinit:jptend)=jetvl(jptinit:jptend)*voltotsub/newvoltot

  call convert_from_density(jetms,jetch,jetvl)

  if(do_evaporation)then

    do i=jptinit,jptend
      if(i==totjptend)then
        if(i>jptinit)then
          tempmod0=lengthpath*(jetptc(i)-jetptc(i-1))
        else
          tempmod0=resolution
        endif
      else
        tempmod0=lengthpath*(jetptc(i+1)-jetptc(i))
      endif
      jetve(i)=tempmod0*Pi*(jetce(i))**2.d0
    enddo

    i=0
    do ipoint=jptend+1,totjptend
      i=i+1
      jetve(ipoint) = keepvolev(i)
    enddo

    newvoltotev=0.d0
    do i=jptinit,jptend
      newvoltotev=newvoltotev+jetve(i)
    enddo

    jetve(jptinit:jptend)=jetve(jptinit:jptend)*voltotevsub/newvoltotev

  endif

  return

 end subroutine reconstruct_refinement_state_host

#ifdef _OPENACC
 subroutine accelerator_reconstruct_refinement_state(jptinit,jptend, &
  totjptend,keepvolsub,keepvolevsub,voltotsub,voltotevsub,do_evaporation)

!***********************************************************************
!
!     OpenACC reconstruction of the accepted-event bead volume and
!     evaporation-volume from the Akima-interpolated cross-section
!     radius, their reference-volume conservation rescale, and the
!     density-to-quantity conversion of mass and charge. Every step is
!     an independent per-bead update or a single reduction over the
!     freshly interpolated mesh, so this reuses the same per-event data
!     window as the Akima kernels without adding any host/device round
!     trip beyond the one already required by the accepted-event state
!     transfer. Only the ordinary (non-variable-mass) density-to-mass
!     branch is implemented; the caller excludes the developer-gated
!     variable-mass path from this device kernel entirely.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: jptinit,jptend,totjptend
  double precision, allocatable, intent(in) :: keepvolsub(:),keepvolevsub(:)
  double precision, intent(in) :: voltotsub,voltotevsub
  logical, intent(in) :: do_evaporation

  integer :: i,ntail
  double precision :: newvoltot,newvoltotev

  ntail=totjptend-jptend

! jetvl/jetms/jetch (like jetxx/jetyy/... in accelerator_eom3_stage) can
! already be present on the device through the persistent topology mapping
! whenever this event did not need capacity growth; when it did, that old
! mapping was already released and the new one is not yet rebound. "create"
! is a safe reference-count increment in the first case and a fresh mapping
! in the second, exactly as fit_akima_accelerator already relies on for the
! same arrays. As in that routine, the structured region names each whole
! allocatable array (never an explicit subrange in the data clause itself,
! which NVHPC rejects as "partially present" the first time an array is
! mapped); only the update/compute clauses below use the exact ranges
! actually needed.
!$acc data create(jetptc,jetcr,jetvl,jetms,jetch)
!$acc update device(jetptc(jptinit:jptend+1),jetcr(jptinit:jptend), &
!$acc& jetms(jptinit:totjptend),jetch(jptinit:totjptend))

! See reconstruct_refinement_state_host for why the true jet endpoint
! (jptend==totjptend, no preserved tail) must not read jetptc(totjptend+1).
!$acc parallel loop present(jetptc,jetcr,jetvl)
  do i=jptinit,jptend
    if(i==totjptend)then
      if(i>jptinit)then
        jetvl(i)=lengthpath*(jetptc(i)-jetptc(i-1))*Pi*jetcr(i)**2.d0
      else
        jetvl(i)=resolution*Pi*jetcr(i)**2.d0
      endif
    else
      jetvl(i)=lengthpath*(jetptc(i+1)-jetptc(i))*Pi*jetcr(i)**2.d0
    endif
  enddo
!$acc end parallel loop

  if(ntail>0)then
!$acc parallel loop present(jetvl) copyin(keepvolsub(1:ntail))
    do i=1,ntail
      jetvl(jptend+i)=keepvolsub(i)
    enddo
!$acc end parallel loop
  endif

  newvoltot=0.d0
!$acc parallel loop present(jetvl) reduction(+:newvoltot)
  do i=jptinit,jptend
    newvoltot=newvoltot+jetvl(i)
  enddo
!$acc end parallel loop

!$acc parallel loop present(jetvl)
  do i=jptinit,jptend
    jetvl(i)=jetvl(i)*voltotsub/newvoltot
  enddo
!$acc end parallel loop

!$acc parallel loop present(jetms,jetch,jetvl)
  do i=jptinit,totjptend
    jetms(i)=jetms(i)*jetvl(i)
    jetch(i)=jetch(i)*jetvl(i)
  enddo
!$acc end parallel loop

!$acc update self(jetvl(jptinit:totjptend),jetms(jptinit:totjptend), &
!$acc& jetch(jptinit:totjptend))
!$acc end data

  if(do_evaporation)then

!$acc data create(jetptc,jetce,jetve)
!$acc update device(jetptc(jptinit:jptend+1),jetce(jptinit:jptend))

!$acc parallel loop present(jetptc,jetce,jetve)
    do i=jptinit,jptend
      if(i==totjptend)then
        if(i>jptinit)then
          jetve(i)=lengthpath*(jetptc(i)-jetptc(i-1))*Pi*jetce(i)**2.d0
        else
          jetve(i)=resolution*Pi*jetce(i)**2.d0
        endif
      else
        jetve(i)=lengthpath*(jetptc(i+1)-jetptc(i))*Pi*jetce(i)**2.d0
      endif
    enddo
!$acc end parallel loop

    if(ntail>0)then
!$acc parallel loop present(jetve) copyin(keepvolevsub(1:ntail))
      do i=1,ntail
        jetve(jptend+i)=keepvolevsub(i)
      enddo
!$acc end parallel loop
    endif

    newvoltotev=0.d0
!$acc parallel loop present(jetve) reduction(+:newvoltotev)
    do i=jptinit,jptend
      newvoltotev=newvoltotev+jetve(i)
    enddo
!$acc end parallel loop

!$acc parallel loop present(jetve)
    do i=jptinit,jptend
      jetve(i)=jetve(i)*voltotevsub/newvoltotev
    enddo
!$acc end parallel loop

!$acc update self(jetve(jptinit:totjptend))
!$acc end data

  endif

  return

 end subroutine accelerator_reconstruct_refinement_state

#ifdef JETSPIN_COMPARE_REFINEMENT_ASSEMBLY
 subroutine compare_refinement_assembly_device(jptinit,jptend,totjptend, &
  voltotsub,voltotevsub,do_evaporation)

!***********************************************************************
!
!     Development-only oracle: back up the pre-reconstruction inputs,
!     evaluate the trusted host reconstruction into reference buffers,
!     restore those same inputs, then run the device path so it remains
!     the real, authoritative state -- exactly the JETSPIN_COMPARE_AKIMA
!     pattern applied to this assembly step. Both paths therefore start
!     from an identical pre-reconstruction state and the device kernel
!     is invoked exactly once.
!
!***********************************************************************

  implicit none

  integer, intent(in) :: jptinit,jptend,totjptend
  double precision, intent(in) :: voltotsub,voltotevsub
  logical, intent(in) :: do_evaporation

  double precision, allocatable :: vl_bak(:),ms_bak(:),ch_bak(:),ve_bak(:)
  double precision, allocatable :: vl_ref(:),ms_ref(:),ch_ref(:),ve_ref(:)
  double precision :: vl_abs,vl_rel,ms_abs,ms_rel,ch_abs,ch_rel
  double precision :: ve_abs,ve_rel,scale

  allocate(vl_bak(jptinit:totjptend),ms_bak(jptinit:totjptend))
  allocate(ch_bak(jptinit:totjptend))
  vl_bak=jetvl(jptinit:totjptend)
  ms_bak=jetms(jptinit:totjptend)
  ch_bak=jetch(jptinit:totjptend)
  if(do_evaporation)then
    allocate(ve_bak(jptinit:totjptend))
    ve_bak=jetve(jptinit:totjptend)
  endif

  call reconstruct_refinement_state_host(jptinit,jptend,totjptend, &
   voltotsub,voltotevsub,do_evaporation)

  allocate(vl_ref(jptinit:totjptend),ms_ref(jptinit:totjptend))
  allocate(ch_ref(jptinit:totjptend))
  vl_ref=jetvl(jptinit:totjptend)
  ms_ref=jetms(jptinit:totjptend)
  ch_ref=jetch(jptinit:totjptend)
  if(do_evaporation)then
    allocate(ve_ref(jptinit:totjptend))
    ve_ref=jetve(jptinit:totjptend)
  endif

! Restore the exact pre-reconstruction inputs (densities, zeroed volumes)
! that the host reference above consumed, so the device kernel starts
! from the same state instead of from the host's already-converted
! mass/charge quantities.
  jetvl(jptinit:totjptend)=vl_bak
  jetms(jptinit:totjptend)=ms_bak
  jetch(jptinit:totjptend)=ch_bak
  if(do_evaporation)jetve(jptinit:totjptend)=ve_bak

  call accelerator_reconstruct_refinement_state(jptinit,jptend,totjptend, &
   keepvol,keepvolev,voltotsub,voltotevsub,do_evaporation)

  vl_abs=maxval(dabs(jetvl(jptinit:totjptend)-vl_ref))
  scale=max(maxval(dabs(vl_ref)),1.d-30)
  vl_rel=vl_abs/scale
  ms_abs=maxval(dabs(jetms(jptinit:totjptend)-ms_ref))
  scale=max(maxval(dabs(ms_ref)),1.d-30)
  ms_rel=ms_abs/scale
  ch_abs=maxval(dabs(jetch(jptinit:totjptend)-ch_ref))
  scale=max(maxval(dabs(ch_ref)),1.d-30)
  ch_rel=ch_abs/scale
  ve_abs=0.d0
  ve_rel=0.d0
  if(do_evaporation)then
    ve_abs=maxval(dabs(jetve(jptinit:totjptend)-ve_ref))
    scale=max(maxval(dabs(ve_ref)),1.d-30)
    ve_rel=ve_abs/scale
  endif

  if(idrank==0)write(6,'(a,6(a,es12.4))') &
   'Refinement assembly device comparison:', &
   ' volume_max_abs=',vl_abs,' volume_max_rel=',vl_rel, &
   ' mass_max_abs=',ms_abs,' mass_max_rel=',ms_rel, &
   ' charge_max_abs=',ch_abs,' charge_max_rel=',ch_rel
  if(do_evaporation .and. idrank==0)write(6,'(a,2(a,es12.4))') &
   'Refinement assembly evaporation comparison:', &
   ' evap_volume_max_abs=',ve_abs,' evap_volume_max_rel=',ve_rel

  deallocate(vl_bak,ms_bak,ch_bak,vl_ref,ms_ref,ch_ref)
  if(do_evaporation)deallocate(ve_bak,ve_ref)

  return

 end subroutine compare_refinement_assembly_device
#endif
#endif

 subroutine fit_jetpt_spline()
 
!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the array jetpt which 
!     parametrizes the nanojet if the spline fitting is performed
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
  
  if(mydoallocate)then
    deallocate(jetpt)
    allocate(jetpt(0:mxnpjet))
  endif
  jetpt(0:npjet)=jetptc(0:npjet)
  jetpt(npjet+1:mxnpjet)=0.d0
  
  
  return
  
 end subroutine fit_jetpt_spline
 
 subroutine store_backup(timesub)

!***********************************************************************
!     
!     JETSPIN subroutine for storing all the data on the service 
!     backup arrays if it is requested
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(inout) :: timesub
  
  inpjetbackup=inpjet
  npjetbackup=npjet
  timebackup=timesub
  
  call allocate_backup(timesub)
  
  if(idrank==0)then
  
  jetxxbak(0:npjet)=jetxx(0:npjet)
  jetyybak(0:npjet)=jetyy(0:npjet)
  jetzzbak(0:npjet)=jetzz(0:npjet)
  jetvxbak(0:npjet)=jetvx(0:npjet)
  jetvybak(0:npjet)=jetvy(0:npjet)
  jetvzbak(0:npjet)=jetvz(0:npjet)
  jetstbak(0:npjet)=jetst(0:npjet)
  jetmsbak(0:npjet)=jetms(0:npjet)
  jetchbak(0:npjet)=jetch(0:npjet)
  jetvlbak(0:npjet)=jetvl(0:npjet)
  if(levaporation)jetvebak(0:npjet)=jetve(0:npjet)
  
  endif
  
  return
  
 end subroutine store_backup
 
 subroutine restore_backup(timesub)

!***********************************************************************
!     
!     JETSPIN subroutine for restoring all the data from the service 
!     backup arrays if it is requested
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(inout) :: timesub
  
  inpjet=inpjetbackup
  npjet=npjetbackup
  timesub=timebackup
  
  if(idrank==0)then
  
  jetxx(0:npjet)=jetxxbak(0:npjet)
  jetyy(0:npjet)=jetyybak(0:npjet)
  jetzz(0:npjet)=jetzzbak(0:npjet)
  jetvx(0:npjet)=jetvxbak(0:npjet)
  jetvy(0:npjet)=jetvybak(0:npjet)
  jetvz(0:npjet)=jetvzbak(0:npjet)
  jetst(0:npjet)=jetstbak(0:npjet)
  jetms(0:npjet)=jetmsbak(0:npjet)
  jetch(0:npjet)=jetchbak(0:npjet)
  jetvl(0:npjet)=jetvlbak(0:npjet)
  if(levaporation)jetve(0:npjet)=jetvebak(0:npjet)
  
  
  endif
  
  call bcast_world_darr(jetxx,npjet+1)
  call bcast_world_darr(jetyy,npjet+1)
  call bcast_world_darr(jetzz,npjet+1)
  call bcast_world_darr(jetvx,npjet+1)
  call bcast_world_darr(jetvy,npjet+1)
  call bcast_world_darr(jetvz,npjet+1)
  call bcast_world_darr(jetst,npjet+1)
  call bcast_world_darr(jetms,npjet+1)
  call bcast_world_darr(jetch,npjet+1)
  call bcast_world_darr(jetvl,npjet+1)
  if(levaporation)call bcast_world_darr(jetve,npjet+1)
  
  return
  
 end subroutine restore_backup
 
 subroutine allocate_backup(timesub)
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the service backup arrays
!     which are used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!***********************************************************************
 
  implicit none
  
  double precision, intent(inout) :: timesub
  
  logical, save :: mioallocfirst=.true.
  
  if(idrank==0)then
    if(mxnpjet>nbackup)then
      nbackup=mxnpjet
      if(.not.mioallocfirst)then
        deallocate(jetxxbak,jetyybak,jetzzbak)
        deallocate(jetvxbak,jetvybak,jetvzbak)
        deallocate(jetstbak,jetmsbak,jetchbak)
        deallocate(jetvlbak)
        if(levaporation)deallocate(jetvebak)
      endif
      allocate(jetxxbak(0:mxnpjet),jetyybak(0:mxnpjet), &
       jetzzbak(0:mxnpjet))
      allocate(jetvxbak(0:mxnpjet),jetvybak(0:mxnpjet), &
       jetvzbak(0:mxnpjet))
      allocate(jetstbak(0:mxnpjet),jetmsbak(0:mxnpjet), &
       jetchbak(0:mxnpjet))
      allocate(jetvlbak(0:mxnpjet))
      if(levaporation)allocate(jetvebak(0:mxnpjet))
    endif
  endif
  
  mioallocfirst=.false.
  
  return
  
 end subroutine allocate_backup
 
 subroutine allocate_jetlthr()
 
!***********************************************************************
!     
!     JETSPIN subroutine for allocating the service backup arrays
!     which are used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  logical, save :: mioallocfirst=.true.
  
  if(mxnpjet>njetlthr)then
    njetlthr=mxnpjet
    if(.not.mioallocfirst)then
      deallocate(jetlthr)
    endif
    allocate(jetlthr(0:mxnpjet))
  endif
  
  mioallocfirst=.false.
  
  return
  
 end subroutine allocate_jetlthr
 
 subroutine allocate_array_myindex(jmiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array myindex 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jmiomax
  
  if(nmyindex/=0)then
    if(jmiomax>nmyindex)then
      deallocate(myindex)
      nmyindex=jmiomax+10
      allocate(myindex(nmyindex,6))
    endif
  else
    nmyindex=jmiomax+10
    allocate(myindex(nmyindex,6))
  endif
  
  return
  
 end subroutine allocate_array_myindex
 
 subroutine allocate_array_keepvol(jmiomax)

!***********************************************************************
!     
!     JETSPIN subroutine for reallocating the service array myindex 
!     which is used within this module
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2017
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: jmiomax
  
  if(nkeepvol/=0)then
    if(jmiomax>nkeepvol)then
      deallocate(keepvol)
      if(levaporation)deallocate(keepvolev)
      nkeepvol=jmiomax+10
      allocate(keepvol(nkeepvol))
      if(levaporation)allocate(keepvolev(nkeepvol))
    endif
  else
    nkeepvol=jmiomax+10
    allocate(keepvol(nkeepvol))
    if(levaporation)allocate(keepvolev(nkeepvol))
  endif
  
  return
  
 end subroutine allocate_array_keepvol
 
 end module dynamic_refinement_mod
 
