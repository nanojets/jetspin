
 program traj2dcd
 
!***********************************************************************
!
!     tool code for converting the trajectory bynary file obtained
!     a JETSPIN run into a xyz and a dcd file which can be opened
!     in a visualization code (e.g. VMD)
!
!     compiling command: gfortran -O2 -o main.x traj2dcd.f90
!     running command  : ./main.x
!     notes  : the file traj.dat should be in the same directory of
!              the executive file main.x
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2016
!     
!***********************************************************************

  implicit none
  
  logical :: lexist
  
  character(len=*), parameter :: nfilein='traj.dat'
  integer, parameter :: myinp=15
  integer, parameter :: jincrement=100
  logical, parameter :: unscale=.true.
  
  logical :: lallocated=.false.
  logical :: lredo
  integer :: istep,mynframes
  
  integer :: maxnjetbd,maxnjetbd2,njetbd,njetbd2
  integer :: nmaxjet
  integer :: inpjet,npjet
  integer :: ncutoff,integrator,units,doreorder,inpjetread,npjetread
  integer :: addread,remread
  integer :: ncounterevel,ncounterevelrel,ncountergeom
  integer :: reprinttime,ncounterivel,ncounterlpath,irefinementdone
  integer :: nmulstep,nmulstepdone,naddtrack,nremtrack
  integer :: maxnatms,natms,sprintdat,systype
  logical :: linserting,linserted
  double precision :: surfacet,airdragamp(3),aird,airv,icrossec
  double precision :: resolution,findex
  double precision :: ilength,imassa,icharge,istress,ivelocity,pfreq
  double precision :: pampl,mu,G,yieldstress,h,V0,tstep,consistency
  double precision :: massscale,chargescale,lengthscale,tao
  double precision :: imassadev,lencorrmassa,ichargedev,lencorrcharge
  double precision :: addedcharge,addedmass,removedcharge
  double precision :: countericurr,meanlpath,timesub
  double precision :: q,v,fve,fvere,Hg,Lrg,Gr,ks,ksre,att,Li,lire
  double precision :: counterecurr,meanicurr,meanecurr,removedmass
  double precision :: counterimass,counteremass,meanimass,meanemass
  double precision :: counterevel,counterevelrel,meanevel,meanevelrel
  double precision :: counterelen,counterecross,meanelen,meanecross
  double precision :: counterivel,meanivel,counterlpath,meancputime
  double precision :: velext,corr,radcorr,lenprobmassa,massratio
  double precision :: lenthresholdbead,oldgaussn,refbeadstartfit
  double precision :: refinementevery,refinementstart
  double precision :: refinementthreshold
  
  integer, allocatable, dimension(:) :: jetin
  logical, allocatable, dimension(:) :: jetbd
  double precision, allocatable, dimension(:) :: jetxx,jetyy,jetzz
  double precision, allocatable, dimension(:) :: jetst,jetvx,jetvy,jetvz
  double precision, allocatable, dimension(:) :: jetms,jetch
  double precision, allocatable, dimension(:) :: jetvl
  
  logical ::lrefinement,ltagbeads,lrefinementthreshold, &
         lrefbeadstart,llenthresholdbead, &
         lrefinementevery,lrefinementstart,lmassavariable,limassadev, &
         llencorrmassa,lfirstmass
         
  integer :: refinementcons,typemass
  
  
  inquire(file=trim(nfilein),exist=lexist)
  
  if(.not. lexist)then
    write(6,'(a)')'ERROR:'
    write(6,'(a)')'input file traj.dat not found!'
    stop
  endif
  
  open(unit=myinp,file=trim(nfilein),status='old',action='read', &
   form='unformatted')
  
  call open_dat_file(myinp,trim(nfilein))
    
  lredo=.true.
  call read_dat_frame(myinp,lredo)
  
  maxnatms=0
  maxnjetbd=0
  maxnjetbd2=0
  istep=0
  do
    
    call read_dat_frame(myinp,lredo)
    
    if(.not.lredo)write(6,*)istep,lredo
    
    if(.not.lredo)exit
    
    istep=istep+1
    maxnatms=max(maxnatms,natms)
    call compute_maxnjetbd()
    
    if(mod(istep,1000)==0)write(6,'(4i8)')istep,natms,maxnatms
    
  enddo
  
  call close_dat_file(myinp)
  
  mynframes=istep
  
  write(6,'(a,i10)')'Max number of beads = ',maxnatms
  write(6,'(a,i10)')'in frames = ',mynframes
  
  write(6,'(a)')'Writing xyz file'
  
  call write_xyz_dcd('trajout.xyz')
  
  call open_dat_file(myinp,trim(nfilein))
    
  lredo=.true.
  call read_dat_frame(myinp,lredo,.true.)
  
  write(6,'(a)')'Writing dcd file'
  call write_dcd_header('trajout.dcd',222,maxnatms,mynframes, &
   mynframes,0,1)
  
  do istep=1,mynframes
    
    jetxx(:)=0.d0
    jetyy(:)=0.d0
    jetzz(:)=0.d0
    call read_dat_frame(myinp,lredo,.true.)
    
    if(.not.lredo)write(6,*)istep,lredo
    
    if(.not.lredo)exit
    
    call write_dcd_frame(222,jetxx,jetyy,jetzz,h,h,h,1.d0,maxnatms)
    
  enddo
    
  call close_dat_file(myinp) 
  
  write(6,'(a)')'Program closed!' 
  
  stop
  
  
 contains
 
 subroutine compute_maxnjetbd()
 
  implicit none
  
  integer :: i
  
  if(sprintdat/=6)return
  
  
  njetbd=0
  njetbd2=0
  do i=inpjet,npjet
    if(jetbd(i))then
      njetbd=njetbd+1
    else
      njetbd2=njetbd2+1
    endif
  enddo
  maxnjetbd=max(maxnjetbd,njetbd)
  maxnjetbd2=max(maxnjetbd2,njetbd2)
  
  return
  
 end subroutine compute_maxnjetbd
 
 subroutine open_dat_file(fileout,filename)
 
  implicit none
  
  integer, intent(in) :: fileout
  character(len=*), intent(in) :: filename
  
  
  open(unit=fileout,file=filename,form='unformatted',status='old',action='read')
  
  return
  
 end subroutine open_dat_file
 
 subroutine read_dat_frame(fileout,lstart,lsecond)
  
  implicit none
  
  integer, intent(in) :: fileout
  logical, intent(out) :: lstart
  logical, intent(in), optional :: lsecond
  
  integer :: i,j,mioind
  real(4) :: rtemp(10) 
  
  logical :: lredo
  double precision :: dtemp(12)
  integer :: itemp(12)
  
  lredo=.true.
!  do while(lredo)
  read(fileout,end=120)lstart
  
  if(.not.lstart)then
    lredo=.false.
    return
  endif
  read(fileout)natms,sprintdat,systype,timesub,linserted
  inpjet=0
  npjet=natms-1
  if(present(lsecond))then
    if(lsecond)then
      call allocate_jet(maxnatms)
    else
      call allocate_jet(npjet)
    endif
  else
    call allocate_jet(npjet)
  endif
  
  
  if(sprintdat==0)then
    read(fileout)ncutoff,systype,integrator,units
    read(fileout)resolution,ilength,imassa,icharge,icrossec,istress,ivelocity
    read(fileout)pfreq,pampl,airdragamp(1),airdragamp(2),airdragamp(3),aird,airv
    read(fileout)mu,G,yieldstress,h,V0,surfacet,tstep,consistency,findex
    read(fileout)massscale,chargescale,lengthscale,tao
    read(fileout)imassadev,lencorrmassa,ichargedev,lencorrcharge,velext
    
    read(fileout)mioind
    
    if(mioind>=1)then
      read(fileout)q,v,fve,fvere,Hg,Lrg,Gr,ks,ksre,att,Li,lire
    endif
    
    if(mioind>=2)then
      read(fileout)naddtrack,nremtrack,linserting
    endif
  
    if(mioind>=3)then
        read(fileout)addedcharge,removedcharge,countericurr, &
         counterecurr,meanicurr,meanecurr,addedmass,removedmass, &
         counterimass,counteremass,meanimass,meanemass
    endif
    
    if(mioind>=4)then
        read(fileout)counterevel,counterevelrel,meanevel,meanevelrel, &
         counterelen,counterecross,meanelen,meanecross,meancputime, &
         counterivel,meanivel,counterlpath
    endif
    
    if(mioind>=5)then
        read(fileout)meanlpath,(dtemp(i),i=1,11)
    endif
    
    if(mioind>=6)then
          read(fileout)ncounterevel,ncounterevelrel,ncountergeom, &
       reprinttime,ncounterivel,ncounterlpath,irefinementdone, &
         nmulstep,nmulstepdone,(itemp(i),i=1,3)
    endif
    
    if(mioind>=7)then
        read(fileout)lrefinement,ltagbeads,lrefinementthreshold, &
         lrefbeadstart,llenthresholdbead, &
         lrefinementevery,lrefinementstart,lmassavariable,limassadev, &
         llencorrmassa,lmassavariable,lfirstmass
    endif
    
    if(mioind>=8)then
          read(fileout)refinementcons,typemass,(itemp(i),i=1,10)
    endif
    
    if(mioind>=9)then
        read(fileout)refinementthreshold,refbeadstartfit,lencorrmassa, &
         lenprobmassa,massratio,imassadev,oldgaussn,corr,radcorr, &
         lenthresholdbead,refinementevery,refinementstart
    endif
    
  elseif(sprintdat==1)then
    lredo=.false.
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,3)
        jetxx(i)=dble(rtemp(1))
        jetst(i)=dble(rtemp(2))
        jetvx(i)=dble(rtemp(3))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
      endif
    case default
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,7)
        jetxx(i)=dble(rtemp(1))
        jetyy(i)=dble(rtemp(2))
        jetzz(i)=dble(rtemp(3))
        jetst(i)=dble(rtemp(4))
        jetvx(i)=dble(rtemp(5))
        jetvy(i)=dble(rtemp(6))
        jetvz(i)=dble(rtemp(7))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetyy(:)=jetyy(:)*lengthscale
        jetzz(:)=jetzz(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetvy(:)=jetvy(:)*lengthscale/tao
        jetvz(:)=jetvz(:)*lengthscale/tao
      endif
    end select
  elseif(sprintdat==2)then
    lredo=.false.
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,5)
        jetxx(i)=dble(rtemp(1))
        jetst(i)=dble(rtemp(2))
        jetvx(i)=dble(rtemp(3))
        jetms(i)=dble(rtemp(4))
        jetch(i)=dble(rtemp(5))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
      endif
    case default
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,9)
        jetxx(i)=dble(rtemp(1))
        jetyy(i)=dble(rtemp(2))
        jetzz(i)=dble(rtemp(3))
        jetst(i)=dble(rtemp(4))
        jetvx(i)=dble(rtemp(5))
        jetvy(i)=dble(rtemp(6))
        jetvz(i)=dble(rtemp(7))
        jetms(i)=dble(rtemp(8))
        jetch(i)=dble(rtemp(9))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetyy(:)=jetyy(:)*lengthscale
        jetzz(:)=jetzz(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetvy(:)=jetvy(:)*lengthscale/tao
        jetvz(:)=jetvz(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
      endif
    end select
  elseif(sprintdat==3)then
    lredo=.false.
    read(fileout)doreorder,inpjetread,npjetread
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,5),itemp(1)
        jetxx(i)=dble(rtemp(1))
        jetst(i)=dble(rtemp(2))
        jetvx(i)=dble(rtemp(3))
        jetms(i)=dble(rtemp(4))
        jetch(i)=dble(rtemp(5))
        jetin(i)=itemp(1)
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
      endif
    case default
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,9),itemp(1)
        jetxx(i)=dble(rtemp(1))
        jetyy(i)=dble(rtemp(2))
        jetzz(i)=dble(rtemp(3))
        jetst(i)=dble(rtemp(4))
        jetvx(i)=dble(rtemp(5))
        jetvy(i)=dble(rtemp(6))
        jetvz(i)=dble(rtemp(7))
        jetms(i)=dble(rtemp(8))
        jetch(i)=dble(rtemp(9))
        jetin(i)=itemp(1)
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetyy(:)=jetyy(:)*lengthscale
        jetzz(:)=jetzz(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetvy(:)=jetvy(:)*lengthscale/tao
        jetvz(:)=jetvz(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
      endif
    end select
  elseif(sprintdat==4)then
    lredo=.false.
    read(fileout)doreorder,addread,remread 
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,5),itemp(1)
        jetxx(i)=dble(rtemp(1))
        jetst(i)=dble(rtemp(2))
        jetvx(i)=dble(rtemp(3))
        jetms(i)=dble(rtemp(4))
        jetch(i)=dble(rtemp(5))
        jetin(i)=itemp(1)
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
      endif
    case default
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,9),itemp(1)
        jetxx(i)=dble(rtemp(1))
        jetyy(i)=dble(rtemp(2))
        jetzz(i)=dble(rtemp(3))
        jetst(i)=dble(rtemp(4))
        jetvx(i)=dble(rtemp(5))
        jetvy(i)=dble(rtemp(6))
        jetvz(i)=dble(rtemp(7))
        jetms(i)=dble(rtemp(8))
        jetch(i)=dble(rtemp(9))
        jetin(i)=itemp(1)
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetyy(:)=jetyy(:)*lengthscale
        jetzz(:)=jetzz(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetvy(:)=jetvy(:)*lengthscale/tao
        jetvz(:)=jetvz(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
      endif
    end select
  elseif(sprintdat==5)then
    lredo=.false.
    read(fileout)doreorder,addread,remread
    select case(systype)
    case (1:2)
      jetxx(:)=0.d0
      jetst(:)=0.d0
      jetvx(:)=0.d0
      jetms(:)=0.d0
      jetch(:)=0.d0
      jetvl(:)=0.d0
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,6)
        jetxx(i)=dble(rtemp(1))
        jetst(i)=dble(rtemp(2))
        jetvx(i)=dble(rtemp(3))
        jetms(i)=dble(rtemp(4))
        jetch(i)=dble(rtemp(5))
        jetvl(i)=dble(rtemp(6))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
        jetvl(:)=jetvl(:)*lengthscale**3.d0
      endif
    case default
      jetxx(:)=0.d0
      jetyy(:)=0.d0
      jetzz(:)=0.d0
      jetst(:)=0.d0
      jetvx(:)=0.d0
      jetvy(:)=0.d0
      jetvz(:)=0.d0
      jetms(:)=0.d0
      jetch(:)=0.d0
      jetvl(:)=0.d0
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,10)
        jetxx(i)=dble(rtemp(1))
        jetyy(i)=dble(rtemp(2))
        jetzz(i)=dble(rtemp(3))
        jetst(i)=dble(rtemp(4))
        jetvx(i)=dble(rtemp(5))
        jetvy(i)=dble(rtemp(6))
        jetvz(i)=dble(rtemp(7))
        jetms(i)=dble(rtemp(8))
        jetch(i)=dble(rtemp(9))
        jetvl(i)=dble(rtemp(10))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetyy(:)=jetyy(:)*lengthscale
        jetzz(:)=jetzz(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetvy(:)=jetvy(:)*lengthscale/tao
        jetvz(:)=jetvz(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
        jetvl(:)=jetvl(:)*lengthscale**3.d0
      endif
    end select
  elseif(sprintdat==6)then
    lredo=.false.
    read(fileout)doreorder,addread,remread
    select case(systype)
    case (1:2)
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,6),jetbd(i)
        jetxx(i)=dble(rtemp(1))
        jetst(i)=dble(rtemp(2))
        jetvx(i)=dble(rtemp(3))
        jetms(i)=dble(rtemp(4))
        jetch(i)=dble(rtemp(5))
        jetvl(i)=dble(rtemp(6))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
        jetvl(:)=jetvl(:)*lengthscale**3.d0
      endif
    case default
      do i=inpjet,npjet
        read(fileout)(rtemp(j),j=1,10),jetbd(i)
        jetxx(i)=dble(rtemp(1))
        jetyy(i)=dble(rtemp(2))
        jetzz(i)=dble(rtemp(3))
        jetst(i)=dble(rtemp(4))
        jetvx(i)=dble(rtemp(5))
        jetvy(i)=dble(rtemp(6))
        jetvz(i)=dble(rtemp(7))
        jetms(i)=dble(rtemp(8))
        jetch(i)=dble(rtemp(9))
        jetvl(i)=dble(rtemp(10))
      end do
      if(unscale)then
        jetxx(:)=jetxx(:)*lengthscale
        jetyy(:)=jetyy(:)*lengthscale
        jetzz(:)=jetzz(:)*lengthscale
        jetst(:)=jetst(:)*G
        jetvx(:)=jetvx(:)*lengthscale/tao
        jetvy(:)=jetvy(:)*lengthscale/tao
        jetvz(:)=jetvz(:)*lengthscale/tao
        jetms(:)=jetms(:)*massscale
        jetch(:)=jetch(:)*chargescale
        jetvl(:)=jetvl(:)*lengthscale**3.d0
      endif
    end select
  endif
  
  return
  
 120 continue
  
  lstart=.false.
  lredo=.false.
  
  return
    
 end subroutine read_dat_frame
 
 subroutine close_dat_file(fileout)
 
  implicit none

  integer, intent(in) :: fileout
  
  close(fileout)
  
  return
  
 end subroutine close_dat_file
 
 subroutine write_xyz_dcd(filename)

    character(len=*), intent(in) :: filename
    integer :: j
    integer,parameter :: fileout=122

    open(fileout,file=trim(filename),status='replace')

    write(fileout,*) maxnatms
    write(fileout,*)
    do j=1,maxnatms
       write(fileout,"(a8,3f10.5)")'C       ',0.d0,0.d0,0.d0
    end do
    

   close(fileout)


 end subroutine write_xyz_dcd
 
 subroutine write_xyz_dcd_bd(filename)

    character(len=*), intent(in) :: filename
    
    integer :: j
    integer,parameter :: fileout=122

    open(fileout,file=trim(filename),status='replace')
    
    maxnatms=maxnjetbd+maxnjetbd2

    write(fileout,*) maxnatms
    write(fileout,*)
    do j=1,maxnjetbd2
       write(fileout,"(a8,3f10.5)")'C       ',0.d0,0.d0,0.d0
    end do
    
    do j=1,maxnjetbd
       write(fileout,"(a8,3f10.5)")'N       ',0.d0,0.d0,0.d0
    end do

   close(fileout)


 end subroutine write_xyz_dcd_bd
 
 subroutine write_dcd_header(filename,fileunit,natmssub,nframessub, &
  nstepsub,nstraj,istraj)

    ! Program for writing a NAMD dcd file 
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: fileunit,natmssub,nframessub,nstepsub,nstraj,istraj
    


    character(len=4) :: coor
    character(len=80) :: title1,title2
    real(8) :: dum,delta
    integer :: zer,fin
    integer :: ntitle

    integer,dimension(8) :: values
    

    title1="Written by BOB"
    call date_and_time(VALUES=values)
    write(title2,"('Created ',i2,1x,a9,1x,i4,' at ',i2,':',i2)")values(3), &
      write_month(values(2)),values(1),values(5),values(6)
    coor = 'CORD'
    zer = 0
    fin=24
    ntitle=2
    dum=0.d0
    delta=1.d0/(huge(1.d0)*0.001d0)

    open(fileunit,file=filename,form='unformatted',status='replace')
    write(fileunit) coor,nframessub,nstraj,istraj,nstepsub, &
     zer,zer,zer,zer,zer,delta,zer,zer,zer,zer,zer,zer,zer,zer,fin
    write(fileunit) ntitle,title1,title2
    write(fileunit) natmssub
    

    return


 end subroutine write_dcd_header
 
 subroutine write_dcd_frame(fileunit,x,y,z,a,b,c,rescale,natmssub)

    ! Program for writing a NAMD dcd file 
    
    integer, intent(in) :: fileunit,natmssub
    real(8), intent(in), allocatable, dimension(:) :: x, y, z
    real(8), intent(in) :: a, b, c, rescale


    integer :: j
    real(4),dimension(natmssub) :: x4,y4,z4
    real(8) :: dum,as,bs,cs

    dum=0.d0
    as=a*rescale
    bs=b*rescale
    cs=c*rescale
    
    
    x4(1:natmssub)=real(x(1:natmssub)*rescale,kind=4)
    y4(1:natmssub)=real(y(1:natmssub)*rescale,kind=4)
    z4(1:natmssub)=real(z(1:natmssub)*rescale,kind=4)


    write(fileunit)   as, dum, bs, dum, dum, cs
    write(fileunit)  (x4(j),j=1,natmssub)
    write(fileunit)  (y4(j),j=1,natmssub)
    write(fileunit)  (z4(j),j=1,natmssub)


    return


 end subroutine write_dcd_frame
 
 function write_month(mese)

   ! writes only ever nfreq steps

    integer, intent(in) :: mese
    character(len=9),parameter,dimension(12) :: smesi=(/'  January', &
    ' February','    March','    April', &
    '      May','     June','     July','   August','September', &
    '  October',' November',' December'/)
    character(len=9) :: write_month

    if(mese.lt.1.or.mese.gt.12)then
      write(6,'(a)')'Error - month < 1 or > 12'
      stop
    endif
    write_month=smesi(mese)


    return

 end function write_month
 
 subroutine allocate_jet(npjetsub)
 
  implicit none
  
  integer, intent(in) :: npjetsub
  
  
  if(lallocated)then
    if(npjetsub>=nmaxjet)then
      deallocate(jetbd)
      deallocate(jetxx)
      deallocate(jetyy)
      deallocate(jetzz)
      deallocate(jetst)
      deallocate(jetvx)
      deallocate(jetvy)
      deallocate(jetvz)
      deallocate(jetms)
      deallocate(jetch)
      deallocate(jetvl)
      deallocate(jetin)
      do while(npjetsub>=nmaxjet)
        nmaxjet=nmaxjet+jincrement
      enddo
      allocate(jetbd(0:nmaxjet-1))
      allocate(jetxx(0:nmaxjet-1))
      allocate(jetyy(0:nmaxjet-1))
      allocate(jetzz(0:nmaxjet-1))
      allocate(jetst(0:nmaxjet-1))
      allocate(jetvx(0:nmaxjet-1))
      allocate(jetvy(0:nmaxjet-1))
      allocate(jetvz(0:nmaxjet-1))
      allocate(jetms(0:nmaxjet-1))
      allocate(jetch(0:nmaxjet-1))
      allocate(jetvl(0:nmaxjet-1))
      allocate(jetin(0:nmaxjet-1))
    endif
  else 
    nmaxjet=npjetsub+jincrement
    allocate(jetbd(0:nmaxjet-1))
    allocate(jetxx(0:nmaxjet-1))
    allocate(jetyy(0:nmaxjet-1))
    allocate(jetzz(0:nmaxjet-1))
    allocate(jetst(0:nmaxjet-1))
    allocate(jetvx(0:nmaxjet-1))
    allocate(jetvy(0:nmaxjet-1))
    allocate(jetvz(0:nmaxjet-1))
    allocate(jetms(0:nmaxjet-1))
    allocate(jetch(0:nmaxjet-1)) 
    allocate(jetvl(0:nmaxjet-1))
    allocate(jetin(0:nmaxjet-1)) 
    lallocated=.true.
  endif
  
  return
  
 end subroutine allocate_jet
 
 end program traj2dcd
