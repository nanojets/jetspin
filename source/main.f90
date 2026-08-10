

 program JetSpin
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! JetSpin is a Program for Electrospinning Simulations of Nanofibers
!
! This is an experimental code. The authors accept no responsibility
! for the performance of the code or for the correctness of the results.
!
! The code is licensed under Open Software License v. 3.0 (OSL-3.0). 
! The full text of the licence can be found on the website:
! http://opensource.org/licenses/OSL-3.0 
!
! A brief explanation of this license is available on the website:  
! http://rosenlaw.com/OSL3.0-explained.htm  
!
! The software development process has received funding from the 
! European Research Council under the European Union's Seventh Framework
! Programme (FP/2007-2013)/ERC Grant Agreement n. 306357 (NANO-JETS). 
!
! If results obtained with this code are published, an
! appropriate citation would be:
!
! Marco Lauricella, Ivan Coluzza, Giuseppe Pontrelli,
! Dario Pisignano, Sauro Succi,
! JETSPIN: A specific-purpose open-source software for
! electrospinning simulations of nanofibers,  
! Computer Physics Communications, 197 (2015), pp. 227-238. 
!
! author: M. Lauricella
!
! contributors: I. Coluzza,G. Pontrelli, D. Pisignano, S. Succi
!
!                        JETSPIN VERSION 1.22
!
! (May 2017)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use version_mod,    only : init_world,get_rank_world,get_size_world,&
                       alloc_domain,time_world,wall_time_world, &
                       get_sync_world,finalize_world,idrank
  use utility_mod,    only : init_random_seed
#ifdef _OPENACC
  use accelerator_mod, only : accelerator_prepare, &
                         accelerator_update_host_state, &
                         accelerator_is_persistent, &
                         accelerator_remove_bead, &
                         accelerator_update_device_removed, &
                         accelerator_add_bead
#endif
  use profiling_mod, only : profiling_initialize,profiling_reset, &
                       profiling_start,profiling_stop,profiling_report, &
                       prof_integrator,prof_add_bead,prof_remove_bead, &
                       prof_breakup,prof_statistics,prof_erase_bead, &
                       prof_output,prof_restart
  use nanojet_mod,    only : inpjet,npjet,linserted,linserting,lremove,myseed,systype, &
                       tstep,xyzrescale,set_resolution_length, &
                       allocate_jet,set_initial_jet,add_jetbead, &
                       remove_jetbead,erase_jetbead,lengthscale, &
                       pdbrescale,lreadrest,lKVfluid,levaporation, &
                       jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
                       jetms,jetch,jetvl,jetfr,topology_add_total, &
                       topology_remove_total,nremtrack,naddtrack,h,mxnpjet, &
                       resolution,dresolution,thresolution,ivelocity,istress, &
                       imassa,icharge,ivolume,timedeposition
  use breaking_mod,   only : ckeck_breakup
  use dynamic_refinement_mod, only : refinementthreshold, &
                               set_refinement_threshold, &
                               refbeadstartfit
  use integrator_mod, only : initime,endtime,driver_integrator, &
                            prepare_integrator_random_history
  use integrator_kv_ev_mod, only : driver_integrator_KV_ev
  use statistic_mod,  only : statistic_driver
  use io_mod,         only : iprintdat,iprintxyz,lprintdat,lprintxyz,&
                       maxnumxyz,sprintdat,print_logo,read_input,&
                       print_internal_units, &
                       print_adimensional_parameters,&
                       print_adimensional_groups, &
                       print_legend_observables,&
                       allocate_print,outprint_driver,open_xyz_file, &
                       open_dat_file,write_dat_parameter, &
                       write_xyz_frame,write_dat_frame,&
                       finish_print,close_xyz_file,close_dat_file, &
                       write_xyz_singlefile,lprintxyzsing, &
                       iprintxyzsing,open_datrem_file, &
                       write_datrem_frame,close_datrem_file, &
                       write_restart_file,read_restart_file, &
                       write_pdb_singlefile,lprintpdbsing, &
                       iprintpdbsing,timcls,timjob,nrestartdump
  
  implicit none
  
  integer:: nstep
  integer :: nremoved
  
  double precision :: mytime
  double precision :: itime,ctime,ftime
  double precision :: loop_start_time,loop_end_time,loop_elapsed_time
  
  logical :: ladd,lrem,lremdat,ldorefinment,lrecycle
  logical :: lfullhostoutput
  logical :: ltopologysnapshot
  character(len=32) :: topology_snapshot_env
  
  integer :: i,j,k,atype

! set up the communications 
  call init_world()
  
! determine processor identities
  call get_rank_world()
  call get_size_world()
  
  call alloc_domain()
  
! start clock
  call time_world(itime)
  
! print logo on terminal
  call print_logo(6)
  
! read the input file
  call read_input(600,'input.dat')
  
! set the seed
  call init_random_seed(myseed)
  
! define the resolution length step 
  call set_resolution_length()
  
! define the spline threshold
  call set_refinement_threshold()
  
! allocate arrays of the nanofiber quantities
  call allocate_jet()
  
! set the nanofiber quantities
  call set_initial_jet(initime,endtime,refinementthreshold, &
   refbeadstartfit)
  
! print the internal scaling units
  call print_internal_units(6)
  
! print the dimensionless parameter
  call print_adimensional_parameters(6)
  
! print the dimensionless groups
  call print_adimensional_groups(6)
  
! print the list of the output observables as requested in input file
  call print_legend_observables(6)
  
! allocate service arrays for printing modules
  call allocate_print()
  
! initialize the counter of the integration steps  
  nstep=0
  
! initialize the time of the integration steps  
  mytime=initime
  
! initialize variables which keep trace if a bead is added and/or removed
! after the integration step
  ladd=.false.
  lrem=.false.
  
! read the restart file if requested
  call read_restart_file(134,'restart.dat', &
   nstep,mytime)
  
! open the  output 'statdat.dat' file and print first record on terminal
! and output file
  call outprint_driver(nstep,mytime)
  
! open the XYZ formatted output file 
  call open_xyz_file(lprintxyz,120,'traj.xyz')
  
! open the binary file for hit beads (only for developers) 
  call open_datrem_file(lprintdat,122,'bead.dat')
  
! open the binary file (only for developers) 
  call open_dat_file(lprintdat,130,'traj.dat')
  
! write the input parameters on the binary file (only for developers) 
  if(.not. lreadrest)call write_dat_parameter(lprintdat,130,mytime)
  
! initialize lrecycle 
  lrecycle=.true.
  
!***********************************************************************
!     start the time integration
!***********************************************************************
! Initialize an accelerator runtime before the measured region. This is a
! no-op for CPU builds and excludes one-time device setup from loop timing.
#ifdef _OPENACC
  call accelerator_prepare()
#endif
! Pre-generate fixed stochastic benchmark noise before loop timing. CPU and
! GPU executions use the same step/bead/component-indexed history.
  call prepare_integrator_random_history(tstep)
  call profiling_initialize()
  call profiling_reset()
  call get_sync_world()
  call wall_time_world(loop_start_time)
  topology_snapshot_env=''
  call get_environment_variable('JETSPIN_TOPOLOGY_SNAPSHOT', &
   topology_snapshot_env)
  ltopologysnapshot=trim(topology_snapshot_env)=='1'
  if(ltopologysnapshot .and. idrank==0)then
    open(unit=151,file='topology-state.dat',status='replace',action='write')
  endif

  do while (lrecycle)
  
!   update the counter
    nstep=nstep+1
    
!   check recycle loop
    lrecycle=((dble(nstep)*tstep)<endtime)
    
!   integrate the system
    call profiling_start(prof_integrator)
    if(lKVfluid.and.levaporation)then
      call driver_integrator_KV_ev(mytime,tstep,nstep,ldorefinment)
    else
      call driver_integrator(mytime,tstep,nstep,ldorefinment)
    endif
    call profiling_stop(prof_integrator)
    
!   check if a new bead should be added and/or removed
    call profiling_start(prof_add_bead)
#ifdef _OPENACC
    if(accelerator_is_persistent() .and. linserting)then
      timedeposition=timedeposition+tstep
      call accelerator_add_bead(npjet,mxnpjet,linserted,ladd,resolution, &
       dresolution,thresolution,ivelocity,istress,imassa,icharge,ivolume, &
       jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz, &
       jetms,jetch,jetvl,jetfr)
      if(ladd)then
        naddtrack=naddtrack+1
        topology_add_total=topology_add_total+1
        timedeposition=0.d0
      endif
    else
      call add_jetbead(nstep,mytime,ladd)
    endif
#else
    call add_jetbead(nstep,mytime,ladd)
#endif
    call profiling_stop(prof_add_bead)
    call profiling_start(prof_remove_bead)
#ifdef _OPENACC
    if(accelerator_is_persistent() .and. linserting .and. lremove)then
      call accelerator_remove_bead(inpjet,npjet,h,jetxx,jetyy,jetzz,jetst, &
       jetvx,jetvy,jetvz,jetms,jetch,jetvl,jetfr,nremoved,lrem)
      lremdat=lrem
      if(lrem)then
        nremtrack=nremtrack+nremoved
        topology_remove_total=topology_remove_total+nremoved
      endif
    else
      call remove_jetbead(nstep,nremoved,mytime,lrem,lremdat)
    endif
#else
    call remove_jetbead(nstep,nremoved,mytime,lrem,lremdat)
#endif
    call profiling_stop(prof_remove_bead)
    if(idrank==0 .and. (ladd .or. lrem))then
      write(6,'(a,i0,a,l1,a,i0,a,i0)')'Topology event: step=',nstep, &
       ' add=',ladd,' remove=',merge(nremoved,0,lrem),' active=',npjet-inpjet
    endif
    if(ltopologysnapshot .and. idrank==0 .and. (ladd .or. lrem))then
#ifdef _OPENACC
      if(accelerator_is_persistent())call accelerator_update_host_state( &
       npjet,jetxx,jetyy,jetzz,jetst,jetvx,jetvy,jetvz,nstep)
#endif
      write(151,'(a,4(1x,i0),1x,l1)')'event',nstep,inpjet,npjet, &
       npjet-inpjet,ladd
      do i=inpjet,npjet
        write(151,'(i0,1x,l1,10(1x,es24.16))')i,jetfr(i),jetxx(i), &
         jetyy(i),jetzz(i),jetst(i),jetvx(i),jetvy(i),jetvz(i), &
         jetms(i),jetch(i),jetvl(i)
      enddo
    endif
    
!   check if the filament is breaking up between two beads
    call profiling_start(prof_breakup)
    call ckeck_breakup(mytime)
    call profiling_stop(prof_breakup)

!   Synchronize the complete state only for output formats that consume bead
!   arrays. Normal statistical output transfers only its selected bead.
#ifdef _OPENACC
    if(accelerator_is_persistent())then
      lfullhostoutput=.false.
      if(lprintxyz)lfullhostoutput=mod(nstep,iprintxyz)==0
      if(lprintxyzsing)lfullhostoutput=lfullhostoutput .or. &
       mod(nstep,iprintxyzsing)==0
      if(lprintpdbsing)lfullhostoutput=lfullhostoutput .or. &
       mod(nstep,iprintpdbsing)==0
      if(lprintdat)lfullhostoutput=lfullhostoutput .or. &
       mod(nstep,iprintdat)==0
      lfullhostoutput=lfullhostoutput .or. mod(nstep,nrestartdump)==0
      if(lfullhostoutput)then
        call accelerator_update_host_state(npjet,jetxx,jetyy,jetzz, &
         jetst,jetvx,jetvy,jetvz,nstep)
      endif
    endif
#endif
    
!   compute statistical quanities
    call profiling_start(prof_statistics)
    call statistic_driver(mytime,tstep,nstep,nremoved,ladd,lrem)
    call profiling_stop(prof_statistics)
    
!   print on the binary file the jet bead which have hit the collector 
!   (only for developers)
    call profiling_start(prof_output)
    call write_datrem_frame(lprintdat,122,nstep,mytime,iprintdat, &
     inpjet,npjet,sprintdat,systype,lremdat,nremoved)
    call profiling_stop(prof_output)
    
!   erase the bead beyond the collector if the variable lrem is .true.
    call profiling_start(prof_erase_bead)
    call erase_jetbead(nstep,mytime,ladd,lrem,nremoved)
#ifdef _OPENACC
    if(accelerator_is_persistent() .and. linserting .and. nremoved>0)then
      call accelerator_update_device_removed(inpjet-nremoved,inpjet-1, &
       jetst,jetvx,jetvy,jetvz,jetms,jetch,jetvl)
    endif
#endif
    call profiling_stop(prof_erase_bead)
    
!   print data on terminal and output 'statdat.dat' file
    call profiling_start(prof_output)
    call outprint_driver(nstep,mytime)
    
!   print the jet geometry on the XYZ formatted output file
    call write_xyz_frame(lprintxyz,120,nstep,mytime,iprintxyz, &
     inpjet,npjet,xyzrescale,systype,linserted,maxnumxyz,.false.)
   
!   print the jet geometry on the XYZ formatted output file
    call write_xyz_singlefile('frame',lprintxyzsing,125,nstep,mytime, &
     iprintxyzsing,inpjet,npjet,xyzrescale,systype,linserted,.false.)
     
!   print the jet geometry on the PDB formatted output file
    call write_pdb_singlefile('frame',lprintpdbsing,126,nstep,mytime, &
     iprintpdbsing,inpjet,npjet,pdbrescale,systype,linserted,.false.)
     
!   print the jet geometry on the binary file (only for developers)
    call write_dat_frame(lprintdat,130,nstep,mytime,iprintdat, &
     inpjet,npjet,sprintdat,systype,linserted)
    call profiling_stop(prof_output)
     
!   print restart file
    call profiling_start(prof_restart)
    call write_restart_file(nrestartdump,135,'save.dat',nstep,mytime)
    call profiling_stop(prof_restart)
    
!   cycle time check
    call time_world(ctime)
    
!   check recycle loop
    lrecycle=(lrecycle .and. timjob-ctime>timcls)
    
  enddo

  call get_sync_world()
  call wall_time_world(loop_end_time)
  loop_elapsed_time=loop_end_time-loop_start_time
  if(idrank==0)then
    write(6,'(/,a,f14.6,a)') &
     'Time-integration loop wall time: ',loop_elapsed_time,' s'
    if(loop_elapsed_time>0.d0)then
      write(6,'(a,f14.3,a/)')'Time-integration throughput: ', &
       dble(nstep)/loop_elapsed_time,' steps/s'
    endif
  endif
  call profiling_report(loop_elapsed_time,idrank)
  if(idrank==0)then
    write(6,'(a,i0)')'Topology additions: ',topology_add_total
    write(6,'(a,i0)')'Topology removals: ',topology_remove_total
    write(6,'(a,i0)')'Topology active beads: ',npjet-inpjet
  endif
!***********************************************************************
!     end of the time integration
!***********************************************************************
  
! stop clock
  call time_world(ftime)
  
! print last record on terminal and close the output 'statdat.dat' file
  call finish_print(nstep,mytime,itime,ftime)
  
! print restart file
#ifdef _OPENACC
  call accelerator_update_host_state(npjet,jetxx,jetyy,jetzz,jetst, &
   jetvx,jetvy,jetvz,nstep)
#endif
  call write_restart_file(1,135,'save.dat',nstep,mytime)
  
! close the XYZ formatted output file 
  call close_xyz_file(lprintxyz,120)

! close the binary file for hit beads (only for developers) 
  call close_datrem_file(lprintdat,122)
    
! close the binary file (only for developers) 
  call close_dat_file(lprintdat,130)
  if(ltopologysnapshot .and. idrank==0)close(151)
  
! close the communications
  call finalize_world()
  
  stop

 end program JetSpin
  
