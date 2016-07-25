

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
!                        JETSPIN VERSION 1.21
!
! (March 2016)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use version_mod,    only : init_world,get_rank_world,get_size_world,&
                       alloc_domain,time_world,time_world,finalize_world
  use utility_mod,    only : init_random_seed
  use nanojet_mod,    only : inpjet,npjet,linserted,myseed,systype, &
                       tstep,xyzrescale,set_resolution_length, &
                       allocate_jet,set_initial_jet,add_jetbead, &
                       remove_jetbead,erase_jetbead,lengthscale, &
                       pdbrescale,lreadrest
  use breaking_mod,   only : ckeck_breakup
  use dynamic_refinement_mod, only : refinementthreshold, &
                               set_refinement_threshold, &
                               refbeadstartfit
  use integrator_mod, only : initime,endtime,driver_integrator
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
                       iprintpdbsing
  
  implicit none
  
  integer:: nstep
  integer :: nremoved
  
  double precision :: mytime
  double precision :: itime,ftime
  
  logical :: ladd,lrem,lremdat,ldorefinment
  
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
  
!***********************************************************************
!     start the time integration
!***********************************************************************
  do while ((dble(nstep)*tstep)<endtime)
  
!   update the counter
    nstep=nstep+1
    
!   integrate the system
    call driver_integrator(mytime,tstep,nstep,ldorefinment)
    
!   check if a new bead should be added and/or removed
    call add_jetbead(nstep,mytime,ladd)
    call remove_jetbead(nstep,nremoved,mytime,lrem,lremdat)
    
!   check if the filament is breaking up between two beads
    call ckeck_breakup(mytime)
    
!   compute statistical quanities
    call statistic_driver(mytime,tstep,nstep,nremoved,ladd,lrem)
    
!   print on the binary file the jet bead which have hit the collector 
!   (only for developers)
    call write_datrem_frame(lprintdat,122,nstep,mytime,iprintdat, &
     inpjet,npjet,sprintdat,systype,lremdat,nremoved)
    
!   erase the bead beyond the collector if the variable lrem is .true.
    call erase_jetbead(nstep,mytime,ladd,lrem,nremoved)
    
!   print data on terminal and output 'statdat.dat' file
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
     
!   print restart file
    call write_restart_file(1000,135,'save.dat',nstep,mytime)
    
  enddo
!***********************************************************************
!     end of the time integration
!***********************************************************************
  
! stop clock
  call time_world(ftime)
  
! print last record on terminal and close the output 'statdat.dat' file
  call finish_print(nstep,mytime,itime,ftime)
  
! print restart file
  call write_restart_file(1,135,'save.dat',nstep,mytime)
  
! close the XYZ formatted output file 
  call close_xyz_file(lprintxyz,120)

! close the binary file for hit beads (only for developers) 
  call close_datrem_file(lprintdat,122)
    
! close the binary file (only for developers) 
  call close_dat_file(lprintdat,130)
  
  
! close the communications
  call finalize_world()
  
  stop

 end program JetSpin
  








