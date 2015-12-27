 
 module error_mod
 
!***********************************************************************
!     
!     JETSPIN module containing subroutines which print warning and
!     close the software if an error is occurred
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification May 2015
!     
!***********************************************************************
 
 use version_mod, only : idrank,abort_world
 
 implicit none
 
 private
 
 public :: error,warning
 
 contains
 
 subroutine error(kode)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing error banners and close the
!     program
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in) :: kode
  
  integer,parameter :: outp=6
  character(len=*),parameter :: outf='(a)'
  character(len=*),parameter :: outf2='(2a)'
  
  if(idrank==0)then
     select case (kode)
      case (1)
        write(outp,outf)'ERROR - ktype error.'
      case (2)
        write(outp,outf)'ERROR - ftype is wrong!'
      case (3)
        write(outp,outf) &
        'ERROR - wrong number of variables for the selected function!'
      case (4)
        write(outp,outf)'ERROR - ftype is wrong for this integrator.'
      case (5)
        write(outp,outf)'ERROR - error in reading input file'
        write(outp,outf2)'ERROR - none FINISH directive found ', &
         'in the input file.'
      case (6)
        write(outp,outf)'ERROR - unknown directive in input file.'
      case (7)
        write(outp,outf)'ERROR - incomplete input file.'
      case (8)
        write(outp,outf2)'ERROR - the resolution of jet ', &
         'discretization is too small.'
      case (9)
        write(outp,outf2)'ERROR - wrong selection of print dat ', &
        'style in input file.'
      case (10)
        write(outp,outf)'ERROR - input file named input.dat not found.'
      case (11)
        write(outp,outf)'ERROR - input file opened with error.'
      case (12)
        write(outp,outf2)'ERROR - requested a spline fitting beyond ', &
         'the spline boundaries.'
      case (13)
        write(outp,outf2)'ERROR - in allocating jetptc in ', &
         'allocate_array_jetptc_mod.'
      case (14)
        write(outp,outf2)'ERROR - numerical instability.', &
         ' Please check the time step.'
      case (15)
        write(outp,outf)'ERROR - restart.dat file not found!'
      case (16)
        write(outp,outf2)'ERROR - dynamic refinement threshold too small!', &
         ' Please increase the dynamic refinement threshold.'
      case default
        write(outp,'(a,i18)')'unknown ERROR! code = ',kode
    end select
  endif
  
  call abort_world()
  stop
  
 end subroutine error
 
 subroutine warning(kode,ddata,wstring)
 
!***********************************************************************
!     
!     JETSPIN subroutine for printing warning banners
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2015
!     
!***********************************************************************
  
  implicit none
 
  integer, intent(in) :: kode
  double precision, intent(in), optional :: ddata
  character(len=*), intent(in), optional :: wstring
  
  integer,parameter :: outp=6
  character(len=*),parameter :: outf='(a)'
  character(len=10) :: r_char
  
  if(idrank/=0)return
  
  select case (kode)
    case (1)
      write(outp,'(/,a,/,a,g20.10,/)')"WARNING - the 'resolution' value is not a submultiple of the 'initial length'", &
      "WARNING - the remaider is equal to",ddata
    case (2)
      write(outp,'(/,a,/)')"WARNING - the 'resolution' value is redefined as a submultiple of the 'initial length'"
    case (3)
      write(outp,'(/,a,g20.10,/)') &
      "WARNING - the new 'resolution' value is equal to ",ddata
    case (4)
      write(outp,'(/,a,/)')"WARNING - 'systype' not specified in the input file"
    case (5)
      write(outp,'(/,a,/)')"WARNING - 'integrator' not specified in the input file"
    case (6)
      write(outp,'(/,a,/)')"WARNING - 'printlist' not specified in the input file"
    case (7)
      write(outp,'(/,a,/)')"WARNING - 'resolution' or 'points' not specified in the input file"
    case (8)
      write(outp,'(/,a,/)')"WARNING - 'initial lenght' not specified in the input file"
    case (9)
      write(outp,'(/,a,/)')"WARNING - 'mass density' not specified in the input file"
    case (10)
      write(outp,'(/,a,/)')"WARNING - 'charge density' not specified in the input file"
    case (11)
      write(outp,'(/,a,/)')"WARNING - 'nozzle cross section' not specified in the input file"
    case (12)
      write(outp,'(/,a,/)')"WARNING - 'nozzle stress' not specified in the input file"
    case (13)
      write(outp,'(/,a,/)')"WARNING - 'nozzle velocity' not specified in the input file"
    case (14)
      write(outp,'(/,a,/)')"WARNING - 'viscosity' not specified in the input file"
    case (15)
      write(outp,'(/,a,/)')"WARNING - 'elastic modulus' not specified in the input file"
    case (16)
      write(outp,'(/,a,/)')"WARNING - 'collector distance' not specified in the input file"
    case (17)
      write(outp,'(/,a,/)')"WARNING - 'external potential' not specified in the input file"
    case (18)
      write(outp,'(/,a,/)')"WARNING - 'surface tension' not specified in the input file"
    case (19)
      write(outp,'(/,a,/)')"WARNING - 'timestep' not specified in the input file"
    case (20)
      write(outp,'(/,a,/)')"WARNING - 'final time' of integration not specified in the input file"
    case (21)
      write (r_char,'(i10)')nint(ddata)
      if(nint(ddata)>1)then
      write(outp,'(/,3a,/)')"The jet is initially discretized in ",trim(adjustl(r_char)),' beads'
      else
      write(outp,'(/,3a,/)')"The jet is initially discretized in ",trim(adjustl(r_char)),' bead'
      endif
    case (22)
      write(outp,'(/,a)')"WARNING - 'resolution' and 'points' cannot be specified at the same time"
      write(outp,'(a,/)')"WARNING - only one of two should be used"
    case (23)
      write(outp,'(/,a,g20.10,a,/)') &
      "The 'resolution' value is equal to ",ddata, ' cm'
    case (24)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'nozzle stress' set equal to ",ddata,' g cm^-1 s^-2'
    case (25)
      write(outp,'(/,a,i18,/)')"WARNING - the number of beads is automatically set equal to ",nint(ddata)
    case (26)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'nozzle velocity' set equal to ",ddata," cm s^-1"
    case (27)
      write(outp,'(/,a,/)')"WARNING - 'yield stress' not specified in the input file"
    case (28)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'yield stress' set equal to ",ddata," g cm^-1 s^-2"
    case (29)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'initial lenght' is automatically set equal to ",ddata," cm"
    case (30)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'cutoff' is automatically set equal to ",ddata," cm"
    case (31)
      write(outp,'(/,a,/)')"WARNING - 'perturbation frequency' not specified in the input file"
    case (32)
      write(outp,'(/,a,/)')"WARNING - 'perturbation amplitude' not specified in the input file"
    case (33)
      write(outp,'(/,a,/)')"WARNING - 'airdrag amplitude' not specified in the input file"
    case (34)
      write(outp,'(/,a,i18,/)')"WARNING - 'airdrag yes' is not compatible with the actual 'integrator' value",nint(ddata)
    case (35)
      write(outp,'(/,a,i18,/)')"WARNING - 'airdrag yes' is not compatible with the actual 'systype' value",nint(ddata)
    case (36)
      write(outp,'(/,a,/)')"WARNING - 'inserting mode' not specified in the input file"
    case (37)
      write(outp,'(/,a)')"WARNING - there is the following unknown directive in the input file:"
      write(outp,'(3a,/)')"'",trim(wstring),"'"
    case (38)
      write(outp,'(/,a,/)')"WARNING - 'airdrag airdensity' not specified in the input file"
    case (39)
      write(outp,'(/,a,/)')"WARNING - 'airdrag airviscosity' not specified in the input file"
    case (40)
      write(outp,'(/,a,/)')"WARNING - 'HBfluid consistency' not specified in the input file"
    case (41)
      write(outp,'(/,a,/)')"WARNING - 'HBfluid index' not specified in the input file"
    case (42)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'initial lenght' is ",ddata," cm"
    case (43)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'nozzle cross' is ",ddata," cm"
    case (44)
      write(outp,'(/,a,/)')"WARNING - 'initial lenght' should be greater than 'nozzle cross'"
    case (45)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'perturbation amplitude' is ",ddata," cm"
    case (46)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'initial lenght' is ",ddata," cm"
    case (47)
      write(outp,'(/,a,/)')"WARNING - 'perturbation amplitude' should be lower than 'initial lenght'"
    case (48)
      write(outp,'(/,a)')"WARNING - 'print list' is not correctly specified"
      write(outp,'(a,/)')"WARNING - 'print list' should be specified as 'print list [keys]'"
    case (49)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'airdrag airvelocity' is automatically set equal to ",ddata," cm/s"
    case (50)
      write(outp,'(/,a,/)')"WARNING - 'variable mass deviation' not specified in the input file"
    case (51)
      write(outp,'(/,a,/)')"WARNING - 'variable mass correlation' not specified in the input file"
    case (52)
      write(outp,'(/,a,/)')"WARNING - 'variable charge deviation' not specified in the input file"
    case (53)
      write(outp,'(/,a,/)')"WARNING - 'variable charge correlation' not specified in the input file"
    case (54)
      write(outp,'(/,a)')"WARNING - 'variable mass deviation' should not be negative"
      write(outp,'(a,g20.10,/)')"WARNING - the actual value of 'variable mass deviation' is",ddata
    case (55)
      write(outp,'(/,a)')"WARNING - 'variable mass correlation' should not be negative"
      write(outp,'(a,g20.10,/)')"WARNING - the actual value of 'variable mass correlation' is",ddata
    case (56)
      write(outp,'(/,a)')"WARNING - 'variable charge deviation' should not be negative"
      write(outp,'(a,g20.10,/)')"WARNING - the actual value of 'variable charge deviation' is",ddata
    case (57)
      write(outp,'(/,a)')"WARNING - 'variable charge correlation' should not be negative"
      write(outp,'(a,g20.10,/)')"WARNING - the actual value of 'variable charge correlation' is",ddata
    case (58)
      write(outp,'(/,a)')"WARNING - 'printstat list' is not correctly specified"
      write(outp,'(a,/)')"WARNING - 'printstat list' should be specified as 'printstat list [keys]'"
    case (59)
      write(outp,'(/,a,g20.10,a,/)') &
      "WARNING - The 'print time' value is automatically set equal to ",ddata, ' s'
    case (60)
      write(outp,'(/,a)') &
      "WARNING - possible [keys] to be used are reported in the following table"
      write(outp,'(a)') &
      "WARNING - each [key] should be separated by a space character (e.g. t x y z rc)"
      write(outp,'(a)') &
      "WARNING - *********************************************************************"
      write(outp,'(a)') &
      "WARNING - * [keys]               * [meanings]                                 *"
      write(outp,'(a)') &
      "WARNING - *********************************************************************"
      write(outp,'(a)') &
      "WARNING - * t                    * unscaled time                              *"
      write(outp,'(a)') &
      "WARNING - * ts                   * scaled time                                *"
      write(outp,'(a)') &
      "WARNING - * x,y,z                * unscaled coordinates                       *"
      write(outp,'(a)') &
      "WARNING - * xs,ys,zs             * scaled coordinates                         *"
      write(outp,'(a)') &
      "WARNING - * st                   * unscaled stress                            *"
      write(outp,'(a)') &
      "WARNING - * sts                  * scaled stress                              *"
      write(outp,'(a)') &
      "WARNING - * vx,vy,vz             * unscaled velocities                        *"
      write(outp,'(a)') &
      "WARNING - * vxs,vys,vzs          * scaled velocities                          *"
      write(outp,'(a)') &
      "WARNING - * yz                   * unscaled normal distance from the x axis   *"
      write(outp,'(a)') &
      "WARNING - * yzs                  * scaled normal distance from the x axis     *"
      write(outp,'(a)') &
      "WARNING - * mass                 * last inserted mass at the nozzle           *"
      write(outp,'(a)') &
      "WARNING - * q                    * last inserted charge at the nozzle         *"
      write(outp,'(a)') &
      "WARNING - * cpu                  * time for every print interval              *"
      write(outp,'(a)') &
      "WARNING - * cpur                 * remaining time to the end                  *"
      write(outp,'(a)') &
      "WARNING - * cpue                 * elapsed time                               *"
      write(outp,'(a)') &
      "WARNING - * n                    * bead number of jet discretization          *"
      write(outp,'(a)') &
      "WARNING - * f                    * index of the first bead                    *"
      write(outp,'(a)') &
      "WARNING - * l                    * index of the last bead                     *"
      write(outp,'(a)') &
      "WARNING - * curn                 * current at the nozzle                      *"
      write(outp,'(a)') &
      "WARNING - * curc                 * current at the collector                   *"
      write(outp,'(a)') &
      "WARNING - * vn                   * velocity modulus of jet at the nozzle      *"
      write(outp,'(a)') &
      "WARNING - * vc                   * velocity modulus of jet at the collector   *"
      write(outp,'(a)') &
      "WARNING - * svc                  * strain velocity at the collector           *"
      write(outp,'(a)') &
      "WARNING - * mfn                  * mass flux at the nozzle                    *"
      write(outp,'(a)') &
      "WARNING - * mfc                  * mass flux at the collector                 *"
      write(outp,'(a)') &
      "WARNING - * rc                   * radius of jet at the collector             *"
      write(outp,'(a)') &
      "WARNING - * rrr                  * radius reduction ratio of jet              *"
      write(outp,'(a)') &
      "WARNING - * lp                   * length path of jet                         *"
      write(outp,'(a)') &
      "WARNING - * rlp                  * length path of jet / collector distance    *"
      write(outp,'(a)') &
      "WARNING - * nref                 * number of performed dynamic refinements    *"
      write(outp,'(a)') &
      "WARNING - * angl                 * instantaneus angular aperture              *"
      write(outp,'(a,/)') &
      "WARNING - *********************************************************************"
    case (61)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,a)')"WARNING - A line of input file is not correctly specified:"
      write(outp,'(4a,/)')'line ',trim(adjustl(r_char)),' : ',wstring
    case (62)
      write(outp,'(/,a)')"WARNING - The 'resolution' value is lower than 'nozzle cross'"
      write(outp,'(a,/)')"WARNING - 'resolution' should be greater than 'nozzle cross'"
    case (63)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,a)')"WARNING - wrong selection of 'print dat style' in input file"
      write(outp,'(2a,/)')"WARNING - the actual selection is : ",trim(adjustl(r_char))
    case (64)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'surface tension'  set equal to ",ddata," g s^-2"
    case (65)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'external potential' set equal to ",ddata," statV"
    case (66)
      write(outp,'(/,a,g20.10,a,/)')"WARNING - 'dynamic refinement threshold' set equal to ",ddata," cm"
    case (67)
      write (r_char,'(i10)')nint(ddata)
      write(outp,'(/,2a,/)')"WARNING - numerical instability at nstep ",trim(adjustl(r_char))
    case (68)
      write(outp,'(/,a,g20.10,a,/)') &
      "WARNING - 'dynamic refinement every' value automatically set equal to ",ddata, ' s'
    case (69)
      write(outp,'(/,2a,g20.10,a,/)')"WARNING - JETSPIN is restarting from 'restart.dat' file", &
      "WARNING - The simulation restarts at the simulation time ",ddata, ' s'
    case (70)
      write(outp,'(/,a,/)') &
       "WARNING - 'variable mass typemass 3' is availabe only in refinement mode activated"
    case default
      write(outp,'(/,a,i8,/)')"unknown WARNING! code = ",kode
  end select

  return

 end subroutine warning

 end module error_mod


