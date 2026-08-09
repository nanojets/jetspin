module profiling_mod

 implicit none
 private

 integer, parameter, public :: prof_integrator=1
 integer, parameter, public :: prof_coulomb=2
 integer, parameter, public :: prof_add_bead=3
 integer, parameter, public :: prof_remove_bead=4
 integer, parameter, public :: prof_breakup=5
 integer, parameter, public :: prof_statistics=6
 integer, parameter, public :: prof_erase_bead=7
 integer, parameter, public :: prof_output=8
 integer, parameter, public :: prof_restart=9
 integer, parameter :: nprofile=9

 logical, public, save :: profiling_enabled=.false.
 double precision, save :: elapsed(nprofile)=0.d0
 integer, save :: starts(nprofile)=0
 integer, save :: calls(nprofile)=0
 integer, save :: clock_rate=0

 public :: profiling_initialize
 public :: profiling_reset
 public :: profiling_start
 public :: profiling_stop
 public :: profiling_report

contains

 subroutine profiling_initialize()
  implicit none
  character(len=32) :: value
  integer :: status,length

  value=''
  call get_environment_variable('JETSPIN_PROFILE',value,length,status)
  if(status==0 .and. length>0)then
    select case(value(1:length))
      case('1','yes','YES','true','TRUE','on','ON')
        profiling_enabled=.true.
      case default
        profiling_enabled=.false.
    end select
  endif
  call system_clock(count_rate=clock_rate)
  call profiling_reset()
 end subroutine profiling_initialize

 subroutine profiling_reset()
  implicit none
  elapsed(:)=0.d0
  starts(:)=0
  calls(:)=0
 end subroutine profiling_reset

 subroutine profiling_start(timer_id)
  implicit none
  integer, intent(in) :: timer_id
  if(.not.profiling_enabled)return
  call system_clock(starts(timer_id))
 end subroutine profiling_start

 subroutine profiling_stop(timer_id)
  implicit none
  integer, intent(in) :: timer_id
  integer :: finish
  if(.not.profiling_enabled)return
  call system_clock(finish)
  elapsed(timer_id)=elapsed(timer_id)+ &
   dble(finish-starts(timer_id))/dble(clock_rate)
  calls(timer_id)=calls(timer_id)+1
 end subroutine profiling_stop

 subroutine profiling_report(loop_time,idrank)
  implicit none
  double precision, intent(in) :: loop_time
  integer, intent(in) :: idrank
  integer :: i
  character(len=24), parameter :: labels(nprofile) = &
   [character(len=24) :: &
   'Integrator total        ', 'Coulomb (nested)       ', &
   'Add bead               ', 'Remove bead            ', &
   'Breakup check          ', 'Statistics             ', &
   'Erase bead             ', 'Scheduled output       ', &
   'Restart output         ']

  if((.not.profiling_enabled) .or. idrank/=0)return
  write(6,'(/,a)')'Subroutine timing profile (wall clock)'
  write(6,'(a)')'--------------------------------------'
  write(6,'(a24,2x,a10,2x,a8,2x,a10)')'Region','Seconds','Loop %','Calls'
  do i=1,nprofile
    if(loop_time>0.d0)then
      write(6,'(a24,2x,f10.6,2x,f8.2,2x,i10)')labels(i),elapsed(i), &
       100.d0*elapsed(i)/loop_time,calls(i)
    else
      write(6,'(a24,2x,f10.6,2x,f8.2,2x,i10)')labels(i),elapsed(i), &
       0.d0,calls(i)
    endif
  enddo
  write(6,'(a)') &
   'Note: Coulomb is included within Integrator total and is not additive.'
 end subroutine profiling_report

end module profiling_mod
