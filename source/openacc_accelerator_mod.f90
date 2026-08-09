module accelerator_mod

 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.true.

 public :: accelerator_prepare

contains

 subroutine accelerator_prepare()
  implicit none
!$acc init
  return
 end subroutine accelerator_prepare

end module accelerator_mod
