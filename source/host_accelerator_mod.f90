module accelerator_mod

 implicit none
 private

 logical, parameter, public :: accelerator_enabled=.false.

 public :: accelerator_prepare

contains

 subroutine accelerator_prepare()
  implicit none
  return
 end subroutine accelerator_prepare

end module accelerator_mod
