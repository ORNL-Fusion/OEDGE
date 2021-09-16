module mod_sol22pcx
  !use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /sol22pcx/ lasts,lastpcx,lastsrc
  ! save /sol22pcx/
  !
  real*8,public :: lasts,lastpcx,lastsrc

  !public :: allocate_mod_sol22pcx,deallocate_mod_sol22pcx

contains

  subroutine allocate_mod_sol22pcx
    !use mod_params
    !use allocate_arrays
    implicit none
    integer :: ierr

    !call pr_trace('mod_sol22pcx','ALLOCATE')


  end subroutine allocate_mod_sol22pcx


  subroutine deallocate_mod_sol22pcx
    implicit none

    !call pr_trace('mod_sol22pcx','DEALLOCATE')


  end subroutine deallocate_mod_sol22pcx

end module mod_sol22pcx
