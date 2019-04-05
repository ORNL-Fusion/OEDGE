module mod_sol22phelpi
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /sol22phelpi/ lasts,lastphelp,lastsrc
  ! save /sol22phelpi/
  !
  real*8,public :: lasts,lastphelp,lastsrc

  public :: allocate_mod_sol22phelpi,deallocate_mod_sol22phelpi

contains

  subroutine allocate_mod_sol22phelpi
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_sol22phelpi','ALLOCATE')


  end subroutine allocate_mod_sol22phelpi


  subroutine deallocate_mod_sol22phelpi
    implicit none

    call pr_trace('mod_sol22phelpi','DEALLOCATE')


  end subroutine deallocate_mod_sol22phelpi

end module mod_sol22phelpi