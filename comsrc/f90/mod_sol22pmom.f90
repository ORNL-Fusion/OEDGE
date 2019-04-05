module mod_sol22pmom
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /sol22pmom/ lasts,lastsmom,lastsrc,lastv,lastte
  ! save /sol22pmom/
  !
  real*8 ,public :: lasts,lastsmom,lastsrc,lastv,lastte

  public :: allocate_mod_sol22pmom,deallocate_mod_sol22pmom

contains

  subroutine allocate_mod_sol22pmom
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_sol22pmom','ALLOCATE')


  end subroutine allocate_mod_sol22pmom


  subroutine deallocate_mod_sol22pmom
    implicit none

    call pr_trace('mod_sol22pmom','DEALLOCATE')


  end subroutine deallocate_mod_sol22pmom

end module mod_sol22pmom