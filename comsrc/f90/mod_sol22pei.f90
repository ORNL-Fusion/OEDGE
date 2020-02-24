module mod_sol22pei
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /sol22pei/ lastpei,lasts
  ! save /sol22pei/
  !
  real*8,public :: lasts,lastpei

  public :: allocate_mod_sol22pei,deallocate_mod_sol22pei

contains

  subroutine allocate_mod_sol22pei
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_sol22pei','ALLOCATE')


  end subroutine allocate_mod_sol22pei


  subroutine deallocate_mod_sol22pei
    implicit none

    call pr_trace('mod_sol22pei','DEALLOCATE')


  end subroutine deallocate_mod_sol22pei

end module mod_sol22pei