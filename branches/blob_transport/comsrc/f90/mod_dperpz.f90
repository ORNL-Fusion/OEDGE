module mod_dperpz
  use debug_options
  implicit none

  !
  !     common block for the dperpz transport option - a parallel displacement
  !     resulting from a perpendicular step along the paramagnetic axis.
  !
  !     -*-fortran-*-
  ! common /dperpz/ base_dperpz_step,dperpz_opt
  !
  ! save /dperpz/
  integer,public :: dperpz_opt
  !
  real,public :: base_dperpz_step

  public :: allocate_mod_dperpz,deallocate_mod_dperpz

contains

  subroutine allocate_mod_dperpz
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dperpz','ALLOCATE')


  end subroutine allocate_mod_dperpz


  subroutine deallocate_mod_dperpz
    implicit none

    call pr_trace('mod_dperpz','DEALLOCATE')


  end subroutine deallocate_mod_dperpz

end module mod_dperpz