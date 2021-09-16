module mod_crand
  use debug_options
  implicit none

  !     -*-fortran-*-
  ! common /crand/ ranv
  ! save /crand/
  real,public,allocatable :: ranv(:)

  public :: allocate_mod_crand,deallocate_mod_crand

contains

  subroutine allocate_mod_crand
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_crand','ALLOCATE')

    call allocate_array(ranv,1000*isect,'ranv',ierr)

  end subroutine allocate_mod_crand


  subroutine deallocate_mod_crand
    implicit none

    call pr_trace('mod_crand','DEALLOCATE')

    if (allocated(ranv)) deallocate(ranv)

  end subroutine deallocate_mod_crand

end module mod_crand