module mod_comprt
  use debug_options
  implicit none

  ! [dummy]

  public :: allocate_mod_comprt,deallocate_mod_comprt

contains

  subroutine allocate_mod_comprt
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_comprt','ALLOCATE')


  end subroutine allocate_mod_comprt


  subroutine deallocate_mod_comprt
    implicit none

    call pr_trace('mod_comprt','DEALLOCATE')


  end subroutine deallocate_mod_comprt

end module mod_comprt