module mod_comxs
  use debug_options
  implicit none

  ! [dummy]

  public :: allocate_mod_comxs,deallocate_mod_comxs

contains

  subroutine allocate_mod_comxs
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_comxs','ALLOCATE')


  end subroutine allocate_mod_comxs


  subroutine deallocate_mod_comxs
    implicit none

    call pr_trace('mod_comxs','DEALLOCATE')


  end subroutine deallocate_mod_comxs

end module mod_comxs