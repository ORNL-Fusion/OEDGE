module mod_ccona
  use debug_options
  implicit none

  ! [dummy]

  public :: allocate_mod_ccona,deallocate_mod_ccona

contains

  subroutine allocate_mod_ccona
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_ccona','ALLOCATE')


  end subroutine allocate_mod_ccona


  subroutine deallocate_mod_ccona
    implicit none

    call pr_trace('mod_ccona','DEALLOCATE')


  end subroutine deallocate_mod_ccona

end module mod_ccona