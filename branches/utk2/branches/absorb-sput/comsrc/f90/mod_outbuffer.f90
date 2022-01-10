module mod_outbuffer
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /outbuffer/ comment
  ! save /outbuffer/
  !
  character*80,public :: comment

  public :: allocate_mod_outbuffer,deallocate_mod_outbuffer

contains

  subroutine allocate_mod_outbuffer
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_outbuffer','ALLOCATE')


  end subroutine allocate_mod_outbuffer


  subroutine deallocate_mod_outbuffer
    implicit none

    call pr_trace('mod_outbuffer','DEALLOCATE')


  end subroutine deallocate_mod_outbuffer

end module mod_outbuffer