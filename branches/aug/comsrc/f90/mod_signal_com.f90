module mod_signal_com
  use debug_options
  implicit none

  !
  !     signal_com: contains information shared between the signal
  !     writing routines in out
  !
  !     -*-fortran-*-
  integer,parameter ,public :: signal_unit = 27

  public :: allocate_mod_signal_com,deallocate_mod_signal_com

contains

  subroutine allocate_mod_signal_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_signal_com','ALLOCATE')


  end subroutine allocate_mod_signal_com


  subroutine deallocate_mod_signal_com
    implicit none

    call pr_trace('mod_signal_com','DEALLOCATE')


  end subroutine deallocate_mod_signal_com

end module mod_signal_com