module mod_plot_switches
  use debug_options
  implicit none

  !
  ! jdemod
  !
  !     the following common contains a flag to be sent to
  !     the contour drawing routines to instruct them to use
  !     .ge. and .le. for the contour boundaries instead of
  !     .gt. and .le. - this will allow the lower boundary to
  !     be properly drawn for contour plots.
  !
  !     -*-fortran-*-
  ! common /cont_data/ first_contour
  ! save /cont_data/
  !
  ! jdemod
  !
  logical,public :: first_contour

  public :: allocate_mod_plot_switches,deallocate_mod_plot_switches

contains

  subroutine allocate_mod_plot_switches
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_plot_switches','ALLOCATE')


  end subroutine allocate_mod_plot_switches


  subroutine deallocate_mod_plot_switches
    implicit none

    call pr_trace('mod_plot_switches','DEALLOCATE')


  end subroutine deallocate_mod_plot_switches

end module mod_plot_switches