module mod_colours
  use debug_options
  implicit none

  !
  !     colour definitions for use in out
  !
  !
  !     -*-fortran-*-
  integer,public :: maxncols,maxmarkers
  parameter(maxncols=255,maxmarkers=8)
  !
  
  ! common /colours/ colour,ncols,defcol,icol,start_col
  !
  ! save /colours/
  !
  !     set up selected markers for plotting
  !
  integer,public :: ncols,defcol,icol,col,start_col
  integer,public,allocatable :: colour(:)
  !
  integer,public,allocatable :: marker(:)
  !
  !data marker /232,250,224,225,227,248,228,229/

  public :: allocate_mod_colours,deallocate_mod_colours

contains

  subroutine allocate_mod_colours
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_colours','ALLOCATE')

    call allocate_array(colour,maxncols,'colour',ierr)
    call allocate_array(marker,maxmarkers,'marker',ierr)

    marker = [232,250,224,225,227,248,228,229]

  end subroutine allocate_mod_colours


  subroutine deallocate_mod_colours
    implicit none

    call pr_trace('mod_colours','DEALLOCATE')

    if (allocated(colour)) deallocate(colour)
    if (allocated(marker)) deallocate(marker)

  end subroutine deallocate_mod_colours

end module mod_colours
