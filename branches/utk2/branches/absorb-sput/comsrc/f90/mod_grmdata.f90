module mod_grmdata
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  integer,public :: ix,iy
  real,public :: xwid,ywid,xpt,ypt
  real,public,allocatable :: xpts(:),ypts(:)
  data    ypts /0.80,0.55,0.30,0.05/
  data    xpts /0.05, 0.55/
  data    xwid /0.40/
  !
  !     draw box - boxindex is 1 to 8
  !
  !     iy is 1 to 4
  !     ix is 1 to 2
  !
  data    ywid /0.15/
  ix = 1
  !
  if (float(boxindex/2).eq.float(boxindex)/2.0) ix = 2
  !
  iy =  (boxindex-1)/2 + 1
  xpt = xpts(ix)
  ypt = ypts(iy)

  public :: allocate_mod_grmdata,deallocate_mod_grmdata

contains

  subroutine allocate_mod_grmdata
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_grmdata','ALLOCATE')

    call allocate_array(xpts,2,'xpts',ierr)
    call allocate_array(ypts,4,'ypts',ierr)

  end subroutine allocate_mod_grmdata


  subroutine deallocate_mod_grmdata
    implicit none

    call pr_trace('mod_grmdata','DEALLOCATE')

    if (allocated(xpts)) deallocate(xpts)
    if (allocated(ypts)) deallocate(ypts)

  end subroutine deallocate_mod_grmdata

end module mod_grmdata