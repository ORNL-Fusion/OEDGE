module mod_local_baffles
  use debug_options
  implicit none

  !
  !     this common block contains baffle related data local to
  !     the baffle routines
  !
  !     -*-fortran-*-
  ! common /local_baffles/ rbufxl,zbufxl,nbufmxl,nbufxl
  !
  ! save /local_baffles/
  integer,public :: nbufmxl
  integer,public,allocatable :: nbufxl(:)
  !
  real*8,public,allocatable :: rbufxl(:,:),zbufxl(:,:)

  public :: allocate_mod_local_baffles,deallocate_mod_local_baffles

contains

  subroutine allocate_mod_local_baffles
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_local_baffles','ALLOCATE')

    call allocate_array(nbufxl,mbufx+2,'nbufxl',ierr)
    call allocate_array(rbufxl,mbufx+2,mves,'rbufxl',ierr)
    call allocate_array(zbufxl,mbufx+2,mves,'zbufxl',ierr)

  end subroutine allocate_mod_local_baffles


  subroutine deallocate_mod_local_baffles
    implicit none

    call pr_trace('mod_local_baffles','DEALLOCATE')

    if (allocated(nbufxl)) deallocate(nbufxl)
    if (allocated(rbufxl)) deallocate(rbufxl)
    if (allocated(zbufxl)) deallocate(zbufxl)

  end subroutine deallocate_mod_local_baffles

end module mod_local_baffles