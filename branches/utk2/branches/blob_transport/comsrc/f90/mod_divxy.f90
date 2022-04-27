module mod_divxy
  use debug_options
  implicit none

  !
  !     this block contains the divimp declarations of ikxys ...
  !     the definitions between out and divimp were allowed to
  !     diverge so that the memory overhead that these variables
  !     incur (which can be substantial for high-resolution arrays)
  !     will not impact both programs. it is important to ensure
  !     that maxgxs and maxgys are greater than or equal to
  !     the corresponding maxnxs and maxnys.
  !
  !     -*-fortran-*-
  ! common /divxy/ ikxys,irxys,ifxys
  !
  ! save /divxy/
  integer,public,allocatable :: ikxys(:,:),irxys(:,:),ifxys(:,:)

  public :: allocate_mod_divxy,deallocate_mod_divxy

contains

  subroutine allocate_mod_divxy
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_divxy','ALLOCATE')

    call allocate_array(ikxys,maxnxs,maxnys,'ikxys',ierr)
    call allocate_array(irxys,maxnxs,maxnys,'irxys',ierr)
    call allocate_array(ifxys,maxnxs,maxnys,'ifxys',ierr)

  end subroutine allocate_mod_divxy


  subroutine deallocate_mod_divxy
    implicit none

    call pr_trace('mod_divxy','DEALLOCATE')

    if (allocated(ikxys)) deallocate(ikxys)
    if (allocated(irxys)) deallocate(irxys)
    if (allocated(ifxys)) deallocate(ifxys)

  end subroutine deallocate_mod_divxy

end module mod_divxy