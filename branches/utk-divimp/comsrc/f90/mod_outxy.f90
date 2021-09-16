module mod_outxy
  use debug_options
  implicit none

  !
  !     this common block contains the subsidiary xy array
  !     declarations for use in the out program only. in order
  !     to reduce the memory requirements for divimp these are
  !     set up in a separate common block.
  !
  !     -*-fortran-*-
  ! common /grxy/ ikxys,irxys,ifxys
  !
  ! save /grxy/
  integer,public,allocatable :: ikxys(:,:),irxys(:,:),ifxys(:,:)
  
  
  
  
  
  

  public :: allocate_mod_outxy,deallocate_mod_outxy

contains

  subroutine allocate_mod_outxy
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_outxy','ALLOCATE')

    call allocate_array(ikxys,maxixs,maxiys,'ikxys',ierr)
    call allocate_array(irxys,maxixs,maxiys,'irxys',ierr)
    call allocate_array(ifxys,maxixs,maxiys,'ifxys',ierr)

  end subroutine allocate_mod_outxy


  subroutine deallocate_mod_outxy
    implicit none

    call pr_trace('mod_outxy','DEALLOCATE')

    if (allocated(ikxys)) deallocate(ikxys)
    if (allocated(irxys)) deallocate(irxys)
    if (allocated(ifxys)) deallocate(ifxys)

  end subroutine deallocate_mod_outxy

end module mod_outxy