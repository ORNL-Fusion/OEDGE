module mod_comhr
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /hrinfo/  soltehr,soltihr,solnehr,solcorhr,cind,maxind
  !
  ! save /hrinfo/
  real,public,allocatable :: soltehr(:),soltihr(:),solnehr(:),solcorhr(:)
  integer,public :: cind,maxind

  public :: allocate_mod_comhr,deallocate_mod_comhr

contains

  subroutine allocate_mod_comhr
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_comhr','ALLOCATE')

    call allocate_array(soltehr,0,'soltehr',msolpt*maxnks+msolpt+1,ierr)
    call allocate_array(soltihr,0,'soltihr',msolpt*maxnks+msolpt+1,ierr)
    call allocate_array(solnehr,0,'solnehr',msolpt*maxnks+msolpt+1,ierr)
    call allocate_array(solcorhr,0,'solcorhr',msolpt*maxnks+msolpt+1,ierr)

  end subroutine allocate_mod_comhr


  subroutine deallocate_mod_comhr
    implicit none

    call pr_trace('mod_comhr','DEALLOCATE')

    if (allocated(soltehr)) deallocate(soltehr)
    if (allocated(soltihr)) deallocate(soltihr)
    if (allocated(solnehr)) deallocate(solnehr)
    if (allocated(solcorhr)) deallocate(solcorhr)

  end subroutine deallocate_mod_comhr

end module mod_comhr