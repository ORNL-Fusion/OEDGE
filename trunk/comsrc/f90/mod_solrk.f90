module mod_solrk
  use debug_options
  implicit none

  !
  !     define coefficients
  !
  !     -*-fortran-*-
  ! common /solrk/  bij,ai,ci,cip
  !
  ! save /solrk/
  !
  real*8,public,allocatable :: bij(:,:),ai(:),ci(:),cip(:)

  public :: allocate_mod_solrk,deallocate_mod_solrk

contains

  subroutine allocate_mod_solrk
    !use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_solrk','ALLOCATE')

    call allocate_array(bij,6,5,'bij',ierr)
    call allocate_array(ai,6,'ai',ierr)
    call allocate_array(ci,6,'ci',ierr)
    call allocate_array(cip,6,'cip',ierr)

  end subroutine allocate_mod_solrk


  subroutine deallocate_mod_solrk
    implicit none

    call pr_trace('mod_solrk','DEALLOCATE')

    if (allocated(bij)) deallocate(bij)
    if (allocated(ai)) deallocate(ai)
    if (allocated(ci)) deallocate(ci)
    if (allocated(cip)) deallocate(cip)

  end subroutine deallocate_mod_solrk

end module mod_solrk
