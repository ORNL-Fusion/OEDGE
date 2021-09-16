module mod_cnoco
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /cnoco/  ptes,pnes,pdvols,pnzs,pradis,prates
  !
  ! save /cnoco/
  real,public,allocatable :: ptes(:),pnes(:),pdvols(:)
  real,public,allocatable :: pnzs(:,:,:),pradis(:,:,:,:)
  real,public,allocatable :: prates(:,:,:,:)

  public :: allocate_mod_cnoco,deallocate_mod_cnoco

contains

  subroutine allocate_mod_cnoco
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cnoco','ALLOCATE')

    call allocate_array(ptes,maxnks,'ptes',ierr)
    call allocate_array(pnes,maxnks,'pnes',ierr)
    call allocate_array(pdvols,maxnks,'pdvols',ierr)
    call allocate_array(pnzs,29,3,maxnks,'pnzs',ierr)
    call allocate_array(pradis,1,14,1,29,1,3,1,maxnks,'pradis',ierr)
    call allocate_array(prates,1,4,1,28,1,3,1,maxnks,'prates',ierr)

  end subroutine allocate_mod_cnoco


  subroutine deallocate_mod_cnoco
    implicit none

    call pr_trace('mod_cnoco','DEALLOCATE')

    if (allocated(ptes)) deallocate(ptes)
    if (allocated(pnes)) deallocate(pnes)
    if (allocated(pdvols)) deallocate(pdvols)
    if (allocated(pnzs)) deallocate(pnzs)
    if (allocated(pradis)) deallocate(pradis)
    if (allocated(prates)) deallocate(prates)

  end subroutine deallocate_mod_cnoco

end module mod_cnoco