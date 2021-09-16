module mod_cnoco


  !COMMON /CNOCO/  PTES,PNES,PDVOLS,PNZS,PRADIS,PRATES                       
  !REAL PTES(MAXNXS),PNES(MAXNXS),PDVOLS(MAXNXS)                             
  !REAL PNZS(29,3,MAXNXS),PRADIS(14,29,3,MAXNXS)                             
  !REAL PRATES(4,28,3,MAXNXS)                                                

  implicit none


  real,public,allocatable :: &
       PTES(:),PNES(:),PDVOLS(:),&
       PNZS(:,:,:),PRADIS(:,:,:,:),&
       PRATES(:,:,:,:)                                                

  private


  public :: allocate_mod_cnoco, deallocate_mod_cnoco


contains

  subroutine allocate_mod_cnoco
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call allocate_array(PTES,maxnxs,'PTES',ierr)
    call allocate_array(PNES,maxnxs,'PNES',ierr)
    call allocate_array(PDVOLS,maxnxs,'PDVOLS',ierr)
    call allocate_array(PNZS,29,3,maxnxs,'PNZS',ierr)
    call allocate_array(PRADIS,1,14,1,29,1,3,1,maxnxs,'PRADIS',ierr)
    call allocate_array(PRATES,1,4,1,28,1,3,1,maxnxs,'PRATES',ierr)

  end subroutine allocate_mod_cnoco


  subroutine deallocate_mod_cnoco
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(PTES)
    deallocate(PNES)
    deallocate(PDVOLS)
    deallocate(PNZS)
    deallocate(PRADIS)
    deallocate(PRATES)

  end subroutine deallocate_mod_cnoco



end module mod_cnoco
