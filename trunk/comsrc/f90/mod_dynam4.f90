module mod_dynam4
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /dynam4/ imode,nts,cstmax,dwelts,dwelfs,ctimes,lims,walks
  ! save /dynam4/
  integer,public :: imode,nts
  real,public :: cstmax
  real,public,allocatable :: dwelts(:),dwelfs(:),ctimes(:,:),lims(:,:,:,:),walks(:,:)

  public :: allocate_mod_dynam4,deallocate_mod_dynam4

contains

  subroutine allocate_mod_dynam4
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam4','ALLOCATE')

    call allocate_array(dwelts,-1,'dwelts',maxizs,ierr)
    call allocate_array(dwelfs,maxnts,'dwelfs',ierr)
    call allocate_array(ctimes,0,maxnts+1,-1,maxizs,'ctimes',ierr)
    call allocate_array(lims,1,maxnks,1,maxnrs,-1,maxizs,1,maxnts,'lims',ierr)
    call allocate_array(walks,maxnws,2,'walks',ierr)

  end subroutine allocate_mod_dynam4


  subroutine deallocate_mod_dynam4
    implicit none

    call pr_trace('mod_dynam4','DEALLOCATE')

    if (allocated(dwelts)) deallocate(dwelts)
    if (allocated(dwelfs)) deallocate(dwelfs)
    if (allocated(ctimes)) deallocate(ctimes)
    if (allocated(lims)) deallocate(lims)
    if (allocated(walks)) deallocate(walks)

  end subroutine deallocate_mod_dynam4

end module mod_dynam4