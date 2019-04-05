module mod_dynam2
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /dynam2/ sdlims,sdts,plrps
  real,public,allocatable :: sdlims(:,:,:)
  real,public,allocatable :: sdts(:,:,:)
  real,public,allocatable :: plrps(:,:,:)
  !
  ! save /dynam2/
  ! common /dynam2a/ schisq1,schisq2,schisq3,schisq4,schisq5
  ! save /dynam2a/
  !
  real,public,allocatable :: schisq1(:),schisq2(:),schisq3(:),schisq4(:),schisq5(:)

  public :: allocate_mod_dynam2,deallocate_mod_dynam2

contains

  subroutine allocate_mod_dynam2
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam2','ALLOCATE')

    call allocate_array(sdlims,1,maxnks,1,maxnrs,-1,maxizs,'sdlims',ierr)
    call allocate_array(sdts,1,maxnks,1,maxnrs,-1,maxizs,'sdts',ierr)
    call allocate_array(plrps,1,maxnks,1,maxnrs,-1,maxplrp,'plrps',ierr)
    call allocate_array(schisq1,maxpiniter,'schisq1',ierr)
    call allocate_array(schisq2,maxpiniter,'schisq2',ierr)
    call allocate_array(schisq3,maxpiniter,'schisq3',ierr)
    call allocate_array(schisq4,maxpiniter,'schisq4',ierr)
    call allocate_array(schisq5,maxpiniter,'schisq5',ierr)

  end subroutine allocate_mod_dynam2


  subroutine deallocate_mod_dynam2
    implicit none

    call pr_trace('mod_dynam2','DEALLOCATE')

    if (allocated(sdlims)) deallocate(sdlims)
    if (allocated(sdts)) deallocate(sdts)
    if (allocated(plrps)) deallocate(plrps)
    if (allocated(schisq1)) deallocate(schisq1)
    if (allocated(schisq2)) deallocate(schisq2)
    if (allocated(schisq3)) deallocate(schisq3)
    if (allocated(schisq4)) deallocate(schisq4)
    if (allocated(schisq5)) deallocate(schisq5)

  end subroutine deallocate_mod_dynam2

end module mod_dynam2