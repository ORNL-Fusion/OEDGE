module mod_clocal
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /clocal/ lfps,lfss,lfts,lllfps,ltolds
  !
  ! save /clocal/
  real,public,allocatable :: ltolds(:,:,:),lfps(:,:,:),lfss(:,:,:),lfts(:,:,:),lllfps(:,:,:)

  public :: allocate_mod_clocal,deallocate_mod_clocal

contains

  subroutine allocate_mod_clocal
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_clocal','ALLOCATE')

    call allocate_array(ltolds,maxnks,maxnrs,maxizs,'ltolds',ierr)
    call allocate_array(lfps,maxnks,maxnrs,maxizs,'lfps',ierr)
    call allocate_array(lfss,maxnks,maxnrs,maxizs,'lfss',ierr)
    call allocate_array(lfts,maxnks,maxnrs,maxizs,'lfts',ierr)
    call allocate_array(lllfps,maxnks,maxnrs,maxizs,'lllfps',ierr)

  end subroutine allocate_mod_clocal


  subroutine deallocate_mod_clocal
    implicit none

    call pr_trace('mod_clocal','DEALLOCATE')

    if (allocated(ltolds)) deallocate(ltolds)
    if (allocated(lfps)) deallocate(lfps)
    if (allocated(lfss)) deallocate(lfss)
    if (allocated(lfts)) deallocate(lfts)
    if (allocated(lllfps)) deallocate(lllfps)

  end subroutine deallocate_mod_clocal

end module mod_clocal