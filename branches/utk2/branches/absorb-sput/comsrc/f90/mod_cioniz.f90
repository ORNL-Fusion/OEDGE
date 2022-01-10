module mod_cioniz
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /cioniz/ c215a,c350a,c350b,kfizs,kfrcs,kpchs,kprcs,kfcxs,knhs,kfps,kfss,&
  !     kfts,kkkfps,kplos,kfssmod,kpmtc,kpizs,lmfp_mtc
  !
  ! save /cioniz/
  real,public :: c215a,c350a,c350b
  real,public,allocatable :: kfizs(:,:,:),kfrcs(:,:,:),kpchs(:,:,:),kprcs(:,:,:),kfcxs(:,:,:),&
       knhs(:,:),kfps(:,:,:),kfss(:,:,:),kfts(:,:,:),kkkfps(:,:,:),kplos(:,:,:),&
       kfssmod(:,:),kpmtc(:,:),kpizs(:,:),lmfp_mtc(:,:)

  public :: allocate_mod_cioniz,deallocate_mod_cioniz

contains

  subroutine allocate_mod_cioniz
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cioniz','ALLOCATE')

    call allocate_array(kfizs,1,maxnks,1,maxnrs,0,maxizs,'kfizs',ierr)
    call allocate_array(kfrcs,1,maxnks,1,maxnrs,0,maxizs,'kfrcs',ierr)
    call allocate_array(kpchs,1,maxnks,1,maxnrs,0,maxizs,'kpchs',ierr)
    call allocate_array(kprcs,1,maxnks,1,maxnrs,0,maxizs,'kprcs',ierr)
    call allocate_array(kfcxs,1,maxnks,1,maxnrs,0,maxizs,'kfcxs',ierr)
    call allocate_array(knhs,maxnks,maxnrs,'knhs',ierr)
    call allocate_array(kfps,maxnks,maxnrs,maxizs,'kfps',ierr)
    call allocate_array(kfss,maxnks,maxnrs,maxizs,'kfss',ierr)
    call allocate_array(kfts,maxnks,maxnrs,maxizs,'kfts',ierr)
    call allocate_array(kkkfps,maxnks,maxnrs,maxizs,'kkkfps',ierr)
    call allocate_array(kplos,maxnks,maxnrs,maxizs,'kplos',ierr)
    call allocate_array(kfssmod,maxnks,maxnrs,'kfssmod',ierr)
    call allocate_array(kpmtc,maxnks,maxnrs,'kpmtc',ierr)
    call allocate_array(kpizs,maxnks,maxnrs,'kpizs',ierr)
    call allocate_array(lmfp_mtc,maxnks,maxnrs,'lmfp_mtc',ierr)

  end subroutine allocate_mod_cioniz


  subroutine deallocate_mod_cioniz
    implicit none

    call pr_trace('mod_cioniz','DEALLOCATE')

    if (allocated(kfizs)) deallocate(kfizs)
    if (allocated(kfrcs)) deallocate(kfrcs)
    if (allocated(kpchs)) deallocate(kpchs)
    if (allocated(kprcs)) deallocate(kprcs)
    if (allocated(kfcxs)) deallocate(kfcxs)
    if (allocated(knhs)) deallocate(knhs)
    if (allocated(kfps)) deallocate(kfps)
    if (allocated(kfss)) deallocate(kfss)
    if (allocated(kfts)) deallocate(kfts)
    if (allocated(kkkfps)) deallocate(kkkfps)
    if (allocated(kplos)) deallocate(kplos)
    if (allocated(kfssmod)) deallocate(kfssmod)
    if (allocated(kpmtc)) deallocate(kpmtc)
    if (allocated(kpizs)) deallocate(kpizs)
    if (allocated(lmfp_mtc)) deallocate(lmfp_mtc)

  end subroutine deallocate_mod_cioniz

end module mod_cioniz