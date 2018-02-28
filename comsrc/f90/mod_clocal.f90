module mod_clocal

  use debug_options

  implicit none

  private

  !include 'clocal'

!     -*-Fortran-*-
!                                                                       
!      COMMON /CLOCAL/ LFPS,LFSS,LFTS,LLLFPS,LTOLDS                      
!      save /clocal/
!
!      REAL LTOLDS(MAXNKS,MAXNRS,MAXIZS),                                
!     >  LFPS (MAXNKS,MAXNRS,MAXIZS)    ,LFSS  (MAXNKS,MAXNRS,MAXIZS),   
!     >  LFTS (MAXNKS,MAXNRS,MAXIZS)    ,LLLFPS(MAXNKS,MAXNRS,MAXIZS)    
!


  real, allocatable,  public :: ltolds(:,:,:),lfps(:,:,:),lfss(:,:,:),lfts(:,:,:),lllfps(:,:,:)

  public :: allocate_clocal,deallocate_clocal

contains

  subroutine allocate_clocal
    use global_parameters
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('MOD_CLOCAL','ALLOCATE')

    call allocate_array(ltolds,maxnks,maxnrs,maxizs,'LTOLDS',ierr)
    call allocate_array(lfps,maxnks,maxnrs,maxizs,'LFPS',ierr)
    call allocate_array(lfss,maxnks,maxnrs,maxizs,'LFSS',ierr)
    call allocate_array(lfts,maxnks,maxnrs,maxizs,'LFTS',ierr)
    call allocate_array(lllfps,maxnks,maxnrs,maxizs,'LLLFPS',ierr)

  end subroutine allocate_clocal


  subroutine deallocate_clocal
    implicit none

    call pr_trace('MOD_CLOCAL','DEALLOCATE')

    deallocate(ltolds)
    deallocate(lfps)
    deallocate(lfss)
    deallocate(lfts)
    deallocate(lllfps)

  end subroutine deallocate_clocal


end module mod_clocal
