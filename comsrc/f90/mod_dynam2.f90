module mod_dynam2

  use debug_options

  implicit none

  private

  !include 'dynam2'

  !     -*-Fortran-*-
  !                                                                       
  !      COMMON /DYNAM2/ SDLIMS,SDTS,PLRPS
  !      REAL SDLIMS(MAXNKS,MAXNRS,-1:MAXIZS)                              
  !      REAL SDTS(MAXNKS,MAXNRS,-1:MAXIZS)                                
  !      real PLRPS(MAXNKS,MAXNRS,-1:MAXPLRP)
  !      save /DYNAM2/
  !
  !      common /dynam2a/ SCHISQ1,SCHISQ2,SCHISQ3,SCHISQ4,SCHISQ5 
  !      save /dynam2a/
  !      real SCHISQ1(maxpiniter),SCHISQ2(maxpiniter),
  !     >     SCHISQ3(maxpiniter),
  !     >     SCHISQ4(maxpiniter),SCHISQ5(maxpiniter)
  !


  real, allocatable, public :: sdlims(:,:,:),sdts(:,:,:),plrps(:,:,:)
  real, allocatable, public :: schisq1(:),schisq2(:),schisq3(:),schisq4(:),schisq5(:)

  public :: allocate_dynam2,deallocate_dynam2

contains

  subroutine allocate_dynam2
    use global_parameters
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('MOD_DYNAM2','ALLOCATE')

    call allocate_array(sdlims,1,maxnks,1,maxnrs,-1,maxizs,'SDLIMS',ierr)
    call allocate_array(sdts,1,maxnks,1,maxnrs,-1,maxizs,'SDTS',ierr)
    call allocate_array(plrps,1,maxnks,1,maxnrs,-1,maxplrp,'PLRPS',ierr)

    call allocate_array(schisq1,maxpiniter,'SCHISQ1',ierr)
    call allocate_array(schisq2,maxpiniter,'SCHISQ2',ierr)
    call allocate_array(schisq3,maxpiniter,'SCHISQ3',ierr)
    call allocate_array(schisq4,maxpiniter,'SCHISQ4',ierr)
    call allocate_array(schisq5,maxpiniter,'SCHISQ5',ierr)


  end subroutine allocate_dynam2


  subroutine deallocate_dynam2
    implicit none

    call pr_trace('MOD_DYNAM2','DEALLOCATE')

    deallocate(sdlims)
    deallocate(sdts)
    deallocate(plrps)

    deallocate(schisq1)
    deallocate(schisq2)
    deallocate(schisq3)
    deallocate(schisq4)
    deallocate(schisq5)

  end subroutine deallocate_dynam2


end module mod_dynam2
