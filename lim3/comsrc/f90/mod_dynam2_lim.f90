module mod_dynam2


  !      COMMON /DYNAM2/ SDLIMS,SDLIM3,SDTS                                       !      
  !      REAL SDLIMS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                             ! 
  !      REAL SDLIM3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)              ! 
  !      REAL SDTS(MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                              

  implicit none

  REAL, allocatable, public ::  SDLIMS(:,:,:), SDLIM3(:,:,:,:), SDTS(:,:,:)                              

  public:: allocate_mod_dynam2,deallocate_mod_dynam2

  private

contains

  subroutine allocate_mod_dynam2
    use mod_params
    use allocate_arrays

    implicit none
    integer :: ierr

    call allocate_array(sdlims,1,maxnxs,-maxnys,maxnys,-1,maxizs,'SDLIMS',ierr)
    call allocate_array(sdlim3,1,maxnxs,-maxy3d,maxy3d,-1,maxizs,-maxnps,maxnps,'SDLIM3',ierr)
    call allocate_array(sdts,1,maxnxs,-maxnys,maxnys,1,maxizs,'SDTS',ierr)

  end subroutine allocate_mod_dynam2


  subroutine deallocate_mod_dynam2
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(sdlims)
    deallocate(sdlim3)
    deallocate(sdts)

  end subroutine deallocate_mod_dynam2

end module mod_dynam2
