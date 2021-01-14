module mod_dynam2


  !      COMMON /DYNAM2/ SDLIMS,SDLIM3,SDTS                                       !      
  !      REAL SDLIMS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                             ! 
  !      REAL SDLIM3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)              ! 
  !      REAL SDTS(MAXNXS,-MAXNYS:MAXNYS,MAXIZS)                              

  implicit none

  real, public :: qtim, fsrate
  
  REAL, allocatable, public ::  SDLIMS(:,:,:), SDLIM3(:,:,:,:), SDTS(:,:,:)
  real, allocatable, public ::  sdvs(:,:,:), sdtimp(:,:,:)
  
  real, allocatable, public :: vtig_array(:,:,:)
  
  logical,public :: debugv = .false.

  

  public:: allocate_mod_dynam2,deallocate_mod_dynam2,allocate_debugv,deallocate_debugv

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

  subroutine allocate_debugv
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr
    
    
    if (debugv) then

       call allocate_array(sdvs,1,maxnxs,-maxnys,maxnys,1,maxizs,'sdvs',ierr)
       call allocate_array(sdtimp,1,maxnxs,-maxnys,maxnys,1,maxizs,'sdtimp',ierr)
       call allocate_array(vtig_array,1,maxnxs,-maxnys,maxnys,1,maxpzone,'vtig_array',ierr)

    endif
    

  end subroutine allocate_debugv

  subroutine deallocate_debugv
    implicit none

    if (debugv) then

       if (allocated(sdvs)) deallocate(sdvs)
       if (allocated(sdtimp)) deallocate(sdtimp)
       if (allocated(vtig_array)) deallocate(vtig_array)

    endif 
   

  end subroutine deallocate_debugv


  
end module mod_dynam2
