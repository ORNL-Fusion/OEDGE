module mod_crand

  use mod_params

!      REAL RANV(1000*ISECT)                                                     
!      COMMON /CRAND/ RANV                                                       


  implicit none
  private

      REAL,public:: RANV(1000*ISECT)                                                     

  
  public :: allocate_mod_crand, deallocate_mod_crand


contains

  subroutine allocate_mod_crand
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_crand


  subroutine deallocate_mod_crand
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_crand



end module mod_crand
