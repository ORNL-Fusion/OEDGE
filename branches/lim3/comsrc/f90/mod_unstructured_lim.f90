module mod_unstructured

  use mod_params

!c
!c     This file will contain declarations for some unstructured 
!c     input values. If it becomes too unwieldy this file will 
!c     be split into separate common blocks for different 
!c     unstructured input values. 
!c



  implicit none
  private


  
 
  public :: allocate_mod_unstructured, deallocate_mod_unstructured


contains

  subroutine allocate_mod_unstructured
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_unstructured


  subroutine deallocate_mod_unstructured
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_unstructured



end module mod_unstructured
