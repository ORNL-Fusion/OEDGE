module mod_global_options

  use mod_params

!c       -*-Fortran-*-
!c
!c	Common block for any options applying globally
!c
!
!	common /global_options/ cprint
!	integer cprint
!c
  

  implicit none
  private

	integer,public:: cprint

  
  public :: allocate_mod_global_options, deallocate_mod_global_options


contains

  subroutine allocate_mod_global_options
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_global_options


  subroutine deallocate_mod_global_options
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_global_options



end module mod_global_options
