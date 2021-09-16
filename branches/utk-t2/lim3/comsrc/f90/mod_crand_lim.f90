module mod_crand
  
  use mod_params
  
  !      real ranv(1000*isect)
  !      common /crand/ ranv
  
  
  implicit none
  private
  
  real,public,allocatable:: ranv(:)
  
  
  public :: allocate_mod_crand, deallocate_mod_crand
  
  
contains
  
  subroutine allocate_mod_crand
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(dtev  ,maxnxs,'dtev',ierr)


    call allocate_array(ranv,1000*isect,'ranv',ierr)

  end subroutine allocate_mod_crand
  
  
  subroutine deallocate_mod_crand
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

    if (allocated(ranv)) deallocate(ranv)

  end subroutine deallocate_mod_crand
  
  
  
end module mod_crand
