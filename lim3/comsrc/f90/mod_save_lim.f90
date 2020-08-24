module mod_save
  
  use mod_params
  
  
  !c
  !      common /save/ r,i,l
  !      real          r(9,maxput)
  !      integer       i(9,maxput)
  !      logical       l(3,maxput)
  
  implicit none
  private
  
  real,public,allocatable:: r(:,:)
  integer,public,allocatable:: i(:,:)
  logical,public,allocatable:: l(:,:)
  
  public :: allocate_mod_save, deallocate_mod_save
  
  
contains
  
  subroutine allocate_mod_save
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(dtev  ,maxnxs,'dtev',ierr)


    call allocate_array(r,9,maxput,'r',ierr)
    call allocate_array(i,9,maxput,'i',ierr)
    call allocate_array(l,3,maxput,'l',ierr)

  end subroutine allocate_mod_save
  
  
  subroutine deallocate_mod_save
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

    if (allocated(r)) deallocate(r)
    if (allocated(i)) deallocate(i)
    if (allocated(l)) deallocate(l)

  end subroutine deallocate_mod_save
  
  
  
end module mod_save
