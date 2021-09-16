module mod_dynam5
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /dynam5/ cnvmf , cirng0 , cirng1 , cj0 , cj1 ,cvmf0 , cvmf1  , cvmf2
  ! save /dynam5/
  integer,public :: cnvmf
  integer,public,allocatable :: cirng0(:),cirng1(:),cj0(:),cj1(:)
  real,public,allocatable :: cvmf0(:),cvmf1(:),cvmf2(:)

  public :: allocate_mod_dynam5,deallocate_mod_dynam5

contains

  subroutine allocate_mod_dynam5
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam5','ALLOCATE')

    call allocate_array(cirng0,maxvmf,'cirng0',ierr)
    call allocate_array(cirng1,maxvmf,'cirng1',ierr)
    call allocate_array(cj0,maxvmf,'cj0',ierr)
    call allocate_array(cj1,maxvmf,'cj1',ierr)
    call allocate_array(cvmf0,maxvmf,'cvmf0',ierr)
    call allocate_array(cvmf1,maxvmf,'cvmf1',ierr)
    call allocate_array(cvmf2,maxvmf,'cvmf2',ierr)

  end subroutine allocate_mod_dynam5


  subroutine deallocate_mod_dynam5
    implicit none

    call pr_trace('mod_dynam5','DEALLOCATE')

    if (allocated(cirng0)) deallocate(cirng0)
    if (allocated(cirng1)) deallocate(cirng1)
    if (allocated(cj0)) deallocate(cj0)
    if (allocated(cj1)) deallocate(cj1)
    if (allocated(cvmf0)) deallocate(cvmf0)
    if (allocated(cvmf1)) deallocate(cvmf1)
    if (allocated(cvmf2)) deallocate(cvmf2)

  end subroutine deallocate_mod_dynam5

end module mod_dynam5