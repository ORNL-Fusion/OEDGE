module mod_cyield
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /cyield/ ceth,cetf,cq,cidata,ntars,flux_frac
  ! save /cyield/
  real,public :: flux_frac
  real,public,allocatable :: ceth(:,:),cetf(:,:),cq(:,:)
  integer,public :: ntars
  logical,public,allocatable :: cidata(:,:)

  public :: allocate_mod_cyield,deallocate_mod_cyield

contains

  subroutine allocate_mod_cyield
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cyield','ALLOCATE')

    call allocate_array(ceth,7,12,'ceth',ierr)
    call allocate_array(cetf,7,12,'cetf',ierr)
    call allocate_array(cq,7,12,'cq',ierr)
    call allocate_array(cidata,7,12,'cidata',ierr)

  end subroutine allocate_mod_cyield


  subroutine deallocate_mod_cyield
    implicit none

    call pr_trace('mod_cyield','DEALLOCATE')

    if (allocated(ceth)) deallocate(ceth)
    if (allocated(cetf)) deallocate(cetf)
    if (allocated(cq)) deallocate(cq)
    if (allocated(cidata)) deallocate(cidata)

  end subroutine deallocate_mod_cyield

end module mod_cyield