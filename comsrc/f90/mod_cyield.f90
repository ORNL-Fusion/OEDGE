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

    ! sazmod - changed from 7 to 8 to accomodate SiC options
    call allocate_array(ceth,8,21,'ceth',ierr)
    call allocate_array(cetf,8,21,'cetf',ierr)
    call allocate_array(cq,8,21,'cq',ierr)
    call allocate_array(cidata,8,21,'cidata',ierr)

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
