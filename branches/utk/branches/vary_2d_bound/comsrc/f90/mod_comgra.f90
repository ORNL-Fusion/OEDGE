module mod_comgra
  use debug_options
  implicit none

  !
  !     common block for graphics in out
  !
  !
  !
  !
  !     -*-fortran-*-
  ! common /comgra/ cxmin,cxmax,cymin,cymax,iplots,gra_nplots,ispot,thickness
  !
  ! save /comgra/
  real,public :: cxmin,cxmax,cymin,cymax
  integer,public :: iplots,gra_nplots,ispot
  !
  integer ,public :: thickness

  public :: allocate_mod_comgra,deallocate_mod_comgra

contains

  subroutine allocate_mod_comgra
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_comgra','ALLOCATE')


  end subroutine allocate_mod_comgra


  subroutine deallocate_mod_comgra
    implicit none

    call pr_trace('mod_comgra','DEALLOCATE')


  end subroutine deallocate_mod_comgra

end module mod_comgra