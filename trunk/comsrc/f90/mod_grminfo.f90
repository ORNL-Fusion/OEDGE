module mod_grminfo
  use debug_options
  implicit none

  !
  !     common block information for grm plots
  !
  !     -*-fortran-*-
  ! common /grminfo/ pageplots,textsize,axistextsize
  !
  ! save /grminfo/
  !
  integer,public :: pageplots,textsize,axistextsize

  public :: allocate_mod_grminfo,deallocate_mod_grminfo

contains

  subroutine allocate_mod_grminfo
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_grminfo','ALLOCATE')


  end subroutine allocate_mod_grminfo


  subroutine deallocate_mod_grminfo
    implicit none

    call pr_trace('mod_grminfo','DEALLOCATE')


  end subroutine deallocate_mod_grminfo

end module mod_grminfo