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
  integer,public :: pageplots,textsize,axistextsize,write_grm_data,iout_grm
  real, public :: absfac_grm_copy

  
  public :: allocate_mod_grminfo,deallocate_mod_grminfo

contains

  subroutine allocate_mod_grminfo
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_grminfo','ALLOCATE')
    write_grm_data = 0
    iout_grm = 26
    

  end subroutine allocate_mod_grminfo


  subroutine deallocate_mod_grminfo
    implicit none

    call pr_trace('mod_grminfo','DEALLOCATE')


  end subroutine deallocate_mod_grminfo

end module mod_grminfo

