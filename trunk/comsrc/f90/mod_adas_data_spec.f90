module mod_adas_data_spec
  use debug_options
  implicit none

  !
  !     this common block exists to allow the load_divdata_array
  !     rotuine to be reused without requiring the adas data specification
  !     to be reread from the input file.
  !
  !     cadas_switch is initialized to zero which deactivates this option
  !
  !     -*-fortran-*-
  ! common /adas_spec/ cisele,ciselr,ciselx,ciseld,cadasyr,           &cadas_switch,&
  !     cadasid,cadasex
  integer,public :: cisele,ciselr,ciselx,ciseld,cadasyr,cadas_switch
  character,public :: cadasid*80,cadasex*3
  
  ! save /adas_spec/

  public :: allocate_mod_adas_data_spec,deallocate_mod_adas_data_spec

contains

  subroutine allocate_mod_adas_data_spec
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_adas_data_spec','ALLOCATE')


  end subroutine allocate_mod_adas_data_spec


  subroutine deallocate_mod_adas_data_spec
    implicit none

    call pr_trace('mod_adas_data_spec','DEALLOCATE')


  end subroutine deallocate_mod_adas_data_spec

end module mod_adas_data_spec