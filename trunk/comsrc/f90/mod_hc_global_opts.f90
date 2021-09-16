module mod_hc_global_opts
  use debug_options
  implicit none

  !
  ! jdemod
  !     all of the hc related variables have been included in modules -
  !     even the actual switch to turn the option on and off - the result
  !     of this is that hc related modules and code need to be included in
  !     a number of places where they are not required - to solve this i am
  !     creating global divimp versions of some of the options so the module
  !     inlusions are not required and the hc code impact is limited to the hc
  !     files.
  !
  !
  !     -*-fortran-*-
  ! common /global_hc_data/ global_hc_follow_option
  !
  ! save /global_hc_data/
  !
  integer,public :: global_hc_follow_option

  public :: allocate_mod_hc_global_opts,deallocate_mod_hc_global_opts

contains

  subroutine allocate_mod_hc_global_opts
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_hc_global_opts','ALLOCATE')


  end subroutine allocate_mod_hc_global_opts


  subroutine deallocate_mod_hc_global_opts
    implicit none

    call pr_trace('mod_hc_global_opts','DEALLOCATE')


  end subroutine deallocate_mod_hc_global_opts

end module mod_hc_global_opts