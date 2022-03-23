module mod_solparams
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! slmod begin
  !...  needed for grid_d3d_10d:
  integer,public :: maxcases,mxspts=500
  !
  !      parameter (maxcases=10,mxspts=100)
  ! slmod end
  !
  ! jdemod - change mxspts from a parameter to a variable to allow code to change it
  !          this is needed in LIM solvers but making it a variable with an initialized value
  !          should not break code elsewhere
  !  
  !parameter (maxcases=10,mxspts=500)
  parameter (maxcases=10)

  integer,public :: maxiter
  real*8,public :: econv,mconv,eps
  parameter (econv=1.602192e-19,mconv=1.672614e-27,eps=1.0e-8)
  
  parameter (maxiter=1000)

  public :: allocate_mod_solparams,deallocate_mod_solparams

contains

  subroutine allocate_mod_solparams
    !use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_solparams','ALLOCATE')


  end subroutine allocate_mod_solparams


  subroutine deallocate_mod_solparams
    implicit none

    call pr_trace('mod_solparams','DEALLOCATE')


  end subroutine deallocate_mod_solparams

end module mod_solparams
