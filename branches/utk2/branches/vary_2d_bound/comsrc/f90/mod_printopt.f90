module mod_printopt
  use debug_options
  implicit none

  !
  !     this common block contains options specifically related to
  !     formatting the divimp printouts
  !
  !     -*-fortran-*-
  ! common /printopts/ prsol21,prsol22,prsol28,coment,inner,outer
  !
  ! save /printopts/
  logical,public :: prsol21,prsol22,prsol28
  !
  character,public :: coment*120,inner*5,outer*5,sp*24,sa*2,s1*6,s2*12
  data sp /'                      '/
  data sa /'  '/
  data s1 /'      '/
  !
  data s2 /'          '/

  public :: allocate_mod_printopt,deallocate_mod_printopt

contains

  subroutine allocate_mod_printopt
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_printopt','ALLOCATE')


  end subroutine allocate_mod_printopt


  subroutine deallocate_mod_printopt
    implicit none

    call pr_trace('mod_printopt','DEALLOCATE')


  end subroutine deallocate_mod_printopt

end module mod_printopt