module mod_reader
  use debug_options
  implicit none

  !     -*-fortran-*-
  !
  ! common /readr1/ buffer,buffer1
  ! save /readr1/
  ! common /readr2/ ibuf
  
  ! save /readr2/
  character*1024,public :: buffer
  character*80,public :: buffer1
  character*7,parameter,public:: buff_format = '(a1024)'
  integer,public :: ibuf

  public :: allocate_mod_reader,deallocate_mod_reader

contains

  subroutine allocate_mod_reader
    !use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_reader','ALLOCATE')


  end subroutine allocate_mod_reader


  subroutine deallocate_mod_reader
    implicit none

    call pr_trace('mod_reader','DEALLOCATE')


  end subroutine deallocate_mod_reader

end module mod_reader
