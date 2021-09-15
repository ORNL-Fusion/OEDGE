program sol22_test
  use debug_options
  use mod_reader

  use mod_io_units
  !use mod_sol22
  !use mod_solswitch
  !use mod_solcommon
  !use mod_sol22_input
  use mod_calcsol_interface
  use mod_allocate_sol22_storage
  
  implicit none
  !integer ierr

  real*8 :: n0, te0, ti0,ringlen
  real*8 :: crmb,rizb
  integer :: npts

  call set_unit_numbers(in_stderr=0,in_stdin=5,in_stdout=6,&
                        in_stddbg=6,in_datunit=7,in_echout=9)

  call allocate_dynamic_input
  call allocate_sol22_storage
  
  
  !
  !      Set hard-coded global trace debugging options
  call init_trace(0,.true.)
  !      call init_trace(0,.false.)
  call pr_trace('SOL22_TEST','BEGIN EXECUTION')

!  open(stdin,file='sol22_input.txt',form='formatted',status='old')
!  rewind(stdin)
  
!10 read(stdin,buff_format,iostat=ierr,end=9999,err=9999) buffer
!  write(stderr,*) ':',trim(buffer),':'
!  goto 10

!  9999 rewind(stdin)

  

  n0 = 1.0e19
  te0 = 10.0
  ti0 = 10.0
  ringlen = 100
  npts = 400
  crmb = 2.0
  rizb = 1.0

  call calcsol_interface (n0,te0,ti0,ringlen,npts,crmb,rizb)

  call pr_trace('SOL22_TEST:','EXECUTION COMPLETE')

  call deallocate_sol22_storage

end program sol22_test


subroutine allocate_dynamic_input
  use mod_solswitch
  implicit none

    call allocate_mod_solswitch_input

end subroutine allocate_dynamic_input
