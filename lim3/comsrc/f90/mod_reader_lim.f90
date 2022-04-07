module mod_reader

  !use mod_params

!      !     -*-Fortran-*-
!      !                                                                       
!      COMMON /READR1/ BUFFER,BUFFER1                                    
!      COMMON /READR2/ IBUF                                              
!      CHARACTER*1024 BUFFER                                               
!      !CHARACTER*512 BUFFER                                               
!      CHARACTER*80 BUFFER1                                              
!      character*7,parameter :: buff_format = '(A1024)'
!      INTEGER IBUF                                                      



  implicit none
  private

      CHARACTER*1024,public:: BUFFER                                               
      CHARACTER*80,public::  BUFFER1                                              
      character*7,parameter,public :: buff_format = '(A1024)'
      INTEGER,public::  IBUF                                                      
 
  public :: allocate_mod_reader, deallocate_mod_reader


contains

  subroutine allocate_mod_reader
    !use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_reader


  subroutine deallocate_mod_reader
    !use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_reader



end module mod_reader
